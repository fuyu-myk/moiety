use std::cell::Cell;
use std::collections::BTreeMap;

use crate::assembly::{Derivation, Failure, assemble, assemble_salt};
use crate::classify::{Choices, classify};
use crate::diagnostic::{Code, Diagnostic, Label};
use crate::dictionary::Dictionary;
use crate::element::Element;
use crate::lexer::{EdgeKind, Lattice, RegionKind, lex};
use crate::molecule::{BondOrder, Molecule};
use crate::normalize::{Normalized, normalize};
use crate::parser::{SyntaxError, Tail, Tree, parse};
use crate::reading::{
    Ambiguity, AmbiguityClass, Direction, Interpretation, Provenance, Reading, Underdetermination,
};
use crate::span::Span;
use crate::suggest::suggest;

/// Version of the grammar this crate implements.
pub const GRAMMAR_VERSION: u32 = 1;
/// Version of the ranking policy this crate implements.
pub const RANKING_POLICY_VERSION: u32 = 2;
/// Longest formula, in bytes, that [`interpret`] accepts.
pub const MAX_INPUT_LEN: usize = 256;

const PATH_LIMIT: usize = 128;

/// How hypervalent nitrogen is represented. See the grammar, section 6.10.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum ChargeConvention {
    /// `[N+](=O)[O-]`.
    #[default]
    ChargeSeparated,
    /// `N(=O)=O`.
    Hypervalent,
}

impl ChargeConvention {
    /// Name as used in output.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            ChargeConvention::ChargeSeparated => "charge-separated",
            ChargeConvention::Hypervalent => "hypervalent",
        }
    }
}

/// Settings for [`interpret`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Options {
    /// Representation of hypervalent nitrogen.
    pub charge_convention: ChargeConvention,
    /// Most readings returned, never fewer than two. Later readings are dropped with warning CF0303.
    pub max_readings: usize,
    /// Whether invalid input triggers a search for a nearby formula that reads.
    pub suggestions: bool,
    /// Most search steps spent on one component before giving up with CF0003. A
    /// step is a search state or an atom copied by a group count.
    pub search_limit: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            charge_convention: ChargeConvention::default(),
            max_readings: 32,
            suggestions: true,
            search_limit: 200_000,
        }
    }
}

/// Reads a condensed formula and returns every reading that survives valence resolution.
#[must_use]
pub fn interpret(input: &str, dictionary: &Dictionary, options: &Options) -> Interpretation {
    interpret_within(input, dictionary, options, None)
}

/// [`interpret`], drawing search states from `shared` instead of a fresh budget per component.
pub(crate) fn interpret_within(
    input: &str,
    dictionary: &Dictionary,
    options: &Options,
    shared: Option<&Cell<usize>>,
) -> Interpretation {
    let normalized = normalize(input);
    let whole = Span::new(0, input.len());
    if input.len() > MAX_INPUT_LEN {
        return invalid(Diagnostic::new(
            Code::TooLong,
            format!("formula is longer than {MAX_INPUT_LEN} bytes"),
            Label::new(whole, "too long to read"),
        ));
    }
    let text = normalized.text();
    let mut components = Vec::new();
    let mut start = 0;
    for (i, c) in text.char_indices() {
        if c == '.' {
            components.push(Span::new(start, i));
            start = i + 1;
        }
    }
    components.push(Span::new(start, text.len()));
    if let Some(empty) = components.iter().find(|span| span.is_empty()) {
        let at = Span::new(
            empty.start.saturating_sub(1),
            (empty.end + 1).min(text.len()),
        );
        return invalid(Diagnostic::new(
            Code::Empty,
            "empty component",
            Label::new(
                if text.is_empty() {
                    whole
                } else {
                    normalized.to_original(at)
                },
                "nothing to read here",
            ),
        ));
    }
    let mut outcomes = Vec::new();
    for span in components {
        match component(input, &normalized, span, dictionary, options, shared) {
            Outcome::Readings(readings) => outcomes.push(readings),
            Outcome::Underdetermined(span, reason) => {
                return Interpretation::Underdetermined { span, reason };
            }
            Outcome::Invalid(diagnostic) => return invalid(*diagnostic),
        }
    }
    combine(outcomes, &normalized, options)
}

fn invalid(diagnostic: Diagnostic) -> Interpretation {
    Interpretation::Invalid {
        diagnostics: vec![diagnostic],
    }
}

enum Outcome {
    Readings(Vec<Partial>),
    Underdetermined(Span, Underdetermination),
    Invalid(Box<Diagnostic>),
}

/// A reading of one component, before components are combined.
struct Partial {
    molecule: Molecule,
    atom_spans: Vec<Span>,
    bond_spans: Vec<Span>,
    smiles: String,
    score: i32,
    direction: Direction,
    ambiguities: Vec<Ambiguity>,
    classes: Vec<AmbiguityClass>,
    warnings: Vec<Diagnostic>,
}

/// One derivation with its tree and total score.
pub(crate) struct Candidate<'a> {
    pub tree: &'a Tree,
    pub derivation: Derivation,
    pub molecule: Molecule,
    pub score: i32,
    pub conventions: Vec<Span>,
}

fn component(
    input: &str,
    normalized: &Normalized,
    span: Span,
    dictionary: &Dictionary,
    options: &Options,
    shared: Option<&Cell<usize>>,
) -> Outcome {
    let original = |s: Span| normalized.to_original(s);
    let lattice = match lex(normalized.text(), span, dictionary) {
        Ok(lattice) => lattice,
        Err(uncovered) if uncovered.looks_like_abbreviation => {
            return Outcome::Underdetermined(
                original(uncovered.span),
                Underdetermination::UnknownAbbreviation,
            );
        }
        Err(uncovered) => {
            return Outcome::Invalid(Box::new(with_suggestion(
                Diagnostic::new(
                    Code::UnrecognizedCharacter,
                    "unrecognized character",
                    Label::new(
                        original(uncovered.span),
                        "no element, abbreviation or symbol starts here",
                    ),
                ),
                input,
                dictionary,
                options,
            )));
        }
    };
    let (paths, paths_truncated) = lattice.paths(PATH_LIMIT);
    let mut trees = Vec::new();
    let mut syntax: Option<SyntaxError> = None;
    for path in &paths {
        match parse(&lattice, path) {
            Ok(parsed) => trees.extend(parsed),
            Err(error) => {
                if syntax
                    .as_ref()
                    .is_none_or(|s| error.span.start > s.span.start)
                {
                    syntax = Some(error);
                }
            }
        }
    }
    if trees.is_empty() {
        let error = syntax.expect("a lattice has at least one path");
        return Outcome::Invalid(Box::new(with_suggestion(
            Diagnostic::new(
                error.code,
                title(error.code),
                Label::new(original(error.span), error.message),
            ),
            input,
            dictionary,
            options,
        )));
    }
    let lexical: Vec<i32> = trees
        .iter()
        .map(|t| lexical_score(&lattice, normalized.text(), t, dictionary))
        .collect();
    let mut candidates = Vec::new();
    let mut failure: Option<(Failure, i32)> = None;
    let mut underdetermined = None;
    let mut truncated = paths_truncated;
    let mut exhausted = false;
    let local = Cell::new(options.search_limit);
    let budget = shared.unwrap_or(&local);
    for salt in [false, true] {
        for (tree, &score) in trees.iter().zip(&lexical) {
            let assembly = if salt {
                assemble_salt(tree, dictionary, budget)
            } else {
                assemble(tree, dictionary, budget)
            };
            truncated |= assembly.truncated || assembly.exhausted;
            exhausted |= assembly.exhausted;
            if let Some(found) = assembly.failure {
                let better = failure.as_ref().is_none_or(|(current, current_score)| {
                    (found.progress, score) > (current.progress, *current_score)
                });
                if better {
                    failure = Some((found, score));
                }
            }
            if underdetermined.is_none() {
                underdetermined = assembly.underdetermined;
            }
            for derivation in assembly.derivations {
                candidates.push(candidate(tree, derivation, score, options));
            }
        }
        if !candidates.is_empty() {
            break;
        }
    }
    if candidates.is_empty() {
        if exhausted {
            return Outcome::Invalid(Box::new(Diagnostic::new(
                Code::TooComplex,
                "formula needs more search than the limit allows",
                Label::new(original(span), "too many alternative readings to explore"),
            )));
        }
        if trees.iter().any(Tree::is_formula_shaped) {
            return Outcome::Underdetermined(original(span), Underdetermination::MolecularFormula);
        }
        if let Some((at, reason)) = underdetermined {
            return Outcome::Underdetermined(original(at), reason);
        }
        let failure = failure.map_or_else(
            || {
                Failure::new(
                    Code::OpenValence,
                    "no reading satisfies every valence",
                    span,
                    "no reading".into(),
                    span.end,
                )
            },
            |(failure, _)| failure,
        );
        return Outcome::Invalid(Box::new(with_suggestion(
            to_diagnostic(&failure, normalized),
            input,
            dictionary,
            options,
        )));
    }
    Outcome::Readings(rank(candidates, &lattice, normalized, truncated, options))
}

fn candidate<'a>(
    tree: &'a Tree,
    derivation: Derivation,
    lexical: i32,
    options: &Options,
) -> Candidate<'a> {
    let mut molecule = derivation.molecule.clone();
    let changed = apply_convention(&mut molecule, options.charge_convention);
    let conventions = changed.iter().map(|&i| derivation.atom_spans[i]).collect();
    let single = molecule.clone();
    for _ in 1..tree.multiplier {
        molecule.extend(&single);
    }
    Candidate {
        tree,
        score: lexical + derivation.score,
        derivation,
        molecule,
        conventions,
    }
}

fn rank(
    candidates: Vec<Candidate<'_>>,
    lattice: &Lattice,
    normalized: &Normalized,
    truncated: bool,
    options: &Options,
) -> Vec<Partial> {
    let mut groups: BTreeMap<String, Vec<Candidate<'_>>> = BTreeMap::new();
    for candidate in candidates {
        groups
            .entry(candidate.molecule.to_smiles())
            .or_default()
            .push(candidate);
    }
    let mut groups: Vec<(String, Vec<Candidate<'_>>)> = groups.into_iter().collect();
    for (_, members) in &mut groups {
        members.sort_by_key(|c| (-c.score, c.tree.path.len()));
    }
    groups.sort_by(|(a_smiles, a), (b_smiles, b)| {
        b[0].score
            .cmp(&a[0].score)
            .then_with(|| a_smiles.cmp(b_smiles))
    });
    let dropped = groups.len() > reading_limit(options);
    groups.truncate(reading_limit(options));
    let choices: Vec<Choices> = groups
        .iter()
        .map(|(_, members)| Choices::of(members, lattice))
        .collect();
    let mut ambiguities: Vec<Vec<Ambiguity>> = vec![Vec::new(); groups.len()];
    let mut classes: Vec<Vec<AmbiguityClass>> = vec![Vec::new(); groups.len()];
    for first in 0..groups.len() {
        for second in first + 1..groups.len() {
            let sites = classify(
                &choices[first],
                &choices[second],
                &groups[first].1[0],
                &groups[second].1[0],
                lattice,
                normalized.text(),
            );
            for site in sites {
                push_class(&mut classes[first], site.class);
                push_class(&mut classes[second], site.class);
                let span = normalized.to_original(site.span);
                push_unique(
                    &mut ambiguities[first],
                    Ambiguity {
                        class: site.class,
                        span,
                        choice: site.top,
                    },
                );
                push_unique(
                    &mut ambiguities[second],
                    Ambiguity {
                        class: site.class,
                        span,
                        choice: site.other,
                    },
                );
            }
        }
    }
    groups
        .into_iter()
        .zip(ambiguities)
        .zip(classes)
        .map(|(((smiles, members), mut ambiguities), classes)| {
            let best = members.into_iter().next().expect("groups are non-empty");
            for &span in &best.conventions {
                push_unique(
                    &mut ambiguities,
                    Ambiguity {
                        class: AmbiguityClass::ChargeRepresentation,
                        span: normalized.to_original(span),
                        choice: options.charge_convention.as_str().to_string(),
                    },
                );
            }
            let mut warnings: Vec<Diagnostic> = best
                .tree
                .warnings
                .iter()
                .map(|&(code, span)| warning(code, normalized.to_original(span)))
                .collect();
            if truncated || dropped {
                warnings.push(warning(
                    Code::ReadingLimit,
                    normalized.to_original(best.tree.span),
                ));
            }
            let d = &best.derivation;
            let direction = match (d.forward, d.backward) {
                (_, 0) => Direction::LeftToRight,
                (0, _) => Direction::RightToLeft,
                _ => Direction::Mixed,
            };
            let copies = best.tree.multiplier as usize;
            Partial {
                atom_spans: d.atom_spans.repeat(copies),
                bond_spans: d.bond_spans.repeat(copies),
                molecule: best.molecule,
                smiles,
                score: best.score,
                direction,
                ambiguities,
                classes,
                warnings,
            }
        })
        .collect()
}

fn push_class(list: &mut Vec<AmbiguityClass>, class: AmbiguityClass) {
    if !list.contains(&class) {
        list.push(class);
    }
}

fn push_unique(list: &mut Vec<Ambiguity>, ambiguity: Ambiguity) {
    if !list
        .iter()
        .any(|a| a.class == ambiguity.class && a.span == ambiguity.span)
    {
        list.push(ambiguity);
    }
}

fn warning(code: Code, span: Span) -> Diagnostic {
    let (message, label) = match code {
        Code::StereoIgnored => (
            "stereo descriptor ignored",
            "stereochemistry is out of scope",
        ),
        Code::IsotopeIgnored => (
            "isotope label read as hydrogen",
            "isotopes are out of scope",
        ),
        Code::MixedBrackets => (
            "group closed with a different bracket",
            "brackets of both kinds",
        ),
        _ => ("reading limit reached", "later readings were dropped"),
    };
    Diagnostic::new(code, message, Label::new(span, label))
}

fn combine(
    outcomes: Vec<Vec<Partial>>,
    normalized: &Normalized,
    options: &Options,
) -> Interpretation {
    let mut combined: Vec<Partial> = vec![Partial {
        molecule: Molecule::new(),
        atom_spans: Vec::new(),
        bond_spans: Vec::new(),
        smiles: String::new(),
        score: 0,
        direction: Direction::LeftToRight,
        ambiguities: Vec::new(),
        classes: Vec::new(),
        warnings: Vec::new(),
    }];
    let mut dropped = false;
    for readings in outcomes {
        let mut next = Vec::new();
        for base in &combined {
            for reading in &readings {
                let mut molecule = base.molecule.clone();
                molecule.extend(&reading.molecule);
                let mut ambiguities = base.ambiguities.clone();
                for a in &reading.ambiguities {
                    push_unique(&mut ambiguities, a.clone());
                }
                let direction = if base.molecule.is_empty() || base.direction == reading.direction {
                    reading.direction
                } else {
                    Direction::Mixed
                };
                next.push(Partial {
                    smiles: String::new(),
                    atom_spans: [base.atom_spans.clone(), reading.atom_spans.clone()].concat(),
                    bond_spans: [base.bond_spans.clone(), reading.bond_spans.clone()].concat(),
                    molecule,
                    score: base.score + reading.score,
                    direction,
                    ambiguities,
                    classes: [base.classes.clone(), reading.classes.clone()].concat(),
                    warnings: [base.warnings.clone(), reading.warnings.clone()].concat(),
                });
            }
        }
        if next.len() > reading_limit(options) {
            dropped = true;
            next.sort_by_key(|p| -p.score);
            next.truncate(reading_limit(options));
        }
        combined = next;
    }
    for partial in &mut combined {
        partial.smiles = partial.molecule.to_smiles();
    }
    combined.sort_by(|a, b| b.score.cmp(&a.score).then_with(|| a.smiles.cmp(&b.smiles)));
    let combined = merge_duplicates(combined);
    let mut classes: Vec<AmbiguityClass> = combined
        .iter()
        .flat_map(|p| p.classes.iter().copied())
        .collect();
    classes.sort_unstable();
    classes.dedup();
    let readings: Vec<Reading> = combined
        .into_iter()
        .enumerate()
        .map(|(i, mut p)| {
            if dropped && !p.warnings.iter().any(|w| w.code == Code::ReadingLimit) {
                p.warnings
                    .push(warning(Code::ReadingLimit, Span::new(0, 0)));
            }
            Reading {
                molecule: p.molecule,
                smiles: p.smiles,
                rank: u32::try_from(i + 1).unwrap_or(u32::MAX),
                score: p.score,
                direction: p.direction,
                ambiguities: p.ambiguities,
                provenance: Provenance::new(
                    p.atom_spans,
                    p.bond_spans,
                    normalized.folds().to_vec(),
                ),
                warnings: p.warnings,
            }
        })
        .collect();
    if readings.len() == 1 {
        Interpretation::Unique(readings.into_iter().next().expect("one reading"))
    } else {
        Interpretation::Ambiguous { readings, classes }
    }
}

/// Most readings returned: the option, but never fewer than two, so an ambiguous
/// result always shows an alternative to its first reading.
fn reading_limit(options: &Options) -> usize {
    options.max_readings.max(2)
}

/// Folds readings with the same SMILES into the best-scoring one, which must come first.
fn merge_duplicates(sorted: Vec<Partial>) -> Vec<Partial> {
    let mut merged: Vec<Partial> = Vec::new();
    for partial in sorted {
        match merged.iter_mut().find(|m| m.smiles == partial.smiles) {
            Some(kept) => {
                for ambiguity in partial.ambiguities {
                    push_unique(&mut kept.ambiguities, ambiguity);
                }
                for class in partial.classes {
                    push_class(&mut kept.classes, class);
                }
            }
            None => merged.push(partial),
        }
    }
    for partial in &mut merged {
        let mut seen = Vec::new();
        partial.warnings.retain(|w| {
            let key = (w.code, w.primary.span);
            let fresh = !seen.contains(&key);
            seen.push(key);
            fresh
        });
    }
    merged
}

fn lexical_score(lattice: &Lattice, text: &str, tree: &Tree, dictionary: &Dictionary) -> i32 {
    let edges = lattice.edges();
    let mut priors: i32 = tree.path.iter().map(|&i| edges[i].prior).sum();
    if let [.., count, sign] = tree.path[..] {
        let (count, sign) = (edges[count], edges[sign]);
        if matches!(count.kind, EdgeKind::Count(_)) && matches!(sign.kind, EdgeKind::Sign(_)) {
            let after_hydrogen = text[..count.span.start].ends_with('H');
            priors -= i32::from(tree.count_as_charge.is_some() == after_hydrogen);
        }
    }
    let dash = lattice
        .regions()
        .iter()
        .find(|r| r.kind == RegionKind::Dash)
        .map(|r| r.span);
    let Some(dash) = dash else {
        return priors;
    };
    let as_charge = matches!(tree.tail, Tail::Charge(..));
    let prefers_marker = last_symbol_is_carbon(&text[lattice.span().start..dash.start], dictionary);
    priors - i32::from(as_charge == prefers_marker)
}

/// Whether the atom a trailing `-` would sit on is carbon: the last element symbol
/// before it, skipping hydrogens, counts and whole groups, or the attachment atom
/// of a dictionary spelling that ends there.
fn last_symbol_is_carbon(prefix: &str, dictionary: &Dictionary) -> bool {
    let mut rest = prefix;
    loop {
        rest = rest.trim_end_matches(|c: char| c.is_ascii_digit() || c == 'H');
        let Some(last) = rest.chars().last() else {
            return false;
        };
        if !matches!(last, ')' | ']') {
            break;
        }
        let mut depth = 0usize;
        let open = rest.char_indices().rev().find(|&(_, c)| {
            match c {
                ')' | ']' => depth += 1,
                '(' | '[' => depth = depth.saturating_sub(1),
                _ => {}
            }
            depth == 0
        });
        match open {
            Some((at, _)) => rest = &rest[..at],
            None => return false,
        }
    }
    let spelled = dictionary.entries().iter().find_map(|entry| {
        entry
            .spellings()
            .iter()
            .any(|spelling| rest.ends_with(spelling.as_str()))
            .then(|| entry.attachments().first())
            .flatten()
            .and_then(|site| entry.fragment().atoms()[site.atom].element())
    });
    if let Some(element) = spelled {
        return element == Element::C;
    }
    rest.rfind(|c: char| c.is_ascii_uppercase())
        .is_some_and(|start| &rest[start..] == "C")
}

/// Applies the charge convention and returns the indices of nitrogen atoms it changed.
fn apply_convention(molecule: &mut Molecule, convention: ChargeConvention) -> Vec<usize> {
    let adjacency = molecule.adjacency();
    let mut changed = Vec::new();
    for n in 0..molecule.atoms().len() {
        let atom = molecule.atoms()[n];
        if atom.element() != Some(Element::N) || atom.aromatic {
            continue;
        }
        let total = molecule.bond_valence(n) + atom.hydrogens;
        let (want_charge, from, to, other_charge) = match convention {
            ChargeConvention::ChargeSeparated if atom.charge == 0 && total == 5 => {
                (1, BondOrder::Double, BondOrder::Single, -1)
            }
            ChargeConvention::Hypervalent if atom.charge == 1 && total == 4 => {
                (0, BondOrder::Single, BondOrder::Double, 0)
            }
            _ => continue,
        };
        let oxygen = adjacency[n].iter().copied().find(|&(o, bond)| {
            let other = molecule.atoms()[o];
            other.element() == Some(Element::O)
                && adjacency[o].len() == 1
                && other.hydrogens == 0
                && molecule.bonds()[bond].order == from
                && other.charge == if other_charge == -1 { 0 } else { -1 }
        });
        if let Some((o, bond)) = oxygen {
            molecule.set_bond_order(bond, to);
            molecule.atoms_mut()[n].charge = want_charge;
            molecule.atoms_mut()[o].charge = other_charge;
            changed.push(n);
        }
    }
    changed
}

fn title(code: Code) -> &'static str {
    match code {
        Code::Empty => "empty component",
        Code::UnclosedGroup => "unclosed group",
        Code::UnmatchedClose => "unmatched closing bracket",
        Code::EmptyGroup => "empty group",
        Code::CountOutOfRange => "count or charge out of range",
        _ => "unexpected token",
    }
}

fn to_diagnostic(failure: &Failure, normalized: &Normalized) -> Diagnostic {
    let mut diagnostic = Diagnostic::new(
        failure.code,
        failure.message.clone(),
        Label::new(
            normalized.to_original(failure.primary.0),
            failure.primary.1.clone(),
        ),
    );
    for (span, label) in &failure.secondary {
        diagnostic =
            diagnostic.with_label(Label::new(normalized.to_original(*span), label.clone()));
    }
    diagnostic
}

fn with_suggestion(
    diagnostic: Diagnostic,
    input: &str,
    dictionary: &Dictionary,
    options: &Options,
) -> Diagnostic {
    if !options.suggestions {
        return diagnostic;
    }
    match suggest(input, diagnostic.primary.span, dictionary, options) {
        Some(fixed) => diagnostic.with_note(format!("did you mean {fixed}?")),
        None => diagnostic,
    }
}
