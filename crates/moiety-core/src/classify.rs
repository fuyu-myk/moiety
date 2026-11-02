//! Classifies what separates two readings into the ambiguity taxonomy (grammar, section 7.2).

use std::collections::{BTreeMap, BTreeSet};

use crate::interpret::Candidate;
use crate::lexer::{EdgeKind, Lattice, RegionKind};
use crate::molecule::Molecule;
use crate::reading::AmbiguityClass;
use crate::span::Span;

/// Every choice made by any derivation of one molecule.
pub(crate) struct Choices {
    regions: Vec<BTreeSet<Vec<usize>>>,
    count_as_charge: BTreeSet<Option<Span>>,
    forks: BTreeMap<Span, BTreeSet<&'static str>>,
    charges: BTreeSet<Vec<Span>>,
}

impl Choices {
    pub fn of(members: &[Candidate<'_>], lattice: &Lattice) -> Self {
        let mut regions = vec![BTreeSet::new(); lattice.regions().len()];
        let mut count_as_charge = BTreeSet::new();
        let mut forks: BTreeMap<Span, BTreeSet<&'static str>> = BTreeMap::new();
        let mut charges = BTreeSet::new();
        for member in members {
            for (r, region) in lattice.regions().iter().enumerate() {
                let inside: Vec<usize> = member
                    .tree
                    .path
                    .iter()
                    .copied()
                    .filter(|&e| {
                        let span = lattice.edges()[e].span;
                        region.span.start <= span.start && span.end <= region.span.end
                    })
                    .collect();
                regions[r].insert(inside);
            }
            count_as_charge.insert(member.tree.count_as_charge);
            for &(span, choice) in &member.derivation.forks {
                forks.entry(span).or_default().insert(choice);
            }
            charges.insert(member.derivation.charge_sites.clone());
        }
        Self {
            regions,
            count_as_charge,
            forks,
            charges,
        }
    }
}

/// One ambiguity separating two readings.
pub(crate) struct Site {
    pub class: AmbiguityClass,
    pub span: Span,
    pub top: String,
    pub other: String,
}

/// What separates `other` from `top`, by the first rule of section 7.2 that applies.
pub(crate) fn classify(
    top: &Choices,
    other: &Choices,
    top_best: &Candidate<'_>,
    other_best: &Candidate<'_>,
    lattice: &Lattice,
    text: &str,
) -> Vec<Site> {
    let lexical: Vec<Site> = lattice
        .regions()
        .iter()
        .enumerate()
        .filter(|(r, _)| top.regions[*r].is_disjoint(&other.regions[*r]))
        .map(|(r, region)| Site {
            class: match region.kind {
                RegionKind::Segmentation => AmbiguityClass::SymbolSegmentation,
                RegionKind::Shadowing => AmbiguityClass::AbbreviationShadowing,
                RegionKind::Dash => AmbiguityClass::ChargeVsAttachment,
            },
            span: region.span,
            top: describe(first(&top.regions[r]), lattice, text),
            other: describe(first(&other.regions[r]), lattice, text),
        })
        .collect();
    if !lexical.is_empty() {
        return lexical;
    }
    if top.count_as_charge.is_disjoint(&other.count_as_charge) {
        let span = top
            .count_as_charge
            .iter()
            .chain(&other.count_as_charge)
            .find_map(|s| *s)
            .unwrap_or(top_best.tree.span);
        let name = |set: &BTreeSet<Option<Span>>| {
            if set.contains(&None) { "subscript" } else { "charge magnitude" }.to_string()
        };
        return vec![Site {
            class: AmbiguityClass::ChargeRepresentation,
            span,
            top: name(&top.count_as_charge),
            other: name(&other.count_as_charge),
        }];
    }
    let forks: Vec<Site> = top
        .forks
        .iter()
        .filter_map(|(span, mine)| {
            let theirs = other.forks.get(span)?;
            mine.is_disjoint(theirs).then(|| Site {
                class: AmbiguityClass::BranchPointFork,
                span: *span,
                top: join(mine),
                other: join(theirs),
            })
        })
        .collect();
    if !forks.is_empty() {
        return forks;
    }
    if top.charges.is_disjoint(&other.charges) {
        let span = first(&top.charges).first().copied().unwrap_or(top_best.tree.span);
        return vec![Site {
            class: AmbiguityClass::ChargeRepresentation,
            span,
            top: charge_text(first(&top.charges), text),
            other: charge_text(first(&other.charges), text),
        }];
    }
    vec![topology(top_best, other_best, text)]
}

fn first<T: Default + Clone>(set: &BTreeSet<T>) -> T {
    set.iter().next().cloned().unwrap_or_default()
}

fn join(set: &BTreeSet<&'static str>) -> String {
    set.iter().copied().collect::<Vec<_>>().join(" or ")
}

fn charge_text(sites: Vec<Span>, text: &str) -> String {
    let atoms: Vec<&str> = sites.iter().map(|s| &text[s.start..s.end]).collect();
    format!("charge on {}", atoms.join(", "))
}

fn describe(edges: Vec<usize>, lattice: &Lattice, text: &str) -> String {
    edges
        .iter()
        .map(|&e| {
            let edge = lattice.edges()[e];
            let written = &text[edge.span.start..edge.span.end];
            match edge.kind {
                EdgeKind::Element { element, .. } => element.symbol().to_string(),
                EdgeKind::Abbreviation { .. } => format!("abbreviation {written}"),
                EdgeKind::Sign(_) | EdgeKind::Charge(_) => "charge".to_string(),
                EdgeKind::Marker(_) => "attachment".to_string(),
                _ => written.to_string(),
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}

/// The first atom, by input position, whose valence differs between the two readings.
fn topology(top: &Candidate<'_>, other: &Candidate<'_>, text: &str) -> Site {
    let signatures = |c: &Candidate<'_>| -> BTreeMap<Span, (String, u8)> {
        let molecule: &Molecule = &c.derivation.molecule;
        let mut map = BTreeMap::new();
        for (i, atom) in molecule.atoms().iter().enumerate() {
            let Some(element) = atom.element() else { continue };
            let valence = molecule.bond_valence(i) + atom.hydrogens;
            map.entry(c.derivation.atom_spans[i])
                .or_insert((element.symbol().to_string(), valence));
        }
        map
    };
    let (mine, theirs) = (signatures(top), signatures(other));
    for (span, (symbol, valence)) in &mine {
        if let Some((_, their_valence)) = theirs.get(span).filter(|t| t.1 != *valence) {
            return Site {
                class: AmbiguityClass::BranchPointFork,
                span: *span,
                top: format!("{symbol} with valence {valence}"),
                other: format!("{symbol} with valence {their_valence}"),
            };
        }
    }
    let span = top.tree.span;
    Site {
        class: AmbiguityClass::BranchPointFork,
        span,
        top: text[span.start..span.end].to_string(),
        other: text[span.start..span.end].to_string(),
    }
}
