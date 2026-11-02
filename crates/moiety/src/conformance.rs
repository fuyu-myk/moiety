//! The conformance runner.
//!
//! The suite is a versioned TOML file of inputs with their expected
//! interpretations (see `conformance/cases.toml` and the grammar, section 10).
//! The runner scores an implementation against it and against a round-trip
//! corpus of molecules written by the inverse writer and read back.

use std::fmt;

use moiety_core::{
    AmbiguityClass, Code, Dictionary, Interpretation, Molecule, Options, Severity, Span, interpret,
    smiles,
};
use serde::Deserialize;

use crate::writer::condense;

/// The interpretation kind a case expects.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Expected {
    /// One reading.
    Unique,
    /// Two or more readings.
    Ambiguous,
    /// No structure encoded.
    Underdetermined,
    /// No reading.
    Invalid,
}

impl Expected {
    /// Name as used in the suite file.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Expected::Unique => "unique",
            Expected::Ambiguous => "ambiguous",
            Expected::Underdetermined => "underdetermined",
            Expected::Invalid => "invalid",
        }
    }
}

/// One conformance case.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Case {
    /// Unique identifier.
    pub id: String,
    /// `curated` or `adversarial`.
    pub corpus: String,
    /// The formula.
    pub input: String,
    /// Expected interpretation kind.
    pub expect: Expected,
    /// Expected readings in rank order, as SMILES.
    #[serde(default)]
    pub smiles: Vec<String>,
    /// Expected ambiguity classes.
    #[serde(default)]
    pub classes: Vec<String>,
    /// Expected underdetermination reason.
    #[serde(default)]
    pub reason: Option<String>,
    /// Expected diagnostic code.
    #[serde(default)]
    pub code: Option<String>,
    /// Expected primary span, as byte offsets in the input.
    #[serde(default)]
    pub span: Option<[usize; 2]>,
    /// Expected warning codes.
    #[serde(default)]
    pub warnings: Vec<String>,
    /// Free-form explanation.
    #[serde(default)]
    pub note: Option<String>,
}

/// Header of a suite file.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SuiteInfo {
    /// Suite name.
    pub name: String,
    /// Suite version.
    pub version: String,
    /// Grammar version the suite targets.
    pub grammar: u32,
    /// Ranking policy version the suite targets.
    pub ranking: u32,
}

/// A parsed conformance suite.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Suite {
    /// Header.
    pub suite: SuiteInfo,
    /// Cases.
    #[serde(rename = "case", default)]
    pub cases: Vec<Case>,
}

/// A suite file that could not be read.
#[derive(Debug, thiserror::Error)]
pub enum SuiteError {
    /// The file is not valid TOML or does not have the suite layout.
    #[error("conformance suite is malformed: {0}")]
    Toml(#[from] toml::de::Error),
    /// An expected SMILES string does not parse.
    #[error("case {id}: expected SMILES {smiles:?} does not parse: {error}")]
    Smiles {
        /// Case identifier.
        id: String,
        /// The SMILES.
        smiles: String,
        /// The parse error.
        error: smiles::SmilesError,
    },
    /// Two cases share an identifier.
    #[error("case id {0} is used twice")]
    DuplicateId(String),
}

/// Reads a suite from TOML text and checks every expected SMILES.
///
/// # Errors
///
/// Fails on malformed TOML, unparseable expected SMILES, or duplicate case ids.
pub fn load_suite(text: &str) -> Result<Suite, SuiteError> {
    let suite: Suite = toml::from_str(text)?;
    let mut ids = std::collections::BTreeSet::new();
    for case in &suite.cases {
        if !ids.insert(case.id.as_str()) {
            return Err(SuiteError::DuplicateId(case.id.clone()));
        }
        for expected in &case.smiles {
            smiles::parse(expected).map_err(|error| SuiteError::Smiles {
                id: case.id.clone(),
                smiles: expected.clone(),
                error,
            })?;
        }
    }
    Ok(suite)
}

/// A count of hits out of a total.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Ratio {
    /// Cases that met the condition.
    pub hits: usize,
    /// Cases the metric applies to.
    pub total: usize,
}

impl Ratio {
    fn record(&mut self, hit: bool) {
        self.total += 1;
        self.hits += usize::from(hit);
    }

    /// Hits as a fraction of the total, or `None` when the metric applies to no case.
    #[must_use]
    pub fn rate(self) -> Option<f64> {
        (self.total > 0).then(|| self.hits as f64 / self.total as f64)
    }
}

impl fmt::Display for Ratio {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.rate() {
            Some(rate) => write!(f, "{:.1}% ({}/{})", rate * 100.0, self.hits, self.total),
            None => f.write_str("n/a (0 cases)"),
        }
    }
}

/// The result of one case.
#[derive(Clone, Debug)]
pub struct Outcome {
    /// Case identifier.
    pub id: String,
    /// The case input.
    pub input: String,
    /// What the case expected.
    pub expected: Expected,
    /// The kind actually returned.
    pub actual: &'static str,
    /// Every way the result differs from the expectation. Empty when the case passes.
    pub problems: Vec<String>,
}

impl Outcome {
    /// Whether the result met every expectation.
    #[must_use]
    pub fn passed(&self) -> bool {
        self.problems.is_empty()
    }
}

/// Scores over a suite run. See the grammar, section 10.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Metrics {
    /// Cases that met every expectation.
    pub passed: Ratio,
    /// Expected-ambiguous cases returned as unique. Must be zero.
    pub false_unique: Ratio,
    /// Expected-ambiguous cases returned ambiguous with every expected class.
    pub class_recall: Ratio,
    /// Expected-ambiguous cases whose rank 1 reading is the expected first reading.
    pub top1: Ratio,
    /// Expected-invalid cases with the expected code and primary span.
    pub span_accuracy: Ratio,
    /// Expected-unique cases returned unique with the expected reading.
    pub unique_accuracy: Ratio,
    /// All cases returned ambiguous.
    pub ambiguity_rate: Ratio,
}

/// Every outcome of a suite run with its metrics.
#[derive(Clone, Debug)]
pub struct Report {
    /// Per-case outcomes, in suite order.
    pub outcomes: Vec<Outcome>,
    /// Aggregate scores.
    pub metrics: Metrics,
}

/// Runs every case of `suite`.
#[must_use]
pub fn run(suite: &Suite, dictionary: &Dictionary, options: &Options) -> Report {
    let mut metrics = Metrics::default();
    let outcomes: Vec<Outcome> = suite
        .cases
        .iter()
        .map(|case| {
            let result = interpret(&case.input, dictionary, options);
            let outcome = check(case, &result, &mut metrics);
            metrics.passed.record(outcome.passed());
            metrics
                .ambiguity_rate
                .record(matches!(result, Interpretation::Ambiguous { .. }));
            outcome
        })
        .collect();
    Report { outcomes, metrics }
}

fn canonical(smiles_text: &str) -> String {
    smiles::parse(smiles_text).map_or_else(|_| smiles_text.to_string(), |m| m.to_smiles())
}

fn check(case: &Case, result: &Interpretation, metrics: &mut Metrics) -> Outcome {
    let mut problems = Vec::new();
    if result.kind() != case.expect.as_str() {
        problems.push(format!(
            "expected {}, got {}",
            case.expect.as_str(),
            result.kind()
        ));
    }
    let expected: Vec<String> = case.smiles.iter().map(|s| canonical(s)).collect();
    let got: Vec<String> = result.readings().iter().map(|r| r.smiles.clone()).collect();
    match case.expect {
        Expected::Unique => {
            metrics
                .unique_accuracy
                .record(matches!(result, Interpretation::Unique(_)) && got == expected);
        }
        Expected::Ambiguous => {
            metrics
                .false_unique
                .record(matches!(result, Interpretation::Unique(_)));
            metrics.top1.record(got.first() == expected.first());
            let classes: Vec<&str> = match result {
                Interpretation::Ambiguous { classes, .. } => {
                    classes.iter().map(|c| c.as_str()).collect()
                }
                _ => Vec::new(),
            };
            let recalled = matches!(result, Interpretation::Ambiguous { .. })
                && case.classes.iter().all(|c| classes.contains(&c.as_str()));
            metrics.class_recall.record(recalled);
            if !case.classes.iter().all(|c| classes.contains(&c.as_str())) {
                problems.push(format!(
                    "expected classes {:?}, got {classes:?}",
                    case.classes
                ));
            }
            for class in &case.classes {
                if AmbiguityClass::parse(class).is_none() {
                    problems.push(format!("unknown class {class:?} in the suite"));
                }
            }
        }
        Expected::Underdetermined => {
            if let Interpretation::Underdetermined { span, reason } = result {
                if case.reason.as_deref().is_some_and(|r| r != reason.as_str()) {
                    problems.push(format!("expected reason {:?}, got {}", case.reason, reason));
                }
                if let Some(expected_span) = case.span {
                    compare_span(expected_span, *span, &mut problems);
                }
            }
        }
        Expected::Invalid => {
            let hit = match result {
                Interpretation::Invalid { diagnostics } => {
                    let first = &diagnostics[0];
                    let code_ok = case
                        .code
                        .as_deref()
                        .is_none_or(|c| Code::parse(c) == Some(first.code));
                    if !code_ok {
                        problems.push(format!("expected code {:?}, got {}", case.code, first.code));
                    }
                    let before = problems.len();
                    if let Some(expected_span) = case.span {
                        compare_span(expected_span, first.primary.span, &mut problems);
                    }
                    code_ok && problems.len() == before
                }
                _ => false,
            };
            metrics.span_accuracy.record(hit);
        }
    }
    if matches!(case.expect, Expected::Unique | Expected::Ambiguous)
        && !expected.is_empty()
        && got != expected
    {
        problems.push(format!("expected readings {expected:?}, got {got:?}"));
    }
    let warnings: Vec<&str> = result
        .readings()
        .first()
        .map(|r| {
            r.warnings
                .iter()
                .filter(|w| w.severity() == Severity::Warning)
                .map(|w| w.code.as_str())
                .collect()
        })
        .unwrap_or_default();
    for code in &case.warnings {
        if !warnings.contains(&code.as_str()) {
            problems.push(format!("expected warning {code}, got {warnings:?}"));
        }
    }
    Outcome {
        id: case.id.clone(),
        input: case.input.clone(),
        expected: case.expect,
        actual: result.kind(),
        problems,
    }
}

fn compare_span(expected: [usize; 2], actual: Span, problems: &mut Vec<String>) {
    if Span::new(expected[0], expected[1]) != actual {
        problems.push(format!(
            "expected span {}..{}, got {actual}",
            expected[0], expected[1]
        ));
    }
}

/// Acceptance thresholds for v1.0, from the grammar, section 10.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Thresholds {
    /// Highest allowed false-unique rate.
    pub false_unique: f64,
    /// Lowest allowed round-trip recovery.
    pub round_trip: f64,
    /// Lowest allowed ambiguity class recall.
    pub class_recall: f64,
    /// Lowest allowed top-1 accuracy on ambiguous cases.
    pub top1: f64,
    /// Lowest allowed diagnostic span accuracy.
    pub span_accuracy: f64,
}

impl Default for Thresholds {
    fn default() -> Self {
        Self {
            false_unique: 0.0,
            round_trip: 0.95,
            class_recall: 0.98,
            top1: 0.85,
            span_accuracy: 0.90,
        }
    }
}

/// Lists every threshold the metrics miss. Empty when the gate passes.
#[must_use]
pub fn gate(
    metrics: &Metrics,
    round_trip: Option<&RoundTrip>,
    thresholds: &Thresholds,
) -> Vec<String> {
    let mut misses = Vec::new();
    let mut at_most = |name: &str, ratio: Ratio, limit: f64| {
        if ratio.rate().is_some_and(|rate| rate > limit) {
            misses.push(format!("{name} {ratio} is above {:.1}%", limit * 100.0));
        }
    };
    at_most(
        "false-unique rate",
        metrics.false_unique,
        thresholds.false_unique,
    );
    let mut at_least = |name: &str, ratio: Ratio, limit: f64| {
        if ratio.rate().is_some_and(|rate| rate < limit) {
            misses.push(format!("{name} {ratio} is below {:.1}%", limit * 100.0));
        }
    };
    at_least(
        "ambiguity class recall",
        metrics.class_recall,
        thresholds.class_recall,
    );
    at_least("top-1 accuracy", metrics.top1, thresholds.top1);
    at_least(
        "diagnostic span accuracy",
        metrics.span_accuracy,
        thresholds.span_accuracy,
    );
    if let Some(round_trip) = round_trip {
        at_least(
            "round-trip recovery",
            round_trip.recovered,
            thresholds.round_trip,
        );
    }
    misses
}

/// One molecule that did not survive the round trip.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RoundTripMiss {
    /// Canonical SMILES of the original molecule.
    pub smiles: String,
    /// What the inverse writer produced.
    pub written: String,
    /// What reading it back produced: the rank 1 SMILES, or the result kind.
    pub read: String,
}

/// Results of writing molecules as condensed formulas and reading them back.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RoundTrip {
    /// Molecules read back to the same canonical SMILES at rank 1.
    pub recovered: Ratio,
    /// Molecules the writer could not express, excluded from `recovered`.
    pub unwritable: usize,
    /// Molecules read back as ambiguous, whatever their rank 1 reading.
    pub ambiguous: usize,
    /// Every molecule that was not recovered.
    pub misses: Vec<RoundTripMiss>,
}

/// Writes each molecule with the inverse writer, reads it back, and compares canonical SMILES.
#[must_use]
pub fn round_trip(molecules: &[Molecule], dictionary: &Dictionary, options: &Options) -> RoundTrip {
    let mut result = RoundTrip::default();
    for molecule in molecules {
        let Ok(written) = condense(molecule, dictionary) else {
            result.unwritable += 1;
            continue;
        };
        let original = molecule.to_smiles();
        let interpretation = interpret(&written, dictionary, options);
        result.ambiguous += usize::from(matches!(interpretation, Interpretation::Ambiguous { .. }));
        let read = interpretation
            .top()
            .map_or_else(|| interpretation.kind().to_string(), |r| r.smiles.clone());
        let hit = read == original;
        result.recovered.record(hit);
        if !hit {
            result.misses.push(RoundTripMiss {
                smiles: original,
                written,
                read,
            });
        }
    }
    result
}

/// Reads a SMILES list: one molecule per line, SMILES first, `#` starts a comment.
///
/// # Errors
///
/// Returns the line number and error of the first line that does not parse.
pub fn read_smiles_list(text: &str) -> Result<Vec<Molecule>, (usize, smiles::SmilesError)> {
    text.lines()
        .enumerate()
        .filter_map(|(i, line)| {
            let content = line.split('#').next().unwrap_or_default().trim();
            let first = content.split_whitespace().next()?;
            Some(smiles::parse(first).map_err(|e| (i + 1, e)))
        })
        .collect()
}

/// A deterministic corpus of `count` random acyclic organic molecules.
///
/// Carbon forms the backbone; oxygen, nitrogen, sulfur and halogens attach to carbon,
/// and carbon may also attach to oxygen, nitrogen and sulfur, giving ethers, amines
/// and thioethers. Some carbon-carbon and carbon-heteroatom bonds are multiple.
#[must_use]
pub fn generated_corpus(count: usize, seed: u64) -> Vec<Molecule> {
    let mut rng = XorShift(seed.max(1));
    let mut seen = std::collections::BTreeSet::new();
    let mut corpus = Vec::new();
    let mut attempts = 0;
    while corpus.len() < count && attempts < count * 50 {
        attempts += 1;
        let molecule = random_molecule(&mut rng);
        if seen.insert(molecule.to_smiles()) {
            corpus.push(molecule);
        }
    }
    corpus
}

struct XorShift(u64);

impl XorShift {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

fn random_molecule(rng: &mut XorShift) -> Molecule {
    use moiety_core::{Atom, Element};
    const CHOICES: [(Element, u64); 7] = [
        (Element::C, 62),
        (Element::O, 14),
        (Element::N, 10),
        (Element::CL, 4),
        (Element::F, 4),
        (Element::S, 3),
        (Element::BR, 3),
    ];
    let size = 2 + rng.below(9);
    let mut molecule = Molecule::new();
    let mut free: Vec<u8> = Vec::new();
    molecule.add_atom(Atom::new(Element::C));
    free.push(4);
    for _ in 1..size {
        let mut roll = rng.below(100);
        let element = CHOICES
            .iter()
            .find(|(_, weight)| {
                let hit = roll < *weight;
                roll = roll.saturating_sub(*weight);
                hit
            })
            .map_or(Element::C, |(e, _)| *e);
        let valence = element.valences()[0];
        let hosts: Vec<usize> = (0..free.len())
            .filter(|&i| free[i] > 0)
            .filter(|&i| element == Element::C || molecule.atoms()[i].element() == Some(Element::C))
            .collect();
        if hosts.is_empty() {
            break;
        }
        let host = hosts[usize::try_from(rng.below(hosts.len() as u64)).unwrap_or(0)];
        let carbon_host = molecule.atoms()[host].element() == Some(Element::C);
        let max_order = free[host].min(valence).min(if carbon_host { 3 } else { 1 });
        let multiple =
            matches!(element, Element::C | Element::O | Element::N) && rng.below(100) < 18;
        let order = if multiple && max_order >= 2 {
            2 + u8::from(max_order == 3 && rng.below(3) == 0)
        } else {
            1
        };
        let order = order.min(max_order).max(1);
        let atom = molecule.add_atom(Atom::new(element));
        free.push(valence - order);
        free[host] -= order;
        molecule.add_bond(
            host,
            atom,
            moiety_core::BondOrder::from_count(order).expect("orders are 1 to 3"),
        );
    }
    for (i, &left) in free.iter().enumerate() {
        molecule.atoms_mut()[i].hydrogens = left;
    }
    molecule
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::default_dictionary;

    const SUITE: &str = r#"
        [suite]
        name = "t"
        version = "1"
        grammar = 1
        ranking = 1

        [[case]]
        id = "a"
        corpus = "curated"
        input = "CH3CH2OH"
        expect = "unique"
        smiles = ["OCC"]

        [[case]]
        id = "b"
        corpus = "adversarial"
        input = "COCl2"
        expect = "ambiguous"
        smiles = ["ClC(Cl)=O", "Cl[Co]Cl"]
        classes = ["symbol-segmentation"]

        [[case]]
        id = "c"
        corpus = "adversarial"
        input = "CH3CH3OH"
        expect = "invalid"
        code = "CF0204"
        span = [3, 6]
    "#;

    #[test]
    fn scores_a_small_suite() {
        let suite = load_suite(SUITE).unwrap();
        let report = run(&suite, default_dictionary(), &Options::default());
        assert!(
            report.outcomes.iter().all(Outcome::passed),
            "{:?}",
            report.outcomes
        );
        assert_eq!(report.metrics.false_unique, Ratio { hits: 0, total: 1 });
        assert_eq!(report.metrics.span_accuracy, Ratio { hits: 1, total: 1 });
        assert!(gate(&report.metrics, None, &Thresholds::default()).is_empty());
    }

    #[test]
    fn rejects_bad_suites() {
        let duplicate = SUITE.replace("id = \"b\"", "id = \"a\"");
        assert!(matches!(
            load_suite(&duplicate),
            Err(SuiteError::DuplicateId(_))
        ));
        let bad_smiles = SUITE.replace("\"OCC\"", "\"C(\"");
        assert!(matches!(
            load_suite(&bad_smiles),
            Err(SuiteError::Smiles { .. })
        ));
    }

    #[test]
    fn generated_corpus_is_deterministic_and_valid() {
        let first = generated_corpus(50, 7);
        let second = generated_corpus(50, 7);
        assert_eq!(first.len(), 50);
        assert_eq!(first, second);
        for molecule in &first {
            assert!(smiles::parse(&molecule.to_smiles()).is_ok());
        }
    }

    #[test]
    fn reads_smiles_lists() {
        let molecules = read_smiles_list("CCO ethanol\n# comment\n\nc1ccccc1 # benzene\n").unwrap();
        assert_eq!(molecules.len(), 2);
        assert_eq!(read_smiles_list("CC\nC(\n").unwrap_err().0, 2);
    }
}
