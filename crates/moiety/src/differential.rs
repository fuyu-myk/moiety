//! Differential testing against another condensed-formula implementation.
//!
//! The other implementation's results for a set of substituent labels are read
//! from TOML (see `conformance/tools/molscribe_labels.py`). moiety reads each
//! label with an attachment marker on the side the label bonds from, and the
//! two results are compared by canonical SMILES. Disagreements are cases to
//! adjudicate, not failures: a recorded verdict says which side is right.

use std::collections::BTreeMap;
use std::fmt;

use moiety_core::{Dictionary, Interpretation, Options, interpret, smiles};
use serde::{Deserialize, Serialize};

/// Which end of a label bonds to the parent structure.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Deserialize, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Side {
    /// The label is written from its attachment, as in `CH2OH`.
    Left,
    /// The label is written toward its attachment, as in `HOOC`.
    Right,
}

/// The other implementation's results, as written by the harness.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct External {
    /// The implementation and version that produced the results.
    pub source: String,
    /// The toolkit used to canonicalize them.
    pub toolkit: String,
    /// One result per label.
    #[serde(rename = "label", default)]
    pub labels: Vec<ExternalLabel>,
}

/// One label as the other implementation read it.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExternalLabel {
    /// The label.
    pub text: String,
    /// Which end bonds to the parent.
    pub side: Side,
    /// How the label was expanded, such as `abbreviation` or `formula`.
    #[serde(default)]
    pub route: Option<String>,
    /// SMILES with `*` on the attachment atom, when the label was expanded.
    #[serde(default)]
    pub smiles: Option<String>,
    /// Why the label was not expanded.
    #[serde(default)]
    pub error: Option<String>,
}

/// A recorded judgement on one disagreement.
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Adjudication {
    /// The label.
    pub text: String,
    /// Which end bonds to the parent.
    pub side: Side,
    /// Which side is right.
    pub verdict: Verdict,
    /// Why.
    pub note: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct AdjudicationFile {
    #[serde(rename = "adjudication", default)]
    adjudications: Vec<Adjudication>,
}

/// Which implementation a disagreement is decided for.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Verdict {
    /// moiety's result is right.
    Moiety,
    /// The other implementation's result is right, so moiety has a bug.
    Other,
    /// Both results are defensible readings of the label.
    Both,
}

/// How moiety's result relates to the other implementation's.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Agreement {
    /// The other result is moiety's rank 1 reading.
    Same,
    /// The other result is one of moiety's lower-ranked readings.
    Alternative,
    /// Both read the label, to different molecules.
    Different,
    /// moiety reports no reading; the other implementation gives one.
    MoietyRefuses,
    /// moiety reads the label; the other implementation does not.
    OtherFails,
    /// Neither reads the label.
    BothFail,
}

impl Agreement {
    /// Whether the outcome needs a recorded verdict.
    #[must_use]
    pub const fn needs_verdict(self) -> bool {
        matches!(
            self,
            Agreement::Alternative | Agreement::Different | Agreement::MoietyRefuses
        )
    }

    /// Name as used in reports.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Agreement::Same => "same",
            Agreement::Alternative => "alternative",
            Agreement::Different => "different",
            Agreement::MoietyRefuses => "moiety-refuses",
            Agreement::OtherFails => "other-fails",
            Agreement::BothFail => "both-fail",
        }
    }
}

impl fmt::Display for Agreement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// One label compared.
#[derive(Clone, Debug, Serialize)]
pub struct Comparison {
    /// The label.
    pub text: String,
    /// Which end bonds to the parent.
    pub side: Side,
    /// What moiety was given: the label with an attachment marker.
    pub input: String,
    /// moiety's result kind.
    pub moiety: &'static str,
    /// moiety's readings in rank order, as canonical SMILES.
    pub readings: Vec<String>,
    /// The other implementation's result, as canonical SMILES.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub other: Option<String>,
    /// Why the other implementation gave no result.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub other_error: Option<String>,
    /// How the results relate.
    pub agreement: Agreement,
    /// The recorded verdict, if any.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub verdict: Option<Verdict>,
    /// The recorded reason, if any.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub note: Option<String>,
}

impl Comparison {
    /// Whether the label needs a verdict and has none.
    #[must_use]
    pub fn unadjudicated(&self) -> bool {
        self.agreement.needs_verdict() && self.verdict.is_none()
    }
}

/// A file that could not be read.
#[derive(Debug, thiserror::Error)]
pub enum DifferentialError {
    /// Malformed TOML.
    #[error("malformed differential file: {0}")]
    Toml(#[from] toml::de::Error),
}

/// Reads the other implementation's results.
///
/// # Errors
///
/// Fails on malformed TOML.
pub fn load_external(text: &str) -> Result<External, DifferentialError> {
    Ok(toml::from_str(text)?)
}

/// Reads recorded verdicts.
///
/// # Errors
///
/// Fails on malformed TOML.
pub fn load_adjudications(text: &str) -> Result<Vec<Adjudication>, DifferentialError> {
    Ok(toml::from_str::<AdjudicationFile>(text)?.adjudications)
}

/// The moiety input for a label: the label with an attachment marker on its bonding side.
#[must_use]
pub fn moiety_input(text: &str, side: Side) -> String {
    match side {
        Side::Left => format!("-{text}"),
        Side::Right => format!("{text}-"),
    }
}

/// Compares moiety with the other implementation on every label.
#[must_use]
pub fn compare(
    external: &External,
    adjudications: &[Adjudication],
    dictionary: &Dictionary,
    options: &Options,
) -> Vec<Comparison> {
    let verdicts: BTreeMap<(&str, Side), &Adjudication> = adjudications
        .iter()
        .map(|a| ((a.text.as_str(), a.side), a))
        .collect();
    external
        .labels
        .iter()
        .map(|label| {
            let input = moiety_input(&label.text, label.side);
            let result = interpret(&input, dictionary, options);
            let readings: Vec<String> =
                result.readings().iter().map(|r| r.smiles.clone()).collect();
            let other = label
                .smiles
                .as_deref()
                .map(|s| smiles::parse(s).map_or_else(|_| s.to_string(), |m| m.to_smiles()));
            let agreement = agreement(&result, &readings, other.as_deref());
            let recorded = verdicts.get(&(label.text.as_str(), label.side));
            Comparison {
                text: label.text.clone(),
                side: label.side,
                input,
                moiety: result.kind(),
                readings,
                other,
                other_error: label.error.clone(),
                agreement,
                verdict: recorded.map(|a| a.verdict),
                note: recorded.map(|a| a.note.clone()),
            }
        })
        .collect()
}

fn agreement(result: &Interpretation, readings: &[String], other: Option<&str>) -> Agreement {
    match (result.readings().is_empty(), other) {
        (true, None) => Agreement::BothFail,
        (true, Some(_)) => Agreement::MoietyRefuses,
        (false, None) => Agreement::OtherFails,
        (false, Some(other)) if readings[0] == other => Agreement::Same,
        (false, Some(other)) if readings.iter().any(|r| r == other) => Agreement::Alternative,
        (false, Some(_)) => Agreement::Different,
    }
}

/// Count of comparisons per agreement, and of those still needing a verdict.
#[must_use]
pub fn summary(comparisons: &[Comparison]) -> (BTreeMap<Agreement, usize>, usize) {
    let mut counts = BTreeMap::new();
    for comparison in comparisons {
        *counts.entry(comparison.agreement).or_insert(0) += 1;
    }
    let open = comparisons.iter().filter(|c| c.unadjudicated()).count();
    (counts, open)
}

#[derive(Serialize)]
struct Report<'a> {
    source: &'a str,
    toolkit: &'a str,
    counts: BTreeMap<String, usize>,
    unadjudicated: usize,
    #[serde(rename = "label")]
    labels: Vec<&'a Comparison>,
}

/// The comparisons as a TOML report, listing every label that is not `same`.
///
/// # Panics
///
/// Never in practice: the report contains only strings, counts and lists.
#[must_use]
pub fn report(external: &External, comparisons: &[Comparison]) -> String {
    let (counts, unadjudicated) = summary(comparisons);
    let report = Report {
        source: &external.source,
        toolkit: &external.toolkit,
        counts: counts
            .into_iter()
            .map(|(agreement, n)| (agreement.as_str().to_string(), n))
            .collect(),
        unadjudicated,
        labels: comparisons
            .iter()
            .filter(|c| c.agreement != Agreement::Same)
            .collect(),
    };
    toml::to_string(&report).expect("the report serializes to TOML")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::default_dictionary;

    const EXTERNAL: &str = r#"
        source = "test"
        toolkit = "none"

        [[label]]
        text = "CH2OH"
        side = "left"
        smiles = "*CO"

        [[label]]
        text = "HOOC"
        side = "right"
        smiles = "*C(=O)O"

        [[label]]
        text = "C3H7"
        side = "left"
        smiles = "*CCC"

        [[label]]
        text = "CH2"
        side = "left"
        error = "no expansion"
    "#;

    #[test]
    fn classifies_agreement() {
        let external = load_external(EXTERNAL).unwrap();
        let verdicts = load_adjudications(
            "[[adjudication]]\ntext = \"C3H7\"\nside = \"left\"\nverdict = \"moiety\"\nnote = \"isomers\"\n",
        )
        .unwrap();
        let comparisons = compare(
            &external,
            &verdicts,
            default_dictionary(),
            &Options::default(),
        );
        let agreements: Vec<Agreement> = comparisons.iter().map(|c| c.agreement).collect();
        assert_eq!(
            agreements,
            [
                Agreement::Same,
                Agreement::Same,
                Agreement::MoietyRefuses,
                Agreement::BothFail
            ]
        );
        assert_eq!(comparisons[1].input, "HOOC-");
        assert_eq!(comparisons[2].verdict, Some(Verdict::Moiety));
        assert_eq!(summary(&comparisons).1, 0);
        let report = report(&external, &comparisons);
        assert!(report.contains("moiety-refuses"));
        assert!(!report.contains("text = \"CH2OH\""));
    }
}
