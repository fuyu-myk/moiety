//! Checks the recorded differential results against MolScribe: every disagreement
//! has a verdict, none is decided against moiety, and the committed report is current.

use moiety::differential::{Verdict, compare, load_adjudications, load_external, report, summary};
use moiety::{Options, default_dictionary};

fn differential_file(name: &str) -> String {
    let path = format!(
        "{}/../../conformance/differential/{name}",
        env!("CARGO_MANIFEST_DIR")
    );
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"))
}

#[test]
fn every_disagreement_with_molscribe_has_a_verdict() {
    let external = load_external(&differential_file("molscribe.toml")).expect("results parse");
    let verdicts =
        load_adjudications(&differential_file("adjudications.toml")).expect("verdicts parse");
    let comparisons = compare(
        &external,
        &verdicts,
        default_dictionary(),
        &Options::default(),
    );
    let open: Vec<String> = comparisons
        .iter()
        .filter(|c| c.unadjudicated())
        .map(|c| {
            format!(
                "{} {} moiety {:?} other {:?}",
                c.input, c.agreement, c.readings, c.other
            )
        })
        .collect();
    assert!(
        open.is_empty(),
        "disagreements without a verdict:\n{}",
        open.join("\n")
    );
    let against: Vec<&str> = comparisons
        .iter()
        .filter(|c| c.verdict == Some(Verdict::Other))
        .map(|c| c.input.as_str())
        .collect();
    assert!(
        against.is_empty(),
        "verdicts decided against moiety: {against:?}"
    );
    assert_eq!(summary(&comparisons).1, 0);
    assert_eq!(
        report(&external, &comparisons),
        differential_file("report.toml"),
        "report.toml is stale; regenerate it with `moiety differential conformance/differential/molscribe.toml -o conformance/differential/report.toml`"
    );
}
