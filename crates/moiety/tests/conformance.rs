//! Runs the published conformance suite and round-trip corpora against the v1.0 gate.

use moiety::conformance::{
    Thresholds, gate, generated_corpus, load_suite, read_smiles_list, round_trip, run,
};
use moiety::{Options, default_dictionary};

fn repository_file(name: &str) -> String {
    let path = format!("{}/../../conformance/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"))
}

#[test]
fn suite_passes_every_case_and_the_gate() {
    let suite = load_suite(&repository_file("cases.toml")).expect("the suite parses");
    let report = run(&suite, default_dictionary(), &Options::default());
    let failures: Vec<String> = report
        .outcomes
        .iter()
        .filter(|o| !o.passed())
        .map(|o| format!("{} ({}): {}", o.id, o.input, o.problems.join("; ")))
        .collect();
    assert!(
        failures.is_empty(),
        "failing cases:\n{}",
        failures.join("\n")
    );
    let misses = gate(&report.metrics, None, &Thresholds::default());
    assert!(misses.is_empty(), "{misses:?}");
}

#[test]
fn curated_corpus_round_trips() {
    let molecules = read_smiles_list(&repository_file("roundtrip.smi")).expect("the corpus parses");
    let result = round_trip(&molecules, default_dictionary(), &Options::default());
    let misses: Vec<String> = result
        .misses
        .iter()
        .map(|m| format!("{} written {} read {}", m.smiles, m.written, m.read))
        .collect();
    assert_eq!(result.unwritable, 0);
    assert!(
        misses.is_empty(),
        "round-trip misses:\n{}",
        misses.join("\n")
    );
}

#[test]
fn generated_corpus_meets_the_round_trip_threshold() {
    let molecules = generated_corpus(500, 2026);
    let result = round_trip(&molecules, default_dictionary(), &Options::default());
    let misses = gate(&Default::default(), Some(&result), &Thresholds::default());
    let sample: Vec<String> = result
        .misses
        .iter()
        .take(20)
        .map(|m| format!("{} written {} read {}", m.smiles, m.written, m.read))
        .collect();
    assert!(misses.is_empty(), "{misses:?}\n{}", sample.join("\n"));
}
