//! Property tests: any input yields a well-formed interpretation without panicking.

use moiety::{Interpretation, Options, default_dictionary, interpret_with, smiles};
use proptest::prelude::*;

fn formula() -> impl Strategy<Value = String> {
    let pieces = prop::sample::select(vec![
        "C", "H", "O", "N", "S", "Cl", "Br", "F", "P", "Co", "Na", "Ts", "Ph", "Me", "Et", "Ac",
        "COOH", "HOOC", "CO2Et", "NO2", "2", "3", "4", "(", ")", "[", "]", "-", "=", "#", "+", ".",
        "^2-", "Xy", "D", "$", "₂", "·",
    ]);
    prop::collection::vec(pieces, 0..12).prop_map(|parts| parts.concat())
}

fn check(input: &str, result: &Interpretation) -> Result<(), TestCaseError> {
    let within = |span: moiety::Span| span.start <= span.end && span.end <= input.len();
    match result {
        Interpretation::Unique(_) | Interpretation::Ambiguous { .. } => {
            for (i, reading) in result.readings().iter().enumerate() {
                prop_assert_eq!(reading.rank as usize, i + 1);
                let reparsed = smiles::parse(&reading.smiles);
                prop_assert!(
                    reparsed.is_ok(),
                    "{} gave unreadable SMILES {}",
                    input,
                    reading.smiles
                );
                prop_assert_eq!(reparsed.unwrap().to_smiles(), reading.smiles.clone());
                prop_assert_eq!(
                    reading.provenance.atoms().len(),
                    reading.molecule.atoms().len()
                );
                prop_assert_eq!(
                    reading.provenance.bonds().len(),
                    reading.molecule.bonds().len()
                );
                for atom in 0..reading.molecule.atoms().len() {
                    prop_assert!(within(reading.provenance.atom_source(atom).unwrap()));
                }
            }
            if let Interpretation::Ambiguous { readings, classes } = result {
                prop_assert!(readings.len() >= 2);
                let mut distinct: Vec<&str> = readings.iter().map(|r| r.smiles.as_str()).collect();
                distinct.sort_unstable();
                distinct.dedup();
                prop_assert_eq!(
                    distinct.len(),
                    readings.len(),
                    "{} repeats a reading",
                    input
                );
                prop_assert!(!classes.is_empty(), "{} is ambiguous with no class", input);
            }
        }
        Interpretation::Underdetermined { span, .. } => prop_assert!(within(*span)),
        Interpretation::Invalid { diagnostics } => {
            prop_assert!(!diagnostics.is_empty());
            for diagnostic in diagnostics {
                prop_assert!(
                    within(diagnostic.primary.span),
                    "{} {:?}",
                    input,
                    diagnostic
                );
            }
        }
    }
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(400))]

    #[test]
    fn random_formulas_never_panic(input in formula()) {
        let options = Options { suggestions: false, ..Options::default() };
        let result = interpret_with(&input, default_dictionary(), &options);
        check(&input, &result)?;
    }

    #[test]
    fn suggestions_never_panic(input in formula().prop_filter("short", |s| s.len() <= 14)) {
        let result = interpret_with(&input, default_dictionary(), &Options::default());
        check(&input, &result)?;
    }

    #[test]
    fn interpretation_is_deterministic(input in formula()) {
        let options = Options { suggestions: false, ..Options::default() };
        let first = interpret_with(&input, default_dictionary(), &options);
        let second = interpret_with(&input, default_dictionary(), &options);
        prop_assert_eq!(first, second);
    }
}
