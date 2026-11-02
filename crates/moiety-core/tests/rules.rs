//! The assembly rules of the grammar, checked through the public API with a small dictionary.

use moiety_core::{
    AmbiguityClass, ChargeConvention, Code, Dictionary, Direction, Entry, Interpretation, Options,
    Span, Underdetermination, interpret,
};

fn dictionary() -> Dictionary {
    let mut d = Dictionary::new("rules", "1");
    for (name, smiles, spellings, reversed) in [
        ("phenyl", "*c1ccccc1", &["Ph"][..], &[][..]),
        ("carboxyl", "*C(=O)O", &["COOH"][..], &["HOOC"][..]),
        ("ethylene", "*CC*", &["C2H4"][..], &[][..]),
        ("water", "O", &["Aq"][..], &[][..]),
    ] {
        let entry = Entry::new(
            name,
            smiles,
            spellings.iter().map(ToString::to_string).collect(),
            reversed.iter().map(ToString::to_string).collect(),
        )
        .unwrap();
        d.insert(entry).unwrap();
    }
    d
}

fn read(input: &str) -> Interpretation {
    interpret(input, &dictionary(), &Options::default())
}

fn canon(smiles: &str) -> String {
    moiety_core::smiles::parse(smiles).unwrap().to_smiles()
}

fn unique(input: &str) -> String {
    match read(input) {
        Interpretation::Unique(reading) => reading.smiles,
        other => panic!("{input}: expected one reading, got {other:?}"),
    }
}

fn reads_as(input: &str, smiles: &str) {
    assert_eq!(unique(input), canon(smiles), "{input}");
}

fn error_code(input: &str) -> Code {
    match read(input) {
        Interpretation::Invalid { diagnostics } => diagnostics[0].code,
        other => panic!("{input}: expected invalid, got {}", other.kind()),
    }
}

#[test]
fn r1_separates_pendants_from_chain_segments() {
    reads_as("(CH3)2CHOH", "CC(C)O");
    reads_as("CH3(CH2)3CH3", "CCCCC");
    reads_as("CH2(CH)2CH2", "C=CC=C");
    assert_eq!(error_code("CH3(CCl3)2"), Code::HostFull);
}

#[test]
fn r2_anchors_leading_pendants_on_the_first_atom() {
    reads_as("(CH3)3COH", "CC(C)(C)O");
    reads_as("H2NCH3", "CN");
    reads_as("Cl2", "ClCl");
    assert_eq!(error_code("(CH3)3OH"), Code::NoAnchor);
}

#[test]
fn r3_solves_bond_orders_with_topology() {
    reads_as("CH3COOH", "CC(=O)O");
    reads_as("HCCH", "C#C");
    reads_as("CH2CCH2", "C=C=C");
    assert_eq!(error_code("CH3=CH3"), Code::BondNotPossible);
}

#[test]
fn r4_rejects_open_valence_and_unplaced_charge() {
    assert_eq!(error_code("CH3CH2"), Code::OpenValence);
    assert_eq!(error_code("CH4+"), Code::ChargeNotPlaced);
    reads_as("-CH2CH3", "*CC");
}

#[test]
fn r5_reads_salts_with_metals_first() {
    reads_as("K2SO4", "[K+].[K+].[O-]S(=O)(=O)[O-]");
    reads_as("NaCl", "[Na]Cl");
}

#[test]
fn hydrogens_bond_to_the_unit_before_them() {
    reads_as("H2NOH", "NO");
    reads_as("CH3CO2H", "CC(=O)O");
}

#[test]
fn host_search_stops_at_hydrogenated_atoms() {
    reads_as("NCCH2CN", "N#CCC#N");
    assert_eq!(error_code("CH3CH3OH"), Code::NoHost);
}

#[test]
fn transition_metals_avoid_carbon_and_metals_avoid_each_other() {
    reads_as("CH3COOH", "CC(=O)O");
    assert_eq!(error_code("CH3CoCH3"), Code::BondNotPossible);
    reads_as("CS2", "S=C=S");
    match read("COCl2") {
        Interpretation::Ambiguous { readings, classes } => {
            assert_eq!(readings[0].smiles, canon("O=C(Cl)Cl"));
            assert_eq!(classes, vec![AmbiguityClass::SymbolSegmentation]);
        }
        other => panic!("expected ambiguous, got {}", other.kind()),
    }
}

#[test]
fn charge_convention_is_switchable() {
    let separated = Options::default();
    let hypervalent = Options {
        charge_convention: ChargeConvention::Hypervalent,
        ..Options::default()
    };
    let top = |options: &Options| {
        interpret("CH3NO2", &dictionary(), options)
            .top()
            .unwrap()
            .smiles
            .clone()
    };
    assert_eq!(top(&separated), canon("C[N+](=O)[O-]"));
    assert_eq!(top(&hypervalent), canon("CN(=O)=O"));
}

#[test]
fn dictionary_fragments_attach_by_their_sites() {
    reads_as("PhCOOH", "OC(=O)c1ccccc1");
    reads_as("HOOCPh", "OC(=O)c1ccccc1");
    reads_as("ClC2H4Cl", "ClCCCl");
    reads_as("Aq", "O");
    assert_eq!(error_code("AqCH3"), Code::MoleculeAsGroup);
}

#[test]
fn direction_is_recorded() {
    let direction = |input: &str| read(input).top().unwrap().direction;
    assert_eq!(direction("CH3CH2OH"), Direction::LeftToRight);
    assert_eq!(direction("(CH3)3N"), Direction::RightToLeft);
    assert_eq!(direction("(CH3)2CHOH"), Direction::Mixed);
    assert_eq!(direction("HOOCCH3"), Direction::RightToLeft);
}

#[test]
fn provenance_maps_atoms_to_original_spans() {
    let result = read("CH₃CH₂OH");
    let reading = result.top().unwrap();
    let sources: Vec<Span> = (0..3)
        .map(|i| reading.provenance.atom_source(i).unwrap())
        .collect();
    assert_eq!(
        sources,
        vec![Span::new(0, 5), Span::new(5, 10), Span::new(10, 12)]
    );
    assert_eq!(reading.provenance.atoms_in(Span::new(6, 7)), vec![1]);
}

#[test]
fn underdetermined_inputs_name_their_reason() {
    let reason = |input: &str| match read(input) {
        Interpretation::Underdetermined { reason, .. } => reason,
        other => panic!("{input}: expected underdetermined, got {}", other.kind()),
    };
    assert_eq!(reason("C6H6"), Underdetermination::MolecularFormula);
    assert_eq!(reason("CH3C6H4CH3"), Underdetermination::UnknownRing);
    assert_eq!(reason("CH3Qq"), Underdetermination::UnknownAbbreviation);
}

#[test]
fn search_limit_reports_too_complex() {
    let tiny = Options {
        search_limit: 3,
        suggestions: false,
        ..Options::default()
    };
    match interpret("CH3CH2CH2CH2OH", &dictionary(), &tiny) {
        Interpretation::Invalid { diagnostics } => {
            assert_eq!(diagnostics[0].code, Code::TooComplex)
        }
        other => panic!("expected invalid, got {}", other.kind()),
    }
}

#[test]
fn reading_limit_never_hides_an_ambiguity() {
    let one = Options {
        max_readings: 1,
        ..Options::default()
    };
    match interpret("COCl2", &dictionary(), &one) {
        Interpretation::Ambiguous { readings, .. } => assert_eq!(readings.len(), 2),
        other => panic!("expected ambiguous, got {}", other.kind()),
    }
}

#[test]
fn caret_charges_do_not_take_a_count() {
    assert_ne!(
        read("Fe2^+").top().map(|r| r.smiles.clone()),
        Some(canon("[Fe+2]"))
    );
    assert_eq!(error_code("C^12-"), Code::CountOutOfRange);
}

#[test]
fn repeated_groups_spend_the_search_budget() {
    let nested = format!("{}CH2{}", "(".repeat(18), ")2".repeat(18));
    assert_eq!(error_code(&nested), Code::TooComplex);
    reads_as("CH3(CH2)100CH3", &format!("C{}C", "C".repeat(100)));
}

#[test]
fn pathological_inputs_stop_at_the_search_limit() {
    let quick = Options {
        suggestions: false,
        ..Options::default()
    };
    for input in [
        "CH3SSSSSSSS100",
        "MnMnMnMnMnMnMnMnMnMnMnMnMnMnMnMnMnMnMnMnO",
    ] {
        let result = interpret(input, &dictionary(), &quick);
        assert!(
            matches!(&result, Interpretation::Invalid { diagnostics } if diagnostics[0].code == Code::TooComplex)
                || !matches!(result, Interpretation::Invalid { .. }),
            "{input}: {}",
            result.kind()
        );
    }
}
