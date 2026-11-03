use super::*;
use crate::{default_dictionary, interpret};

fn write(smiles: &str) -> String {
    let molecule = moiety_core::smiles::parse(smiles).unwrap();
    condense(&molecule, default_dictionary()).unwrap()
}

fn round_trips(smiles: &str) {
    let molecule = moiety_core::smiles::parse(smiles).unwrap();
    let written = condense(&molecule, default_dictionary()).unwrap();
    let top = interpret(&written).top().map(|r| r.smiles.clone());
    assert_eq!(
        top,
        Some(molecule.to_smiles()),
        "{smiles} written as {written}"
    );
}

#[test]
fn writes_natural_condensed_formulas() {
    assert_eq!(write("CCO"), "CH3CH2OH");
    assert_eq!(write("CC(=O)O"), "CH3COOH");
    assert_eq!(write("CC(C)C"), "CH3CH(CH3)CH3");
    assert_eq!(write("OCc1ccccc1"), "BnOH");
    assert_eq!(write("C=CC"), "CH2=CHCH3");
    assert_eq!(write("O.O.O"), "3H2O");
    assert_eq!(write("[Na+].[Cl-]"), "Na+.Cl-");
    assert_eq!(write("c1ccc(-c2ccccn2)cc1"), "Ph-2-C5H4N");
    assert_eq!(write("Cc1ccccn1.Cc1ccccn1"), "2-C5H4NCH3.2-C5H4NCH3");
    assert_eq!(write("*C.*C"), "-CH3.-CH3");
}

#[test]
fn round_trips_common_molecules() {
    for smiles in [
        "CCCCCC",
        "CC(C)(C)O",
        "OC(=O)CCC(=O)O",
        "CC(=O)OCC",
        "C[N+](=O)[O-]",
        "Cc1ccc(S(=O)(=O)Cl)cc1",
        "FC(F)(F)C(=O)O",
        "CN(C)C",
        "C#CC",
        "*CC(=O)O",
        "OCC(O)CO",
        "CC(C)CC(C)(C)C",
        "O=CCOC(C=O)=O",
        "c1ccc(-c2ccccn2)cc1",
        "OCc1cccnc1",
        "CC(=O)c1cccs1",
        "c1ccc2cc(C(F)(F)F)ccc2c1",
        "Cc1ccccn1.Cc1ccccn1",
        "*c1ccccn1.*c1ccccn1",
    ] {
        round_trips(smiles);
    }
}

#[test]
fn writes_benzene_rings_with_locants() {
    assert_eq!(write("CCc1ccc(C=O)cc1"), "CH3CH2C6H4-4-CHO");
    assert_eq!(write("Clc1ccc(Cl)cc1"), "1,4-Cl2C6H4");
    assert_eq!(write("c1ccccc1"), "C6H6");
    assert_eq!(write("*c1ccc(*)cc1"), "-1,4-C6H4-");
    assert_eq!(write("*c1cc(C(F)(F)F)cc(C(F)(F)F)c1"), "-C6H3-3,5-(CF3)2");
}

#[test]
fn round_trips_substituted_benzene_rings() {
    for smiles in [
        "Cc1ccc(C=O)cc1",
        "COc1ccc(C=O)cc1",
        "OC(=O)c1ccc(Cl)cc1Cl",
        "Clc1ccc(Cl)cc1",
        "Clc1ccccc1Cl",
        "Cc1cc(C)cc(C)c1",
        "c1ccccc1",
        "*c1ccc(*)cc1",
        "*Cc1ccc(*)cc1",
        "*c1ccc(OC)cc1",
        "*c1cc(C(F)(F)F)cc(C(F)(F)F)c1",
        "CC1(C)OB(c2ccc(Cl)cc2)OC1(C)C",
        "Nc1ccc(N)cc1",
        "Oc1ccc(-c2ccc(O)cc2)cc1",
        "Fc1ccc(P(c2ccc(F)cc2)c2ccc(F)cc2)cc1",
        "O=C(O)c1ccccc1C(=O)O",
        "OC(=O)c1cc(Br)cc(C(=O)O)c1",
        "C[N+](C)(C)c1ccc(OC)cc1",
        "COc1ccc(CCN)cc1OC",
        "Cc1cc(Cl)c(O)c(Cl)c1",
        "CCOC(=O)c1ccc(N)cc1",
        "O=[N+]([O-])c1ccc(Cl)cc1",
        "Cc1cc(ccc1CC=N)N",
        "Cc1cc(cc(c1SC=N)Cl)O",
        "CC(=O)Oc1ccc(C)cc1C=C",
        "CCC(CC(=O)Br)COC",
    ] {
        round_trips(smiles);
    }
}

#[test]
fn writes_template_rings_with_locants() {
    assert_eq!(write("CN(C)c1ccncc1"), "4-C5H4NN(CH3)CH3");
    assert_eq!(write("Cc1cccc(C)n1"), "CH3-2,6-C5H3NCH3");
    assert_eq!(write("Brc1cccc(Br)n1"), "2,6-Br2C5H3N");
    assert_eq!(write("Clc1ccsc1"), "3-C4H3SCl");
    for smiles in [
        "CN(C)c1ccncc1",
        "Cc1cccc(C)n1",
        "OCc1ccccn1",
        "OCc1cccc(CO)n1",
        "c1ccc(P(c2ccccn2)c2ccccn2)nc1",
        "Cc1ccc(C)s1",
        "Cc1ccco1",
        "Clc1ccsc1",
        "*c1ccc(Br)cn1",
        "OCc1cccc2ccccc12",
        "Oc1ccc2ccccc2c1",
        "C(=NNc1ccccc1)c1ccccn1",
        "C[Si](C)(C)Cc1ccccn1",
        "Brc1cccc(Br)n1",
        "CCCCCC(CCCCCC)c1ccc(O)cn1",
        "CCCCCC(CCCCCC)c1c(O)ccc2ccccc12",
        "OC(=O)c1ccc(o1)-c1ccc(cc1)[N+](=O)[O-]",
        "Oc1ccc(Cl)cc1Cc1ccc(C)cn1",
    ] {
        round_trips(smiles);
    }
}

#[test]
fn writes_rings_between_attachment_points() {
    assert_eq!(write("*c1ccc(-c2ccc(*)cc2)cc1"), "-1,4-C6H4-1,4-C6H4-");
    round_trips("*c1ccc(C(=O)c2ccc(*)cc2)cc1");
}

#[test]
fn never_writes_a_formula_that_reads_back_differently() {
    for smiles in [
        "[O-]c1ccc(SC)cc1",
        "[O-]c1ccc(Cl)cc1Cl",
        "OC(=O)c1ccc([O-])cc1[O-]",
        "*c1ccc(C*)cc1[NH3+]",
        "*c1c(Cl)cc(*)cc1Cl",
    ] {
        let molecule = moiety_core::smiles::parse(smiles).unwrap();
        if let Ok(written) = condense(&molecule, default_dictionary()) {
            let top = interpret(&written).top().map(|r| r.smiles.clone());
            assert_eq!(
                top,
                Some(molecule.to_smiles()),
                "{smiles} written as {written}"
            );
        }
    }
    let biphenylene = moiety_core::smiles::parse("Clc1ccc2-c3ccc(Cl)cc3-c2c1").unwrap();
    assert!(matches!(
        condense(&biphenylene, default_dictionary()),
        Err(WriteError::UnknownRing(_))
    ));
}

#[test]
fn rejects_unknown_rings_and_inner_attachments() {
    let ring = moiety_core::smiles::parse("C1CCCC1").unwrap();
    assert!(matches!(
        condense(&ring, default_dictionary()),
        Err(WriteError::UnknownRing(_))
    ));
    assert_eq!(write("*CC(*)C"), "-CH2CH(CH3)-");
    let inner = moiety_core::smiles::parse("*CC(*)C*").unwrap();
    assert!(matches!(
        condense(&inner, default_dictionary()),
        Err(WriteError::InnerAttachment(_))
    ));
}
