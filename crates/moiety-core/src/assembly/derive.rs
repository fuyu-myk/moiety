use super::Derivation;
use super::state::State;
use crate::element::Element;
use crate::molecule::{Atom, Bond, BondOrder, Molecule};

pub(crate) fn pentavalent_nitrogen_without_oxo(state: &State) -> Option<usize> {
    (0..state.atoms.len()).find(|&i| {
        let atom = &state.atoms[i];
        atom.element_kind() == Some(Element::N)
            && atom.charge == 0
            && atom.valence == 5
            && !state.bonds.iter().any(|b| {
                b.order == BondOrder::Double
                    && b.other(i)
                        .is_some_and(|o| state.atoms[o].element_kind() == Some(Element::O))
            })
    })
}

/// Builds the final molecule, folding hydrogens into implicit counts.
pub(crate) fn derive(state: State) -> Derivation {
    let mut neighbours = vec![Vec::new(); state.atoms.len()];
    for (index, bond) in state.bonds.iter().enumerate() {
        neighbours[bond.a].push((bond.b, index));
        neighbours[bond.b].push((bond.a, index));
    }
    let collapsible = |i: usize| {
        let atom = &state.atoms[i];
        atom.element_kind() == Some(Element::H)
            && atom.charge == 0
            && neighbours[i].len() == 1
            && state.atoms[neighbours[i][0].0]
                .element_kind()
                .is_some_and(|e| e != Element::H)
    };
    let mut index = vec![None; state.atoms.len()];
    let mut atoms = Vec::new();
    let mut atom_spans = Vec::new();
    for (i, atom) in state.atoms.iter().enumerate() {
        if collapsible(i) {
            continue;
        }
        index[i] = Some(atoms.len());
        atoms.push(Atom {
            kind: atom.kind,
            charge: atom.charge,
            hydrogens: atom.hydrogens,
            aromatic: atom.aromatic,
        });
        atom_spans.push(atom.group);
    }
    for i in 0..state.atoms.len() {
        if collapsible(i) {
            let heavy = index[neighbours[i][0].0].expect("heavy atoms are kept");
            atoms[heavy].hydrogens += 1;
        }
    }
    let mut bonds = Vec::new();
    let mut bond_spans = Vec::new();
    for bond in &state.bonds {
        if let (Some(a), Some(b)) = (index[bond.a], index[bond.b]) {
            bonds.push(Bond {
                a,
                b,
                order: bond.order,
            });
            bond_spans.push(
                bond.token
                    .unwrap_or_else(|| state.atoms[bond.a].group.cover(state.atoms[bond.b].group)),
            );
        }
    }
    let charge_sites = state
        .atoms
        .iter()
        .filter(|a| a.budget_charge)
        .map(|a| a.span)
        .collect();
    Derivation {
        molecule: Molecule::from_parts(atoms, bonds),
        atom_spans,
        bond_spans,
        score: state.score,
        forks: state.forks,
        forward: state.forward,
        backward: state.backward,
        charge_sites,
    }
}
