//! Aromaticity for the ring systems moiety can express: six-membered rings of
//! carbon and nitrogen, alone or fused.

use crate::element::Element;
use crate::molecule::{BondOrder, Molecule};

/// Marks Kekulé six-membered carbon and nitrogen rings as aromatic.
///
/// A candidate ring qualifies when each of its atoms has exactly one double bond
/// and that double bond stays inside the union of candidate rings. This covers
/// benzene, pyridine and their fused systems, and leaves quinones alone.
pub fn perceive_aromaticity(molecule: &mut Molecule) {
    let adjacency = molecule.adjacency();
    let candidates: Vec<Vec<usize>> = six_rings(molecule, &adjacency)
        .into_iter()
        .filter(|ring| {
            ring.iter().all(|&atom| {
                let a = molecule.atoms()[atom];
                matches!(a.element(), Some(Element::C | Element::N))
            })
        })
        .collect();
    if candidates.is_empty() {
        return;
    }
    let mut in_system = vec![false; molecule.atoms().len()];
    for &atom in candidates.iter().flatten() {
        in_system[atom] = true;
    }
    let qualifies = |atom: usize| {
        let doubles: Vec<usize> = adjacency[atom]
            .iter()
            .filter(|&&(_, bond)| molecule.bonds()[bond].order == BondOrder::Double)
            .map(|&(other, _)| other)
            .collect();
        (doubles.len() == 1 && in_system[doubles[0]])
            || (molecule.atoms()[atom].aromatic && doubles.is_empty())
    };
    let aromatic_rings: Vec<&Vec<usize>> = candidates
        .iter()
        .filter(|ring| ring.iter().all(|&atom| qualifies(atom)))
        .collect();
    for ring in aromatic_rings {
        for (i, &atom) in ring.iter().enumerate() {
            let next = ring[(i + 1) % ring.len()];
            molecule.atoms_mut()[atom].aromatic = true;
            if let Some(&(_, bond)) = adjacency[atom].iter().find(|&&(other, _)| other == next) {
                molecule.bonds_mut()[bond].order = BondOrder::Aromatic;
            }
        }
    }
}

/// A copy of `molecule` with aromatic bonds replaced by alternating single and double bonds.
///
/// Returns `None` when no Kekulé structure exists.
#[must_use]
pub fn kekulize(molecule: &Molecule) -> Option<Molecule> {
    let adjacency = molecule.adjacency();
    let atoms = molecule.atoms();
    let needs_double: Vec<bool> = (0..atoms.len())
        .map(|i| {
            let atom = atoms[i];
            if !atom.aromatic {
                return false;
            }
            let aromatic_bonds = adjacency[i]
                .iter()
                .filter(|&&(_, b)| molecule.bonds()[b].order == BondOrder::Aromatic)
                .count();
            let other_valence: u8 = adjacency[i]
                .iter()
                .filter(|&&(_, b)| molecule.bonds()[b].order != BondOrder::Aromatic)
                .map(|&(_, b)| molecule.bonds()[b].order.valence())
                .sum();
            let used =
                u8::try_from(aromatic_bonds).unwrap_or(u8::MAX) + other_valence + atom.hydrogens;
            let target = atom
                .element()
                .map(|e| e.charged_valences(atom.charge))
                .and_then(|valences| valences.into_iter().find(|&v| v > used))
                .unwrap_or(0);
            target == used + 1
        })
        .collect();
    let mut matched: Vec<Option<usize>> = vec![None; atoms.len()];
    let order: Vec<usize> = (0..atoms.len()).filter(|&i| needs_double[i]).collect();
    if !match_from(0, &order, molecule, &adjacency, &needs_double, &mut matched) {
        return None;
    }
    let mut result = molecule.clone();
    for bond in result.bonds_mut() {
        if bond.order == BondOrder::Aromatic {
            bond.order = if matched[bond.a] == Some(bond.b) {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
        }
    }
    for atom in result.atoms_mut() {
        atom.aromatic = false;
    }
    Some(result)
}

fn match_from(
    index: usize,
    order: &[usize],
    molecule: &Molecule,
    adjacency: &[Vec<(usize, usize)>],
    needs_double: &[bool],
    matched: &mut [Option<usize>],
) -> bool {
    let Some(&atom) = order.get(index) else {
        return true;
    };
    if matched[atom].is_some() {
        return match_from(index + 1, order, molecule, adjacency, needs_double, matched);
    }
    for &(other, bond) in &adjacency[atom] {
        if molecule.bonds()[bond].order == BondOrder::Aromatic
            && needs_double[other]
            && matched[other].is_none()
        {
            matched[atom] = Some(other);
            matched[other] = Some(atom);
            if match_from(index + 1, order, molecule, adjacency, needs_double, matched) {
                return true;
            }
            matched[atom] = None;
            matched[other] = None;
        }
    }
    false
}

/// Every simple cycle of exactly six atoms, each listed once in ring order.
pub(crate) fn six_rings(molecule: &Molecule, adjacency: &[Vec<(usize, usize)>]) -> Vec<Vec<usize>> {
    let mut rings = Vec::new();
    let mut path = Vec::with_capacity(6);
    for start in 0..molecule.atoms().len() {
        path.clear();
        path.push(start);
        extend_ring(start, adjacency, &mut path, &mut rings);
    }
    rings
}

fn extend_ring(
    start: usize,
    adjacency: &[Vec<(usize, usize)>],
    path: &mut Vec<usize>,
    rings: &mut Vec<Vec<usize>>,
) {
    let last = *path.last().expect("path starts non-empty");
    for &(next, _) in &adjacency[last] {
        if path.len() == 6 {
            if next == start && path[1] < path[5] {
                rings.push(path.clone());
            }
            continue;
        }
        if next > start && !path.contains(&next) {
            path.push(next);
            extend_ring(start, adjacency, path, rings);
            path.pop();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles;

    #[test]
    fn kekulizes_benzene_and_pyridine() {
        for input in [
            "c1ccccc1",
            "c1ccncc1",
            "c1ccc2ccccc2c1",
            "Cc1ccc(S(=O)(=O)Cl)cc1",
        ] {
            let molecule = smiles::parse(input).unwrap();
            let kekule = kekulize(&molecule).unwrap();
            assert!(
                kekule
                    .bonds()
                    .iter()
                    .all(|b| b.order != BondOrder::Aromatic)
            );
            let reread = smiles::parse(&smiles::write(&kekule)).unwrap();
            assert_eq!(smiles::write(&reread), smiles::write(&molecule));
        }
    }

    #[test]
    fn leaves_quinone_alone() {
        let molecule = smiles::parse("O=C1C=CC(=O)C=C1").unwrap();
        assert!(molecule.atoms().iter().all(|a| !a.aromatic));
    }

    #[test]
    fn finds_each_six_ring_once() {
        let molecule = smiles::parse("C1CCCCC1").unwrap();
        assert_eq!(six_rings(&molecule, &molecule.adjacency()).len(), 1);
    }
}
