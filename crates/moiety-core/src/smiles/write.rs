use std::collections::{BTreeMap, BTreeSet};

use super::{aromatic_symbol, default_valence};
use crate::molecule::{Atom, AtomKind, BondOrder, Molecule};

/// Writes canonical SMILES.
///
/// Atoms are ranked by iterative refinement of their local invariants, ties are
/// broken one class at a time, and each component is written by a depth-first
/// walk that visits neighbours in rank order. Components are sorted and joined
/// with `.`.
#[must_use]
pub fn write(molecule: &Molecule) -> String {
    let adjacency = molecule.adjacency();
    let ranks = canonical_ranks(molecule, &adjacency);
    let mut parts: Vec<String> = molecule
        .components()
        .into_iter()
        .map(|component| {
            let start = *component
                .iter()
                .min_by_key(|&&atom| ranks[atom])
                .expect("components are non-empty");
            Writer::new(molecule, &adjacency, &ranks).component(start)
        })
        .collect();
    parts.sort();
    parts.join(".")
}

/// Canonical rank of every atom: equal ranks only for symmetry-equivalent atoms.
pub(crate) fn canonical_ranks(
    molecule: &Molecule,
    adjacency: &[Vec<(usize, usize)>],
) -> Vec<usize> {
    let atoms = molecule.atoms();
    let invariants: Vec<(u8, bool, usize, u8, i8)> = atoms
        .iter()
        .enumerate()
        .map(|(i, atom)| {
            let z = match atom.kind {
                AtomKind::Element(element) => element.atomic_number(),
                AtomKind::Attachment => 0,
            };
            (
                z,
                atom.aromatic,
                adjacency[i].len(),
                atom.hydrogens,
                atom.charge,
            )
        })
        .collect();
    let mut ranks = dense_rank(&invariants);
    refine(molecule, adjacency, &mut ranks);
    while let Some(tied) = lowest_tie(&ranks) {
        for rank in &mut ranks {
            *rank = *rank * 2 + 1;
        }
        ranks[tied] -= 1;
        ranks = dense_rank(&ranks);
        refine(molecule, adjacency, &mut ranks);
    }
    ranks
}

fn refine(molecule: &Molecule, adjacency: &[Vec<(usize, usize)>], ranks: &mut Vec<usize>) {
    let mut classes = distinct(ranks);
    loop {
        let keys: Vec<(usize, Vec<(usize, BondOrder)>)> = (0..ranks.len())
            .map(|i| {
                let mut neighbours: Vec<(usize, BondOrder)> = adjacency[i]
                    .iter()
                    .map(|&(other, bond)| (ranks[other], molecule.bonds()[bond].order))
                    .collect();
                neighbours.sort_unstable();
                (ranks[i], neighbours)
            })
            .collect();
        let next = dense_rank(&keys);
        let next_classes = distinct(&next);
        *ranks = next;
        if next_classes == classes {
            return;
        }
        classes = next_classes;
    }
}

fn lowest_tie(ranks: &[usize]) -> Option<usize> {
    let mut counts = BTreeMap::new();
    for &rank in ranks {
        *counts.entry(rank).or_insert(0usize) += 1;
    }
    let (tied_rank, _) = counts.into_iter().find(|&(_, count)| count > 1)?;
    ranks.iter().position(|&rank| rank == tied_rank)
}

fn distinct(ranks: &[usize]) -> usize {
    ranks.iter().collect::<BTreeSet<_>>().len()
}

fn dense_rank<T: Ord>(keys: &[T]) -> Vec<usize> {
    let sorted: Vec<&T> = keys.iter().collect::<BTreeSet<_>>().into_iter().collect();
    keys.iter()
        .map(|key| {
            sorted
                .binary_search(&key)
                .expect("every key is in the sorted set")
        })
        .collect()
}

struct Writer<'a> {
    molecule: &'a Molecule,
    adjacency: &'a [Vec<(usize, usize)>],
    ranks: &'a [usize],
    visited: Vec<bool>,
    children: Vec<Vec<(usize, usize)>>,
    rings: Vec<Vec<(usize, usize)>>,
    ring_seen: Vec<bool>,
    written: Vec<bool>,
    digits: Vec<Option<usize>>,
}

impl<'a> Writer<'a> {
    fn new(
        molecule: &'a Molecule,
        adjacency: &'a [Vec<(usize, usize)>],
        ranks: &'a [usize],
    ) -> Self {
        let atoms = molecule.atoms().len();
        Self {
            molecule,
            adjacency,
            ranks,
            visited: vec![false; atoms],
            children: vec![Vec::new(); atoms],
            rings: vec![Vec::new(); atoms],
            ring_seen: vec![false; molecule.bonds().len()],
            written: vec![false; atoms],
            digits: vec![None; molecule.bonds().len()],
        }
    }

    fn component(mut self, start: usize) -> String {
        self.walk(start, None);
        let mut out = String::new();
        let mut open = Vec::new();
        self.emit(start, &mut out, &mut open);
        out
    }

    fn walk(&mut self, atom: usize, parent: Option<usize>) {
        self.visited[atom] = true;
        let mut neighbours = self.adjacency[atom].clone();
        neighbours.sort_by_key(|&(other, _)| self.ranks[other]);
        for (other, bond) in neighbours {
            if Some(bond) == parent {
                continue;
            }
            if self.visited[other] {
                if !self.ring_seen[bond] {
                    self.ring_seen[bond] = true;
                    self.rings[other].push((atom, bond));
                    self.rings[atom].push((other, bond));
                }
            } else {
                self.children[atom].push((other, bond));
                self.walk(other, Some(bond));
            }
        }
    }

    fn emit(&mut self, atom: usize, out: &mut String, open: &mut Vec<bool>) {
        self.written[atom] = true;
        out.push_str(&self.atom_text(atom));
        let mut rings = self.rings[atom].clone();
        rings.sort_by_key(|&(other, _)| (!self.written[other], self.ranks[other]));
        for (other, bond) in rings {
            if let Some(digit) = self.digits[bond] {
                open[digit] = false;
                push_digit(out, digit);
            } else {
                let digit = open.iter().position(|used| !used).unwrap_or_else(|| {
                    open.push(false);
                    open.len() - 1
                });
                open[digit] = true;
                self.digits[bond] = Some(digit);
                out.push_str(self.bond_text(bond, atom, other));
                push_digit(out, digit);
            }
        }
        let children = self.children[atom].clone();
        let last = children.len().saturating_sub(1);
        for (i, (child, bond)) in children.into_iter().enumerate() {
            if i < last {
                out.push('(');
            }
            out.push_str(self.bond_text(bond, atom, child));
            self.emit(child, out, open);
            if i < last {
                out.push(')');
            }
        }
    }

    fn bond_text(&self, bond: usize, a: usize, b: usize) -> &'static str {
        let atoms = self.molecule.atoms();
        let both_aromatic = atoms[a].aromatic && atoms[b].aromatic;
        match self.molecule.bonds()[bond].order {
            BondOrder::Single if both_aromatic => "-",
            BondOrder::Single => "",
            BondOrder::Double => "=",
            BondOrder::Triple => "#",
            BondOrder::Aromatic if both_aromatic => "",
            BondOrder::Aromatic => ":",
        }
    }

    fn atom_text(&self, index: usize) -> String {
        let atom = self.molecule.atoms()[index];
        let element = match atom.kind {
            AtomKind::Attachment => return "*".into(),
            AtomKind::Element(element) => element,
        };
        let symbol = element.symbol();
        let used = self.molecule.bond_valence(index) + u8::from(atom.aromatic);
        let bare = atom.charge == 0
            && (!atom.aromatic || matches!(symbol, "B" | "C" | "N" | "O" | "P" | "S"))
            && default_valence(symbol, used).is_some_and(|v| v - used == atom.hydrogens);
        let written_symbol = if atom.aromatic && aromatic_symbol(symbol) {
            symbol.to_ascii_lowercase()
        } else {
            symbol.to_string()
        };
        if bare {
            return written_symbol;
        }
        bracket(&written_symbol, &atom)
    }
}

fn bracket(symbol: &str, atom: &Atom) -> String {
    let mut text = format!("[{symbol}");
    match atom.hydrogens {
        0 => {}
        1 => text.push('H'),
        n => text.push_str(&format!("H{n}")),
    }
    match atom.charge {
        0 => {}
        1 => text.push('+'),
        -1 => text.push('-'),
        c if c > 0 => text.push_str(&format!("+{c}")),
        c => text.push_str(&format!("-{}", -i16::from(c))),
    }
    text.push(']');
    text
}

fn push_digit(out: &mut String, digit: usize) {
    let label = digit + 1;
    if label < 10 {
        out.push_str(&label.to_string());
    } else {
        out.push_str(&format!("%{label}"));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse;

    fn canon(smiles: &str) -> String {
        write(&parse(smiles).unwrap())
    }

    #[test]
    fn equivalent_spellings_agree() {
        assert_eq!(canon("OCC"), canon("CCO"));
        assert_eq!(canon("C(C)(C)C(=O)O"), canon("CC(C)C(O)=O"));
        assert_eq!(canon("OC(=O)C(C)C"), canon("CC(C)C(=O)O"));
        assert_eq!(canon("C1=CC=CC=C1C"), canon("Cc1ccccc1"));
        assert_eq!(
            canon("[O-]S(=O)(=O)[O-].[Cu+2]"),
            canon("[Cu+2].O=S(=O)([O-])[O-]")
        );
        assert_eq!(
            canon("Cc1ccc(cc1)S(Cl)(=O)=O"),
            canon("Cc1ccc(S(=O)(=O)Cl)cc1")
        );
    }

    #[test]
    fn different_molecules_differ() {
        assert_ne!(canon("CCO"), canon("COC"));
        assert_ne!(canon("CC(C)CC"), canon("CCCCC"));
        assert_ne!(canon("C=CC=C"), canon("CC=CC"));
    }

    #[test]
    fn writes_brackets_only_when_needed() {
        assert_eq!(canon("[CH3][CH2][OH]"), canon("CCO"));
        assert_eq!(canon("CC(=O)[O-]").matches('[').count(), 1);
        assert!(canon("[NH4+]").starts_with("[NH4+]"));
        assert_eq!(canon("*CC"), canon("CC*"));
    }

    fn permuted(molecule: &Molecule, order: &[usize]) -> Molecule {
        let mut position = vec![0; order.len()];
        for (new, &old) in order.iter().enumerate() {
            position[old] = new;
        }
        let atoms = order.iter().map(|&old| molecule.atoms()[old]).collect();
        let bonds = molecule
            .bonds()
            .iter()
            .rev()
            .map(|bond| crate::molecule::Bond {
                a: position[bond.b],
                b: position[bond.a],
                order: bond.order,
            })
            .collect();
        Molecule::from_parts(atoms, bonds)
    }

    proptest::proptest! {
        #[test]
        fn atom_order_does_not_change_output(
            seed in proptest::sample::select(vec![
                "CC(C)(C)C(=O)OCC",
                "OCC(O)CO",
                "Cc1ccc(S(=O)(=O)Cl)cc1",
                "C1CC2CCC1CC2",
                "c1ccc2ccccc2c1",
                "CC(C)(C)OC(=O)N(C)CC(=O)[O-].[Na+]",
                "*OC(C)=O",
            ]),
            shuffle in proptest::collection::vec(proptest::num::usize::ANY, 32),
        ) {
            let molecule = parse(seed).unwrap();
            let mut order: Vec<usize> = (0..molecule.atoms().len()).collect();
            for (i, r) in shuffle.iter().enumerate().take(order.len()) {
                let j = r % order.len();
                order.swap(i, j);
            }
            proptest::prop_assert_eq!(write(&permuted(&molecule, &order)), write(&molecule));
        }
    }

    #[test]
    fn output_reads_back_to_itself() {
        for input in [
            "CC(C)(C)O",
            "c1ccc2ccccc2c1",
            "C1CC2CCC1CC2",
            "O=C(Cl)Cl",
            "[Na+].[Cl-]",
            "C[N+](=O)[O-]",
            "*OC(C)=O",
        ] {
            let once = canon(input);
            assert_eq!(canon(&once), once, "{input}");
        }
    }
}
