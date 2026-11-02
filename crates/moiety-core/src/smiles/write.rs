use std::cell::Cell;
use std::collections::{BTreeMap, BTreeSet};

use super::{aromatic_symbol, default_valence};
use crate::molecule::{Atom, AtomKind, BondOrder, Molecule};

/// Most leaves explored when breaking symmetry ties in one component.
const TIE_BREAK_LIMIT: usize = 128;

/// Writes canonical SMILES.
///
/// Each component is canonicalized on its own, and identical components only
/// once. Atoms are ranked by iterative refinement of their local invariants.
/// Where refinement leaves a tie, every member of the lowest tied class is
/// tried in turn and the smallest resulting string is kept, so the output does
/// not depend on atom order. Each string is written by a depth-first walk that
/// visits neighbours in rank order. Components are sorted and joined with `.`.
#[must_use]
pub fn write(molecule: &Molecule) -> String {
    write_within(molecule, &Cell::new(usize::MAX)).expect("unlimited work never runs out")
}

/// [`write`], spending one unit of `work` for every atom that refinement or
/// writing visits. `None` once `work` runs out.
pub(crate) fn write_within(molecule: &Molecule, work: &Cell<usize>) -> Option<String> {
    let mut written: Vec<(Molecule, String)> = Vec::new();
    let mut parts: Vec<String> = Vec::new();
    for component in molecule.components() {
        let part = submolecule(molecule, &component);
        let text = match written.iter().find(|(seen, _)| *seen == part) {
            Some((_, text)) => text.clone(),
            None => {
                let text = canonical_component(&part, work)?;
                written.push((part, text.clone()));
                text
            }
        };
        parts.push(text);
    }
    parts.sort();
    Some(parts.join("."))
}

fn spend(work: &Cell<usize>, n: usize) -> bool {
    let left = work.get().checked_sub(n);
    work.set(left.unwrap_or(0));
    left.is_some()
}

fn submolecule(molecule: &Molecule, atoms: &[usize]) -> Molecule {
    let mut index = vec![None; molecule.atoms().len()];
    let mut part = Molecule::new();
    for &atom in atoms {
        index[atom] = Some(part.add_atom(molecule.atoms()[atom]));
    }
    for bond in molecule.bonds() {
        if let (Some(a), Some(b)) = (index[bond.a], index[bond.b]) {
            part.add_bond(a, b, bond.order);
        }
    }
    part
}

fn canonical_component(molecule: &Molecule, work: &Cell<usize>) -> Option<String> {
    let adjacency = molecule.adjacency();
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
    refine(molecule, &adjacency, &mut ranks, work)?;
    let mut best = None;
    let mut leaves = 0;
    search(molecule, &adjacency, ranks, &mut best, &mut leaves, work)?;
    best
}

/// Breaks the lowest remaining tie every possible way and keeps the smallest string.
fn search(
    molecule: &Molecule,
    adjacency: &[Vec<(usize, usize)>],
    ranks: Vec<usize>,
    best: &mut Option<String>,
    leaves: &mut usize,
    work: &Cell<usize>,
) -> Option<()> {
    let Some(tied) = lowest_tie(&ranks) else {
        *leaves += 1;
        if !spend(work, ranks.len()) {
            return None;
        }
        let start = (0..ranks.len())
            .min_by_key(|&atom| (adjacency[atom].len(), ranks[atom]))
            .expect("components are non-empty");
        let written = Writer::new(molecule, adjacency, &ranks).component(start);
        if best.as_ref().is_none_or(|current| written < *current) {
            *best = Some(written);
        }
        return Some(());
    };
    for member in representatives(tied, molecule, adjacency) {
        if *leaves >= TIE_BREAK_LIMIT && best.is_some() {
            return Some(());
        }
        let mut split: Vec<usize> = ranks.iter().map(|&rank| rank * 2 + 1).collect();
        split[member] -= 1;
        let mut split = dense_rank(&split);
        refine(molecule, adjacency, &mut split, work)?;
        search(molecule, adjacency, split, best, leaves, work)?;
    }
    Some(())
}

/// The members of a tied class worth trying. Swapping two terminal atoms on the
/// same neighbour by the same bond is a symmetry of the molecule, so the first
/// of them stands for the rest.
fn representatives(
    tied: Vec<usize>,
    molecule: &Molecule,
    adjacency: &[Vec<(usize, usize)>],
) -> Vec<usize> {
    let mut seen: Vec<(usize, BondOrder)> = Vec::new();
    tied.into_iter()
        .filter(|&atom| match adjacency[atom].as_slice() {
            [(neighbour, bond)] => {
                let key = (*neighbour, molecule.bonds()[*bond].order);
                let first = !seen.contains(&key);
                seen.push(key);
                first
            }
            _ => true,
        })
        .collect()
}

fn refine(
    molecule: &Molecule,
    adjacency: &[Vec<(usize, usize)>],
    ranks: &mut Vec<usize>,
    work: &Cell<usize>,
) -> Option<()> {
    let n = ranks.len();
    let mut classes = ranks.iter().max().map_or(0, |&rank| rank + 1);
    let mut neighbours: Vec<(usize, BondOrder)> = Vec::new();
    let mut starts: Vec<usize> = Vec::with_capacity(n + 1);
    let mut order: Vec<usize> = (0..n).collect();
    loop {
        if !spend(work, n) {
            return None;
        }
        neighbours.clear();
        starts.clear();
        for links in adjacency {
            let from = neighbours.len();
            starts.push(from);
            neighbours.extend(
                links
                    .iter()
                    .map(|&(other, bond)| (ranks[other], molecule.bonds()[bond].order)),
            );
            neighbours[from..].sort_unstable();
        }
        starts.push(neighbours.len());
        let next = {
            let key = |i: usize| (ranks[i], &neighbours[starts[i]..starts[i + 1]]);
            order.sort_unstable_by(|&a, &b| key(a).cmp(&key(b)));
            let mut next = vec![0; n];
            let mut rank = 0;
            for pair in order.windows(2) {
                if key(pair[0]) != key(pair[1]) {
                    rank += 1;
                }
                next[pair[1]] = rank;
            }
            next
        };
        let next_classes = next.iter().max().map_or(0, |&rank| rank + 1);
        *ranks = next;
        if next_classes == classes {
            return Some(());
        }
        classes = next_classes;
    }
}

/// Every atom of the lowest rank shared by more than one atom.
fn lowest_tie(ranks: &[usize]) -> Option<Vec<usize>> {
    let mut counts = BTreeMap::new();
    for &rank in ranks {
        *counts.entry(rank).or_insert(0usize) += 1;
    }
    let (tied_rank, _) = counts.into_iter().find(|&(_, count)| count > 1)?;
    Some(
        (0..ranks.len())
            .filter(|&atom| ranks[atom] == tied_rank)
            .collect(),
    )
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
    fn regular_graphs_that_refinement_cannot_split_agree() {
        assert_eq!(canon("C12C3C1C1C4C(C23)C14"), canon("C12C3C1C3C1C3C2C31"));
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
                "C12C3C1C1C4C(C23)C14",
                "O.O.O.[Cu+2].[O-]S(=O)(=O)[O-]",
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
