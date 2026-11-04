//! Ring templates: ring systems that locants number, such as pyridine, with
//! their numbering (grammar, sections 6.12 and 9). Benzene is built in.

use std::sync::OnceLock;

use crate::dictionary::DictionaryError;
use crate::element::Element;
use crate::molecule::{AtomKind, Molecule};
use crate::smiles;

/// Most positions a template may number, the most a one-digit locant can reach.
pub const MAX_POSITIONS: usize = 8;

/// A ring system that locants number, with its parent molecule and numbering.
///
/// The ring atoms are numbered from 1 in SMILES order, skipping atoms shared by
/// two rings, and a template ring is written as the parent's formula in Hill
/// order with the hydrogens left after substitution, as `C5H4N` for pyridine.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RingTemplate {
    name: String,
    smiles: String,
    parent: Molecule,
    positions: Vec<usize>,
    formula: Vec<(Element, usize)>,
}

impl RingTemplate {
    /// Defines a template from the SMILES of its parent ring system.
    ///
    /// # Errors
    ///
    /// Fails when the SMILES does not parse, has an attachment point, or is not one
    /// connected ring system with no bond outside its rings, when a numbered atom
    /// holds more than one hydrogen, when an atom shared by two rings holds one,
    /// when no numbered atom holds one, when the parent has no carbon, or when it
    /// numbers more than [`MAX_POSITIONS`] atoms.
    pub fn new(name: &str, smiles: &str) -> Result<Self, DictionaryError> {
        let bad = |reason: &str| DictionaryError::BadRing {
            ring: name.to_string(),
            reason: reason.to_string(),
        };
        let parent = smiles::parse(smiles).map_err(|error| bad(&error.to_string()))?;
        if parent
            .atoms()
            .iter()
            .any(|a| a.kind == AtomKind::Attachment)
        {
            return Err(bad("a ring template has no attachment points"));
        }
        if parent.components().len() != 1 || !every_bond_in_a_ring(&parent) {
            return Err(bad(
                "a ring template is one ring system with no bonds outside its rings",
            ));
        }
        let adjacency = parent.adjacency();
        let positions: Vec<usize> = (0..parent.atoms().len())
            .filter(|&atom| adjacency[atom].len() == 2)
            .collect();
        if positions.len() > MAX_POSITIONS {
            return Err(bad("a ring template numbers at most 8 positions"));
        }
        if positions
            .iter()
            .any(|&atom| parent.atoms()[atom].hydrogens > 1)
        {
            return Err(bad("a numbered atom may hold at most one hydrogen"));
        }
        if (0..parent.atoms().len())
            .any(|atom| adjacency[atom].len() != 2 && parent.atoms()[atom].hydrogens > 0)
        {
            return Err(bad("an atom shared by two rings may not hold a hydrogen"));
        }
        if positions
            .iter()
            .all(|&atom| parent.atoms()[atom].hydrogens == 0)
        {
            return Err(bad("a ring template needs a numbered atom with a hydrogen"));
        }
        let mut formula: Vec<(Element, usize)> = Vec::new();
        for atom in parent.atoms() {
            let AtomKind::Element(element) = atom.kind else {
                continue;
            };
            match formula.iter_mut().find(|(e, _)| *e == element) {
                Some((_, count)) => *count += 1,
                None => formula.push((element, 1)),
            }
        }
        formula.sort_by_key(|(element, _)| (*element != Element::C, element.symbol()));
        if formula
            .first()
            .is_none_or(|(element, _)| *element != Element::C)
        {
            return Err(bad(
                "a ring template needs carbon, which its formula starts with",
            ));
        }
        Ok(Self {
            name: name.to_string(),
            smiles: smiles.to_string(),
            parent,
            positions,
            formula,
        })
    }

    /// The built-in benzene template.
    #[must_use]
    pub fn benzene() -> &'static RingTemplate {
        static BENZENE: OnceLock<RingTemplate> = OnceLock::new();
        BENZENE.get_or_init(|| Self::new("benzene", "c1ccccc1").expect("benzene is a ring"))
    }

    /// Template name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The parent ring system as written.
    #[must_use]
    pub fn smiles(&self) -> &str {
        &self.smiles
    }

    /// The parent ring system.
    #[must_use]
    pub fn parent(&self) -> &Molecule {
        &self.parent
    }

    /// How many positions the template numbers.
    #[must_use]
    pub fn positions(&self) -> usize {
        self.positions.len()
    }

    /// The parent atom at `position`, counting from 1.
    #[must_use]
    pub fn atom(&self, position: u8) -> Option<usize> {
        let index = usize::from(position).checked_sub(1)?;
        self.positions.get(index).copied()
    }

    /// The position of a parent atom, if it is numbered.
    #[must_use]
    pub fn position(&self, atom: usize) -> Option<u8> {
        let index = self.positions.iter().position(|&a| a == atom)?;
        u8::try_from(index + 1).ok()
    }

    /// Whether a substituent can take `position`: the parent atom there holds a hydrogen.
    #[must_use]
    pub fn substitutable(&self, position: u8) -> bool {
        self.atom(position)
            .is_some_and(|atom| self.parent.atoms()[atom].hydrogens == 1)
    }

    /// The hydrogens of the parent, all on numbered positions.
    #[must_use]
    pub fn hydrogens(&self) -> usize {
        self.positions
            .iter()
            .map(|&atom| usize::from(self.parent.atoms()[atom].hydrogens))
            .sum()
    }

    /// The heavy elements of the parent and their counts, carbon first and then
    /// alphabetical, as the formula writes them around the hydrogen count.
    #[must_use]
    pub fn formula(&self) -> &[(Element, usize)] {
        &self.formula
    }

    /// The ring written with `hydrogens` hydrogens left, as `C5H4N`.
    #[must_use]
    pub fn written(&self, hydrogens: usize) -> String {
        let count = |n: usize| if n == 1 { String::new() } else { n.to_string() };
        let mut text = String::new();
        for (index, &(element, n)) in self.formula.iter().enumerate() {
            text.push_str(element.symbol());
            text.push_str(&count(n));
            if index == 0 && hydrogens > 0 {
                text.push('H');
                text.push_str(&count(hydrogens));
            }
        }
        text
    }
}

/// Whether every bond of `molecule` lies on a cycle.
fn every_bond_in_a_ring(molecule: &Molecule) -> bool {
    let adjacency = molecule.adjacency();
    molecule.bonds().iter().enumerate().all(|(index, bond)| {
        let mut seen = vec![false; molecule.atoms().len()];
        let mut stack = vec![bond.a];
        seen[bond.a] = true;
        while let Some(atom) = stack.pop() {
            for &(next, via) in &adjacency[atom] {
                if via != index && !seen[next] {
                    seen[next] = true;
                    stack.push(next);
                }
            }
        }
        seen[bond.b]
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numbers_parent_rings() {
        let pyridine = RingTemplate::new("pyridine", "n1ccccc1").unwrap();
        assert_eq!(pyridine.positions(), 6);
        assert!(!pyridine.substitutable(1));
        assert!(pyridine.substitutable(2));
        assert_eq!(pyridine.hydrogens(), 5);
        assert_eq!(pyridine.written(4), "C5H4N");
        assert_eq!(pyridine.written(1), "C5HN");
        let naphthalene = RingTemplate::new("naphthalene", "c1cccc2ccccc12").unwrap();
        assert_eq!(naphthalene.positions(), 8);
        assert_eq!(naphthalene.written(7), "C10H7");
        assert_eq!(RingTemplate::benzene().written(4), "C6H4");
        assert_eq!(
            RingTemplate::new("thiophene", "s1cccc1")
                .unwrap()
                .written(3),
            "C4H3S"
        );
    }

    #[test]
    fn rejects_what_is_not_a_ring_system() {
        for smiles in [
            "*c1ccccc1",
            "CCO",
            "c1ccccc1-c1ccccc1",
            "C1CCCCC1",
            "n1nnnn1",
            "C12C=CC(C=C1)C=C2",
        ] {
            assert!(
                matches!(
                    RingTemplate::new("x", smiles),
                    Err(DictionaryError::BadRing { .. })
                ),
                "{smiles}"
            );
        }
    }
}
