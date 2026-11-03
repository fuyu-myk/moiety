use std::collections::HashMap;

use super::{SmilesError, aromatic_symbol, implicit_hydrogens};
use crate::aromatic::perceive_aromaticity;
use crate::element::Element;
use crate::molecule::{Atom, BondOrder, Molecule};

const MAX_SMILES_CHARGE: i8 = 15;

/// Reads a SMILES string into a molecule.
///
/// Kekulé six-membered rings are perceived as aromatic, so `C1=CC=CC=C1` and
/// `c1ccccc1` read as the same molecule.
///
/// # Errors
///
/// Returns a [`SmilesError`] naming the first byte that does not fit the SMILES grammar.
pub fn parse(smiles: &str) -> Result<Molecule, SmilesError> {
    let mut reader = Reader {
        text: smiles.as_bytes(),
        pos: 0,
        molecule: Molecule::new(),
        implicit: Vec::new(),
        rings: HashMap::new(),
    };
    reader.read()?;
    let mut molecule = reader.finish()?;
    perceive_aromaticity(&mut molecule);
    Ok(molecule)
}

struct Reader<'a> {
    text: &'a [u8],
    pos: usize,
    molecule: Molecule,
    implicit: Vec<bool>,
    rings: HashMap<u32, (usize, Option<BondOrder>, usize)>,
}

impl Reader<'_> {
    fn error(&self, message: impl Into<String>) -> SmilesError {
        SmilesError {
            position: self.pos,
            message: message.into(),
        }
    }

    fn peek(&self) -> Option<u8> {
        self.text.get(self.pos).copied()
    }

    fn read(&mut self) -> Result<(), SmilesError> {
        let mut previous: Option<usize> = None;
        let mut branches: Vec<Option<usize>> = Vec::new();
        let mut bond: Option<BondOrder> = None;
        let mut expect_atom = true;
        while let Some(byte) = self.peek() {
            match byte {
                b'(' => {
                    if previous.is_none() {
                        return Err(self.error("branch opened before any atom"));
                    }
                    branches.push(previous);
                    self.pos += 1;
                }
                b')' => {
                    previous = branches.pop().ok_or_else(|| self.error("unmatched ')'"))?;
                    self.pos += 1;
                }
                b'-' | b'/' | b'\\' => {
                    bond = Some(BondOrder::Single);
                    self.pos += 1;
                }
                b'=' => {
                    bond = Some(BondOrder::Double);
                    self.pos += 1;
                }
                b'#' => {
                    bond = Some(BondOrder::Triple);
                    self.pos += 1;
                }
                b':' => {
                    bond = Some(BondOrder::Aromatic);
                    self.pos += 1;
                }
                b'.' => {
                    if bond.is_some() {
                        return Err(self.error("bond before '.'"));
                    }
                    previous = None;
                    expect_atom = true;
                    self.pos += 1;
                }
                b'0'..=b'9' | b'%' => {
                    let at = previous.ok_or_else(|| self.error("ring closure before any atom"))?;
                    let label = self.ring_label()?;
                    self.close_ring(label, at, bond.take())?;
                }
                _ => {
                    let atom = self.atom()?;
                    if let Some(from) = previous {
                        let order = bond
                            .take()
                            .unwrap_or_else(|| self.implied_order(from, atom));
                        self.molecule.add_bond(from, atom, order);
                    } else if bond.is_some() {
                        return Err(self.error("bond with no atom before it"));
                    }
                    previous = Some(atom);
                    expect_atom = false;
                }
            }
        }
        if !branches.is_empty() {
            return Err(self.error("unclosed '('"));
        }
        if bond.is_some() || (expect_atom && !self.molecule.is_empty()) {
            return Err(self.error("SMILES ends with a dangling bond or '.'"));
        }
        if let Some((_, (_, _, position))) = self.rings.iter().next() {
            return Err(SmilesError {
                position: *position,
                message: "unclosed ring".into(),
            });
        }
        if self.molecule.is_empty() {
            return Err(self.error("empty SMILES"));
        }
        Ok(())
    }

    fn implied_order(&self, a: usize, b: usize) -> BondOrder {
        let atoms = self.molecule.atoms();
        if atoms[a].aromatic && atoms[b].aromatic {
            BondOrder::Aromatic
        } else {
            BondOrder::Single
        }
    }

    fn ring_label(&mut self) -> Result<u32, SmilesError> {
        if self.peek() == Some(b'%') {
            self.pos += 1;
            let digits = self.text.get(self.pos..self.pos + 2).unwrap_or_default();
            if digits.len() != 2 || !digits.iter().all(u8::is_ascii_digit) {
                return Err(self.error("'%' must be followed by two digits"));
            }
            self.pos += 2;
            return Ok(u32::from(digits[0] - b'0') * 10 + u32::from(digits[1] - b'0'));
        }
        let digit = self.text[self.pos] - b'0';
        self.pos += 1;
        Ok(u32::from(digit))
    }

    fn close_ring(
        &mut self,
        label: u32,
        at: usize,
        bond: Option<BondOrder>,
    ) -> Result<(), SmilesError> {
        match self.rings.remove(&label) {
            None => {
                self.rings.insert(label, (at, bond, self.pos));
                Ok(())
            }
            Some((other, opening, _)) => {
                if other == at {
                    return Err(self.error("ring closure bonds an atom to itself"));
                }
                let order = match (opening, bond) {
                    (Some(a), Some(b)) if a != b => {
                        return Err(self.error("ring closure bond orders disagree"));
                    }
                    (Some(order), _) | (None, Some(order)) => order,
                    (None, None) => self.implied_order(other, at),
                };
                self.molecule.add_bond(other, at, order);
                Ok(())
            }
        }
    }

    fn atom(&mut self) -> Result<usize, SmilesError> {
        let byte = self.peek().ok_or_else(|| self.error("expected an atom"))?;
        if byte == b'[' {
            return self.bracket_atom();
        }
        if byte == b'*' {
            self.pos += 1;
            self.implicit.push(false);
            return Ok(self.molecule.add_atom(Atom::attachment()));
        }
        let two = self.text.get(self.pos..self.pos + 2);
        let (symbol, aromatic, len) = match (byte, two) {
            (_, Some(b"Cl")) => ("Cl", false, 2),
            (_, Some(b"Br")) => ("Br", false, 2),
            (b'B', _) => ("B", false, 1),
            (b'C', _) => ("C", false, 1),
            (b'N', _) => ("N", false, 1),
            (b'O', _) => ("O", false, 1),
            (b'P', _) => ("P", false, 1),
            (b'S', _) => ("S", false, 1),
            (b'F', _) => ("F", false, 1),
            (b'I', _) => ("I", false, 1),
            (b'b', _) => ("B", true, 1),
            (b'c', _) => ("C", true, 1),
            (b'n', _) => ("N", true, 1),
            (b'o', _) => ("O", true, 1),
            (b'p', _) => ("P", true, 1),
            (b's', _) => ("S", true, 1),
            _ => return Err(self.error(format!("unexpected '{}'", char::from(byte)))),
        };
        self.pos += len;
        let element = Element::from_symbol(symbol).expect("organic subset symbols are elements");
        self.implicit.push(true);
        Ok(self.molecule.add_atom(Atom {
            aromatic,
            ..Atom::new(element)
        }))
    }

    fn bracket_atom(&mut self) -> Result<usize, SmilesError> {
        self.pos += 1;
        while self.peek().is_some_and(|b| b.is_ascii_digit()) {
            self.pos += 1;
        }
        let atom = if self.peek() == Some(b'*') {
            self.pos += 1;
            Atom::attachment()
        } else {
            let (element, aromatic) = self.bracket_symbol()?;
            Atom {
                aromatic,
                ..Atom::new(element)
            }
        };
        while matches!(self.peek(), Some(b'@')) {
            self.pos += 1;
        }
        let mut hydrogens = 0;
        if self.peek() == Some(b'H') {
            self.pos += 1;
            hydrogens = self.number("hydrogen count")?.unwrap_or(1);
        }
        let charge = self.charge()?;
        if self.peek() == Some(b':') {
            self.pos += 1;
            while self.peek().is_some_and(|b| b.is_ascii_digit()) {
                self.pos += 1;
            }
        }
        if self.peek() != Some(b']') {
            return Err(self.error("expected ']'"));
        }
        self.pos += 1;
        self.implicit.push(false);
        Ok(self.molecule.add_atom(Atom {
            charge,
            hydrogens,
            ..atom
        }))
    }

    fn bracket_symbol(&mut self) -> Result<(Element, bool), SmilesError> {
        let rest = &self.text[self.pos..];
        let first = *rest
            .first()
            .ok_or_else(|| self.error("expected an element symbol"))?;
        if first.is_ascii_lowercase() {
            for len in [2, 1] {
                let Some(candidate) = rest.get(..len) else {
                    continue;
                };
                let mut symbol = String::from_utf8_lossy(candidate).into_owned();
                symbol[..1].make_ascii_uppercase();
                if aromatic_symbol(&symbol) {
                    if let Some(element) = Element::from_symbol(&symbol) {
                        self.pos += len;
                        return Ok((element, true));
                    }
                }
            }
            return Err(self.error("unknown aromatic symbol"));
        }
        if let Some(candidate) = rest.get(..2) {
            if candidate[1].is_ascii_lowercase() {
                if let Some(element) = Element::from_symbol(&String::from_utf8_lossy(candidate)) {
                    self.pos += 2;
                    return Ok((element, false));
                }
            }
        }
        let symbol = String::from_utf8_lossy(&rest[..1]).into_owned();
        let element =
            Element::from_symbol(&symbol).ok_or_else(|| self.error("unknown element symbol"))?;
        self.pos += 1;
        Ok((element, false))
    }

    /// Reads a run of digits as a number below 100, or `None` when there are no digits.
    fn number(&mut self, what: &str) -> Result<Option<u8>, SmilesError> {
        let start = self.pos;
        while self.peek().is_some_and(|b| b.is_ascii_digit()) {
            self.pos += 1;
        }
        if start == self.pos {
            return Ok(None);
        }
        std::str::from_utf8(&self.text[start..self.pos])
            .ok()
            .and_then(|digits| digits.parse::<u8>().ok())
            .filter(|&n| n < 100)
            .map(Some)
            .ok_or_else(|| SmilesError {
                position: start,
                message: format!("{what} out of range"),
            })
    }

    fn charge(&mut self) -> Result<i8, SmilesError> {
        let sign: i8 = match self.peek() {
            Some(b'+') => 1,
            Some(b'-') => -1,
            _ => return Ok(0),
        };
        let sign_byte = self.text[self.pos];
        self.pos += 1;
        let mut magnitude: i8 = 1;
        while self.peek() == Some(sign_byte) {
            self.pos += 1;
            magnitude = magnitude
                .checked_add(1)
                .filter(|&m| m <= MAX_SMILES_CHARGE)
                .ok_or_else(|| self.error("charge out of range"))?;
        }
        if magnitude == 1 {
            if let Some(n) = self.number("charge")? {
                magnitude = i8::try_from(n)
                    .ok()
                    .filter(|&m| m <= MAX_SMILES_CHARGE)
                    .ok_or_else(|| self.error("charge out of range"))?;
            }
        }
        Ok(sign * magnitude)
    }

    fn finish(mut self) -> Result<Molecule, SmilesError> {
        for index in 0..self.molecule.atoms().len() {
            if !self.implicit[index] {
                continue;
            }
            let atom = self.molecule.atoms()[index];
            let bonds = self.molecule.bond_valence(index);
            let element = atom.element().expect("implicit atoms are elements");
            self.molecule.atoms_mut()[index].hydrogens =
                implicit_hydrogens(element.symbol(), bonds, atom.aromatic)
                    .expect("implicit atoms are in the organic subset");
        }
        Ok(self.molecule)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counts(smiles: &str) -> (usize, usize) {
        let molecule = parse(smiles).unwrap();
        (molecule.atoms().len(), molecule.bonds().len())
    }

    #[test]
    fn reads_chains_branches_and_rings() {
        assert_eq!(counts("CCO"), (3, 2));
        assert_eq!(counts("CC(C)C(=O)O"), (6, 5));
        assert_eq!(counts("c1ccccc1"), (6, 6));
        assert_eq!(counts("C1CC%10CC1.C%10"), (6, 6));
    }

    #[test]
    fn computes_implicit_hydrogens() {
        let molecule = parse("CC(=O)[O-]").unwrap();
        let hydrogens: Vec<u8> = molecule.atoms().iter().map(|a| a.hydrogens).collect();
        assert_eq!(hydrogens, vec![3, 0, 0, 0]);
        assert_eq!(molecule.atoms()[3].charge, -1);
        let benzene = parse("c1ccccc1").unwrap();
        assert!(benzene.atoms().iter().all(|a| a.hydrogens == 1));
        for (smiles, hydrogens) in [
            ("s1cccc1", vec![0, 1, 1, 1, 1]),
            ("o1cccc1", vec![0, 1, 1, 1, 1]),
            ("Cn1cccc1", vec![3, 0, 1, 1, 1, 1]),
        ] {
            let counts: Vec<u8> = parse(smiles)
                .unwrap()
                .atoms()
                .iter()
                .map(|a| a.hydrogens)
                .collect();
            assert_eq!(counts, hydrogens, "{smiles}");
        }
    }

    #[test]
    fn reads_bracket_atoms() {
        let molecule = parse("[NH4+].[Cu+2].[13CH3-]").unwrap();
        let atoms = molecule.atoms();
        assert_eq!((atoms[0].hydrogens, atoms[0].charge), (4, 1));
        assert_eq!(atoms[1].charge, 2);
        assert_eq!((atoms[2].hydrogens, atoms[2].charge), (3, -1));
    }

    #[test]
    fn perceives_kekule_benzene() {
        let molecule = parse("C1=CC=CC=C1").unwrap();
        assert!(molecule.atoms().iter().all(|a| a.aromatic));
    }

    #[test]
    fn rejects_out_of_range_numbers() {
        let many_plus = format!("[C{}]", "+".repeat(128));
        assert_eq!(
            parse(&many_plus).unwrap_err().message,
            "charge out of range"
        );
        assert_eq!(parse("[C+300]").unwrap_err().message, "charge out of range");
        assert_eq!(
            parse("[CH300]").unwrap_err().message,
            "hydrogen count out of range"
        );
        assert_eq!(parse("[C+16]").unwrap_err().message, "charge out of range");
    }

    #[test]
    fn rejects_malformed_input() {
        assert_eq!(parse("C(C").unwrap_err().message, "unclosed '('");
        assert_eq!(parse("C1CC").unwrap_err().message, "unclosed ring");
        assert_eq!(parse("CX").unwrap_err().position, 1);
        assert!(parse("").is_err());
        assert!(parse("C=").is_err());
    }
}
