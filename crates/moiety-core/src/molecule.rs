use crate::element::Element;
use crate::smiles;

/// What an atom position holds.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum AtomKind {
    /// A real atom.
    Element(Element),
    /// An open valence, written `*` in SMILES and `-` at either end of a condensed formula.
    Attachment,
}

/// One atom of a [`Molecule`]. Hydrogens bonded to a single heavy atom are stored as a count.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Atom {
    /// Element, or an attachment point.
    pub kind: AtomKind,
    /// Formal charge.
    pub charge: i8,
    /// Number of attached hydrogens not stored as atoms.
    pub hydrogens: u8,
    /// Whether the atom belongs to an aromatic ring.
    pub aromatic: bool,
}

impl Atom {
    /// A neutral, non-aromatic atom of `element` with no hydrogens.
    #[must_use]
    pub const fn new(element: Element) -> Self {
        Self {
            kind: AtomKind::Element(element),
            charge: 0,
            hydrogens: 0,
            aromatic: false,
        }
    }

    /// An attachment point.
    #[must_use]
    pub const fn attachment() -> Self {
        Self {
            kind: AtomKind::Attachment,
            charge: 0,
            hydrogens: 0,
            aromatic: false,
        }
    }

    /// The element, or `None` for an attachment point.
    #[must_use]
    pub const fn element(&self) -> Option<Element> {
        match self.kind {
            AtomKind::Element(element) => Some(element),
            AtomKind::Attachment => None,
        }
    }
}

/// Order of a bond.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum BondOrder {
    /// Single bond.
    Single,
    /// Double bond.
    Double,
    /// Triple bond.
    Triple,
    /// Bond inside an aromatic ring.
    Aromatic,
}

impl BondOrder {
    /// The order for a bond of `n` electron pairs, if `n` is 1, 2 or 3.
    #[must_use]
    pub const fn from_count(n: u8) -> Option<BondOrder> {
        match n {
            1 => Some(BondOrder::Single),
            2 => Some(BondOrder::Double),
            3 => Some(BondOrder::Triple),
            _ => None,
        }
    }

    /// Valence consumed at each end. An aromatic bond counts 1; the aromatic atom adds the rest.
    #[must_use]
    pub const fn valence(self) -> u8 {
        match self {
            BondOrder::Single | BondOrder::Aromatic => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
        }
    }
}

/// A bond between two atoms of a [`Molecule`], by atom index.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Bond {
    /// First atom.
    pub a: usize,
    /// Second atom.
    pub b: usize,
    /// Bond order.
    pub order: BondOrder,
}

impl Bond {
    /// The atom at the other end from `atom`, if `atom` is an endpoint.
    #[must_use]
    pub const fn other(&self, atom: usize) -> Option<usize> {
        if self.a == atom {
            Some(self.b)
        } else if self.b == atom {
            Some(self.a)
        } else {
            None
        }
    }
}

/// A molecular graph with implicit hydrogens.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Molecule {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
}

impl Molecule {
    /// An empty molecule.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Builds a molecule from atoms and bonds.
    ///
    /// # Panics
    ///
    /// Panics if a bond refers to an atom index out of range or bonds an atom to itself.
    #[must_use]
    pub fn from_parts(atoms: Vec<Atom>, bonds: Vec<Bond>) -> Self {
        for bond in &bonds {
            assert!(
                bond.a < atoms.len() && bond.b < atoms.len() && bond.a != bond.b,
                "bond {bond:?} does not join two distinct atoms"
            );
        }
        Self { atoms, bonds }
    }

    /// Appends an atom and returns its index.
    pub fn add_atom(&mut self, atom: Atom) -> usize {
        self.atoms.push(atom);
        self.atoms.len() - 1
    }

    /// Appends a bond between atoms `a` and `b` and returns its index.
    ///
    /// # Panics
    ///
    /// Panics if either index is out of range or `a == b`.
    pub fn add_bond(&mut self, a: usize, b: usize, order: BondOrder) -> usize {
        assert!(
            a < self.atoms.len() && b < self.atoms.len() && a != b,
            "bond {a}-{b} does not join two distinct atoms"
        );
        self.bonds.push(Bond { a, b, order });
        self.bonds.len() - 1
    }

    /// The atoms, by index.
    #[must_use]
    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }

    /// Mutable access to the atoms, by index.
    pub fn atoms_mut(&mut self) -> &mut [Atom] {
        &mut self.atoms
    }

    /// The bonds, by index.
    #[must_use]
    pub fn bonds(&self) -> &[Bond] {
        &self.bonds
    }

    /// Mutable access to the bonds, by index.
    pub fn bonds_mut(&mut self) -> &mut [Bond] {
        &mut self.bonds
    }

    /// Whether the molecule has no atoms.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// For each atom, its neighbours as `(atom index, bond index)` pairs in bond order.
    #[must_use]
    pub fn adjacency(&self) -> Vec<Vec<(usize, usize)>> {
        let mut adjacency = vec![Vec::new(); self.atoms.len()];
        for (index, bond) in self.bonds.iter().enumerate() {
            adjacency[bond.a].push((bond.b, index));
            adjacency[bond.b].push((bond.a, index));
        }
        adjacency
    }

    /// Connected components, each a sorted list of atom indices, ordered by their first atom.
    #[must_use]
    pub fn components(&self) -> Vec<Vec<usize>> {
        let adjacency = self.adjacency();
        let mut seen = vec![false; self.atoms.len()];
        let mut components = Vec::new();
        for start in 0..self.atoms.len() {
            if seen[start] {
                continue;
            }
            seen[start] = true;
            let mut stack = vec![start];
            let mut members = Vec::new();
            while let Some(atom) = stack.pop() {
                members.push(atom);
                for &(next, _) in &adjacency[atom] {
                    if !seen[next] {
                        seen[next] = true;
                        stack.push(next);
                    }
                }
            }
            members.sort_unstable();
            components.push(members);
        }
        components
    }

    /// Sum of bond valences at `atom`, counting aromatic bonds as 1 each.
    #[must_use]
    pub fn bond_valence(&self, atom: usize) -> u8 {
        self.bonds
            .iter()
            .filter(|bond| bond.a == atom || bond.b == atom)
            .map(|bond| bond.order.valence())
            .sum()
    }

    /// Appends every atom and bond of `other`, returning the index offset applied to its atoms.
    pub fn extend(&mut self, other: &Molecule) -> usize {
        let offset = self.atoms.len();
        self.atoms.extend_from_slice(&other.atoms);
        self.bonds.extend(other.bonds.iter().map(|bond| Bond {
            a: bond.a + offset,
            b: bond.b + offset,
            order: bond.order,
        }));
        offset
    }

    /// Net formal charge.
    #[must_use]
    pub fn charge(&self) -> i32 {
        self.atoms.iter().map(|atom| i32::from(atom.charge)).sum()
    }

    /// Canonical SMILES. Two molecules have the same canonical SMILES exactly when they are the same graph.
    #[must_use]
    pub fn to_smiles(&self) -> String {
        smiles::write(self)
    }
}
