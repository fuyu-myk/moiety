//! Stereochemistry read from SMILES: the order in which each atom's neighbours
//! are written, which `@` and `@@` refer to, and the direction of bonds written
//! `/` and `\`, which fixes the configuration of the double bonds between them.

use crate::molecule::{BondOrder, DoubleBondStereo, Molecule, Tetrahedral};

/// One entry in the order of an atom's neighbours.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Slot {
    Atom(usize),
    Hydrogen,
    /// A ring bond opened here and not yet closed.
    Pending,
}

/// What the reader has seen of the stereochemistry so far.
#[derive(Debug, Default)]
pub(super) struct Recorder {
    neighbours: Vec<Vec<Slot>>,
    anticlockwise: Vec<Option<bool>>,
    preceded: Vec<bool>,
    /// Bonds written with `/` (`true`) or `\` (`false`), from the atom written first.
    directions: Vec<(usize, usize, bool)>,
}

impl Recorder {
    /// A new atom, with its `@` (`Some(true)`) or `@@` (`Some(false)`) mark and
    /// whether its bracket gives it a hydrogen.
    pub(super) fn atom(&mut self, anticlockwise: Option<bool>, hydrogen: bool) {
        let first = if hydrogen && anticlockwise.is_some() {
            vec![Slot::Hydrogen]
        } else {
            Vec::new()
        };
        self.neighbours.push(first);
        self.anticlockwise.push(anticlockwise);
        self.preceded.push(false);
    }

    /// A bond written from `from` to the atom `to` just after it, with its
    /// direction if it was written `/` or `\`.
    pub(super) fn bond(&mut self, from: usize, to: usize, direction: Option<bool>) {
        self.neighbours[from].push(Slot::Atom(to));
        self.neighbours[to].insert(0, Slot::Atom(from));
        self.preceded[to] = true;
        if let Some(slash) = direction {
            self.directions.push((from, to, slash));
        }
    }

    /// A ring bond opened at `atom`, returning the slot that its partner fills.
    pub(super) fn ring_opened(&mut self, atom: usize) -> usize {
        self.neighbours[atom].push(Slot::Pending);
        self.neighbours[atom].len() - 1
    }

    /// The ring bond opened at `opener`, in `slot`, closed at `closer`, with the
    /// directions written at its opening and closing digits. A direction at a
    /// digit reads as if the partner were written just after that atom.
    pub(super) fn ring_closed(
        &mut self,
        opener: usize,
        slot: usize,
        closer: usize,
        (opened, closed): (Option<bool>, Option<bool>),
    ) {
        self.neighbours[opener][slot] = Slot::Atom(closer);
        self.neighbours[closer].push(Slot::Atom(opener));
        if let Some(slash) = opened {
            self.directions.push((opener, closer, slash));
        }
        if let Some(slash) = closed {
            self.directions.push((closer, opener, slash));
        }
    }

    /// Records in `molecule` every stereocentre with four neighbours and every
    /// double bond with a directional bond on each side.
    pub(super) fn finish(self, molecule: &mut Molecule) {
        for (centre, mark) in self.anticlockwise.iter().enumerate() {
            let Some(anticlockwise) = *mark else {
                continue;
            };
            let slots = &self.neighbours[centre];
            if slots.contains(&Slot::Pending) {
                continue;
            }
            let mut neighbours: Vec<Option<usize>> = slots
                .iter()
                .map(|slot| match slot {
                    Slot::Atom(atom) => Some(*atom),
                    _ => None,
                })
                .collect();
            let lone_pair = molecule.atoms()[centre]
                .element()
                .is_some_and(|e| ["N", "P", "As", "S", "Se", "Te"].contains(&e.symbol()));
            if neighbours.len() == 3 && !slots.contains(&Slot::Hydrogen) && lone_pair {
                // A lone pair stands where an implicit hydrogen would.
                neighbours.insert(usize::from(self.preceded[centre]), None);
            }
            if let Ok(neighbours) = <[Option<usize>; 4]>::try_from(neighbours) {
                molecule.add_tetrahedral(Tetrahedral {
                    centre,
                    neighbours,
                    anticlockwise,
                });
            }
        }
        let bonds: Vec<(usize, usize, usize)> = molecule
            .bonds()
            .iter()
            .enumerate()
            .filter(|(_, bond)| bond.order == BondOrder::Double)
            .map(|(index, bond)| (index, bond.a, bond.b))
            .collect();
        for (bond, a, b) in bonds {
            if let (Some((na, below_a)), Some((nb, below_b))) = (self.side(a, b), self.side(b, a)) {
                molecule.add_double_bond_stereo(DoubleBondStereo {
                    bond,
                    neighbours: (na, nb),
                    cis: below_a == below_b,
                });
            }
        }
    }

    /// The first neighbour of `end`, other than `across`, joined to it by a
    /// directional bond, and whether it lies below `end`: `n/end` and `end\n`
    /// put `n` below, `n\end` and `end/n` above.
    fn side(&self, end: usize, across: usize) -> Option<(usize, bool)> {
        self.directions.iter().find_map(|&(from, to, slash)| {
            if to == end && from != across {
                Some((from, slash))
            } else if from == end && to != across {
                Some((to, !slash))
            } else {
                None
            }
        })
    }
}
