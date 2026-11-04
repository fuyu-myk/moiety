//! Ring perception: a cycle basis of smallest rings, grouped into ring systems.
//!
//! Every bond that lies on a cycle contributes the shortest cycle through it.
//! These candidates are taken smallest first and kept when they are linearly
//! independent of the rings already kept (over GF(2), as sets of bonds). Rings
//! that share an atom belong to the same ring system.

use moiety_core::Molecule;

use crate::graph::{Adjacency, bond_between, shortest_path};

/// The rings of a molecule and the ring systems they form.
#[derive(Clone, Debug, Default)]
pub(crate) struct RingInfo {
    rings: Vec<Vec<usize>>,
    systems: Vec<RingSystem>,
    atom_system: Vec<Option<usize>>,
    ring_bonds: Vec<bool>,
    corners: Vec<(usize, usize, usize)>,
}

/// Rings that share atoms, laid out together as one rigid unit.
#[derive(Clone, Debug, Default)]
pub(crate) struct RingSystem {
    /// Indices into [`RingInfo::rings`].
    pub(crate) rings: Vec<usize>,
    /// Every atom of the system, sorted.
    pub(crate) atoms: Vec<usize>,
}

impl RingInfo {
    /// Perceives the rings of `molecule`.
    pub(crate) fn new(molecule: &Molecule, adjacency: &Adjacency) -> Self {
        let mut ring_bonds = vec![false; molecule.bonds().len()];
        let mut candidates: Vec<Vec<usize>> = Vec::new();
        for (index, bond) in molecule.bonds().iter().enumerate() {
            if let Some(path) = shortest_path(adjacency, bond.a, bond.b, Some(index)) {
                ring_bonds[index] = true;
                candidates.push(canonical(path));
            }
        }
        candidates.sort_by_cached_key(|ring| {
            let mut atoms = ring.clone();
            atoms.sort_unstable();
            (ring.len(), atoms)
        });
        candidates.dedup();
        let rings = independent(candidates, adjacency, molecule.bonds().len());
        let systems = group(&rings, molecule.atoms().len());
        let mut atom_system = vec![None; molecule.atoms().len()];
        for (index, system) in systems.iter().enumerate() {
            for &atom in &system.atoms {
                atom_system[atom] = Some(index);
            }
        }
        let corners = fused_corners(&rings, &systems);
        Self {
            rings,
            systems,
            atom_system,
            ring_bonds,
            corners,
        }
    }

    /// The rings, each listed once as atom indices in ring order.
    pub(crate) fn rings(&self) -> &[Vec<usize>] {
        &self.rings
    }

    /// The ring systems.
    pub(crate) fn systems(&self) -> &[RingSystem] {
        &self.systems
    }

    /// The ring system containing `atom`, if the atom is in a ring.
    pub(crate) fn system_of(&self, atom: usize) -> Option<usize> {
        self.atom_system.get(atom).copied().flatten()
    }

    /// The corners where a ring fills the angle at an atom, sorted, each as
    /// `(atom, lower neighbour, higher neighbour)`. Only rings that share at
    /// most two atoms with every other ring count: a bridged ring may be drawn
    /// in perspective, where its smallest rings are not the faces drawn.
    pub(crate) fn corners(&self) -> &[(usize, usize, usize)] {
        &self.corners
    }

    /// Whether a ring runs from `a` through `atom` to `b`, filling the angle
    /// between them, as [`RingInfo::corners`] counts them.
    pub(crate) fn fills(&self, a: usize, atom: usize, b: usize) -> bool {
        corner_between(&self.corners, a, atom, b)
    }

    /// Whether bond `index` lies on a cycle.
    pub(crate) fn is_ring_bond(&self, index: usize) -> bool {
        self.ring_bonds.get(index).copied().unwrap_or(false)
    }
}

/// Whether the sorted `corners` hold the corner at `atom` between `a` and `b`.
pub(crate) fn corner_between(
    corners: &[(usize, usize, usize)],
    a: usize,
    atom: usize,
    b: usize,
) -> bool {
    corners.binary_search(&(atom, a.min(b), a.max(b))).is_ok()
}

/// The corners of the rings that share at most two atoms with every other
/// ring of their system, sorted, as [`RingInfo::corners`] describes them.
fn fused_corners(rings: &[Vec<usize>], systems: &[RingSystem]) -> Vec<(usize, usize, usize)> {
    let shared = |a: &[usize], b: &[usize]| a.iter().filter(|atom| b.contains(atom)).count();
    let mut corners: Vec<(usize, usize, usize)> = systems
        .iter()
        .flat_map(|system| {
            system.rings.iter().filter(|&&i| {
                system
                    .rings
                    .iter()
                    .all(|&j| i == j || shared(&rings[i], &rings[j]) <= 2)
            })
        })
        .flat_map(|&i| {
            let ring = &rings[i];
            let n = ring.len();
            (0..n).map(move |k| {
                let (before, after) = (ring[(k + n - 1) % n], ring[(k + 1) % n]);
                (ring[k], before.min(after), before.max(after))
            })
        })
        .collect();
    corners.sort_unstable();
    corners
}

/// The ring rotated to start at its lowest atom, running toward the lower of that atom's two ring neighbours.
fn canonical(mut ring: Vec<usize>) -> Vec<usize> {
    let lowest = ring
        .iter()
        .enumerate()
        .min_by_key(|&(_, &atom)| atom)
        .map_or(0, |(position, _)| position);
    ring.rotate_left(lowest);
    if ring.len() > 2 && ring[1] > ring[ring.len() - 1] {
        ring[1..].reverse();
    }
    ring
}

/// The candidates, in order, that are independent of the ones kept before them.
fn independent(
    candidates: Vec<Vec<usize>>,
    adjacency: &Adjacency,
    bond_count: usize,
) -> Vec<Vec<usize>> {
    let words = bond_count.div_ceil(64);
    let mut basis: Vec<(usize, Vec<u64>)> = Vec::new();
    let mut kept = Vec::new();
    for ring in candidates {
        let mut bits = vec![0u64; words];
        for (i, &atom) in ring.iter().enumerate() {
            let next = ring[(i + 1) % ring.len()];
            if let Some(bond) = bond_between(adjacency, atom, next) {
                bits[bond / 64] ^= 1 << (bond % 64);
            }
        }
        for (pivot, row) in &basis {
            if bits[pivot / 64] & (1 << (pivot % 64)) != 0 {
                bits.iter_mut().zip(row).for_each(|(bit, &r)| *bit ^= r);
            }
        }
        let Some(pivot) = bits
            .iter()
            .enumerate()
            .find(|&(_, &word)| word != 0)
            .map(|(i, &word)| i * 64 + word.trailing_zeros() as usize)
        else {
            continue;
        };
        basis.push((pivot, bits));
        kept.push(ring);
    }
    kept
}

/// Groups rings that share an atom into systems, in order of their first ring.
fn group(rings: &[Vec<usize>], atom_count: usize) -> Vec<RingSystem> {
    let mut parent: Vec<usize> = (0..rings.len()).collect();
    let mut owner: Vec<Option<usize>> = vec![None; atom_count];
    for (index, ring) in rings.iter().enumerate() {
        for &atom in ring {
            match owner[atom] {
                Some(other) => union(&mut parent, other, index),
                None => owner[atom] = Some(index),
            }
        }
    }
    let mut systems: Vec<RingSystem> = Vec::new();
    let mut slot: Vec<Option<usize>> = vec![None; rings.len()];
    for (index, ring) in rings.iter().enumerate() {
        let root = find(&mut parent, index);
        let position = *slot[root].get_or_insert_with(|| {
            systems.push(RingSystem::default());
            systems.len() - 1
        });
        let system = &mut systems[position];
        system.rings.push(index);
        system.atoms.extend_from_slice(ring);
    }
    for system in &mut systems {
        system.atoms.sort_unstable();
        system.atoms.dedup();
    }
    systems
}

fn find(parent: &mut [usize], mut index: usize) -> usize {
    while parent[index] != index {
        parent[index] = parent[parent[index]];
        index = parent[index];
    }
    index
}

fn union(parent: &mut [usize], a: usize, b: usize) {
    let (ra, rb) = (find(parent, a), find(parent, b));
    if ra != rb {
        parent[ra.max(rb)] = ra.min(rb);
    }
}

#[cfg(test)]
mod tests {
    use moiety_core::smiles;

    use super::*;

    fn perceive(input: &str) -> RingInfo {
        let molecule = smiles::parse(input).unwrap();
        RingInfo::new(&molecule, &molecule.adjacency())
    }

    #[test]
    fn finds_a_single_ring_in_order() {
        let info = perceive("c1ccccc1C");
        assert_eq!(info.rings(), [vec![0, 1, 2, 3, 4, 5]]);
        assert_eq!(info.systems().len(), 1);
        assert_eq!(info.system_of(6), None);
        assert!(info.is_ring_bond(0));
        assert!(!info.is_ring_bond(6));
    }

    #[test]
    fn fills_the_corners_of_fused_rings_only() {
        let decalin = perceive("C1CCC2CCCCC2C1");
        assert!(decalin.fills(2, 3, 8) && decalin.fills(8, 3, 4));
        assert!(!decalin.fills(2, 3, 4) && !decalin.fills(0, 3, 4));
        let norbornane = perceive("C1CC2CCC1C2");
        assert!(!norbornane.fills(1, 2, 3) && !norbornane.fills(1, 2, 6));
        assert_eq!(norbornane.corners().len(), 0);
    }

    #[test]
    fn splits_naphthalene_into_two_hexagons() {
        let info = perceive("c1ccc2ccccc2c1");
        assert_eq!(info.rings().len(), 2);
        assert!(info.rings().iter().all(|ring| ring.len() == 6));
        assert_eq!(info.systems().len(), 1);
        assert_eq!(info.systems()[0].atoms.len(), 10);
    }

    #[test]
    fn keeps_separate_rings_in_separate_systems() {
        let info = perceive("c1ccccc1C(c1ccccc1)c1ccccc1");
        assert_eq!(info.rings().len(), 3);
        assert_eq!(info.systems().len(), 3);
    }

    #[test]
    fn finds_enough_rings_for_cage_compounds() {
        let cubane = perceive("C12C3C4C1C5C2C3C45");
        assert_eq!(cubane.rings().len(), 5);
        assert!(cubane.rings().iter().all(|ring| ring.len() == 4));
        let spiro = perceive("C1CCC2(CC1)CCCC2");
        assert_eq!(spiro.rings().len(), 2);
        assert_eq!(spiro.systems().len(), 1);
    }
}
