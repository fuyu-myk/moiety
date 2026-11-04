//! Hubs: acyclic atoms that hold three or more alike branches, as the
//! phosphorus of `PPh3` or the boron of `B(C6F5)3`.
//!
//! A component with a hub is drawn outward from it, with its branches evenly
//! spaced and the first pointing up. Branches that crowd each other, and that
//! no flip or turn can part without moving the hub's bonds, are moved apart by
//! stretching every bond of the hub equally, which keeps the angles.

use std::collections::VecDeque;

use moiety_core::{AtomKind, BondOrder, Molecule};

use super::geometry::{EPSILON, Point};
use crate::graph::Adjacency;
use crate::rings::RingInfo;

/// Unbonded atoms of two branches closer than this, in bond lengths, crowd each other.
const CROWDED: f64 = 0.6;

/// Lengths tried for crowded hub bonds, shortest first, in bond lengths.
const STRETCHES: [f64; 2] = [1.25, 1.5];

/// Upper bound on the atom pairs compared while stretching one hub.
const MAX_PAIRS: usize = 4_000_000;

/// A branch's atoms, each by kind, charge, hydrogens, aromaticity, degree and
/// bonds from the branch's first atom, sorted, with its sorted bond orders.
/// Alike branches have equal signatures.
type Signature = (Vec<(usize, AtomKind, i8, u8, bool, usize)>, Vec<BondOrder>);

/// The hub of `component`, if it has one.
///
/// A hub is an atom outside every ring with at least three alike branches of two
/// or more atoms, and no branch holding half the component. Only one atom of a
/// component can have no branch that large, so a component has at most one hub.
pub(super) fn find(
    molecule: &Molecule,
    adjacency: &Adjacency,
    rings: &RingInfo,
    component: &[usize],
) -> Option<usize> {
    let limit = component.len().div_ceil(2);
    component.iter().copied().find_map(|atom| {
        if adjacency[atom].len() < 3 || rings.system_of(atom).is_some() {
            return None;
        }
        let branches = branches(adjacency, atom, limit)?;
        let (count, _) = alike(molecule, adjacency, &branches);
        Some((count >= 3).then_some(atom))
    })?
}

/// The neighbours of `hub` in drawing order: the alike ones in input order, with
/// the others spread evenly among them starting from the first.
pub(super) fn order(molecule: &Molecule, adjacency: &Adjacency, hub: usize) -> Vec<usize> {
    let branches = branches(adjacency, hub, usize::MAX).expect("no branch is over the limit");
    let (_, mask) = alike(molecule, adjacency, &branches);
    let (common, odd): (Vec<usize>, Vec<usize>) = (0..branches.len()).partition(|&i| mask[i]);
    let mut slots: Vec<Option<usize>> = vec![None; branches.len()];
    for (k, &i) in odd.iter().enumerate() {
        slots[k * branches.len() / odd.len()] = Some(i);
    }
    let mut common = common.into_iter();
    slots
        .into_iter()
        .map(|slot| slot.or_else(|| common.next()))
        .map(|slot| branches[slot.expect("every slot is filled")][0])
        .collect()
}

/// Whether `after` leaves no two branches of `hub` crowding each other with the
/// hub and its neighbours where they are in `before`.
pub(super) fn kept(before: &[Point], after: &[Point], adjacency: &Adjacency, hub: usize) -> bool {
    let unmoved = |atom: usize| before[atom].distance(after[atom]) < EPSILON;
    let branches = branches(adjacency, hub, usize::MAX).expect("no branch is over the limit");
    unmoved(hub) && branches.iter().all(|branch| unmoved(branch[0])) && !crowded(after, &branches)
}

/// Stretches every bond of `hub` to the shortest of [`STRETCHES`] that leaves no
/// two of its branches crowding each other, and returns whether it did. The bonds
/// stay as they are when no stretch is enough or none is needed.
pub(super) fn stretch(positions: &mut [Point], adjacency: &Adjacency, hub: usize) -> bool {
    let branches = branches(adjacency, hub, usize::MAX).expect("no branch is over the limit");
    let size: usize = branches.iter().map(Vec::len).sum();
    if size.saturating_mul(size) > MAX_PAIRS || !crowded(positions, &branches) {
        return false;
    }
    let origin = positions[hub];
    let directions: Vec<Point> = branches
        .iter()
        .map(|branch| {
            (positions[branch[0]] - origin)
                .normalized()
                .unwrap_or_default()
        })
        .collect();
    let mut stretched = positions.to_vec();
    let mut length = 1.0;
    for target in STRETCHES {
        for (branch, &direction) in branches.iter().zip(&directions) {
            for &atom in branch {
                stretched[atom] = stretched[atom] + direction * (target - length);
            }
        }
        length = target;
        if !crowded(&stretched, &branches) {
            positions.copy_from_slice(&stretched);
            return true;
        }
    }
    false
}

/// The atoms of each branch of `hub`, its neighbour first, or `None` when one
/// branch reaches `limit` atoms.
fn branches(adjacency: &Adjacency, hub: usize, limit: usize) -> Option<Vec<Vec<usize>>> {
    adjacency[hub]
        .iter()
        .map(|&(start, _)| side_within(adjacency, start, hub, limit))
        .collect()
}

/// The atoms reachable from `start` without passing `blocked`, breadth first, or
/// `None` once they reach `limit`.
fn side_within(
    adjacency: &Adjacency,
    start: usize,
    blocked: usize,
    limit: usize,
) -> Option<Vec<usize>> {
    let mut seen = vec![false; adjacency.len()];
    seen[start] = true;
    seen[blocked] = true;
    let mut members = vec![start];
    let mut queue = VecDeque::from([start]);
    while let Some(atom) = queue.pop_front() {
        for &(next, _) in &adjacency[atom] {
            if !seen[next] {
                seen[next] = true;
                members.push(next);
                queue.push_back(next);
            }
        }
        if members.len() >= limit {
            return None;
        }
    }
    Some(members)
}

/// The size of the largest set of alike branches of two or more atoms, and which
/// branches belong to it.
fn alike(
    molecule: &Molecule,
    adjacency: &Adjacency,
    branches: &[Vec<usize>],
) -> (usize, Vec<bool>) {
    let signatures: Vec<Option<Signature>> = branches
        .iter()
        .map(|branch| (branch.len() >= 2).then(|| signature(molecule, adjacency, branch)))
        .collect();
    let count = |s: &Signature| signatures.iter().flatten().filter(|t| *t == s).count();
    let best = signatures
        .iter()
        .flatten()
        .max_by_key(|s| count(s))
        .cloned();
    let Some(best) = best else {
        return (0, vec![false; branches.len()]);
    };
    let mask = signatures
        .iter()
        .map(|s| s.as_ref() == Some(&best))
        .collect();
    (count(&best), mask)
}

fn signature(molecule: &Molecule, adjacency: &Adjacency, branch: &[usize]) -> Signature {
    let depths = depths(adjacency, branch);
    let mut atoms: Vec<_> = branch
        .iter()
        .zip(depths)
        .map(|(&index, depth)| {
            let atom = molecule.atoms()[index];
            let degree = adjacency[index].len();
            (
                depth,
                atom.kind,
                atom.charge,
                atom.hydrogens,
                atom.aromatic,
                degree,
            )
        })
        .collect();
    atoms.sort_unstable();
    let mut bonds: Vec<BondOrder> = branch
        .iter()
        .flat_map(|&a| adjacency[a].iter().map(|&(_, bond)| bond))
        .map(|bond| molecule.bonds()[bond].order)
        .collect();
    bonds.sort_unstable();
    (atoms, bonds)
}

/// The number of bonds from the first atom of `branch` to each of its atoms,
/// in branch order, counting only paths inside the branch.
fn depths(adjacency: &Adjacency, branch: &[usize]) -> Vec<usize> {
    let mut depth = vec![usize::MAX; adjacency.len()];
    let mut member = vec![false; adjacency.len()];
    for &atom in branch {
        member[atom] = true;
    }
    depth[branch[0]] = 0;
    let mut queue = VecDeque::from([branch[0]]);
    while let Some(atom) = queue.pop_front() {
        for &(next, _) in &adjacency[atom] {
            if member[next] && depth[next] == usize::MAX {
                depth[next] = depth[atom] + 1;
                queue.push_back(next);
            }
        }
    }
    branch.iter().map(|&atom| depth[atom]).collect()
}

/// Whether two atoms of different branches sit closer than [`CROWDED`].
fn crowded(positions: &[Point], branches: &[Vec<usize>]) -> bool {
    branches.iter().enumerate().any(|(i, first)| {
        branches[i + 1..].iter().any(|second| {
            first.iter().any(|&a| {
                second
                    .iter()
                    .any(|&b| positions[a].distance(positions[b]) < CROWDED - EPSILON)
            })
        })
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hub_of(smiles: &str) -> Option<usize> {
        let molecule = moiety_core::smiles::parse(smiles).unwrap();
        let adjacency = molecule.adjacency();
        let rings = RingInfo::new(&molecule, &adjacency);
        let component: Vec<usize> = (0..molecule.atoms().len()).collect();
        find(&molecule, &adjacency, &rings, &component)
    }

    #[test]
    fn finds_atoms_with_three_alike_branches() {
        assert_eq!(hub_of("P(c1ccccc1)(c1ccccc1)c1ccccc1"), Some(0));
        assert_eq!(hub_of("O=P(OC)(OC)OC"), Some(1));
        assert_eq!(hub_of("CC(C)(C)C(=O)Cl"), None, "methyls are too small");
        assert_eq!(hub_of("OC(=O)CC(c1ccccc1)c1ccccc1"), None, "two are alike");
        assert_eq!(
            hub_of("CCCCCCCCCCCCCCCCCCCCCCCCC[Si](CC)(CC)CC"),
            None,
            "a branch holds half the molecule"
        );
        assert_eq!(
            hub_of("P(c1ccccc1C)(c1cccc(C)c1)c1ccc(C)cc1"),
            None,
            "tolyl isomers differ"
        );
        assert_eq!(
            hub_of("P(CCC(C)C)(CC(C)CC)CCC(C)C"),
            None,
            "pentyl isomers differ"
        );
    }

    #[test]
    fn orders_unlike_branches_apart() {
        let molecule = moiety_core::smiles::parse("P(c1ccccc1)(c1ccccc1)(c1ccccc1)=O").unwrap();
        let adjacency = molecule.adjacency();
        assert_eq!(order(&molecule, &adjacency, 0), [19, 1, 7, 13]);
        let molecule = moiety_core::smiles::parse("ClP(Cl)(c1ccccc1)(c1ccccc1)c1ccccc1").unwrap();
        let adjacency = molecule.adjacency();
        assert_eq!(order(&molecule, &adjacency, 1), [0, 3, 2, 9, 15]);
    }
}
