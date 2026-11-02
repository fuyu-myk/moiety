//! Clash removal: acyclic branches are flipped or turned about single bonds
//! until no two unbonded atoms sit closer than [`CLASH`] bond lengths, or no
//! move helps.

use std::f64::consts::{FRAC_PI_3, FRAC_PI_6};

use moiety_core::{BondOrder, Molecule};

use super::geometry::Point;
use crate::graph::{self, Adjacency, bond_between, shortest_path};
use crate::rings::RingInfo;

/// Unbonded atoms closer than this, in bond lengths, clash.
const CLASH: f64 = 0.6;

/// Upper bound on the moves made for one component.
const MAX_MOVES: usize = 64;

/// Turns tried about a branch point, after the flip.
const TURNS: [f64; 4] = [FRAC_PI_6, -FRAC_PI_6, FRAC_PI_3, -FRAC_PI_3];

/// Penalty added to a turn when comparing it with a flip, which keeps every bond angle.
const TURN_COST: f64 = 0.02;

/// How a branch moves about its pivot.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Motion {
    /// Mirror the branch across the bond from the pivot.
    Flip,
    /// Rotate the branch about the pivot by this many radians.
    Turn(f64),
}

/// One candidate move of the atoms `branch` about `pivot`, along the bond to `axis`.
struct Move {
    branch: Vec<usize>,
    pivot: usize,
    axis: usize,
    motion: Motion,
}

/// The read-only context of a clash search.
struct Scene<'a> {
    molecule: &'a Molecule,
    adjacency: &'a Adjacency,
    rings: &'a RingInfo,
    component: &'a [usize],
}

/// Moves branches of `component` until none of its atoms clash or no move reduces the clashes.
pub(super) fn resolve(
    positions: &mut [Point],
    component: &[usize],
    molecule: &Molecule,
    adjacency: &Adjacency,
    rings: &RingInfo,
) {
    let scene = Scene {
        molecule,
        adjacency,
        rings,
        component,
    };
    for _ in 0..MAX_MOVES {
        let clashes = scene.clashes(positions);
        let Some(best) = clashes
            .iter()
            .find_map(|&(a, b)| scene.best_move(positions, a, b))
        else {
            return;
        };
        let moved = best.apply(positions);
        for (&atom, point) in best.branch.iter().zip(moved) {
            positions[atom] = point;
        }
    }
}

impl Move {
    /// New positions of the branch atoms, in branch order.
    fn apply(&self, positions: &[Point]) -> Vec<Point> {
        let pivot = positions[self.pivot];
        let direction = (positions[self.axis] - pivot)
            .normalized()
            .unwrap_or(Point::new(1.0, 0.0));
        self.branch
            .iter()
            .map(|&atom| match self.motion {
                Motion::Flip => positions[atom].reflected(pivot, direction),
                Motion::Turn(angle) => pivot + (positions[atom] - pivot).rotated(angle),
            })
            .collect()
    }
}

impl Scene<'_> {
    /// Every clashing pair of unbonded atoms, closest first.
    fn clashes(&self, positions: &[Point]) -> Vec<(usize, usize)> {
        let mut pairs: Vec<(f64, usize, usize)> = Vec::new();
        for (i, &a) in self.component.iter().enumerate() {
            for &b in &self.component[i + 1..] {
                let distance = positions[a].distance(positions[b]);
                if distance < CLASH && bond_between(self.adjacency, a, b).is_none() {
                    pairs.push((distance, a, b));
                }
            }
        }
        pairs.sort_by(|x, y| x.0.total_cmp(&y.0).then((x.1, x.2).cmp(&(y.1, y.2))));
        pairs.into_iter().map(|(_, a, b)| (a, b)).collect()
    }

    /// The move about a bond between `a` and `b` that reduces the clash penalty most, if any does.
    fn best_move(&self, positions: &[Point], a: usize, b: usize) -> Option<Move> {
        let path = shortest_path(self.adjacency, a, b, None)?;
        let mut best: Option<(f64, Move)> = None;
        for pair in path.windows(2) {
            for candidate in self.moves_about(pair[0], pair[1]) {
                let change = self.penalty_change(positions, &candidate);
                let cost = match candidate.motion {
                    Motion::Flip => change,
                    Motion::Turn(_) => change + TURN_COST,
                };
                if change < -1e-9 && best.as_ref().is_none_or(|(c, _)| cost < *c) {
                    best = Some((cost, candidate));
                }
            }
        }
        best.map(|(_, candidate)| candidate)
    }

    /// Candidate moves of the smaller side of the bond `x`-`y`, if it is an acyclic single bond.
    fn moves_about(&self, x: usize, y: usize) -> Vec<Move> {
        let Some(bond) = bond_between(self.adjacency, x, y) else {
            return Vec::new();
        };
        if self.rings.is_ring_bond(bond) || self.molecule.bonds()[bond].order != BondOrder::Single {
            return Vec::new();
        }
        let y_side = graph::side(self.adjacency, y, x);
        let x_side = graph::side(self.adjacency, x, y);
        let (branch, pivot, axis) = if x_side.len() < y_side.len() {
            (x_side, y, x)
        } else {
            (y_side, x, y)
        };
        let turns = if graph::is_linear(self.molecule, self.adjacency, pivot) {
            &[][..]
        } else {
            &TURNS[..]
        };
        std::iter::once(Motion::Flip)
            .chain(turns.iter().map(|&angle| Motion::Turn(angle)))
            .map(|motion| Move {
                branch: branch.clone(),
                pivot,
                axis,
                motion,
            })
            .collect()
    }

    /// Change in the total clash penalty if `candidate` were applied.
    fn penalty_change(&self, positions: &[Point], candidate: &Move) -> f64 {
        let mut moving = vec![false; positions.len()];
        for &atom in &candidate.branch {
            moving[atom] = true;
        }
        let moved = candidate.apply(positions);
        let mut change = 0.0;
        for (&atom, &after) in candidate.branch.iter().zip(&moved) {
            for &other in self.component.iter().filter(|&&o| !moving[o]) {
                if bond_between(self.adjacency, atom, other).is_none() {
                    change += penalty(after.distance(positions[other]))
                        - penalty(positions[atom].distance(positions[other]));
                }
            }
        }
        change
    }
}

fn penalty(distance: f64) -> f64 {
    (CLASH - distance).max(0.0).powi(2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flips_a_branch_off_a_clash() {
        let molecule = moiety_core::smiles::parse("CCCC").unwrap();
        let adjacency = molecule.adjacency();
        let rings = RingInfo::new(&molecule, &adjacency);
        let mut positions = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(1.0, 1.0),
            Point::new(0.2, 0.4),
        ];
        assert!(positions[0].distance(positions[3]) < CLASH);
        resolve(&mut positions, &[0, 1, 2, 3], &molecule, &adjacency, &rings);
        assert!(positions[0].distance(positions[3]) >= CLASH);
        assert!((positions[2].distance(positions[3]) - 1.0).abs() < 1e-9);
        assert!(
            (positions[3].x - 1.8).abs() < 1e-9,
            "the branch is flipped, not turned"
        );
    }

    #[test]
    fn penalises_only_close_pairs() {
        assert_eq!(penalty(1.0), 0.0);
        assert!(penalty(0.1) > penalty(0.5));
    }
}
