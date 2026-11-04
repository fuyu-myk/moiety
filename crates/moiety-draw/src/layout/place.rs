//! Placement of connected components, one atom at a time outward from a start.
//!
//! Ring systems are placed whole from their rigid template. Chains grow as
//! zig-zags: each atom turns the opposite way from the atom before it, so bonds
//! alternate sides at 120 degrees. Branches spread evenly around their atom,
//! substituents on rings take the widest free angle, and atoms with a triple
//! bond or two double bonds are straight. A component with a hub is placed
//! outward from it, its branches evenly spaced and the first pointing up.

use std::collections::VecDeque;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_6, PI, TAU};

use moiety_core::Molecule;

use super::geometry::{EPSILON, Point, largest_gap};
use super::hub;
use super::ring_system::layout_system;
use crate::graph::{self, Adjacency};
use crate::rings::RingInfo;

/// Direction of the first bond of a chain: up and to the right.
const FIRST_BOND: f64 = -FRAC_PI_6;

/// Places atoms component by component.
pub(super) struct Placer<'a> {
    molecule: &'a Molecule,
    adjacency: &'a Adjacency,
    rings: &'a RingInfo,
    positions: Vec<Option<Point>>,
    turns: Vec<i8>,
    queue: VecDeque<usize>,
}

impl<'a> Placer<'a> {
    /// A placer with no atom placed yet.
    pub(super) fn new(
        molecule: &'a Molecule,
        adjacency: &'a Adjacency,
        rings: &'a RingInfo,
    ) -> Self {
        let count = molecule.atoms().len();
        Self {
            molecule,
            adjacency,
            rings,
            positions: vec![None; count],
            turns: vec![0; count],
            queue: VecDeque::new(),
        }
    }

    /// Places every atom of `component` around the origin, and returns its hub if it has one.
    pub(super) fn place(&mut self, component: &[usize]) -> Option<usize> {
        let hub = hub::find(self.molecule, self.adjacency, self.rings, component);
        match (hub, self.start_system(component)) {
            (Some(hub), _) => self.place_hub(hub),
            (None, Some(system)) => self.place_start_system(system),
            (None, None) => {
                let root = self.root(component);
                self.positions[root] = Some(Point::default());
                self.queue.push_back(root);
            }
        }
        while let Some(atom) = self.queue.pop_front() {
            self.expand(atom);
        }
        hub
    }

    /// The position of every atom.
    ///
    /// # Panics
    ///
    /// Panics if an atom was never placed, which cannot happen once every component has been placed.
    pub(super) fn finish(self) -> Vec<Point> {
        self.positions
            .into_iter()
            .map(|point| point.expect("placing a component reaches every atom in it"))
            .collect()
    }

    /// The largest ring system of the component, the one with the lowest atom on ties.
    fn start_system(&self, component: &[usize]) -> Option<usize> {
        let mut best: Option<usize> = None;
        for system in component
            .iter()
            .filter_map(|&atom| self.rings.system_of(atom))
        {
            let size = |s: usize| self.rings.systems()[s].atoms.len();
            if best.is_none_or(|b| size(system) > size(b)) {
                best = Some(system);
            }
        }
        best
    }

    /// One end of a longest path, the lower-numbered end, so chains read in input order.
    fn root(&self, component: &[usize]) -> usize {
        let one_end = graph::farthest(self.adjacency, component[0]);
        let other_end = graph::farthest(self.adjacency, one_end);
        one_end.min(other_end)
    }

    /// Places `hub` at the origin and its neighbours evenly around it, the first straight up.
    fn place_hub(&mut self, hub: usize) {
        self.positions[hub] = Some(Point::default());
        let order = hub::order(self.molecule, self.adjacency, hub);
        let step = TAU / order.len() as f64;
        let slots = order
            .into_iter()
            .enumerate()
            .map(|(j, next)| (next, -FRAC_PI_2 + step * j as f64))
            .collect();
        self.place_slots(hub, slots);
    }

    fn place_start_system(&mut self, system: usize) {
        let atoms = &self.rings.systems()[system].atoms;
        let local = layout_system(self.rings, &self.rings.systems()[system], self.adjacency);
        let single_ring = self.rings.systems()[system].rings.len() == 1;
        let exit = atoms.iter().copied().find(|&atom| {
            self.adjacency[atom]
                .iter()
                .any(|&(next, _)| self.rings.system_of(next) != Some(system))
        });
        let rotation = match exit {
            Some(atom) if single_ring => FIRST_BOND - self.outward(atom, &local),
            _ => 0.0,
        };
        for &atom in atoms {
            self.positions[atom] = local[atom].map(|p| p.rotated(rotation));
            self.queue.push_back(atom);
        }
    }

    fn place_entered_system(&mut self, system: usize, entry: usize, target: Point, from: Point) {
        let local = layout_system(self.rings, &self.rings.systems()[system], self.adjacency);
        let anchor = local[entry].expect("the entry atom belongs to the system");
        let rotation = (from - target).angle() - self.outward(entry, &local);
        for &atom in &self.rings.systems()[system].atoms {
            self.positions[atom] = local[atom].map(|p| target + (p - anchor).rotated(rotation));
            self.queue.push_back(atom);
        }
    }

    /// Direction from `atom` into the widest free angle between its neighbours placed in `local`.
    fn outward(&self, atom: usize, local: &[Option<Point>]) -> f64 {
        let origin = local[atom].unwrap_or_default();
        let (start, width) = largest_gap(
            self.adjacency[atom]
                .iter()
                .filter_map(|&(next, _)| local[next])
                .map(|p| (p - origin).angle()),
        );
        start + width / 2.0
    }

    /// Places every unplaced neighbour of `atom`.
    fn expand(&mut self, atom: usize) {
        let pending: Vec<usize> = self
            .neighbours(atom)
            .filter(|&n| self.positions[n].is_none())
            .collect();
        if pending.is_empty() {
            return;
        }
        let placed: Vec<usize> = self
            .neighbours(atom)
            .filter(|&n| self.positions[n].is_some())
            .collect();
        let slots = match placed.as_slice() {
            [] => self.root_slots(&pending),
            [parent] if self.rings.system_of(atom).is_none() => {
                self.chain_slots(atom, *parent, &pending)
            }
            _ => self.gap_slots(atom, &placed, &pending),
        };
        self.place_slots(atom, slots);
    }

    /// Places each neighbour of `atom` one bond away in the direction given with it.
    fn place_slots(&mut self, atom: usize, slots: Vec<(usize, f64)>) {
        let origin = self.position(atom);
        for (next, angle) in slots {
            let target = origin + Point::from_angle(angle);
            match self.rings.system_of(next) {
                Some(system) => self.place_entered_system(system, next, target, origin),
                None => {
                    self.positions[next] = Some(target);
                    self.queue.push_back(next);
                }
            }
        }
    }

    fn root_slots(&self, pending: &[usize]) -> Vec<(usize, f64)> {
        let first = match pending {
            [only] if self.rings.system_of(*only).is_none() && self.prefers_axis(*only) => 0.0,
            _ => FIRST_BOND,
        };
        let step = TAU / pending.len() as f64;
        pending
            .iter()
            .enumerate()
            .map(|(j, &next)| (next, first + step * j as f64))
            .collect()
    }

    /// Whether a chain should reach `atom` along the horizontal: it is straight, or a four-way branch.
    fn prefers_axis(&self, atom: usize) -> bool {
        self.is_linear(atom) || self.adjacency[atom].len() >= 4
    }

    fn chain_slots(&mut self, atom: usize, parent: usize, pending: &[usize]) -> Vec<(usize, f64)> {
        let incoming = (self.position(atom) - self.position(parent)).angle();
        let turn = match self.turns[parent] {
            0 => horizontal_turn(incoming),
            previous => -previous,
        };
        let sign = f64::from(turn);
        let (angles, own_turn) = match pending.len() + 1 {
            2 if self.is_linear(atom) => (vec![incoming], self.turns[parent]),
            2 => (vec![incoming + sign * FRAC_PI_3], turn),
            3 => (
                vec![incoming + sign * FRAC_PI_3, incoming - sign * FRAC_PI_3],
                turn,
            ),
            4 => (
                vec![
                    incoming,
                    incoming - sign * FRAC_PI_2,
                    incoming + sign * FRAC_PI_2,
                ],
                -turn,
            ),
            degree => (spread(incoming, degree), -turn),
        };
        self.turns[atom] = own_turn;
        self.main_first(atom, pending)
            .into_iter()
            .zip(angles)
            .collect()
    }

    fn gap_slots(&self, atom: usize, placed: &[usize], pending: &[usize]) -> Vec<(usize, f64)> {
        let origin = self.position(atom);
        let (start, width) =
            largest_gap(placed.iter().map(|&n| (self.position(n) - origin).angle()));
        let parts = (pending.len() + 1) as f64;
        pending
            .iter()
            .enumerate()
            .map(|(j, &next)| (next, start + width * (j + 1) as f64 / parts))
            .collect()
    }

    /// `pending` with the neighbour leading to the most atoms first, so it continues the main chain.
    fn main_first(&self, atom: usize, pending: &[usize]) -> Vec<usize> {
        let mut ordered: Vec<(usize, usize)> = pending
            .iter()
            .map(|&next| (graph::side(self.adjacency, next, atom).len(), next))
            .collect();
        ordered.sort_by_key(|&(size, _)| std::cmp::Reverse(size));
        ordered.into_iter().map(|(_, next)| next).collect()
    }

    fn is_linear(&self, atom: usize) -> bool {
        graph::is_linear(self.molecule, self.adjacency, atom)
    }

    fn neighbours(&self, atom: usize) -> impl Iterator<Item = usize> + '_ {
        self.adjacency[atom].iter().map(|&(next, _)| next)
    }

    fn position(&self, atom: usize) -> Point {
        self.positions[atom].expect("only placed atoms are expanded or measured")
    }
}

/// The turn, `1` or `-1`, that keeps a chain arriving along `incoming` closest to horizontal.
fn horizontal_turn(incoming: f64) -> i8 {
    let score = |angle: f64| (angle.sin().abs(), -angle.cos());
    let (left, right) = (score(incoming + FRAC_PI_3), score(incoming - FRAC_PI_3));
    if (left.0 - right.0).abs() > EPSILON {
        if left.0 < right.0 { 1 } else { -1 }
    } else if left.1 > right.1 + EPSILON {
        -1
    } else {
        1
    }
}

/// Directions for the `degree - 1` children of an atom reached along `incoming`,
/// evenly spaced, the one nearest to straight ahead first.
fn spread(incoming: f64, degree: usize) -> Vec<f64> {
    let step = TAU / degree as f64;
    let mut slots: Vec<f64> = (1..degree)
        .map(|j| incoming + PI + step * j as f64)
        .collect();
    let ahead = slots
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| deviation(**a, incoming).total_cmp(&deviation(**b, incoming)))
        .map_or(0, |(i, _)| i);
    let main = slots.remove(ahead);
    slots.insert(0, main);
    slots
}

/// Absolute angle between two directions, in `[0, pi]`.
fn deviation(a: f64, b: f64) -> f64 {
    let difference = (a - b).rem_euclid(TAU);
    difference.min(TAU - difference)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn turns_toward_the_horizontal() {
        assert_eq!(horizontal_turn(-FRAC_PI_6), 1);
        assert_eq!(horizontal_turn(FRAC_PI_6), -1);
        assert_eq!(horizontal_turn(PI - FRAC_PI_6), 1);
        assert_eq!(horizontal_turn(0.0), 1);
    }

    #[test]
    fn spreads_five_neighbours_evenly() {
        let slots = spread(0.0, 5);
        assert_eq!(slots.len(), 4);
        assert!(deviation(slots[0], 0.0) < 0.7);
        assert!((deviation(PI, 0.0) - PI).abs() < EPSILON);
    }
}
