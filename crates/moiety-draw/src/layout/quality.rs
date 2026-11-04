//! How legible a drawing of a ring system is, counting the substituents that
//! will hang from it, so that candidate drawings can be compared.

use super::geometry::{EPSILON, Point, largest_gap};
use crate::graph::Adjacency;
use crate::rings::RingSystem;

/// Atoms of a ring system closer than this, in bond lengths, are crowded.
pub(super) const CLEARANCE: f64 = 0.5;

/// Atoms closer than this, in bond lengths, are piled on one another.
const PILED: f64 = 0.4;

/// A bond that misses unit length by more than this is stretched.
const STRETCH: f64 = 0.25;

/// Most bonds in a ring system whose drawings are scored; larger systems keep
/// their first drawing, since scoring compares every pair of bonds.
const MAX_SCORED_BONDS: usize = 600;

/// Bond lengths a drawing may use at all, in bond lengths.
const LENGTHS: std::ops::RangeInclusive<f64> = 0.55..=1.6;

/// How legible a drawing is, worst first: atoms piled on one another and bonds
/// outside [`LENGTHS`], then bonds that cross, then atoms closer than
/// [`CLEARANCE`], then how far bonds miss unit length beyond [`STRETCH`].
/// Smaller is better.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub(super) struct Quality {
    broken: usize,
    crossings: usize,
    crowded: usize,
    stretch: f64,
}

impl Quality {
    /// The piled atoms, bad bond lengths and crossing bonds, the faults no
    /// placement should accept.
    pub(super) fn severe(&self) -> (usize, usize) {
        (self.broken, self.crossings)
    }

    /// Whether nothing is piled, crossed, crowded or stretched.
    pub(super) fn is_clean(&self) -> bool {
        self.broken == 0 && self.crossings == 0 && self.crowded == 0 && self.stretch < EPSILON
    }
}

/// A ring system's bonds, and for each of its atoms the neighbours inside it and
/// the number outside it, which become stubs one bond long in the atom's widest
/// free angle once its neighbours inside are placed.
pub(super) struct Frame {
    atoms: Vec<usize>,
    bonds: Vec<(usize, usize)>,
    inside: Vec<Vec<usize>>,
    outside: Vec<usize>,
}

impl Frame {
    pub(super) fn new(system: &RingSystem, adjacency: &Adjacency) -> Self {
        let member = |atom: usize| system.atoms.binary_search(&atom).is_ok();
        let mut inside = vec![Vec::new(); adjacency.len()];
        let mut outside = vec![0; adjacency.len()];
        let mut bonds = Vec::new();
        for &a in &system.atoms {
            for &(b, _) in &adjacency[a] {
                if member(b) {
                    inside[a].push(b);
                    if a < b {
                        bonds.push((a, b));
                    }
                } else {
                    outside[a] += 1;
                }
            }
        }
        Self {
            atoms: system.atoms.clone(),
            bonds,
            inside,
            outside,
        }
    }

    /// Whether the system is small enough for its drawings to be compared.
    pub(super) fn scored(&self) -> bool {
        self.bonds.len() <= MAX_SCORED_BONDS
    }

    /// The quality of the placed part of `local`, stubs included.
    pub(super) fn quality(&self, local: &[Option<Point>]) -> Quality {
        let mut points: Vec<Point> = self.atoms.iter().filter_map(|&a| local[a]).collect();
        let mut segments: Vec<(usize, usize, Point, Point)> = self
            .bonds
            .iter()
            .filter_map(|&(a, b)| Some((a, b, local[a]?, local[b]?)))
            .collect();
        let stretch = segments
            .iter()
            .map(|&(_, _, p, q)| ((p.distance(q) - 1.0).abs() - STRETCH).max(0.0))
            .sum();
        let bad_lengths = segments
            .iter()
            .filter(|&&(_, _, p, q)| !LENGTHS.contains(&p.distance(q)))
            .count();
        let mut stub = usize::MAX;
        for &atom in &self.atoms {
            for end in self.stubs(atom, local) {
                points.push(end);
                segments.push((atom, stub, local[atom].unwrap_or_default(), end));
                stub -= 1;
            }
        }
        let near = |limit: f64| {
            points
                .iter()
                .enumerate()
                .map(|(i, p)| {
                    points[i + 1..]
                        .iter()
                        .filter(|q| p.distance(**q) < limit)
                        .count()
                })
                .sum()
        };
        let crossings = segments
            .iter()
            .enumerate()
            .map(|(i, &(a, b, p, q))| {
                segments[i + 1..]
                    .iter()
                    .filter(|&&(c, d, r, s)| {
                        ![c, d].contains(&a) && ![c, d].contains(&b) && cross(p, q, r, s)
                    })
                    .count()
            })
            .sum();
        Quality {
            broken: near(PILED) + bad_lengths,
            crossings,
            crowded: near(CLEARANCE),
            stretch,
        }
    }

    /// Ends of the stubs of `atom`, once it and its neighbours inside are placed.
    fn stubs(&self, atom: usize, local: &[Option<Point>]) -> Vec<Point> {
        let count = self.outside[atom];
        let Some(origin) = local[atom].filter(|_| count > 0) else {
            return Vec::new();
        };
        let directions: Option<Vec<f64>> = self.inside[atom]
            .iter()
            .map(|&n| local[n].map(|p| (p - origin).angle()))
            .collect();
        let Some(directions) = directions else {
            return Vec::new();
        };
        let (start, width) = largest_gap(directions);
        (1..=count)
            .map(|j| origin + Point::from_angle(start + width * j as f64 / (count + 1) as f64))
            .collect()
    }
}

/// Whether the segments `p`-`q` and `r`-`s` cross at a point inside both.
fn cross(p: Point, q: Point, r: Point, s: Point) -> bool {
    let side = |o: Point, a: Point, b: Point| (a - o).cross(b - o);
    let apart = |x: f64, y: f64| (x > EPSILON && y < -EPSILON) || (x < -EPSILON && y > EPSILON);
    apart(side(r, s, p), side(r, s, q)) && apart(side(p, q, r), side(p, q, s))
}
