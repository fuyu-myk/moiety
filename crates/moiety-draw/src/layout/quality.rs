//! How legible a drawing of a ring system is, counting the substituents that
//! will hang from it, so that candidate drawings can be compared.

use moiety_core::Molecule;

use super::geometry::{EPSILON, Point, largest_gap};
use crate::graph::{self, Adjacency};
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

/// Largest branch size counted toward the cost of pointing an atom into its ring.
const BRANCH_CAP: usize = 20;

/// Cost of putting an atom where another ring of the system cannot lie outside
/// its ring: the end of a shared path on an inward corner, or the middle of one
/// on an outward corner.
const FUSED: u32 = 400;

/// Bond lengths a drawing may use at all, in bond lengths.
const LENGTHS: std::ops::RangeInclusive<f64> = 0.55..=1.6;

/// How legible a drawing is, worst first: atoms piled on one another and bonds
/// outside [`LENGTHS`], then bonds that cross, then atoms closer than
/// [`CLEARANCE`] to an atom or a stub, then how far bonds miss unit length
/// beyond [`STRETCH`], and last stubs closer than [`CLEARANCE`] to each other,
/// which the placement of substituents can still ease. Smaller is better.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub(super) struct Quality {
    broken: usize,
    crossings: usize,
    crowded: usize,
    stretch: f64,
    cramped: usize,
}

impl Quality {
    /// Whether this drawing is better than `other` leaving stubs aside.
    pub(super) fn beats(&self, other: &Self) -> bool {
        let key = |q: &Self| (q.broken, q.crossings, q.crowded, q.stretch);
        key(self) < key(other)
    }

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
    branches: Vec<usize>,
    /// Stereo double bonds between two atoms of the system: the atoms, a
    /// neighbour of each and whether those neighbours are cis.
    double_bonds: Vec<(usize, usize, usize, usize, bool)>,
}

impl Frame {
    pub(super) fn new(system: &RingSystem, molecule: &Molecule, adjacency: &Adjacency) -> Self {
        let member = |atom: usize| system.atoms.binary_search(&atom).is_ok();
        let double_bonds = molecule
            .double_bond_stereo()
            .iter()
            .map(|d| {
                let bond = molecule.bonds()[d.bond];
                (bond.a, bond.b, d.neighbours.0, d.neighbours.1, d.cis)
            })
            .filter(|&(a, b, ..)| member(a) && member(b))
            .collect();
        let mut inside = vec![Vec::new(); adjacency.len()];
        let mut outside = vec![0; adjacency.len()];
        let mut branches = vec![0; adjacency.len()];
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
                    branches[a] += graph::side(adjacency, b, a).len().min(BRANCH_CAP);
                }
            }
        }
        Self {
            atoms: system.atoms.clone(),
            bonds,
            inside,
            outside,
            branches,
            double_bonds,
        }
    }

    /// The stereo double bonds along `ring`, each as the ring position of its
    /// first atom and whether the ring runs trans across it, the ring neighbours
    /// of its atoms lying on opposite sides.
    pub(super) fn ring_double_bonds(&self, ring: &[usize]) -> Vec<(usize, bool)> {
        let n = ring.len();
        (0..n)
            .filter_map(|i| {
                let (a, b) = (ring[i], ring[(i + 1) % n]);
                let (before, after) = (ring[(i + n - 1) % n], ring[(i + 2) % n]);
                self.double_bonds.iter().find_map(|&(x, y, nx, ny, cis)| {
                    let (na, nb) = if (x, y) == (a, b) {
                        (nx, ny)
                    } else if (x, y) == (b, a) {
                        (ny, nx)
                    } else {
                        return None;
                    };
                    // A neighbour off the ring stands opposite the ring neighbour.
                    let ring_cis = cis ^ (na != before) ^ (nb != after);
                    Some((i, !ring_cis))
                })
            })
            .collect()
    }

    /// For each atom of `ring`, in ring order, what it costs on an inward corner
    /// and on an outward corner of the ring's outline, so that the system's other
    /// `rings` and the substituents lie outside it.
    ///
    /// Where another ring shares a path with this one, the path's ends belong on
    /// outward corners and its middle atoms on inward ones, the way a ring fused
    /// on an edge or hung from three atoms sits outside a macrocycle. An atom
    /// bonded to another ring's atom outside the path also belongs outward, and
    /// a substituted atom costs twenty on an inward corner for each atom of its
    /// substituents, up to [`BRANCH_CAP`] per substituent, so that a large
    /// substituent weighs as much as another ring.
    pub(super) fn corner_costs(&self, ring: &[usize], rings: &[&[usize]]) -> Vec<(u32, u32)> {
        let n = ring.len();
        let mut costs: Vec<(u32, u32)> = ring
            .iter()
            .map(|&atom| {
                let bridging = self.inside[atom].iter().any(|n| !ring.contains(n));
                let branch = 20 * u32::try_from(self.branches[atom]).unwrap_or(u32::MAX / 100);
                (if bridging { FUSED } else { branch }, 0)
            })
            .collect();
        for other in rings.iter().filter(|other| **other != ring) {
            let shared: Vec<bool> = ring.iter().map(|a| other.contains(a)).collect();
            for i in (0..n).filter(|&i| shared[i]) {
                let (before, after) = (shared[(i + n - 1) % n], shared[(i + 1) % n]);
                if before && after {
                    costs[i].1 = costs[i].1.saturating_add(FUSED);
                } else if before || after {
                    costs[i].0 = costs[i].0.saturating_add(FUSED);
                }
            }
        }
        costs
    }

    /// The atom of `ring` with the largest substituents, if any has substituents.
    pub(super) fn heaviest(&self, ring: &[usize]) -> Option<usize> {
        (0..ring.len())
            .filter(|&i| self.branches[ring[i]] > 0)
            .max_by_key(|&i| (self.branches[ring[i]], std::cmp::Reverse(i)))
    }

    /// The atoms of the ring system.
    pub(super) fn atoms(&self) -> &[usize] {
        &self.atoms
    }

    /// The bonds of the ring system.
    pub(super) fn bonds(&self) -> &[(usize, usize)] {
        &self.bonds
    }

    /// The neighbours of `atom` inside the ring system.
    pub(super) fn neighbours(&self, atom: usize) -> &[usize] {
        &self.inside[atom]
    }

    /// Whether the system is small enough for its drawings to be compared.
    pub(super) fn scored(&self) -> bool {
        self.bonds.len() <= MAX_SCORED_BONDS
    }

    /// The quality of the placed part of `local`, stubs included.
    pub(super) fn quality(&self, local: &[Option<Point>]) -> Quality {
        let mut points: Vec<Point> = self.atoms.iter().filter_map(|&a| local[a]).collect();
        let atoms = points.len();
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
        // Pairs closer than `limit`; a pair of stubs counts only when `stubs` is set.
        // Pairs closer than `limit` that `counted` accepts, given whether each is an atom.
        let near = |limit: f64, counted: fn(bool, bool) -> bool| {
            points
                .iter()
                .enumerate()
                .map(|(i, p)| {
                    (i + 1..points.len())
                        .filter(|&j| counted(i < atoms, j < atoms) && p.distance(points[j]) < limit)
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
            broken: near(PILED, |a, b| a || b) + bad_lengths,
            crossings,
            crowded: near(CLEARANCE, |a, b| a || b),
            stretch,
            cramped: near(CLEARANCE, |a, b| !a && !b),
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
