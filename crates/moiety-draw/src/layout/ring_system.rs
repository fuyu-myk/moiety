//! Rigid layouts of ring systems in local coordinates.
//!
//! The first ring is a regular polygon. Every further ring is built onto the
//! atoms already placed: a ring fused on an edge becomes the regular polygon on
//! the far side of that edge, a spiro ring a polygon pointing away from the
//! system, and the unplaced atoms of any other ring an arc of unit chords
//! between the placed atoms that bound them. Systems that this cannot draw
//! without atoms piling up, such as cages, are relaxed by stress majorization.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use super::geometry::{EPSILON, Point, centroid, circumradius};
use super::relax::relax;
use crate::graph::Adjacency;
use crate::rings::{RingInfo, RingSystem};

/// Atoms of a ring system closer than this, in bond lengths, are crowded.
const CLEARANCE: f64 = 0.5;

/// A template bond that misses unit length by more than this is stretched.
const STRETCH: f64 = 0.25;

/// Local positions of the atoms of `system`, indexed by atom; atoms outside the system are `None`.
pub(super) fn layout_system(
    info: &RingInfo,
    system: &RingSystem,
    adjacency: &Adjacency,
) -> Vec<Option<Point>> {
    let mut local = template(info, system, adjacency.len());
    let distance = |a: usize, b: usize| match (local[a], local[b]) {
        (Some(p), Some(q)) => p.distance(q),
        _ => f64::INFINITY,
    };
    let crowded = system.atoms.iter().enumerate().any(|(i, &a)| {
        system.atoms[i + 1..]
            .iter()
            .any(|&b| distance(a, b) < CLEARANCE)
    });
    let stretched = system.atoms.iter().any(|&a| {
        adjacency[a]
            .iter()
            .any(|&(b, _)| local[b].is_some() && (distance(a, b) - 1.0).abs() > STRETCH)
    });
    if crowded || stretched {
        relax(&system.atoms, adjacency, &mut local);
    }
    local
}

/// Local positions from ring templates alone.
fn template(info: &RingInfo, system: &RingSystem, atom_count: usize) -> Vec<Option<Point>> {
    let rings: Vec<&[usize]> = system
        .rings
        .iter()
        .map(|&index| info.rings()[index].as_slice())
        .collect();
    let mut local = vec![None; atom_count];
    let mut done = vec![false; rings.len()];
    let first = first_ring(&rings);
    place_first(&rings, first, &mut local);
    done[first] = true;
    while let Some(next) = next_ring(&rings, &done, &local) {
        place_ring(rings[next], &mut local);
        done[next] = true;
    }
    local
}

/// The ring fused to the most others, then the largest, then the earliest.
fn first_ring(rings: &[&[usize]]) -> usize {
    (0..rings.len())
        .max_by_key(|&i| {
            let fused = (0..rings.len())
                .filter(|&j| j != i && shared(rings[i], rings[j]) >= 2)
                .count();
            (fused, rings[i].len(), std::cmp::Reverse(i))
        })
        .expect("a ring system has at least one ring")
}

fn shared(a: &[usize], b: &[usize]) -> usize {
    a.iter().filter(|atom| b.contains(atom)).count()
}

/// Places `rings[first]` centred on the origin, turned so that its first fused
/// edge is vertical on the right, or pointing up when nothing is fused to it.
fn place_first(rings: &[&[usize]], first: usize, local: &mut [Option<Point>]) {
    let ring = rings[first];
    let n = ring.len();
    let step = TAU / n as f64;
    let fused_edge = (0..n).find(|&i| {
        let (a, b) = (ring[i], ring[(i + 1) % n]);
        rings
            .iter()
            .enumerate()
            .any(|(j, other)| j != first && other.contains(&a) && other.contains(&b))
    });
    let base = match fused_edge {
        Some(i) => -(i as f64 + 0.5) * step,
        None => upright_base(n),
    };
    let radius = circumradius(n);
    for (i, &atom) in ring.iter().enumerate() {
        local[atom] = Some(Point::from_angle(base + i as f64 * step) * radius);
    }
}

/// Angle of the first vertex of a free-standing ring: a vertex on top, or a flat top for squares and octagons.
fn upright_base(sides: usize) -> f64 {
    let flat = if sides % 4 == 0 {
        PI / sides as f64
    } else {
        0.0
    };
    flat - FRAC_PI_2
}

/// The unplaced ring with the most atoms already placed, the earliest on ties.
fn next_ring(rings: &[&[usize]], done: &[bool], local: &[Option<Point>]) -> Option<usize> {
    (0..rings.len()).filter(|&i| !done[i]).max_by_key(|&i| {
        let placed = rings[i].iter().filter(|&&a| local[a].is_some()).count();
        (placed, std::cmp::Reverse(i))
    })
}

fn place_ring(ring: &[usize], local: &mut [Option<Point>]) {
    let placed: Vec<usize> = (0..ring.len())
        .filter(|&i| local[ring[i]].is_some())
        .collect();
    let system_centre = centroid(local.iter().flatten().copied()).unwrap_or_default();
    match placed.as_slice() {
        [] => place_apart(ring, local),
        [only] => place_spiro(ring, *only, system_centre, local),
        _ => {
            let reference = if placed.len() > 2 {
                centroid(placed.iter().filter_map(|&i| local[ring[i]])).unwrap_or_default()
            } else {
                system_centre
            };
            for (start, end) in gaps(ring.len(), &placed) {
                fill_gap(ring, start, end, reference, local);
            }
        }
    }
}

/// Every run of unplaced ring positions as `(placed position before, placed position after)`.
fn gaps(n: usize, placed: &[usize]) -> Vec<(usize, usize)> {
    placed
        .iter()
        .zip(placed.iter().skip(1).chain(placed.first()))
        .filter(|&(&from, &to)| (to + n - from) % n != 1)
        .map(|(&from, &to)| (from, to))
        .collect()
}

/// Places the atoms strictly between ring positions `start` and `end` on an arc.
///
/// The arc bulges away from `reference` unless that brings it closer than
/// [`CLEARANCE`] to atoms already placed and the other side does not, which is
/// how a bridge across a ring ends up inside it.
fn fill_gap(
    ring: &[usize],
    start: usize,
    end: usize,
    reference: Point,
    local: &mut [Option<Point>],
) {
    let n = ring.len();
    let count = (end + n - start) % n - 1;
    let (Some(from), Some(to)) = (local[ring[start]], local[ring[end]]) else {
        return;
    };
    let away = away_from(from, to, reference);
    let clearance = |points: &[Point]| {
        points
            .iter()
            .flat_map(|p| local.iter().flatten().map(move |q| p.distance(*q)))
            .fold(f64::INFINITY, f64::min)
    };
    let outer = arc(from, to, count, away);
    let points = if clearance(&outer) >= CLEARANCE {
        outer
    } else {
        let inner = arc(from, to, count, -away);
        if clearance(&inner) > clearance(&outer) {
            inner
        } else {
            outer
        }
    };
    for (offset, point) in points.into_iter().enumerate() {
        local[ring[(start + 1 + offset) % n]] = Some(point);
    }
}

/// The unit normal of the chord from `from` to `to` on the side away from `reference`.
fn away_from(from: Point, to: Point, reference: Point) -> Point {
    let normal = (to - from)
        .normalized()
        .map_or(Point::new(0.0, -1.0), Point::perpendicular);
    if (reference - (from + to) * 0.5).dot(normal) > EPSILON {
        -normal
    } else {
        normal
    }
}

/// `count` points joining `from` to `to` with `count + 1` equal chords, as close
/// to length 1 as the distance allows, bulging toward the unit vector `normal`.
fn arc(from: Point, to: Point, count: usize, normal: Point) -> Vec<Point> {
    let chords = (count + 1) as f64;
    let middle = (from + to) * 0.5;
    let chord = to - from;
    let length = chord.length();
    if length < EPSILON || length >= chords - EPSILON {
        let along = if length < EPSILON {
            normal
        } else {
            chord * chords.recip()
        };
        return (1..=count).map(|j| from + along * j as f64).collect();
    }
    let sweep = solve_sweep(chords, length);
    let radius = 0.5 / (sweep / (2.0 * chords)).sin();
    let centre = middle - normal * (radius * (sweep / 2.0).cos());
    let (start, end) = (from - centre, to - centre);
    let miss = |sign: f64| (start.rotated(sign * sweep) - end).length();
    let sign = if miss(1.0) < miss(-1.0) { 1.0 } else { -1.0 };
    (1..=count)
        .map(|j| centre + start.rotated(sign * sweep * j as f64 / chords))
        .collect()
}

/// The central angle in `(0, 2 pi)` over which `chords` unit chords span a straight distance `length`.
fn solve_sweep(chords: f64, length: f64) -> f64 {
    let span = |sweep: f64| (sweep / 2.0).sin() / (sweep / (2.0 * chords)).sin();
    let (mut low, mut high) = (0.0, TAU);
    for _ in 0..80 {
        let mid = 0.5 * (low + high);
        if span(mid) > length {
            low = mid;
        } else {
            high = mid;
        }
    }
    0.5 * (low + high)
}

/// Places a ring sharing one atom with the system as a polygon pointing away from the system.
fn place_spiro(
    ring: &[usize],
    shared_at: usize,
    system_centre: Point,
    local: &mut [Option<Point>],
) {
    let n = ring.len();
    let anchor = local[ring[shared_at]].unwrap_or_default();
    let outward = (anchor - system_centre)
        .normalized()
        .unwrap_or(Point::new(1.0, 0.0));
    let radius = circumradius(n);
    let centre = anchor + outward * radius;
    let base = (-outward).angle();
    for j in 1..n {
        let point = centre + Point::from_angle(base + TAU * j as f64 / n as f64) * radius;
        local[ring[(shared_at + j) % n]] = Some(point);
    }
}

/// Places a ring that shares nothing with what is placed to the right of it.
fn place_apart(ring: &[usize], local: &mut [Option<Point>]) {
    let n = ring.len();
    let radius = circumradius(n);
    let right = local
        .iter()
        .flatten()
        .map(|point| point.x)
        .fold(f64::NEG_INFINITY, f64::max);
    let centre = Point::new(
        if right.is_finite() {
            right + radius + 1.0
        } else {
            0.0
        },
        0.0,
    );
    let base = upright_base(n);
    for (i, &atom) in ring.iter().enumerate() {
        local[atom] = Some(centre + Point::from_angle(base + TAU * i as f64 / n as f64) * radius);
    }
}

#[cfg(test)]
mod tests {
    use moiety_core::smiles;

    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-6
    }

    fn system_layout(input: &str) -> (moiety_core::Molecule, Vec<Option<Point>>) {
        let molecule = smiles::parse(input).unwrap();
        let adjacency = molecule.adjacency();
        let info = RingInfo::new(&molecule, &adjacency);
        let local = layout_system(&info, &info.systems()[0], &adjacency);
        (molecule, local)
    }

    fn bond_lengths(molecule: &moiety_core::Molecule, local: &[Option<Point>]) -> Vec<f64> {
        molecule
            .bonds()
            .iter()
            .filter_map(|b| Some(local[b.a]?.distance(local[b.b]?)))
            .collect()
    }

    fn closest_pair(local: &[Option<Point>]) -> f64 {
        let points: Vec<Point> = local.iter().flatten().copied().collect();
        points
            .iter()
            .enumerate()
            .flat_map(|(i, p)| points[i + 1..].iter().map(move |q| p.distance(*q)))
            .fold(f64::INFINITY, f64::min)
    }

    #[test]
    fn arcs_rebuild_a_regular_hexagon_on_an_edge() {
        let (top, bottom) = (Point::new(0.0, -0.5), Point::new(0.0, 0.5));
        let centre = Point::new(3f64.sqrt() / 2.0, 0.0);
        for (from, to) in [(top, bottom), (bottom, top)] {
            let points = arc(from, to, 4, Point::new(1.0, 0.0));
            assert_eq!(points.len(), 4);
            assert!(points.iter().all(|p| close(p.distance(centre), 1.0)));
            assert!(close(points[0].distance(from), 1.0));
            assert!(points.windows(2).all(|w| close(w[0].distance(w[1]), 1.0)));
            assert!(close(points[3].distance(to), 1.0));
        }
        assert_eq!(
            away_from(top, bottom, Point::new(-1.0, 0.0)),
            Point::new(1.0, 0.0)
        );
    }

    #[test]
    fn arcs_fall_back_to_a_line_when_the_gap_is_too_wide() {
        let points = arc(
            Point::new(0.0, 0.0),
            Point::new(5.0, 0.0),
            1,
            Point::new(0.0, 1.0),
        );
        assert_eq!(points, vec![Point::new(2.5, 0.0)]);
    }

    #[test]
    fn fuses_rings_on_unit_bonds() {
        for (input, size) in [
            ("c1ccc2ccccc2c1", 10),
            ("c1cc2ccc3cccc4ccc(c1)c2c34", 16),
            ("CC12CCC3C(CCC4CC(O)CCC34C)C1CCC2O", 17),
        ] {
            let (molecule, local) = system_layout(input);
            assert_eq!(local.iter().flatten().count(), size, "{input}");
            let lengths = bond_lengths(&molecule, &local);
            assert!(
                lengths.iter().all(|&l| close(l, 1.0)),
                "{input}: {lengths:?}"
            );
            assert!(closest_pair(&local) > 0.9, "{input}");
        }
    }

    #[test]
    fn draws_a_bridge_inside_its_ring() {
        let (molecule, local) = system_layout("C1CC2CCC1C2");
        assert!(
            bond_lengths(&molecule, &local)
                .iter()
                .all(|&l| close(l, 1.0))
        );
        assert!(closest_pair(&local) >= CLEARANCE);
    }

    #[test]
    fn relaxes_cages_into_separate_atoms() {
        // A cube drawn flat cannot keep every pair of atoms half a bond apart.
        const LEGIBLE: f64 = 0.35;
        for input in [
            "C12C3C4C1C5C2C3C45",
            "C1C2CC3CC1CC(C2)C3",
            "C1CCC2(CC1)CCCC2",
        ] {
            let (molecule, local) = system_layout(input);
            assert!(
                local
                    .iter()
                    .flatten()
                    .all(|p| p.x.is_finite() && p.y.is_finite())
            );
            assert_eq!(
                local.iter().flatten().count(),
                molecule.atoms().len(),
                "{input}"
            );
            assert!(
                closest_pair(&local) >= LEGIBLE,
                "{input}: {}",
                closest_pair(&local)
            );
            let lengths = bond_lengths(&molecule, &local);
            assert!(
                lengths.iter().all(|&l| (0.6..1.4).contains(&l)),
                "{input}: {lengths:?}"
            );
        }
    }
}
