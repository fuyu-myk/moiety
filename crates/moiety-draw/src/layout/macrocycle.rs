//! Large rings drawn the way chemists draw macrocycles: as a closed zig-zag on
//! a hexagonal lattice, a large hexagon whose sides are zig-zags, with the atoms
//! that carry other rings, then substituents, on its outward corners.
//!
//! Walking such an outline counterclockwise, every atom turns left or right by
//! 60 degrees. A side of the large hexagon is a left turn followed by `p` pairs
//! of right and left turns, so an outline of `m` atoms has `(m + 6) / 2` outward
//! corners. It closes when the zig-zag counts of opposite sides balance. A ring
//! with an odd number of atoms is drawn on the outline one larger, without one
//! corner, and its bonds are evened out.

use std::f64::consts::FRAC_PI_3;

use super::geometry::Point;

/// The smallest ring drawn as a lattice outline; smaller rings are regular polygons.
pub(super) const MIN_SIZE: usize = 12;

/// The largest ring drawn as a lattice outline, which bounds the outlines tried.
const MAX_SIZE: usize = 46;

/// Outlines tried for one ring, the roundest first.
const SHAPES: usize = 8;

/// Passes that even out the bonds of an odd ring.
const EVENING: usize = 200;

/// Positions for the atoms of a ring of `inward.len()` atoms, in ring order, on
/// the outline whose inward corners cost least, `inward[i]` being the cost of
/// atom `i` on one; or `None` for a ring outside [`MIN_SIZE`] to [`MAX_SIZE`].
pub(super) fn outline(inward: &[u32]) -> Option<Vec<Point>> {
    let n = inward.len();
    if !(MIN_SIZE..=MAX_SIZE).contains(&n) {
        return None;
    }
    let m = n + n % 2;
    let mut best: Option<Choice> = None;
    for turns in shapes(m) {
        let removals: Vec<Option<usize>> = if n == m {
            vec![None]
        } else {
            corners(&turns).into_iter().map(Some).collect()
        };
        for removed in removals {
            let kept: Vec<i8> = (0..m)
                .filter(|&v| Some(v) != removed)
                .map(|v| turns[v])
                .collect();
            for start in 0..n {
                for forward in [true, false] {
                    let cost = (0..n)
                        .filter(|&i| kept[i] < 0)
                        .map(|i| inward[atom_at(i, start, forward, n)])
                        .sum::<u32>();
                    if best.as_ref().is_none_or(|b| cost < b.cost) {
                        best = Some(Choice {
                            cost,
                            turns: turns.clone(),
                            removed,
                            start,
                            forward,
                        });
                    }
                }
            }
        }
    }
    let Choice {
        turns,
        removed,
        start,
        forward,
        ..
    } = best?;
    let lattice = walk(&turns);
    let vertices: Vec<Point> = (0..m)
        .filter(|&v| Some(v) != removed)
        .map(|v| lattice[v])
        .collect();
    let mut points = vec![Point::default(); n];
    for (i, &point) in vertices.iter().enumerate() {
        points[atom_at(i, start, forward, n)] = point;
    }
    if removed.is_some() {
        even_out(&mut points);
    }
    Some(points)
}

/// An outline and how the ring sits on it: the vertex left out for an odd ring,
/// and the ring position on the first vertex and the direction it runs.
struct Choice {
    cost: u32,
    turns: Vec<i8>,
    removed: Option<usize>,
    start: usize,
    forward: bool,
}

/// The ring position of the atom on outline vertex `i`.
fn atom_at(i: usize, start: usize, forward: bool, n: usize) -> usize {
    if forward {
        (start + i) % n
    } else {
        (start + n - i % n) % n
    }
}

/// Turn sequences of closed outlines of `m` atoms, the roundest first.
fn shapes(m: usize) -> Vec<Vec<i8>> {
    let total = (m - 6) / 2;
    let mut found: Vec<[usize; 6]> = Vec::new();
    for a in 0..=total {
        for b in 0..=total - a {
            for c in 0..=total - a - b {
                for d in 0..=total - a - b - c {
                    let e_f = total - a - b - c - d;
                    for e in 0..=e_f {
                        let f = e_f - e;
                        if a + b == d + e && b + c == e + f {
                            found.push([a, b, c, d, e, f]);
                        }
                    }
                }
            }
        }
    }
    let spread = |p: &[usize; 6]| p.iter().max().unwrap_or(&0) - p.iter().min().unwrap_or(&0);
    found.sort_by_key(|p| (spread(p), std::cmp::Reverse(p[0] + p[3])));
    found
        .into_iter()
        .take(SHAPES)
        .map(|p| {
            p.iter()
                .flat_map(|&zigzags| {
                    std::iter::once(1).chain(std::iter::repeat_n([-1, 1], zigzags).flatten())
                })
                .collect()
        })
        .collect()
}

/// The vertices where a side of the large hexagon starts.
fn corners(turns: &[i8]) -> Vec<usize> {
    (0..turns.len())
        .filter(|&v| turns[v] > 0 && turns[(v + turns.len() - 1) % turns.len()] > 0)
        .collect()
}

/// Lattice positions of an outline walked counterclockwise from the origin,
/// where `turns[v]` is the turn at vertex `v`, between its two bonds.
fn walk(turns: &[i8]) -> Vec<Point> {
    let m = turns.len();
    let mut points = Vec::with_capacity(m);
    let (mut at, mut heading) = (Point::default(), 0i32);
    for v in 0..m {
        points.push(at);
        at = at + Point::from_angle(-f64::from(heading) * FRAC_PI_3);
        heading += i32::from(turns[(v + 1) % m]);
    }
    points
}

/// Nudges the atoms of a closed ring until every bond is close to unit length.
fn even_out(points: &mut [Point]) {
    let n = points.len();
    for _ in 0..EVENING {
        for i in 0..n {
            let j = (i + 1) % n;
            let delta = points[j] - points[i];
            let length = delta.length();
            if length > 1e-9 {
                let correction = delta * ((length - 1.0) / (2.0 * length));
                points[i] = points[i] + correction;
                points[j] = points[j] - correction;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bonds(points: &[Point]) -> Vec<f64> {
        (0..points.len())
            .map(|i| points[i].distance(points[(i + 1) % points.len()]))
            .collect()
    }

    #[test]
    fn closes_even_and_odd_rings_with_unit_bonds() {
        for n in [12, 14, 18, 23, 31, 46] {
            let points = outline(&vec![0; n]).unwrap();
            assert_eq!(points.len(), n);
            assert!(
                bonds(&points).iter().all(|&b| (b - 1.0).abs() < 0.05),
                "{n}"
            );
            for i in 0..n {
                for j in i + 1..n {
                    assert!(points[i].distance(points[j]) > 0.8, "{n}: {i} and {j}");
                }
            }
        }
        assert!(outline(&[0; 11]).is_none());
        assert!(outline(&[0; 47]).is_none());
    }

    #[test]
    fn puts_marked_atoms_on_outward_corners() {
        let mut inward = vec![0; 14];
        for i in [0, 3, 5, 9, 12] {
            inward[i] = 10;
        }
        let points = outline(&inward).unwrap();
        let centre = points.iter().fold(Point::default(), |s, &p| s + p) * (1.0 / 14.0);
        for i in [0, 3, 5, 9, 12] {
            let (before, after) = (points[(i + 13) % 14], points[(i + 1) % 14]);
            let mid = (before + after) * 0.5;
            assert!(
                points[i].distance(centre) > mid.distance(centre),
                "atom {i} points outward"
            );
        }
    }
}
