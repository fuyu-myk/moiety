//! Large rings drawn the way chemists draw macrocycles: as a closed zig-zag on
//! a hexagonal lattice, a large hexagon whose sides are zig-zags, with the atoms
//! that carry other rings, then substituents, on its outward corners.
//!
//! Walking such an outline counterclockwise, every atom turns left or right by
//! 60 degrees. A side of the large hexagon is a left turn followed by `p` pairs
//! of right and left turns, so an outline of `m` atoms has `(m + 6) / 2` outward
//! corners. It closes when the zig-zag counts of opposite sides balance. A ring
//! with an odd number of atoms is drawn on the outline one larger without one
//! corner, whose two neighbours are joined by one longer bond, placed where
//! neither carries a substituent if it can be. The outline is turned so that
//! its bonds run at 30 and 90 degrees, as chemists draw zig-zags, and so that
//! it is wider than tall.

use std::f64::consts::{FRAC_PI_3, FRAC_PI_6};

use super::geometry::Point;

/// The smallest ring drawn as a lattice outline; smaller rings are regular polygons.
pub(super) const MIN_SIZE: usize = 12;

/// The largest ring drawn as a lattice outline, which bounds the outlines tried.
const MAX_SIZE: usize = 46;

/// Outlines tried for one ring, the roundest first.
const SHAPES: usize = 24;

/// Cost of each step away from the roundest outline, which breaks near ties
/// toward round outlines without outweighing a substituent.
const ROUNDNESS: u32 = 4;

/// Passes that even out the three bonds around the corner an odd ring leaves out.
const EVENING: usize = 100;

/// Positions for the atoms of a ring of `costs.len()` atoms, in ring order, on
/// the outline where they cost least, `costs[i]` being what atom `i` costs on an
/// inward corner and on an outward one; the atom `top` ends up highest when that
/// keeps the outline as wide; or `None` for a ring
/// outside [`MIN_SIZE`] to [`MAX_SIZE`].
pub(super) fn outline(costs: &[(u32, u32)], top: Option<usize>) -> Option<Vec<Point>> {
    let n = costs.len();
    if !(MIN_SIZE..=MAX_SIZE).contains(&n) {
        return None;
    }
    let m = n + n % 2;
    let mut best: Option<Choice> = None;
    for (rank, turns) in (0u32..).zip(shapes(m)) {
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
                    let at = |i: usize| costs[atom_at(i, start, forward, n)];
                    let defect = removed.map_or(0, |v| {
                        let (before, after) = beside(v, m);
                        let (a, b) = (at(before), at(after));
                        (a.0 + a.1 + b.0 + b.1) / 2
                    });
                    let cost = defect
                        + rank * ROUNDNESS
                        + (0..n)
                            .map(|i| if kept[i] < 0 { at(i).0 } else { at(i).1 })
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
    if let Some(v) = removed {
        let (before, after) = beside(v, m);
        let (a, b) = (
            atom_at(before, start, forward, n),
            atom_at(after, start, forward, n),
        );
        let first = if (a + 1) % n == b { a } else { b };
        even_out(&mut points, first);
    }
    Some(upright(points, top))
}

/// The outline positions, among the vertices kept, on either side of the
/// removed vertex `v` of an outline of `m` vertices.
fn beside(v: usize, m: usize) -> (usize, usize) {
    let kept = |x: usize| if x < v { x } else { x - 1 };
    (kept((v + m - 1) % m), kept((v + 1) % m))
}

/// The outline turned to bonds at 30 and 90 degrees, and then by the multiple of
/// 60 degrees that makes it widest for its height, with the atom `top` highest
/// among the turns that are as wide.
fn upright(points: Vec<Point>, top: Option<usize>) -> Vec<Point> {
    let score = |angle: f64| {
        let turned: Vec<Point> = points.iter().map(|p| p.rotated(angle)).collect();
        let span = |f: fn(&Point) -> f64| {
            let values = turned.iter().map(f);
            values.clone().fold(f64::NEG_INFINITY, f64::max) - values.fold(f64::INFINITY, f64::min)
        };
        let ratio = span(|p| p.x) / span(|p| p.y).max(1e-9);
        // Rounded so that turns as wide as each other tie; y grows downward.
        let height = top.map_or(0.0, |i| -turned[i].y);
        ((ratio * 1e6).round(), height)
    };
    let best = (0..6)
        .map(|k| FRAC_PI_6 + FRAC_PI_3 * f64::from(k))
        .map(|angle| (angle, score(angle)))
        .fold(
            (FRAC_PI_6, (f64::NEG_INFINITY, f64::NEG_INFINITY)),
            |best, next| if next.1 > best.1 { next } else { best },
        )
        .0;
    points.iter().map(|p| p.rotated(best)).collect()
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

/// Moves the ring atoms at `first` and the one after it, joined by the long bond
/// across a removed corner, until it and the bonds on either side of them are
/// as close to unit length as they can be. Every other atom stays on the lattice.
fn even_out(points: &mut [Point], first: usize) {
    let n = points.len();
    let (p, q) = (first, (first + 1) % n);
    let bonds = [((p + n - 1) % n, p), (p, q), (q, (q + 1) % n)];
    let movable = |atom: usize| atom == p || atom == q;
    for _ in 0..EVENING {
        for &(i, j) in &bonds {
            let delta = points[j] - points[i];
            let length = delta.length();
            if length < 1e-9 {
                continue;
            }
            let error = delta * ((length - 1.0) / length);
            match (movable(i), movable(j)) {
                (true, true) => {
                    points[i] = points[i] + error * 0.5;
                    points[j] = points[j] - error * 0.5;
                }
                (true, false) => points[i] = points[i] + error,
                (false, true) => points[j] = points[j] - error,
                (false, false) => {}
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
            let points = outline(&vec![(0, 0); n], None).unwrap();
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
        assert!(outline(&[(0, 0); 11], None).is_none());
        assert!(outline(&[(0, 0); 47], None).is_none());
    }

    #[test]
    fn puts_marked_atoms_on_outward_corners() {
        let mut costs = vec![(0, 0); 14];
        for i in [0, 3, 5, 9, 12] {
            costs[i] = (10, 0);
        }
        let points = outline(&costs, None).unwrap();
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
