//! Fixed drawings for ring systems that ring templates and bridges cannot draw
//! the way chemists do, with coordinates from PubChem's two-dimensional
//! depictions, which are in the public domain: porphine, the morphine skeleton
//! and norbornane in perspective.
//!
//! A drawing applies to any ring system with the same skeleton, whatever its
//! elements and bond orders, so porphine's also draws chlorins and the
//! norbornane drawing camphor. Of the ways a skeleton maps onto the system, the
//! one where substituents fit best is kept.

use super::geometry::Point;
use super::quality::Frame;

/// Most mappings of one skeleton onto a ring system that are compared.
const MAX_MAPPINGS: usize = 64;

/// A ring skeleton and its drawing, one unit per mean bond length, `y` downward.
struct Template {
    points: &'static [(f64, f64)],
    bonds: &'static [(usize, usize)],
}

const TEMPLATES: [Template; 3] = [
    Template {
        points: &[
            (1.5061, 0.0139),
            (0.0044, 1.3834),
            (-0.0291, -1.3556),
            (-1.4991, 0.0139),
            (2.1072, 0.7489),
            (0.8063, 1.9845),
            (2.0737, -0.7880),
            (-0.7976, 2.0180),
            (1.8075, 1.7166),
            (3.0431, 0.4810),
            (0.5050, 2.9522),
            (3.0431, -0.5201),
            (-0.4962, 2.9522),
            (1.7740, -1.7222),
            (-1.7653, 1.7166),
            (0.8063, -2.0236),
            (-2.0666, 0.7824),
            (-0.8310, -1.9901),
            (-2.0666, -0.7545),
            (-1.7987, -1.7222),
            (0.4715, -2.9578),
            (-3.0343, 0.5145),
            (-0.5297, -2.9578),
            (-3.0343, -0.4866),
        ],
        bonds: &[
            (0, 4),
            (0, 6),
            (1, 5),
            (1, 7),
            (2, 15),
            (2, 17),
            (3, 16),
            (3, 18),
            (4, 8),
            (4, 9),
            (5, 8),
            (5, 10),
            (6, 11),
            (6, 13),
            (7, 12),
            (7, 14),
            (9, 11),
            (10, 12),
            (13, 15),
            (14, 16),
            (15, 20),
            (16, 21),
            (17, 19),
            (17, 22),
            (18, 19),
            (18, 23),
            (20, 22),
            (21, 23),
        ],
    },
    Template {
        points: &[
            (-1.8338, -0.0780),
            (2.0955, 0.9258),
            (-0.3755, 0.4502),
            (0.4481, 0.9258),
            (1.2718, 0.4502),
            (-1.1992, 0.9258),
            (0.1740, -0.2471),
            (-0.3755, -0.5007),
            (1.2718, -0.5007),
            (1.5266, -0.2471),
            (0.4481, -0.9763),
            (-1.1992, 1.8769),
            (0.4481, 1.8769),
            (-1.1992, -0.9763),
            (-0.3755, 2.3524),
            (0.4481, -1.9274),
            (-1.1992, -1.9274),
            (-0.3755, -2.4029),
        ],
        bonds: &[
            (0, 5),
            (0, 13),
            (1, 4),
            (1, 9),
            (2, 3),
            (2, 5),
            (2, 6),
            (2, 7),
            (3, 4),
            (3, 12),
            (4, 8),
            (5, 11),
            (6, 9),
            (7, 10),
            (7, 13),
            (8, 10),
            (10, 15),
            (11, 14),
            (12, 14),
            (13, 16),
            (15, 17),
            (16, 17),
        ],
    },
    Template {
        points: &[
            (0.0242, 0.2542),
            (0.2001, -0.4022),
            (-0.3004, -1.0698),
            (-1.1408, 0.0179),
            (-0.9084, 0.6564),
            (1.1506, -0.0564),
            (0.9747, 0.5999),
        ],
        bonds: &[
            (0, 2),
            (0, 4),
            (0, 6),
            (1, 2),
            (1, 3),
            (1, 5),
            (3, 4),
            (5, 6),
        ],
    },
];

/// Positions for the ring system of `frame` from the first skeleton that matches
/// it, mapped so that substituents fit best, if one does. A template's own
/// geometry is trusted, so the bridge of the morphine skeleton may pass close
/// to the atoms behind it, as it does in the usual drawing.
pub(super) fn matched(frame: &Frame, atom_count: usize) -> Option<Vec<Option<Point>>> {
    let atoms = frame.atoms();
    TEMPLATES
        .iter()
        .filter(|t| t.points.len() == atoms.len() && t.bonds.len() == frame.bonds().len())
        .find_map(|template| {
            mappings(template, frame)
                .into_iter()
                .map(|mapping| {
                    let mut local = vec![None; atom_count];
                    for (vertex, &atom) in mapping.iter().enumerate() {
                        let (x, y) = template.points[vertex];
                        local[atom] = Some(Point::new(x, y));
                    }
                    let quality = frame.quality(&local);
                    (quality, local)
                })
                .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(_, local)| local)
        })
}

/// Up to [`MAX_MAPPINGS`] maps from the template's vertices onto the system's
/// atoms that carry every template bond onto a system bond.
fn mappings(template: &Template, frame: &Frame) -> Vec<Vec<usize>> {
    let n = template.points.len();
    let mut neighbours = vec![Vec::new(); n];
    for &(a, b) in template.bonds {
        neighbours[a].push(b);
        neighbours[b].push(a);
    }
    let mut order = vec![0];
    let mut seen = vec![false; n];
    seen[0] = true;
    let mut cursor = 0;
    while let Some(&vertex) = order.get(cursor) {
        cursor += 1;
        for &next in &neighbours[vertex] {
            if !seen[next] {
                seen[next] = true;
                order.push(next);
            }
        }
    }
    if order.len() != n {
        return Vec::new();
    }
    let mut found = Vec::new();
    let mut map = vec![usize::MAX; n];
    extend(&order, 0, &neighbours, frame, &mut map, &mut found);
    found
}

fn extend(
    order: &[usize],
    depth: usize,
    neighbours: &[Vec<usize>],
    frame: &Frame,
    map: &mut [usize],
    found: &mut Vec<Vec<usize>>,
) {
    if found.len() >= MAX_MAPPINGS {
        return;
    }
    let Some(&vertex) = order.get(depth) else {
        found.push(map.to_vec());
        return;
    };
    let candidates: Vec<usize> = match neighbours[vertex].iter().find(|&&v| map[v] != usize::MAX) {
        Some(&placed) => frame.neighbours(map[placed]).to_vec(),
        None => frame.atoms().to_vec(),
    };
    for atom in candidates {
        let fits = !map.contains(&atom)
            && frame.neighbours(atom).len() == neighbours[vertex].len()
            && neighbours[vertex]
                .iter()
                .filter(|&&v| map[v] != usize::MAX)
                .all(|&v| frame.neighbours(atom).contains(&map[v]));
        if fits {
            map[vertex] = atom;
            extend(order, depth + 1, neighbours, frame, map, found);
            map[vertex] = usize::MAX;
        }
    }
}
