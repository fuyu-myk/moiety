//! Fixed drawings for ring systems that ring templates and bridges cannot draw
//! the way chemists do, with coordinates from PubChem's two-dimensional
//! depictions, which are in the public domain: porphine, the morphine skeleton,
//! norbornane in perspective, the taxane core of paclitaxel and artemisinin.
//! The steroid nucleus is drawn in PubChem's orientation, ring A at the lower
//! left and ring D at the upper right, with its rings made regular.
//!
//! A drawing applies to any ring system with the same skeleton, whatever its
//! elements and bond orders, so porphine's also draws chlorins, the norbornane
//! drawing camphor and the steroid drawing estradiol. Of the ways a skeleton
//! maps onto the system, the one where substituents fit best is kept.

use super::geometry::Point;
use super::quality::Frame;

/// Most mappings of one skeleton onto a ring system that are compared.
const MAX_MAPPINGS: usize = 64;

/// A ring skeleton and its drawing, one unit per mean bond length, `y` downward.
struct Template {
    points: &'static [(f64, f64)],
    bonds: &'static [(usize, usize)],
}

const TEMPLATES: [Template; 6] = [
    // Porphine.
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
    // The morphine skeleton.
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
    // Norbornane, in perspective.
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
    // Gonane, the steroid nucleus.
    Template {
        points: &[
            (0.612673627, 0.205882353),
            (-0.253351777, -0.294117647),
            (1.478699031, -0.294117647),
            (-1.119377181, 0.205882353),
            (1.478699031, -1.294117647),
            (-1.119377181, 1.205882353),
            (-0.253351777, -1.294117647),
            (0.612673627, 1.205882353),
            (0.612673627, -1.794117647),
            (-0.253351777, 1.705882353),
            (2.429755547, 0.014899347),
            (2.429755547, -1.603134641),
            (-1.985402584, -0.294117647),
            (3.017540799, -0.794117647),
            (-1.985402584, 1.705882353),
            (-2.851427988, 0.205882353),
            (-2.851427988, 1.205882353),
        ],
        bonds: &[
            (0, 1),
            (0, 2),
            (0, 7),
            (1, 3),
            (1, 6),
            (2, 4),
            (2, 10),
            (3, 5),
            (3, 12),
            (4, 8),
            (4, 11),
            (5, 9),
            (5, 14),
            (6, 8),
            (7, 9),
            (10, 13),
            (11, 13),
            (12, 15),
            (14, 16),
            (15, 16),
        ],
    },
    // The taxane core of paclitaxel, with its oxetane.
    Template {
        points: &[
            (-0.1013, 3.0478),
            (-0.0746, 0.5908),
            (0.0746, 1.6028),
            (-0.9818, 0.2150),
            (-0.7322, 2.2318),
            (0.8327, 0.2150),
            (1.2086, -0.6922),
            (-1.8029, 0.8251),
            (-1.6772, 1.8404),
            (0.8327, -1.5996),
            (0.7121, 2.4136),
            (-1.3576, -0.6922),
            (1.8696, -1.4206),
            (-0.0746, -1.9754),
            (-0.9818, -1.5996),
            (1.5997, -2.3665),
            (0.6538, -2.6364),
        ],
        bonds: &[
            (0, 4),
            (0, 10),
            (1, 2),
            (1, 3),
            (1, 5),
            (2, 4),
            (2, 10),
            (3, 7),
            (3, 11),
            (4, 8),
            (5, 6),
            (6, 9),
            (6, 12),
            (7, 8),
            (9, 13),
            (11, 14),
            (12, 15),
            (13, 14),
            (13, 16),
            (15, 16),
        ],
    },
    // Artemisinin.
    Template {
        points: &[
            (-0.3526, -0.4803),
            (1.1264, -1.1732),
            (1.8269, 0.2641),
            (-0.1291, -1.4598),
            (-0.0023, 0.2421),
            (-0.9074, -0.1938),
            (0.0921, 1.2843),
            (-1.7810, 0.3823),
            (0.9029, -0.1938),
            (-0.7680, 1.8805),
            (-1.1309, -1.1732),
            (-1.7109, 1.4264),
            (0.9747, 1.8174),
            (0.5000, -1.9587),
            (-0.5046, -1.9587),
            (1.8636, 1.2947),
        ],
        bonds: &[
            (0, 3),
            (4, 0),
            (1, 8),
            (1, 13),
            (2, 8),
            (2, 15),
            (3, 13),
            (4, 5),
            (4, 6),
            (4, 8),
            (5, 7),
            (5, 10),
            (6, 9),
            (6, 12),
            (7, 11),
            (9, 11),
            (10, 14),
            (12, 15),
            (13, 14),
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

#[cfg(test)]
mod tests {
    use moiety_core::smiles;

    use super::*;
    use crate::rings::RingInfo;

    fn largest_system_matches(input: &str) -> bool {
        let molecule = smiles::parse(input).unwrap();
        let adjacency = molecule.adjacency();
        let info = RingInfo::new(&molecule, &adjacency);
        let system = info
            .systems()
            .iter()
            .max_by_key(|system| system.atoms.len())
            .unwrap();
        let frame = Frame::new(&info, system, &molecule, &adjacency);
        matched(&frame, adjacency.len()).is_some()
    }

    #[test]
    fn matches_each_skeleton_whatever_its_elements_and_bonds() {
        for (name, input) in [
            (
                "porphine",
                "C1=CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2",
            ),
            ("morphine", "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O"),
            ("camphor", "CC1(C)C2CCC1(C)C(=O)C2"),
            ("androstane", "CC12CCCCC1CCC1C2CCC2(C)CCCC12"),
            (
                "estradiol",
                "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O",
            ),
            (
                "docetaxel",
                "CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O",
            ),
            (
                "artemether",
                "C[C@@H]1CC[C@H]2[C@H]([C@H](O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)OC)C",
            ),
        ] {
            assert!(largest_system_matches(input), "{name}");
        }
        for input in ["c1ccc2ccccc2c1", "C1CCC2CC3CC4CCCC4CC3CC2C1"] {
            assert!(!largest_system_matches(input), "{input}");
        }
    }
}
