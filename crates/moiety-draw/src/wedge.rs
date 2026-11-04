//! Wedge and hash bonds that show the configuration of tetrahedral stereocentres.
//!
//! Each stereocentre gets one single bond, narrow at the centre: preferably to
//! a heteroatom, then to an atom with few neighbours, outside every ring and
//! not itself a stereocentre. The bond is a wedge when lifting its far end
//! toward the viewer, with the hidden hydrogen or lone pair completing the
//! tetrahedron, turns the neighbours the way `@` or `@@` says, and a hash when
//! pushing it away does. The turn is the signed volume of all four neighbours,
//! which does not depend on which of them comes first. A centre with four
//! drawn neighbours is never wedged on a bond whose three others lie within a
//! half-plane, where the drawing would be ambiguous.
//!
//! A centre where fused rings meet, with three ring bonds and a hydrogen,
//! shows the hydrogen instead, as steroids are drawn: on a short wedge or hash
//! out of the rings, or, where that is crowded, into one of them, upright where
//! it can be. Bridged rings keep their wedged ring bond, since they may be
//! drawn in perspective.

use moiety_core::{BondOrder, Element, Molecule, Tetrahedral};

use crate::graph::{Adjacency, bond_between};
use crate::layout::{EPSILON, Point, largest_gap};
use crate::rings::RingInfo;

/// Length of the bond to a hydrogen drawn at a stereocentre, in bond lengths.
pub(crate) const HYDROGEN_BOND: f64 = 0.7;

/// Distance a drawn hydrogen keeps from every other atom where it can, in bond lengths.
const HYDROGEN_CLEARANCE: f64 = 0.5;

/// Distance drawn hydrogens keep from each other where they can, in bond
/// lengths, wider than [`HYDROGEN_CLEARANCE`] since both carry labels.
const HYDROGEN_SPACING: f64 = 0.9;

/// How a stereo bond is drawn.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Wedge {
    /// A filled triangle: the far end is toward the viewer.
    Solid,
    /// A hashed triangle: the far end is away from the viewer.
    Hashed,
}

/// A wedged bond: the atom at its narrow end and how it is drawn.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Wedged {
    pub narrow: usize,
    pub wedge: Wedge,
}

/// A hydrogen drawn at a stereocentre, on a bond wedged from the centre.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Hydrogen {
    pub centre: usize,
    pub position: Point,
    pub wedge: Wedge,
}

/// The stereo marks of a drawing.
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct Marks {
    /// The wedge, if any, of every bond.
    pub bonds: Vec<Option<Wedged>>,
    /// Hydrogens to draw at stereocentres where fused rings meet.
    pub hydrogens: Vec<Hydrogen>,
}

/// The stereo marks of `molecule` drawn at `positions`, one unit per bond
/// length with `y` downward.
pub(crate) fn wedges(molecule: &Molecule, positions: &[Point]) -> Marks {
    let mut marks = Marks {
        bonds: vec![None; molecule.bonds().len()],
        hydrogens: Vec::new(),
    };
    if molecule.tetrahedral().is_empty() {
        return marks;
    }
    let adjacency = molecule.adjacency();
    let rings = RingInfo::new(molecule, &adjacency);
    let centres: Vec<usize> = molecule.tetrahedral().iter().map(|t| t.centre).collect();
    marks.hydrogens = fusion_hydrogens(molecule, positions, &adjacency, &rings);
    for stereo in molecule.tetrahedral() {
        if marks.hydrogens.iter().any(|h| h.centre == stereo.centre) {
            continue;
        }
        let mut candidates: Vec<(usize, usize)> = stereo
            .neighbours
            .iter()
            .flatten()
            .filter_map(|&n| Some((n, bond_between(&adjacency, stereo.centre, n)?)))
            .filter(|&(_, bond)| {
                molecule.bonds()[bond].order == BondOrder::Single && marks.bonds[bond].is_none()
            })
            .collect();
        candidates.sort_by_key(|&(n, bond)| {
            let carbon = molecule.atoms()[n].element() == Some(Element::C);
            (
                centres.contains(&n),
                rings.is_ring_bond(bond),
                carbon,
                adjacency[n].len(),
                n,
            )
        });
        let chosen = candidates.into_iter().find_map(|(n, bond)| {
            let wedge = orientation(stereo, positions, n, &adjacency)?;
            Some((bond, wedge))
        });
        if let Some((bond, wedge)) = chosen {
            marks.bonds[bond] = Some(Wedged {
                narrow: stereo.centre,
                wedge,
            });
        }
    }
    marks
}

/// A stereocentre where fused rings meet, and the directions its hydrogen
/// may take, best first: the angles no ring fills, widest first, then those
/// a ring fills, nearest to vertical first, each clear of the other atoms.
struct Site<'a> {
    stereo: &'a Tetrahedral,
    options: Vec<Point>,
    /// The widest angle at the centre that no ring fills.
    open: f64,
}

impl<'a> Site<'a> {
    /// The site of `stereo`, if its centre has three ring bonds, a hydrogen,
    /// and a ring of its own filling one of its angles.
    fn new(
        stereo: &'a Tetrahedral,
        molecule: &Molecule,
        positions: &[Point],
        adjacency: &Adjacency,
        rings: &RingInfo,
    ) -> Option<Self> {
        let centre = stereo.centre;
        let links = &adjacency[centre];
        let fused = molecule.atoms()[centre].hydrogens > 0
            && links.len() == 3
            && links.iter().all(|&(_, bond)| rings.is_ring_bond(bond))
            && stereo.neighbours.iter().filter(|n| n.is_none()).count() == 1;
        if !fused {
            return None;
        }
        let origin = positions[centre];
        let gaps = gaps(centre, links, positions, rings);
        if !gaps.iter().any(|gap| gap.0) {
            return None;
        }
        let clear = |d: Point| {
            positions.iter().enumerate().all(|(i, p)| {
                i == centre || p.distance(origin + d * HYDROGEN_BOND) >= HYDROGEN_CLEARANCE
            })
        };
        let mut options: Vec<Point> = gaps.iter().map(|gap| gap.1).filter(|&d| clear(d)).collect();
        if options.is_empty() {
            options.push(gaps[0].1);
        }
        let open = gaps
            .iter()
            .filter(|gap| !gap.0)
            .map(|gap| gap.2)
            .fold(0.0, f64::max);
        Some(Self {
            stereo,
            options,
            open,
        })
    }

    fn at(&self, positions: &[Point], option: usize) -> Point {
        positions[self.stereo.centre] + self.options[option] * HYDROGEN_BOND
    }

    /// The hydrogen drawn along option `option`, if its wedge can be told.
    fn hydrogen(&self, positions: &[Point], option: usize) -> Option<Hydrogen> {
        let origin = positions[self.stereo.centre];
        let direction = self.options[option];
        let wedge = sense(self.stereo, |slot, height| match slot {
            Some(n) => {
                let d = positions[n] - origin;
                normalized([d.x, -d.y, 0.0])
            }
            None => normalized([direction.x, -direction.y, height]),
        })?;
        Some(Hydrogen {
            centre: self.stereo.centre,
            position: self.at(positions, option),
            wedge,
        })
    }
}

/// The angles between the three neighbours of `centre`, as whether a ring
/// fills it, the direction through its middle and its width: those no ring
/// fills first, widest first, then the rest, nearest to vertical first.
fn gaps(
    centre: usize,
    links: &[(usize, usize)],
    positions: &[Point],
    rings: &RingInfo,
) -> Vec<(bool, Point, f64)> {
    let origin = positions[centre];
    let mut around: Vec<(usize, f64)> = links
        .iter()
        .map(|&(n, _)| (n, (positions[n] - origin).angle()))
        .collect();
    around.sort_by(|a, b| a.1.total_cmp(&b.1));
    let count = around.len();
    let mut gaps: Vec<(bool, Point, f64)> = (0..count)
        .map(|i| {
            let ((a, from), (b, to)) = (around[i], around[(i + 1) % count]);
            let width = (to - from).rem_euclid(std::f64::consts::TAU);
            let middle = Point::from_angle(from + width / 2.0);
            (rings.fills(a, centre, b), middle, width)
        })
        .collect();
    let rank = |&(filled, middle, width): &(bool, Point, f64)| {
        if filled { -middle.y.abs() } else { -width }
    };
    gaps.sort_by(|x, y| x.0.cmp(&y.0).then(rank(x).total_cmp(&rank(y))));
    gaps
}

/// The hydrogens to draw where fused rings meet, each at its site's first
/// option, except that of two hydrogens closer than [`HYDROGEN_SPACING`], the
/// one with the narrower open angle moves to its next option clear of the rest.
fn fusion_hydrogens(
    molecule: &Molecule,
    positions: &[Point],
    adjacency: &Adjacency,
    rings: &RingInfo,
) -> Vec<Hydrogen> {
    let sites: Vec<Site> = molecule
        .tetrahedral()
        .iter()
        .filter_map(|stereo| Site::new(stereo, molecule, positions, adjacency, rings))
        .collect();
    let mut chosen = vec![0; sites.len()];
    for i in 0..sites.len() {
        for j in i + 1..sites.len() {
            let apart = sites[i]
                .at(positions, chosen[i])
                .distance(sites[j].at(positions, chosen[j]));
            if apart >= HYDROGEN_SPACING {
                continue;
            }
            let mover = if sites[i].open < sites[j].open { i } else { j };
            let clear = |option: usize| {
                let at = sites[mover].at(positions, option);
                (0..sites.len()).all(|other| {
                    other == mover
                        || at.distance(sites[other].at(positions, chosen[other]))
                            >= HYDROGEN_SPACING
                })
            };
            if let Some(option) =
                (0..sites[mover].options.len()).find(|&o| o != chosen[mover] && clear(o))
            {
                chosen[mover] = option;
            }
        }
    }
    sites
        .iter()
        .zip(chosen)
        .filter_map(|(site, option)| site.hydrogen(positions, option))
        .collect()
}

/// Whether the bond from the centre to `lifted` shows `stereo` as a wedge or a
/// hash, or `None` when the drawing is too flat to tell.
fn orientation(
    stereo: &Tetrahedral,
    positions: &[Point],
    lifted: usize,
    adjacency: &Adjacency,
) -> Option<Wedge> {
    let centre = positions[stereo.centre];
    if adjacency[stereo.centre].len() == 4 {
        let others = adjacency[stereo.centre]
            .iter()
            .filter(|&&(n, _)| n != lifted)
            .map(|&(n, _)| (positions[n] - centre).angle());
        if largest_gap(others).1 >= std::f64::consts::PI - EPSILON {
            return None;
        }
    }
    let vector = |n: usize, height: f64| {
        let d = positions[n] - centre;
        let z = if n == lifted { height } else { 0.0 };
        normalized([d.x, -d.y, z])
    };
    sense(stereo, |slot, height| match slot {
        Some(n) => vector(n, height),
        None => {
            let explicit: Vec<[f64; 3]> = adjacency[stereo.centre]
                .iter()
                .map(|&(n, _)| vector(n, height))
                .collect::<Option<_>>()?;
            normalized(explicit.iter().fold([0.0; 3], |sum, v| {
                [sum[0] - v[0], sum[1] - v[1], sum[2] - v[2]]
            }))
        }
    })
}

/// A wedge when raising one neighbour toward the viewer turns the neighbours
/// of `stereo` the way it says, a hash when lowering it does, or `None` when
/// neither does clearly. `vector` gives the unit vector to each neighbour, or
/// to the hidden hydrogen or lone pair for `None`, with the raised neighbour at
/// the given height, in coordinates with `y` up and `z` toward the viewer.
fn sense(
    stereo: &Tetrahedral,
    vector: impl Fn(Option<usize>, f64) -> Option<[f64; 3]>,
) -> Option<Wedge> {
    let volume = |height: f64| -> Option<f64> {
        let ordered: Vec<[f64; 3]> = stereo
            .neighbours
            .iter()
            .map(|&n| vector(n, height))
            .collect::<Option<_>>()?;
        let from = |v: [f64; 3]| {
            [
                v[0] - ordered[0][0],
                v[1] - ordered[0][1],
                v[2] - ordered[0][2],
            ]
        };
        Some(determinant(
            from(ordered[1]),
            from(ordered[2]),
            from(ordered[3]),
        ))
    };
    // Seen from the first neighbour, the others turn anticlockwise when the
    // volume they span from it is negative.
    let wanted = |volume: f64| volume.abs() > 1e-6 && (volume < 0.0) == stereo.anticlockwise;
    if volume(1.0).is_some_and(wanted) {
        Some(Wedge::Solid)
    } else if volume(-1.0).is_some_and(wanted) {
        Some(Wedge::Hashed)
    } else {
        None
    }
}

fn normalized(v: [f64; 3]) -> Option<[f64; 3]> {
    let length = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    (length > 1e-9).then(|| [v[0] / length, v[1] / length, v[2] / length])
}

/// The signed volume of `a`, `b` and `c`, in coordinates with `y` up and `z`
/// toward the viewer.
fn determinant(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
}

#[cfg(test)]
mod tests {
    use moiety_core::smiles;

    use super::*;

    /// (R)-butan-2-ol drawn as a zig-zag with the ethyl up-left, the methyl
    /// up-right and the hydroxyl down: a solid wedge to the hydroxyl, worked
    /// out by hand from the CIP rules.
    #[test]
    fn wedges_the_hydroxyl_of_r_butanol() {
        let positions = [
            Point::new(0.0, 0.0),
            Point::new(0.866, -0.5),
            Point::new(1.732, 0.0),
            Point::new(2.598, -0.5),
            Point::new(1.732, 1.0),
        ];
        let oxygen_bond = |molecule: &Molecule| {
            molecule
                .bonds()
                .iter()
                .position(|b| (b.a, b.b) == (2, 4))
                .unwrap()
        };
        for (input, wedge) in [
            ("CC[C@@H](C)O", Wedge::Solid),
            ("CC[C@H](C)O", Wedge::Hashed),
        ] {
            let molecule = smiles::parse(input).unwrap();
            let wedged = wedges(&molecule, &positions).bonds;
            assert_eq!(
                wedged[oxygen_bond(&molecule)],
                Some(Wedged { narrow: 2, wedge }),
                "{input}"
            );
            assert_eq!(wedged.iter().flatten().count(), 1);
        }
    }

    /// The same centre with its neighbours listed in another order of the same
    /// handedness gets the same wedge.
    #[test]
    fn does_not_depend_on_the_order_of_neighbours() {
        let parsed = smiles::parse("O[C@]12CCCCC(O1)CCCCCCCC2").unwrap();
        let stereo = parsed.tetrahedral()[0].clone();
        assert_eq!(stereo.neighbours, [Some(0), Some(7), Some(15), Some(2)]);
        let mut reordered = Molecule::from_parts(parsed.atoms().to_vec(), parsed.bonds().to_vec());
        reordered.add_tetrahedral(Tetrahedral {
            neighbours: [Some(0), Some(15), Some(2), Some(7)],
            ..stereo
        });
        let drawn = crate::layout::layout(&parsed);
        let first = wedges(&parsed, drawn.positions());
        assert!(first.bonds.iter().flatten().count() == 1);
        assert_eq!(first, wedges(&reordered, drawn.positions()));
    }

    #[test]
    fn draws_the_hydrogen_at_a_ring_fusion() {
        for (input, alike) in [
            ("C1CC[C@H]2CCCC[C@@H]2C1", false),
            ("C1CC[C@H]2CCCC[C@H]2C1", true),
        ] {
            let molecule = smiles::parse(input).unwrap();
            let drawn = crate::layout::layout(&molecule);
            let p = drawn.positions();
            let marks = wedges(&molecule, p);
            assert!(marks.bonds.iter().all(Option::is_none), "{input}");
            let kinds: Vec<Wedge> = marks.hydrogens.iter().map(|h| h.wedge).collect();
            assert_eq!(kinds.len(), 2, "{input}");
            assert_eq!(kinds[0] == kinds[1], alike, "{input}");
            for hydrogen in &marks.hydrogens {
                let length = hydrogen.position.distance(p[hydrogen.centre]);
                assert!((length - HYDROGEN_BOND).abs() < 1e-9);
                assert!(
                    p.iter()
                        .enumerate()
                        .all(|(i, q)| i == hydrogen.centre || q.distance(hydrogen.position) >= 0.5)
                );
            }
        }
        let camphor = smiles::parse("CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2").unwrap();
        let drawn = crate::layout::layout(&camphor);
        assert!(wedges(&camphor, drawn.positions()).hydrogens.is_empty());
    }

    /// Testosterone's hydrogens at C9, C8 and C14 stand upright, as in the usual
    /// drawing: C8's up into ring C, the others down, none crowding another.
    #[test]
    fn stands_steroid_hydrogens_upright_and_apart() {
        let molecule =
            smiles::parse("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C").unwrap();
        let drawn = crate::layout::layout(&molecule);
        let p = drawn.positions();
        let marks = wedges(&molecule, p);
        let rise = |centre: usize| {
            let hydrogen = marks.hydrogens.iter().find(|h| h.centre == centre).unwrap();
            (p[centre].y - hydrogen.position.y) / HYDROGEN_BOND
        };
        assert!(rise(5) > 0.95, "C8 points up");
        assert!(rise(4) < -0.95 && rise(6) < -0.9, "C9 and C14 point down");
        for (i, a) in marks.hydrogens.iter().enumerate() {
            for b in &marks.hydrogens[i + 1..] {
                assert!(a.position.distance(b.position) >= HYDROGEN_SPACING);
            }
        }
    }

    #[test]
    fn gives_each_centre_its_own_bond() {
        let molecule = smiles::parse("C[C@H](O)[C@@H](O)C").unwrap();
        let positions = crate::layout::layout(&molecule);
        let wedged = wedges(&molecule, positions.positions()).bonds;
        let narrow: Vec<usize> = wedged.iter().flatten().map(|w| w.narrow).collect();
        assert_eq!(narrow.len(), 2);
        assert!(narrow.contains(&1) && narrow.contains(&3));
        assert!(
            wedges(&smiles::parse("CC(O)CC").unwrap(), &[Point::default(); 5])
                .bonds
                .iter()
                .all(Option::is_none)
        );
    }
}
