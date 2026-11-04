//! Wedge and hash bonds that show the configuration of tetrahedral stereocentres.
//!
//! Each stereocentre gets one single bond, narrow at the centre: preferably to
//! a heteroatom, then to an atom with few neighbours, outside every ring and
//! not itself a stereocentre. The bond is a wedge when lifting its far end
//! toward the viewer, with the hidden hydrogen or lone pair completing the
//! tetrahedron, turns the neighbours the way `@` or `@@` says, and a hash when
//! pushing it away does.

use moiety_core::{BondOrder, Element, Molecule, Tetrahedral};

use crate::graph::{Adjacency, bond_between};
use crate::layout::Point;
use crate::rings::RingInfo;

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

/// The wedge, if any, of every bond of `molecule` drawn at `positions`, one
/// unit per bond length with `y` downward.
pub(crate) fn wedges(molecule: &Molecule, positions: &[Point]) -> Vec<Option<Wedged>> {
    let mut wedged = vec![None; molecule.bonds().len()];
    if molecule.tetrahedral().is_empty() {
        return wedged;
    }
    let adjacency = molecule.adjacency();
    let rings = RingInfo::new(molecule, &adjacency);
    let centres: Vec<usize> = molecule.tetrahedral().iter().map(|t| t.centre).collect();
    for stereo in molecule.tetrahedral() {
        let mut candidates: Vec<(usize, usize)> = stereo
            .neighbours
            .iter()
            .flatten()
            .filter_map(|&n| Some((n, bond_between(&adjacency, stereo.centre, n)?)))
            .filter(|&(_, bond)| {
                molecule.bonds()[bond].order == BondOrder::Single && wedged[bond].is_none()
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
            wedged[bond] = Some(Wedged {
                narrow: stereo.centre,
                wedge,
            });
        }
    }
    wedged
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
    let turn = |height: f64| -> Option<f64> {
        let vector = |n: usize| {
            let d = positions[n] - centre;
            let z = if n == lifted { height } else { 0.0 };
            normalized([d.x, -d.y, z])
        };
        let explicit: Vec<[f64; 3]> = adjacency[stereo.centre]
            .iter()
            .map(|&(n, _)| vector(n))
            .collect::<Option<_>>()?;
        let hidden = normalized(explicit.iter().fold([0.0; 3], |sum, v| {
            [sum[0] - v[0], sum[1] - v[1], sum[2] - v[2]]
        }));
        let ordered: Vec<[f64; 3]> = stereo
            .neighbours
            .iter()
            .map(|n| match n {
                Some(n) => vector(*n),
                None => hidden,
            })
            .collect::<Option<_>>()?;
        Some(determinant(ordered[0], ordered[1], ordered[2]))
    };
    let wanted = |volume: f64| volume.abs() > 1e-6 && (volume > 0.0) == stereo.anticlockwise;
    if turn(1.0).is_some_and(wanted) {
        Some(Wedge::Solid)
    } else if turn(-1.0).is_some_and(wanted) {
        Some(Wedge::Hashed)
    } else {
        None
    }
}

fn normalized(v: [f64; 3]) -> Option<[f64; 3]> {
    let length = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    (length > 1e-9).then(|| [v[0] / length, v[1] / length, v[2] / length])
}

/// The signed volume of `a`, `b` and `c`, positive when they turn anticlockwise
/// seen from the side `a` points to, in coordinates with `y` up and `z` toward the viewer.
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
            let wedged = wedges(&molecule, &positions);
            assert_eq!(
                wedged[oxygen_bond(&molecule)],
                Some(Wedged { narrow: 2, wedge }),
                "{input}"
            );
            assert_eq!(wedged.iter().flatten().count(), 1);
        }
    }

    #[test]
    fn gives_each_centre_its_own_bond() {
        let molecule = smiles::parse("C[C@H](O)[C@@H](O)C").unwrap();
        let positions = crate::layout::layout(&molecule);
        let wedged = wedges(&molecule, positions.positions());
        let narrow: Vec<usize> = wedged.iter().flatten().map(|w| w.narrow).collect();
        assert_eq!(narrow.len(), 2);
        assert!(narrow.contains(&1) && narrow.contains(&3));
        assert!(
            wedges(&smiles::parse("CC(O)CC").unwrap(), &[Point::default(); 5])
                .iter()
                .all(Option::is_none)
        );
    }
}
