//! Double bonds outside rings drawn with the configuration their source gave.
//!
//! Placement draws every chain as a zig-zag, which makes most double bonds in a
//! chain trans. Where the molecule says otherwise, one side of the bond is
//! mirrored across it, which keeps every other bond length and angle. Clash
//! removal afterwards only flips and turns branches rigidly, so it keeps the
//! configuration.

use moiety_core::Molecule;

use super::geometry::Point;
use crate::graph::{self, Adjacency};
use crate::rings::RingInfo;

/// Mirrors one side of each stereo double bond outside a ring whose drawn
/// configuration differs from the molecule's.
pub(super) fn settle(
    positions: &mut [Point],
    molecule: &Molecule,
    adjacency: &Adjacency,
    rings: &RingInfo,
) {
    for stereo in molecule.double_bond_stereo() {
        if rings.is_ring_bond(stereo.bond) {
            continue;
        }
        let bond = molecule.bonds()[stereo.bond];
        let (a, b) = (positions[bond.a], positions[bond.b]);
        let side = |n: usize| (b - a).cross(positions[n] - a);
        let (first, second) = (side(stereo.neighbours.0), side(stereo.neighbours.1));
        if first.abs() < 1e-9 || second.abs() < 1e-9 || (first * second > 0.0) == stereo.cis {
            continue;
        }
        let Some(direction) = (b - a).normalized() else {
            continue;
        };
        for atom in graph::side(adjacency, bond.b, bond.a) {
            if atom != bond.b {
                positions[atom] = positions[atom].reflected(a, direction);
            }
        }
    }
}
