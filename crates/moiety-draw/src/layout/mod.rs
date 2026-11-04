//! Two-dimensional coordinates for skeletal structure drawings.
//!
//! [`layout`] is a template layout, not a general embedding engine. It knows
//! the shapes that condensed formulas produce and draws them the way a chemist
//! would:
//!
//! - chains zig-zag at 120 degrees, alternating sides, along the longest path;
//! - branches spread evenly: three neighbours at 120 degrees, four at 90;
//! - an atom outside the rings with three or more alike branches, such as the
//!   phosphorus of `PPh3`, is drawn first, its branches evenly spaced with one
//!   pointing up; when they crowd each other its bonds stretch equally, up to
//!   1.5 bond lengths, rather than bend;
//! - atoms with a triple bond or two double bonds (`O=C=O`, allenes) are straight;
//! - double bonds keep the cis or trans configuration their source gave, in
//!   chains and in macrocycles;
//! - rings are regular polygons with substituents pointing away from the ring
//!   centre, and rings fused on an edge lie on opposite sides of it;
//! - porphyrins, the morphine skeleton and norbornanes take fixed drawings
//!   from PubChem's depictions, so that they look the way they are usually drawn;
//! - a ring of 12 to 46 atoms is a closed zig-zag on a hexagonal lattice, the
//!   way macrocycles are drawn, with the atoms that carry other rings or the
//!   largest substituents on its outward corners;
//! - a bridge across a ring is an arc or a shallow bow inside it, chosen so that
//!   no bonds cross and no atom, or substituent, lands on another; a ring system
//!   that still draws badly is started from another ring, then relaxed, and the
//!   most legible drawing is kept;
//! - a final pass flips or turns branches about single bonds where unbonded
//!   atoms come too close, stretching a crowded hub's bonds only when no flip or
//!   turn parts its branches without bending it;
//! - disconnected components sit side by side, left to right, in input order.
//!
//! The output is deterministic: the same molecule always gets the same coordinates.

mod double_bonds;
mod geometry;
mod hub;
mod macrocycle;
mod overlap;
mod place;
mod quality;
mod relax;
mod ring_system;
mod templates;

use moiety_core::Molecule;

pub use geometry::Point;
pub(crate) use geometry::{Bounds, EPSILON, centroid};

use crate::rings::RingInfo;
use place::Placer;

/// Horizontal space between the atoms of neighbouring components, in bond lengths.
pub const COMPONENT_GAP: f64 = 2.0;

/// Coordinates for every atom of a molecule, one unit per bond length, `y` downward.
///
/// Most bonds are one unit long. A crowded hub's bonds may stretch, a bridge
/// across a ring may be shorter, and a template keeps its source's lengths, all
/// between 0.55 and 1.6.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Layout {
    positions: Vec<Point>,
    rings: Vec<Vec<usize>>,
}

impl Layout {
    /// The position of every atom, by atom index.
    #[must_use]
    pub fn positions(&self) -> &[Point] {
        &self.positions
    }

    /// The position of atom `index`, if the molecule has that atom.
    #[must_use]
    pub fn position(&self, index: usize) -> Option<Point> {
        self.positions.get(index).copied()
    }

    /// The smallest rings of the molecule, each listed once as atom indices in ring order.
    #[must_use]
    pub fn rings(&self) -> &[Vec<usize>] {
        &self.rings
    }

    /// The centre of ring `index`, the mean position of its atoms.
    #[must_use]
    pub fn ring_centre(&self, index: usize) -> Option<Point> {
        let ring = self.rings.get(index)?;
        centroid(ring.iter().filter_map(|&atom| self.position(atom)))
    }
}

/// Lays out `molecule` for a skeletal drawing.
///
/// Bond lengths are as [`Layout`] describes. Rings come from a smallest-ring
/// basis, so fused and bridged systems are accepted: macrocycles, cages and
/// bridged polycycles get the drawings the module documentation lists, and any
/// other shape an approximate but complete one.
#[must_use]
pub fn layout(molecule: &Molecule) -> Layout {
    let adjacency = molecule.adjacency();
    let rings = RingInfo::new(molecule, &adjacency);
    let components = molecule.components();
    let mut placer = Placer::new(molecule, &adjacency, &rings);
    let hubs: Vec<Option<usize>> = components
        .iter()
        .map(|component| placer.place(component))
        .collect();
    let mut positions = placer.finish();
    double_bonds::settle(&mut positions, molecule, &adjacency, &rings);
    for (component, hub) in components.iter().zip(hubs) {
        overlap::resolve(&mut positions, component, molecule, &adjacency, &rings, hub);
    }
    arrange(&mut positions, &components);
    Layout {
        positions,
        rings: rings.rings().to_vec(),
    }
}

/// Moves the components into a row, left to right, each centred on `y = 0`.
fn arrange(positions: &mut [Point], components: &[Vec<usize>]) {
    let mut cursor = 0.0;
    for component in components {
        let mut bounds = Bounds::default();
        component
            .iter()
            .for_each(|&atom| bounds.include(positions[atom]));
        let Some((min, max)) = bounds.extent() else {
            continue;
        };
        let shift = Point::new(cursor - min.x, -(min.y + max.y) / 2.0);
        for &atom in component {
            positions[atom] = positions[atom] + shift;
        }
        cursor += max.x - min.x + COMPONENT_GAP;
    }
}

#[cfg(test)]
mod tests;
