use moiety_core::{Molecule, smiles};

use super::*;
use crate::graph;

const TOLERANCE: f64 = 1e-6;

fn read(input: &str) -> Molecule {
    moiety::interpret(input)
        .top()
        .unwrap_or_else(|| panic!("{input} has a reading"))
        .molecule
        .clone()
}

/// Lays out `molecule` and checks the invariants every layout must keep.
fn checked(name: &str, molecule: &Molecule) -> Layout {
    let layout = layout(molecule);
    let points = layout.positions();
    assert_eq!(points.len(), molecule.atoms().len(), "{name}");
    assert!(
        points.iter().all(|p| p.x.is_finite() && p.y.is_finite()),
        "{name}"
    );
    for bond in molecule.bonds() {
        let length = points[bond.a].distance(points[bond.b]);
        assert!(
            (length - 1.0).abs() < TOLERANCE,
            "{name}: bond {bond:?} has length {length}"
        );
    }
    for (i, a) in points.iter().enumerate() {
        for (j, b) in points.iter().enumerate().skip(i + 1) {
            assert!(
                a.distance(*b) >= 0.5,
                "{name}: atoms {i} and {j} are too close"
            );
        }
    }
    for (index, ring) in layout.rings().iter().enumerate() {
        assert_regular(name, &layout, index, ring);
    }
    assert_components_apart(name, molecule, &layout);
    layout
}

fn assert_regular(name: &str, layout: &Layout, index: usize, ring: &[usize]) {
    let centre = layout.ring_centre(index).unwrap();
    let points = layout.positions();
    let radius = points[ring[0]].distance(centre);
    for (i, &atom) in ring.iter().enumerate() {
        let next = ring[(i + 1) % ring.len()];
        assert!(
            (points[atom].distance(centre) - radius).abs() < TOLERANCE,
            "{name}: ring {index}"
        );
        assert!(
            (points[atom].distance(points[next]) - 1.0).abs() < TOLERANCE,
            "{name}: ring {index}"
        );
    }
}

fn assert_components_apart(name: &str, molecule: &Molecule, layout: &Layout) {
    let spans: Vec<(f64, f64)> = molecule
        .components()
        .iter()
        .map(|component| {
            let xs = component.iter().map(|&a| layout.positions()[a].x);
            let min = xs.clone().fold(f64::INFINITY, f64::min);
            (min, xs.fold(f64::NEG_INFINITY, f64::max))
        })
        .collect();
    for pair in spans.windows(2) {
        assert!(
            pair[0].1 + 1.0 < pair[1].0,
            "{name}: components overlap: {spans:?}"
        );
    }
}

/// Angle at `centre` between the bonds to `a` and `b`, in degrees.
fn angle(layout: &Layout, centre: usize, a: usize, b: usize) -> f64 {
    let p = layout.positions();
    let (u, v) = (p[a] - p[centre], p[b] - p[centre]);
    u.cross(v).atan2(u.dot(v)).abs().to_degrees()
}

/// Every angle at every acyclic atom: 120 for two or three neighbours, 90 or 180 for four, 180 when straight.
fn assert_chain_angles(name: &str, molecule: &Molecule, layout: &Layout) {
    let adjacency = molecule.adjacency();
    let in_ring = |atom: usize| layout.rings().iter().any(|ring| ring.contains(&atom));
    for atom in (0..molecule.atoms().len()).filter(|&a| !in_ring(a)) {
        let neighbours: Vec<usize> = adjacency[atom].iter().map(|&(n, _)| n).collect();
        let allowed: &[f64] = match neighbours.len() {
            2 if graph::is_linear(molecule, &adjacency, atom) => &[180.0],
            2 | 3 => &[120.0],
            4 => &[90.0, 180.0],
            _ => continue,
        };
        for (i, &a) in neighbours.iter().enumerate() {
            for &b in &neighbours[i + 1..] {
                let measured = angle(layout, atom, a, b);
                assert!(
                    allowed.iter().any(|&want| (measured - want).abs() < 1e-4),
                    "{name}: angle {a}-{atom}-{b} is {measured}"
                );
            }
        }
    }
}

#[test]
fn hexane_is_a_horizontal_zig_zag() {
    let molecule = read("CH3(CH2)4CH3");
    let layout = checked("hexane", &molecule);
    assert_chain_angles("hexane", &molecule, &layout);
    let p = layout.positions();
    let turns: Vec<f64> = (1..5)
        .map(|i| (p[i] - p[i - 1]).cross(p[i + 1] - p[i]).signum())
        .collect();
    assert!(
        turns.windows(2).all(|w| w[0] == -w[1]),
        "turns alternate: {turns:?}"
    );
    assert!(p[5].x - p[0].x > 4.0, "the chain runs left to right");
    assert!(
        p.iter().all(|q| q.y.abs() < 0.5),
        "the chain stays horizontal"
    );
}

#[test]
fn tert_butanol_branches_at_right_angles() {
    let molecule = read("(CH3)3COH");
    let layout = checked("tert-butanol", &molecule);
    assert_chain_angles("tert-butanol", &molecule, &layout);
}

#[test]
fn benzyl_alcohol_points_its_substituent_away_from_the_ring() {
    let molecule = read("PhCH2OH");
    let layout = checked("benzyl alcohol", &molecule);
    assert_eq!(layout.rings().len(), 1);
    assert_chain_angles("benzyl alcohol", &molecule, &layout);
    let ring = &layout.rings()[0];
    let attached = *ring.iter().find(|&&a| a == 0).unwrap();
    let centre = layout.ring_centre(0).unwrap();
    let p = layout.positions();
    let outward = (p[attached] - centre).normalized().unwrap();
    let bond = (p[6] - p[attached]).normalized().unwrap();
    assert!((outward.dot(bond) - 1.0).abs() < TOLERANCE);
}

#[test]
fn trityl_alcohol_keeps_three_regular_rings_apart() {
    let molecule = read("TrOH");
    let layout = checked("trityl alcohol", &molecule);
    assert_eq!(layout.rings().len(), 3);
    assert_chain_angles("trityl alcohol", &molecule, &layout);
}

#[test]
fn hydrated_salt_components_sit_side_by_side() {
    let molecule = read("CuSO4·5H2O");
    assert_eq!(molecule.components().len(), 7);
    let layout = checked("copper sulfate pentahydrate", &molecule);
    assert_chain_angles("copper sulfate pentahydrate", &molecule, &layout);
    let xs: Vec<f64> = molecule
        .components()
        .iter()
        .map(|c| layout.positions()[c[0]].x)
        .collect();
    assert!(
        xs.windows(2).all(|w| w[0] < w[1]),
        "components keep input order: {xs:?}"
    );
}

#[test]
fn triple_bonds_and_cumulated_double_bonds_are_straight() {
    for (input, centre) in [("HC≡CCH3", 1), ("CO2", 0), ("CH2=C=CH2", 1)] {
        let molecule = read(input);
        let layout = checked(input, &molecule);
        let neighbours: Vec<usize> = molecule.adjacency()[centre]
            .iter()
            .map(|&(n, _)| n)
            .collect();
        let measured = angle(&layout, centre, neighbours[0], neighbours[1]);
        assert!((measured - 180.0).abs() < 1e-6, "{input}: {measured}");
        assert!(
            layout.positions().iter().all(|p| p.y.abs() < TOLERANCE),
            "{input} lies flat"
        );
    }
}

#[test]
fn acetate_has_a_trigonal_carboxyl_carbon() {
    for reading in moiety::interpret("CH3COO-").readings() {
        let layout = checked(&reading.smiles, &reading.molecule);
        assert_chain_angles(&reading.smiles, &reading.molecule, &layout);
    }
}

#[test]
fn naphthalene_rings_lie_on_opposite_sides_of_the_shared_edge() {
    let molecule = smiles::parse("c1ccc2ccccc2c1").unwrap();
    let layout = checked("naphthalene", &molecule);
    assert_eq!(layout.rings().len(), 2);
    let shared: Vec<usize> = layout.rings()[0]
        .iter()
        .copied()
        .filter(|a| layout.rings()[1].contains(a))
        .collect();
    assert_eq!(shared.len(), 2);
    let p = layout.positions();
    let edge = p[shared[1]] - p[shared[0]];
    let sides: Vec<f64> = (0..2)
        .map(|r| edge.cross(layout.ring_centre(r).unwrap() - p[shared[0]]))
        .collect();
    assert!(
        sides[0] * sides[1] < 0.0,
        "ring centres straddle the shared edge"
    );
    let centres = layout
        .ring_centre(0)
        .unwrap()
        .distance(layout.ring_centre(1).unwrap());
    assert!((centres - 3f64.sqrt()).abs() < TOLERANCE);
}

#[test]
fn other_ring_shapes_get_a_complete_drawing() {
    for input in [
        "C1CC1CO",
        "C1CCCCO1",
        "c1cc2ccc3cccc4ccc(c1)c2c34",
        "C1CCC2(CC1)CCCC2",
        "C1CC2CCC1C2",
        "C12C3C4C1C5C2C3C45",
        "C1CCCCCCCCCCC1",
        "c1ccccc1-c1ccccc1",
        "[Na+].[Cl-]",
        "FS(F)(F)(F)(F)F",
        "*CC*",
    ] {
        let molecule = smiles::parse(input).unwrap();
        let layout = layout(&molecule);
        assert_eq!(layout.positions().len(), molecule.atoms().len(), "{input}");
        assert!(
            layout
                .positions()
                .iter()
                .all(|p| p.x.is_finite() && p.y.is_finite()),
            "{input}"
        );
    }
    for input in [
        "C1CC1CO",
        "C1CCCCO1",
        "c1cc2ccc3cccc4ccc(c1)c2c34",
        "C1CCC2(CC1)CCCC2",
    ] {
        checked(input, &smiles::parse(input).unwrap());
    }
    let hexafluoride = smiles::parse("FS(F)(F)(F)(F)F").unwrap();
    checked("sulfur hexafluoride", &hexafluoride);
}

#[test]
fn layout_is_deterministic_and_handles_empty_molecules() {
    let molecule = read("TrOH");
    assert_eq!(layout(&molecule), layout(&molecule));
    assert_eq!(layout(&Molecule::new()), Layout::default());
    let methane = layout(&read("CH4"));
    assert_eq!(methane.positions(), [Point::new(0.0, 0.0)]);
}
