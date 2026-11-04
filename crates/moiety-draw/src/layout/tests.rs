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

#[test]
fn lays_out_large_readings_with_unit_bonds() {
    let molecule = read("-(BNO-O2NCHBNCH2)66H");
    let layout = layout(&molecule);
    let points = layout.positions();
    assert_eq!(points.len(), 661);
    assert!(points.iter().all(|p| p.x.is_finite() && p.y.is_finite()));
    for bond in molecule.bonds() {
        assert!((points[bond.a].distance(points[bond.b]) - 1.0).abs() < TOLERANCE);
    }
}

/// The direction of the bond from `from` to `to`, in degrees clockwise from the positive `x` axis.
fn heading(layout: &Layout, from: usize, to: usize) -> f64 {
    (layout.positions()[to] - layout.positions()[from])
        .angle()
        .to_degrees()
}

fn neighbours(molecule: &Molecule, atom: usize) -> Vec<usize> {
    molecule.adjacency()[atom].iter().map(|&(n, _)| n).collect()
}

fn element_index(molecule: &Molecule, symbol: &str) -> usize {
    molecule
        .atoms()
        .iter()
        .position(|atom| atom.element().is_some_and(|e| e.symbol() == symbol))
        .unwrap_or_else(|| panic!("the molecule has {symbol}"))
}

#[test]
fn hubs_spread_alike_branches_evenly_with_one_pointing_up() {
    for (input, centre) in [
        ("P(C6H5)3", "P"),
        ("N(CH2CH2OH)3", "N"),
        ("Al(OiPr)3", "Al"),
    ] {
        let molecule = read(input);
        let layout = checked(input, &molecule);
        assert_chain_angles(input, &molecule, &layout);
        let hub = element_index(&molecule, centre);
        let headings: Vec<f64> = neighbours(&molecule, hub)
            .into_iter()
            .map(|n| heading(&layout, hub, n))
            .collect();
        assert!(
            headings.iter().any(|h| (h + 90.0).abs() < 1e-6),
            "{input}: one branch points up, {headings:?}"
        );
    }
}

#[test]
fn four_way_hubs_are_crosses_with_the_odd_branch_up() {
    for (input, centre, up) in [("Ph3PO", "P", Some("O")), ("C(CH2OH)4", "C", None)] {
        let molecule = read(input);
        let layout = checked(input, &molecule);
        assert_chain_angles(input, &molecule, &layout);
        let hub = if centre == "C" {
            (0..molecule.atoms().len())
                .find(|&a| neighbours(&molecule, a).len() == 4)
                .unwrap()
        } else {
            element_index(&molecule, centre)
        };
        for n in neighbours(&molecule, hub) {
            let h = heading(&layout, hub, n).rem_euclid(90.0);
            assert!(h.min(90.0 - h) < 1e-6, "{input}: bonds lie on the axes");
        }
        if let Some(symbol) = up {
            let odd = element_index(&molecule, symbol);
            assert!((heading(&layout, hub, odd) + 90.0).abs() < 1e-6, "{input}");
        }
    }
}

#[test]
fn crowded_hubs_stretch_their_bonds_and_keep_even_angles() {
    let molecule = read("B(C6F5)3");
    let layout = layout(&molecule);
    let points = layout.positions();
    let hub = element_index(&molecule, "B");
    let arms = neighbours(&molecule, hub);
    let lengths: Vec<f64> = arms
        .iter()
        .map(|&n| points[hub].distance(points[n]))
        .collect();
    assert!(
        lengths
            .iter()
            .all(|&l| l > 1.0 + TOLERANCE && l <= 1.5 + TOLERANCE)
    );
    assert!(lengths.iter().all(|&l| (l - lengths[0]).abs() < TOLERANCE));
    for bond in molecule.bonds().iter().filter(|b| b.a != hub && b.b != hub) {
        assert!((points[bond.a].distance(points[bond.b]) - 1.0).abs() < TOLERANCE);
    }
    for (i, &a) in arms.iter().enumerate() {
        for &b in &arms[i + 1..] {
            assert!((angle(&layout, hub, a, b) - 120.0).abs() < 1e-6);
        }
    }
    for (i, a) in points.iter().enumerate() {
        for b in &points[i + 1..] {
            assert!(a.distance(*b) >= 0.6, "no fluorines collide");
        }
    }
    for (index, ring) in layout.rings().iter().enumerate() {
        assert_regular("B(C6F5)3", &layout, index, ring);
    }
}

#[test]
fn groups_without_three_alike_branches_are_not_hubs() {
    let molecule = read("(CH3)3CCOCl");
    let layout = checked("pivaloyl chloride", &molecule);
    let tertiary = (0..molecule.atoms().len())
        .find(|&a| neighbours(&molecule, a).len() == 4)
        .unwrap();
    let acyl = neighbours(&molecule, tertiary)
        .into_iter()
        .find(|&n| neighbours(&molecule, n).len() == 3)
        .unwrap();
    assert!(
        heading(&layout, tertiary, acyl).abs() < 1e-6,
        "the chain still runs to the right"
    );
}

#[test]
fn hubs_that_turns_can_part_keep_unit_bonds() {
    for (name, molecule) in [
        ("trityl tert-butane", read("(C6H5)3CC(CH3)3")),
        (
            "tert-butyltriphenylsilane",
            smiles::parse("CC(C)(C)[Si](c1ccccc1)(c1ccccc1)c1ccccc1").unwrap(),
        ),
        (
            "TIPS triflate",
            smiles::parse("FC(F)(F)S(=O)(=O)O[Si](C(C)C)(C(C)C)C(C)C").unwrap(),
        ),
    ] {
        checked(name, &molecule);
    }
}

/// Pairs of bonds with no atom in common whose segments cross.
fn crossings(molecule: &Molecule, layout: &Layout) -> Vec<(usize, usize)> {
    let p = layout.positions();
    let side = |o: Point, a: Point, b: Point| (a - o).cross(b - o);
    let bonds = molecule.bonds();
    let mut found = Vec::new();
    for (i, x) in bonds.iter().enumerate() {
        for (j, y) in bonds.iter().enumerate().skip(i + 1) {
            if [x.a, x.b].iter().any(|a| *a == y.a || *a == y.b) {
                continue;
            }
            let (a, b, c, d) = (p[x.a], p[x.b], p[y.a], p[y.b]);
            let apart = |s: f64, t: f64| (s > 1e-9 && t < -1e-9) || (s < -1e-9 && t > 1e-9);
            if apart(side(c, d, a), side(c, d, b)) && apart(side(a, b, c), side(a, b, d)) {
                found.push((i, j));
            }
        }
    }
    found
}

#[test]
fn bridged_ring_systems_draw_without_crossings() {
    for (name, input) in [
        ("bicyclo[2.2.2]octane", "C1CC2CCC1CC2"),
        ("quinuclidine", "C1CN2CCC1CC2"),
        ("tropine", "CN1C2CCC1CC(O)C2"),
        ("cocaine", "CN1C2CCC1C(C(C2)OC(=O)c1ccccc1)C(=O)OC"),
        ("quinine", "COc1ccc2nccc(C(O)C3CC4CCN3CC4C=C)c2c1"),
    ] {
        let molecule = smiles::parse(input).unwrap();
        let layout = layout(&molecule);
        let p = layout.positions();
        assert_eq!(crossings(&molecule, &layout), [], "{name}: bonds cross");
        let adjacency = molecule.adjacency();
        for i in 0..p.len() {
            for j in i + 1..p.len() {
                if graph::bond_between(&adjacency, i, j).is_none() {
                    assert!(
                        p[i].distance(p[j]) >= 0.5,
                        "{name}: atoms {i} and {j} are too close"
                    );
                }
            }
        }
        for bond in molecule.bonds() {
            let length = p[bond.a].distance(p[bond.b]);
            assert!(
                (0.55..=1.6).contains(&length),
                "{name}: bond {bond:?} has length {length}"
            );
        }
    }
}

#[test]
fn templates_draw_porphine_morphine_and_norbornane_as_chemists_do() {
    let porphine = smiles::parse("C1=CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2").unwrap();
    let drawn = layout(&porphine);
    assert_eq!(crossings(&porphine, &drawn), []);
    let p = drawn.positions();
    let centre = super::centroid(p.iter().copied()).unwrap();
    let mut by_distance: Vec<usize> = (0..p.len()).collect();
    by_distance.sort_by(|&a, &b| p[a].distance(centre).total_cmp(&p[b].distance(centre)));
    let inner: Vec<bool> = by_distance[..4]
        .iter()
        .map(|&a| {
            porphine.atoms()[a]
                .element()
                .is_some_and(|e| e.symbol() == "N")
        })
        .collect();
    assert_eq!(inner, [true; 4], "the nitrogens point to the centre");

    let morphine = smiles::parse("CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O").unwrap();
    assert!(crossings(&morphine, &layout(&morphine)).len() <= 1);

    let camphor = smiles::parse("CC1(C)C2CCC1(C)C(=O)C2").unwrap();
    let drawn = layout(&camphor);
    let longest = camphor
        .bonds()
        .iter()
        .map(|b| drawn.positions()[b.a].distance(drawn.positions()[b.b]))
        .fold(0.0, f64::max);
    assert!(longest > 1.2, "the bridge is drawn in perspective");
}

#[test]
fn steroids_are_drawn_the_usual_way_round() {
    for (name, input) in [
        (
            "cholesterol",
            "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
        ),
        (
            "testosterone",
            "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C",
        ),
        (
            "testosterone from the A ring",
            "C[C@@]12CCC(=O)C=C1CC[C@@H]1[C@@H]2CC[C@]2(C)[C@@H](O)CC[C@@H]12",
        ),
        (
            "estradiol",
            "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O",
        ),
    ] {
        let molecule = smiles::parse(input).unwrap();
        let drawn = layout(&molecule);
        let centres: Vec<(usize, Point)> = drawn
            .rings()
            .iter()
            .enumerate()
            .map(|(i, ring)| (ring.len(), drawn.ring_centre(i).unwrap()))
            .collect();
        let (_, d) = centres.iter().find(|(size, _)| *size == 5).unwrap();
        for (_, centre) in centres.iter().filter(|(size, _)| *size == 6) {
            assert!(centre.x < d.x, "{name}: ring D is on the right");
            assert!(
                centre.y > d.y - 1.0,
                "{name}: ring D is not below the others"
            );
        }
        for (ring_atom, methyl) in angular_methyls(&molecule) {
            let angle = heading(&drawn, ring_atom, methyl);
            assert!(
                (angle + 90.0).abs() < 10.0,
                "{name}: the angular methyl {methyl} points up, not {angle}"
            );
        }
    }
}

/// Each methyl on an atom with three other neighbours, with that atom.
fn angular_methyls(molecule: &Molecule) -> Vec<(usize, usize)> {
    let adjacency = molecule.adjacency();
    adjacency
        .iter()
        .enumerate()
        .filter_map(|(atom, links)| match links.as_slice() {
            &[(hub, _)] if adjacency[hub].len() == 4 => Some((hub, atom)),
            _ => None,
        })
        .filter(|&(_, atom)| {
            molecule.atoms()[atom]
                .element()
                .is_some_and(|e| e.symbol() == "C")
        })
        .collect()
}

/// Paclitaxel takes PubChem's drawing, where the methyls of the bridging C15
/// cross the bonds beside it, as they do in the usual drawing; artemisinin
/// takes one where nothing crosses.
#[test]
fn taxanes_and_artemisinin_take_their_usual_drawings() {
    let paclitaxel = smiles::parse(
        "CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
    )
    .unwrap();
    let drawn = layout(&paclitaxel);
    let methyls = [15, 16];
    let crossed = crossings(&paclitaxel, &drawn);
    assert!(!crossed.is_empty() && crossed.len() <= 2);
    for (x, y) in crossed {
        let bonds = [paclitaxel.bonds()[x], paclitaxel.bonds()[y]];
        assert!(
            bonds
                .iter()
                .any(|b| methyls.contains(&b.a) || methyls.contains(&b.b)),
            "only the C15 methyls cross: {bonds:?}"
        );
    }
    let artemisinin =
        smiles::parse("C[C@@H]1CC[C@H]2[C@H](C(=O)O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)C").unwrap();
    assert_eq!(crossings(&artemisinin, &layout(&artemisinin)), []);
}

#[test]
fn double_bonds_keep_the_configuration_their_smiles_gives() {
    for input in [
        "F/C=C\\F",
        "F/C=C/F",
        "CCCC/C=C\\CCCC",
        "OC(=O)CCCCCCC/C=C\\CCCCCCCC",
        "C/C=C/C=C\\C=C/C",
        "Cc1ccccc1/C=C\\c1ccccc1",
        "C(C(C)C)C(=C/C(C)C)/C(=C\\C(C)C)CC(C)C",
        "CC1(C)/C=C/C(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C1(C)C",
    ] {
        let molecule = smiles::parse(input).unwrap();
        let drawn = layout(&molecule);
        let p = drawn.positions();
        for stereo in molecule.double_bond_stereo() {
            let bond = molecule.bonds()[stereo.bond];
            let side = |n: usize| (p[bond.b] - p[bond.a]).cross(p[n] - p[bond.a]);
            let (a, b) = (side(stereo.neighbours.0), side(stereo.neighbours.1));
            assert_eq!(a * b > 0.0, stereo.cis, "{input}: bond {}", stereo.bond);
        }
    }
}
