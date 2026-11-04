use moiety_core::{Atom, Bond, BondOrder, Element, Molecule, smiles};

use super::*;

fn top(input: &str) -> Reading {
    moiety::interpret(input)
        .top()
        .unwrap_or_else(|| panic!("{input} has a reading"))
        .clone()
}

/// Checks that every tag is closed in order and no markup character is left raw.
fn assert_well_formed(svg: &str) {
    let mut open: Vec<&str> = Vec::new();
    let mut rest = svg;
    while let Some(start) = rest.find('<') {
        assert!(
            !rest[..start].contains('>'),
            "raw > in text before {}",
            &rest[start..]
        );
        let end = start + rest[start..].find('>').expect("every tag is closed");
        let tag = &rest[start + 1..end];
        assert!(!tag.contains('<'), "raw < inside <{tag}>");
        assert_eq!(
            tag.matches('"').count() % 2,
            0,
            "unbalanced quotes in <{tag}>"
        );
        if let Some(name) = tag.strip_prefix('/') {
            assert_eq!(open.pop(), Some(name), "mismatched </{name}>");
        } else if !tag.ends_with('/') {
            open.push(tag.split_whitespace().next().expect("tags have names"));
        }
        rest = &rest[end + 1..];
    }
    assert!(open.is_empty(), "unclosed tags: {open:?}");
    assert!(!rest.contains('>'));
}

fn assert_complete(name: &str, molecule: &Molecule, svg: &str) {
    assert_well_formed(svg);
    for index in 0..molecule.atoms().len() {
        assert!(
            svg.contains(&format!("data-atom=\"{index}\"")),
            "{name}: atom {index}"
        );
    }
    for index in 0..molecule.bonds().len() {
        assert!(
            svg.contains(&format!("data-bond=\"{index}\"")),
            "{name}: bond {index}"
        );
    }
}

#[test]
fn drawings_are_well_formed_and_complete() {
    for input in [
        "CH3(CH2)4CH3",
        "(CH3)3COH",
        "PhCH2OH",
        "TrOH",
        "CuSO4·5H2O",
        "HC≡CCH3",
        "CO2",
        "CH3COO-",
        "CH3NO2",
        "-CH2CH3",
        "CH4",
    ] {
        for reading in moiety::interpret(input).readings() {
            let svg = render_reading_svg(reading, &SvgOptions::default());
            assert_complete(input, &reading.molecule, &svg);
            assert!(svg.starts_with("<svg xmlns=\"http://www.w3.org/2000/svg\""));
        }
    }
    for input in [
        "c1ccc2ccccc2c1",
        "C12C3C4C1C5C2C3C45",
        "[Na+].[Cl-]",
        "FS(F)(F)(F)(F)F",
    ] {
        let molecule = smiles::parse(input).unwrap();
        assert_complete(
            input,
            &molecule,
            &render_svg(&molecule, &SvgOptions::default()),
        );
    }
}

#[test]
fn readings_carry_their_input_spans() {
    let reading = top("CH3CH2OH");
    let svg = render_reading_svg(&reading, &SvgOptions::default());
    assert!(svg.contains("data-atom=\"0\" data-span=\"0-3\""));
    assert!(svg.contains("data-atom=\"2\" data-span=\"6-8\""));
    assert!(svg.contains("data-bond=\"1\" data-span=\"3-8\""));
    assert!(!render_svg(&reading.molecule, &SvgOptions::default()).contains("data-span"));
}

#[test]
fn titles_and_labels_are_escaped() {
    let options = SvgOptions {
        title: Some("<script>alert(\"x\") & 'y'</script>".to_string()),
        ..SvgOptions::default()
    };
    let svg = render_svg(&top("CH3CH2OH").molecule, &options);
    assert!(svg.contains(
        "<title>&lt;script&gt;alert(&quot;x&quot;) &amp; &#39;y&#39;&lt;/script&gt;</title>"
    ));
    assert!(!svg.contains("<script"));
    assert_well_formed(&svg);
}

#[test]
fn only_heteroatoms_and_lone_carbons_are_labeled() {
    let ethane = render_svg(&top("CH3CH3").molecule, &SvgOptions::default());
    assert!(!ethane.contains("<text"));
    let methane = render_svg(&top("CH4").molecule, &SvgOptions::default());
    assert!(methane.contains(">C</text>"));
    assert!(methane.contains(">H<tspan dy=\"4.8\" font-size=\"11.2\">4</tspan></text>"));
    let water = render_svg(&smiles::parse("O").unwrap(), &SvgOptions::default());
    assert!(
        water.contains("text-anchor=\"end\">H<tspan"),
        "H2O puts its hydrogens first"
    );
}

#[test]
fn hydrogens_sit_on_the_outside() {
    let ethanol = render_svg(&smiles::parse("CCO").unwrap(), &SvgOptions::default());
    assert!(
        ethanol.contains("text-anchor=\"start\">H</text>"),
        "OH: {ethanol}"
    );
    let reversed = render_svg(&smiles::parse("OCC").unwrap(), &SvgOptions::default());
    assert!(
        reversed.contains("text-anchor=\"end\">H</text>"),
        "HO: {reversed}"
    );
}

#[test]
fn charges_are_superscripts() {
    let acetate = render_svg(&top("CH3COO-").molecule, &SvgOptions::default());
    assert!(acetate.contains("font-size=\"11.2\">\u{2212}</tspan>"));
    let copper = render_svg(&top("CuSO4").molecule, &SvgOptions::default());
    assert!(copper.contains(">2+</tspan>"));
}

#[test]
fn multiple_bonds_get_parallel_strokes() {
    let strokes = |svg: &str, order: &str| -> Vec<usize> {
        svg.split("<g class=\"bond ")
            .skip(1)
            .filter(|group| group.starts_with(order))
            .map(|group| {
                group
                    .split("</g>")
                    .next()
                    .unwrap()
                    .matches("<line x1")
                    .count()
            })
            .collect()
    };
    let benzyl = render_svg(&top("PhCH2OH").molecule, &SvgOptions::default());
    assert_eq!(strokes(&benzyl, "double"), [2, 2, 2]);
    assert!(strokes(&benzyl, "aromatic").is_empty());
    let propyne = render_svg(&top("HC≡CCH3").molecule, &SvgOptions::default());
    assert_eq!(strokes(&propyne, "triple"), [3]);
    let dioxide = render_svg(&top("CO2").molecule, &SvgOptions::default());
    assert_eq!(strokes(&dioxide, "double"), [2, 2]);
}

#[test]
fn rings_without_a_kekule_structure_get_a_circle() {
    let atom = Atom {
        hydrogens: 1,
        aromatic: true,
        ..Atom::new(Element::C)
    };
    let bonds = (0..5)
        .map(|i| Bond {
            a: i,
            b: (i + 1) % 5,
            order: BondOrder::Aromatic,
        })
        .collect();
    let molecule = Molecule::from_parts(vec![atom; 5], bonds);
    let svg = render_svg(&molecule, &SvgOptions::default());
    assert!(svg.contains("<circle class=\"aromatic-ring\""));
    assert_complete("odd aromatic ring", &molecule, &svg);
}

#[test]
fn attachment_points_are_wavy_lines() {
    let svg = render_reading_svg(&top("-CH2CH3"), &SvgOptions::default());
    assert!(svg.contains("<path class=\"attachment\""));
    assert!(!svg.contains("<text"));
}

#[test]
fn falls_back_to_default_sizes_and_draws_empty_molecules() {
    let options = SvgOptions {
        bond_length: -1.0,
        font_size: f64::NAN,
        line_width: 0.0,
        padding: f64::INFINITY,
        title: None,
    };
    let svg = render_svg(&Molecule::new(), &options);
    assert!(svg.contains("viewBox=\"-12 -12 24 24\""));
    assert!(svg.contains("<title>Skeletal structure of </title>"));
    assert_well_formed(&svg);
}

#[test]
fn view_box_fits_the_drawing() {
    let svg = render_svg(&top("CH3(CH2)4CH3").molecule, &SvgOptions::default());
    let view_box: Vec<f64> = svg
        .split("viewBox=\"")
        .nth(1)
        .and_then(|rest| rest.split('"').next())
        .unwrap()
        .split(' ')
        .map(|n| n.parse().unwrap())
        .collect();
    assert!((view_box[2] - (5.0 * 40.0 * 3f64.sqrt() / 2.0 + 24.0)).abs() < 1.0);
    assert!((view_box[3] - (20.0 + 24.0)).abs() < 1.0);
    assert_eq!(num(1.0), "1");
    assert_eq!(num(-0.001), "0");
    assert_eq!(num(2.345_67), "2.35");
}

/// The strokes of every bond of the given order, as (x1, y1, x2, y2), keyed by bond index.
fn strokes_of(svg: &str, order: &str) -> Vec<(usize, Vec<[f64; 4]>)> {
    let mut found = Vec::new();
    let marker = format!("<g class=\"bond {order}\" data-bond=\"");
    let mut rest = svg;
    while let Some(start) = rest.find(&marker) {
        rest = &rest[start + marker.len()..];
        let index: usize = rest[..rest.find('"').expect("closing quote")]
            .parse()
            .expect("bond index");
        let group = &rest[..rest.find("</g>").expect("group end")];
        let lines = group
            .split("<line ")
            .skip(1)
            .filter(|line| !line.starts_with("class=\"hit\""))
            .map(|line| {
                let value = |name: &str| -> f64 {
                    let at = line.find(&format!("{name}=\"")).expect("coordinate") + name.len() + 2;
                    line[at..at + line[at..].find('"').expect("closing quote")]
                        .parse()
                        .expect("number")
                };
                [value("x1"), value("y1"), value("x2"), value("y2")]
            })
            .collect();
        found.push((index, lines));
    }
    found
}

/// Positions of a stroke's two ends along the bond axis from atom `a` to atom `b`.
fn along_axis(stroke: [f64; 4], a: Point, b: Point) -> (f64, f64) {
    let axis = (b - a).normalized().expect("bonded atoms are apart");
    let project = |x: f64, y: f64| (Point::new(x, y) - a).dot(axis);
    (project(stroke[0], stroke[1]), project(stroke[2], stroke[3]))
}

#[test]
fn parallel_strokes_end_level_at_labels() {
    for input in [
        "TsCl",
        "CH3SO2Cl",
        "CH3S(O)CH3",
        "CH3CH2CN",
        "NCCH2CN",
        "PhNO2",
    ] {
        let reading = top(input);
        let molecule = &reading.molecule;
        let options = SvgOptions::default();
        let svg = render_svg(molecule, &options);
        let positions = crate::layout::layout(molecule);
        for order in ["double", "triple"] {
            for (index, lines) in strokes_of(&svg, order) {
                let bond = molecule.bonds()[index];
                let labelled = |atom: usize| molecule.atoms()[atom].element() != Some(Element::C);
                let placed = |atom: usize| positions.position(atom).expect("every atom is placed");
                let (a, b) = (placed(bond.a), placed(bond.b));
                let scale = options.bond_length;
                let ends: Vec<(f64, f64)> = lines
                    .iter()
                    .map(|&line| along_axis(line, a * scale, b * scale))
                    .collect();
                for pair in ends.windows(2) {
                    if labelled(bond.a) {
                        assert!(
                            (pair[0].0 - pair[1].0).abs() < 0.01,
                            "{input} bond {index} start {ends:?}"
                        );
                    }
                    if labelled(bond.b) {
                        assert!(
                            (pair[0].1 - pair[1].1).abs() < 0.01,
                            "{input} bond {index} end {ends:?}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn bonds_meet_double_bonds_without_gaps() {
    for input in [
        "CH2=CHCH=CH2",
        "CH2=CHCH3",
        "CH3COCH3",
        "CH3COOH",
        "(CH3)2C=CH2",
    ] {
        let svg = render_reading_svg(&top(input), &SvgOptions::default());
        let singles: Vec<[f64; 4]> = strokes_of(&svg, "single")
            .into_iter()
            .flat_map(|(_, lines)| lines)
            .collect();
        let doubles: Vec<[f64; 4]> = strokes_of(&svg, "double")
            .into_iter()
            .flat_map(|(_, lines)| lines)
            .collect();
        let on_segment = |p: Point, l: &[f64; 4]| {
            let (a, b) = (Point::new(l[0], l[1]), Point::new(l[2], l[3]));
            let t = ((p - a).dot(b - a) / (b - a).dot(b - a)).clamp(0.0, 1.0);
            p.distance(a + (b - a) * t) < 0.01
        };
        let ends = |l: &[f64; 4]| [Point::new(l[0], l[1]), Point::new(l[2], l[3])];
        let joined = doubles.iter().any(|double| {
            ends(double)
                .iter()
                .any(|&p| singles.iter().any(|single| on_segment(p, single)))
        });
        assert!(
            joined,
            "{input}: no single bond meets a double-bond stroke\n{svg}"
        );
    }
}

#[test]
fn draws_one_wedge_or_hash_per_stereocentre() {
    let molecule = moiety_core::smiles::parse("C[C@@H](C(=O)O)N").unwrap();
    let svg = render_svg(&molecule, &SvgOptions::default());
    let wedges = svg.matches("bond single wedge").count() + svg.matches("bond single hash").count();
    assert_eq!(wedges, 1);
    let mirror = render_svg(
        &moiety_core::smiles::parse("C[C@H](C(=O)O)N").unwrap(),
        &SvgOptions::default(),
    );
    assert_ne!(
        svg.contains("bond single wedge"),
        mirror.contains("bond single wedge"),
        "the enantiomers differ"
    );
    let flat = render_svg(
        &moiety_core::smiles::parse("CC(C(=O)O)N").unwrap(),
        &SvgOptions::default(),
    );
    assert!(!flat.contains("wedge") && !flat.contains("hash"));
}

#[test]
fn draws_the_hydrogens_where_fused_rings_meet() {
    let molecule = moiety_core::smiles::parse("C1CC[C@H]2CCCC[C@@H]2C1").unwrap();
    let svg = render_svg(&molecule, &SvgOptions::default());
    assert_well_formed(&svg);
    assert_eq!(svg.matches(">H</text>").count(), 2, "{svg}");
    assert_eq!(svg.matches("bond single wedge").count(), 1);
    assert_eq!(svg.matches("bond single hash").count(), 1);
    assert!(
        svg.contains("data-atom=\"11\""),
        "the hydrogens follow the atoms"
    );
}

#[test]
fn small_drawings_wedge_a_ring_bond_instead() {
    let molecule = moiety_core::smiles::parse("C1CC[C@H]2CCCC[C@@H]2C1").unwrap();
    let options = SvgOptions {
        bond_length: 20.0,
        ..SvgOptions::default()
    };
    let svg = render_svg(&molecule, &options);
    assert!(!svg.contains(">H</text>"), "{svg}");
    let wedges = svg.matches("bond single wedge").count() + svg.matches("bond single hash").count();
    assert_eq!(wedges, 2);
}
