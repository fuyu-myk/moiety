//! SVG drawings of skeletal structures.
//!
//! The drawing follows the usual skeletal conventions:
//!
//! - carbon atoms are bare vertices, unless they stand alone (`CH4`) or carry a charge;
//! - other atoms are labeled with their implicit hydrogens on the side away
//!   from their bonds (`OH`, `HO`, `NH2`) and their charge as a superscript;
//! - bonds stop short of labels;
//! - a double bond in a ring or chain gets a shorter second line on the inner
//!   side, a terminal one such as `C=O` a centred pair, and a triple bond three lines;
//! - aromatic rings are drawn as a Kekulé structure, or with a circle inside
//!   when no Kekulé structure exists;
//! - an attachment point is a short wavy line across the end of its bond.
//!
//! Every atom is drawn as `<g class="atom" data-atom="i">` and every bond as
//! `<g class="bond single" data-bond="j">` (or `double`, `triple`, `aromatic`).
//! A drawing made from a [`Reading`] adds `data-span="start-end"` to both: the
//! byte span of the original input the atom or bond came from. All strokes and
//! fills use `currentColor`, so the drawing takes the text colour of the page
//! it is placed in and reads on light and dark backgrounds alike.

mod atom;
mod bond;
mod label;

use moiety_core::{BondOrder, Molecule, Provenance, Reading, Span, aromatic};

use crate::layout::{Bounds, Layout, Point, layout};
use crate::xml::escape;
use label::Label;

/// Sizes and text of an SVG drawing.
#[derive(Clone, Debug, PartialEq)]
pub struct SvgOptions {
    /// Length of one bond, in pixels.
    pub bond_length: f64,
    /// Font size of atom labels, in pixels.
    pub font_size: f64,
    /// Width of bond lines, in pixels.
    pub line_width: f64,
    /// Empty margin around the drawing, in pixels.
    pub padding: f64,
    /// Accessible title of the drawing. `None` names the molecule by its SMILES.
    pub title: Option<String>,
}

impl Default for SvgOptions {
    fn default() -> Self {
        Self {
            bond_length: 40.0,
            font_size: 16.0,
            line_width: 1.6,
            padding: 12.0,
            title: None,
        }
    }
}

/// Draws `molecule` as a standalone SVG document.
///
/// Sizes in `options` that are not positive and finite fall back to their defaults.
#[must_use]
pub fn render_svg(molecule: &Molecule, options: &SvgOptions) -> String {
    render(molecule, None, options)
}

/// Draws the molecule of `reading`, tagging every atom and bond with the
/// `data-span` of the input it came from.
#[must_use]
pub fn render_reading_svg(reading: &Reading, options: &SvgOptions) -> String {
    render(&reading.molecule, Some(&reading.provenance), options)
}

/// Sizes in pixels, each positive and finite.
#[derive(Clone, Copy, Debug)]
struct Metrics {
    bond: f64,
    font: f64,
    line: f64,
    padding: f64,
}

impl Metrics {
    fn new(options: &SvgOptions) -> Self {
        let defaults = SvgOptions::default();
        let pick = |value: f64, fallback: f64| {
            if value.is_finite() && value > 0.0 {
                value
            } else {
                fallback
            }
        };
        Self {
            bond: pick(options.bond_length, defaults.bond_length),
            font: pick(options.font_size, defaults.font_size),
            line: pick(options.line_width, defaults.line_width),
            padding: pick(options.padding, defaults.padding),
        }
    }
}

/// Everything the atom and bond painters need, in pixel coordinates.
struct Scene {
    molecule: Molecule,
    adjacency: Vec<Vec<(usize, usize)>>,
    points: Vec<Point>,
    labels: Vec<Option<Label>>,
    rings: Vec<Vec<usize>>,
    ring_centres: Vec<Point>,
    circled_rings: Vec<usize>,
    metrics: Metrics,
}

impl Scene {
    fn new(molecule: &Molecule, metrics: Metrics) -> Self {
        let geometry = layout(molecule);
        let (drawable, circled_rings) = kekule(molecule, &geometry);
        let adjacency = drawable.adjacency();
        let points: Vec<Point> = geometry
            .positions()
            .iter()
            .map(|&p| p * metrics.bond)
            .collect();
        let labels = (0..points.len())
            .map(|index| Label::for_atom(&drawable, &adjacency, &points, index))
            .collect();
        let ring_centres = (0..geometry.rings().len())
            .map(|index| geometry.ring_centre(index).unwrap_or_default() * metrics.bond)
            .collect();
        Self {
            molecule: drawable,
            adjacency,
            points,
            labels,
            rings: geometry.rings().to_vec(),
            ring_centres,
            circled_rings,
            metrics,
        }
    }
}

/// The molecule with aromatic bonds kekulized, and the rings that still need a
/// circle because no Kekulé structure exists.
fn kekule(molecule: &Molecule, geometry: &Layout) -> (Molecule, Vec<usize>) {
    let is_aromatic = |bond: &moiety_core::Bond| bond.order == BondOrder::Aromatic;
    if !molecule.bonds().iter().any(is_aromatic) {
        return (molecule.clone(), Vec::new());
    }
    if let Some(kekule) = aromatic::kekulize(molecule) {
        return (kekule, Vec::new());
    }
    let adjacency = molecule.adjacency();
    let circled = geometry
        .rings()
        .iter()
        .enumerate()
        .filter(|(_, ring)| {
            ring.iter().enumerate().all(|(i, &atom)| {
                let next = ring[(i + 1) % ring.len()];
                crate::graph::bond_between(&adjacency, atom, next)
                    .is_some_and(|bond| is_aromatic(&molecule.bonds()[bond]))
            })
        })
        .map(|(index, _)| index)
        .collect();
    (molecule.clone(), circled)
}

fn render(molecule: &Molecule, provenance: Option<&Provenance>, options: &SvgOptions) -> String {
    let metrics = Metrics::new(options);
    let scene = Scene::new(molecule, metrics);
    let mut bounds = Bounds::default();
    let mut bonds = String::new();
    for index in 0..scene.molecule.bonds().len() {
        let span = span_attribute(provenance.and_then(|p| p.bond_source(index)));
        bond::render(&scene, index, &span, &mut bonds, &mut bounds);
    }
    bond::render_circles(&scene, &mut bonds, &mut bounds);
    let mut atoms = String::new();
    for index in 0..scene.points.len() {
        let span = span_attribute(provenance.and_then(|p| p.atom_source(index)));
        atom::render(&scene, index, &span, &mut atoms, &mut bounds);
    }
    let title = options
        .title
        .clone()
        .unwrap_or_else(|| format!("Skeletal structure of {}", molecule.to_smiles()));
    document(&title, metrics, bounds, &bonds, &atoms)
}

fn span_attribute(span: Option<Span>) -> String {
    span.map(|s| format!(" data-span=\"{}-{}\"", s.start, s.end))
        .unwrap_or_default()
}

fn document(title: &str, metrics: Metrics, bounds: Bounds, bonds: &str, atoms: &str) -> String {
    let (min, max) = bounds.extent().unwrap_or_default();
    let origin = min - Point::new(metrics.padding, metrics.padding);
    let size = max - min + Point::new(2.0 * metrics.padding, 2.0 * metrics.padding);
    format!(
        concat!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" class=\"moiety-structure\" role=\"img\" ",
            "viewBox=\"{x} {y} {w} {h}\" width=\"{w}\" height=\"{h}\">",
            "<title>{title}</title>",
            "<g class=\"bonds\" fill=\"none\" stroke=\"currentColor\" stroke-width=\"{line}\" ",
            "stroke-linecap=\"round\">{bonds}</g>",
            "<g class=\"atoms\" fill=\"currentColor\" font-family=\"Helvetica, Arial, sans-serif\" ",
            "font-size=\"{font}\">{atoms}</g>",
            "</svg>"
        ),
        x = num(origin.x),
        y = num(origin.y),
        w = num(size.x),
        h = num(size.y),
        title = escape(title),
        line = num(metrics.line),
        font = num(metrics.font),
        bonds = bonds,
        atoms = atoms,
    )
}

/// `value` with at most two decimals and no trailing zeros, as used in SVG attributes.
fn num(value: f64) -> String {
    if !value.is_finite() {
        return "0".to_string();
    }
    let text = format!("{value:.2}");
    let trimmed = text.trim_end_matches('0').trim_end_matches('.');
    if trimmed == "-0" {
        "0".to_string()
    } else {
        trimmed.to_string()
    }
}

#[cfg(test)]
mod tests;
