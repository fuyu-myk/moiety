//! Atom labels: an element symbol centred on the atom, its implicit hydrogens
//! on the side away from its bonds, and its charge as a superscript.

use moiety_core::{AtomKind, Element, Molecule};

use super::num;
use crate::graph::Adjacency;
use crate::layout::{Bounds, Point};
use crate::xml::escape;

/// Baseline of label text below the atom position, in font sizes, so capitals sit centred on it.
const BASELINE: f64 = 0.35;

/// Height of a capital letter, in font sizes.
const CAP_HEIGHT: f64 = 0.72;

/// Clearance kept between label text and the bonds that stop at it, in font sizes.
const MARGIN: f64 = 0.12;

/// Where a label puts the hydrogens of its atom.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Side {
    Right,
    Left,
    Above,
    Below,
}

impl Side {
    /// The side across the symbol from this one.
    fn opposite(self) -> Side {
        match self {
            Side::Right => Side::Left,
            Side::Left => Side::Right,
            Side::Above => Side::Below,
            Side::Below => Side::Above,
        }
    }
}

/// Vertical position of a piece of label text.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Level {
    Base,
    Subscript,
    Superscript,
}

impl Level {
    /// Baseline shift in font sizes, positive downward.
    fn shift(self) -> f64 {
        match self {
            Level::Base => 0.0,
            Level::Subscript => 0.3,
            Level::Superscript => -0.45,
        }
    }

    /// Font size relative to the label font.
    fn scale(self) -> f64 {
        match self {
            Level::Base => 1.0,
            Level::Subscript | Level::Superscript => 0.7,
        }
    }
}

/// Horizontal alignment of a run of text.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Anchor {
    Start,
    Middle,
    End,
}

/// A rectangle relative to the atom position, as `(lowest corner, highest corner)`.
pub(super) type Rect = (Point, Point);

/// One `<text>` element of a label, anchored on its baseline relative to the atom position.
struct Run {
    at: Point,
    anchor: Anchor,
    pieces: Vec<(String, Level)>,
}

impl Run {
    fn width(&self, font: f64) -> f64 {
        self.pieces
            .iter()
            .map(|(text, level)| text_width(text) * level.scale() * font)
            .sum()
    }

    /// The box the run's glyphs cover, grown by `margin` on every side.
    fn rect(&self, font: f64, margin: f64) -> Rect {
        let width = self.width(font);
        let left = match self.anchor {
            Anchor::Start => self.at.x,
            Anchor::Middle => self.at.x - width / 2.0,
            Anchor::End => self.at.x - width,
        };
        let (top, bottom) = self.pieces.iter().fold(
            (f64::INFINITY, f64::NEG_INFINITY),
            |(top, bottom), (_, level)| {
                let baseline = self.at.y + level.shift() * font;
                (
                    top.min(baseline - CAP_HEIGHT * level.scale() * font),
                    bottom.max(baseline),
                )
            },
        );
        (
            Point::new(left - margin, top - margin),
            Point::new(left + width + margin, bottom + margin),
        )
    }

    fn write(&self, centre: Point, font: f64, out: &mut String) {
        let anchor = match self.anchor {
            Anchor::Start => "start",
            Anchor::Middle => "middle",
            Anchor::End => "end",
        };
        let mut shift = 0.0;
        let mut body = String::new();
        for (text, level) in &self.pieces {
            let target = level.shift() * font;
            let dy = target - shift;
            shift = target;
            if *level == Level::Base && dy.abs() < 1e-9 {
                body.push_str(&escape(text));
                continue;
            }
            let size = match level {
                Level::Base => String::new(),
                _ => format!(" font-size=\"{}\"", num(level.scale() * font)),
            };
            body.push_str(&format!(
                "<tspan dy=\"{}\"{size}>{}</tspan>",
                num(dy),
                escape(text)
            ));
        }
        let at = centre + self.at;
        out.push_str(&format!(
            "<text x=\"{}\" y=\"{}\" text-anchor=\"{anchor}\">{body}</text>",
            num(at.x),
            num(at.y)
        ));
    }
}

/// The text of a labeled atom.
#[derive(Clone, Debug, PartialEq)]
pub(super) struct Label {
    symbol: String,
    hydrogens: u8,
    charge: i8,
    side: Side,
}

impl Label {
    /// The label of atom `index`, or `None` for a bare carbon vertex or a bonded attachment point.
    pub(super) fn for_atom(
        molecule: &Molecule,
        adjacency: &Adjacency,
        points: &[Point],
        index: usize,
    ) -> Option<Label> {
        let atom = molecule.atoms()[index];
        let bonded = !adjacency[index].is_empty();
        let symbol = match atom.kind {
            AtomKind::Attachment if bonded => return None,
            AtomKind::Attachment => "*",
            AtomKind::Element(Element::C) if bonded && atom.charge == 0 => return None,
            AtomKind::Element(element) => element.symbol(),
        };
        let side = match atom.kind {
            AtomKind::Element(element) if !bonded => free_side(element),
            _ => bonded_side(&directions(adjacency, points, index)),
        };
        Some(Label {
            symbol: symbol.to_string(),
            hydrogens: atom.hydrogens,
            charge: atom.charge,
            side,
        })
    }

    /// Boxes around every part of the label, relative to the atom, that bonds must stop short of.
    pub(super) fn keep_out(&self, font: f64) -> Vec<Rect> {
        self.runs(font)
            .iter()
            .map(|run| run.rect(font, MARGIN * font))
            .collect()
    }

    /// Radius of a circle around the element symbol.
    pub(super) fn radius(&self, font: f64) -> f64 {
        let (low, high) = self.runs(font)[0].rect(font, MARGIN * font);
        [low.x, low.y, high.x, high.y]
            .into_iter()
            .map(f64::abs)
            .fold(0.0, f64::max)
    }

    /// Writes the label text centred on `centre`.
    pub(super) fn render(&self, centre: Point, font: f64, out: &mut String, bounds: &mut Bounds) {
        for run in self.runs(font) {
            let (low, high) = run.rect(font, 0.0);
            bounds.include(centre + low);
            bounds.include(centre + high);
            run.write(centre, font, out);
        }
    }

    /// The text runs of the label, the element symbol first.
    fn runs(&self, font: f64) -> Vec<Run> {
        let baseline = BASELINE * font;
        let half = text_width(&self.symbol) * font / 2.0;
        let hydrogens = self.hydrogen_pieces();
        let charge: Vec<(String, Level)> = self.charge_piece().into_iter().collect();
        let run = |x: f64, y: f64, anchor: Anchor, pieces: Vec<(String, Level)>| Run {
            at: Point::new(x, y),
            anchor,
            pieces,
        };
        let mut runs = vec![run(
            0.0,
            baseline,
            Anchor::Middle,
            vec![(self.symbol.clone(), Level::Base)],
        )];
        match self.side {
            Side::Right => {
                let pieces = hydrogens.into_iter().chain(charge).collect();
                runs.push(run(half, baseline, Anchor::Start, pieces));
            }
            Side::Left => {
                runs.push(run(-half, baseline, Anchor::End, hydrogens));
                runs.push(run(half, baseline, Anchor::Start, charge));
            }
            Side::Above | Side::Below => {
                let step = if self.side == Side::Above {
                    -font
                } else {
                    font
                };
                runs.push(run(0.0, baseline + step, Anchor::Middle, hydrogens));
                runs.push(run(half, baseline, Anchor::Start, charge));
            }
        }
        runs.retain(|run| !run.pieces.is_empty());
        runs
    }

    fn hydrogen_pieces(&self) -> Vec<(String, Level)> {
        match self.hydrogens {
            0 => Vec::new(),
            1 => vec![("H".to_string(), Level::Base)],
            n => vec![
                ("H".to_string(), Level::Base),
                (n.to_string(), Level::Subscript),
            ],
        }
    }

    fn charge_piece(&self) -> Option<(String, Level)> {
        let sign = match self.charge {
            0 => return None,
            c if c > 0 => "+",
            _ => "\u{2212}",
        };
        let magnitude = self.charge.unsigned_abs();
        let text = if magnitude == 1 {
            sign.to_string()
        } else {
            format!("{magnitude}{sign}")
        };
        Some((text, Level::Superscript))
    }
}

/// Hydrogen side for an atom with no bonds: written first for the chalcogens
/// and halogens (`H2O`, `HCl`), last otherwise (`CH4`, `NH3`).
fn free_side(element: Element) -> Side {
    match element.atomic_number() {
        8 | 9 | 16 | 17 | 34 | 35 | 52 | 53 | 84 | 85 => Side::Left,
        _ => Side::Right,
    }
}

/// The label of every atom, as [`Label::for_atom`] gives it, except that an
/// atom whose hydrogens would cover another label moves them to the opposite
/// side, left for right or above for below, when its bonds leave that side
/// free and nothing is covered there; a charge stays where it is. Atoms are
/// settled in order, each against the labels settled before it and the
/// preferred labels of the rest.
pub(super) fn labels(
    molecule: &Molecule,
    adjacency: &Adjacency,
    points: &[Point],
    font: f64,
) -> Vec<Option<Label>> {
    let mut settled: Vec<Option<Label>> = (0..points.len())
        .map(|index| Label::for_atom(molecule, adjacency, points, index))
        .collect();
    let mut boxes = Boxes::new(&settled, points, font);
    for index in 0..points.len() {
        let Some(label) = settled[index].clone().filter(|l| l.hydrogens > 0) else {
            continue;
        };
        if adjacency[index].is_empty() {
            continue;
        }
        let clear = |side: Side| {
            let candidate = Label {
                side,
                ..label.clone()
            };
            placed(&candidate, points[index], font)
                .skip(1)
                .all(|rect| boxes.clear(points, index, rect))
        };
        let opposite = label.side.opposite();
        let free = free_sides(&directions(adjacency, points, index));
        if !clear(label.side) && free.contains(&opposite) && clear(opposite) {
            let moved = Label {
                side: opposite,
                ..label
            };
            boxes.rects[index] = placed(&moved, points[index], font).collect();
            settled[index] = Some(moved);
        }
    }
    settled
}

/// The boxes of the labels of a drawing, with the labelled atoms ordered by
/// `x` so that those near an atom are found without comparing every pair.
struct Boxes {
    rects: Vec<Vec<Rect>>,
    by_x: Vec<usize>,
    /// Labels whose atoms are further apart than this in `x` cannot touch.
    reach: f64,
}

impl Boxes {
    fn new(labels: &[Option<Label>], points: &[Point], font: f64) -> Self {
        let rects: Vec<Vec<Rect>> = labels
            .iter()
            .enumerate()
            .map(|(i, label)| {
                label
                    .as_ref()
                    .map_or_else(Vec::new, |l| placed(l, points[i], font).collect())
            })
            .collect();
        // No run reaches further from its atom than the whole label is wide.
        let reach = 2.0
            * rects
                .iter()
                .map(|boxes| boxes.iter().map(|(low, high)| high.x - low.x).sum::<f64>())
                .fold(0.0, f64::max);
        let mut by_x: Vec<usize> = (0..points.len())
            .filter(|&i| !rects[i].is_empty())
            .collect();
        by_x.sort_by(|&a, &b| points[a].x.total_cmp(&points[b].x));
        Self { rects, by_x, reach }
    }

    /// Whether `rect` overlaps no label but that of atom `index`.
    fn clear(&self, points: &[Point], index: usize, rect: Rect) -> bool {
        let x = points[index].x;
        let first = self.by_x.partition_point(|&i| points[i].x < x - self.reach);
        self.by_x[first..]
            .iter()
            .take_while(|&&i| points[i].x <= x + self.reach)
            .filter(|&&i| i != index)
            .flat_map(|&i| &self.rects[i])
            .all(|&other| !overlaps(rect, other))
    }
}

/// The boxes [`Label::keep_out`] gives for `label`, placed at `at`.
fn placed(label: &Label, at: Point, font: f64) -> impl Iterator<Item = Rect> {
    label
        .keep_out(font)
        .into_iter()
        .map(move |(low, high)| (at + low, at + high))
}

/// Unit vectors from atom `index` along each of its bonds.
fn directions(adjacency: &Adjacency, points: &[Point], index: usize) -> Vec<Point> {
    adjacency[index]
        .iter()
        .filter_map(|&(next, _)| (points[next] - points[index]).normalized())
        .collect()
}

fn overlaps(a: Rect, b: Rect) -> bool {
    a.0.x < b.1.x && b.0.x < a.1.x && a.0.y < b.1.y && b.0.y < a.1.y
}

/// Hydrogen side for an atom with bonds along the unit vectors `directions`:
/// the first of its free sides, or the right if none is free.
fn bonded_side(directions: &[Point]) -> Side {
    free_sides(directions)
        .first()
        .copied()
        .unwrap_or(Side::Right)
}

/// The sides no bond along the unit vectors `directions` points into: right,
/// left, above and below, in that order of preference.
fn free_sides(directions: &[Point]) -> Vec<Side> {
    const REACH: f64 = 0.3;
    let blocks = |side: Side, d: Point| match side {
        Side::Right => d.x > REACH,
        Side::Left => d.x < -REACH,
        Side::Above => d.y < -REACH,
        Side::Below => d.y > REACH,
    };
    [Side::Right, Side::Left, Side::Above, Side::Below]
        .into_iter()
        .filter(|&side| !directions.iter().any(|&d| blocks(side, d)))
        .collect()
}

/// Approximate advance width of `text` in a Helvetica-like face, in font sizes.
fn text_width(text: &str) -> f64 {
    text.chars()
        .map(|c| match c {
            'I' | 'i' | 'j' | 'l' => 0.28,
            'f' | 'r' | 't' => 0.33,
            'M' | 'W' | 'm' | 'w' => 0.85,
            'A'..='Z' => 0.7,
            'a'..='z' => 0.55,
            _ => 0.58,
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn label(hydrogens: u8, charge: i8, side: Side) -> Label {
        Label {
            symbol: "N".to_string(),
            hydrogens,
            charge,
            side,
        }
    }

    #[test]
    fn hydrogens_go_where_no_bond_points() {
        let from = |angles: &[f64]| {
            let directions: Vec<Point> = angles
                .iter()
                .map(|a| Point::from_angle(a.to_radians()))
                .collect();
            bonded_side(&directions)
        };
        assert_eq!(from(&[180.0]), Side::Right);
        assert_eq!(from(&[-30.0]), Side::Left);
        assert_eq!(from(&[30.0, 150.0]), Side::Above);
        assert_eq!(from(&[-30.0, -150.0]), Side::Below);
        assert_eq!(from(&[-90.0]), Side::Right);
        assert_eq!(free_side(Element::O), Side::Left);
        assert_eq!(free_side(Element::N), Side::Right);
    }

    #[test]
    fn moves_hydrogens_off_a_neighbouring_label() {
        let molecule = moiety_core::smiles::parse("OCCO").unwrap();
        let adjacency = molecule.adjacency();
        let sides = |points: &[Point]| -> Vec<Option<Side>> {
            labels(&molecule, &adjacency, points, 16.0)
                .iter()
                .map(|label| label.as_ref().map(|l| l.side))
                .collect()
        };
        let close = [
            Point::new(0.0, 0.0),
            Point::new(0.0, -40.0),
            Point::new(70.0, -40.0),
            Point::new(36.0, 0.0),
        ];
        assert_eq!(
            sides(&close),
            [Some(Side::Left), None, None, Some(Side::Left)]
        );
        let apart = [
            Point::new(0.0, 0.0),
            Point::new(0.0, -40.0),
            Point::new(90.0, -40.0),
            Point::new(56.0, 0.0),
        ];
        assert_eq!(
            sides(&apart),
            [Some(Side::Right), None, None, Some(Side::Left)]
        );
    }

    #[test]
    fn writes_charges_as_superscripts() {
        assert_eq!(label(0, 0, Side::Right).charge_piece(), None);
        assert_eq!(
            label(0, -1, Side::Right).charge_piece().unwrap().0,
            "\u{2212}"
        );
        assert_eq!(label(0, 2, Side::Right).charge_piece().unwrap().0, "2+");
        let mut out = String::new();
        label(3, 1, Side::Right).render(Point::default(), 10.0, &mut out, &mut Bounds::default());
        assert_eq!(
            out,
            concat!(
                "<text x=\"0\" y=\"3.5\" text-anchor=\"middle\">N</text>",
                "<text x=\"3.5\" y=\"3.5\" text-anchor=\"start\">H<tspan dy=\"3\" font-size=\"7\">3</tspan>",
                "<tspan dy=\"-7.5\" font-size=\"7\">+</tspan></text>"
            )
        );
    }

    #[test]
    fn keeps_bonds_clear_of_every_run() {
        let plain = label(0, 0, Side::Right).keep_out(10.0);
        assert_eq!(plain.len(), 1);
        let charged = label(2, 1, Side::Left).keep_out(10.0);
        assert_eq!(charged.len(), 3);
        let (low, high) = charged[2];
        assert!(
            low.x > 0.0 && high.y < plain[0].1.y,
            "the charge sits upper right"
        );
        assert!((label(0, 0, Side::Right).radius(10.0) - 4.9).abs() < 1e-9);
    }
}
