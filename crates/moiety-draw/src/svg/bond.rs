//! Bond groups: one, two or three strokes per bond, stopped short of labels.

use moiety_core::{Bond, BondOrder};

use super::label::Rect;
use super::{Scene, num};
use crate::layout::{Bounds, EPSILON, Point};
use crate::xml::escape;

/// Distance between the strokes of a double bond, in bond lengths.
const DOUBLE_GAP: f64 = 0.18;

/// Distance between neighbouring strokes of a triple bond, in bond lengths.
const TRIPLE_GAP: f64 = 0.15;

/// How much shorter the inner stroke of a double bond is at each bare vertex, in bond lengths.
const INNER_INSET: f64 = 0.14;

/// Longest a mitred stroke may extend or shrink at a vertex, in bond lengths.
const MAX_MITER: f64 = 0.5;

/// Strokes shorter than this, in pixels, are left out.
const MIN_STROKE: f64 = 1.0;

/// Radius of the circle drawn in an aromatic ring that has no Kekulé structure, relative to the ring's apothem.
const CIRCLE_RATIO: f64 = 0.62;

/// One stroke of a bond: the line between the atom centres moved by `offset`,
/// cut back by `inset` pixels at each end that has no label.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Stroke {
    offset: Point,
    inset: f64,
}

impl Stroke {
    const CENTRE: Stroke = Stroke {
        offset: Point::new(0.0, 0.0),
        inset: 0.0,
    };

    fn shifted(offset: Point) -> Stroke {
        Stroke { offset, inset: 0.0 }
    }
}

/// Writes the group of bond `index`: a wide invisible hit line, its strokes and a title.
pub(super) fn render(
    scene: &Scene,
    index: usize,
    span: &str,
    out: &mut String,
    bounds: &mut Bounds,
) {
    let bond = scene.molecule.bonds()[index];
    let (from, to) = (scene.points[bond.a], scene.points[bond.b]);
    let order = order_name(bond.order);
    let title = format!("bond {index}: {order}, atoms {} and {}", bond.a, bond.b);
    out.push_str(&format!(
        "<g class=\"bond {order}\" data-bond=\"{index}\"{span}><title>{}</title>",
        escape(&title)
    ));
    if let Some(along) = (to - from).normalized() {
        let hit = (4.0 * scene.metrics.line).max(8.0);
        out.push_str(&format!(
            "<line class=\"hit\" {} stroke-width=\"{}\" stroke-opacity=\"0\"/>",
            ends(from, to),
            num(hit)
        ));
        for (start, end) in stroke_ends(scene, &bond, along, &strokes(scene, &bond, along)) {
            bounds.include(start);
            bounds.include(end);
            out.push_str(&format!("<line {}/>", ends(start, end)));
        }
    }
    out.push_str("</g>");
}

/// Writes a circle inside every ring that is aromatic but has no Kekulé structure.
pub(super) fn render_circles(scene: &Scene, out: &mut String, bounds: &mut Bounds) {
    for &ring_index in &scene.circled_rings {
        let ring = &scene.rings[ring_index];
        let centre = scene.ring_centres[ring_index];
        let apothem = ring
            .iter()
            .enumerate()
            .map(|(i, &atom)| {
                let next = ring[(i + 1) % ring.len()];
                ((scene.points[atom] + scene.points[next]) * 0.5).distance(centre)
            })
            .sum::<f64>()
            / ring.len() as f64;
        let radius = CIRCLE_RATIO * apothem;
        bounds.include_box(centre, radius, radius);
        out.push_str(&format!(
            "<circle class=\"aromatic-ring\" cx=\"{}\" cy=\"{}\" r=\"{}\"/>",
            num(centre.x),
            num(centre.y),
            num(radius)
        ));
    }
}

fn order_name(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Single => "single",
        BondOrder::Double => "double",
        BondOrder::Triple => "triple",
        BondOrder::Aromatic => "aromatic",
    }
}

fn ends(start: Point, end: Point) -> String {
    format!(
        "x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\"",
        num(start.x),
        num(start.y),
        num(end.x),
        num(end.y)
    )
}

fn strokes(scene: &Scene, bond: &Bond, along: Point) -> Vec<Stroke> {
    let length = scene.metrics.bond;
    let across = along.perpendicular();
    match bond.order {
        BondOrder::Single | BondOrder::Aromatic => vec![Stroke::CENTRE],
        BondOrder::Triple => {
            let gap = across * (TRIPLE_GAP * length);
            vec![Stroke::CENTRE, Stroke::shifted(gap), Stroke::shifted(-gap)]
        }
        BondOrder::Double => match inner_side(scene, bond, along) {
            Some(side) => vec![
                Stroke::CENTRE,
                Stroke {
                    offset: side * (DOUBLE_GAP * length),
                    inset: INNER_INSET * length,
                },
            ],
            None => {
                let half = across * (DOUBLE_GAP * length / 2.0);
                vec![Stroke::shifted(half), Stroke::shifted(-half)]
            }
        },
    }
}

/// The unit normal toward the side that the second stroke of a double bond
/// goes on: the ring centre, or the side with more neighbours. `None` centres
/// the pair, as for a terminal `C=O`, `CH2=CH2`, or a bond with balanced
/// neighbours. A bare terminal carbon, as in `CH2=CH-`, keeps the main stroke
/// on the axis so the next bond joins it.
fn inner_side(scene: &Scene, bond: &Bond, along: Point) -> Option<Point> {
    let across = along.perpendicular();
    let toward = |side: f64| if side >= 0.0 { across } else { -across };
    if let Some(centre) = ring_centre(scene, bond.a, bond.b) {
        return Some(toward((centre - scene.points[bond.a]).dot(across)));
    }
    let terminal = |atom: usize| scene.adjacency[atom].len() == 1;
    let labelled_terminal = |atom: usize| terminal(atom) && scene.labels[atom].is_some();
    if (terminal(bond.a) && terminal(bond.b))
        || labelled_terminal(bond.a)
        || labelled_terminal(bond.b)
    {
        return None;
    }
    let balance: i32 = [(bond.a, bond.b), (bond.b, bond.a)]
        .into_iter()
        .flat_map(|(atom, other)| {
            scene.adjacency[atom]
                .iter()
                .filter(move |&&(next, _)| next != other)
                .map(move |&(next, _)| along.cross(scene.points[next] - scene.points[atom]))
        })
        .map(|side| i32::from(side > EPSILON) - i32::from(side < -EPSILON))
        .sum();
    (balance != 0).then(|| toward(f64::from(balance)))
}

/// Centre of the smallest ring in which `a` and `b` are neighbours.
fn ring_centre(scene: &Scene, a: usize, b: usize) -> Option<Point> {
    scene
        .rings
        .iter()
        .enumerate()
        .filter(|(_, ring)| {
            ring.iter().enumerate().any(|(i, &atom)| {
                let next = ring[(i + 1) % ring.len()];
                (atom, next) == (a, b) || (atom, next) == (b, a)
            })
        })
        .min_by_key(|(i, ring)| (ring.len(), *i))
        .map(|(i, _)| scene.ring_centres[i])
}

/// The visible ends of every stroke that is long enough to draw.
fn stroke_ends(
    scene: &Scene,
    bond: &Bond,
    along: Point,
    strokes: &[Stroke],
) -> Vec<(Point, Point)> {
    let (from, to) = (scene.points[bond.a], scene.points[bond.b]);
    let starts = end_cuts(scene, bond.a, bond.b, along, strokes);
    let finishes = end_cuts(scene, bond.b, bond.a, -along, strokes);
    strokes
        .iter()
        .zip(starts.into_iter().zip(finishes))
        .filter(|&(_, (start, finish))| from.distance(to) - start - finish >= MIN_STROKE)
        .map(|(stroke, (start, finish))| {
            (
                from + stroke.offset + along * start,
                to + stroke.offset - along * finish,
            )
        })
        .collect()
}

/// How far each stroke is cut back from `atom`, measured along `direction` away
/// from it. At a label every stroke takes the same cut, the largest any of them
/// needs, so parallel strokes end level whatever the bond's angle. At a bare
/// vertex an inner stroke keeps its inset and an offset stroke is mitred to the
/// neighbouring bond on its side, so the strokes meet that bond exactly.
fn end_cuts(
    scene: &Scene,
    atom: usize,
    other: usize,
    direction: Point,
    strokes: &[Stroke],
) -> Vec<f64> {
    if let Some(label) = &scene.labels[atom] {
        let rects = label.keep_out(scene.metrics.font);
        let cut = strokes
            .iter()
            .flat_map(|stroke| {
                rects
                    .iter()
                    .filter_map(move |&rect| exit(stroke.offset, direction, rect))
            })
            .fold(0.0, f64::max);
        return vec![cut; strokes.len()];
    }
    strokes
        .iter()
        .map(|stroke| {
            if stroke.inset > 0.0 {
                stroke.inset
            } else {
                miter(scene, atom, other, direction, stroke.offset)
            }
        })
        .collect()
}

/// Where a stroke offset from the bond axis meets the nearest neighbouring bond
/// on its own side, as a distance along `direction`. Negative values extend the
/// stroke past the atom. Zero when no neighbour lies on that side or the stroke
/// runs along the axis.
fn miter(scene: &Scene, atom: usize, other: usize, direction: Point, offset: Point) -> f64 {
    let side = direction.cross(offset);
    if side.abs() < EPSILON {
        return 0.0;
    }
    let limit = MAX_MITER * scene.metrics.bond;
    scene.adjacency[atom]
        .iter()
        .filter(|&&(next, _)| next != other)
        .filter_map(|&(next, _)| (scene.points[next] - scene.points[atom]).normalized())
        .filter(|&toward| direction.cross(toward) * side > EPSILON)
        .map(|toward| -offset.cross(toward) / direction.cross(toward))
        .filter(|cut| cut.abs() <= limit)
        .min_by(|a, b| a.abs().total_cmp(&b.abs()))
        .unwrap_or(0.0)
}

/// Distance along the ray from `origin` in `direction` to where it leaves
/// `rect`, if the ray passes through the rectangle ahead of its origin.
fn exit(origin: Point, direction: Point, (low, high): Rect) -> Option<f64> {
    let mut enter = f64::NEG_INFINITY;
    let mut leave = f64::INFINITY;
    for (start, step, min, max) in [
        (origin.x, direction.x, low.x, high.x),
        (origin.y, direction.y, low.y, high.y),
    ] {
        if step.abs() < EPSILON {
            if start < min || start > max {
                return None;
            }
        } else {
            let (a, b) = ((min - start) / step, (max - start) / step);
            enter = enter.max(a.min(b));
            leave = leave.min(a.max(b));
        }
    }
    (enter <= leave && leave > 0.0).then_some(leave)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn leaves_a_label_box_along_the_bond() {
        let rect = (Point::new(-10.0, -5.0), Point::new(10.0, 5.0));
        let right = Point::new(1.0, 0.0);
        assert!((exit(Point::default(), right, rect).unwrap() - 10.0).abs() < 1e-9);
        assert!((exit(Point::default(), Point::new(0.0, -1.0), rect).unwrap() - 5.0).abs() < 1e-9);
        assert!((exit(Point::new(0.0, 4.0), right, rect).unwrap() - 10.0).abs() < 1e-9);
        assert_eq!(exit(Point::new(20.0, 0.0), right, rect), None);
        assert_eq!(exit(Point::new(0.0, 8.0), right, rect), None);
        let ahead = (Point::new(12.0, -1.0), Point::new(16.0, 1.0));
        assert!((exit(Point::default(), right, ahead).unwrap() - 16.0).abs() < 1e-9);
    }
}
