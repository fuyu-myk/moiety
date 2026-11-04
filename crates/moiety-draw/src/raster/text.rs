//! Label text as glyph outlines, set the way an SVG renderer sets [`Text`].

use super::glyphs::{GLYPHS, UNITS_PER_EM};
use crate::layout::Point;
use crate::svg::shape::{Anchor, Text};

/// One step of a glyph outline, in font units with `y` growing upward.
#[derive(Clone, Copy, Debug)]
pub(super) enum Step {
    /// Starts a contour.
    M(i16, i16),
    /// A straight edge to a point.
    L(i16, i16),
    /// A quadratic curve through a control point to an end point.
    Q(i16, i16, i16, i16),
}

/// Advance width, in ems, of a character the glyph table lacks.
const MISSING_ADVANCE: f64 = 0.58;

/// Straight pieces each quadratic curve is split into.
const CURVE_PIECES: usize = 6;

/// The advance width and outline of `character`.
fn glyph(character: char) -> Option<(u16, &'static [Step])> {
    GLYPHS
        .binary_search_by(|(c, _, _)| c.cmp(&character))
        .ok()
        .map(|index| (GLYPHS[index].1, GLYPHS[index].2))
}

/// Advance width of `text` set at `size` pixels.
fn advance(text: &str, size: f64) -> f64 {
    text.chars()
        .map(|c| glyph(c).map_or(MISSING_ADVANCE, |(a, _)| f64::from(a) / UNITS_PER_EM))
        .sum::<f64>()
        * size
}

/// The outlines of every glyph of `text`, in pixels: one list of closed
/// contours per glyph, holes included.
pub(super) fn outlines(text: &Text) -> Vec<Vec<Vec<Point>>> {
    let width: f64 = text
        .pieces
        .iter()
        .map(|(piece, level)| advance(piece, level.scale() * text.font))
        .sum();
    let mut pen = text.at.x
        - match text.anchor {
            Anchor::Start => 0.0,
            Anchor::Middle => width / 2.0,
            Anchor::End => width,
        };
    let mut glyphs = Vec::new();
    for (piece, level) in &text.pieces {
        let size = level.scale() * text.font;
        let baseline = text.at.y + level.shift() * text.font;
        for character in piece.chars() {
            let Some((advance, steps)) = glyph(character) else {
                pen += MISSING_ADVANCE * size;
                continue;
            };
            let scale = size / UNITS_PER_EM;
            let place = |x: i16, y: i16| {
                Point::new(pen + f64::from(x) * scale, baseline - f64::from(y) * scale)
            };
            glyphs.push(contours(steps, place));
            pen += f64::from(advance) * scale;
        }
    }
    glyphs
}

/// The contours of an outline as polygons, its points placed by `place`.
fn contours(steps: &[Step], place: impl Fn(i16, i16) -> Point) -> Vec<Vec<Point>> {
    let mut contours: Vec<Vec<Point>> = Vec::new();
    for step in steps {
        match *step {
            Step::M(x, y) => contours.push(vec![place(x, y)]),
            Step::L(x, y) => {
                if let Some(contour) = contours.last_mut() {
                    contour.push(place(x, y));
                }
            }
            Step::Q(cx, cy, x, y) => {
                if let Some(contour) = contours.last_mut() {
                    let from = contour.last().copied().unwrap_or_default();
                    let (control, to) = (place(cx, cy), place(x, y));
                    contour.extend(curve(from, control, to));
                }
            }
        }
    }
    contours
}

/// Points along the quadratic curve from `from` through `control` to `to`,
/// after `from`.
pub(super) fn curve(from: Point, control: Point, to: Point) -> impl Iterator<Item = Point> {
    (1..=CURVE_PIECES).map(move |k| {
        let t = k as f64 / CURVE_PIECES as f64;
        let u = 1.0 - t;
        from * (u * u) + control * (2.0 * u * t) + to * (t * t)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::svg::shape::Level;

    fn text(anchor: Anchor, pieces: &[(&str, Level)]) -> Text {
        Text {
            at: Point::new(100.0, 50.0),
            anchor,
            pieces: pieces.iter().map(|&(t, l)| (t.to_string(), l)).collect(),
            font: 16.0,
        }
    }

    fn extent(glyphs: &[Vec<Vec<Point>>]) -> (Point, Point) {
        glyphs.iter().flatten().flatten().fold(
            (
                Point::new(f64::MAX, f64::MAX),
                Point::new(f64::MIN, f64::MIN),
            ),
            |(low, high), p| {
                (
                    Point::new(low.x.min(p.x), low.y.min(p.y)),
                    Point::new(high.x.max(p.x), high.y.max(p.y)),
                )
            },
        )
    }

    #[test]
    fn holds_every_label_character_in_order() {
        assert!(GLYPHS.windows(2).all(|pair| pair[0].0 < pair[1].0));
        let needed = ('0'..='9')
            .chain('A'..='Z')
            .chain('a'..='z')
            .chain(['+', '*', '\u{2212}']);
        assert!(needed.into_iter().all(|c| glyph(c).is_some()));
        assert!(glyph('é').is_none());
    }

    #[test]
    fn sets_capitals_on_the_baseline_at_cap_height() {
        let glyphs = outlines(&text(Anchor::Middle, &[("N", Level::Base)]));
        assert_eq!(glyphs.len(), 1);
        let (low, high) = extent(&glyphs);
        assert!((high.y - 50.0).abs() < 0.01, "baseline {}", high.y);
        let cap = 50.0 - low.y;
        assert!((cap / 16.0 - 0.727).abs() < 0.01, "cap height {cap}");
        assert!(
            (low.x + high.x - 200.0).abs() < 1.0,
            "centred {} {}",
            low.x,
            high.x
        );
    }

    #[test]
    fn anchors_runs_and_shrinks_raised_pieces() {
        let start = extent(&outlines(&text(Anchor::Start, &[("H", Level::Base)])));
        let end = extent(&outlines(&text(Anchor::End, &[("H", Level::Base)])));
        assert!(start.0.x > 100.0 && end.1.x < 100.0);
        let charged = outlines(&text(
            Anchor::Start,
            &[
                ("H", Level::Base),
                ("2", Level::Subscript),
                ("+", Level::Superscript),
            ],
        ));
        assert_eq!(charged.len(), 3);
        let (_, two) = extent(&charged[1..2]);
        let (plus, _) = extent(&charged[2..3]);
        assert!(two.y > 50.0, "subscript sits below the baseline");
        assert!(plus.y < 50.0 - 0.6 * 16.0, "superscript is raised");
        assert!(
            extent(&outlines(&text(Anchor::Start, &[("é", Level::Base)])))
                .0
                .x
                > 1e9
        );
    }
}
