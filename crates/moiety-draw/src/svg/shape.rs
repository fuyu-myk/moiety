//! The marks a drawing is made of, in pixel coordinates, and how each is
//! written as SVG.

use super::num;
use crate::layout::Point;
use crate::xml::escape;

/// One mark of a drawing.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Shape {
    /// A stroked straight line.
    Line(Point, Point),
    /// A filled triangle: the solid wedge of a stereo bond.
    Wedge([Point; 3]),
    /// Stroked lines across a bond, widening away from its stereocentre:
    /// the hashed wedge of a stereo bond.
    Hash(Vec<(Point, Point)>),
    /// A stroked circle: the ring drawn inside an aromatic ring.
    Ring { centre: Point, radius: f64 },
    /// A stroked path of quadratic curves from `start`, each given as
    /// `(control, end)`: the wavy mark of an attachment point.
    Wave {
        start: Point,
        curves: Vec<(Point, Point)>,
    },
    /// The label of the atom at `centre`, as runs of text, the element symbol first.
    Label { centre: Point, runs: Vec<Text> },
}

impl Shape {
    /// Writes the shape as an SVG element, stroking with lines `line_width` wide.
    pub(crate) fn write(&self, line_width: f64, out: &mut String) {
        match self {
            Shape::Line(start, end) => out.push_str(&format!("<line {}/>", ends(*start, *end))),
            Shape::Hash(lines) => {
                for (start, end) in lines {
                    out.push_str(&format!("<line {}/>", ends(*start, *end)));
                }
            }
            Shape::Wedge(corners) => {
                let corners = corners
                    .iter()
                    .map(|p| format!("{},{}", num(p.x), num(p.y)))
                    .collect::<Vec<_>>()
                    .join(" ");
                out.push_str(&format!(
                    "<polygon points=\"{corners}\" fill=\"currentColor\" stroke-linejoin=\"round\"/>"
                ));
            }
            Shape::Ring { centre, radius } => out.push_str(&format!(
                "<circle class=\"aromatic-ring\" cx=\"{}\" cy=\"{}\" r=\"{}\"/>",
                num(centre.x),
                num(centre.y),
                num(*radius)
            )),
            Shape::Wave { start, curves } => {
                let mut path = format!("M{} {}", num(start.x), num(start.y));
                for (control, end) in curves {
                    path.push_str(&format!(
                        " Q{} {} {} {}",
                        num(control.x),
                        num(control.y),
                        num(end.x),
                        num(end.y)
                    ));
                }
                out.push_str(&format!(
                    "<path class=\"attachment\" d=\"{path}\" fill=\"none\" stroke=\"currentColor\" \
                     stroke-width=\"{}\" stroke-linecap=\"round\"/>",
                    num(line_width)
                ));
            }
            Shape::Label { runs, .. } => runs.iter().for_each(|run| run.write(out)),
        }
    }
}

/// The `x1`, `y1`, `x2` and `y2` attributes of a line from `start` to `end`.
pub(crate) fn ends(start: Point, end: Point) -> String {
    format!(
        "x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\"",
        num(start.x),
        num(start.y),
        num(end.x),
        num(end.y)
    )
}

/// Vertical position of a piece of label text.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Level {
    Base,
    Subscript,
    Superscript,
}

impl Level {
    /// Baseline shift in font sizes, positive downward.
    pub(crate) fn shift(self) -> f64 {
        match self {
            Level::Base => 0.0,
            Level::Subscript => 0.3,
            Level::Superscript => -0.45,
        }
    }

    /// Font size relative to the label font.
    pub(crate) fn scale(self) -> f64 {
        match self {
            Level::Base => 1.0,
            Level::Subscript | Level::Superscript => 0.7,
        }
    }
}

/// Horizontal alignment of a run of text.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Anchor {
    Start,
    Middle,
    End,
}

/// A run of label text on one baseline, its pieces raised or lowered by their level.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Text {
    /// Where the baseline meets the anchor.
    pub(crate) at: Point,
    /// Which end of the run, or its middle, sits at `at`.
    pub(crate) anchor: Anchor,
    /// The text of the run, piece by piece.
    pub(crate) pieces: Vec<(String, Level)>,
    /// Font size of base text, in pixels.
    pub(crate) font: f64,
}

impl Text {
    fn write(&self, out: &mut String) {
        let anchor = match self.anchor {
            Anchor::Start => "start",
            Anchor::Middle => "middle",
            Anchor::End => "end",
        };
        let mut shift = 0.0;
        let mut body = String::new();
        for (text, level) in &self.pieces {
            let target = level.shift() * self.font;
            let dy = target - shift;
            shift = target;
            if *level == Level::Base && dy.abs() < 1e-9 {
                body.push_str(&escape(text));
                continue;
            }
            let size = match level {
                Level::Base => String::new(),
                _ => format!(" font-size=\"{}\"", num(level.scale() * self.font)),
            };
            body.push_str(&format!(
                "<tspan dy=\"{}\"{size}>{}</tspan>",
                num(dy),
                escape(text)
            ));
        }
        out.push_str(&format!(
            "<text x=\"{}\" y=\"{}\" text-anchor=\"{anchor}\">{body}</text>",
            num(self.at.x),
            num(self.at.y)
        ));
    }
}
