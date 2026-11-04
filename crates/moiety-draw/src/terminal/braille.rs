//! Drawings in Unicode text: bonds as braille dots, two across and four down
//! in each character cell, and labels as the characters themselves.

use super::{CELL_HEIGHT, CELL_WIDTH, TEXT_ZOOM};
use crate::layout::Point;
use crate::svg::Picture;
use crate::svg::shape::{Anchor, Level, Shape, Text};

/// Dots across and down one character cell.
const DOTS_ACROSS: usize = 2;
const DOTS_DOWN: usize = 4;

/// The bit of each dot of a braille character, by column and then row.
const BITS: [[u32; DOTS_DOWN]; DOTS_ACROSS] = [[0x01, 0x02, 0x04, 0x40], [0x08, 0x10, 0x20, 0x80]];

/// `picture` as lines of text at most `columns` characters wide.
pub(super) fn render(picture: &Picture, columns: usize) -> String {
    let budget = columns.max(1);
    let mut fit = Fit::new(picture, budget);
    if fit.columns > budget {
        let overhang = fit.columns - fit.drawing;
        fit = Fit::new(picture, budget.saturating_sub(overhang).max(1));
    }
    let mut grid = Grid::new(fit.columns.min(budget), fit.rows);
    for shape in &picture.shapes {
        grid.shape(shape, &fit.frame);
    }
    for label in &fit.labels {
        grid.label(label);
    }
    grid.text()
}

/// A drawing fitted to a grid: its frame, its labels placed in cells, and
/// the columns and rows the grid needs to hold both.
struct Fit {
    frame: Frame,
    labels: Vec<Vec<Placed>>,
    /// Columns the bonds span.
    drawing: usize,
    columns: usize,
    rows: usize,
}

impl Fit {
    /// The fit that makes the bonds of `picture` at most `columns` wide.
    fn new(picture: &Picture, columns: usize) -> Self {
        let (min, max) = picture.bounds.extent().unwrap_or_default();
        let margin = Point::new(picture.padding, picture.padding);
        let size = max - min + margin * 2.0;
        let zoom = (columns as f64 * CELL_WIDTH / size.x.max(1.0)).min(TEXT_ZOOM);
        let cell = Point::new(CELL_WIDTH / zoom, CELL_HEIGHT / zoom);
        let mut frame = Frame {
            origin: min - margin,
            cell,
        };
        let mut labels: Vec<Vec<Placed>> = picture
            .shapes
            .iter()
            .filter_map(|shape| match shape {
                Shape::Label { centre, runs } => Some(placements(frame.cell(*centre), runs)),
                _ => None,
            })
            .collect();
        let drawing = (
            (size.x / cell.x).ceil().max(1.0),
            (size.y / cell.y).ceil().max(1.0),
        );
        let (low, high) = labels.iter().flatten().fold(
            (Point::default(), Point::new(drawing.0, drawing.1)),
            |(low, high), placed| {
                let end = placed.column + placed.text.len() as f64;
                (
                    Point::new(low.x.min(placed.column), low.y.min(placed.row)),
                    Point::new(high.x.max(end), high.y.max(placed.row + 1.0)),
                )
            },
        );
        frame.origin = frame.origin + Point::new(low.x * cell.x, low.y * cell.y);
        for placed in labels.iter_mut().flatten() {
            placed.column -= low.x;
            placed.row -= low.y;
        }
        Self {
            frame,
            labels,
            drawing: drawing.0 as usize,
            columns: (high.x - low.x) as usize,
            rows: (high.y - low.y) as usize,
        }
    }
}

/// A run of label characters placed on the grid, from its first cell.
struct Placed {
    column: f64,
    row: f64,
    text: Vec<char>,
}

/// The runs of a label whose atom lies in the cell `(column, row)`: the
/// symbol on that cell, hydrogens and charges in the cells beside, above or
/// below it.
fn placements((column, row): (f64, f64), runs: &[Text]) -> Vec<Placed> {
    let Some((symbol, rest)) = runs.split_first() else {
        return Vec::new();
    };
    let text = letters(symbol);
    let left = column - ((text.len().max(1) - 1) / 2) as f64;
    let mut right = left + text.len() as f64;
    let mut placed = vec![Placed {
        column: left,
        row,
        text,
    }];
    for run in rest {
        let text = letters(run);
        let width = text.len() as f64;
        let shift = run.at.y - symbol.at.y;
        let line = if shift.abs() < 1e-6 {
            row
        } else {
            row + shift.signum()
        };
        let start = match run.anchor {
            Anchor::Start if line == row => right,
            Anchor::Start => left,
            Anchor::End => left - width,
            Anchor::Middle => column - ((text.len().max(1) - 1) / 2) as f64,
        };
        if line == row && run.anchor == Anchor::Start {
            right += width;
        }
        placed.push(Placed {
            column: start,
            row: line,
            text,
        });
    }
    placed
}

/// Where a drawing lies on the grid: the drawing point at its top left
/// corner, and the size of a cell in drawing pixels.
struct Frame {
    origin: Point,
    cell: Point,
}

impl Frame {
    /// The point `p` of the drawing, in dots from the top left of the grid.
    fn dots(&self, p: Point) -> Point {
        let offset = p - self.origin;
        Point::new(
            offset.x / self.cell.x * DOTS_ACROSS as f64,
            offset.y / self.cell.y * DOTS_DOWN as f64,
        )
    }

    /// The column and row of the cell that holds the point `p` of the drawing.
    fn cell(&self, p: Point) -> (f64, f64) {
        let offset = p - self.origin;
        (
            (offset.x / self.cell.x).floor(),
            (offset.y / self.cell.y).floor(),
        )
    }
}

/// Character cells, each with its braille dots or a label character.
struct Grid {
    columns: usize,
    rows: usize,
    dots: Vec<u32>,
    letters: Vec<Option<char>>,
}

impl Grid {
    fn new(columns: usize, rows: usize) -> Self {
        Self {
            columns,
            rows,
            dots: vec![0; columns * rows],
            letters: vec![None; columns * rows],
        }
    }

    /// Draws the dots of `shape`, placed on the grid by `frame`. Labels are
    /// written by [`Grid::label`] instead.
    fn shape(&mut self, shape: &Shape, frame: &Frame) {
        match shape {
            Shape::Line(a, b) => self.line(frame.dots(*a), frame.dots(*b)),
            Shape::Hash(lines) => {
                for &(a, b) in lines {
                    self.dot(frame.dots((a + b) * 0.5));
                }
            }
            Shape::Wedge(corners) => self.triangle(corners.map(|p| frame.dots(p))),
            Shape::Ring { centre, radius } => {
                let steps = (radius / frame.cell.x * 8.0).ceil().max(12.0) as usize;
                let points: Vec<Point> = (0..=steps)
                    .map(|k| {
                        let angle = std::f64::consts::TAU * k as f64 / steps as f64;
                        frame.dots(*centre + Point::from_angle(angle) * *radius)
                    })
                    .collect();
                points.windows(2).for_each(|w| self.line(w[0], w[1]));
            }
            Shape::Wave { start, curves } => {
                let mut from = frame.dots(*start);
                for &(_, end) in curves {
                    let to = frame.dots(end);
                    self.line(from, to);
                    from = to;
                }
            }
            Shape::Label { .. } => {}
        }
    }

    /// Sets the dot at `p`, in dots from the top left, if it is on the grid.
    fn dot(&mut self, p: Point) {
        if p.x < 0.0 || p.y < 0.0 {
            return;
        }
        let (x, y) = (p.x as usize, p.y as usize);
        let (column, row) = (x / DOTS_ACROSS, y / DOTS_DOWN);
        if column < self.columns && row < self.rows {
            self.dots[row * self.columns + column] |= BITS[x % DOTS_ACROSS][y % DOTS_DOWN];
        }
    }

    fn line(&mut self, a: Point, b: Point) {
        let steps = ((b - a).length() * 2.0).ceil().max(1.0) as usize;
        for k in 0..=steps {
            self.dot(a + (b - a) * (k as f64 / steps as f64));
        }
    }

    fn triangle(&mut self, corners: [Point; 3]) {
        let low = corners.iter().fold(Point::new(f64::MAX, f64::MAX), |l, p| {
            Point::new(l.x.min(p.x), l.y.min(p.y))
        });
        let high = corners.iter().fold(Point::new(f64::MIN, f64::MIN), |h, p| {
            Point::new(h.x.max(p.x), h.y.max(p.y))
        });
        let side = |a: Point, b: Point, p: Point| (b - a).cross(p - a);
        for y in low.y.floor().max(0.0) as usize..=high.y.ceil().max(0.0) as usize {
            for x in low.x.floor().max(0.0) as usize..=high.x.ceil().max(0.0) as usize {
                let p = Point::new(x as f64 + 0.5, y as f64 + 0.5);
                let signs = [0, 1, 2].map(|i| side(corners[i], corners[(i + 1) % 3], p));
                if signs.iter().all(|&s| s >= 0.0) || signs.iter().all(|&s| s <= 0.0) {
                    self.dot(p);
                }
            }
        }
        for i in 0..3 {
            self.line(corners[i], corners[(i + 1) % 3]);
        }
    }

    /// Writes the runs of a label, moved across together as far as it takes
    /// to keep them on the grid.
    fn label(&mut self, runs: &[Placed]) {
        let low = runs.iter().map(|p| p.column).fold(f64::INFINITY, f64::min);
        let high = runs
            .iter()
            .map(|p| p.column + p.text.len() as f64)
            .fold(f64::NEG_INFINITY, f64::max);
        let over = (high - self.columns as f64).max(0.0);
        let shift = if low - over < 0.0 { -low } else { -over };
        for run in runs {
            self.write(run.column + shift, run.row, &run.text);
        }
    }

    fn write(&mut self, column: f64, row: f64, text: &[char]) {
        if row < 0.0 || row >= self.rows as f64 {
            return;
        }
        for (k, &letter) in text.iter().enumerate() {
            let at = column + k as f64;
            if at >= 0.0 && at < self.columns as f64 {
                self.letters[row as usize * self.columns + at as usize] = Some(letter);
            }
        }
    }

    /// The grid as lines of text without trailing spaces, from the first
    /// line that holds anything to the last.
    fn text(&self) -> String {
        let filled = |row: usize| {
            let cells = row * self.columns..(row + 1) * self.columns;
            cells.clone().any(|cell| self.dots[cell] != 0)
                || self.letters[cells].iter().any(Option::is_some)
        };
        let first = (0..self.rows).find(|&row| filled(row)).unwrap_or(0);
        let last = (0..self.rows).rfind(|&row| filled(row)).unwrap_or(0);
        let mut out = String::new();
        for row in first..=last {
            let line: String = (0..self.columns)
                .map(|column| {
                    let cell = row * self.columns + column;
                    match (self.letters[cell], self.dots[cell]) {
                        (Some(letter), _) => letter,
                        (None, 0) => ' ',
                        (None, bits) => char::from_u32(0x2800 + bits).unwrap_or(' '),
                    }
                })
                .collect();
            out.push_str(line.trim_end());
            out.push('\n');
        }
        out
    }
}

/// The characters of a run, with subscripts and superscripts in their Unicode forms.
fn letters(run: &Text) -> Vec<char> {
    run.pieces
        .iter()
        .flat_map(|(text, level)| text.chars().map(move |c| raised(c, *level)))
        .collect()
}

fn raised(character: char, level: Level) -> char {
    const SUBSCRIPTS: &str = "₀₁₂₃₄₅₆₇₈₉";
    const SUPERSCRIPTS: &str = "⁰¹²³⁴⁵⁶⁷⁸⁹";
    let digit = |set: &str| {
        character
            .to_digit(10)
            .and_then(|d| set.chars().nth(d as usize))
    };
    match (level, character) {
        (Level::Base, _) => character,
        (Level::Superscript, '+') => '⁺',
        (Level::Superscript, '\u{2212}' | '-') => '⁻',
        (Level::Subscript, _) => digit(SUBSCRIPTS).unwrap_or(character),
        (Level::Superscript, _) => digit(SUPERSCRIPTS).unwrap_or(character),
    }
}
