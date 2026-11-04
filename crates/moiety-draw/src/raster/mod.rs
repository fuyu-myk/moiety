//! Pixel drawings of skeletal structures.
//!
//! The drawing is the one [`render_svg`](crate::render_svg) writes, painted
//! with antialiasing into RGBA pixels at the same size in pixels: lines with
//! round ends, filled wedges, and labels set in Inter, whose glyph outlines
//! ship with this crate under the SIL Open Font License 1.1.
//!
//! ```
//! let molecule = moiety::smiles::parse("c1ccccc1O").expect("phenol");
//! let image = moiety_draw::render_image(&molecule, &moiety_draw::RasterOptions::default());
//! assert_eq!(image.pixels().len(), (image.width() * image.height() * 4) as usize);
//! assert!(image.to_png().starts_with(b"\x89PNG"));
//! ```

mod fill;
mod glyphs;
mod png;
mod text;

use moiety_core::Molecule;

use crate::layout::Point;
use crate::svg::shape::Shape;
use crate::svg::{SvgOptions, picture};
use fill::Coverage;
use text::{Step, curve, outlines};

/// Most pixels a drawing may hold; larger drawings are scaled down to fit.
const MAX_PIXELS: f64 = (1 << 24) as f64;

/// Sizes and colours of a pixel drawing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RasterOptions {
    /// Length of one bond, in pixels.
    pub bond_length: f64,
    /// Font size of atom labels, in pixels.
    pub font_size: f64,
    /// Width of bond lines, in pixels.
    pub line_width: f64,
    /// Empty margin around the drawing, in pixels.
    pub padding: f64,
    /// Colour of bonds and labels, as red, green, blue and alpha.
    pub ink: [u8; 4],
    /// Colour behind the drawing, as red, green, blue and alpha.
    pub background: [u8; 4],
}

impl Default for RasterOptions {
    /// The sizes of [`SvgOptions::default`], in black on a transparent background.
    fn default() -> Self {
        let svg = SvgOptions::default();
        Self {
            bond_length: svg.bond_length,
            font_size: svg.font_size,
            line_width: svg.line_width,
            padding: svg.padding,
            ink: [0, 0, 0, 255],
            background: [0, 0, 0, 0],
        }
    }
}

/// An RGBA image, row by row from the top, without premultiplied alpha.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Image {
    width: u32,
    height: u32,
    pixels: Vec<u8>,
}

impl Image {
    /// Width in pixels.
    #[must_use]
    pub fn width(&self) -> u32 {
        self.width
    }

    /// Height in pixels.
    #[must_use]
    pub fn height(&self) -> u32 {
        self.height
    }

    /// Red, green, blue and alpha of every pixel, row by row from the top.
    #[must_use]
    pub fn pixels(&self) -> &[u8] {
        &self.pixels
    }

    /// The image as a PNG file.
    #[must_use]
    pub fn to_png(&self) -> Vec<u8> {
        png::encode(self.width, self.height, &self.pixels)
    }
}

/// Draws `molecule` as pixels.
///
/// Sizes in `options` that are not positive and finite fall back to their
/// defaults. A drawing that would hold more than 2²⁴ pixels is scaled down to
/// that many.
#[must_use]
pub fn render_image(molecule: &Molecule, options: &RasterOptions) -> Image {
    Mask::new(molecule, options, 1.0, f64::INFINITY).to_image(options.ink, options.background)
}

/// Draws `molecule` as a PNG file, as [`render_image`] paints it.
#[must_use]
pub fn render_png(molecule: &Molecule, options: &RasterOptions) -> Vec<u8> {
    render_image(molecule, options).to_png()
}

/// `ink` laid over `background` with the share `alpha` of its own opacity.
pub(crate) fn blend(background: [u8; 4], ink: [u8; 4], alpha: f32) -> [u8; 4] {
    let top = alpha * f32::from(ink[3]) / 255.0;
    let bottom = f32::from(background[3]) / 255.0 * (1.0 - top);
    let total = top + bottom;
    if total <= 0.0 {
        return [0; 4];
    }
    let channel = |i: usize| {
        ((f32::from(ink[i]) * top + f32::from(background[i]) * bottom) / total).round() as u8
    };
    [
        channel(0),
        channel(1),
        channel(2),
        (total * 255.0).round() as u8,
    ]
}

/// How much of each pixel the ink of a drawing covers.
pub(crate) struct Mask {
    pub(crate) width: usize,
    pub(crate) height: usize,
    /// Coverage from 0 to 1, row by row from the top.
    pub(crate) alpha: Vec<f32>,
}

impl Mask {
    /// The ink of the drawing of `molecule`, enlarged `zoom` times, or less
    /// so that it is at most `max_width` pixels wide.
    pub(crate) fn new(
        molecule: &Molecule,
        options: &RasterOptions,
        zoom: f64,
        max_width: f64,
    ) -> Self {
        let svg = SvgOptions {
            bond_length: options.bond_length,
            font_size: options.font_size,
            line_width: options.line_width,
            padding: options.padding,
            title: None,
        };
        let picture = picture(molecule, &svg);
        let (min, max) = picture.bounds.extent().unwrap_or_default();
        let margin = Point::new(picture.padding, picture.padding);
        let origin = min - margin;
        let size = max - min + margin * 2.0;
        let fit = MAX_PIXELS.sqrt() / (size.x.max(1.0).sqrt() * size.y.max(1.0).sqrt());
        let scale = zoom.min(fit).min(max_width / size.x.max(1.0));
        if !(scale.is_finite() && scale > 0.0) {
            return Mask::empty();
        }
        let round =
            |f: fn(f64) -> f64| ((f(size.x * scale)).max(1.0), (f(size.y * scale)).max(1.0));
        let (mut width, mut height) = round(f64::ceil);
        if width * height > MAX_PIXELS {
            (width, height) = round(f64::floor);
        }
        let (width, height) = (width as usize, height as usize);
        let mut painter = Painter {
            mask: Mask {
                width,
                height,
                alpha: vec![0.0; width * height],
            },
            origin,
            scale,
            half_line: picture.line_width * scale / 2.0,
        };
        for shape in &picture.shapes {
            painter.shape(shape);
        }
        painter.mask
    }

    /// A single uncovered pixel, for a drawing too large to measure.
    fn empty() -> Self {
        Mask {
            width: 1,
            height: 1,
            alpha: vec![0.0],
        }
    }

    /// The mask as `ink` laid over `background`.
    pub(crate) fn to_image(&self, ink: [u8; 4], background: [u8; 4]) -> Image {
        Image {
            width: self.width as u32,
            height: self.height as u32,
            pixels: {
                let mut pixels = Vec::with_capacity(self.alpha.len() * 4);
                for &alpha in &self.alpha {
                    pixels.extend(blend(background, ink, alpha));
                }
                pixels
            },
        }
    }

    fn darken(&mut self, x: usize, y: usize, alpha: f32) {
        let cell = &mut self.alpha[y * self.width + x];
        *cell = cell.max(alpha);
    }
}

/// Paints shapes given in drawing coordinates into a mask.
struct Painter {
    mask: Mask,
    origin: Point,
    scale: f64,
    half_line: f64,
}

impl Painter {
    fn to_pixels(&self, point: Point) -> Point {
        (point - self.origin) * self.scale
    }

    fn shape(&mut self, shape: &Shape) {
        match shape {
            Shape::Line(start, end) => self.stroke(*start, *end),
            Shape::Hash(lines) => {
                for &(start, end) in lines {
                    self.stroke(start, end);
                }
            }
            Shape::Wedge(corners) => {
                self.fill(&[corners.iter().map(|&p| self.to_pixels(p)).collect()]);
                for (i, &corner) in corners.iter().enumerate() {
                    self.stroke(corner, corners[(i + 1) % 3]);
                }
            }
            Shape::Ring { centre, radius } => self.ring(*centre, *radius),
            Shape::Wave { start, curves } => {
                let mut from = *start;
                for &(control, to) in curves {
                    for point in curve(from, control, to) {
                        self.stroke(from, point);
                        from = point;
                    }
                }
            }
            Shape::Label { runs, .. } => {
                for run in runs {
                    for glyph in outlines(run) {
                        let pixels: Vec<Vec<Point>> = glyph
                            .iter()
                            .map(|contour| contour.iter().map(|&p| self.to_pixels(p)).collect())
                            .collect();
                        self.fill(&pixels);
                    }
                }
            }
        }
    }

    /// Paints the pixels of the box from `low` to `high` with the coverage
    /// `coverage` gives at each pixel centre.
    fn each_pixel(&mut self, low: Point, high: Point, coverage: impl Fn(Point) -> f64) {
        let (width, height) = (self.mask.width as f64, self.mask.height as f64);
        let x0 = low.x.floor().clamp(0.0, width) as usize;
        let x1 = high.x.ceil().clamp(0.0, width) as usize;
        let y0 = low.y.floor().clamp(0.0, height) as usize;
        let y1 = high.y.ceil().clamp(0.0, height) as usize;
        for y in y0..y1 {
            for x in x0..x1 {
                let alpha = coverage(Point::new(x as f64 + 0.5, y as f64 + 0.5));
                if alpha > 0.0 {
                    self.mask.darken(x, y, alpha.min(1.0) as f32);
                }
            }
        }
    }

    /// A line of the drawing's line width from `start` to `end`, with round ends.
    fn stroke(&mut self, start: Point, end: Point) {
        let (a, b) = (self.to_pixels(start), self.to_pixels(end));
        let reach = self.half_line + 1.0;
        let low = Point::new(a.x.min(b.x) - reach, a.y.min(b.y) - reach);
        let high = Point::new(a.x.max(b.x) + reach, a.y.max(b.y) + reach);
        let half = self.half_line;
        self.each_pixel(low, high, |p| half + 0.5 - segment_distance(p, a, b));
    }

    /// A circle of the drawing's line width.
    fn ring(&mut self, centre: Point, radius: f64) {
        let centre = self.to_pixels(centre);
        let radius = radius * self.scale;
        let reach = Point::new(radius + self.half_line + 1.0, radius + self.half_line + 1.0);
        let half = self.half_line;
        self.each_pixel(centre - reach, centre + reach, |p| {
            half + 0.5 - (p.distance(centre) - radius).abs()
        });
    }

    /// The area inside `contours`, given in pixels, by the nonzero rule.
    fn fill(&mut self, contours: &[Vec<Point>]) {
        let points = contours.iter().flatten();
        let (low, high) = points.fold(
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
        );
        if low.x > high.x {
            return;
        }
        let x0 = low.x.floor().max(0.0);
        let y0 = low.y.floor().max(0.0);
        let x1 = high.x.ceil().min(self.mask.width as f64);
        let y1 = high.y.ceil().min(self.mask.height as f64);
        if x1 <= x0 || y1 <= y0 {
            return;
        }
        let corner = Point::new(x0 - 1.0, y0);
        let (width, height) = ((x1 - x0) as usize + 2, (y1 - y0) as usize);
        let mut coverage = Coverage::new(width, height);
        for contour in contours {
            let local: Vec<Point> = contour.iter().map(|&p| p - corner).collect();
            coverage.polygon(&local);
        }
        for (i, alpha) in coverage.finish().into_iter().enumerate() {
            let (x, y) = (i % width, i / width);
            let column = (corner.x + x as f64) as isize;
            if alpha > 0.0 && column >= 0 && (column as usize) < self.mask.width {
                self.mask.darken(column as usize, y0 as usize + y, alpha);
            }
        }
    }
}

/// Distance from `p` to the segment from `a` to `b`.
fn segment_distance(p: Point, a: Point, b: Point) -> f64 {
    let along = b - a;
    let length = along.dot(along);
    let t = if length > 0.0 {
        ((p - a).dot(along) / length).clamp(0.0, 1.0)
    } else {
        0.0
    };
    p.distance(a + along * t)
}

#[cfg(test)]
mod tests;
