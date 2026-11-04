//! Exact-area coverage of filled outlines on a pixel grid.
//!
//! Every edge adds, to each pixel it crosses, the signed share of that pixel
//! that lies to its right; a running sum along each row then gives the
//! coverage. Outlines may hold holes, as letters do, under the nonzero rule.

use crate::layout::Point;

/// Edges of closed outlines accumulated over a grid of pixels.
pub(super) struct Coverage {
    width: usize,
    height: usize,
    cells: Vec<f32>,
}

impl Coverage {
    /// An empty grid `width` by `height` pixels.
    pub(super) fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            cells: vec![0.0; width * height + 2],
        }
    }

    /// Adds every edge of the closed polygon through `points`.
    pub(super) fn polygon(&mut self, points: &[Point]) {
        for (i, &a) in points.iter().enumerate() {
            self.edge(a, points[(i + 1) % points.len()]);
        }
    }

    /// Adds the edge from `a` to `b`, in pixels. The part above or below the
    /// grid is left out, and a part beyond its left or right side runs along
    /// that side instead, which leaves the coverage of the grid unchanged.
    pub(super) fn edge(&mut self, a: Point, b: Point) {
        let finite = |p: Point| p.x.is_finite() && p.y.is_finite();
        if !finite(a) || !finite(b) {
            return;
        }
        let height = self.height as f64;
        let at = |t: f64| a + (b - a) * t;
        let crossing = |from: f64, to: f64, value: f64| {
            ((from - value) * (to - value) < 0.0).then(|| (value - from) / (to - from))
        };
        let (mut low, mut high) = (0.0_f64, 1.0_f64);
        for bound in [0.0, height] {
            if let Some(t) = crossing(a.y, b.y, bound) {
                let entering = (b.y - a.y) * (if bound == 0.0 { 1.0 } else { -1.0 }) > 0.0;
                if entering {
                    low = low.max(t)
                } else {
                    high = high.min(t)
                }
            }
        }
        let inside = |y: f64| (0.0..=height).contains(&y);
        if low >= high || !inside(at(0.5 * (low + high)).y) {
            return;
        }
        let right = self.width as f64 - 1e-3;
        let mut cuts = [low, high, low, low];
        for (k, side) in [0.0, right].into_iter().enumerate() {
            if let Some(t) = crossing(a.x, b.x, side).filter(|t| (low..=high).contains(t)) {
                cuts[2 + k] = t;
            }
        }
        cuts.sort_by(f64::total_cmp);
        let clamp = |p: Point| Point::new(p.x.clamp(0.0, right), p.y.clamp(0.0, height));
        for pair in cuts.windows(2) {
            self.rows(clamp(at(pair[0])), clamp(at(pair[1])));
        }
    }

    /// Adds an edge that lies within the grid.
    fn rows(&mut self, a: Point, b: Point) {
        if (a.y - b.y).abs() < 1e-9 {
            return;
        }
        let (sign, top, bottom) = if a.y < b.y { (1.0, a, b) } else { (-1.0, b, a) };
        let slope = (bottom.x - top.x) / (bottom.y - top.y);
        let mut x = top.x;
        for row in (top.y.floor() as usize)..(bottom.y.ceil() as usize).min(self.height) {
            let dy = (row as f64 + 1.0).min(bottom.y) - (row as f64).max(top.y);
            let next = x + slope * dy;
            self.span(row, x, next, sign * dy);
            x = next;
        }
    }

    /// Adds, along `row`, an edge that runs from `x0` to `x1` while falling `dy`.
    fn span(&mut self, row: usize, x0: f64, x1: f64, dy: f64) {
        let start = row * self.width;
        let (left, right) = if x0 < x1 { (x0, x1) } else { (x1, x0) };
        let first = left.floor();
        let last = right.ceil();
        let column = first as usize;
        if last <= first + 1.0 {
            let middle = 0.5 * (x0 + x1) - first;
            self.cells[start + column] += (dy * (1.0 - middle)) as f32;
            self.cells[start + column + 1] += (dy * middle) as f32;
            return;
        }
        let step = (right - left).recip();
        let head = left - first;
        let tail = right - last + 1.0;
        let a0 = 0.5 * step * (1.0 - head) * (1.0 - head);
        let am = 0.5 * step * tail * tail;
        let end = last as usize;
        self.cells[start + column] += (dy * a0) as f32;
        if end == column + 2 {
            self.cells[start + column + 1] += (dy * (1.0 - a0 - am)) as f32;
        } else {
            let a1 = step * (1.5 - head);
            self.cells[start + column + 1] += (dy * (a1 - a0)) as f32;
            for cell in column + 2..end - 1 {
                self.cells[start + cell] += (dy * step) as f32;
            }
            let a2 = a1 + (end - column - 3) as f64 * step;
            self.cells[start + end - 1] += (dy * (1.0 - a2 - am)) as f32;
        }
        self.cells[start + end] += (dy * am) as f32;
    }

    /// The coverage of each pixel, from 0 to 1, row by row.
    pub(super) fn finish(mut self) -> Vec<f32> {
        let mut sum = 0.0;
        self.cells.truncate(self.width * self.height);
        for cell in &mut self.cells {
            sum += *cell;
            *cell = sum.abs().min(1.0);
        }
        self.cells
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn area(width: usize, height: usize, polygons: &[&[(f64, f64)]]) -> (f32, Vec<f32>) {
        let mut coverage = Coverage::new(width, height);
        for polygon in polygons {
            let points: Vec<Point> = polygon.iter().map(|&(x, y)| Point::new(x, y)).collect();
            coverage.polygon(&points);
        }
        let cells = coverage.finish();
        (cells.iter().sum(), cells)
    }

    fn close(a: f32, b: f32) -> bool {
        (a - b).abs() < 1e-3
    }

    #[test]
    fn covers_a_rectangle_exactly() {
        let (total, cells) = area(10, 6, &[&[(2.5, 1.0), (7.5, 1.0), (7.5, 4.0), (2.5, 4.0)]]);
        assert!(close(total, 15.0), "{total}");
        assert!(close(cells[10 + 2], 0.5) && close(cells[10 + 3], 1.0));
        assert!(close(cells[10 + 7], 0.5) && close(cells[10 + 8], 0.0));
        assert!(close(cells[0], 0.0) && close(cells[50 + 3], 0.0));
    }

    #[test]
    fn covers_slanted_edges_by_area() {
        let (total, _) = area(12, 12, &[&[(1.0, 1.0), (11.0, 1.0), (1.0, 9.0)]]);
        assert!(close(total, 40.0), "{total}");
        let (reversed, _) = area(12, 12, &[&[(1.0, 9.0), (11.0, 1.0), (1.0, 1.0)]]);
        assert!(close(reversed, 40.0), "{reversed}");
    }

    #[test]
    fn leaves_holes_open() {
        let outer: &[(f64, f64)] = &[(1.0, 1.0), (9.0, 1.0), (9.0, 9.0), (1.0, 9.0)];
        let inner: &[(f64, f64)] = &[(3.0, 3.0), (3.0, 7.0), (7.0, 7.0), (7.0, 3.0)];
        let (total, cells) = area(10, 10, &[outer, inner]);
        assert!(close(total, 48.0), "{total}");
        assert!(close(cells[50 + 5], 0.0) && close(cells[50 + 2], 1.0));
    }

    #[test]
    fn clips_edges_that_cross_the_border() {
        let (left, _) = area(10, 10, &[&[(-4.0, 0.0), (4.0, 8.0), (-4.0, 8.0)]]);
        assert!(close(left, 8.0), "{left}");
        let (top, _) = area(10, 10, &[&[(0.0, -4.0), (8.0, 4.0), (0.0, 4.0)]]);
        assert!(close(top, 24.0), "{top}");
        let (lost, _) = area(4, 4, &[&[(1.0, 1.0), (f64::NAN, 2.0), (2.0, 3.0)]]);
        assert!(lost.is_finite());
    }

    #[test]
    fn clamps_outlines_to_the_grid() {
        let (total, _) = area(
            4,
            4,
            &[&[(-3.0, -3.0), (9.0, -3.0), (9.0, 9.0), (-3.0, 9.0)]],
        );
        assert!(total <= 16.0 + 1e-3 && total > 15.0, "{total}");
    }
}
