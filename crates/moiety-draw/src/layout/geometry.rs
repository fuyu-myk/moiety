//! Points, vectors and bounding boxes in the plane.

use std::f64::consts::TAU;
use std::ops::{Add, Mul, Neg, Sub};

/// Lengths below this are treated as zero.
pub(crate) const EPSILON: f64 = 1e-9;

/// How much narrower than the widest gap, in radians, an open gap may be and
/// still take a substituent: 15 degrees.
const NEAR_TIE: f64 = 0.26;

/// A point, or a vector, in the plane.
///
/// In a [`Layout`](super::Layout) one unit is one bond length. `x` grows to the
/// right and `y` grows downward, as in SVG, so a positive angle turns clockwise
/// on screen.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Point {
    /// Horizontal coordinate.
    pub x: f64,
    /// Vertical coordinate, growing downward.
    pub y: f64,
}

impl Point {
    /// The point `(x, y)`.
    #[must_use]
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// The unit vector at `angle` radians from the positive `x` axis.
    #[must_use]
    pub fn from_angle(angle: f64) -> Self {
        let (sin, cos) = angle.sin_cos();
        Self::new(cos, sin)
    }

    /// Length of the vector.
    #[must_use]
    pub fn length(self) -> f64 {
        self.x.hypot(self.y)
    }

    /// Distance to `other`.
    #[must_use]
    pub fn distance(self, other: Point) -> f64 {
        (self - other).length()
    }

    /// Angle of the vector from the positive `x` axis, in radians.
    #[must_use]
    pub fn angle(self) -> f64 {
        self.y.atan2(self.x)
    }

    /// The vector rotated by `angle` radians.
    #[must_use]
    pub fn rotated(self, angle: f64) -> Self {
        let (sin, cos) = angle.sin_cos();
        Self::new(self.x * cos - self.y * sin, self.x * sin + self.y * cos)
    }

    /// The vector scaled to length 1, or `None` when it is too short to have a direction.
    #[must_use]
    pub fn normalized(self) -> Option<Self> {
        let length = self.length();
        (length > EPSILON).then(|| self * length.recip())
    }

    /// The vector turned a quarter turn, so that `v.cross(v.perpendicular())` is positive.
    #[must_use]
    pub fn perpendicular(self) -> Self {
        Self::new(-self.y, self.x)
    }

    /// Dot product.
    #[must_use]
    pub fn dot(self, other: Point) -> f64 {
        self.x * other.x + self.y * other.y
    }

    /// The `z` component of the cross product: positive when `other` lies on the
    /// side that [`perpendicular`](Self::perpendicular) points to.
    #[must_use]
    pub fn cross(self, other: Point) -> f64 {
        self.x * other.y - self.y * other.x
    }

    /// Mirror image across the line through `origin` along the unit vector `direction`.
    #[must_use]
    pub fn reflected(self, origin: Point, direction: Point) -> Self {
        let offset = self - origin;
        origin + direction * (2.0 * offset.dot(direction)) - offset
    }
}

impl Add for Point {
    type Output = Point;

    fn add(self, other: Point) -> Point {
        Point::new(self.x + other.x, self.y + other.y)
    }
}

impl Sub for Point {
    type Output = Point;

    fn sub(self, other: Point) -> Point {
        Point::new(self.x - other.x, self.y - other.y)
    }
}

impl Mul<f64> for Point {
    type Output = Point;

    fn mul(self, factor: f64) -> Point {
        Point::new(self.x * factor, self.y * factor)
    }
}

impl Neg for Point {
    type Output = Point;

    fn neg(self) -> Point {
        Point::new(-self.x, -self.y)
    }
}

/// Mean of the points, or `None` for no points.
pub(crate) fn centroid(points: impl IntoIterator<Item = Point>) -> Option<Point> {
    let (sum, count) = points
        .into_iter()
        .fold((Point::default(), 0.0_f64), |(sum, count), point| {
            (sum + point, count + 1.0)
        });
    (count > 0.0).then(|| sum * count.recip())
}

/// The widest angular gap between the given directions, as `(start angle, width)`.
///
/// One direction leaves a full turn starting at that direction; none leaves a
/// full turn starting at angle zero.
pub(crate) fn largest_gap(directions: impl IntoIterator<Item = f64>) -> (f64, f64) {
    let directions: Vec<(usize, f64)> = directions.into_iter().enumerate().collect();
    open_gap(&directions, |_, _| false)
}

/// The widest angular gap between the directions to `neighbours`, given as
/// `(neighbour, angle)`, as `(start angle, width)`. A gap that a ring fills,
/// between two neighbours for which `filled` holds, gives way to an open gap
/// less than [`NEAR_TIE`] narrower, so that substituents point out of rings
/// whose angles are nearly even.
pub(crate) fn open_gap(
    neighbours: &[(usize, f64)],
    filled: impl Fn(usize, usize) -> bool,
) -> (f64, f64) {
    let mut sorted: Vec<(usize, f64)> = neighbours
        .iter()
        .map(|&(n, angle)| (n, angle.rem_euclid(TAU)))
        .collect();
    sorted.sort_by(|x, y| x.1.total_cmp(&y.1));
    let count = sorted.len();
    let gaps: Vec<(f64, f64, bool)> = (0..count)
        .map(|i| {
            let ((a, from), (b, to)) = (sorted[i], sorted[(i + 1) % count]);
            let width = if i + 1 < count {
                to - from
            } else {
                to + TAU - from
            };
            (from, width, count > 1 && filled(a, b))
        })
        .collect();
    let widest = |open_only: bool| {
        gaps.iter().filter(|gap| !(open_only && gap.2)).fold(
            None,
            |best: Option<(f64, f64)>, &(start, width, _)| match best {
                Some(best) if width <= best.1 + EPSILON => Some(best),
                _ => Some((start, width)),
            },
        )
    };
    let Some(any) = widest(false) else {
        return (0.0, TAU);
    };
    match widest(true) {
        Some(open) if open.1 > any.1 - NEAR_TIE => open,
        _ => any,
    }
}

/// Circumradius of a regular polygon with `sides` sides of length 1.
pub(crate) fn circumradius(sides: usize) -> f64 {
    0.5 / (std::f64::consts::PI / sides as f64).sin()
}

/// An axis-aligned box that grows to cover what it is shown.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct Bounds {
    extent: Option<(Point, Point)>,
}

impl Bounds {
    /// Grows the box to cover `point`.
    pub(crate) fn include(&mut self, point: Point) {
        let (min, max) = self.extent.get_or_insert((point, point));
        *min = Point::new(min.x.min(point.x), min.y.min(point.y));
        *max = Point::new(max.x.max(point.x), max.y.max(point.y));
    }

    /// Grows the box to cover the box centred at `centre` with half extents `half_width` and `half_height`.
    pub(crate) fn include_box(&mut self, centre: Point, half_width: f64, half_height: f64) {
        self.include(centre - Point::new(half_width, half_height));
        self.include(centre + Point::new(half_width, half_height));
    }

    /// The lowest and highest corners, or `None` when nothing was included.
    pub(crate) fn extent(&self) -> Option<(Point, Point)> {
        self.extent
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_PI_2, PI};

    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn rotates_and_reflects() {
        let v = Point::new(1.0, 0.0).rotated(FRAC_PI_2);
        assert!(close(v.x, 0.0) && close(v.y, 1.0));
        let r = Point::new(1.0, 1.0).reflected(Point::default(), Point::new(1.0, 0.0));
        assert!(close(r.x, 1.0) && close(r.y, -1.0));
        assert!(Point::new(1.0, 0.0).cross(Point::new(1.0, 0.0).perpendicular()) > 0.0);
        assert_eq!(Point::default().normalized(), None);
    }

    #[test]
    fn finds_the_widest_gap() {
        let (start, width) = largest_gap([0.0, FRAC_PI_2]);
        assert!(close(start, FRAC_PI_2) && close(width, 1.5 * PI));
        let (start, width) = largest_gap([PI]);
        assert!(close(start, PI) && close(width, TAU));
        assert_eq!(largest_gap([]), (0.0, TAU));
    }

    #[test]
    fn prefers_a_gap_no_ring_fills_when_it_is_nearly_the_widest() {
        let neighbours = [(1, 0.0), (2, 2.15), (3, 4.2)];
        let filled = |a: usize, b: usize| a + b != 4;
        let (start, width) = open_gap(&neighbours, filled);
        assert!(close(start, 4.2) && close(width, TAU - 4.2));
        let (start, _) = open_gap(&neighbours, |_, _| false);
        assert!(close(start, 0.0));
        let lopsided = [(1, 0.0), (2, 2.6), (3, 4.2)];
        let (start, _) = open_gap(&lopsided, filled);
        assert!(close(start, 0.0), "a much wider filled gap still wins");
        assert_eq!(open_gap(&[(7, PI)], filled), (PI, TAU));
    }

    #[test]
    fn measures_polygons_and_boxes() {
        assert!(close(circumradius(6), 1.0));
        let mut bounds = Bounds::default();
        assert_eq!(bounds.extent(), None);
        bounds.include_box(Point::new(1.0, 1.0), 0.5, 0.25);
        bounds.include(Point::new(-1.0, 0.0));
        assert_eq!(
            bounds.extent(),
            Some((Point::new(-1.0, 0.0), Point::new(1.5, 1.25)))
        );
        assert_eq!(
            centroid([Point::new(0.0, 0.0), Point::new(2.0, 2.0)]),
            Some(Point::new(1.0, 1.0))
        );
    }
}
