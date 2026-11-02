use std::fmt;

/// A half-open byte range `[start, end)` into a formula string.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Span {
    /// First byte of the range.
    pub start: usize,
    /// One past the last byte of the range.
    pub end: usize,
}

impl Span {
    /// Creates the span `[start, end)`.
    #[must_use]
    pub const fn new(start: usize, end: usize) -> Self {
        Self { start, end }
    }

    /// Number of bytes covered.
    #[must_use]
    pub const fn len(self) -> usize {
        self.end - self.start
    }

    /// Whether the span covers no bytes.
    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.start == self.end
    }

    /// The smallest span covering both `self` and `other`.
    #[must_use]
    pub fn cover(self, other: Span) -> Span {
        Span::new(self.start.min(other.start), self.end.max(other.end))
    }

    /// Whether `offset` lies inside the span.
    #[must_use]
    pub const fn contains(self, offset: usize) -> bool {
        self.start <= offset && offset < self.end
    }

    /// Whether the two spans share at least one byte.
    #[must_use]
    pub const fn overlaps(self, other: Span) -> bool {
        self.start < other.end && other.start < self.end
    }

    /// The same span moved right by `offset` bytes.
    #[must_use]
    pub const fn shift(self, offset: usize) -> Span {
        Span::new(self.start + offset, self.end + offset)
    }
}

impl fmt::Display for Span {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}..{}", self.start, self.end)
    }
}
