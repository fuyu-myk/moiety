use std::fmt;

use crate::span::Span;

/// How serious a diagnostic is.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Severity {
    /// The formula has no reading.
    Error,
    /// The formula was read, but something in it was skipped or limited.
    Warning,
}

/// Stable diagnostic codes. See the grammar specification, section 8.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Code {
    /// CF0001: empty formula or empty component.
    Empty,
    /// CF0002: formula longer than the parser accepts.
    TooLong,
    /// CF0003: formula needs more search than the limit allows.
    TooComplex,
    /// CF0004: reading with more atoms than the limit allows.
    TooLarge,
    /// CF0101: character that no token can start with.
    UnrecognizedCharacter,
    /// CF0102: group opened and never closed.
    UnclosedGroup,
    /// CF0103: close bracket with no matching open.
    UnmatchedClose,
    /// CF0104: group opened with one bracket kind and closed with the other.
    MixedBrackets,
    /// CF0105: token where the grammar does not allow one.
    UnexpectedToken,
    /// CF0106: empty group.
    EmptyGroup,
    /// CF0107: count of zero or above the limit.
    CountOutOfRange,
    /// CF0201: a bond the two atoms cannot form.
    BondNotPossible,
    /// CF0202: valence left open at the end of a component.
    OpenValence,
    /// CF0203: leading groups with no anchor that has room for them.
    NoAnchor,
    /// CF0204: unit after a terminal group, with nothing to bond to.
    NoHost,
    /// CF0205: written charge that cannot be placed.
    ChargeNotPlaced,
    /// CF0206: repeated group with no free valence.
    GroupWithoutValence,
    /// CF0207: complete molecule used as a group.
    MoleculeAsGroup,
    /// CF0208: more pendants than the host atom has room for.
    HostFull,
    /// CF0301: stereo descriptor ignored.
    StereoIgnored,
    /// CF0302: isotope label read as hydrogen.
    IsotopeIgnored,
    /// CF0303: reading limit reached.
    ReadingLimit,
}

impl Code {
    /// The code as written in output, such as `CF0204`.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Code::Empty => "CF0001",
            Code::TooLong => "CF0002",
            Code::TooComplex => "CF0003",
            Code::TooLarge => "CF0004",
            Code::UnrecognizedCharacter => "CF0101",
            Code::UnclosedGroup => "CF0102",
            Code::UnmatchedClose => "CF0103",
            Code::MixedBrackets => "CF0104",
            Code::UnexpectedToken => "CF0105",
            Code::EmptyGroup => "CF0106",
            Code::CountOutOfRange => "CF0107",
            Code::BondNotPossible => "CF0201",
            Code::OpenValence => "CF0202",
            Code::NoAnchor => "CF0203",
            Code::NoHost => "CF0204",
            Code::ChargeNotPlaced => "CF0205",
            Code::GroupWithoutValence => "CF0206",
            Code::MoleculeAsGroup => "CF0207",
            Code::HostFull => "CF0208",
            Code::StereoIgnored => "CF0301",
            Code::IsotopeIgnored => "CF0302",
            Code::ReadingLimit => "CF0303",
        }
    }

    /// The code written as `CF` and four digits.
    #[must_use]
    pub fn parse(code: &str) -> Option<Code> {
        ALL_CODES.iter().copied().find(|c| c.as_str() == code)
    }

    /// Whether the code is an error or a warning.
    #[must_use]
    pub const fn severity(self) -> Severity {
        match self {
            Code::MixedBrackets
            | Code::StereoIgnored
            | Code::IsotopeIgnored
            | Code::ReadingLimit => Severity::Warning,
            _ => Severity::Error,
        }
    }
}

const ALL_CODES: [Code; 22] = [
    Code::Empty,
    Code::TooLong,
    Code::TooComplex,
    Code::TooLarge,
    Code::UnrecognizedCharacter,
    Code::UnclosedGroup,
    Code::UnmatchedClose,
    Code::MixedBrackets,
    Code::UnexpectedToken,
    Code::EmptyGroup,
    Code::CountOutOfRange,
    Code::BondNotPossible,
    Code::OpenValence,
    Code::NoAnchor,
    Code::NoHost,
    Code::ChargeNotPlaced,
    Code::GroupWithoutValence,
    Code::MoleculeAsGroup,
    Code::HostFull,
    Code::StereoIgnored,
    Code::IsotopeIgnored,
    Code::ReadingLimit,
];

impl fmt::Display for Code {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// A span with an explanation attached.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Label {
    /// Byte span in the original input.
    pub span: Span,
    /// What this span contributes to the problem.
    pub message: String,
}

impl Label {
    /// Creates a label.
    #[must_use]
    pub fn new(span: Span, message: impl Into<String>) -> Self {
        Self {
            span,
            message: message.into(),
        }
    }
}

/// A problem found in a formula, located by spans of the original input.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Diagnostic {
    /// Stable code.
    pub code: Code,
    /// One-line summary.
    pub message: String,
    /// The span the problem is about.
    pub primary: Label,
    /// Further spans that explain the conflict.
    pub secondary: Vec<Label>,
    /// Free-form notes, including any suggested correction.
    pub notes: Vec<String>,
}

impl Diagnostic {
    /// Creates a diagnostic with no secondary labels or notes.
    #[must_use]
    pub fn new(code: Code, message: impl Into<String>, primary: Label) -> Self {
        Self {
            code,
            message: message.into(),
            primary,
            secondary: Vec::new(),
            notes: Vec::new(),
        }
    }

    /// Adds a secondary label.
    #[must_use]
    pub fn with_label(mut self, label: Label) -> Self {
        self.secondary.push(label);
        self
    }

    /// Adds a note.
    #[must_use]
    pub fn with_note(mut self, note: impl Into<String>) -> Self {
        self.notes.push(note.into());
        self
    }

    /// Severity of the code.
    #[must_use]
    pub const fn severity(&self) -> Severity {
        self.code.severity()
    }

    /// Renders the diagnostic against the original input, in the style of rustc.
    ///
    /// ```text
    /// error[CF0204]: chain cannot continue past a terminal group
    ///   --> CH3CH3OH
    ///          ^^^ this group has no free valence left
    ///             ^^ but more units follow it
    /// ```
    #[must_use]
    pub fn render(&self, source: &str) -> String {
        const GUTTER: &str = "  --> ";
        let severity = match self.severity() {
            Severity::Error => "error",
            Severity::Warning => "warning",
        };
        let mut out = format!(
            "{severity}[{}]: {}\n{GUTTER}{source}\n",
            self.code, self.message
        );
        for label in std::iter::once(&self.primary).chain(&self.secondary) {
            let column = column(source, label.span.start);
            let width = column_width(source, label.span).max(1);
            out.push_str(&format!(
                "{}{} {}\n",
                " ".repeat(GUTTER.len() + column),
                "^".repeat(width),
                label.message
            ));
        }
        for note in &self.notes {
            out.push_str(&format!("  note: {note}\n"));
        }
        out
    }
}

fn column(source: &str, offset: usize) -> usize {
    source
        .get(..offset.min(source.len()))
        .map_or(offset, |prefix| prefix.chars().count())
}

fn column_width(source: &str, span: Span) -> usize {
    let start = span.start.min(source.len());
    let end = span.end.clamp(start, source.len());
    source
        .get(start..end)
        .map_or(span.len(), |text| text.chars().count())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn renders_labels_under_their_spans() {
        let diagnostic = Diagnostic::new(
            Code::NoHost,
            "chain cannot continue past a terminal group",
            Label::new(Span::new(3, 6), "this group has no free valence left"),
        )
        .with_label(Label::new(Span::new(6, 8), "but more units follow it"))
        .with_note("did you mean CH3CH2OH?");
        let expected = "error[CF0204]: chain cannot continue past a terminal group\n  --> CH3CH3OH\n         ^^^ this group has no free valence left\n            ^^ but more units follow it\n  note: did you mean CH3CH2OH?\n";
        assert_eq!(diagnostic.render("CH3CH3OH"), expected);
    }

    #[test]
    fn counts_columns_in_characters() {
        let diagnostic = Diagnostic::new(
            Code::UnrecognizedCharacter,
            "x",
            Label::new(Span::new(5, 6), "here"),
        );
        assert!(diagnostic.render("CH₃$").contains("\n         ^ here"));
    }

    #[test]
    fn codes_round_trip() {
        for code in ALL_CODES {
            assert_eq!(Code::parse(code.as_str()), Some(code));
        }
    }
}
