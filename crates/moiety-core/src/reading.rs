use std::fmt;

use crate::diagnostic::Diagnostic;
use crate::molecule::Molecule;
use crate::normalize::{self, Fold};
use crate::span::Span;

/// The result of reading a formula. Ambiguous input has a different shape from
/// resolved input, so a caller cannot mistake one for the other.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Interpretation {
    /// Exactly one molecule survives.
    Unique(Reading),
    /// Two or more molecules survive, ranked by the ranking policy.
    Ambiguous {
        /// Every surviving reading, rank 1 first.
        readings: Vec<Reading>,
        /// The classes of every ambiguity separating the readings.
        classes: Vec<AmbiguityClass>,
    },
    /// The formula is well formed but does not encode a structure.
    Underdetermined {
        /// Span of the underdetermined part, in the original input.
        span: Span,
        /// Why no structure can be derived.
        reason: Underdetermination,
    },
    /// The formula has no reading.
    Invalid {
        /// Errors, most relevant first.
        diagnostics: Vec<Diagnostic>,
    },
}

impl Interpretation {
    /// Every reading, rank 1 first. Empty for underdetermined and invalid input.
    #[must_use]
    pub fn readings(&self) -> &[Reading] {
        match self {
            Interpretation::Unique(reading) => std::slice::from_ref(reading),
            Interpretation::Ambiguous { readings, .. } => readings,
            _ => &[],
        }
    }

    /// The rank 1 reading, with its ambiguities still attached.
    #[must_use]
    pub fn top(&self) -> Option<&Reading> {
        self.readings().first()
    }

    /// The kind of result: `unique`, `ambiguous`, `underdetermined` or `invalid`.
    #[must_use]
    pub fn kind(&self) -> &'static str {
        match self {
            Interpretation::Unique(_) => "unique",
            Interpretation::Ambiguous { .. } => "ambiguous",
            Interpretation::Underdetermined { .. } => "underdetermined",
            Interpretation::Invalid { .. } => "invalid",
        }
    }
}

/// One way to read a formula.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Reading {
    /// The molecular graph.
    pub molecule: Molecule,
    /// Canonical SMILES of the molecule.
    pub smiles: String,
    /// Position in the ranking, starting at 1.
    pub rank: u32,
    /// Ranking score under the ranking policy. Higher ranks first.
    pub score: i32,
    /// The direction the formula was written in.
    pub direction: Direction,
    /// The ambiguities this reading resolved, with the choice it made at each.
    pub ambiguities: Vec<Ambiguity>,
    /// Where every atom and bond came from in the input.
    pub provenance: Provenance,
    /// Warnings about parts of the input that were skipped or limited.
    pub warnings: Vec<Diagnostic>,
}

/// The direction a reading's bonds were written in. See the grammar, section 7.4.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Direction {
    /// Every unit bonds to a host on its left.
    LeftToRight,
    /// Every unit bonds to an anchor on its right.
    RightToLeft,
    /// Both.
    Mixed,
}

impl Direction {
    /// Name as used in output: `left-to-right`, `right-to-left` or `mixed`.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Direction::LeftToRight => "left-to-right",
            Direction::RightToLeft => "right-to-left",
            Direction::Mixed => "mixed",
        }
    }
}

/// The classes of the ambiguity taxonomy.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum AmbiguityClass {
    /// `C` + `O` against `Co`.
    SymbolSegmentation,
    /// An abbreviation against an element with the same spelling, such as tosyl against tennessine.
    AbbreviationShadowing,
    /// Spellings of the same molecule. Never separates readings.
    OrthographicVariant,
    /// Where a written charge sits, or which convention represents a hypervalent atom.
    ChargeRepresentation,
    /// The direction a group was written in.
    Direction,
    /// A unit read as a pendant in one reading and as a chain link in another.
    BranchPointFork,
    /// A trailing `-` read as a charge or as an open valence.
    ChargeVsAttachment,
    /// Connectivity the formula does not encode.
    TopologicalUnderdetermination,
    /// Notation outside the grammar, skipped with a warning.
    OutOfScope,
    /// Components parsed independently.
    MultiComponent,
}

impl AmbiguityClass {
    /// Kebab-case name, as used in the conformance suite.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            AmbiguityClass::SymbolSegmentation => "symbol-segmentation",
            AmbiguityClass::AbbreviationShadowing => "abbreviation-shadowing",
            AmbiguityClass::OrthographicVariant => "orthographic-variant",
            AmbiguityClass::ChargeRepresentation => "charge-representation",
            AmbiguityClass::Direction => "direction",
            AmbiguityClass::BranchPointFork => "branch-point-fork",
            AmbiguityClass::ChargeVsAttachment => "charge-vs-attachment",
            AmbiguityClass::TopologicalUnderdetermination => "topological-underdetermination",
            AmbiguityClass::OutOfScope => "out-of-scope",
            AmbiguityClass::MultiComponent => "multi-component",
        }
    }

    /// The class with kebab-case name `name`.
    #[must_use]
    pub fn parse(name: &str) -> Option<AmbiguityClass> {
        [
            AmbiguityClass::SymbolSegmentation,
            AmbiguityClass::AbbreviationShadowing,
            AmbiguityClass::OrthographicVariant,
            AmbiguityClass::ChargeRepresentation,
            AmbiguityClass::Direction,
            AmbiguityClass::BranchPointFork,
            AmbiguityClass::ChargeVsAttachment,
            AmbiguityClass::TopologicalUnderdetermination,
            AmbiguityClass::OutOfScope,
            AmbiguityClass::MultiComponent,
        ]
        .into_iter()
        .find(|class| class.as_str() == name)
    }
}

impl fmt::Display for AmbiguityClass {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// One ambiguity a reading resolved, and how.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Ambiguity {
    /// Class from the taxonomy.
    pub class: AmbiguityClass,
    /// Span in the original input.
    pub span: Span,
    /// The choice this reading made, such as `C O` or `Co`.
    pub choice: String,
}

/// Why a well-formed formula yields no structure.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Underdetermination {
    /// A molecular formula such as `C6H6`, which encodes no connectivity.
    MolecularFormula,
    /// A ring the dictionary cannot expand.
    UnknownRing,
    /// A token that matches no element and no dictionary entry.
    UnknownAbbreviation,
}

impl Underdetermination {
    /// Kebab-case name, as used in the conformance suite.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Underdetermination::MolecularFormula => "molecular-formula",
            Underdetermination::UnknownRing => "unknown-ring",
            Underdetermination::UnknownAbbreviation => "unknown-abbreviation",
        }
    }
}

impl fmt::Display for Underdetermination {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Byte spans from every atom and bond of a reading back to the input.
///
/// Spans are stored in normalized coordinates. The fold table maps them back
/// to the original input, so `CH₃` still underlines correctly.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Provenance {
    atoms: Vec<Span>,
    bonds: Vec<Span>,
    folds: Vec<Fold>,
}

impl Provenance {
    /// Builds provenance from per-atom and per-bond normalized spans and the fold table.
    #[must_use]
    pub fn new(atoms: Vec<Span>, bonds: Vec<Span>, folds: Vec<Fold>) -> Self {
        Self {
            atoms,
            bonds,
            folds,
        }
    }

    /// Normalized span of every atom, by atom index.
    #[must_use]
    pub fn atoms(&self) -> &[Span] {
        &self.atoms
    }

    /// Normalized span of every bond, by bond index.
    #[must_use]
    pub fn bonds(&self) -> &[Span] {
        &self.bonds
    }

    /// The normalization folds, mapping normalized spans to original spans.
    #[must_use]
    pub fn folds(&self) -> &[Fold] {
        &self.folds
    }

    /// Maps a normalized span to the original input.
    #[must_use]
    pub fn to_original(&self, span: Span) -> Span {
        normalize::to_original(&self.folds, span)
    }

    /// Original-input span of atom `index`.
    #[must_use]
    pub fn atom_source(&self, index: usize) -> Option<Span> {
        self.atoms.get(index).map(|&span| self.to_original(span))
    }

    /// Original-input span of bond `index`.
    #[must_use]
    pub fn bond_source(&self, index: usize) -> Option<Span> {
        self.bonds.get(index).map(|&span| self.to_original(span))
    }

    /// Atoms whose original span overlaps `span` of the original input.
    #[must_use]
    pub fn atoms_in(&self, span: Span) -> Vec<usize> {
        (0..self.atoms.len())
            .filter(|&i| self.atom_source(i).is_some_and(|s| s.overlaps(span)))
            .collect()
    }
}
