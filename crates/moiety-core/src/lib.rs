//! moiety-core reads condensed organic formulas such as `(CH3)2CHCOOH` into
//! molecular graphs.
//!
//! It never guesses. Every reading that survives valence resolution is returned,
//! ranked by a declared policy, and an ambiguous input has a different result
//! shape from a resolved one. The rules are specified in `docs/grammar.md`.
//!
//! ```
//! use moiety_core::{Dictionary, Interpretation, Options, interpret};
//!
//! let result = interpret("(CH3)2CHCOOH", &Dictionary::default(), &Options::default());
//! let Interpretation::Unique(reading) = result else { panic!("expected one reading") };
//! assert_eq!(reading.smiles, "CC(C)C(=O)O");
//! ```

pub mod aromatic;
mod assembly;
mod classify;
mod diagnostic;
pub mod dictionary;
mod element;
mod interpret;
pub mod lexer;
mod molecule;
pub mod normalize;
mod parser;
mod reading;
pub mod smiles;
mod span;
mod suggest;
mod template;

pub use diagnostic::{Code, Diagnostic, Label, Severity};
pub use dictionary::{Dictionary, DictionaryError, Entry};
pub use element::{Category, Element};
pub use interpret::{
    ChargeConvention, GRAMMAR_VERSION, MAX_ATOMS, MAX_INPUT_LEN, Options, RANKING_POLICY_VERSION,
    interpret,
};
pub use molecule::{Atom, AtomKind, Bond, BondOrder, DoubleBondStereo, Molecule, Tetrahedral};
pub use reading::{
    Ambiguity, AmbiguityClass, Direction, Interpretation, Provenance, Reading, Underdetermination,
};
pub use span::Span;
pub use template::{MAX_POSITIONS, RingTemplate};
