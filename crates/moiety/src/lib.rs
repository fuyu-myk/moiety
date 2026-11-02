//! moiety reads condensed organic formulas such as `CH3CH2OH`, `(CH3)2CHCOOH` or
//! `EtO2CCH2CO2Et` into molecular graphs and canonical SMILES.
//!
//! It reports every reading that survives valence resolution and never returns a
//! single best guess for an ambiguous formula. This crate adds the bundled
//! abbreviation dictionary, an inverse writer, and the conformance runner to the
//! dependency-free [`moiety_core`].
//!
//! ```
//! use moiety::Interpretation;
//!
//! match moiety::interpret("COCl2") {
//!     Interpretation::Ambiguous { readings, classes } => {
//!         assert_eq!(readings[0].smiles, "O=C(Cl)Cl");
//!         assert_eq!(classes[0].as_str(), "symbol-segmentation");
//!     }
//!     other => panic!("expected two readings, got {}", other.kind()),
//! }
//! ```

pub mod conformance;
pub mod dictionary;
pub mod writer;

pub use dictionary::{LoadError, default_dictionary};
pub use moiety_core::*;

/// Reads a formula with the bundled dictionary and default options.
#[must_use]
pub fn interpret(input: &str) -> Interpretation {
    moiety_core::interpret(input, default_dictionary(), &Options::default())
}

/// Reads a formula with a chosen dictionary and options.
#[must_use]
pub fn interpret_with(input: &str, dictionary: &Dictionary, options: &Options) -> Interpretation {
    moiety_core::interpret(input, dictionary, options)
}
