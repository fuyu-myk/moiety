//! moiety-core reads condensed organic formulas into molecular graphs and reports
//! every reading that survives valence resolution.

pub mod aromatic;
mod element;
mod molecule;
pub mod smiles;
mod span;

pub use element::{Category, Element};
pub use molecule::{Atom, AtomKind, Bond, BondOrder, Molecule};
pub use span::Span;
