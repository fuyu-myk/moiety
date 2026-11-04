//! Reading and canonical writing of SMILES.
//!
//! The reader accepts the organic subset, bracket atoms, branches, ring closures,
//! aromatic lowercase atoms and `*` attachment points. Stereo marks and atom
//! classes are read and dropped, since moiety does not represent them. The writer
//! always produces canonical SMILES with aromatic atoms in lowercase.

mod read;
mod stereo;
mod write;

use std::fmt;

pub use read::parse;
pub use write::write;
pub(crate) use write::write_within;

/// A SMILES string that could not be read.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SmilesError {
    /// Byte offset of the problem.
    pub position: usize,
    /// What went wrong.
    pub message: String,
}

impl fmt::Display for SmilesError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} at byte {}", self.message, self.position)
    }
}

impl std::error::Error for SmilesError {}

/// Lowest SMILES default valence of an organic-subset element at least `used`, if it has one.
pub(crate) fn default_valence(symbol: &str, used: u8) -> Option<u8> {
    let valences: &[u8] = match symbol {
        "B" => &[3],
        "C" => &[4],
        "N" | "P" => &[3, 5],
        "O" => &[2],
        "S" => &[2, 4, 6],
        "F" | "Cl" | "Br" | "I" => &[1],
        _ => return None,
    };
    Some(
        valences
            .iter()
            .copied()
            .find(|&v| v >= used)
            .unwrap_or(used),
    )
}

/// The implicit hydrogens of an organic-subset atom whose bonds use `bonds` of its
/// valence: the lowest default valence that fits, less one when the atom is
/// aromatic, as OpenSMILES counts them. Aromatic `c` gets one, `n`, `o` and `s` none.
pub(crate) fn implicit_hydrogens(symbol: &str, bonds: u8, aromatic: bool) -> Option<u8> {
    let valence = default_valence(symbol, bonds)?;
    Some(
        valence
            .saturating_sub(bonds)
            .saturating_sub(u8::from(aromatic)),
    )
}

/// Whether `symbol` may be written in lowercase as an aromatic atom.
pub(crate) fn aromatic_symbol(symbol: &str) -> bool {
    matches!(symbol, "B" | "C" | "N" | "O" | "P" | "S" | "Se" | "As")
}
