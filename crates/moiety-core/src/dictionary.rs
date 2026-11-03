//! Abbreviation dictionaries: named fragments with attachment points and the spellings that name them.

use std::fmt;

use crate::molecule::{AtomKind, BondOrder, Molecule};
use crate::smiles::{self, SmilesError};
use crate::template::RingTemplate;

/// Where an abbreviation bonds to the rest of the formula.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Attachment {
    /// Index of the attachment atom in [`Entry::fragment`].
    pub atom: usize,
    /// Order of the bond the attachment makes.
    pub order: BondOrder,
}

/// One abbreviation: a fragment with zero, one or two attachment points and the spellings that name it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Entry {
    name: String,
    smiles: String,
    fragment: Molecule,
    attachments: Vec<Attachment>,
    spellings: Vec<String>,
    reversed: Vec<String>,
}

impl Entry {
    /// Defines an abbreviation from SMILES in which each `*` marks an attachment point.
    ///
    /// `spellings` are read with the attachment on the left, as in `COOH`.
    /// `reversed` are read with the attachment on the right, as in `HOOC`.
    ///
    /// # Errors
    ///
    /// Fails when the SMILES does not parse, an attachment point does not have exactly
    /// one neighbour, there are more than two attachment points, no spelling is given,
    /// or a spelling uses characters other than ASCII letters, digits and `-`, or
    /// starts with neither a letter nor a locant such as the `2-` of `2-C5H4N`.
    pub fn new(
        name: &str,
        smiles: &str,
        spellings: Vec<String>,
        reversed: Vec<String>,
    ) -> Result<Entry, DictionaryError> {
        let parsed = smiles::parse(smiles).map_err(|error| DictionaryError::Smiles {
            entry: name.to_string(),
            error,
        })?;
        if spellings.is_empty() && reversed.is_empty() {
            return Err(DictionaryError::NoSpellings {
                entry: name.to_string(),
            });
        }
        if let Some(bad) = spellings
            .iter()
            .chain(&reversed)
            .find(|s| !valid_spelling(s))
        {
            return Err(DictionaryError::BadSpelling {
                entry: name.to_string(),
                spelling: bad.clone(),
            });
        }
        let (fragment, attachments) =
            split_attachments(&parsed).ok_or_else(|| DictionaryError::BadAttachment {
                entry: name.to_string(),
            })?;
        Ok(Entry {
            name: name.to_string(),
            smiles: smiles.to_string(),
            fragment,
            attachments,
            spellings,
            reversed,
        })
    }

    /// Descriptive name, such as `carboxyl`.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The defining SMILES, with `*` attachment points.
    #[must_use]
    pub fn smiles(&self) -> &str {
        &self.smiles
    }

    /// The fragment without its attachment points.
    #[must_use]
    pub fn fragment(&self) -> &Molecule {
        &self.fragment
    }

    /// Attachment points, in the order their `*` appear in the SMILES.
    #[must_use]
    pub fn attachments(&self) -> &[Attachment] {
        &self.attachments
    }

    /// Spellings read with the attachment on the left.
    #[must_use]
    pub fn spellings(&self) -> &[String] {
        &self.spellings
    }

    /// Spellings read with the attachment on the right.
    #[must_use]
    pub fn reversed(&self) -> &[String] {
        &self.reversed
    }
}

fn valid_spelling(spelling: &str) -> bool {
    let digits = spelling.bytes().take_while(u8::is_ascii_digit).count();
    let word = match digits {
        0 => Some(spelling),
        1 | 2 if !spelling.starts_with('0') => spelling[digits..].strip_prefix('-'),
        _ => None,
    };
    word.is_some_and(|word| word.starts_with(|c: char| c.is_ascii_alphabetic()))
        && spelling
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || c == '-')
}

fn split_attachments(parsed: &Molecule) -> Option<(Molecule, Vec<Attachment>)> {
    let adjacency = parsed.adjacency();
    let mut attachments = Vec::new();
    let mut index = vec![None; parsed.atoms().len()];
    let mut fragment = Molecule::new();
    for (i, atom) in parsed.atoms().iter().enumerate() {
        if atom.kind != AtomKind::Attachment {
            index[i] = Some(fragment.add_atom(*atom));
        }
    }
    for (i, atom) in parsed.atoms().iter().enumerate() {
        if atom.kind != AtomKind::Attachment {
            continue;
        }
        let [(neighbour, bond)] = adjacency[i][..] else {
            return None;
        };
        attachments.push(Attachment {
            atom: index[neighbour]?,
            order: parsed.bonds()[bond].order,
        });
    }
    if attachments.len() > 2 || fragment.is_empty() {
        return None;
    }
    for bond in parsed.bonds() {
        if let (Some(a), Some(b)) = (index[bond.a], index[bond.b]) {
            fragment.add_bond(a, b, bond.order);
        }
    }
    Some((fragment, attachments))
}

/// A named, versioned set of abbreviations and ring templates.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Dictionary {
    name: String,
    version: String,
    entries: Vec<Entry>,
    rings: Vec<RingTemplate>,
}

/// A dictionary spelling found in the input.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SpellingMatch {
    /// Index of the entry in [`Dictionary::entries`].
    pub entry: usize,
    /// Length of the spelling in bytes.
    pub len: usize,
    /// Whether the spelling is one of the entry's reversed spellings.
    pub reversed: bool,
}

impl Dictionary {
    /// An empty dictionary.
    #[must_use]
    pub fn new(name: &str, version: &str) -> Self {
        Self {
            name: name.to_string(),
            version: version.to_string(),
            entries: Vec::new(),
            rings: Vec::new(),
        }
    }

    /// Dictionary name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Dictionary version.
    #[must_use]
    pub fn version(&self) -> &str {
        &self.version
    }

    /// All entries.
    #[must_use]
    pub fn entries(&self) -> &[Entry] {
        &self.entries
    }

    /// The entry named `name`.
    #[must_use]
    pub fn get(&self, name: &str) -> Option<&Entry> {
        self.entries.iter().find(|entry| entry.name == name)
    }

    /// The ring templates, besides the built-in benzene.
    #[must_use]
    pub fn rings(&self) -> &[RingTemplate] {
        &self.rings
    }

    /// The template numbered `id`: 0 is benzene, and `i + 1` is [`Dictionary::rings`]`[i]`.
    pub(crate) fn template(&self, id: usize) -> &RingTemplate {
        match id {
            0 => RingTemplate::benzene(),
            _ => &self.rings[id - 1],
        }
    }

    /// How many templates there are, benzene included.
    pub(crate) fn templates(&self) -> usize {
        self.rings.len() + 1
    }

    /// Adds a ring template, replacing any template with the same name, and returns
    /// the one it replaced.
    ///
    /// # Errors
    ///
    /// Fails when another template, or benzene, has the same formula.
    pub fn insert_ring(
        &mut self,
        ring: RingTemplate,
    ) -> Result<Option<RingTemplate>, DictionaryError> {
        let others = std::iter::once(RingTemplate::benzene()).chain(
            self.rings
                .iter()
                .filter(|other| other.name() != ring.name()),
        );
        for other in others {
            if other.formula() == ring.formula() {
                return Err(DictionaryError::DuplicateRing {
                    formula: ring.written(1),
                    first: other.name().to_string(),
                    second: ring.name().to_string(),
                });
            }
        }
        match self
            .rings
            .iter()
            .position(|other| other.name() == ring.name())
        {
            Some(index) => Ok(Some(std::mem::replace(&mut self.rings[index], ring))),
            None => {
                self.rings.push(ring);
                Ok(None)
            }
        }
    }

    /// Adds an entry, replacing any entry with the same name, and returns the one it replaced.
    ///
    /// # Errors
    ///
    /// Fails when one of the entry's spellings already names a different entry.
    pub fn insert(&mut self, entry: Entry) -> Result<Option<Entry>, DictionaryError> {
        for spelling in entry.spellings.iter().chain(&entry.reversed) {
            if let Some(owner) = self.entries.iter().find(|other| {
                other.name != entry.name
                    && other
                        .spellings
                        .iter()
                        .chain(&other.reversed)
                        .any(|s| s == spelling)
            }) {
                return Err(DictionaryError::DuplicateSpelling {
                    spelling: spelling.clone(),
                    first: owner.name.clone(),
                    second: entry.name.clone(),
                });
            }
        }
        match self
            .entries
            .iter()
            .position(|other| other.name == entry.name)
        {
            Some(index) => Ok(Some(std::mem::replace(&mut self.entries[index], entry))),
            None => {
                self.entries.push(entry);
                Ok(None)
            }
        }
    }

    /// Every spelling that starts at byte `pos` of `text`.
    #[must_use]
    pub fn matches_at(&self, text: &str, pos: usize) -> Vec<SpellingMatch> {
        let rest = &text[pos..];
        let mut found = Vec::new();
        for (index, entry) in self.entries.iter().enumerate() {
            let forward = entry.spellings.iter().map(|s| (s, false));
            let backward = entry.reversed.iter().map(|s| (s, true));
            for (spelling, reversed) in forward.chain(backward) {
                if rest.starts_with(spelling.as_str()) {
                    found.push(SpellingMatch {
                        entry: index,
                        len: spelling.len(),
                        reversed,
                    });
                }
            }
        }
        found
    }
}

/// A dictionary entry that could not be defined or added.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DictionaryError {
    /// The entry's SMILES does not parse.
    Smiles {
        /// Entry name.
        entry: String,
        /// The SMILES error.
        error: SmilesError,
    },
    /// The entry has no spellings.
    NoSpellings {
        /// Entry name.
        entry: String,
    },
    /// A spelling is empty or uses characters other than ASCII letters, digits and `-`.
    BadSpelling {
        /// Entry name.
        entry: String,
        /// The offending spelling.
        spelling: String,
    },
    /// An attachment point does not have exactly one neighbour, or there are more than two.
    BadAttachment {
        /// Entry name.
        entry: String,
    },
    /// Two entries claim the same spelling.
    DuplicateSpelling {
        /// The spelling.
        spelling: String,
        /// Entry that already owns it.
        first: String,
        /// Entry that tried to claim it.
        second: String,
    },
    /// A ring template is not a ring system that locants can number.
    BadRing {
        /// Template name.
        ring: String,
        /// What is wrong with it.
        reason: String,
    },
    /// Two ring templates have the same formula.
    DuplicateRing {
        /// The formula, written with one hydrogen.
        formula: String,
        /// Template that already has it.
        first: String,
        /// Template that tried to take it.
        second: String,
    },
}

impl fmt::Display for DictionaryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Smiles { entry, error } => write!(f, "entry {entry}: {error}"),
            Self::NoSpellings { entry } => write!(f, "entry {entry} has no spellings"),
            Self::BadSpelling { entry, spelling } => {
                write!(f, "entry {entry}: spelling {spelling:?} is not allowed")
            }
            Self::BadAttachment { entry } => write!(
                f,
                "entry {entry}: each '*' needs exactly one neighbour, and at most two are allowed"
            ),
            Self::DuplicateSpelling {
                spelling,
                first,
                second,
            } => write!(
                f,
                "spelling {spelling:?} is claimed by both {first} and {second}"
            ),
            Self::BadRing { ring, reason } => write!(f, "ring {ring}: {reason}"),
            Self::DuplicateRing {
                formula,
                first,
                second,
            } => write!(f, "rings {first} and {second} are both written {formula}"),
        }
    }
}

impl std::error::Error for DictionaryError {}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(name: &str, smiles: &str, spellings: &[&str], reversed: &[&str]) -> Entry {
        Entry::new(
            name,
            smiles,
            spellings.iter().map(ToString::to_string).collect(),
            reversed.iter().map(ToString::to_string).collect(),
        )
        .unwrap()
    }

    #[test]
    fn splits_attachment_points() {
        let carboxyl = entry("carboxyl", "*C(=O)O", &["COOH"], &["HOOC"]);
        assert_eq!(carboxyl.fragment().atoms().len(), 3);
        assert_eq!(
            carboxyl.attachments(),
            &[Attachment {
                atom: 0,
                order: BondOrder::Single
            }]
        );
        let linker = entry("ethylene", "*CC*", &["C2H4"], &[]);
        assert_eq!(linker.attachments().len(), 2);
    }

    #[test]
    fn rejects_bad_entries() {
        let bad = |smiles: &str, spellings: Vec<String>| Entry::new("x", smiles, spellings, vec![]);
        assert!(matches!(
            bad("*C(", vec!["X".into()]),
            Err(DictionaryError::Smiles { .. })
        ));
        assert!(matches!(
            bad("*C", vec![]),
            Err(DictionaryError::NoSpellings { .. })
        ));
        for spelling in ["X Y", "2C5H4N", "2-", "-C5H4N", "2-4C", "0-Py", "123-Py"] {
            assert!(
                matches!(
                    bad("*C", vec![spelling.into()]),
                    Err(DictionaryError::BadSpelling { .. })
                ),
                "{spelling}"
            );
        }
        assert!(bad("*c1ccccn1", vec!["2-C5H4N".into()]).is_ok());
        assert!(matches!(
            bad("*C(*)(*)", vec!["X".into()]),
            Err(DictionaryError::BadAttachment { .. })
        ));
    }

    #[test]
    fn finds_spellings_and_rejects_duplicates() {
        let mut dictionary = Dictionary::new("test", "1");
        dictionary
            .insert(entry("methyl", "*C", &["Me"], &[]))
            .unwrap();
        dictionary
            .insert(entry("methoxy", "*OC", &["OMe"], &["MeO"]))
            .unwrap();
        let found = dictionary.matches_at("CH3OMe", 3);
        assert_eq!(
            found,
            vec![SpellingMatch {
                entry: 1,
                len: 3,
                reversed: false
            }]
        );
        assert!(matches!(
            dictionary.insert(entry("other", "*C", &["Me"], &[])),
            Err(DictionaryError::DuplicateSpelling { .. })
        ));
        assert!(
            dictionary
                .insert(entry("methyl", "*C", &["Me", "CH3"], &[]))
                .unwrap()
                .is_some()
        );
    }
}
