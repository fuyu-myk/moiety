//! Loading abbreviation dictionaries from TOML.

use std::sync::OnceLock;

use moiety_core::{Dictionary, DictionaryError, Entry};
use serde::Deserialize;

const DEFAULT: &str = include_str!("../data/dictionary.toml");

/// A dictionary file that could not be loaded.
#[derive(Debug, thiserror::Error)]
pub enum LoadError {
    /// The file is not valid TOML or does not have the dictionary layout.
    #[error("dictionary file is malformed: {0}")]
    Toml(#[from] toml::de::Error),
    /// An entry could not be defined or added.
    #[error("dictionary entry rejected: {0}")]
    Entry(#[from] DictionaryError),
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct File {
    dictionary: Header,
    #[serde(default)]
    entry: Vec<EntryFile>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Header {
    name: String,
    version: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct EntryFile {
    name: String,
    smiles: String,
    #[serde(default)]
    spellings: Vec<String>,
    #[serde(default)]
    reversed: Vec<String>,
}

/// The dictionary that ships with moiety.
///
/// # Panics
///
/// Panics only if the bundled dictionary file is broken, which the test suite rules out.
#[must_use]
pub fn default_dictionary() -> &'static Dictionary {
    static DICTIONARY: OnceLock<Dictionary> = OnceLock::new();
    DICTIONARY.get_or_init(|| load(DEFAULT).expect("the bundled dictionary is valid"))
}

/// Reads a dictionary from TOML text.
///
/// # Errors
///
/// Fails when the text is not a dictionary file or an entry is invalid.
pub fn load(text: &str) -> Result<Dictionary, LoadError> {
    let file: File = toml::from_str(text)?;
    let mut dictionary = Dictionary::new(&file.dictionary.name, &file.dictionary.version);
    add_entries(&mut dictionary, file.entry)?;
    Ok(dictionary)
}

/// Adds the entries of a TOML dictionary to `base`, replacing entries with the same name.
///
/// The result keeps the name and version of `base`.
///
/// # Errors
///
/// Fails when the text is not a dictionary file or an entry is invalid.
pub fn extend(base: &Dictionary, text: &str) -> Result<Dictionary, LoadError> {
    let file: File = toml::from_str(text)?;
    let mut dictionary = base.clone();
    add_entries(&mut dictionary, file.entry)?;
    Ok(dictionary)
}

fn add_entries(dictionary: &mut Dictionary, entries: Vec<EntryFile>) -> Result<(), LoadError> {
    for entry in entries {
        dictionary.insert(Entry::new(
            &entry.name,
            &entry.smiles,
            entry.spellings,
            entry.reversed,
        )?)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bundled_dictionary_loads() {
        let dictionary = default_dictionary();
        assert_eq!(dictionary.name(), "moiety-default");
        assert!(dictionary.get("phenyl").is_some());
        assert!(dictionary.entries().len() > 40);
    }

    #[test]
    fn extends_and_overrides() {
        let extra = r#"
            [dictionary]
            name = "lab"
            version = "1"

            [[entry]]
            name = "methyl"
            smiles = "*C"
            spellings = ["Me", "Mt"]

            [[entry]]
            name = "pyridin-2-yl"
            smiles = "*c1ccccn1"
            spellings = ["2Pyr"]
        "#;
        assert!(matches!(
            extend(default_dictionary(), extra),
            Err(LoadError::Entry(_))
        ));
        let fixed = extra.replace("2Pyr", "2-Pyr");
        let dictionary = extend(default_dictionary(), &fixed).unwrap();
        assert_eq!(dictionary.name(), "moiety-default");
        assert_eq!(dictionary.get("methyl").unwrap().spellings(), ["Me", "Mt"]);
        assert!(dictionary.get("pyridin-2-yl").is_some());
    }

    #[test]
    fn rejects_malformed_files() {
        assert!(matches!(
            load("[dictionary]\nname = 1"),
            Err(LoadError::Toml(_))
        ));
        let bad = "[dictionary]\nname = \"x\"\nversion = \"1\"\n[[entry]]\nname = \"x\"\nsmiles = \"*C(\"\nspellings = [\"X\"]\n";
        assert!(matches!(
            load(bad),
            Err(LoadError::Entry(DictionaryError::Smiles { .. }))
        ));
    }
}
