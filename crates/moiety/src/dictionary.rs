//! Loading abbreviation dictionaries and ring templates from TOML.

use std::sync::OnceLock;

use moiety_core::{Dictionary, DictionaryError, Entry, RingTemplate};
use serde::Deserialize;

const DEFAULT: &str = include_str!("../data/dictionary.toml");

/// A dictionary file that could not be loaded.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum LoadError {
    /// The file is not valid TOML or does not have the dictionary layout.
    #[error("dictionary file is malformed")]
    Toml(#[from] toml::de::Error),
    /// An entry or ring template could not be defined or added.
    #[error("dictionary entry or ring template rejected")]
    Entry(#[from] DictionaryError),
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct File {
    dictionary: Header,
    #[serde(default)]
    entry: Vec<EntryFile>,
    #[serde(default)]
    ring: Vec<RingFile>,
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

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RingFile {
    name: String,
    smiles: String,
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
/// Fails when the text is not a dictionary file or an entry or ring template is invalid.
pub fn load(text: &str) -> Result<Dictionary, LoadError> {
    let file: File = toml::from_str(text)?;
    let mut dictionary = Dictionary::new(&file.dictionary.name, &file.dictionary.version);
    add_file(&mut dictionary, file)?;
    Ok(dictionary)
}

/// Adds the entries and ring templates of a TOML dictionary to `base`, replacing
/// those with the same name. A ring template also replaces one of `base` with
/// the same formula, so a dictionary can bring pyrazole in place of the bundled
/// imidazole.
///
/// The result keeps the name and version of `base`.
///
/// # Errors
///
/// Fails when the text is not a dictionary file or an entry or ring template is invalid.
pub fn extend(base: &Dictionary, text: &str) -> Result<Dictionary, LoadError> {
    let file: File = toml::from_str(text)?;
    let mut dictionary = base.clone();
    add_file(&mut dictionary, file)?;
    Ok(dictionary)
}

fn add_file(dictionary: &mut Dictionary, file: File) -> Result<(), LoadError> {
    for entry in file.entry {
        dictionary.insert(Entry::new(
            &entry.name,
            &entry.smiles,
            entry.spellings,
            entry.reversed,
        )?)?;
    }
    let mut added: Vec<String> = Vec::new();
    for ring in file.ring {
        let ring = RingTemplate::new(&ring.name, &ring.smiles)?;
        let replaced = dictionary
            .rings()
            .iter()
            .find(|other| {
                other.formula() == ring.formula() && !added.contains(&other.name().to_string())
            })
            .map(|other| other.name().to_string());
        if let Some(name) = replaced {
            dictionary.remove_ring(&name);
        }
        added.push(ring.name().to_string());
        dictionary.insert_ring(ring)?;
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
    fn rejects_bad_ring_templates() {
        let header = "[dictionary]\nname = \"x\"\nversion = \"1\"\n";
        let ring = |name: &str, smiles: &str| {
            format!("[[ring]]\nname = \"{name}\"\nsmiles = \"{smiles}\"\n")
        };
        let unparsed = format!("{header}{}", ring("broken", "c1cc("));
        assert!(matches!(
            load(&unparsed),
            Err(LoadError::Entry(DictionaryError::BadRing { ring, .. })) if ring == "broken"
        ));
        let isomers = format!(
            "{header}{}{}",
            ring("pyrimidine", "n1cnccc1"),
            ring("pyrazine", "n1ccncc1")
        );
        let Err(LoadError::Entry(error)) = load(&isomers) else {
            panic!("isomeric templates share a formula");
        };
        assert_eq!(
            error.to_string(),
            "rings pyrimidine and pyrazine share the formula C4H4N2"
        );
    }

    #[test]
    fn ring_templates_replace_bundled_ones_with_their_formula() {
        let pyrazole = "[dictionary]\nname = \"x\"\nversion = \"1\"\n[[ring]]\nname = \"pyrazole\"\nsmiles = \"[nH]1nccc1\"\n";
        let dictionary = extend(default_dictionary(), pyrazole).unwrap();
        let names: Vec<&str> = dictionary.rings().iter().map(|r| r.name()).collect();
        assert!(names.contains(&"pyrazole"));
        assert!(!names.contains(&"imidazole"));
        let benzene = pyrazole.replace("[nH]1nccc1", "c1ccccc1");
        assert!(matches!(
            extend(default_dictionary(), &benzene),
            Err(LoadError::Entry(DictionaryError::DuplicateRing { .. }))
        ));
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
