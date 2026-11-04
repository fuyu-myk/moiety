use std::path::PathBuf;

use anyhow::{Context, Result, bail};

pub const USAGE: &str = "\
usage: moiety <command> [options] [arguments]

commands:
  parse <formula>...       print every reading of each formula
  write <smiles>...        write molecules as condensed formulas
  draw <formula>           render the skeletal structure as SVG
  html <formula>           render an interactive page with every reading
  conformance              run the conformance suite and round-trip corpora
  differential <results>   compare with another implementation's label results
  dictionary               list the abbreviation dictionary

options:
  --dictionary <file>      add the entries of a TOML dictionary to the default one
  --hypervalent            write nitro groups as N(=O)=O instead of [N+](=O)[O-]
  --reading <n>            draw reading n instead of reading 1
  --smiles                 draw a SMILES string instead of a formula
  -o, --output <file>      write output to a file instead of standard output
  --cases <file>           conformance suite (default conformance/cases.toml)
  --corpus <file>          round-trip SMILES list (default conformance/roundtrip.smi)
  --generated <n>          also round-trip n generated molecules (default 500)
  --seed <n>               seed for generated molecules (default 2026)
  --adjudications <file>   recorded verdicts (default conformance/differential/adjudications.toml)
  -h, --help               show this help
  -V, --version            show the version
";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Command {
    Parse,
    Write,
    Draw,
    Html,
    Conformance,
    Differential,
    Dictionary,
    Help,
    Version,
}

#[derive(Debug)]
pub struct Args {
    pub command: Command,
    pub inputs: Vec<String>,
    pub dictionary: Option<PathBuf>,
    pub hypervalent: bool,
    pub reading: usize,
    pub smiles: bool,
    pub output: Option<PathBuf>,
    pub cases: PathBuf,
    pub corpus: PathBuf,
    pub generated: usize,
    pub seed: u64,
    pub adjudications: PathBuf,
}

pub fn parse(argv: &[String]) -> Result<Args> {
    let mut iter = argv.iter();
    let command = match iter.next().map(String::as_str) {
        None | Some("-h" | "--help" | "help") => Command::Help,
        Some("-V" | "--version" | "version") => Command::Version,
        Some("parse") => Command::Parse,
        Some("write") => Command::Write,
        Some("draw") => Command::Draw,
        Some("html") => Command::Html,
        Some("conformance") => Command::Conformance,
        Some("differential") => Command::Differential,
        Some("dictionary") => Command::Dictionary,
        Some(other) => bail!("unknown command {other:?}\n\n{USAGE}"),
    };
    let mut args = Args {
        command,
        inputs: Vec::new(),
        dictionary: None,
        hypervalent: false,
        reading: 1,
        smiles: false,
        output: None,
        cases: PathBuf::from("conformance/cases.toml"),
        corpus: PathBuf::from("conformance/roundtrip.smi"),
        generated: 500,
        seed: 2026,
        adjudications: PathBuf::from("conformance/differential/adjudications.toml"),
    };
    while let Some(arg) = iter.next() {
        let mut value = |name: &str| {
            iter.next()
                .cloned()
                .with_context(|| format!("{name} needs a value"))
        };
        match arg.as_str() {
            "-h" | "--help" => args.command = Command::Help,
            "-V" | "--version" => args.command = Command::Version,
            "--dictionary" => args.dictionary = Some(value(arg)?.into()),
            "--hypervalent" => args.hypervalent = true,
            "--reading" => args.reading = number(&value(arg)?, arg)?,
            "--smiles" => args.smiles = true,
            "-o" | "--output" => args.output = Some(value(arg)?.into()),
            "--cases" => args.cases = value(arg)?.into(),
            "--corpus" => args.corpus = value(arg)?.into(),
            "--generated" => args.generated = number(&value(arg)?, arg)?,
            "--seed" => args.seed = number(&value(arg)?, arg)?,
            "--adjudications" => args.adjudications = value(arg)?.into(),
            "--" => args.inputs.extend(iter.by_ref().cloned()),
            flag if flag.starts_with("--") => bail!("unknown option {flag:?}"),
            input => args.inputs.push(input.to_string()),
        }
    }
    if matches!(
        args.command,
        Command::Parse | Command::Write | Command::Draw | Command::Html
    ) && args.inputs.is_empty()
    {
        bail!("this command needs at least one input\n\n{USAGE}");
    }
    if args.reading == 0 {
        bail!("--reading counts from 1");
    }
    Ok(args)
}

fn number<T: std::str::FromStr>(text: &str, flag: &str) -> Result<T> {
    text.parse()
        .ok()
        .with_context(|| format!("{flag} expects a number, got {text:?}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(list: &[&str]) -> Result<Args> {
        parse(&list.iter().map(ToString::to_string).collect::<Vec<_>>())
    }

    #[test]
    fn parses_commands_and_options() {
        let parsed = args(&["parse", "CH3COO-", "--hypervalent", "--", "-CH3"]).unwrap();
        assert_eq!(parsed.command, Command::Parse);
        assert_eq!(parsed.inputs, ["CH3COO-", "-CH3"]);
        assert!(parsed.hypervalent);
        let parsed = args(&["draw", "PhOH", "--reading", "2", "-o", "out.svg"]).unwrap();
        assert_eq!(parsed.reading, 2);
        assert_eq!(parsed.output, Some(PathBuf::from("out.svg")));
        assert_eq!(args(&[]).unwrap().command, Command::Help);
        assert_eq!(args(&["--version"]).unwrap().command, Command::Version);
        assert_eq!(args(&["parse", "-V"]).unwrap().command, Command::Version);
    }

    #[test]
    fn rejects_bad_arguments() {
        assert!(args(&["frobnicate"]).is_err());
        assert!(args(&["parse"]).is_err());
        assert!(args(&["parse", "CO", "--reading", "x"]).is_err());
        assert!(args(&["parse", "CO", "--reading", "0"]).is_err());
        assert!(args(&["parse", "CO", "--nope"]).is_err());
        assert!(args(&["parse", "CO", "--dictionary"]).is_err());
    }
}
