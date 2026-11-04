use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use moiety_draw::terminal::{Protocol, Theme};

pub const USAGE: &str = "\
usage: moiety <command> [options] [arguments]

commands:
  parse <formula>...       print every reading of each formula
  write <smiles>...        write molecules as condensed formulas
  draw <formula>           draw the skeletal structure, in the terminal or as SVG
  html <formula>           render an interactive page with every reading
  conformance              run the conformance suite and round-trip corpora
  differential <results>   compare with another implementation's label results
  dictionary               list the abbreviation dictionary

options:
  --dictionary <file>      add the entries of a TOML dictionary to the default one
  --hypervalent            write nitro groups as N(=O)=O instead of [N+](=O)[O-]
  --reading <n>            draw reading n instead of reading 1
  --smiles                 draw a SMILES string instead of a formula
  --format <format>        draw as svg, png, text, kitty, iterm2 or sixel; by default
                           draw shows the structure when writing to a terminal,
                           writes PNG for -o <file>.png, and SVG otherwise
  --theme <dark|light>     colours of terminal images and PNG files
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

/// What `draw` writes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Svg,
    Png,
    Terminal(Protocol),
}

#[derive(Debug)]
pub struct Args {
    pub command: Command,
    pub inputs: Vec<String>,
    pub dictionary: Option<PathBuf>,
    pub hypervalent: bool,
    pub reading: usize,
    pub smiles: bool,
    pub format: Option<Format>,
    pub theme: Option<Theme>,
    pub output: Option<PathBuf>,
    pub cases: PathBuf,
    pub corpus: PathBuf,
    pub generated: usize,
    pub seed: u64,
    pub adjudications: PathBuf,
}

impl Args {
    /// The defaults of every option, for `command`.
    fn new(command: Command) -> Self {
        Self {
            command,
            inputs: Vec::new(),
            dictionary: None,
            hypervalent: false,
            reading: 1,
            smiles: false,
            format: None,
            theme: None,
            output: None,
            cases: PathBuf::from("conformance/cases.toml"),
            corpus: PathBuf::from("conformance/roundtrip.smi"),
            generated: 500,
            seed: 2026,
            adjudications: PathBuf::from("conformance/differential/adjudications.toml"),
        }
    }
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
    let mut args = Args::new(command);
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
            "--format" => args.format = Some(format(&value(arg)?)?),
            "--theme" => args.theme = Some(theme(&value(arg)?)?),
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
    check(&args)?;
    Ok(args)
}

/// Fails on options that do not fit together or with the command.
fn check(args: &Args) -> Result<()> {
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
    if args.smiles && !matches!(args.command, Command::Draw | Command::Html) {
        bail!("--smiles applies to draw");
    }
    if (args.format.is_some() || args.theme.is_some()) && args.command != Command::Draw {
        bail!("--format and --theme apply to draw");
    }
    Ok(())
}

fn format(name: &str) -> Result<Format> {
    Ok(match name {
        "svg" => Format::Svg,
        "png" => Format::Png,
        "text" => Format::Terminal(Protocol::Unicode),
        "kitty" => Format::Terminal(Protocol::Kitty),
        "iterm2" => Format::Terminal(Protocol::Iterm2),
        "sixel" => Format::Terminal(Protocol::Sixel),
        _ => bail!("--format is one of svg, png, text, kitty, iterm2 or sixel, not {name:?}"),
    })
}

fn theme(name: &str) -> Result<Theme> {
    Ok(match name {
        "dark" => Theme::Dark,
        "light" => Theme::Light,
        _ => bail!("--theme is dark or light, not {name:?}"),
    })
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
        let parsed = args(&["draw", "CO", "--format", "kitty", "--theme", "light"]).unwrap();
        assert_eq!(parsed.format, Some(Format::Terminal(Protocol::Kitty)));
        assert_eq!(parsed.theme, Some(Theme::Light));
        assert_eq!(
            args(&["draw", "CO", "--format", "text"]).unwrap().format,
            Some(Format::Terminal(Protocol::Unicode))
        );
    }

    #[test]
    fn rejects_bad_arguments() {
        assert!(args(&["frobnicate"]).is_err());
        assert!(args(&["parse"]).is_err());
        assert!(args(&["parse", "CO", "--reading", "x"]).is_err());
        assert!(args(&["parse", "CO", "--reading", "0"]).is_err());
        assert!(args(&["parse", "CO", "--nope"]).is_err());
        assert!(args(&["parse", "CO", "--dictionary"]).is_err());
        assert!(args(&["parse", "CO", "--smiles"]).is_err());
        assert!(args(&["draw", "CO", "--format", "bmp"]).is_err());
        assert!(args(&["draw", "CO", "--theme", "blue"]).is_err());
        assert!(args(&["html", "CO", "--format", "png"]).is_err());
        assert!(args(&["parse", "CO", "--theme", "dark"]).is_err());
    }
}
