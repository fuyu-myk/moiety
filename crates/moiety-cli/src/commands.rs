use std::io::{IsTerminal, Write};
use std::process::ExitCode;

use anyhow::{Context, Result, bail};
use moiety::conformance::{
    Report, RoundTrip, Thresholds, gate, generated_corpus, load_suite, read_smiles_list,
    round_trip, run as run_suite,
};
use moiety::writer::condense;
use moiety::{
    ChargeConvention, Dictionary, Interpretation, Molecule, Options, Reading, default_dictionary,
    smiles,
};
use moiety_draw::terminal::{TerminalOptions, Theme, detect, detect_theme, render_terminal};
use moiety_draw::{
    RasterOptions, SvgOptions, render_html, render_png, render_reading_svg, render_svg,
};

use crate::args::{Args, Command, Format, USAGE};

pub fn run(args: Args) -> Result<ExitCode> {
    let dictionary = load_dictionary(&args)?;
    let options = Options {
        charge_convention: if args.hypervalent {
            ChargeConvention::Hypervalent
        } else {
            ChargeConvention::ChargeSeparated
        },
        ..Options::default()
    };
    match args.command {
        Command::Help => {
            print!("{USAGE}");
            Ok(ExitCode::SUCCESS)
        }
        Command::Version => {
            println!("moiety {}", env!("CARGO_PKG_VERSION"));
            Ok(ExitCode::SUCCESS)
        }
        Command::Parse => parse(&args, &dictionary, &options),
        Command::Write => write(&args, &dictionary),
        Command::Draw | Command::Html => draw(&args, &dictionary, &options),
        Command::Conformance => conformance(&args, &dictionary, &options),
        Command::Differential => differential(&args, &dictionary, &options),
        Command::Dictionary => {
            emit(&args, &list_dictionary(&dictionary))?;
            Ok(ExitCode::SUCCESS)
        }
    }
}

fn load_dictionary(args: &Args) -> Result<Dictionary> {
    let Some(path) = &args.dictionary else {
        return Ok(default_dictionary().clone());
    };
    let text =
        std::fs::read_to_string(path).with_context(|| format!("cannot read {}", path.display()))?;
    moiety::dictionary::extend(default_dictionary(), &text)
        .with_context(|| format!("in {}", path.display()))
}

fn emit_bytes(args: &Args, bytes: &[u8]) -> Result<()> {
    match &args.output {
        Some(path) => {
            std::fs::write(path, bytes).with_context(|| format!("cannot write {}", path.display()))
        }
        None => {
            let mut stdout = std::io::stdout().lock();
            stdout
                .write_all(bytes)
                .and_then(|()| stdout.flush())
                .context("cannot write to standard output")
        }
    }
}

fn emit(args: &Args, text: &str) -> Result<()> {
    emit_bytes(args, text.as_bytes())
}

fn parse(args: &Args, dictionary: &Dictionary, options: &Options) -> Result<ExitCode> {
    let mut out = String::new();
    let mut all_read = true;
    for input in &args.inputs {
        let result = moiety::interpret_with(input, dictionary, options);
        all_read &= !matches!(result, Interpretation::Invalid { .. });
        out.push_str(&describe(input, &result));
    }
    emit(args, &out)?;
    Ok(if all_read {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    })
}

/// Plain-text description of an interpretation, as printed by `moiety parse`.
pub fn describe(input: &str, result: &Interpretation) -> String {
    let mut out = format!("{input}\n");
    match result {
        Interpretation::Unique(reading) => {
            line(&mut out, &format!("  unique: {}", reading_summary(reading)));
            readings_detail(&mut out, input, std::slice::from_ref(reading));
        }
        Interpretation::Ambiguous { readings, classes } => {
            let names: Vec<&str> = classes.iter().map(|c| c.as_str()).collect();
            line(
                &mut out,
                &format!(
                    "  ambiguous: {} readings ({})",
                    readings.len(),
                    names.join(", ")
                ),
            );
            for reading in readings {
                line(
                    &mut out,
                    &format!("  {}. {}", reading.rank, reading_summary(reading)),
                );
            }
            readings_detail(&mut out, input, readings);
        }
        Interpretation::Underdetermined { span, reason } => {
            let text = input.get(span.start..span.end).unwrap_or_default();
            line(
                &mut out,
                &format!("  underdetermined: {reason} at {span} ({text})"),
            );
        }
        Interpretation::Invalid { diagnostics } => {
            for diagnostic in diagnostics {
                for text in diagnostic.render(input).lines() {
                    line(&mut out, &format!("  {text}"));
                }
            }
        }
    }
    out
}

fn line(out: &mut String, text: &str) {
    out.push_str(text);
    out.push('\n');
}

fn reading_summary(reading: &Reading) -> String {
    format!(
        "{}  (score {}, {})",
        reading.smiles,
        reading.score,
        reading.direction.as_str()
    )
}

fn readings_detail(out: &mut String, input: &str, readings: &[Reading]) {
    for reading in readings {
        let prefix = if readings.len() > 1 {
            format!("     {}: ", reading.rank)
        } else {
            "     ".to_string()
        };
        for ambiguity in &reading.ambiguities {
            let text = input
                .get(ambiguity.span.start..ambiguity.span.end)
                .unwrap_or_default();
            line(
                out,
                &format!(
                    "{prefix}{} at {text:?}: {}",
                    ambiguity.class, ambiguity.choice
                ),
            );
        }
        for warning in &reading.warnings {
            for text in warning.render(input).lines() {
                line(out, &format!("     {text}"));
            }
        }
    }
}

fn write(args: &Args, dictionary: &Dictionary) -> Result<ExitCode> {
    let mut out = String::new();
    let mut all_written = true;
    for input in &args.inputs {
        let molecule =
            smiles::parse(input).with_context(|| format!("cannot read SMILES {input:?}"))?;
        match condense(&molecule, dictionary) {
            Ok(formula) => line(&mut out, &format!("{input}\t{formula}")),
            Err(error) => {
                all_written = false;
                line(&mut out, &format!("{input}\terror: {error}"));
            }
        }
    }
    emit(args, &out)?;
    Ok(if all_written {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    })
}

fn draw(args: &Args, dictionary: &Dictionary, options: &Options) -> Result<ExitCode> {
    let [input] = &args.inputs[..] else {
        bail!(
            "{} draws one formula at a time",
            if args.command == Command::Html {
                "html"
            } else {
                "draw"
            }
        );
    };
    let svg_options = SvgOptions {
        title: Some(input.clone()),
        ..SvgOptions::default()
    };
    if args.smiles {
        if args.command == Command::Html {
            bail!("html shows the readings of a formula; draw a SMILES string with draw --smiles");
        }
        let molecule =
            smiles::parse(input).with_context(|| format!("cannot read SMILES {input:?}"))?;
        return show(args, &molecule, || render_svg(&molecule, &svg_options));
    }
    let result = moiety::interpret_with(input, dictionary, options);
    let readings = result.readings();
    if readings.is_empty() {
        eprint!("{}", describe(input, &result));
        return Ok(ExitCode::FAILURE);
    }
    if args.command == Command::Html {
        emit(args, &render_html(input, readings, &svg_options))?;
        return Ok(ExitCode::SUCCESS);
    }
    let reading = readings.get(args.reading - 1).with_context(|| {
        format!(
            "{input} has {} reading(s), not {}",
            readings.len(),
            args.reading
        )
    })?;
    show(args, &reading.molecule, || {
        render_reading_svg(reading, &svg_options)
    })
}

/// Writes the drawing of `molecule` in the format `args` asks for, or the one
/// that suits where it goes, using `svg` for SVG.
fn show(args: &Args, molecule: &Molecule, svg: impl FnOnce() -> String) -> Result<ExitCode> {
    let getenv = |key: &str| std::env::var(key).ok();
    let terminal = args.output.is_none() && std::io::stdout().is_terminal();
    match output_format(args, terminal, getenv) {
        Format::Svg => emit(args, &svg())?,
        Format::Png => {
            if terminal {
                bail!("PNG output needs -o <file> or a pipe");
            }
            let theme = args.theme.unwrap_or(Theme::Light);
            emit_bytes(args, &render_png(molecule, &png_options(theme)))?;
        }
        Format::Terminal(protocol) => {
            let options = TerminalOptions {
                protocol,
                columns: terminal_columns(getenv),
                theme: args.theme.unwrap_or_else(|| detect_theme(getenv)),
            };
            emit(args, &render_terminal(molecule, &options))?;
        }
    }
    Ok(ExitCode::SUCCESS)
}

/// The format `args` asks for, or else PNG for a `.png` output file, the
/// terminal's own protocol when writing to a `terminal`, and SVG otherwise.
fn output_format(args: &Args, terminal: bool, getenv: impl Fn(&str) -> Option<String>) -> Format {
    if let Some(format) = args.format {
        return format;
    }
    let png = args.output.as_ref().is_some_and(|path| {
        path.extension()
            .is_some_and(|ext| ext.eq_ignore_ascii_case("png"))
    });
    if png {
        Format::Png
    } else if terminal {
        Format::Terminal(detect(getenv))
    } else {
        Format::Svg
    }
}

/// Sizes twice those of SVG, so that PNG files stay sharp on dense screens,
/// in the colours of `theme` on a transparent background.
fn png_options(theme: Theme) -> RasterOptions {
    let svg = SvgOptions::default();
    RasterOptions {
        bond_length: 2.0 * svg.bond_length,
        font_size: 2.0 * svg.font_size,
        line_width: 2.0 * svg.line_width,
        padding: 2.0 * svg.padding,
        ink: if theme == Theme::Dark {
            [255, 255, 255, 255]
        } else {
            [0, 0, 0, 255]
        },
        background: [0, 0, 0, 0],
    }
}

/// Columns of the terminal: `COLUMNS` if set, else what `stty` reports for the
/// controlling terminal, else 80.
fn terminal_columns(getenv: impl Fn(&str) -> Option<String>) -> usize {
    let stty = || {
        let tty = std::fs::File::open("/dev/tty").ok()?;
        let out = std::process::Command::new("stty")
            .arg("size")
            .stdin(tty)
            .output()
            .ok()?;
        String::from_utf8(out.stdout).ok()
    };
    getenv("COLUMNS")
        .and_then(|columns| columns.trim().parse().ok())
        .or_else(|| stty().and_then(|size| size.split_whitespace().nth(1)?.parse().ok()))
        .filter(|&columns: &usize| columns > 0)
        .unwrap_or(80)
}

fn conformance(args: &Args, dictionary: &Dictionary, options: &Options) -> Result<ExitCode> {
    let text = std::fs::read_to_string(&args.cases)
        .with_context(|| format!("cannot read {}", args.cases.display()))?;
    let suite = load_suite(&text).with_context(|| format!("in {}", args.cases.display()))?;
    let report = run_suite(&suite, dictionary, options);
    let corpus_text = std::fs::read_to_string(&args.corpus)
        .with_context(|| format!("cannot read {}", args.corpus.display()))?;
    let mut molecules = read_smiles_list(&corpus_text).map_err(|(line, error)| {
        anyhow::anyhow!("{} line {line}: {error}", args.corpus.display())
    })?;
    let curated = molecules.len();
    molecules.extend(generated_corpus(args.generated, args.seed));
    let trip = round_trip(&molecules, dictionary, options);
    let misses = gate(&report.metrics, Some(&trip), &Thresholds::default());
    let mut out = conformance_report(
        &suite.suite.name,
        &suite.suite.version,
        &report,
        &trip,
        curated,
        args.generated,
    );
    if misses.is_empty() {
        line(&mut out, "gate: pass");
    } else {
        line(&mut out, "gate: FAIL");
        for miss in &misses {
            line(&mut out, &format!("  {miss}"));
        }
    }
    emit(args, &out)?;
    Ok(if misses.is_empty() {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    })
}

fn differential(args: &Args, dictionary: &Dictionary, options: &Options) -> Result<ExitCode> {
    use moiety::differential::{compare, load_adjudications, load_external, report, summary};
    let [path] = &args.inputs[..] else {
        bail!("differential compares one results file at a time");
    };
    let text = std::fs::read_to_string(path).with_context(|| format!("cannot read {path}"))?;
    let external = load_external(&text).with_context(|| format!("in {path}"))?;
    let verdicts = match std::fs::read_to_string(&args.adjudications) {
        Ok(text) => load_adjudications(&text)
            .with_context(|| format!("in {}", args.adjudications.display()))?,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Vec::new(),
        Err(error) => {
            return Err(error)
                .with_context(|| format!("cannot read {}", args.adjudications.display()));
        }
    };
    let comparisons = compare(&external, &verdicts, dictionary, options);
    let (counts, open) = summary(&comparisons);
    let mut out = format!(
        "moiety {} against {}\n\n",
        env!("CARGO_PKG_VERSION"),
        external.source
    );
    for (agreement, count) in &counts {
        line(&mut out, &format!("{:<16}{count}", agreement.as_str()));
    }
    line(&mut out, &format!("{:<16}{open}", "unadjudicated"));
    for comparison in comparisons.iter().filter(|c| c.unadjudicated()) {
        line(
            &mut out,
            &format!(
                "  {:<14} {:<10} moiety {:?} other {}",
                comparison.input,
                comparison.agreement.as_str(),
                comparison.readings,
                comparison.other.as_deref().unwrap_or("-")
            ),
        );
    }
    match &args.output {
        Some(path) => {
            std::fs::write(path, report(&external, &comparisons))
                .with_context(|| format!("cannot write {}", path.display()))?;
            print!("{out}");
        }
        None => print!("{out}"),
    }
    Ok(if open == 0 {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    })
}

fn conformance_report(
    name: &str,
    version: &str,
    report: &Report,
    trip: &RoundTrip,
    curated: usize,
    generated: usize,
) -> String {
    let m = &report.metrics;
    let mut out = format!(
        "{name} {version}: grammar {}, ranking policy {}\n\n",
        moiety::GRAMMAR_VERSION,
        moiety::RANKING_POLICY_VERSION
    );
    for (label, value) in [
        ("cases passed", m.passed.to_string()),
        ("false-unique rate", m.false_unique.to_string()),
        ("ambiguity class recall", m.class_recall.to_string()),
        ("top-1 accuracy, ambiguous cases", m.top1.to_string()),
        ("diagnostic span accuracy", m.span_accuracy.to_string()),
        ("unique accuracy", m.unique_accuracy.to_string()),
        ("ambiguity rate", m.ambiguity_rate.to_string()),
        (
            "round-trip recovery",
            format!(
                "{} ({curated} curated, {generated} generated)",
                trip.recovered
            ),
        ),
        ("round-trip read back ambiguous", trip.ambiguous.to_string()),
        ("round-trip unwritable", trip.unwritable.to_string()),
    ] {
        line(&mut out, &format!("{label:<34}{value}"));
    }
    let failing: Vec<_> = report.outcomes.iter().filter(|o| !o.passed()).collect();
    if !failing.is_empty() {
        line(&mut out, "\nfailing cases:");
        for outcome in failing {
            line(
                &mut out,
                &format!(
                    "  {} {:?}: {}",
                    outcome.id,
                    outcome.input,
                    outcome.problems.join("; ")
                ),
            );
        }
    }
    if !trip.misses.is_empty() {
        line(&mut out, "\nround-trip misses:");
        for miss in &trip.misses {
            line(
                &mut out,
                &format!(
                    "  {} written {} read {}",
                    miss.smiles, miss.written, miss.read
                ),
            );
        }
    }
    out.push('\n');
    out
}

fn list_dictionary(dictionary: &Dictionary) -> String {
    let mut out = format!("{} {}\n", dictionary.name(), dictionary.version());
    for entry in dictionary.entries() {
        let mut spellings = entry.spellings().join(" ");
        if !entry.reversed().is_empty() {
            spellings.push_str(&format!(" | reversed: {}", entry.reversed().join(" ")));
        }
        line(
            &mut out,
            &format!("  {:<28}{:<36}{spellings}", entry.name(), entry.smiles()),
        );
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use moiety::interpret;

    #[test]
    fn draws_and_pages_a_formula() {
        let out = std::env::temp_dir().join("moiety-draws-a-formula.svg");
        let out = out.to_str().unwrap().to_string();
        let args =
            crate::args::parse(&["draw".into(), "PhCH2OH".into(), "-o".into(), out]).unwrap();
        let options = Options::default();
        assert!(
            matches!(draw(&args, default_dictionary(), &options), Ok(code) if code == ExitCode::SUCCESS)
        );
        let too_far = crate::args::parse(&[
            "draw".into(),
            "PhCH2OH".into(),
            "--reading".into(),
            "2".into(),
        ])
        .unwrap();
        assert!(draw(&too_far, default_dictionary(), &options).is_err());
        let invalid = crate::args::parse(&["html".to_string(), "CH3CH3OH".to_string()]).unwrap();
        assert!(
            matches!(draw(&invalid, default_dictionary(), &options), Ok(code) if code == ExitCode::FAILURE)
        );
    }

    #[test]
    fn draws_smiles() {
        let options = Options::default();
        let run = |list: &[&str]| {
            let argv: Vec<String> = list.iter().map(ToString::to_string).collect();
            draw(
                &crate::args::parse(&argv).unwrap(),
                default_dictionary(),
                &options,
            )
        };
        let out = std::env::temp_dir().join("moiety-draws-smiles.svg");
        let out = out.to_str().unwrap();
        assert!(
            matches!(run(&["draw", "--smiles", "c1ccccc1O", "-o", out]), Ok(code) if code == ExitCode::SUCCESS)
        );
        assert!(std::fs::read_to_string(out).unwrap().contains("<svg"));
        assert!(run(&["draw", "--smiles", "C1CC"]).is_err());
        assert!(run(&["html", "--smiles", "CCO"]).is_err());
    }

    #[test]
    fn draws_png_files_and_text() {
        let options = Options::default();
        let run = |list: &[&str]| {
            let argv: Vec<String> = list.iter().map(ToString::to_string).collect();
            draw(
                &crate::args::parse(&argv).unwrap(),
                default_dictionary(),
                &options,
            )
        };
        let png = std::env::temp_dir().join("moiety-draws-png.png");
        let png = png.to_str().unwrap();
        assert!(matches!(run(&["draw", "PhOH", "-o", png]), Ok(code) if code == ExitCode::SUCCESS));
        assert!(
            std::fs::read(png)
                .unwrap()
                .starts_with(b"\x89PNG\r\n\x1a\n")
        );
        let text = std::env::temp_dir().join("moiety-draws-text.txt");
        let text = text.to_str().unwrap();
        assert!(run(&["draw", "CH3CH2OH", "--format", "text", "-o", text]).is_ok());
        let drawing = std::fs::read_to_string(text).unwrap();
        assert!(
            drawing.contains("OH") || drawing.contains("HO"),
            "{drawing}"
        );
    }

    #[test]
    fn picks_the_format_that_suits_the_output() {
        let parse = |list: &[&str]| {
            let argv: Vec<String> = list.iter().map(ToString::to_string).collect();
            crate::args::parse(&argv).unwrap()
        };
        let kitty = |key: &str| (key == "TERM").then(|| "xterm-kitty".to_string());
        let plain = parse(&["draw", "CO"]);
        assert_eq!(output_format(&plain, false, kitty), Format::Svg);
        assert_eq!(
            output_format(&plain, true, kitty),
            Format::Terminal(moiety_draw::terminal::Protocol::Kitty)
        );
        let file = parse(&["draw", "CO", "-o", "out.PNG"]);
        assert_eq!(output_format(&file, false, kitty), Format::Png);
        assert_eq!(
            output_format(&parse(&["draw", "CO", "-o", "out.svg"]), false, kitty),
            Format::Svg
        );
        let forced = parse(&["draw", "CO", "--format", "sixel", "-o", "out.png"]);
        assert_eq!(
            output_format(&forced, false, kitty),
            Format::Terminal(moiety_draw::terminal::Protocol::Sixel)
        );
    }

    #[test]
    fn reads_the_terminal_width_from_columns() {
        assert_eq!(terminal_columns(|_| Some("132".to_string())), 132);
        let columns = terminal_columns(|_| Some("0".to_string()));
        assert!(columns > 0);
    }

    #[test]
    fn describes_each_result_kind() {
        let unique = describe("CH3CH2OH", &interpret("CH3CH2OH"));
        assert!(unique.contains("unique: CCO"));
        let ambiguous = describe("COCl2", &interpret("COCl2"));
        assert!(ambiguous.contains("ambiguous: 2 readings (symbol-segmentation)"));
        assert!(ambiguous.contains("1. O=C(Cl)Cl"));
        let invalid = describe("CH3CH3OH", &interpret("CH3CH3OH"));
        assert!(invalid.contains("error[CF0204]"));
        let under = describe("C6H12O6", &interpret("C6H12O6"));
        assert!(under.contains("underdetermined: molecular-formula at 0..7 (C6H12O6)"));
    }
}
