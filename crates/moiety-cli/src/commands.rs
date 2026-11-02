use std::process::ExitCode;

use anyhow::{Context, Result, bail};
use moiety::conformance::{
    Report, RoundTrip, Thresholds, gate, generated_corpus, load_suite, read_smiles_list,
    round_trip, run as run_suite,
};
use moiety::writer::condense;
use moiety::{
    ChargeConvention, Dictionary, Interpretation, Options, Reading, default_dictionary, smiles,
};
use moiety_draw::{SvgOptions, render_html, render_reading_svg};

use crate::args::{Args, Command, USAGE};

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
        Command::Parse => parse(&args, &dictionary, &options),
        Command::Write => write(&args, &dictionary),
        Command::Draw | Command::Html => draw(&args, &dictionary, &options),
        Command::Conformance => conformance(&args, &dictionary, &options),
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

fn emit(args: &Args, text: &str) -> Result<()> {
    match &args.output {
        Some(path) => {
            std::fs::write(path, text).with_context(|| format!("cannot write {}", path.display()))
        }
        None => {
            print!("{text}");
            Ok(())
        }
    }
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
    let result = moiety::interpret_with(input, dictionary, options);
    let readings = result.readings();
    if readings.is_empty() {
        eprint!("{}", describe(input, &result));
        return Ok(ExitCode::FAILURE);
    }
    let svg_options = SvgOptions {
        title: Some(input.clone()),
        ..SvgOptions::default()
    };
    let output = if args.command == Command::Html {
        render_html(input, readings, &svg_options)
    } else {
        let reading = readings.get(args.reading - 1).with_context(|| {
            format!(
                "{input} has {} reading(s), not {}",
                readings.len(),
                args.reading
            )
        })?;
        render_reading_svg(reading, &svg_options)
    };
    emit(args, &output)?;
    Ok(ExitCode::SUCCESS)
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
        let args = crate::args::parse(&["draw".to_string(), "PhCH2OH".to_string()]).unwrap();
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
    fn describes_each_result_kind() {
        let unique = describe("CH3CH2OH", &interpret("CH3CH2OH"));
        assert!(unique.contains("unique: CCO"));
        let ambiguous = describe("COCl2", &interpret("COCl2"));
        assert!(ambiguous.contains("ambiguous: 2 readings (symbol-segmentation)"));
        assert!(ambiguous.contains("1. O=C(Cl)Cl"));
        let invalid = describe("CH3CH3OH", &interpret("CH3CH3OH"));
        assert!(invalid.contains("error[CF0204]"));
        let under = describe("C6H6", &interpret("C6H6"));
        assert!(under.contains("underdetermined: molecular-formula at 0..4 (C6H6)"));
    }
}
