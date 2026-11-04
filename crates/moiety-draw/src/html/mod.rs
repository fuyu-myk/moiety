//! A self-contained HTML page that traces every reading of a formula back to its characters.
//!
//! The page shows the formula with one `<span data-offset>` per character and,
//! below it, the drawing of every reading with its rank and SMILES. Pointing at
//! or focusing an atom or bond highlights the characters it came from; pointing
//! at, focusing or selecting characters highlights every atom and bond whose
//! span overlaps them, in all readings at once. Styles and script are inline,
//! nothing is loaded from elsewhere, and colours follow `prefers-color-scheme`.
//!
//! [`render_fragment`] renders the same body for any interpretation, for pages
//! that embed the viewer and replace its content as the formula changes, and
//! [`render_smiles_fragment`] a body for a molecule given as SMILES.

use moiety_core::{Diagnostic, Interpretation, Reading, Span, Underdetermination};

use crate::svg::{SvgOptions, render_reading_svg, render_svg};
use crate::xml::escape;

/// The viewer's stylesheet.
pub const STYLE: &str = include_str!("viewer.css");

/// The viewer's script. It also traces content added to the page after it runs.
pub const SCRIPT: &str = include_str!("viewer.js");

/// Renders the viewer page for `input` and its `readings`, drawn with `options`.
///
/// Each drawing gets the title "Reading N: SMILES" in place of `options.title`.
/// An empty `readings` slice gives a page that says the formula has no reading.
#[must_use]
pub fn render_html(input: &str, readings: &[Reading], options: &SvgOptions) -> String {
    let mut page =
        String::from("<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n");
    page.push_str("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n");
    page.push_str(&format!("<title>{} (moiety)</title>\n", escape(input)));
    page.push_str(&format!(
        "<style>\n{STYLE}</style>\n</head>\n<body>\n<main>\n"
    ));
    page.push_str(&formula(input, &[]));
    page.push_str(&readings_list(input, readings, options));
    page.push_str(&format!(
        "</main>\n<script>\n{SCRIPT}</script>\n</body>\n</html>\n"
    ));
    page
}

/// Renders the body of the viewer for any `interpretation` of `input`.
///
/// Readings are shown as in [`render_html`]. An underdetermined formula says
/// why it encodes no structure, and an invalid one lists its diagnostics; both
/// underline the characters at fault. Every diagnostic label carries the
/// `data-span` of the characters it points at, so the viewer traces it too.
///
/// The fragment has no page around it: the embedding page provides [`STYLE`]
/// and [`SCRIPT`].
#[must_use]
pub fn render_fragment(
    input: &str,
    interpretation: &Interpretation,
    options: &SvgOptions,
) -> String {
    match interpretation {
        Interpretation::Underdetermined { span, reason } => format!(
            "{}<p class=\"summary\">{}</p>\n",
            formula(input, &[*span]),
            underdetermination(input, *span, *reason)
        ),
        Interpretation::Invalid { diagnostics } => {
            let flagged: Vec<Span> = diagnostics.iter().map(|d| d.primary.span).collect();
            format!(
                "{}{}{}",
                formula(input, &flagged),
                summary(0),
                problems(input, diagnostics)
            )
        }
        _ => formula(input, &[]) + &readings_list(input, interpretation.readings(), options),
    }
}

/// Renders the body of the viewer for a molecule given as `smiles`: its drawing
/// and canonical SMILES, or why the string cannot be read. A drawing from SMILES
/// has no characters to trace.
///
/// The fragment has no page around it, as for [`render_fragment`].
#[must_use]
pub fn render_smiles_fragment(smiles: &str, options: &SvgOptions) -> String {
    let heading = format!(
        "<h1 class=\"formula smiles\" aria-label=\"SMILES {0}\">{0}</h1>\n",
        escape(smiles)
    );
    let molecule = match moiety_core::smiles::parse(smiles) {
        Ok(molecule) => molecule,
        Err(error) => {
            return format!(
                "{heading}<p class=\"summary\">This SMILES string cannot be read: {}.</p>\n",
                escape(&error.to_string())
            );
        }
    };
    let canonical = molecule.to_smiles();
    let drawing = render_svg(
        &molecule,
        &SvgOptions {
            title: Some(canonical.clone()),
            ..options.clone()
        },
    );
    format!(
        concat!(
            "{heading}<p class=\"summary\">Drawn from SMILES.</p>\n",
            "<ol class=\"readings\">\n<li class=\"reading\">\n",
            "<h2>Structure</h2>\n",
            "<p class=\"smiles\"><code>{smiles}</code></p>\n",
            "<figure>{drawing}</figure>\n",
            "</li>\n</ol>\n"
        ),
        heading = heading,
        smiles = escape(&canonical),
        drawing = drawing,
    )
}

/// The formula, one span per character, each with its byte range in `input`.
/// Characters inside a `flagged` span get the class `flag`.
fn formula(input: &str, flagged: &[Span]) -> String {
    let characters: String = input
        .char_indices()
        .map(|(offset, c)| {
            let end = offset + c.len_utf8();
            let flag = if flagged.iter().any(|s| s.start < end && offset < s.end) {
                " class=\"flag\""
            } else {
                ""
            };
            format!(
                "<span data-offset=\"{offset}\" data-end=\"{end}\"{flag}>{}</span>",
                escape(c.encode_utf8(&mut [0; 4]))
            )
        })
        .collect();
    format!(
        "<h1 class=\"formula\" aria-label=\"Formula {}\">{characters}</h1>\n",
        escape(input)
    )
}

fn summary(count: usize) -> String {
    let text = match count {
        0 => "No structure can be read from this formula.".to_string(),
        1 => "1 reading. Point at, focus or select characters, atoms or bonds to trace them."
            .to_string(),
        n => format!(
            "{n} readings, best first. Point at, focus or select characters, atoms or bonds to trace them."
        ),
    };
    format!("<p class=\"summary\">{text}</p>\n")
}

fn readings_list(input: &str, readings: &[Reading], options: &SvgOptions) -> String {
    let mut html = summary(readings.len());
    if !readings.is_empty() {
        html.push_str("<ol class=\"readings\">\n");
        for reading in readings {
            html.push_str(&card(input, reading, options));
        }
        html.push_str("</ol>\n");
    }
    html
}

fn quoted(input: &str, span: Span) -> String {
    let text = input.get(span.start..span.end).unwrap_or_default();
    format!("<code>{}</code>", escape(text))
}

fn underdetermination(input: &str, span: Span, reason: Underdetermination) -> String {
    let text = quoted(input, span);
    match reason {
        Underdetermination::MolecularFormula => format!(
            "{text} is a molecular formula, which gives the atoms but not how they are connected."
        ),
        Underdetermination::UnknownRing => {
            format!("{text} is a ring the dictionary cannot expand.")
        }
        Underdetermination::UnknownAbbreviation => {
            format!("{text} is not an element or an abbreviation in the dictionary.")
        }
    }
}

/// Every diagnostic with its code, message, labelled spans and notes.
fn problems(input: &str, diagnostics: &[Diagnostic]) -> String {
    let mut html = String::from("<ul class=\"problems\">\n");
    for diagnostic in diagnostics {
        html.push_str(&format!(
            "<li>\n<p><code>{}</code> {}</p>\n<ul class=\"labels\">\n",
            diagnostic.code,
            escape(&diagnostic.message)
        ));
        for label in std::iter::once(&diagnostic.primary).chain(&diagnostic.secondary) {
            html.push_str(&format!(
                "<li data-span=\"{}-{}\">{} {}</li>\n",
                label.span.start,
                label.span.end,
                quoted(input, label.span),
                escape(&label.message)
            ));
        }
        html.push_str("</ul>\n");
        for note in &diagnostic.notes {
            html.push_str(&format!("<p class=\"note\">Note: {}</p>\n", escape(note)));
        }
        html.push_str("</li>\n");
    }
    html.push_str("</ul>\n");
    html
}

/// Score, direction and the ambiguities a reading resolved, with the choice it made at each.
fn notes(input: &str, reading: &Reading) -> String {
    let mut notes = format!(
        "<p class=\"meta\">Score {}, written {}</p>\n",
        reading.score,
        reading.direction.as_str()
    );
    if !reading.ambiguities.is_empty() {
        notes.push_str("<ul class=\"ambiguities\">\n");
        for ambiguity in &reading.ambiguities {
            let text = input
                .get(ambiguity.span.start..ambiguity.span.end)
                .unwrap_or_default();
            notes.push_str(&format!(
                "<li>{} at <code>{}</code>: {}</li>\n",
                ambiguity.class,
                escape(text),
                escape(&ambiguity.choice)
            ));
        }
        notes.push_str("</ul>\n");
    }
    notes
}

fn card(input: &str, reading: &Reading, options: &SvgOptions) -> String {
    let title = format!("Reading {}: {}", reading.rank, reading.smiles);
    let drawing = render_reading_svg(
        reading,
        &SvgOptions {
            title: Some(title),
            ..options.clone()
        },
    );
    format!(
        concat!(
            "<li class=\"reading\" data-rank=\"{rank}\">\n",
            "<h2>Reading {rank}</h2>\n",
            "<p class=\"smiles\"><code>{smiles}</code></p>\n",
            "{notes}",
            "<figure>{drawing}</figure>\n",
            "</li>\n"
        ),
        rank = reading.rank,
        smiles = escape(&reading.smiles),
        notes = notes(input, reading),
        drawing = drawing,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn page(input: &str) -> String {
        render_html(
            input,
            moiety::interpret(input).readings(),
            &SvgOptions::default(),
        )
    }

    #[test]
    fn renders_smiles_or_says_why_not() {
        let html = render_smiles_fragment("c1ccccc1O<", &SvgOptions::default());
        assert!(html.contains("cannot be read"));
        assert!(html.contains("c1ccccc1O&lt;"));
        let html = render_smiles_fragment("Oc1ccccc1", &SvgOptions::default());
        assert!(html.contains("Drawn from SMILES."));
        assert!(html.contains("<svg"));
        assert_eq!(html.matches("<li class=\"reading\"").count(), 1);
    }

    #[test]
    fn has_one_span_per_input_character() {
        for input in ["CH3CH2OH", "CuSO4·5H2O", "HC≡CCH3", "CH₃COO⁻"] {
            let html = page(input);
            assert_eq!(
                html.matches("data-offset=").count(),
                input.chars().count(),
                "{input}"
            );
        }
        let html = page("CuSO4·5H2O");
        assert!(html.contains("<span data-offset=\"5\" data-end=\"7\">·</span>"));
        assert!(html.contains("<span data-offset=\"7\" data-end=\"8\">5</span>"));
    }

    #[test]
    fn shows_every_reading_with_rank_and_smiles() {
        let html = page("COCl2");
        assert!(html.starts_with("<!DOCTYPE html>"));
        assert_eq!(html.matches("<li class=\"reading\"").count(), 2);
        assert_eq!(html.matches("<svg ").count(), 2);
        assert!(html.contains("<h2>Reading 2</h2>"));
        assert!(html.contains("<code>Cl[Co]Cl</code>"));
        assert!(html.contains("<li>symbol-segmentation at <code>CO</code>: Co</li>"));
        assert!(html.contains("Score -6, written left-to-right"));
        assert!(html.contains("data-span="));
        assert!(html.contains("prefers-color-scheme"));
    }

    fn fragment(input: &str) -> String {
        render_fragment(input, &moiety::interpret(input), &SvgOptions::default())
    }

    #[test]
    fn fragment_draws_every_reading_without_the_page() {
        let html = fragment("COCl2");
        assert!(html.starts_with("<h1 class=\"formula\""));
        assert!(html.contains("2 readings, best first."));
        assert_eq!(html.matches("<li class=\"reading\"").count(), 2);
        assert!(!html.contains("class=\"flag\""));
        for page in ["<!DOCTYPE", "<style", "<script", "<main"] {
            assert!(!html.contains(page), "{page}");
        }
    }

    #[test]
    fn fragment_explains_underdetermined_formulas() {
        let html = fragment("CH3Xq");
        assert!(html.contains("<span data-offset=\"2\" data-end=\"3\">3</span>"));
        assert!(html.contains("<span data-offset=\"3\" data-end=\"4\" class=\"flag\">X</span>"));
        assert!(
            html.contains(
                "<code>Xq</code> is not an element or an abbreviation in the dictionary."
            )
        );
        assert!(!html.contains("<svg"));
        assert!(fragment("C6H12O6").contains(
            "<code>C6H12O6</code> is a molecular formula, which gives the atoms but not how they are connected."
        ));
        assert!(
            fragment("PhC6H4")
                .contains("<code>C6H4</code> is a ring the dictionary cannot expand.")
        );
    }

    #[test]
    fn fragment_lists_diagnostics_with_traceable_labels() {
        let html = fragment("CH3CH3OH");
        assert!(html.contains("No structure can be read from this formula."));
        assert!(html.contains("<code>CF0204</code> chain cannot continue past a terminal group"));
        assert!(html.contains(
            "<li data-span=\"3-6\"><code>CH3</code> this group has no free valence left</li>"
        ));
        assert!(
            html.contains("<li data-span=\"6-8\"><code>OH</code> but more units follow it</li>")
        );
        assert!(html.contains("<p class=\"note\">Note: did you mean CH3CH2OH?</p>"));
        assert_eq!(html.matches("class=\"flag\"").count(), 3);
    }

    #[test]
    fn fragment_escapes_the_input() {
        let html = fragment("<b>");
        assert!(html.contains("<code>CF0101</code>"));
        assert!(html.contains("<code>&lt;</code>"));
        assert!(!html.contains("<b>"));
    }

    #[test]
    fn is_self_contained_and_escaped() {
        let html = render_html("<b>&\"", &[], &SvgOptions::default());
        assert!(html.contains("&lt;b&gt;&amp;&quot;"));
        assert!(!html.contains("<b>"));
        assert!(html.contains("No structure can be read"));
        for external in ["<link", " src=", "@import", "url("] {
            assert!(!page("PhCH2OH").contains(external), "{external}");
        }
    }
}
