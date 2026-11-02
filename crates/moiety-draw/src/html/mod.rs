//! A self-contained HTML page that traces every reading of a formula back to its characters.
//!
//! The page shows the formula with one `<span data-offset>` per character and,
//! below it, the drawing of every reading with its rank and SMILES. Pointing at
//! or focusing an atom or bond highlights the characters it came from; pointing
//! at, focusing or selecting characters highlights every atom and bond whose
//! span overlaps them, in all readings at once. Styles and script are inline,
//! nothing is loaded from elsewhere, and colours follow `prefers-color-scheme`.

use moiety_core::Reading;

use crate::svg::{SvgOptions, render_reading_svg};
use crate::xml::escape;

const STYLE: &str = include_str!("viewer.css");
const SCRIPT: &str = include_str!("viewer.js");

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
    page.push_str(&formula(input));
    page.push_str(&summary(readings.len()));
    if !readings.is_empty() {
        page.push_str("<ol class=\"readings\">\n");
        for reading in readings {
            page.push_str(&card(input, reading, options));
        }
        page.push_str("</ol>\n");
    }
    page.push_str(&format!(
        "</main>\n<script>\n{SCRIPT}</script>\n</body>\n</html>\n"
    ));
    page
}

/// The formula, one span per character, each with its byte range in `input`.
fn formula(input: &str) -> String {
    let characters: String = input
        .char_indices()
        .map(|(offset, c)| {
            format!(
                "<span data-offset=\"{offset}\" data-end=\"{}\">{}</span>",
                offset + c.len_utf8(),
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
