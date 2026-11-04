//! Skeletal structures shown in a terminal: an inline image over the Kitty,
//! iTerm2 or Sixel graphics protocol, or Unicode text that any terminal shows.
//!
//! Images are painted as [`render_image`](crate::render_image) paints them, on
//! a card in the colours of a [`Theme`], and sized in character cells so that
//! a bond is about five columns long. Unicode drawings put bonds in braille
//! dots and labels in plain characters (`OH`, `NH₂`, `N⁺`) in the terminal's
//! own colours.
//!
//! ```
//! use moiety_draw::terminal::{Protocol, TerminalOptions, render_terminal};
//!
//! let molecule = moiety::smiles::parse("CCO").expect("ethanol");
//! let options = TerminalOptions {
//!     protocol: Protocol::Unicode,
//!     ..TerminalOptions::default()
//! };
//! let text = render_terminal(&molecule, &options);
//! assert!(text.contains("OH"));
//! ```

mod braille;
mod sixel;

use moiety_core::Molecule;

use crate::raster::{Mask, RasterOptions};
use crate::svg::{SvgOptions, picture};

/// Width of a character cell, in pixels of a drawing at the default sizes.
const CELL_WIDTH: f64 = 8.0;

/// Height of a character cell, in pixels of a drawing at the default sizes.
const CELL_HEIGHT: f64 = 16.0;

/// How much larger Unicode drawings are than images, so that the strokes of a
/// double bond fall on separate rows of dots.
const TEXT_ZOOM: f64 = 1.5;

/// How many image pixels a terminal that scales images gets per drawing pixel.
const IMAGE_ZOOM: f64 = 2.0;

/// Device pixels assumed per character cell across, for Sixel images, which
/// terminals do not scale.
const SIXEL_CELL: f64 = 10.0;

/// The narrowest cell assumed when keeping a Sixel image within the columns
/// it may span, so that it is not cut off where cells are small.
const SIXEL_NARROW_CELL: f64 = 7.0;

/// Longest payload of one Kitty escape, fixed by the protocol.
const KITTY_CHUNK: usize = 4096;

/// How a terminal can show a drawing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum Protocol {
    /// Kitty's graphics protocol, also spoken by Ghostty.
    Kitty,
    /// iTerm2's inline images, also shown by WezTerm.
    Iterm2,
    /// Sixel graphics, shown by foot, Konsole, Windows Terminal and others.
    Sixel,
    /// Braille dots and plain characters, for any terminal.
    Unicode,
}

/// The colours of the card an image is drawn on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum Theme {
    /// Light ink on a dark card.
    Dark,
    /// Dark ink on a light card.
    Light,
}

impl Theme {
    /// The ink and card colours, as red, green, blue and alpha.
    fn colours(self) -> ([u8; 4], [u8; 4]) {
        match self {
            Theme::Dark => ([230, 230, 230, 255], [24, 24, 24, 255]),
            Theme::Light => ([0, 0, 0, 255], [255, 255, 255, 255]),
        }
    }
}

/// How to show a drawing in a terminal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TerminalOptions {
    /// The protocol the terminal speaks.
    pub protocol: Protocol,
    /// The most columns the drawing may span. Wider drawings are shrunk to
    /// fit, labels included. Sixel images, which terminals do not scale, are
    /// fitted for cells 7 pixels wide.
    pub columns: usize,
    /// The colours of the card an image is drawn on. Unicode drawings take
    /// the terminal's own colours.
    pub theme: Theme,
}

impl Default for TerminalOptions {
    /// Unicode text in 80 columns, on a dark card.
    fn default() -> Self {
        Self {
            protocol: Protocol::Unicode,
            columns: 80,
            theme: Theme::Dark,
        }
    }
}

/// The protocol a terminal speaks, judged from its environment variables, as
/// `getenv` looks them up. Inside tmux, screen or Zellij, which pass no
/// images on, and in terminals not recognised, it is [`Protocol::Unicode`].
///
/// `TERM_PROGRAM` and `TERM` come first. Variables that only name a session,
/// such as `KITTY_WINDOW_ID`, count only when `TERM_PROGRAM` is unset, since
/// programs started from a terminal inherit them: VS Code opened from kitty
/// still has `KITTY_WINDOW_ID`.
pub fn detect(getenv: impl Fn(&str) -> Option<String>) -> Protocol {
    let value = |key: &str| getenv(key).unwrap_or_default();
    let (term, program) = (value("TERM"), value("TERM_PROGRAM"));
    let multiplexed = getenv("TMUX").is_some() || getenv("ZELLIJ").is_some();
    if multiplexed || term.starts_with("screen") || term.starts_with("tmux") {
        return Protocol::Unicode;
    }
    match program.as_str() {
        "ghostty" => return Protocol::Kitty,
        "iTerm.app" | "WezTerm" => return Protocol::Iterm2,
        _ => {}
    }
    if term.contains("kitty") || term.contains("ghostty") {
        return Protocol::Kitty;
    }
    let sixel_terms = ["foot", "mlterm", "contour", "yaft"];
    if sixel_terms.iter().any(|name| term.starts_with(name)) {
        return Protocol::Sixel;
    }
    if !program.is_empty() {
        return Protocol::Unicode;
    }
    if getenv("KITTY_WINDOW_ID").is_some() {
        Protocol::Kitty
    } else if getenv("ITERM_SESSION_ID").is_some() {
        Protocol::Iterm2
    } else if getenv("KONSOLE_VERSION").is_some() || getenv("WT_SESSION").is_some() {
        Protocol::Sixel
    } else {
        Protocol::Unicode
    }
}

/// The theme that matches the terminal's background, judged from
/// `COLORFGBG` as `getenv` looks it up, or [`Theme::Dark`] when it is unset.
pub fn detect_theme(getenv: impl Fn(&str) -> Option<String>) -> Theme {
    let background = getenv("COLORFGBG").and_then(|value| {
        value
            .rsplit(';')
            .next()
            .and_then(|bg| bg.parse::<u8>().ok())
    });
    match background {
        Some(7 | 9..=15) => Theme::Light,
        _ => Theme::Dark,
    }
}

/// Draws `molecule` for a terminal, ending with a newline so that whatever is
/// printed next starts below the drawing.
#[must_use]
pub fn render_terminal(molecule: &Molecule, options: &TerminalOptions) -> String {
    let columns = options.columns.max(1);
    let (ink, card) = options.theme.colours();
    let raster = RasterOptions {
        ink,
        background: card,
        ..RasterOptions::default()
    };
    match options.protocol {
        Protocol::Unicode => braille::render(&picture(molecule, &SvgOptions::default()), columns),
        Protocol::Sixel => {
            let zoom = SIXEL_CELL / CELL_WIDTH;
            let mask = Mask::new(molecule, &raster, zoom, columns as f64 * SIXEL_NARROW_CELL);
            sixel::encode(&mask, ink, card) + "\n"
        }
        Protocol::Kitty | Protocol::Iterm2 => {
            let cell = CELL_WIDTH * IMAGE_ZOOM;
            let mask = Mask::new(molecule, &raster, IMAGE_ZOOM, columns as f64 * cell);
            let span = ((mask.width as f64 / cell).ceil() as usize).clamp(1, columns);
            let png = mask.to_image(ink, card).to_png();
            if options.protocol == Protocol::Kitty {
                kitty(&png, span)
            } else {
                iterm2(&png, span)
            }
        }
    }
}

/// A Kitty escape that sends `png` and shows it `columns` cells wide.
fn kitty(png: &[u8], columns: usize) -> String {
    let payload = base64(png);
    let count = payload.len().div_ceil(KITTY_CHUNK);
    let mut out = String::with_capacity(payload.len() + count * 24 + 64);
    for (i, start) in (0..payload.len()).step_by(KITTY_CHUNK).enumerate() {
        let chunk = &payload[start..(start + KITTY_CHUNK).min(payload.len())];
        let more = u8::from(i + 1 < count);
        let keys = if i == 0 {
            format!("a=T,f=100,q=2,c={columns},")
        } else {
            "q=2,".to_string()
        };
        out.push_str(&format!("\x1b_G{keys}m={more};{chunk}\x1b\\"));
    }
    out.push('\n');
    out
}

/// An iTerm2 escape that shows `png` inline, `columns` cells wide.
fn iterm2(png: &[u8], columns: usize) -> String {
    format!(
        "\x1b]1337;File=inline=1;size={};width={columns};preserveAspectRatio=1:{}\x07\n",
        png.len(),
        base64(png)
    )
}

/// `bytes` in standard base64 with padding.
fn base64(bytes: &[u8]) -> String {
    const ALPHABET: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for group in bytes.chunks(3) {
        let word = group
            .iter()
            .enumerate()
            .fold(0u32, |word, (i, &b)| word | u32::from(b) << (16 - 8 * i));
        for k in 0..4 {
            if k <= group.len() {
                out.push(char::from(ALPHABET[(word >> (18 - 6 * k) & 63) as usize]));
            } else {
                out.push('=');
            }
        }
    }
    out
}

#[cfg(test)]
mod tests;
