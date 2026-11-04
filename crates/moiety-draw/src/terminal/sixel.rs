//! Sixel images: bands six pixels high, each pixel one of a few shades
//! between the background and the ink.

use crate::raster::{Mask, blend};

/// Shades from the background, shade 0, to the ink.
const SHADES: usize = 16;

/// Runs of one sixel at least this long are written as a repeat.
const MIN_REPEAT: usize = 4;

/// `mask` as a Sixel image in `ink` on an opaque `background`.
pub(super) fn encode(mask: &Mask, ink: [u8; 4], background: [u8; 4]) -> String {
    let shades: Vec<u8> = mask
        .alpha
        .iter()
        .map(|&alpha| (alpha * (SHADES - 1) as f32).round() as u8)
        .collect();
    let mut out = format!("\x1bP0;1;0q\"1;1;{};{}", mask.width, mask.height);
    for shade in 0..SHADES {
        let [r, g, b, _] = blend(background, ink, shade as f32 / (SHADES - 1) as f32);
        let percent = |c: u8| (u32::from(c) * 100 + 127) / 255;
        out.push_str(&format!(
            "#{shade};2;{};{};{}",
            percent(r),
            percent(g),
            percent(b)
        ));
    }
    for top in (0..mask.height).step_by(6) {
        let rows = (mask.height - top).min(6);
        let mut first = true;
        for shade in 0..SHADES as u8 {
            let sixels: Vec<u8> = (0..mask.width)
                .map(|x| {
                    (0..rows)
                        .filter(|&dy| shades[(top + dy) * mask.width + x] == shade)
                        .fold(0, |bits, dy| bits | 1 << dy)
                })
                .collect();
            if sixels.iter().all(|&bits| bits == 0) {
                continue;
            }
            if !first {
                out.push('$');
            }
            first = false;
            out.push_str(&format!("#{shade}"));
            write_runs(&sixels, &mut out);
        }
        out.push('-');
    }
    out.push_str("\x1b\\");
    out
}

/// Writes one colour's sixels across a band, repeating long runs.
fn write_runs(sixels: &[u8], out: &mut String) {
    let end = sixels
        .iter()
        .rposition(|&bits| bits != 0)
        .map_or(0, |i| i + 1);
    let mut at = 0;
    while at < end {
        let bits = sixels[at];
        let run = sixels[at..end].iter().take_while(|&&b| b == bits).count();
        let letter = char::from(63 + bits);
        if run >= MIN_REPEAT {
            out.push_str(&format!("!{run}{letter}"));
        } else {
            (0..run).for_each(|_| out.push(letter));
        }
        at += run;
    }
}
