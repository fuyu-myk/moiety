//! Reads and draws a SMILES string exactly as the playground does.

#![no_main]

use libfuzzer_sys::fuzz_target;
use moiety_draw::SvgOptions;

fuzz_target!(|data: &[u8]| {
    if let Ok(input) = std::str::from_utf8(data) {
        let html = moiety_draw::render_smiles_fragment(input, &SvgOptions::default());
        std::hint::black_box(html);
    }
});
