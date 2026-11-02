//! Reads and draws a formula exactly as the playground does.

#![no_main]

use libfuzzer_sys::fuzz_target;
use moiety_draw::SvgOptions;

fuzz_target!(|data: &[u8]| {
    if let Ok(input) = std::str::from_utf8(data) {
        let interpretation = moiety::interpret(input);
        let html = moiety_draw::render_fragment(input, &interpretation, &SvgOptions::default());
        std::hint::black_box(html);
    }
});
