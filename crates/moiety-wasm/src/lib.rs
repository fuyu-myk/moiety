//! The WebAssembly module behind the moiety playground.
//!
//! The module imports nothing and exports its memory and four functions, so a
//! page can call it without generated bindings:
//!
//! 1. `moiety_input(len)` makes room for a UTF-8 input of `len` bytes and
//!    returns its address in memory, where the caller writes the input.
//! 2. `moiety_render()` reads the input as a formula and returns the length in
//!    bytes of the HTML it rendered, [`render`] of the formula.
//! 3. `moiety_render_smiles()` does the same for a SMILES string, with
//!    [`render_smiles`].
//! 4. `moiety_output()` returns the address of that HTML.
//!
//! Memory can grow during a render, so views of it must be made after the call.

use std::cell::RefCell;

use moiety_draw::SvgOptions;

thread_local! {
    static INPUT: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
    static OUTPUT: RefCell<String> = const { RefCell::new(String::new()) };
}

/// The playground's HTML for `input`: [`moiety_draw::html::render_fragment`] of
/// its interpretation with the bundled dictionary and default options.
#[must_use]
pub fn render(input: &str) -> String {
    moiety_draw::render_fragment(input, &moiety::interpret(input), &SvgOptions::default())
}

/// The playground's HTML for `input` read as SMILES:
/// [`moiety_draw::html::render_smiles_fragment`] with default options.
#[must_use]
pub fn render_smiles(input: &str) -> String {
    moiety_draw::render_smiles_fragment(input, &SvgOptions::default())
}

/// Makes room for an input of `len` bytes and returns where to write it.
#[expect(unsafe_code, reason = "the exported names are unique to this module")]
#[unsafe(no_mangle)]
pub extern "C" fn moiety_input(len: usize) -> *mut u8 {
    INPUT.with_borrow_mut(|input| {
        input.clear();
        input.resize(len, 0);
        input.as_mut_ptr()
    })
}

/// Renders the input as a formula and returns the length of the output in bytes.
#[expect(unsafe_code, reason = "the exported names are unique to this module")]
#[unsafe(no_mangle)]
pub extern "C" fn moiety_render() -> usize {
    store(INPUT.with_borrow(|input| render(&String::from_utf8_lossy(input))))
}

/// Renders the input as SMILES and returns the length of the output in bytes.
#[expect(unsafe_code, reason = "the exported names are unique to this module")]
#[unsafe(no_mangle)]
pub extern "C" fn moiety_render_smiles() -> usize {
    store(INPUT.with_borrow(|input| render_smiles(&String::from_utf8_lossy(input))))
}

fn store(html: String) -> usize {
    OUTPUT.with_borrow_mut(|output| {
        *output = html;
        output.len()
    })
}

/// Where the output of the last [`moiety_render`] starts.
#[expect(unsafe_code, reason = "the exported names are unique to this module")]
#[unsafe(no_mangle)]
pub extern "C" fn moiety_output() -> *const u8 {
    OUTPUT.with_borrow(|output| output.as_ptr())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Calls the exports as a page does, writing the input through the buffer
    /// that `moiety_input` hands out.
    fn call(input: &[u8]) -> String {
        let at = moiety_input(input.len());
        INPUT.with_borrow_mut(|buffer| {
            assert_eq!(buffer.as_mut_ptr(), at);
            buffer.copy_from_slice(input);
        });
        let len = moiety_render();
        OUTPUT.with_borrow(|output| {
            assert_eq!(output.as_ptr(), moiety_output());
            assert_eq!(output.len(), len);
            output.clone()
        })
    }

    #[test]
    fn renders_through_the_exports() {
        let html = call(b"EtO2CCH2CO2Et");
        assert_eq!(html, render("EtO2CCH2CO2Et"));
        assert!(html.contains("1 reading."));
        let html = call(b"COCl2");
        assert_eq!(html.matches("<li class=\"reading\"").count(), 2);
        assert!(call(b"CH3CH3OH").contains("<code>CF0204</code>"));
    }

    #[test]
    fn renders_smiles_through_the_exports() {
        let at = moiety_input(9);
        INPUT.with_borrow_mut(|buffer| {
            assert_eq!(buffer.as_mut_ptr(), at);
            buffer.copy_from_slice(b"Oc1ccccc1");
        });
        let len = moiety_render_smiles();
        let html = OUTPUT.with_borrow(|output| output.clone());
        assert_eq!(html.len(), len);
        assert_eq!(html, render_smiles("Oc1ccccc1"));
        assert!(html.contains("Drawn from SMILES."));
    }

    #[test]
    fn reads_invalid_utf8_as_replacement_characters() {
        let html = call(b"CH3\xffOH");
        assert!(html.contains(">\u{fffd}</span>"));
        assert!(html.contains("<code>CF0101</code>"));
    }
}
