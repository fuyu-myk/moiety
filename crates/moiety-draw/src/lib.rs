//! moiety-draw turns the molecular graphs that moiety reads into skeletal
//! (line-angle) structure drawings.
//!
//! Three layers build on each other:
//!
//! - [`layout`](mod@layout) places every atom in the plane, one unit per bond, with
//!   zig-zag chains, regular-polygon rings and components side by side.
//! - [`svg`] draws a layout as a self-contained SVG in the usual skeletal
//!   conventions. Drawn from a [`Reading`](moiety_core::Reading), every atom
//!   and bond carries the byte span of the input it came from.
//! - [`html`] puts the formula and the drawings of all its readings on one
//!   page, where pointing at a character highlights the atoms and bonds it
//!   produced, and pointing at an atom or bond highlights its characters.
//!
//! ```
//! let interpretation = moiety::interpret("CH3CH2OH");
//! let reading = interpretation.top().expect("ethanol has a reading");
//!
//! let svg = moiety_draw::render_reading_svg(reading, &moiety_draw::SvgOptions::default());
//! assert!(svg.starts_with("<svg"));
//! assert!(svg.contains(r#"data-atom="2" data-span="6-8""#));
//!
//! let page = moiety_draw::render_html("CH3CH2OH", interpretation.readings(), &Default::default());
//! assert_eq!(page.matches("data-offset=").count(), 8);
//! ```

mod graph;
pub mod html;
pub mod layout;
mod rings;
pub mod svg;
mod xml;

pub use html::{render_fragment, render_html};
pub use layout::{Layout, Point, layout};
pub use svg::{SvgOptions, render_reading_svg, render_svg};
