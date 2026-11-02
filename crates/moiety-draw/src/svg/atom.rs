//! Atom groups: the label or attachment mark of each atom, a title, and an
//! invisible halo that a viewer can light up.

use moiety_core::AtomKind;

use super::{Scene, num};
use crate::layout::{Bounds, Point};
use crate::xml::escape;

/// Radius of the halo around an unlabeled atom, in bond lengths.
const HALO: f64 = 0.3;

/// Writes the group of atom `index`.
pub(super) fn render(
    scene: &Scene,
    index: usize,
    span: &str,
    out: &mut String,
    bounds: &mut Bounds,
) {
    let centre = scene.points[index];
    let font = scene.metrics.font;
    let label = scene.labels[index].as_ref();
    let radius = label.map_or(HALO * scene.metrics.bond, |l| l.radius(font));
    out.push_str(&format!(
        "<g class=\"atom\" data-atom=\"{index}\"{span}><title>{}</title>\
         <circle class=\"halo\" cx=\"{}\" cy=\"{}\" r=\"{}\" fill-opacity=\"0\"/>",
        escape(&describe(scene, index)),
        num(centre.x),
        num(centre.y),
        num(radius)
    ));
    match label {
        Some(label) => label.render(centre, font, out, bounds),
        None if scene.molecule.atoms()[index].kind == AtomKind::Attachment => {
            wavy_mark(scene, index, out, bounds);
        }
        None => {}
    }
    out.push_str("</g>");
}

fn describe(scene: &Scene, index: usize) -> String {
    let atom = scene.molecule.atoms()[index];
    match atom.element() {
        None => format!("attachment point {index}"),
        Some(element) if atom.charge != 0 => {
            format!(
                "atom {index}: {}, charge {:+}",
                element.symbol(),
                atom.charge
            )
        }
        Some(element) => format!("atom {index}: {}", element.symbol()),
    }
}

/// A short wavy line across the end of the bond to attachment point `index`.
fn wavy_mark(scene: &Scene, index: usize, out: &mut String, bounds: &mut Bounds) {
    const WAVES: usize = 4;
    let centre = scene.points[index];
    let along = scene.adjacency[index]
        .first()
        .and_then(|&(next, _)| (centre - scene.points[next]).normalized())
        .unwrap_or(Point::new(1.0, 0.0));
    let across = along.perpendicular();
    let length = 0.5 * scene.metrics.bond;
    let amplitude = 0.1 * scene.metrics.bond;
    let start = centre - across * (length / 2.0);
    let mut path = format!("M{} {}", num(start.x), num(start.y));
    for wave in 0..WAVES {
        let bulge = if wave % 2 == 0 { amplitude } else { -amplitude };
        let control =
            start + across * (length * (wave as f64 + 0.5) / WAVES as f64) + along * bulge;
        let end = start + across * (length * (wave + 1) as f64 / WAVES as f64);
        path.push_str(&format!(
            " Q{} {} {} {}",
            num(control.x),
            num(control.y),
            num(end.x),
            num(end.y)
        ));
    }
    bounds.include_box(centre, length / 2.0 + amplitude, length / 2.0 + amplitude);
    out.push_str(&format!(
        "<path class=\"attachment\" d=\"{path}\" fill=\"none\" stroke=\"currentColor\" \
         stroke-width=\"{}\" stroke-linecap=\"round\"/>",
        num(scene.metrics.line)
    ));
}
