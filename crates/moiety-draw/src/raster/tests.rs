use moiety_core::smiles;

use super::*;
use crate::render_svg;

fn image(input: &str) -> Image {
    render_image(&smiles::parse(input).unwrap(), &RasterOptions::default())
}

fn alpha(image: &Image, x: f64, y: f64) -> u8 {
    let (x, y) = (x.floor() as usize, y.floor() as usize);
    image.pixels()[(y * image.width() as usize + x) * 4 + 3]
}

/// The `viewBox` of an SVG drawing, as origin and size.
fn view_box(svg: &str) -> (Point, Point) {
    let start = svg.find("viewBox=\"").unwrap() + 9;
    let values: Vec<f64> = svg[start..svg[start..].find('"').unwrap() + start]
        .split(' ')
        .map(|v| v.parse().unwrap())
        .collect();
    (
        Point::new(values[0], values[1]),
        Point::new(values[2], values[3]),
    )
}

#[test]
fn paints_the_drawing_the_svg_shows() {
    let molecule = smiles::parse("CCO").unwrap();
    let image = render_image(&molecule, &RasterOptions::default());
    let (origin, size) = view_box(&render_svg(&molecule, &SvgOptions::default()));
    assert_eq!(image.width(), size.x.ceil() as u32);
    assert_eq!(image.height(), size.y.ceil() as u32);
    let picture = picture(&molecule, &SvgOptions::default());
    let Shape::Line(start, end) = picture.shapes[0] else {
        panic!("the first bond is a plain line");
    };
    let middle = (start + end) * 0.5 - origin;
    let darkest = (-1..=1)
        .flat_map(|dx| (-1..=1).map(move |dy| (f64::from(dx), f64::from(dy))))
        .map(|(dx, dy)| alpha(&image, middle.x + dx, middle.y + dy))
        .max();
    assert!(darkest > Some(200), "bond is inked: {darkest:?}");
    let Some(Shape::Label { centre, .. }) = picture.shapes.last() else {
        panic!("the oxygen is labelled");
    };
    let label = *centre - origin;
    let inked = (-6..=6)
        .flat_map(|dx| (-6..=6).map(move |dy| (dx, dy)))
        .filter(|&(dx, dy)| alpha(&image, label.x + f64::from(dx), label.y + f64::from(dy)) > 128)
        .count();
    assert!(inked > 20, "the O is drawn: {inked} pixels");
    assert_eq!(alpha(&image, 0.0, 0.0), 0, "padding stays clear");
}

#[test]
fn fills_wedges_and_rings() {
    let solid = image("C[C@H](O)F");
    let plain = image("CC(O)F");
    let ink = |image: &Image| {
        image
            .pixels()
            .chunks(4)
            .map(|p| u32::from(p[3]))
            .sum::<u32>()
    };
    assert!(
        ink(&solid) > ink(&plain),
        "a wedge holds more ink than a line"
    );
    let circled = image("c1cccc1");
    assert!(ink(&circled) > ink(&image("C1CCCC1")));
}

#[test]
fn blends_ink_over_the_background() {
    let white = [255, 255, 255, 255];
    let black = [0, 0, 0, 255];
    assert_eq!(blend(white, black, 0.0), white);
    assert_eq!(blend(white, black, 1.0), black);
    assert_eq!(blend(white, black, 0.5), [128, 128, 128, 255]);
    assert_eq!(blend([0; 4], [10, 20, 30, 255], 0.5), [10, 20, 30, 128]);
    assert_eq!(blend([0; 4], [10, 20, 30, 0], 1.0), [0; 4]);
    let image = render_image(
        &smiles::parse("C").unwrap(),
        &RasterOptions {
            background: white,
            ..RasterOptions::default()
        },
    );
    assert_eq!(&image.pixels()[..4], &white);
}

#[test]
fn scales_huge_drawings_down() {
    let chain = "C".repeat(4000);
    let image = image(&chain);
    assert!(f64::from(image.width()) * f64::from(image.height()) <= MAX_PIXELS);
    assert!(image.width() > 1000);
    let odd = RasterOptions {
        bond_length: f64::NAN,
        font_size: -1.0,
        ..RasterOptions::default()
    };
    let molecule = smiles::parse("CCO").unwrap();
    assert_eq!(
        render_image(&molecule, &odd),
        render_image(&molecule, &RasterOptions::default())
    );
    for bond_length in [1e200, f64::MAX] {
        let huge = RasterOptions {
            bond_length,
            ..RasterOptions::default()
        };
        let image = render_image(&molecule, &huge);
        assert!(
            f64::from(image.width()) * f64::from(image.height()) <= MAX_PIXELS,
            "{bond_length} {}x{}",
            image.width(),
            image.height()
        );
    }
}
