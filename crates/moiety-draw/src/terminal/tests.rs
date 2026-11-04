use moiety_core::smiles;

use super::*;

fn env<'a>(pairs: &'a [(&'a str, &'a str)]) -> impl Fn(&str) -> Option<String> + 'a {
    move |key| {
        pairs
            .iter()
            .find(|(k, _)| *k == key)
            .map(|(_, v)| v.to_string())
    }
}

fn show(input: &str, protocol: Protocol, columns: usize) -> String {
    let options = TerminalOptions {
        protocol,
        columns,
        ..TerminalOptions::default()
    };
    render_terminal(&smiles::parse(input).unwrap(), &options)
}

#[test]
fn encodes_base64_as_the_rfc_does() {
    let cases = [
        ("", ""),
        ("f", "Zg=="),
        ("fo", "Zm8="),
        ("foo", "Zm9v"),
        ("foob", "Zm9vYg=="),
        ("fooba", "Zm9vYmE="),
        ("foobar", "Zm9vYmFy"),
    ];
    for (plain, encoded) in cases {
        assert_eq!(base64(plain.as_bytes()), encoded);
    }
    assert_eq!(base64(&[0xfb, 0xff]), "+/8=");
}

#[test]
fn detects_the_protocol_from_the_environment() {
    let cases: [(&[(&str, &str)], Protocol); 11] = [
        (&[("KITTY_WINDOW_ID", "1")], Protocol::Kitty),
        (&[("TERM", "xterm-kitty")], Protocol::Kitty),
        (&[("TERM", "xterm-ghostty")], Protocol::Kitty),
        (&[("TERM_PROGRAM", "iTerm.app")], Protocol::Iterm2),
        (&[("TERM_PROGRAM", "WezTerm")], Protocol::Iterm2),
        (&[("TERM", "foot")], Protocol::Sixel),
        (&[("WT_SESSION", "x")], Protocol::Sixel),
        (&[("KONSOLE_VERSION", "230804")], Protocol::Sixel),
        (&[("TERM", "xterm-256color")], Protocol::Unicode),
        (
            &[("TERM", "xterm-kitty"), ("TMUX", "/tmp/t,1,0")],
            Protocol::Unicode,
        ),
        (
            &[("TERM", "screen-256color"), ("KITTY_WINDOW_ID", "1")],
            Protocol::Unicode,
        ),
    ];
    for (pairs, protocol) in cases {
        assert_eq!(detect(env(pairs)), protocol, "{pairs:?}");
    }
    assert_eq!(detect(env(&[])), Protocol::Unicode);
    let inherited: [&[(&str, &str)]; 4] = [
        &[("TERM_PROGRAM", "vscode"), ("KITTY_WINDOW_ID", "1")],
        &[
            ("TERM_PROGRAM", "Apple_Terminal"),
            ("ITERM_SESSION_ID", "w0"),
        ],
        &[("TERM_PROGRAM", "vscode"), ("WT_SESSION", "x")],
        &[("TERM", "xterm-kitty"), ("ZELLIJ", "0")],
    ];
    for pairs in inherited {
        assert_eq!(detect(env(pairs)), Protocol::Unicode, "{pairs:?}");
    }
    assert_eq!(
        detect(env(&[
            ("TERM_PROGRAM", "WezTerm"),
            ("KITTY_WINDOW_ID", "1")
        ])),
        Protocol::Iterm2
    );
}

#[test]
fn takes_the_theme_from_the_background_colour() {
    assert_eq!(detect_theme(env(&[("COLORFGBG", "0;15")])), Theme::Light);
    assert_eq!(
        detect_theme(env(&[("COLORFGBG", "15;default;0")])),
        Theme::Dark
    );
    assert_eq!(detect_theme(env(&[("COLORFGBG", "0;7")])), Theme::Light);
    assert_eq!(detect_theme(env(&[("COLORFGBG", "junk")])), Theme::Dark);
    assert_eq!(detect_theme(env(&[])), Theme::Dark);
}

/// The `(control, payload)` pairs of a stream of Kitty escapes.
fn kitty_chunks(out: &str) -> Vec<(&str, &str)> {
    out.trim_end()
        .split("\x1b\\")
        .filter(|chunk| !chunk.is_empty())
        .map(|chunk| {
            let body = chunk.strip_prefix("\x1b_G").expect("APC start");
            body.split_once(';').expect("control and payload")
        })
        .collect()
}

#[test]
fn sends_kitty_images_in_bounded_chunks() {
    let out = show("CC(=O)Oc1ccccc1C(=O)O", Protocol::Kitty, 80);
    assert!(out.ends_with("\x1b\\\n"));
    let chunks = kitty_chunks(&out);
    assert!(chunks.len() > 1, "aspirin needs several chunks");
    let (first, _) = chunks[0];
    assert!(first.starts_with("a=T,f=100,q=2,c="), "{first}");
    assert!(first.ends_with(",m=1"));
    assert!(
        chunks[1..chunks.len() - 1]
            .iter()
            .all(|(c, _)| *c == "q=2,m=1")
    );
    assert_eq!(chunks.last().unwrap().0, "q=2,m=0");
    assert!(
        chunks
            .iter()
            .all(|(_, payload)| payload.len() <= KITTY_CHUNK)
    );
    let payload: String = chunks.iter().map(|(_, p)| *p).collect();
    assert!(payload.starts_with(&base64(b"\x89PNG\r\n\x1a\n")[..8]));
    let columns: usize = first["a=T,f=100,q=2,c=".len()..first.len() - 4]
        .parse()
        .unwrap();
    assert!(
        (20..=35).contains(&columns),
        "aspirin spans {columns} columns"
    );
    let narrow = show("CC(=O)Oc1ccccc1C(=O)O", Protocol::Kitty, 6);
    assert!(narrow.starts_with("\x1b_Ga=T,f=100,q=2,c=6,"));
}

#[test]
fn shows_iterm2_images_inline() {
    let out = show("CCO", Protocol::Iterm2, 80);
    let rest = out
        .strip_prefix("\x1b]1337;File=inline=1;size=")
        .expect("OSC 1337 header");
    let (size, rest) = rest.split_once(';').unwrap();
    let (fields, rest) = rest.split_once(':').unwrap();
    assert!(fields.starts_with("width=") && fields.ends_with(";preserveAspectRatio=1"));
    let payload = rest.strip_suffix("\x07\n").expect("BEL and newline");
    assert_eq!(
        payload.len(),
        size.parse::<usize>().unwrap().div_ceil(3) * 4
    );
}

/// The number at the start of `rest`, which moves past it.
fn number(rest: &mut &str) -> usize {
    let digits = rest
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(rest.len());
    let (n, tail) = rest.split_at(digits);
    *rest = tail;
    n.parse().unwrap()
}

/// The shade of every pixel of a Sixel image, and its size.
fn decode_sixel(out: &str) -> (usize, usize, Vec<u8>) {
    let body = out
        .strip_prefix("\x1bP0;1;0q\"1;1;")
        .and_then(|b| b.strip_suffix("\x1b\\\n"))
        .expect("DCS header and ST");
    let (size, mut rest) = body.split_at(body.find('#').unwrap());
    let (width, height) = size.split_once(';').unwrap();
    let (width, height): (usize, usize) = (width.parse().unwrap(), height.parse().unwrap());
    let mut shades = vec![u8::MAX; width * height];
    let (mut x, mut top, mut shade) = (0, 0, 0u8);
    while let Some(c) = rest.chars().next() {
        rest = &rest[c.len_utf8()..];
        match c {
            '#' => {
                shade = number(&mut rest) as u8;
                if rest.starts_with(";2;") {
                    rest = &rest[rest.find('#').unwrap_or(rest.len())..];
                }
            }
            '$' => x = 0,
            '-' => (x, top) = (0, top + 6),
            '!' | '?'..='~' => {
                let count = if c == '!' { number(&mut rest) } else { 1 };
                let letter = if c == '!' {
                    let l = rest.chars().next().unwrap();
                    rest = &rest[1..];
                    l
                } else {
                    c
                };
                let bits = letter as u8 - 63;
                for _ in 0..count {
                    for dy in (0..6).filter(|dy| bits & 1 << dy != 0) {
                        shades[(top + dy) * width + x] = shade;
                    }
                    x += 1;
                }
            }
            _ => panic!("unexpected {c:?}"),
        }
    }
    (width, height, shades)
}

#[test]
fn paints_every_sixel_pixel_its_shade() {
    let molecule = smiles::parse("C[N+](C)(C)CC([O-])=O").unwrap();
    let out = show("C[N+](C)(C)CC([O-])=O", Protocol::Sixel, 80);
    assert!(out.contains("#15;2;90;90;90") && out.contains("#0;2;9;9;9"));
    let (width, height, shades) = decode_sixel(&out);
    let (ink, card) = Theme::Dark.colours();
    let raster = RasterOptions {
        ink,
        background: card,
        ..RasterOptions::default()
    };
    let mask = Mask::new(
        &molecule,
        &raster,
        SIXEL_CELL / CELL_WIDTH,
        80.0 * SIXEL_NARROW_CELL,
    );
    assert_eq!((width, height), (mask.width, mask.height));
    let expected: Vec<u8> = mask
        .alpha
        .iter()
        .map(|&a| (a * 15.0).round() as u8)
        .collect();
    assert_eq!(shades, expected);
    assert!(out.contains('!'), "long runs are repeated");
}

#[test]
fn writes_labels_as_characters_among_braille_bonds() {
    let ethanol = show("CCO", Protocol::Unicode, 80);
    assert!(
        ethanol.contains("OH") || ethanol.contains("HO"),
        "{ethanol}"
    );
    assert!(
        ethanol
            .chars()
            .any(|c| ('\u{2801}'..='\u{28ff}').contains(&c))
    );
    assert!(show("O", Protocol::Unicode, 80).contains("H₂O"));
    assert!(show("[NH4+]", Protocol::Unicode, 80).contains("NH₄⁺"));
    let acetate = show("CC(=O)[O-]", Protocol::Unicode, 80);
    assert!(acetate.contains("O⁻"), "{acetate}");
    let benzene = show("c1ccccc1", Protocol::Unicode, 80);
    assert!(!benzene.chars().any(char::is_alphanumeric), "{benzene}");
    assert!(benzene.ends_with('\n') && !benzene.lines().any(|l| l.ends_with(' ')));
}

#[test]
fn shrinks_unicode_drawings_to_the_columns_given() {
    let chain = "C".repeat(60);
    let wide = show(&chain, Protocol::Unicode, 1000);
    assert!(wide.lines().map(|l| l.chars().count()).max().unwrap() > 200);
    let narrow = show(&chain, Protocol::Unicode, 40);
    assert!(narrow.lines().all(|l| l.chars().count() <= 40), "{narrow}");
    let label = show(
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCO",
        Protocol::Unicode,
        30,
    );
    assert!(label.lines().all(|l| l.chars().count() <= 30));
    assert!(label.contains('O'));
}

#[test]
fn keeps_whole_labels_at_the_edges() {
    let chain = "C".repeat(30);
    let ions = show(&format!("[O-]{chain}[NH3+]"), Protocol::Unicode, 80);
    assert!(ions.contains("NH₃⁺") && ions.contains("O⁻"), "{ions}");
    let amine = show(
        &format!("N{}C(=O)[O-]", "C".repeat(40)),
        Protocol::Unicode,
        80,
    );
    assert!(amine.contains("H₂N"), "{amine}");
    assert!(amine.lines().all(|l| l.chars().count() <= 80));
    let narrow = show("OCCN", Protocol::Unicode, 12);
    assert!(narrow.contains("HO") || narrow.contains("OH"), "{narrow}");
    assert!(narrow.contains("NH₂") || narrow.contains("H₂N"), "{narrow}");
    assert!(narrow.lines().all(|l| l.chars().count() <= 12), "{narrow}");
}
