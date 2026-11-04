"""Write the label glyph outlines that moiety-draw paints pixel drawings with.

Usage:

    python -I glyphs.py <Inter-Regular.ttf> <out.rs>

Reads Inter Regular 3.019 (SIL Open Font License 1.1, see Inter-LICENSE.txt)
and writes the advance width and outline of every character an atom label can
use: the Latin letters, the digits, `+`, `*` and the minus sign. Outlines are
in font units with `y` growing upward, each contour a move followed by lines
and quadratic curves.

Needs fontTools.
"""

import sys

from fontTools.ttLib import TTFont

CHARACTERS = (
    [chr(c) for c in range(ord("0"), ord("9") + 1)]
    + [chr(c) for c in range(ord("A"), ord("Z") + 1)]
    + [chr(c) for c in range(ord("a"), ord("z") + 1)]
    + ["*", "+", "−"]
)


def contours(glyph, glyf):
    coordinates, ends, flags = glyph.getCoordinates(glyf)
    start = 0
    for end in ends:
        points = [
            (int(x), int(y), bool(flags[i] & 1))
            for i, (x, y) in enumerate(coordinates[start : end + 1], start)
        ]
        start = end + 1
        yield points


def steps(points):
    """Moves, lines and quadratic curves for one closed TrueType contour."""
    if not any(on for _, _, on in points):
        x, y = (points[0][0] + points[1][0]) // 2, (points[0][1] + points[1][1]) // 2
        points = [(x, y, True)] + points
    first = next(i for i, p in enumerate(points) if p[2])
    points = points[first:] + points[:first]
    x0, y0, _ = points[0]
    out = [f"M({x0}, {y0})"]
    control = None
    for x, y, on in points[1:] + [points[0]]:
        if on:
            out.append(f"Q({control[0]}, {control[1]}, {x}, {y})" if control else f"L({x}, {y})")
            control = None
        elif control:
            mx, my = (control[0] + x) // 2, (control[1] + y) // 2
            out.append(f"Q({control[0]}, {control[1]}, {mx}, {my})")
            control = (x, y)
        else:
            control = (x, y)
    return out


def main(font_path, out_path):
    font = TTFont(font_path)
    glyf = font["glyf"]
    cmap = font.getBestCmap()
    lines = [
        "//! Outlines of the characters atom labels use, from Inter Regular 3.019.",
        "//!",
        "//! Copyright 2020 The Inter Project Authors (https://github.com/rsms/inter),",
        "//! licensed under the SIL Open Font License 1.1: see `fonts/Inter-LICENSE.txt`.",
        "//! Written by `fonts/glyphs.py`.",
        "",
        "use super::Step::{self, L, M, Q};",
        "",
        "/// Font units per em.",
        f"pub(super) const UNITS_PER_EM: f64 = {font['head'].unitsPerEm}.0;",
        "",
        "/// Each character with its advance width and outline, in font units, sorted by character.",
        "#[rustfmt::skip]",
        "pub(super) static GLYPHS: &[(char, u16, &[Step])] = &[",
    ]
    for character in sorted(CHARACTERS):
        name = cmap[ord(character)]
        advance = font["hmtx"][name][0]
        outline = [step for contour in contours(glyf[name], glyf) for step in steps(contour)]
        escaped = "\\u{2212}" if character == "−" else character
        lines.append(f"    ('{escaped}', {advance}, &[{', '.join(outline)}]),")
    lines.append("];")
    with open(out_path, "w") as out:
        out.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:])
