# moiety

moiety reads condensed and semi-structural organic formulas such as `CH3CH2OH`, `(CH3)2CHCOOH` or `EtO2CCH2CO2Et`, turns them into molecular graphs and canonical SMILES, and draws them as skeletal structures.

It never guesses. Every reading that survives valence resolution is returned, ranked by a declared policy, and an ambiguous formula comes back in a different shape from a resolved one, so a caller cannot mistake one for the other.

```text
$ moiety parse '(CH3)2CHCOOH' COCl2 CH3CH3OH
(CH3)2CHCOOH
  unique: CC(C)C(=O)O  (score 1, mixed)
COCl2
  ambiguous: 2 readings (symbol-segmentation)
  1. O=C(Cl)Cl  (score 0, left-to-right)
  2. Cl[Co]Cl  (score -6, left-to-right)
     1: symbol-segmentation at "CO": C O
     2: symbol-segmentation at "CO": Co
CH3CH3OH
  error[CF0204]: chain cannot continue past a terminal group
    --> CH3CH3OH
           ^^^ this group has no free valence left
              ^^ but more units follow it
    note: did you mean CH3CH2OH?
```

The language is specified in [docs/grammar.md](docs/grammar.md): normalization, the token lattice, the surface grammar, valence resolution (rules R1 to R5), the ambiguity taxonomy, ranking policy 2 and the diagnostic codes.

## Crates

| Crate | Contents |
| --- | --- |
| `moiety-core` | Normalization, lexer, parser, valence resolution, interpretation, SMILES reading and canonical writing. No dependencies. |
| `moiety` | The bundled abbreviation dictionary and its TOML loader, the inverse writer, and the conformance runner. |
| `moiety-draw` | Template layout for the v1 input language, SVG rendering, and an HTML page with bidirectional provenance highlighting. |
| `moiety-cli` | The `moiety` command. |
| `moiety-wasm` | The WebAssembly module behind the browser playground. Not published. |

Canonical SMILES lives in `moiety-core` rather than `moiety`, because deduplicating readings depends on it.

## Using the library

```rust
use moiety::Interpretation;

match moiety::interpret("CH3COO-") {
    Interpretation::Unique(reading) => println!("{}", reading.smiles),
    Interpretation::Ambiguous { readings, classes } => {
        // Acetate, CC([O-])=O, ranks first; *OC(C)=O is the acetoxy group with an open valence.
        for reading in &readings {
            println!("{} {}", reading.rank, reading.smiles);
        }
        println!("{classes:?}");
    }
    Interpretation::Underdetermined { span, reason } => println!("{reason} at {span}"),
    Interpretation::Invalid { diagnostics } => {
        for diagnostic in &diagnostics {
            eprint!("{}", diagnostic.render("CH3COO-"));
        }
    }
}
```

Every `Reading` carries its canonical SMILES, rank, score, the direction it was written in, the ambiguities it resolved with the choice it made at each, and provenance: the byte span each atom and bond came from, mapped back through Unicode folds so `CH₃` still points at the right characters.

`moiety::interpret_with` takes a custom dictionary and `Options` (charge convention, reading limit, search limit, suggestions). Dictionaries are TOML files; see [crates/moiety/data/dictionary.toml](crates/moiety/data/dictionary.toml) for the format and `moiety::dictionary::extend` to add entries to the default set.

## Command line

```text
moiety parse <formula>...     every reading of each formula
moiety write <smiles>...      condensed formulas from SMILES (the inverse writer)
moiety draw <formula>         skeletal structure as SVG
moiety html <formula>         a page with every reading and provenance highlighting
moiety conformance            the conformance suite and round-trip corpora, with the v1.0 gate
moiety differential <file>    compare with another implementation's results on substituent labels
moiety dictionary             the abbreviation dictionary
```

`--dictionary <file>` adds a TOML dictionary to the default one, `--hypervalent` writes nitro groups as `N(=O)=O`, and `-o <file>` writes output to a file.

## Drawing

`moiety-draw` lays out the closed set of structures the v1 language can express (chains, branches, and the rings that come from the dictionary) with templates rather than a general layout engine: zig-zag chains, regular polygons for rings, straight lines through triple bonds and cumulated double bonds, and components side by side. Drawings follow skeletal conventions and use `currentColor`, so they take the colour of the page they sit in.

```text
moiety draw 'BocNHCH2COOH' -o boc-glycine.svg
moiety html 'CH3SOCH3' -o dmso.html
```

The HTML page puts the formula above a card for every reading, with its rank, score, direction and the choice it made at each ambiguity. Pointing at, focusing or selecting characters of the formula highlights the atoms and bonds they produced in every reading, and pointing at an atom or bond highlights its characters.

## Playground

The [playground](https://fuyu-myk.github.io/moiety/) runs moiety in the browser: every reading is drawn as you type, with the same two-way highlighting as the HTML page. Underdetermined and invalid formulas say why they cannot be read, with the characters at fault underlined.

It is a static site in [playground](playground), deployed to GitHub Pages from `main` by [.github/workflows/pages.yml](.github/workflows/pages.yml). `moiety-wasm` compiles the parser and the renderer to a WebAssembly module that exports three functions and needs no generated bindings, and a Web Worker runs it, so a slow formula never blocks typing. To build and serve it locally:

```text
rustup target add wasm32-unknown-unknown
playground/build.sh target/playground
python3 -m http.server --directory target/playground
```

## Conformance

The suite is data, not Rust tests, so any implementation can run it: [conformance/cases.toml](conformance/cases.toml) lists each input with its expected result kind, readings in rank order, ambiguity classes, and for invalid input the diagnostic code and primary span. [conformance/roundtrip.smi](conformance/roundtrip.smi) is a corpus of well-known molecules that the inverse writer turns into condensed formulas and the parser reads back.

Current results (`moiety conformance`):

| Metric | Result | v1.0 target |
| --- | --- | --- |
| Cases passed | 294 / 294 | |
| False-unique rate | 0 / 27 | 0% |
| Ambiguity class recall | 27 / 27 | at least 98% |
| Top-1 accuracy, ambiguous cases | 27 / 27 | at least 85% |
| Diagnostic span accuracy | 23 / 23 | at least 90% |
| Round-trip recovery | 685 / 685 (185 curated, 500 generated) | at least 95% |
| Ambiguity rate | 27 / 294 | reported |

Differential testing against MolScribe's condensed-formula expander covers 246 substituent labels: 182 match exactly, every disagreement has a recorded verdict, and none is decided against moiety. See [conformance/differential](conformance/differential/README.md).

The ambiguous cases are adversarial by design: symbol collisions such as `COCl2` (phosgene or cobalt chloride), shadowed abbreviations such as `TsCl` (tosyl or tennessine), trailing dashes that could be a charge or an open valence, and units that can be a pendant or a chain link, such as the oxygen in `CH3SOCH3`.

## Development

```text
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo run --release -p moiety-cli -- conformance
```

CI runs the same commands, including the conformance gate and a playground build, from [.github/workflows/ci.yml](.github/workflows/ci.yml).

The fuzzing harness in [fuzz](fuzz) reads and draws arbitrary input under memory and time limits, as the playground does. It needs a nightly toolchain and `cargo-fuzz`:

```text
python3 -I fuzz/seed.py
cd fuzz
cargo +nightly fuzz run -s none playground corpus/playground -- -dict=formula.dict -max_len=300 -timeout=3 -rss_limit_mb=512 -malloc_limit_mb=256
```

Changes are listed in [CHANGELOG.md](CHANGELOG.md).

## Status

The parser meets every item of the v1.0 gate: the metrics above hold, the conformance suite is published and runnable by a third party, CI runs it on every push, the inverse writer is complete, and the crates are published on crates.io.

```text
cargo install moiety-cli
```

## License

MIT. See [LICENSE](LICENSE).
