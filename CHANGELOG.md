# Changelog

Notable changes to the moiety crates. Versions follow [Semantic Versioning](https://semver.org/); before 1.0, a minor version can change how formulas are read. The grammar, the ranking policy, the bundled dictionary and the conformance suite carry their own versions, given where they change.

## Unreleased

Grammar version 2, dictionary 1.3.0.

### Added

- Locants that number benzene rings: a list before the ring as in `4-MeOC6H4CHO` and `2,4-Cl2C6H3CO2H`, one per substituent as in `2-Cl-4-FC6H3CHO`, `o-`, `m-` and `p-` as in `p-CH3C6H4OH`, locants on the ring itself as in `o-C6H4(CO2H)2`, locants after the ring as in `C6H4-4-OMe` and `P(C6H4-4-OMe)3`, and locants after their substituent as in `C6H4OMe-4` and `C6H3Me2-2,6`. A counted unit after a numbered ring puts one copy on each open position. A ring with no locant is still underdetermined.
- CF0108 for locants that do not fit their ring.
- `lexer::Locant`, `EdgeKind::Locant` and `EdgeKind::TrailingLocant`.
- Dictionary spellings may start with a locant, as `2-C5H4N` does, and are matched only where a locant may start.
- Dictionary 1.3.0: 2-, 3- and 4-pyridyl (`2-C5H4N`, `2-Py`), 2- and 3-thienyl (`2-C4H3S`), 2- and 3-furyl (`2-C4H3O`), 1- and 2-naphthyl (`1-C10H7`, `1-Naph`), and `Fmoc`.

### Changed

- The writer puts a `-` before a spelling that starts with a locant unless a locant could start there, as in `Ph-2-C5H4N`, and writes out repeated components that start with one instead of counting them.

### Fixed

- The writer counted repeated components that start with an attachment marker, writing `*C.*C` as `2-CH3`, which does not read back. They are now written out: `-CH3.-CH3`.

## 0.3.0 - 2026-11-02

Dictionary 1.2.0, conformance suite 0.3.0.

### Added

- A [browser playground](https://fuyu-myk.github.io/moiety/) that draws every reading as you type, with two-way highlighting, and shows why a formula cannot be read.
- `moiety_draw::render_fragment`, which renders the viewer body for any interpretation. Underdetermined and invalid formulas underline the characters at fault, and diagnostic labels are traced like atoms and bonds. The viewer's `STYLE` and `SCRIPT` are now public.
- `moiety-wasm`, an unpublished WebAssembly module behind the playground.
- `C6H6` reads as benzene, a conventional reading like `C6H5` for phenyl.
- A fuzzing harness in `fuzz/` that reads and draws arbitrary input under memory and time limits.
- CF0004 for a reading with more than 1,000 atoms, counting every copy that a count makes.
- `moiety --version`.

### Changed

- The viewer script traces content that a page adds after it loads.

### Fixed

- Nested group counts such as `((((CH2)2)2)2)` no longer exhaust memory. Every atom a count copies is charged to the search budget, and deep nesting stops with CF0004.
- Fuzzing found formulas of under 20 characters that took up to 6 seconds or 1 GB to read. Every limit below fails with a diagnostic or marks the readings with CF0303; none silently changes a result:
  - The search charges each state as it is made, gives states that cannot complete an allowance of their own, and caps every frontier at 400,000 atoms. Duplicate states are found by hashing instead of copying their keys.
  - Canonical SMILES reuses the result for identical molecules and components, skips symmetric terminal atoms when breaking ties, refines without allocating per atom, and may visit 25 atoms per search step for a whole formula before CF0003.
  - Readings are limited to 1,000 atoms. A formula with any larger reading reports CF0004 instead of returning only its smaller readings.
  - The suggestion search checks whether an edit reads without ranking it, then confirms the one it offers with a full read.
  - The largest real formula needs a sixteenth of the frontier cap and a five-hundredth of the canonical SMILES allowance, and every conformance case reads as before.
- Drawing a reading of several hundred atoms took seconds, nearly all of it spent moving branches apart. Clash removal now stops after comparing 4,000,000 atom pairs in a molecule, which no real formula's drawing comes close to.

## 0.2.0 - 2026-11-02

Ranking policy 2, dictionary 1.1.0, conformance suite 0.2.0.

### Added

- Differential testing against MolScribe's condensed-formula expander on 246 substituent labels, with a recorded verdict for every disagreement, and the `moiety differential` command.
- Pentafluoroethyl (`C2F5`) in the dictionary.

### Changed

- Ranking policy 2: a written charge on carbon costs 2, so `CH3CO-` ranks the acetyl group above the acyl anion.

### Fixed

- A hydrogen written after a run of halogens bonds to their host, so `CF2H` reads as difluoromethyl.
- A written charge can land on an atom inside a repeated group, so `CH3(CH2)3-` is ambiguous instead of wrongly unique.

## 0.1.0 - 2026-11-02

First release: grammar version 1, ranking policy 1, dictionary 1.0.0, conformance suite 0.1.0.

### Added

- `moiety-core`: the condensed-formula lexer, parser, valence resolution and interpretation, with canonical SMILES and no dependencies.
- `moiety`: the bundled abbreviation dictionary, the inverse writer and the conformance runner.
- `moiety-draw`: skeletal structure layout, SVG rendering and the provenance viewer.
- `moiety-cli`: the `moiety` command.
