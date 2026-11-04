# Changelog

Notable changes to the moiety crates. Versions follow [Semantic Versioning](https://semver.org/); before 1.0, a minor version can change how formulas are read. The grammar, the ranking policy, the bundled dictionary and the conformance suite carry their own versions, given where they change.

## Unreleased

Grammar version 4.

### Added

- Template rings take bracketed locants and locants after their substituents, as benzene does. They must fill every substituted position, with nothing before the ring: `C5H3N(2-Cl)(6-Cl)`, `C9H6N(8-OH)`, `C5H3NCl2-2,6`, `C5H3NCl-2-CH3-6` (grammar section 6.12).
- Drawing from SMILES: `moiety draw --smiles`, `moiety_draw::render_smiles_fragment`, the wasm export `moiety_render_smiles`, and a Formula/SMILES switch in the playground with SMILES examples.
- The SMILES reader keeps stereochemistry: tetrahedral centres written `@`, `@@`, `@TH1` or `@TH2`, and double bond configurations from bonds written `/` and `\`, as `Molecule::tetrahedral` and `Molecule::double_bond_stereo` with the new `Tetrahedral` and `DoubleBondStereo`. Canonical SMILES still leaves stereo out.
- Drawings from SMILES show each tetrahedral stereocentre with one wedge or hash, and draw every double bond with the configuration the SMILES gave, in chains and in macrocycles.
- Macrocycles of 12 to 46 atoms are drawn as closed zig-zags on a hexagonal lattice, as chemists draw them, with fused rings and the largest substituents outside and an odd ring's one irregular bond away from substituents: erythromycin, tacrolimus and rapamycin look like their usual drawings.
- Porphyrins, the morphine skeleton and norbornanes take fixed drawings from PubChem's depictions.
- The `smiles` fuzz target.

### Changed

- A bracket that starts with a locant after a template ring now holds bracketed substituents, as it does after benzene, and a phenyl spelling in it is part of the substituent. So `C5H4N(2-OPh)` reads as 2-phenoxypyridine, and `C5H4N(4-ClPh)` is an error, as `C6H4(4-ClPh)` already was, instead of an underdetermined ring.
- Units after a template ring whose brackets fill it, as in `C5H4N(2-Cl)CH3`, get CF0108 at the ring instead of CF0204.
- A locant that lists a position whose atom cannot take a substituent suggests the nearest position with the same element that can: `3-MeC3H3N2` suggests `1-MeC3H3N2`, the nitrogen that holds a hydrogen.
- Bridged ring systems are drawn without crossing bonds where the convention has none: a bridge across a ring is an arc or a shallow bow, chosen with the substituents in mind, and a system that still draws badly is started from another ring or relaxed, whichever is most legible. Bicyclo[2.2.2]octanes, quinuclidines, tropanes, artemisinin and cryptands no longer cross.
- `layout` no longer keeps every bond one unit long: bridges may be shorter and templates follow their source, within 0.55 to 1.6.

## 0.5.0 - 2026-11-03

Grammar version 3, dictionary 1.4.0, conformance suite 0.5.0.

### Added

- Ring templates, which number rings other than benzene from the dictionary (grammar section 6.12):
  - substituted rings as whole molecules: `4-Me2NC5H4N`, `2,6-Me2C5H3N`, `2-Me3SiCH2C5H4N`;
  - attachment locants before the ring: `2-C5H4NCH2OH`, `-5-Br-2-C5H3N`, `HOCH2-2,6-C5H3N-CH2OH`;
  - attachment locants after it: `(C5H4N-2)3P`, `Ph2P(C5H4N-2)`, `-C10H7-1`;
  - substituents after the ring: `C5H4N-2-CH=NNHPh`.
- `[[ring]]` tables in dictionary files, with `RingTemplate`, `Dictionary::rings`, `Dictionary::insert_ring`, `Dictionary::remove_ring`, `MAX_POSITIONS`, and the errors `DictionaryError::BadRing` and `DictionaryError::DuplicateRing`. A file that extends a dictionary replaces the ring template with the same name or the same formula, so it can bring pyrazole in place of the bundled imidazole.
- Dictionary 1.4.0: templates for pyridine, pyrimidine, pyrrole, imidazole, thiazole, thiophene, furan, naphthalene, quinoline and indole, numbered as IUPAC numbers them: `8-HOC9H6N`, `3-C8H6NCH2COOH`, `2-H2NC3H2NS`, `1-MeC3H3N2`.
- When two ring formulas match at one place, the longer is read, so a user template such as thienothiophene's `C6H2S2` wins over benzene's `C6H2`.
- The writer writes template rings with locants: `4-C5H4NN(CH3)CH3`, `2,6-Br2C5H3N`, `CH3-2,6-C5H3NCH3`, `-(5-Br-2-C5H3N)`, `3-C4H3SCl`. A ring with no hydrogen left fails with the new `WriteError::FullRing`, and a template's parent ring on its own, such as pyridine, with `WriteError::BareRing`.

### Changed

- The public error enums `DictionaryError`, `WriteError`, `LoadError`, `SuiteError` and `DifferentialError` are `#[non_exhaustive]`, so a `match` on them needs a wildcard arm.
- Error messages no longer repeat their source error, which the CLI prints after them, and a ring template's SMILES error names the ring.
- A group that opens with a template ring followed by its substituents, as in `P(C5H4N-2-CH3)3`, reports CF0108 at the ring, since nothing is left for the host, instead of CF0206.
- The generated round-trip corpus grows one molecule in eight each from pyridine and from thiophene.
- A locant directly before a benzene ring may follow its substituents, listing the positions they leave: `2-Cl-1,4-C6H3(OH)2` is 2-chlorohydroquinone.
- Locants may list positions up to 8, which naphthalene needs.
- Drawings start from a hub, an atom outside the rings with three or more alike branches, as in `PPh3`, `Ph3PO` or `(MeO)3PO`, and space its branches evenly with one pointing up; four branches make an upright cross. When the branches crowd each other and no flip or turn parts them without bending the hub, its bonds stretch equally, up to 1.5 bond lengths, so `B(C6F5)3` keeps its rings at 120 degrees.

### Fixed

- The SMILES reader gave aromatic sulfur, and aromatic nitrogen with three bonds, an implicit hydrogen they do not have. Canonical SMILES text is unchanged.
- The SMILES reader made a bond between two aromatic atoms aromatic even when it lies on no ring, so the writer could not write biaryls whose SMILES joins the rings without `-`, such as `c1ccc(cc1)c1ccccn1`. Such a bond is now single.
- A ring that opens a group starting with a bond now bonds to the group's host by its lowest listed position, as the grammar says. `CH3CH(-1,2,4-C6H3(Cl)Br)CH3` put the chlorine at 1 and the host at 4 before.

## 0.4.0 - 2026-11-02

Grammar version 2, dictionary 1.3.0, conformance suite 0.4.0.

### Added

- Locants that number benzene rings: a list before the ring as in `4-MeOC6H4CHO` and `2,4-Cl2C6H3CO2H`, one per substituent as in `2-Cl-4-FC6H3CHO`, `o-`, `m-` and `p-` as in `p-CH3C6H4OH`, locants on the ring itself as in `o-C6H4(CO2H)2`, locants after the ring as in `C6H4-4-OMe` and `P(C6H4-4-OMe)3`, locants after their substituents as in `C6H4OMe-4`, `C6H3Me2-2,6` and `C6H3Cl-2-F-4`, bracketed locants as in `C6H4(4-OMe)` and `C6H3(2-Cl)(4-F)`, and a bare locant in a chain as in `HO-p-C6H4-CO2H` and `-p-C6H4-`. `ortho-`, `meta-` and `para-` may be spelled out.
  A counted unit after a numbered ring puts one copy on each open position. A ring with no locant is still underdetermined.
- Locants before a phenyl spelling number a substituted phenyl: `4-ClPh`, `3,5-(CF3)2PhBr`, `(4-FPh)3P`.
- CF0108 for locants that do not fit their ring.
- `lexer::Locant`, `EdgeKind::Locant` and `EdgeKind::TrailingLocant`.
- Dictionary spellings may start with a locant, as `2-C5H4N` does, and are matched only where a locant may start.
- Dictionary 1.3.0:
  - heteroaryl and naphthyl: 2-, 3- and 4-pyridyl (`2-C5H4N`, `2-Py`), 2- and 3-thienyl (`2-C4H3S`), 2- and 3-furyl (`2-C4H3O`), 1- and 2-naphthyl (`1-C10H7`, `1-Naph`);
  - aryl and alkyl: pentafluorophenyl (`C6F5`), `PMP`, `Dipp`, 1- and 2-adamantyl (`Ad`, `2-Ad`);
  - protecting groups: `Fmoc`, `Alloc`, `Troc`, `Teoc`, `TBDPS`, `SEM`, `MEM`;
  - boryl and imido groups: `Bpin`, `Bcat`, `NPhth`.

### Changed

- The writer writes benzene rings the dictionary does not cover with locants: `CH3CH2C6H4-4-CHO`, `1,4-Cl2C6H4`, `-C6H3-3,5-(CF3)2`, `-1,4-C6H4-`, and `C6H6` for benzene. Fused and bridged rings and rings with no hydrogen left stay unwritable.
- The writer reads each formula back before returning it. When leaving out a bond order would read as another molecule, it writes every multiple bond out, as `CH2CH=NH`, and when no formula reads back it returns the new `WriteError::Unreadable` instead of a wrong formula.
- The generated round-trip corpus grows one molecule in four from a benzene ring.
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
