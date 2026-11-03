# moiety grammar, version 3

This document is the normative specification of the condensed-formula language that moiety reads. The reference implementation lives in `crates/moiety-core`, and the conformance suite in `conformance/` is the executable form of the rules below. Where this document and the implementation disagree, the implementation has a bug.

Grammar version: **3**. Ranking policy version: **2**.

Version 3 adds ring templates, which number heteroaromatic rings such as pyridine (section 6.12), and lets a locant directly before a benzene ring follow its substituents (section 6.11). Version 2 added locants that number benzene rings. Version 1 read every `C6H4` between other units as an unknown ring.

## 1. Scope

moiety turns a condensed or semi-structural formula into a molecular graph and reports every reading that survives valence resolution. It never returns a single best guess for an ambiguous input.

In scope:

- Condensed and semi-structural formulas: `CH3CH2OH`, `(CH3)2CHCOOH`, `CH3(CH2)4CH3`, `HOCH2CH2OH`.
- Abbreviations from a declared, versioned dictionary: `Ph`, `Me`, `Et`, `Boc`, `Ts`, `OMe`, `CO2Et`.
- Benzene rings numbered by locants: `4-MeOC6H4CHO`, `2,4-Cl2C6H3CO2H`, `o-C6H4(CO2H)2`.
- Heteroaromatic rings numbered by locants through ring templates: `4-Me2NC5H4N`, `2,6-Me2C5H3N`, `(C5H4N-2)3P`.
- Formal charges where written, attachment markers, and dot-separated multi-component formulas.
- Simple salts whose cation is written first: `CuSO4`, `Na2CO3` (rule R5).

Out of scope, reported rather than guessed:

- Stereochemistry and isotopes. Descriptors are skipped with a warning; `D` and `T` are read as hydrogen with a warning.
- Systematic names.
- Molecular formulas such as `C6H12O6`, which carry no connectivity.
- Rings that the dictionary does not define, apart from benzene and template rings numbered by locants.
- Polymer repeat units such as `[CH2CH2]n`.

## 2. Normalization

Normalization runs before lexing and records every substitution in a fold table, so any span in the normalized string maps back to the original input.

| Input | Becomes |
| --- | --- |
| Subscript digits `₀` to `₉` | ASCII digits |
| A run of superscript digits and signs, such as `²⁻` or `⁺` | `^` followed by the ASCII digits and sign, such as `^2-` |
| `·` `•` `∙` `*` | component separator `.` |
| `≡` | `#` |
| U+2010, U+2011, U+2012, U+2013, U+2014, U+2212 | `-` |
| Whitespace | removed |

All other characters pass through unchanged. Leading and trailing whitespace is trimmed before folding.

## 3. Lexical layer

The input is split into components at `.` first. Each component is lexed into a **token lattice**: a directed acyclic graph over byte positions in which an edge from `i` to `j` is one token that could occupy that span. The lexer never commits to one tokenization. Edges that lie on no complete path through the component are pruned.

### 3.1 Token kinds

| Kind | Matches | Notes |
| --- | --- | --- |
| Element | An uppercase letter optionally followed by a lowercase letter, matching an element symbol exactly | `C`, `Cl`, `Co`, `Si` |
| Folded element | Two uppercase letters `XY` where `Xy` is an element symbol | `CO` for `Co`, `CL` for `Cl`, `NO` for `No`. Lower prior. |
| Isotope hydrogen | `D` or `T` | Read as hydrogen, warning CF0302 |
| Abbreviation | Any spelling in the loaded dictionary, case-exact | May span digits, letters and hyphens: `CO2Et`, `t-Bu`, `C6H5`, `2-C5H4N`. A spelling that starts with a digit matches only where a locant may start |
| Count | A maximal run of digits | Subscript, multiplier or charge magnitude; the grammar decides |
| Open, close | `(` `[` and `)` `]` | Brackets are parenthesis synonyms |
| Bond | `-` `=` `#` between two units | |
| Marker | `-` `=` `#` at the start or end of a component | An open valence of that order |
| Charge | `+`; `-` at the end of a component; `^` followed by an optional count and a sign | |
| Descriptor | At the start of a component: `(E)-`, `(Z)-`, `(R)-`, `(S)-`, any parenthesized list of those letters or `+` `-` `±` separated by commas followed by `-`, and `cis-`, `trans-`, `rac-` | Skipped, warning CF0301 |
| Locant | At the start of a component, after a descriptor, an open bracket or a bond sign: digits separated by commas, or `o`, `m`, `p`, `ortho`, `meta` or `para`, then `-` and more input | `4-`, `2,4-`, `p-`. Section 6.11. Its digits are not also a count. A longer dictionary spelling that starts there, such as `p-Tol`, replaces it |
| Trailing locant | `-`, then digits separated by commas or a word as for a locant, at the end of a component or before a close bracket | The `-4` of `C6H4OMe-4`. Section 6.11. Its digits are not also a count |

A trailing `-` emits both a Marker edge and a Charge edge. Every other character position must be covered by at least one edge.

### 3.2 Uncovered input

When no complete path exists, the lexer reports the first position that no edge can leave:

- If the uncovered text is an uppercase letter followed by lowercase letters, such as `Xy` or `Mtt` without a dictionary entry, the interpretation is `Underdetermined` with reason `UnknownAbbreviation`.
- Otherwise the interpretation is `Invalid` with CF0101.

### 3.3 Lattice priors

Every edge carries an integer prior. Priors only order readings that survive valence resolution. They never remove one.

| Edge | Prior |
| --- | --- |
| Element in the organic set (H B C N O F Si P S Cl Br I) | 0 |
| Any other element | -2 |
| Folded element | -4, plus -2 if outside the organic set |
| Isotope hydrogen | -1 |
| Abbreviation | +1 |
| Descriptor, locant, count, bond, group, charge | 0 |

## 4. Surface grammar

```
formula    = component , { "." , component } ;
component  = [ descriptor ] , [ count ] , [ marker ] , unit-seq , [ trail ] ;
unit-seq   = ( arene | unit | bracketed-ring ) ,
             { [ bond ] , ( unit | bracketed-ring ) | bond , locant , ring } ,
             [ [ bond ] , numbered ]
           | numbered ;
unit       = ( atom | abbreviation | group ) , [ count ] ;
group      = open , [ bond ] , unit-seq , close ;
arene      = substituent , { "-" , substituent } , [ bond ] , ( ring | phenyl )
           | locant , ring ;
numbered   = ring , "-" , substituent , { "-" , substituent }
           | ring , suffixed , { "-" , locant , suffixed } , trailing-locant
           | ring , trailing-locant ;
bracketed-ring = ring , bracketed , { bracketed } ;
suffixed   = [ bond ] , unit , { [ bond ] , unit } ;
bracketed  = open , locant , unit , { [ bond ] , unit } , close ;
phenyl     = abbreviation ;
substituent = locant , unit , { [ bond ] , unit } ;
ring       = "C6H" , [ "2" | "3" | "4" ] | template ;
template   = "C" , count , [ "H" , [ count ] ] , { element , [ count ] } ;
locant     = ( position , { "," , position } | relative ) , "-" ;
relative   = "o" | "m" | "p" | "ortho" | "meta" | "para" ;
trailing-locant = "-" , ( position , { "," , position } | relative ) ;
position   = "1" | "2" | "3" | "4" | "5" | "6" | "7" | "8" ;
trail      = charge | marker ;
charge     = [ count ] , ( "+" | "-" )
           | "^" , [ count ] , ( "+" | "-" ) ;
marker     = "-" | "=" | "#" ;
bond       = "-" | "=" | "#" ;
open       = "(" | "[" ;
close      = ")" | "]" ;
count      = digit , { digit } ;
```

An arene starts the sequence of a group or of a component, after any leading marker, and an arene with no substituents may also follow a bond inside a sequence. A component count cannot come before an arene, since its digits would run into the locant's. Each substituent ends at a `-` followed by a locant, or at the ring: the first carbon with count 6 directly followed by hydrogen with a count of 4 or less. Units after the ring, including later `C6H4`s, are read as usual. A numbered ring whose substituents follow it unbracketed ends its sequence: its last substituent runs to the end of the group or component. A ring with bracketed substituents does not.

The grammar is deliberately permissive. Parsing is a chart parse over the lattice and its result is a parse forest: one tree per combination of lattice path and grammatical choice. The grammatical choices are:

- A trailing `-` is a charge of -1 or a single-bond marker.
- A count directly before a trailing bare sign is the last unit's count or the charge magnitude, so `Fe2+` parses as `Fe` with charge `2+` and as `Fe2` with charge `+`. A `^` charge never takes the count before it: `Fe2^+` is `Fe2` with charge `+`.

Syntax errors:

| Code | Condition |
| --- | --- |
| CF0001 | Empty formula or empty component |
| CF0002 | Formula longer than 256 bytes |
| CF0003 | Formula needs more work than the limits allow. Each component may take 200,000 search steps by default, where a step is a search state or an atom copied by a group count, and discard up to ten states that cannot complete for each step. Writing canonical SMILES for the readings of the whole formula may visit 25 atoms for each step |
| CF0004 | A reading would have more than 1,000 atoms, counting every copy that a count makes and every component but not hydrogens held by other atoms. The formula fails even if smaller readings exist, so a size limit never decides between readings |
| CF0102 | A group is opened and never closed |
| CF0103 | A close bracket has no matching open |
| CF0104 | Warning: a group opened with one bracket kind is closed with the other |
| CF0105 | A token appears where the grammar does not allow it |
| CF0106 | A group is empty |
| CF0107 | A count of zero or above 100, or a charge of zero or above 8 |
| CF0108 | Locants that do not fit their ring (section 6.11) |

## 5. Valence model

### 5.1 Element valences

Each element has a list of allowed neutral valences. The organic set uses the SMILES defaults: H 1, B 3, C 4, N 3 or 5, O 2, F 1, Si 4, P 3 or 5, S 2, 4 or 6, Cl 1, Br 1, I 1, 3 or 5. Other main-group elements use their group valences. Metals use their common oxidation states as valences. The full table is `crates/moiety-core/src/element.rs`.

A charged atom's valence follows its charge:

- For a nonmetal or metalloid with charge `c`, the valence is the lowest valence of the element with atomic number `Z - c`. N+ has valence 4, O- has 1, C- has 3, B- has 4, Cl- has 0.
- For a metal with charge `c > 0`, each neutral valence `v` with `v >= c` gives valence `v - c`. Fe3+ has valence 0. Metals do not take negative charges.

### 5.2 Bond restrictions

- Multiple bonds are allowed only between two atoms that both belong to B C N O Si P S Se As Te Ge Sb Cl Br I Xe.
- No bond joins two metals.
- Transition metals (groups 3 to 11), lanthanides, actinides and every element with atomic number 104 or more bond only to heteroatoms: never to carbon and never to hydrogen. This is what keeps the cobalt reading of `CH3COOH` from surviving while the cobalt chloride reading of `COCl2` does.
- An aromatic atom from a dictionary fragment takes only single bonds from outside its fragment.
- Nitrogen reaches valence 5 only with a double bond to oxygen, which the charge convention then rewrites as `[N+]-[O-]`.
- A neutral atom using more than its lowest valence, such as sulfur in a sulfone or nitrogen in a nitro group, forms multiple bonds only to oxygen, and at most double ones. This keeps the expanded octet for `S=O`, `P=O` and `N=O` and stops it from inventing `S#S` or `N#N` bonds.

### 5.3 Free valence

For an assembled unit `U`, **free(U)** is the sum of its atoms' chosen valences minus twice the sum of its internal bond orders. Hydrogens are atoms like any other, so `CH3` has free 1, `CH2` free 2, `CH` free 3, `C` free 4, `O` free 2 and `OH` free 1.

## 6. Assembly

Assembly turns one parse tree into zero or more molecules. It processes units from left to right, forking wherever a rule offers a choice, and keeps every branch that ends with all valences satisfied.

### 6.1 The host stack

Assembly keeps a **host stack**: the atoms that later units may bond to, in the order they were created. The **host** for a new attachment is the most recently created atom on the stack that still has free valence. The search for it steps back over saturated atoms, which is how the second oxygen of `COOH` reaches the carbon after the first takes a double bond. It does not step back past a saturated atom that carries written hydrogens: a `CH2` or `NH` written mid-chain is a chain link, so nothing after it bonds to an atom before it. A unit that finds no host while atoms already exist fails with CF0204, "chain cannot continue past a terminal group".

Atoms inside a parenthesized branch are not pushed on the outer stack. Atoms created by a subscript are.

### 6.2 Units

**Atom `X`.** Create `X` with each allowed valence in turn. If a host exists, bond it to `X`. Otherwise `X` starts the structure. Push `X`.

**Hydrogen.** A hydrogen, with or without a count, bonds only to an open atom of the unit written immediately before it; it never steps back further. `H2NOH` is hydroxylamine, not an amine oxide with the last hydrogen on nitrogen. The one exception: when that unit is a run of halogens, the hydrogen bonds to the atom they hang from, so `CF2H` is difluoromethyl.

**Subscripted atom `X n`, n of 2 or more.** A count on an atom is always a subscript: R1 applies to groups and abbreviations, not to atoms.

- If `X` is carbon, the span encodes no connectivity and the tree is underdetermined (section 6.6).
- Otherwise create `n` copies of `X`, each bonded to the same host. The host must have room for all of them, or the tree fails with CF0208. Every copy is pushed, so `CO2H` reads as two oxygens on the carbon with the hydrogen on the open one. A disulfide is written `CH3SSCH3`; `CH3S2CH3` asks for two sulfurs on one methyl carbon and fails.

**Abbreviation.** Insert the dictionary fragment. Its attachment atoms are the neighbors of the fragment's `*` atoms, which are removed.

- One attachment: the attachment atom bonds to the host, or starts the structure. Push it.
- Two attachments: the first attachment atom bonds to the host and the second is pushed. A reversed spelling swaps them.
- No attachments: the abbreviation must be the only unit of its component, otherwise CF0207.
- A count `n` of 2 or more repeats the fragment by R1, treating one attachment as free 1 and two as free 2.

**Group `( unit-seq )` with no count.** Assemble the inner sequence on its own and keep only the assemblies with the least free valence, so `(COOH)` is the carboxyl group and not a peroxide with three open valences. Then bond the inner fragment's first open atom to the host as a branch. A bond written just inside the open bracket, as in `C(=O)`, sets the order of that branch bond. The branch's atoms are not pushed.

**Group with count `n` of 2 or more.** Assemble the inner sequence the same way, then apply R1.

### 6.3 R1, repeated groups

For `(G)n` with `n >= 2`:

- free(G) = 0: CF0206.
- free(G) = 1: the copies are **pendants**, each bonded to the same host, which must have room for all of them.
- free(G) = 2: the copies **catenate**. The first copy's first open atom bonds to the host, the last open atom of each copy bonds to the first open atom of the next, and the last copy's atoms are pushed.
- free(G) of 3 or more: fork into both readings.

This rule alone separates `(CH3)2`, two pendant methyls, from `(CH2)4`, a four-carbon chain segment.

### 6.4 R2, anchors for leading pendants

Pendant-type units that appear before any atom exists (a subscripted atom, a group, or a counted abbreviation with one attachment point) are held as **pending**. The first unit that is not pendant-type becomes the **anchor**: every pending copy bonds to it, and the anchor must have free valence for all of them. If it does not, the tree fails with CF0203.

If the component ends while copies are still pending, the first copy becomes the anchor for the others. This reads `H2`, `Cl2`, `O2` and `(COOH)2` as single molecules.

### 6.5 R3, bond orders

Every bond between two units gets an order. An explicit bond token fixes it. Otherwise each order from 1 up to the smaller remaining free valence of the two atoms is tried, capped at 3 and limited to 1 by section 5.2. Topology and bond order are solved together because the search forks on both at once.

### 6.6 Carbon counts

A carbon atom with a count of 2 or more, such as `C2H5` or `C6H12O6`, carries no connectivity. Unless a dictionary spelling covers the span, the tree is `Underdetermined`. The bundled dictionary spells out the conventional readings of a few formulas: `C2H5` is ethyl, `C6H5` is phenyl and `C6H6` is benzene. With a locant that names the attachment, `2-C5H4N` is 2-pyridyl and `1-C10H7` is 1-naphthyl. Formulas with no single conventional reading, such as `C3H7`, stay underdetermined:

- Reason `MolecularFormula` when the span is the whole component, or when its degree of unsaturation is zero.
- Reason `UnknownRing` when it sits between other units and its degree of unsaturation is above zero, as in the `C6H4` of `CH3C6H4CH3`. Locants number such a ring (section 6.11).

The degree of unsaturation of `Cn Hm` with `k` neighbouring units is `(2n + 2 - m - k) / 2`.

A component that is written as a molecular formula and has no reading is also `Underdetermined` with reason `MolecularFormula`. Written as a molecular formula means: atoms only, every element once, no groups, bonds, abbreviations or charges, and either three or more elements or a count on an atom other than hydrogen. `CH4O`, `H2SO4` and `N2O` fall here; `CH3` does not, and is invalid for its open valence.

### 6.7 Markers and charges

- A leading marker creates an attachment atom first, so the first unit bonds to it.
- A trailing marker bonds an attachment atom to the host after the last unit.
- A written charge `q`, between 1 and 8 in magnitude, is placed during assembly. No single atom carries a charge above 4 in magnitude. Atoms may take a charge of the same sign as `q` while charge remains to be placed. Charges larger than 1 in magnitude are only placed on a component of a single atom. All of `q` must be placed.
- Charge still unplaced after the last unit goes on the atoms left open, when there is exactly one open atom per unit of charge and charging each closes its valence. This reaches atoms inside groups, so `CH3(CH2)3-` has a carbanion reading as `CH3CH2-` does.
- Charges inside dictionary fragments are kept and do not count against `q`.

### 6.8 R4, saturation

A branch survives only if, after the last unit and any trailing marker, every atom has zero free valence and all of the written charge is placed. A branch that ends with free valence fails with CF0202.

### 6.9 R5, salts

When no tree of a component yields a molecule, and the component begins with one or more metal atoms followed by at least one nonmetal unit, with no written charge or markers, the component is read again as a salt. Each metal atom becomes a separate cation whose charge is one of its valences. The remainder is assembled as an anion carrying the opposite total charge. `CuSO4` reads as `[Cu+2].[O-]S(=O)(=O)[O-]`.

### 6.10 Conventions applied after assembly

- **Charge representation.** By default a neutral nitrogen with bond order sum 5 and a double-bonded oxygen is rewritten as `[N+]` with a single-bonded `[O-]`, so `CH3NO2` reads `C[N+](=O)[O-]`. The hypervalent convention does the reverse. The convention in force is recorded on every reading it changes.
- **Hydrogens.** A hydrogen bonded to exactly one heavy atom becomes an implicit hydrogen of that atom.
- **Multi-component.** A leading count repeats its component. Components keep their input order in the molecule.

### 6.11 Numbered benzene rings

A locant numbers the ring of its arene (section 4): a unit written `C6H4`, `C6H3`, `C6H2` or `C6H`, read as benzene, or a substituted phenyl as described below. The ring's **substituted positions** are the ones without hydrogen: two for `C6H4`, three for `C6H3`, four for `C6H2` and five for `C6H`.

**Substituents.** Each substituent is assembled like a group (section 6.2), so the `MeO` of `4-MeOC6H4CHO` is a methoxy group and the `O2N` of `4-O2NC6H4OH` a nitro group. A substituent that is one counted unit, such as `Cl2`, `Me3` or `(MeO)2`, puts one copy at each position its locant lists, so the count must equal the number of positions. Any other substituent lists one position. Each copy bonds to the ring by its first open atom. A bond written before the ring, as in `4-MeO-C6H4`, sets the order of the last substituent's bond.

**Position 1.**

- When no substituent lists position 1, it is left for the units after the ring. The listed positions and position 1 must then be exactly the substituted positions. `4-MeOC6H4CHO` bonds the formyl group at position 1, and `2,4-Cl2C6H3CO2H` and `2-Cl-4-FC6H3CHO` put their substituents at 2 and 4.
- When a substituent lists position 1, the substituents must fill every substituted position, and nothing bonds after the ring: `1,4-Cl2C6H4`, `1,3,5-Me3C6H3`.
- A locant directly before the ring, with no substituent between them, lists the positions left for the units around the ring. With no other locant it must list every substituted position: `1,4-C6H4(OH)2`, `1,3,5-C6H3(CO2H)3`. After other substituents it lists the positions they leave: `2-Cl-1,4-C6H3(OH)2` is 2-chlorohydroquinone.

**Ortho, meta and para.** `o-`, `m-` and `p-`, or `ortho-`, `meta-` and `para-`, place the two substituents of a `C6H4` ring at positions 1 and 2, 1 and 3, or 1 and 4, and are the only locant of their ring. A substituent with one copy takes the second position, so `p-MeC6H4OH` is 4-methylphenol. A substituent with two copies, as in `p-Cl2C6H4`, or a locant directly before the ring, as in `p-C6H4(CO2Me)2`, takes both.

**After the ring.** The positions left for later units are pushed on the host stack lowest on top, so the next unit bonds to the lowest one, and the host search steps back to the next once it is filled. A counted unit whose host is one of these positions puts one copy on each open position of that ring, lowest first, instead of all on one atom: `o-C6H4(CO2H)2` is phthalic acid and `o-C6H4Cl2` is 1,2-dichlorobenzene.

**A locant in a chain.** A locant directly before its ring may also follow a bond inside a sequence, or a leading marker: `HO-p-C6H4-CO2H`, `CH3-o-C6H4-CH3`, `-p-C6H4-`. The ring bonds by its lowest listed position to the unit before it, and its other listed positions are left for the units after it. A ring with substituents before it may follow a leading marker too, which then takes position 1: `-4-MeOC6H4`. A group that starts with a bond works the same way: in `CH3CH(-2,5-C5H3NOH)CH3` the ring bonds to the CH by position 2 and the hydroxy group takes position 5.

**Substituted phenyl.** Before a dictionary spelling of phenyl, one whose fragment is `*c1ccccc1` and whose spelling states no hydrogen count, such as `Ph` but not `C6H5`, locants number a substituted phenyl ring as they number `C6H4`: `4-ClPh` is 4-chlorophenyl, `3,5-(CF3)2PhBr` is 1-bromo-3,5-bis(trifluoromethyl)benzene, and `(4-FPh)3P` is tris(4-fluorophenyl)phosphine. Position 1 is the phenyl's own attachment, so the locants must not list it, and the ring has one substituted position for each listed position and one for position 1. The phenyl is numbered only when no numbered ring such as `C6H4`, or a template ring such as `C5H4N` (section 6.12), comes later in the sequence; otherwise it is part of a substituent, as the benzyl group of `4-CH2PhC6H4OH` is.

**Locants after the ring.** A numbered ring (section 4) is written with its substituents after it, each as `-`, a locant and the substituent: `C6H4-4-OMe`, `C6H3-3,5-(CF3)2`, `C6H3-2-Cl-4-F`. Position 1 is then the ring's bond to what comes before it, either the host or, when the ring starts its group or component, the open valence that group or component leaves: `CH2C6H4-4-OMe`, `P(C6H4-4-OMe)3`, `-C6H4-4-OMe`. The locants must not list position 1, and with it they must number exactly the substituted positions. `o-`, `m-` and `p-` place a single substituent as before: `C6H4-p-OMe`. The units of the last substituent continue the sequence from its positions, so a charge or marker can follow them, as in `MeOC6H4-4-NH3+`. A substituent may itself be a numbered ring, which takes every substituent written after it: `-C6H4-4-C6H4-4-OMe` is 4'-methoxybiphenyl-4-yl.

**Locants after their substituent.** A ring may instead be followed by its substituents with each locant after its own substituent: `C6H4OMe-4`, `C6H4Me-p`, `C6H3Me2-2,6`, `C6H3Cl-2-F-4`, `P(C6H4OMe-4)3`. Between two substituents the locant stands between hyphens, as the `-2-` of `C6H3Cl-2-F-4`. Each substituent runs from the ring, or from the locant before it, to its own locant, and the last locant ends the ring's group or component, so no charge or marker follows it. The last locant numbers the last ring written before it, so in `CH3C6H4C6H4OMe-4` the first ring keeps no locant and stays underdetermined. The substituents read as they do written after the ring, so `C6H4OMe-4` is `C6H4-4-OMe`.

**Bracketed locants.** A ring may also be followed by bracketed substituents that each start with their locant: `C6H4(4-OMe)`, `C6H3(2-Cl)(4-F)`, `C6H4(p-Me)`. The locants must not list position 1, which is the ring's bond to the unit before it or, when the ring starts its group or component, the position left for the units after the brackets: `HOCH2C6H4(4-OMe)`, `C6H4(4-OMe)CH2OH`. A bracket that holds a benzene ring such as `C6H4` of its own, as in the `(4-ClC6H4)` of `C6H4(4-ClC6H4)`, is an ordinary group. A phenyl spelling in the bracket is part of the substituent: `C6H4(4-CH2Ph)Br` is 4-benzylbromobenzene.

Bonds to substituents written before their ring are backward under section 7.4. Bonds to substituents written after it, and to later units, are forward. A ring with no locant stays underdetermined under section 6.6, so `CH3C6H4CH3` still names no isomer.

**Errors.** CF0108 reports a locant that:

- has no ring after it in its sequence, or, written after its substituent, none before it;
- lists a position that is not a single digit from 1 to 8, one its ring does not number, or one position twice;
- lists a different number of positions than its substituent has copies;
- with the other locants of its ring, does not number exactly the substituted positions as above;
- is `o-`, `m-` or `p-`, or spelled out, on a ring without exactly two substituted positions, or beside another locant;
- leaves no open position on its ring while more units follow the ring;
- lists position 1 after the ring, in brackets after it, or before a phenyl spelling.

Not covered: a charge or marker after a trailing locant, and bracketed and trailing locants on one ring, as in `C6H3(4-OMe)Cl-2`.

### 6.12 Numbered template rings

A **ring template** in the dictionary (section 9) defines a ring other than benzene, such as pyridine, thiophene, furan or naphthalene, with its numbering. A template ring is written as the template's formula with the hydrogens left after substitution, in Hill order: `C5H4N` is pyridine with one substituent, `C4H2S` thiophene with two, `C10H7` naphthalene with one. Its **substituted positions** are the numbered positions that hold a hydrogen in the parent ring and have lost it: one for `C5H4N`, two for `C5H3N`. Like benzene, a template ring is numbered only when a locant numbers it; `C5H4NCH2OH` stays underdetermined under section 6.6. When the formulas of two rings both match at one place, the longer is read: with a thienothiophene template, `C6H2S2` is that ring, not benzene followed by sulfur.

The rules of section 6.11 apply, with these differences:

- **No implicit attachment.** Benzene leaves position 1 for the units after it; a template ring leaves nothing. Substituents written before it, with no locant directly before the ring, fill every substituted position and make it a complete ring: `4-Me2NC5H4N` is 4-(dimethylamino)pyridine and `2,6-Me2C5H3N` is 2,6-lutidine, while `4-ClC5H4N` is 4-chloropyridine where `4-ClC6H4` is a 4-chlorophenyl group.
- **Attachment locants.** A locant directly before the ring lists its attachments, as it does for benzene: `2-C5H4NCH2OH` is 2-pyridylmethanol, `5-Br-2-C5H3N` a 5-bromo-2-pyridyl group, and `HOCH2-2,6-C5H3N-CH2OH` 2,6-bis(hydroxymethyl)pyridine.
- **A trailing attachment.** A trailing locant directly after the ring, with no substituent between them, lists the attachment to the unit before the ring or, when the ring starts its group or component, the position that group or component leaves open: `Ph2P(C5H4N-2)`, `(C5H4N-2)3P`, `C10H7-1`.
- **Substituents after the ring.** The ring may be followed by its substituents, each as `-`, a locant and the substituent, as benzene may: `C5H4N-2-CH=NNHPh` puts the hydrazone chain at position 2. Nothing may come before the ring then, since no position is left for it.
- **Positions without hydrogen.** A position whose parent atom holds no hydrogen, such as the nitrogen at position 1 of pyridine or the sulfur of thiophene, cannot be listed.

Bracketed locants, locants after their substituents and `o-`, `m-` and `p-` apply to benzene only. CF0108 reports a template ring whose locants break these rules or those of section 6.11 that apply.

## 7. Interpretation

### 7.1 Kinds

All surviving molecules of all trees are deduplicated by canonical SMILES. The derivation with the highest score is kept for each.

| Result | When |
| --- | --- |
| `Unique` | Exactly one distinct molecule |
| `Ambiguous` | Two or more distinct molecules |
| `Underdetermined` | No molecule, and some tree stopped at section 3.2 or 6.6 |
| `Invalid` | No molecule, and no tree was underdetermined |

A component that is invalid makes the whole formula invalid. Otherwise an underdetermined component makes it underdetermined. Otherwise the readings of the components are combined, and their scores add.

### 7.2 Ambiguity classes

Every reading of an `Ambiguous` result is compared with the top-ranked reading, and the difference is classified by the first rule that applies:

| Class | The two readings differ in |
| --- | --- |
| Symbol segmentation | The tokens chosen inside one lattice region, where no abbreviation competes with an element of the same spelling |
| Abbreviation shadowing | The tokens chosen inside a region where an abbreviation and an element share a spelling |
| Charge vs attachment | The role of a trailing `-` |
| Branch-point fork | An R1 fork, or the valence and bond orders that make a unit a pendant in one reading and a chain link in the other |
| Charge representation | Where a written charge was placed |

Three more classes describe outcomes rather than competing readings:

| Class | Meaning |
| --- | --- |
| Orthographic variant | Different spellings that give the same molecule. Never an ambiguity: they deduplicate. |
| Direction | The direction a group was written in. Recorded on the reading (section 7.4). |
| Topological underdetermination | Reported as `Underdetermined`, never ranked. |

A lattice region is a maximal span between two positions that every complete path passes through, containing more than one edge sequence.

### 7.3 Ranking policy, version 2

A reading's score is the sum of:

- the priors of every lattice edge it used (section 3.3);
- 0 for a chain and -1 for a pendant at an R1 fork;
- 0 for the preferred role of a trailing `-` and -1 for the other. Attachment is preferred when the atom the `-` would sit on is carbon: the last element symbol before the `-`, skipping `H`, digits and whole parenthesized groups, or the attachment atom of a dictionary spelling that ends there. Charge is preferred otherwise, so `CH3C(=O)-` prefers the acetyl group and `CH3COO-` the acetate ion;
- 0 for the preferred reading of a count written just before a trailing sign and -1 for the other. The count is preferably a subscript when it follows `H`, as in `CH2-`, and a charge magnitude otherwise, as in `Fe2+`;
- -1 for every single bond the assembly made between two atoms of N, O, S and Se, so that `CH3SOCH3` ranks the sulfoxide above the sulfenate ester;
- -2 for every written charge that lands on carbon, so `CH3CO-` ranks the acetyl group above the acyl anion;
- -1 for a salt reading under R5.

Readings are ordered by score, highest first, then by canonical SMILES. Ranks start at 1. The policy reorders and never deletes.

Version 2 adds the rule for charges on carbon. Version 1 had no such rule.

### 7.4 Direction

Every bond made between two units, other than a bond to a written hydrogen, is **forward** when a unit bonds to a host on its left, and **backward** when a pending unit bonds to an anchor on its right or when the host is the attachment atom of a reversed dictionary spelling. A reading is `LeftToRight` when all such bonds are forward, `RightToLeft` when all are backward, and `Mixed` otherwise. A component with a single unit is `LeftToRight`.

## 8. Diagnostics

Diagnostics follow the rustc model: a code, a message, a primary span and labelled secondary spans. Spans in diagnostics are byte offsets in the original input.

| Code | Meaning |
| --- | --- |
| CF0101 | Unrecognized character |
| CF0201 | Two atoms cannot form the bond: an explicit bond stronger than they can take, a bond section 5.2 forbids, or a pentavalent nitrogen without a double bond to oxygen |
| CF0202 | Valence left open at the end of the component |
| CF0203 | Leading groups have no anchor with enough free valence |
| CF0204 | Chain cannot continue past a terminal group |
| CF0205 | The written charge cannot be placed |
| CF0206 | A repeated group has no free valence |
| CF0207 | A complete molecule is used as a group |
| CF0208 | More pendants than the host atom has room for |
| CF0301 | Warning: stereo descriptor ignored |
| CF0302 | Warning: isotope label read as hydrogen |
| CF0303 | Warning: reading limit reached, later readings dropped |

When every tree fails, the reported diagnostic is the failure that got furthest into the input, preferring higher lexical scores on ties. A single-edit search over the input then looks for a nearby formula that does read; if one exists it is offered as a note and never applied.

```
error[CF0204]: chain cannot continue past a terminal group
  --> CH3CH3OH
         ^^^ this group has no free valence left
            ^^ but more units follow it
  note: did you mean CH3CH2OH?
```

## 9. Dictionary

The dictionary is a TOML file, versioned separately from the grammar:

```toml
[dictionary]
name = "moiety-default"
version = "1.0.0"

[[entry]]
name = "carboxyl"
smiles = "*C(=O)O"
spellings = ["COOH", "CO2H"]
reversed = ["HOOC", "HO2C"]
```

`smiles` defines the fragment. Each `*` is an attachment point. `spellings` are read with the attachment on the left, `reversed` with the attachment on the right. Callers can load their own file or extend the default one.

A spelling starts with a letter, or with a one- or two-digit locant followed by `-` and a spelling, as `2-C5H4N` names 2-pyridyl and `1-C10H7` 1-naphthyl. A spelling that starts with a digit is matched only where a locant may start (section 3.1), so the count of `(CH2)2-C5H4N` is never read as part of one.

A dictionary may also define ring templates (section 6.12):

```toml
[[ring]]
name = "pyridine"
smiles = "n1ccccc1"
```

`smiles` is the parent ring system, aromatic and with no attachment points. Its ring atoms are numbered from 1 in SMILES order, skipping atoms shared by two rings, so `c1cccc2ccccc12` numbers naphthalene 1 to 8. A numbered atom may hold at most one hydrogen, and an atom shared by two rings none. The template's formula is the parent's in Hill order, written with the hydrogens left, so pyridine is written `C5H4N` to `C5HN`. Two templates may not share a formula, and benzene is built in. The bundled dictionary defines pyridine, thiophene, furan and naphthalene.

## 10. Conformance

The suite in `conformance/cases.toml` lists inputs with their expected result kind, the expected readings in rank order as SMILES, the expected ambiguity classes, and for invalid inputs the expected code and primary span. Expected SMILES are compared after canonicalization, so any valid SMILES spelling is accepted.

Differential testing against another implementation is described in `conformance/differential/README.md`.

Metrics reported by the runner:

| Metric | Definition | v1.0 target |
| --- | --- | --- |
| False-unique rate | Expected-ambiguous cases returned as `Unique` | 0% |
| Round-trip recovery | Molecules written by the inverse writer and read back to the same canonical SMILES at rank 1 | at least 95% |
| Ambiguity class recall | Expected-ambiguous cases returned `Ambiguous` with every expected class | at least 98% |
| Top-1 accuracy | Expected-ambiguous cases whose rank-1 reading is the expected first reading | at least 85% |
| Diagnostic span accuracy | Expected-invalid cases with the expected code and primary span | at least 90% |
| Ambiguity rate | Share of all cases returned `Ambiguous` | reported, no target |
