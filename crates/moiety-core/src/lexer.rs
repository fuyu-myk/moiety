//! The token lattice: every token that could occupy every span of a component (grammar, section 3).

use crate::dictionary::Dictionary;
use crate::element::Element;
use crate::molecule::BondOrder;
use crate::span::Span;

/// Largest count the lexer accepts.
pub const MAX_COUNT: u32 = 100;

/// Bracket style of a group.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Bracket {
    /// `(` and `)`.
    Round,
    /// `[` and `]`.
    Square,
}

/// What a lattice edge could be.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum EdgeKind {
    /// An element symbol. `folded` marks an all-uppercase spelling such as `CO` for cobalt.
    Element {
        /// The element.
        element: Element,
        /// Whether the symbol was matched by case folding.
        folded: bool,
    },
    /// `D` or `T`, read as hydrogen.
    IsotopeHydrogen,
    /// A dictionary spelling.
    Abbreviation {
        /// Index of the entry in the dictionary.
        entry: usize,
        /// Whether the spelling is a reversed one.
        reversed: bool,
    },
    /// A run of digits.
    Count(u32),
    /// An opening bracket.
    Open(Bracket),
    /// A closing bracket.
    Close(Bracket),
    /// An explicit bond between two units.
    Bond(BondOrder),
    /// An open valence at either end of a component.
    Marker(BondOrder),
    /// A bare `+` or `-` charge sign, `1` or `-1`.
    Sign(i8),
    /// A complete `^` charge such as `^2-`.
    Charge(i32),
    /// A stereo descriptor such as `(E)-`.
    Descriptor,
}

/// One token that could occupy a span of the input.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Edge {
    /// Span in the normalized formula.
    pub span: Span,
    /// What the token is.
    pub kind: EdgeKind,
    /// Ranking prior. Higher is preferred; priors never remove a reading.
    pub prior: i32,
}

/// How the alternatives inside a lattice region compete.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RegionKind {
    /// Different element segmentations, such as `C` + `O` against `Co`.
    Segmentation,
    /// An abbreviation against an element with the same spelling, such as tosyl against tennessine.
    Shadowing,
    /// A trailing `-` read as a charge or as an attachment marker.
    Dash,
}

/// A stretch of input with more than one tokenization.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Region {
    /// Span in the normalized formula.
    pub span: Span,
    /// How the alternatives compete.
    pub kind: RegionKind,
}

/// The token lattice of one component: every token that could occupy every span.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Lattice {
    span: Span,
    edges: Vec<Edge>,
    regions: Vec<Region>,
}

/// Input that no complete tokenization covers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Uncovered {
    /// The first stretch that no path through the lattice can cross.
    pub span: Span,
    /// Whether the stretch looks like an abbreviation: an uppercase letter and lowercase letters.
    pub looks_like_abbreviation: bool,
}

impl Lattice {
    /// Span of the component in the normalized formula.
    #[must_use]
    pub fn span(&self) -> Span {
        self.span
    }

    /// Every edge that lies on at least one complete path.
    #[must_use]
    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }

    /// Stretches with more than one tokenization.
    #[must_use]
    pub fn regions(&self) -> &[Region] {
        &self.regions
    }

    /// Complete paths through the lattice as lists of edge indices, at most `limit` of them.
    /// The flag is set when paths were dropped.
    #[must_use]
    pub fn paths(&self, limit: usize) -> (Vec<Vec<usize>>, bool) {
        let mut paths = Vec::new();
        let mut current = Vec::new();
        let mut truncated = false;
        self.walk(self.span.start, &mut current, &mut paths, limit, &mut truncated);
        (paths, truncated)
    }

    fn walk(
        &self,
        pos: usize,
        current: &mut Vec<usize>,
        paths: &mut Vec<Vec<usize>>,
        limit: usize,
        truncated: &mut bool,
    ) {
        if pos == self.span.end {
            if paths.len() < limit {
                paths.push(current.clone());
            } else {
                *truncated = true;
            }
            return;
        }
        for (index, edge) in self.edges.iter().enumerate() {
            if edge.span.start == pos {
                current.push(index);
                self.walk(edge.span.end, current, paths, limit, truncated);
                current.pop();
                if *truncated {
                    return;
                }
            }
        }
    }
}

/// Builds the token lattice of one component.
///
/// `text` is the whole normalized formula and `span` the component inside it.
///
/// # Errors
///
/// Returns the first stretch of input that no complete tokenization covers.
pub fn lex(text: &str, span: Span, dictionary: &Dictionary) -> Result<Lattice, Uncovered> {
    let component = &text[span.start..span.end];
    let mut edges = Vec::new();
    for (pos, c) in component.char_indices() {
        let at = |len: usize| Span::new(span.start + pos, span.start + pos + len);
        let rest = &component[pos..];
        let first = pos == 0;
        let last = pos + c.len_utf8() == component.len();
        if first {
            if let Some(len) = descriptor_len(rest) {
                edges.push(edge(at(len), EdgeKind::Descriptor));
            }
        }
        for found in dictionary.matches_at(component, pos) {
            edges.push(edge(
                at(found.len),
                EdgeKind::Abbreviation {
                    entry: found.entry,
                    reversed: found.reversed,
                },
            ));
        }
        match c {
            'A'..='Z' => element_edges(rest, at, &mut edges),
            '0'..='9' => {
                let previous_is_digit = pos > 0 && component.as_bytes()[pos - 1].is_ascii_digit();
                if !previous_is_digit {
                    let len = rest.bytes().take_while(u8::is_ascii_digit).count();
                    let value = rest[..len].parse::<u32>().unwrap_or(u32::MAX);
                    edges.push(edge(at(len), EdgeKind::Count(value)));
                }
            }
            '(' => edges.push(edge(at(1), EdgeKind::Open(Bracket::Round))),
            '[' => edges.push(edge(at(1), EdgeKind::Open(Bracket::Square))),
            ')' => edges.push(edge(at(1), EdgeKind::Close(Bracket::Round))),
            ']' => edges.push(edge(at(1), EdgeKind::Close(Bracket::Square))),
            '-' | '=' | '#' => {
                let order = match c {
                    '-' => BondOrder::Single,
                    '=' => BondOrder::Double,
                    _ => BondOrder::Triple,
                };
                if first || last {
                    edges.push(edge(at(1), EdgeKind::Marker(order)));
                } else {
                    edges.push(edge(at(1), EdgeKind::Bond(order)));
                }
                if c == '-' && last && !first {
                    edges.push(edge(at(1), EdgeKind::Sign(-1)));
                }
            }
            '+' => edges.push(edge(at(1), EdgeKind::Sign(1))),
            '^' => {
                if let Some((len, charge)) = caret_charge(rest) {
                    edges.push(edge(at(len), EdgeKind::Charge(charge)));
                }
            }
            _ => {}
        }
    }
    prune(span, edges, component)
}

fn edge(span: Span, kind: EdgeKind) -> Edge {
    Edge {
        span,
        kind,
        prior: prior(kind),
    }
}

fn prior(kind: EdgeKind) -> i32 {
    match kind {
        EdgeKind::Element { element, folded } => {
            let exotic = if element.is_organic() { 0 } else { -2 };
            exotic + if folded { -4 } else { 0 }
        }
        EdgeKind::IsotopeHydrogen => -1,
        EdgeKind::Abbreviation { .. } => 1,
        _ => 0,
    }
}

fn element_edges(rest: &str, at: impl Fn(usize) -> Span, edges: &mut Vec<Edge>) {
    let bytes = rest.as_bytes();
    if let Some(element) = Element::from_symbol(&rest[..1]) {
        edges.push(edge(at(1), EdgeKind::Element { element, folded: false }));
    }
    if matches!(bytes[0], b'D' | b'T') {
        edges.push(edge(at(1), EdgeKind::IsotopeHydrogen));
    }
    let Some(&second) = bytes.get(1) else {
        return;
    };
    if second.is_ascii_lowercase() {
        if let Some(element) = Element::from_symbol(&rest[..2]) {
            edges.push(edge(at(2), EdgeKind::Element { element, folded: false }));
        }
    } else if second.is_ascii_uppercase() {
        let folded = format!("{}{}", char::from(bytes[0]), char::from(second.to_ascii_lowercase()));
        if let Some(element) = Element::from_symbol(&folded) {
            edges.push(edge(at(2), EdgeKind::Element { element, folded: true }));
        }
    }
}

fn caret_charge(rest: &str) -> Option<(usize, i32)> {
    let digits = rest[1..].bytes().take_while(u8::is_ascii_digit).count();
    let magnitude = if digits == 0 {
        1
    } else {
        rest[1..=digits].parse::<i32>().ok()?
    };
    let sign = match rest.as_bytes().get(1 + digits)? {
        b'+' => 1,
        b'-' => -1,
        _ => return None,
    };
    Some((digits + 2, sign * magnitude))
}

fn descriptor_len(rest: &str) -> Option<usize> {
    for prefix in ["cis-", "trans-", "rac-"] {
        if rest.starts_with(prefix) {
            return Some(prefix.len());
        }
    }
    let inner = rest.strip_prefix('(')?;
    let close = inner.find(')')?;
    let body = &inner[..close];
    let allowed = |c: char| matches!(c, 'E' | 'Z' | 'R' | 'S' | '+' | '-' | '±' | ',');
    if body.is_empty() || !body.chars().all(allowed) {
        return None;
    }
    inner[close + 1..]
        .starts_with('-')
        .then_some(close + 3)
}

fn prune(span: Span, edges: Vec<Edge>, component: &str) -> Result<Lattice, Uncovered> {
    let len = span.len();
    let local = |p: usize| p - span.start;
    let mut forward = vec![false; len + 1];
    forward[0] = true;
    for pos in 0..len {
        if forward[pos] {
            for e in edges.iter().filter(|e| local(e.span.start) == pos) {
                forward[local(e.span.end)] = true;
            }
        }
    }
    if !forward[len] {
        let stuck = (0..=len).rev().find(|&p| forward[p]).unwrap_or(0);
        return Err(uncovered(component, span.start, stuck));
    }
    let mut backward = vec![false; len + 1];
    backward[len] = true;
    for pos in (0..len).rev() {
        backward[pos] = edges
            .iter()
            .any(|e| local(e.span.start) == pos && backward[local(e.span.end)]);
    }
    let edges: Vec<Edge> = edges
        .into_iter()
        .filter(|e| forward[local(e.span.start)] && backward[local(e.span.end)])
        .collect();
    let regions = regions(span, &edges);
    Ok(Lattice {
        span,
        edges,
        regions,
    })
}

fn uncovered(component: &str, offset: usize, stuck: usize) -> Uncovered {
    let rest = &component[stuck..];
    let mut chars = rest.char_indices();
    let first = chars.next().map(|(_, c)| c);
    let word = rest
        .char_indices()
        .skip(1)
        .find(|(_, c)| !c.is_ascii_lowercase())
        .map_or(rest.len(), |(i, _)| i);
    let looks_like_abbreviation = first.is_some_and(|c| c.is_ascii_uppercase()) && word > 1;
    let len = if looks_like_abbreviation {
        word
    } else {
        first.map_or(0, char::len_utf8)
    };
    Uncovered {
        span: Span::new(offset + stuck, offset + stuck + len),
        looks_like_abbreviation,
    }
}

fn regions(span: Span, edges: &[Edge]) -> Vec<Region> {
    let cut = |p: usize| !edges.iter().any(|e| e.span.start < p && p < e.span.end);
    let cuts: Vec<usize> = (span.start..=span.end).filter(|&p| cut(p)).collect();
    cuts.windows(2)
        .filter_map(|pair| {
            let region = Span::new(pair[0], pair[1]);
            let inside: Vec<&Edge> = edges
                .iter()
                .filter(|e| region.start <= e.span.start && e.span.end <= region.end)
                .collect();
            if path_count(&inside, region) < 2 {
                return None;
            }
            let same_span = |a: &Edge, b: &Edge| a.span == b.span;
            let kind = if inside.iter().any(|a| {
                matches!(a.kind, EdgeKind::Abbreviation { .. })
                    && inside
                        .iter()
                        .any(|b| same_span(a, b) && matches!(b.kind, EdgeKind::Element { .. }))
            }) {
                RegionKind::Shadowing
            } else if inside.iter().any(|e| matches!(e.kind, EdgeKind::Sign(-1)))
                && inside.iter().any(|e| matches!(e.kind, EdgeKind::Marker(_)))
            {
                RegionKind::Dash
            } else {
                RegionKind::Segmentation
            };
            Some(Region { span: region, kind })
        })
        .collect()
}

fn path_count(inside: &[&Edge], region: Span) -> usize {
    let mut count = vec![0usize; region.len() + 1];
    count[region.len()] = 1;
    for pos in (region.start..region.end).rev() {
        count[pos - region.start] = inside
            .iter()
            .filter(|e| e.span.start == pos)
            .map(|e| count[e.span.end - region.start])
            .fold(0, usize::saturating_add);
    }
    count[0]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dictionary::Entry;

    fn dictionary() -> Dictionary {
        let mut d = Dictionary::new("test", "1");
        let add = |d: &mut Dictionary, name: &str, smiles: &str, spellings: &[&str]| {
            d.insert(
                Entry::new(
                    name,
                    smiles,
                    spellings.iter().map(ToString::to_string).collect(),
                    vec![],
                )
                .unwrap(),
            )
            .unwrap();
        };
        add(&mut d, "tosyl", "*S(=O)(=O)c1ccc(C)cc1", &["Ts"]);
        add(&mut d, "tert-butyl", "*C(C)(C)C", &["t-Bu"]);
        d
    }

    fn lattice(text: &str) -> Lattice {
        lex(text, Span::new(0, text.len()), &dictionary()).unwrap()
    }

    #[test]
    fn keeps_both_readings_of_co() {
        let l = lattice("COCl2");
        let cobalt = Element::from_symbol("Co").unwrap();
        assert!(l.edges().iter().any(|e| e.kind == EdgeKind::Element { element: cobalt, folded: true }));
        assert_eq!(l.paths(10).0.len(), 2);
        assert_eq!(l.regions(), &[Region { span: Span::new(0, 2), kind: RegionKind::Segmentation }]);
    }

    #[test]
    fn marks_abbreviation_shadowing() {
        let l = lattice("TsCl");
        assert_eq!(l.regions()[0].kind, RegionKind::Shadowing);
        assert_eq!(l.paths(10).0.len(), 2);
    }

    #[test]
    fn trailing_dash_is_charge_or_marker() {
        let l = lattice("CH3COO-");
        assert_eq!(l.regions().last().unwrap().kind, RegionKind::Dash);
        let leading = lattice("-CH3");
        assert!(matches!(leading.edges()[0].kind, EdgeKind::Marker(BondOrder::Single)));
    }

    #[test]
    fn prunes_dead_edges() {
        let l = lattice("Cl");
        assert_eq!(l.edges().len(), 1);
        assert!(l.regions().is_empty());
    }

    #[test]
    fn reads_descriptors_counts_and_charges() {
        let l = lattice("(E)-CH3CH=CHCH3");
        assert_eq!(l.edges()[0].kind, EdgeKind::Descriptor);
        let l = lattice("SO4^2-");
        assert!(l.edges().iter().any(|e| e.kind == EdgeKind::Charge(-2)));
        let l = lattice("CH3(CH2)12CH3");
        assert!(l.edges().iter().any(|e| e.kind == EdgeKind::Count(12)));
        assert!(!l.edges().iter().any(|e| e.kind == EdgeKind::Count(1)));
        assert!(lattice("t-BuOH").paths(10).0.len() == 1);
    }

    #[test]
    fn reports_uncovered_input() {
        let err = lex("CH3Xy", Span::new(0, 5), &dictionary()).unwrap_err();
        assert_eq!(err, Uncovered { span: Span::new(3, 5), looks_like_abbreviation: true });
        let err = lex("C$H4", Span::new(0, 4), &dictionary()).unwrap_err();
        assert_eq!(err, Uncovered { span: Span::new(1, 2), looks_like_abbreviation: false });
    }
}
