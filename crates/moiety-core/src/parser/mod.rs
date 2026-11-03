use crate::diagnostic::Code;
use crate::element::Element;
use crate::lexer::{Bracket, Edge, EdgeKind, Lattice, MAX_CHARGE, MAX_COUNT};
use crate::molecule::BondOrder;
use crate::span::Span;

mod arene;

pub(crate) use arene::{Arene, Substituent};

/// One parse tree of a component, built from one path through its lattice.
#[derive(Clone, Debug)]
pub(crate) struct Tree {
    pub span: Span,
    pub path: Vec<usize>,
    pub multiplier: u32,
    pub lead: Option<(BondOrder, Span)>,
    pub seq: Seq,
    pub tail: Tail,
    pub warnings: Vec<(Code, Span)>,
    pub count_as_charge: Option<Span>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Tail {
    None,
    Marker(BondOrder, Span),
    Charge(i32, Span),
}

#[derive(Clone, Debug, Default)]
pub(crate) struct Seq {
    pub items: Vec<Item>,
}

#[derive(Clone, Debug)]
pub(crate) struct Item {
    pub bond: Option<(BondOrder, Span)>,
    pub unit: Unit,
}

#[derive(Clone, Debug)]
pub(crate) struct Unit {
    pub kind: UnitKind,
    pub count: Option<(u32, Span)>,
    pub span: Span,
}

#[derive(Clone, Debug)]
pub(crate) enum UnitKind {
    Atom {
        element: Element,
        symbol: Span,
    },
    Abbreviation {
        entry: usize,
        reversed: bool,
        symbol: Span,
    },
    Group {
        bond: Option<(BondOrder, Span)>,
        seq: Seq,
    },
    Arene(Arene),
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct SyntaxError {
    pub code: Code,
    pub span: Span,
    pub message: String,
}

impl Unit {
    pub fn count(&self) -> u32 {
        self.count.map_or(1, |(n, _)| n)
    }
}

impl Tree {
    /// Whether the component is written as a molecular formula: atoms only, each
    /// element once, and either three elements or a count on an atom other than hydrogen.
    pub fn is_formula_shaped(&self) -> bool {
        if self.multiplier != 1 || self.lead.is_some() || self.tail != Tail::None {
            return false;
        }
        let mut seen = Vec::new();
        let mut heavy_count = false;
        for item in &self.seq.items {
            let UnitKind::Atom { element, .. } = item.unit.kind else {
                return false;
            };
            if item.bond.is_some() || seen.contains(&element) {
                return false;
            }
            heavy_count |= element != Element::H && item.unit.count() > 1;
            seen.push(element);
        }
        seen.len() >= 3 || heavy_count
    }
}

/// Parses one lattice path into its parse trees.
///
/// A path yields two trees when a count sits directly before a trailing sign,
/// since the count may be the last unit's subscript or the charge magnitude.
pub(crate) fn parse(lattice: &Lattice, path: &[usize]) -> Result<Vec<Tree>, SyntaxError> {
    let tokens: Vec<Edge> = path.iter().map(|&i| lattice.edges()[i]).collect();
    let mut parser = Parser {
        tokens: &tokens,
        pos: 0,
        warnings: Vec::new(),
    };
    let tree = parser.component(lattice.span(), path)?;
    let mut trees = vec![tree.clone()];
    let ends_with_sign = tokens
        .last()
        .is_some_and(|e| matches!(e.kind, EdgeKind::Sign(_)));
    if let (Tail::Charge(sign @ (1 | -1), sign_span), true) = (tree.tail, ends_with_sign) {
        if let Some(last) = tree.seq.items.last() {
            if let Some((count, count_span)) = last.unit.count {
                if count_span.end == sign_span.start && count <= MAX_CHARGE {
                    let mut split = tree;
                    let unit = &mut split.seq.items.last_mut().expect("checked above").unit;
                    unit.count = None;
                    unit.span = Span::new(unit.span.start, count_span.start);
                    let magnitude = i32::try_from(count).unwrap_or(i32::MAX);
                    split.tail = Tail::Charge(sign * magnitude, count_span.cover(sign_span));
                    split.count_as_charge = Some(count_span);
                    trees.push(split);
                }
            }
        }
    }
    Ok(trees)
}

struct Parser<'a> {
    tokens: &'a [Edge],
    pos: usize,
    warnings: Vec<(Code, Span)>,
}

impl Parser<'_> {
    fn peek(&self) -> Option<Edge> {
        self.tokens.get(self.pos).copied()
    }

    fn error(code: Code, span: Span, message: &str) -> SyntaxError {
        SyntaxError {
            code,
            span,
            message: message.to_string(),
        }
    }

    fn component(&mut self, span: Span, path: &[usize]) -> Result<Tree, SyntaxError> {
        if let Some(edge) = self.peek().filter(|e| e.kind == EdgeKind::Descriptor) {
            self.warnings.push((Code::StereoIgnored, edge.span));
            self.pos += 1;
        }
        let mut multiplier = 1;
        if let Some(edge) = self.peek() {
            if let EdgeKind::Count(n) = edge.kind {
                multiplier = Self::checked_count(n, edge.span)?;
                self.pos += 1;
            }
        }
        let mut lead = None;
        if let Some(edge) = self.peek() {
            if let EdgeKind::Marker(order) = edge.kind {
                lead = Some((order, edge.span));
                self.pos += 1;
            }
        }
        let seq = self.seq(lead.is_none())?;
        if seq.items.is_empty() {
            return Err(match self.peek() {
                Some(edge) => Self::unexpected(edge),
                None => Self::error(Code::Empty, span, "component has no atoms"),
            });
        }
        let tail = match self.peek() {
            None => Tail::None,
            Some(edge) => {
                let tail = match edge.kind {
                    EdgeKind::Marker(order) => Tail::Marker(order, edge.span),
                    EdgeKind::Sign(sign) => Tail::Charge(i32::from(sign), edge.span),
                    EdgeKind::Charge(charge) => {
                        let magnitude = charge.unsigned_abs();
                        if !(1..=MAX_CHARGE).contains(&magnitude) {
                            return Err(Self::error(
                                Code::CountOutOfRange,
                                edge.span,
                                &format!("charges must be between 1 and {MAX_CHARGE}"),
                            ));
                        }
                        Tail::Charge(charge, edge.span)
                    }
                    EdgeKind::Close(_) => {
                        return Err(Self::error(
                            Code::UnmatchedClose,
                            edge.span,
                            "this bracket closes a group that was never opened",
                        ));
                    }
                    _ => return Err(Self::unexpected(edge)),
                };
                self.pos += 1;
                if let Some(extra) = self.peek() {
                    return Err(Self::unexpected(extra));
                }
                tail
            }
        };
        Ok(Tree {
            span,
            path: path.to_vec(),
            multiplier,
            lead,
            seq,
            tail,
            warnings: std::mem::take(&mut self.warnings),
            count_as_charge: None,
        })
    }

    fn unexpected(edge: Edge) -> SyntaxError {
        if let EdgeKind::TrailingLocant(_) = edge.kind {
            return arene::mismatch(
                edge.span,
                "this locant has no unnumbered benzene ring before it",
            );
        }
        let what = match edge.kind {
            EdgeKind::Count(_) => "a count must follow an atom, abbreviation or group",
            EdgeKind::Bond(_) | EdgeKind::Marker(_) => "a bond must sit between two units",
            EdgeKind::Sign(_) | EdgeKind::Charge(_) => "a charge must end its component",
            EdgeKind::Descriptor => "a descriptor must start its component",
            EdgeKind::Locant(_) => "a locant must start its component, a group or a substituent",
            _ => "unexpected token",
        };
        Self::error(Code::UnexpectedToken, edge.span, what)
    }

    fn checked_count(n: u32, span: Span) -> Result<u32, SyntaxError> {
        if (1..=MAX_COUNT).contains(&n) {
            Ok(n)
        } else {
            Err(Self::error(
                Code::CountOutOfRange,
                span,
                &format!("counts must be between 1 and {MAX_COUNT}"),
            ))
        }
    }

    /// Parses a unit sequence. `arene_allowed` is whether it may start with locants.
    fn seq(&mut self, arene_allowed: bool) -> Result<Seq, SyntaxError> {
        let mut items = Vec::new();
        if arene_allowed {
            if let Some(unit) = self.arene()? {
                items.push(Item { bond: None, unit });
            }
        }
        loop {
            let start = self.pos;
            let bond = match self.peek() {
                Some(edge) if !items.is_empty() => match edge.kind {
                    EdgeKind::Bond(order) => {
                        self.pos += 1;
                        Some((order, edge.span))
                    }
                    _ => None,
                },
                _ => None,
            };
            let (unit, rest) = match self.numbered()? {
                Some((unit, rest)) => (Some(unit), rest),
                None => (self.unit()?, Vec::new()),
            };
            match unit {
                Some(unit) => {
                    if let [Item { unit: first, .. }] = items.as_slice() {
                        if matches!(&first.kind, UnitKind::Arene(arene) if arene.open.is_empty()) {
                            return Err(arene::mismatch(
                                first.span,
                                "this ring has no open position for the units after it",
                            ));
                        }
                    }
                    items.push(Item { bond, unit });
                    items.extend(rest);
                }
                None => {
                    if bond.is_some() {
                        let blamed = match self.peek() {
                            Some(
                                edge @ Edge {
                                    kind: EdgeKind::Locant(_),
                                    ..
                                },
                            ) => edge,
                            _ => self.tokens[start],
                        };
                        return Err(Self::unexpected(blamed));
                    }
                    return Ok(Seq { items });
                }
            }
        }
    }

    fn unit(&mut self) -> Result<Option<Unit>, SyntaxError> {
        let Some(edge) = self.peek() else {
            return Ok(None);
        };
        let kind = match edge.kind {
            EdgeKind::Element { element, .. } => {
                self.pos += 1;
                UnitKind::Atom {
                    element,
                    symbol: edge.span,
                }
            }
            EdgeKind::IsotopeHydrogen => {
                self.pos += 1;
                self.warnings.push((Code::IsotopeIgnored, edge.span));
                UnitKind::Atom {
                    element: Element::H,
                    symbol: edge.span,
                }
            }
            EdgeKind::Abbreviation { entry, reversed } => {
                self.pos += 1;
                UnitKind::Abbreviation {
                    entry,
                    reversed,
                    symbol: edge.span,
                }
            }
            EdgeKind::Open(bracket) => self.group(edge, bracket)?,
            _ => return Ok(None),
        };
        let end = self.tokens[self.pos - 1].span.end;
        let mut unit = Unit {
            kind,
            count: None,
            span: Span::new(edge.span.start, end),
        };
        if let Some(next) = self.peek() {
            if let EdgeKind::Count(n) = next.kind {
                unit.count = Some((Self::checked_count(n, next.span)?, next.span));
                unit.span.end = next.span.end;
                self.pos += 1;
            }
        }
        Ok(Some(unit))
    }

    fn group(&mut self, open: Edge, bracket: Bracket) -> Result<UnitKind, SyntaxError> {
        self.pos += 1;
        let bond = match self.peek() {
            Some(Edge {
                kind: EdgeKind::Bond(order),
                span,
                ..
            }) => {
                self.pos += 1;
                Some((order, span))
            }
            _ => None,
        };
        let seq = self.seq(true)?;
        match self.peek() {
            Some(edge) => match edge.kind {
                EdgeKind::Close(closing) => {
                    if seq.items.is_empty() {
                        return Err(Self::error(
                            Code::EmptyGroup,
                            open.span.cover(edge.span),
                            "this group is empty",
                        ));
                    }
                    if closing != bracket {
                        self.warnings
                            .push((Code::MixedBrackets, open.span.cover(edge.span)));
                    }
                    self.pos += 1;
                    Ok(UnitKind::Group { bond, seq })
                }
                EdgeKind::TrailingLocant(_) => Err(Self::unexpected(edge)),
                _ if seq.items.is_empty() => Err(Self::unexpected(edge)),
                _ => Err(Self::unclosed(open)),
            },
            None => Err(Self::unclosed(open)),
        }
    }

    fn unclosed(open: Edge) -> SyntaxError {
        Self::error(Code::UnclosedGroup, open.span, "this group is never closed")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dictionary::Dictionary;
    use crate::lexer::lex;

    fn trees(text: &str) -> Result<Vec<Tree>, SyntaxError> {
        let lattice = lex(text, Span::new(0, text.len()), &Dictionary::default()).unwrap();
        let (paths, _) = lattice.paths(64);
        let mut all = Vec::new();
        for path in paths {
            all.extend(parse(&lattice, &path)?);
        }
        Ok(all)
    }

    #[test]
    fn parses_groups_and_counts() {
        let all = trees("(CH3)2CHCOOH").unwrap();
        let tree = &all[0];
        assert_eq!(tree.seq.items.len(), 7);
        assert_eq!(tree.seq.items[0].unit.count(), 2);
        assert!(matches!(
            tree.seq.items[0].unit.kind,
            UnitKind::Group { .. }
        ));
    }

    #[test]
    fn forks_count_before_sign() {
        let all = trees("Fe2+").unwrap();
        assert_eq!(all.len(), 2);
        assert_eq!(all[1].tail, Tail::Charge(2, Span::new(2, 4)));
    }

    #[test]
    fn reports_syntax_errors() {
        let error = trees("CH3(CH2").unwrap_err();
        assert_eq!(
            (error.code, error.span),
            (Code::UnclosedGroup, Span::new(3, 4))
        );
        let error = trees("CH3)2").unwrap_err();
        assert_eq!(
            (error.code, error.span),
            (Code::UnmatchedClose, Span::new(3, 4))
        );
        assert_eq!(trees("C()H4").unwrap_err().code, Code::EmptyGroup);
        assert_eq!(
            trees("CH3(CH2)0CH3").unwrap_err().code,
            Code::CountOutOfRange
        );
        assert_eq!(trees("CH3-=CH3").unwrap_err().code, Code::UnexpectedToken);
        let error = trees("CH3-4-ClC6H4").unwrap_err();
        assert_eq!(
            (error.code, error.span),
            (Code::UnexpectedToken, Span::new(4, 6))
        );
    }

    #[test]
    fn numbers_benzene_rings() {
        let all = trees("2-Cl-4-FC6H3CHO").unwrap();
        let UnitKind::Arene(arene) = &all[0].seq.items[0].unit.kind else {
            panic!("expected a numbered ring");
        };
        let positions: Vec<&[u8]> = arene
            .substituents
            .iter()
            .map(|s| s.positions.as_slice())
            .collect();
        assert_eq!(positions, [&[2][..], &[4][..]]);
        assert_eq!(arene.open, [1]);
        assert_eq!(arene.ring, Span::new(8, 12));
        assert_eq!(all[0].seq.items.len(), 4);
        let all = trees("o-C6H4(OH)2").unwrap();
        let UnitKind::Arene(arene) = &all[0].seq.items[0].unit.kind else {
            panic!("expected a numbered ring");
        };
        assert!(arene.substituents.is_empty());
        assert_eq!(arene.open, [1, 2]);
        let error = trees("2,4-ClC6H3OH").unwrap_err();
        assert_eq!(
            (error.code, error.span),
            (Code::LocantMismatch, Span::new(0, 4))
        );
    }

    #[test]
    fn warns_on_mixed_brackets_and_descriptors() {
        let all = trees("(E)-CH3[CH2)2CH3").unwrap();
        let codes: Vec<Code> = all[0].warnings.iter().map(|w| w.0).collect();
        assert_eq!(codes, vec![Code::StereoIgnored, Code::MixedBrackets]);
    }

    #[test]
    fn recognizes_formula_shape() {
        assert!(trees("C2H6O").unwrap()[0].is_formula_shaped());
        assert!(!trees("CH3OH").unwrap()[0].is_formula_shaped());
    }
}
