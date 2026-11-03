//! Benzene rings numbered by locants and the substituents written before or
//! after them (grammar, section 6.11).

use super::{Item, Parser, Seq, SyntaxError, Unit, UnitKind};
use crate::diagnostic::Code;
use crate::element::Element;
use crate::lexer::{Edge, EdgeKind, Locant};
use crate::molecule::BondOrder;
use crate::span::Span;

/// A benzene ring numbered by locants, with its substituents (grammar, section 6.11).
#[derive(Clone, Debug)]
pub(crate) struct Arene {
    /// The `C6H4` of the ring.
    pub ring: Span,
    pub substituents: Vec<Substituent>,
    /// Positions left for the units after the ring, lowest first. When the
    /// substituents follow the ring, these are the last substituent's, whose
    /// units continue the sequence.
    pub open: Vec<u8>,
    /// Whether the substituents are written after the ring.
    pub after: bool,
}

/// A substituent of a numbered ring, assembled like a group.
#[derive(Clone, Debug)]
pub(crate) struct Substituent {
    pub seq: Seq,
    /// The ring position of each copy.
    pub positions: Vec<u8>,
    /// A bond written between the substituent and the ring.
    pub bond: Option<(BondOrder, Span)>,
    /// The locant written before the substituent.
    pub locant: Span,
    pub span: Span,
}

impl Parser<'_> {
    fn take_locant(&mut self) -> Option<(Span, Locant)> {
        let edge = self.peek()?;
        let EdgeKind::Locant(locant) = edge.kind else {
            return None;
        };
        self.pos += 1;
        Some((edge.span, locant))
    }

    /// Takes a benzene ring written `C6H4`, `C6H3`, `C6H2` or `C6H`, returning its
    /// span and hydrogen count.
    fn ring(&mut self) -> Option<(Span, u8)> {
        let (span, hydrogens, taken) = ring_at(&self.tokens[self.pos..])?;
        self.pos += taken;
        Some((span, hydrogens))
    }

    /// Parses a benzene ring numbered by locants, with the substituents written before it.
    pub(super) fn arene(&mut self) -> Result<Option<Unit>, SyntaxError> {
        let Some(first) = self.take_locant() else {
            return Ok(None);
        };
        let mut written = vec![(first.0, first.1, Vec::<Item>::new())];
        loop {
            let started = written
                .last()
                .is_some_and(|(_, _, items)| !items.is_empty());
            let bond = match self.peek() {
                Some(edge) if started => match edge.kind {
                    EdgeKind::Bond(order) => {
                        self.pos += 1;
                        Some((order, edge))
                    }
                    _ => None,
                },
                _ => None,
            };
            if let Some((ring, hydrogens)) = self.ring() {
                let arene = number(
                    written,
                    bond.map(|(order, edge)| (order, edge.span)),
                    ring,
                    hydrogens,
                    false,
                )?;
                return Ok(Some(Unit {
                    kind: UnitKind::Arene(arene),
                    count: None,
                    span: first.0.cover(ring),
                }));
            }
            if let Some((order, edge)) = bond {
                if let Some((span, locant)) = self.take_locant() {
                    if order != BondOrder::Single {
                        return Err(Self::unexpected(edge));
                    }
                    written.push((span, locant, Vec::new()));
                    continue;
                }
            }
            let Some(unit) = self.unit()? else {
                return Err(match bond {
                    Some((_, edge)) => Self::unexpected(edge),
                    None => {
                        let (span, _, _) = written.last().expect("never empty");
                        mismatch(*span, "this locant has no benzene ring after it")
                    }
                });
            };
            let bond = bond.map(|(order, edge)| (order, edge.span));
            written
                .last_mut()
                .expect("never empty")
                .2
                .push(Item { bond, unit });
        }
    }
}

impl Parser<'_> {
    /// Whether the next tokens are `-` and a locant, which start a substituent
    /// written after its ring.
    fn locant_follows(&self) -> bool {
        matches!(
            self.tokens.get(self.pos..self.pos + 2),
            Some([
                Edge {
                    kind: EdgeKind::Bond(BondOrder::Single),
                    ..
                },
                Edge {
                    kind: EdgeKind::Locant(_),
                    ..
                },
            ])
        )
    }

    /// Parses a benzene ring followed by its substituents, each written `-`, a
    /// locant and the substituent. The last substituent runs to the end of the
    /// sequence, and its units are returned to continue it after the ring.
    pub(super) fn numbered(&mut self) -> Result<Option<(Unit, Vec<Item>)>, SyntaxError> {
        let start = self.pos;
        let Some((ring, hydrogens)) = self.ring() else {
            return Ok(None);
        };
        if !self.locant_follows() {
            if self.trailing_locant_ahead() {
                return self.suffixed(ring, hydrogens).map(Some);
            }
            self.pos = start;
            return Ok(None);
        }
        let mut written = Vec::new();
        while self.locant_follows() {
            self.pos += 1;
            let (span, locant) = self.take_locant().expect("checked above");
            let mut items: Vec<Item> = Vec::new();
            loop {
                if !items.is_empty() && self.locant_follows() {
                    break;
                }
                let bond = match self.peek() {
                    Some(edge) if !items.is_empty() => match edge.kind {
                        EdgeKind::Bond(order) => {
                            self.pos += 1;
                            Some((order, edge))
                        }
                        _ => None,
                    },
                    _ => None,
                };
                let (unit, rest) = match self.numbered()? {
                    Some((unit, rest)) => (Some(unit), rest),
                    None => (self.unit()?, Vec::new()),
                };
                let Some(unit) = unit else {
                    if let Some((_, edge)) = bond {
                        return Err(Self::unexpected(edge));
                    }
                    if items.is_empty() {
                        return Err(mismatch(span, "this locant has no substituent after it"));
                    }
                    break;
                };
                let bond = bond.map(|(order, edge)| (order, edge.span));
                items.push(Item { bond, unit });
                items.extend(rest);
            }
            written.push((span, locant, items));
        }
        let rest = written
            .last()
            .map(|(_, _, items)| items.clone())
            .unwrap_or_default();
        let mut arene = number(written, None, ring, hydrogens, true)?;
        let last = arene.substituents.pop().expect("at least one substituent");
        arene.open = last.positions;
        let unit = Unit {
            kind: UnitKind::Arene(arene),
            count: None,
            span: ring.cover(last.locant),
        };
        Ok(Some((unit, rest)))
    }
}

impl Parser<'_> {
    /// Whether a trailing locant, such as the `-4` of `C6H4OMe-4`, ends the current
    /// sequence with no other ring before it, so that it numbers the ring just taken.
    fn trailing_locant_ahead(&self) -> bool {
        let mut depth = 0usize;
        for (index, edge) in self.tokens.iter().enumerate().skip(self.pos) {
            match edge.kind {
                EdgeKind::Open(_) => depth += 1,
                EdgeKind::Close(_) if depth == 0 => return false,
                EdgeKind::Close(_) => depth -= 1,
                EdgeKind::TrailingLocant(_) if depth == 0 => return true,
                _ if depth == 0 && ring_at(&self.tokens[index..]).is_some() => return false,
                _ => {}
            }
        }
        false
    }

    /// Parses the substituent after a ring up to its trailing locant, as the
    /// `OMe-4` of `C6H4OMe-4`. Its units are returned to continue the sequence.
    fn suffixed(&mut self, ring: Span, hydrogens: u8) -> Result<(Unit, Vec<Item>), SyntaxError> {
        let mut items: Vec<Item> = Vec::new();
        loop {
            let edge = self.peek().expect("a trailing locant lies ahead");
            if let EdgeKind::TrailingLocant(locant) = edge.kind {
                self.pos += 1;
                if items.is_empty() {
                    return Err(mismatch(
                        edge.span,
                        "this locant has no substituent before it",
                    ));
                }
                let written = vec![(edge.span, locant, items.clone())];
                let mut arene = number(written, None, ring, hydrogens, true)?;
                let last = arene.substituents.pop().expect("one substituent");
                arene.open = last.positions;
                let unit = Unit {
                    kind: UnitKind::Arene(arene),
                    count: None,
                    span: ring,
                };
                return Ok((unit, items));
            }
            let bond = match edge.kind {
                EdgeKind::Bond(order) => {
                    self.pos += 1;
                    Some((order, edge.span))
                }
                _ => None,
            };
            let Some(unit) = self.unit()? else {
                let next = self.peek().expect("a trailing locant lies ahead");
                let blamed = match (bond, next.kind) {
                    (Some(_), EdgeKind::Locant(_)) | (None, _) => next,
                    (Some(_), _) => edge,
                };
                return Err(Self::unexpected(blamed));
            };
            items.push(Item { bond, unit });
        }
    }
}

/// The span, hydrogen count and token count of a benzene ring written `C6H4`,
/// `C6H3`, `C6H2` or `C6H` at the start of `tokens`.
fn ring_at(tokens: &[Edge]) -> Option<(Span, u8, usize)> {
    let [
        Edge {
            kind:
                EdgeKind::Element {
                    element: Element::C,
                    folded: false,
                },
            span: carbon,
            ..
        },
        Edge {
            kind: EdgeKind::Count(6),
            ..
        },
        Edge {
            kind:
                EdgeKind::Element {
                    element: Element::H,
                    folded: false,
                },
            span: hydrogen,
            ..
        },
        rest @ ..,
    ] = tokens
    else {
        return None;
    };
    let (end, hydrogens, taken) = match rest.first().map(|e| (e.kind, e.span)) {
        Some((EdgeKind::Count(n @ 2..=4), span)) => (span, n, 4),
        Some((EdgeKind::Count(_), _)) => return None,
        _ => (*hydrogen, 1, 3),
    };
    Some((
        carbon.cover(end),
        u8::try_from(hydrogens).expect("1 to 4"),
        taken,
    ))
}

pub(super) fn mismatch(span: Span, message: &str) -> SyntaxError {
    Parser::error(Code::LocantMismatch, span, message)
}

/// Gives ring positions to a ring's substituents and to the units around it,
/// checking them against the ring's hydrogens. `after` is whether the
/// substituents are written after the ring, where position 1 is never listed.
fn number(
    written: Vec<(Span, Locant, Vec<Item>)>,
    bond: Option<(BondOrder, Span)>,
    ring: Span,
    hydrogens: u8,
    after: bool,
) -> Result<Arene, SyntaxError> {
    let substituted = 6 - usize::from(hydrogens);
    let name = if hydrogens == 1 {
        "C6H".to_string()
    } else {
        format!("C6H{hydrogens}")
    };
    let several = written.len() > 1;
    for (span, locant, _) in &written {
        if matches!(locant, Locant::Relative(_)) && (several || hydrogens != 4) {
            return Err(mismatch(
                *span,
                "o-, m- and p- number the two substituents of a C6H4 ring, with no other locant",
            ));
        }
    }
    let bare = written.iter().position(|(_, _, items)| items.is_empty());
    if let Some(index) = bare {
        let (span, locant, _) = written[index];
        if several {
            return Err(mismatch(
                span,
                "a locant directly before the ring must be its only locant",
            ));
        }
        let open = positions(span, locant, substituted)?;
        if open.len() != substituted {
            return Err(mismatch(
                span,
                &format!(
                    "{name} has {substituted} substituted positions, and this locant lists {}",
                    open.len()
                ),
            ));
        }
        return Ok(Arene {
            ring,
            substituents: Vec::new(),
            open,
            after,
        });
    }
    let last = written.len() - 1;
    let mut taken = 0u8;
    let mut substituents = Vec::new();
    for (index, (span, locant, items)) in written.into_iter().enumerate() {
        let end = items.last().expect("not bare").unit.span;
        let (seq, copies, own) = substituent(items);
        let positions = positions(span, locant, copies)?;
        if positions.len() != copies {
            let noun = |n: usize, one: &str, many: &str| {
                if n == 1 {
                    format!("1 {one}")
                } else {
                    format!("{n} {many}")
                }
            };
            return Err(mismatch(
                span,
                &format!(
                    "this locant lists {} for {}",
                    noun(positions.len(), "position", "positions"),
                    noun(copies, "copy", "copies")
                ),
            ));
        }
        if after && positions.contains(&1) {
            return Err(mismatch(
                span,
                "position 1 is where the ring bonds to what comes before it",
            ));
        }
        for &position in &positions {
            let bit = 1 << (position - 1);
            if taken & bit != 0 {
                return Err(mismatch(
                    span,
                    &format!("position {position} is listed twice"),
                ));
            }
            taken |= bit;
        }
        substituents.push(Substituent {
            seq,
            positions,
            bond: if index == last { bond.or(own) } else { own },
            locant: span,
            span: span.cover(end),
        });
    }
    let open = if taken & 1 == 0 { vec![1] } else { Vec::new() };
    let numbered = taken.count_ones() as usize + open.len();
    if numbered != substituted {
        return Err(mismatch(
            ring,
            &format!(
                "{name} has {substituted} substituted positions, and its locants number {numbered}"
            ),
        ));
    }
    Ok(Arene {
        ring,
        substituents,
        open,
        after,
    })
}

/// The positions a locant lists, ascending. `o-`, `m-` and `p-` give one position
/// for a single copy and two otherwise.
fn positions(span: Span, locant: Locant, copies: usize) -> Result<Vec<u8>, SyntaxError> {
    match locant {
        Locant::Positions(0) => Err(mismatch(
            span,
            "a benzene ring is numbered 1 to 6, each position once",
        )),
        Locant::Positions(bits) => Ok((1..=6).filter(|p| bits & (1 << (p - 1)) != 0).collect()),
        Locant::Relative(partner) if copies == 1 => Ok(vec![partner]),
        Locant::Relative(partner) => Ok(vec![1, partner]),
    }
}

/// A substituent's sequence, how many copies of it are written, and any bond
/// written at the start of its group.
fn substituent(items: Vec<Item>) -> (Seq, usize, Option<(BondOrder, Span)>) {
    if let [item] = items.as_slice() {
        let copies = usize::try_from(item.unit.count()).unwrap_or(usize::MAX);
        if let UnitKind::Group { bond, seq } = &item.unit.kind {
            return (seq.clone(), copies, *bond);
        }
        let unit = Unit {
            count: None,
            ..item.unit.clone()
        };
        return (
            Seq {
                items: vec![Item { bond: None, unit }],
            },
            copies,
            None,
        );
    }
    (Seq { items }, 1, None)
}
