//! Rings whose substituents or attachments are written after them: after the
//! ring as `C6H4-4-OMe`, with each locant after its substituent as `C6H4OMe-4`,
//! bracketed as `C6H4(4-OMe)`, or a template ring's attachment as `C5H4N-2`
//! (grammar, sections 6.11 and 6.12).

use super::arene::{Ring, Written, mismatch, number, ring_at};
use super::{Item, Parser, SyntaxError, Unit, UnitKind};
use crate::diagnostic::Code;
use crate::lexer::{Edge, EdgeKind};
use crate::molecule::BondOrder;

impl Parser<'_> {
    /// Whether the next tokens are `-` and a locant, which start a substituent
    /// written after its ring, or separate two written before their locants.
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

    /// Parses a ring followed by its substituents or, for a template ring, by the
    /// locant of its attachment. Units that continue the sequence after the ring,
    /// those of a last substituent whose locant comes first, are returned with it.
    pub(super) fn numbered(&mut self) -> Result<Option<(Unit, Vec<Item>)>, SyntaxError> {
        let start = self.pos;
        let Some(ring) = self.ring() else {
            return Ok(None);
        };
        let benzene = ring.template == 0;
        if let Some(edge) = self.peek().filter(|_| !benzene) {
            if let EdgeKind::TrailingLocant(locant) = edge.kind {
                self.pos += 1;
                let written = vec![(edge.span, locant, Vec::new())];
                let arene = number(written, None, ring, self.dictionary, false)?;
                let unit = Unit {
                    kind: UnitKind::Arene(arene),
                    count: None,
                    span: ring.span.cover(edge.span),
                };
                return Ok(Some((unit, Vec::new())));
            }
        }
        if self.locant_follows() {
            return self.locants_first(ring).map(Some);
        }
        if self.trailing_locant_ahead() {
            return self.locants_last(ring).map(Some);
        }
        if self.bracketed_follows() {
            return self.bracketed(ring).map(|unit| Some((unit, Vec::new())));
        }
        self.pos = start;
        Ok(None)
    }

    /// Parses substituents written after the ring, each as `-`, a locant and the
    /// substituent, as the `-4-OMe` of `C6H4-4-OMe`. The last runs to the end of
    /// the sequence, and its units are returned to continue it.
    fn locants_first(&mut self, ring: Ring) -> Result<(Unit, Vec<Item>), SyntaxError> {
        let mut written: Vec<Written> = Vec::new();
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
        let (unit, rest) = self.continued(written, ring)?;
        Ok((unit, rest))
    }

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
                _ if depth == 0 && ring_at(&self.tokens[index..], self.dictionary).is_some() => {
                    return false;
                }
                _ => {}
            }
        }
        false
    }

    /// Parses substituents that each end with their locant, as the `Cl-2-F-4` of
    /// `C6H3Cl-2-F-4`. A locant between two substituents is written `-2-`, and a
    /// trailing locant ends the last one, whose units are returned to continue
    /// the sequence.
    fn locants_last(&mut self, ring: Ring) -> Result<(Unit, Vec<Item>), SyntaxError> {
        let mut written: Vec<Written> = Vec::new();
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
                written.push((edge.span, locant, items));
                return self.continued(written, ring);
            }
            if !items.is_empty() && self.locant_follows() {
                self.pos += 1;
                let (span, locant) = self.take_locant().expect("checked above");
                written.push((span, locant, std::mem::take(&mut items)));
                continue;
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

    /// Numbers a ring whose substituents follow it. The last substituent's units
    /// are returned to continue the sequence from its positions.
    fn continued(
        &self,
        written: Vec<Written>,
        ring: Ring,
    ) -> Result<(Unit, Vec<Item>), SyntaxError> {
        let rest = written
            .last()
            .map(|(_, _, items)| items.clone())
            .unwrap_or_default();
        let mut arene = number(written, None, ring, self.dictionary, true)?;
        let last = arene.substituents.pop().expect("at least one substituent");
        arene.open = last.positions;
        let unit = Unit {
            kind: UnitKind::Arene(arene),
            count: None,
            span: ring.span,
        };
        Ok((unit, rest))
    }

    /// Whether a bracketed substituent that starts with its locant follows, as the
    /// `(4-OMe)` of `C6H4(4-OMe)`. A bracket that holds a benzene ring of its own,
    /// as `(4-ClC6H4)` does, is an ordinary group instead.
    fn bracketed_follows(&self) -> bool {
        let tokens = &self.tokens[self.pos..];
        if !matches!(
            tokens,
            [
                Edge {
                    kind: EdgeKind::Open(_),
                    ..
                },
                Edge {
                    kind: EdgeKind::Locant(_),
                    ..
                },
                ..
            ]
        ) {
            return false;
        }
        let mut depth = 0usize;
        for (index, edge) in tokens.iter().enumerate() {
            match edge.kind {
                EdgeKind::Open(_) => depth += 1,
                EdgeKind::Close(_) => {
                    depth -= 1;
                    if depth == 0 {
                        return true;
                    }
                }
                _ if depth == 1 && ring_at(&tokens[index..], self.dictionary).is_some() => {
                    return false;
                }
                _ => {}
            }
        }
        false
    }

    /// Parses bracketed substituents after a ring, each starting with its locant,
    /// as the `(2-Cl)(4-F)` of `C6H3(2-Cl)(4-F)`. Position 1 of benzene is left for
    /// the bond to the unit before the ring or, when there is none, the units after
    /// the brackets. A template ring's brackets fill every substituted position.
    fn bracketed(&mut self, ring: Ring) -> Result<Unit, SyntaxError> {
        let mut written: Vec<Written> = Vec::new();
        let mut end = ring.span;
        while self.bracketed_follows() {
            let open = self.peek().expect("checked above");
            self.pos += 1;
            let (span, locant) = self.take_locant().expect("checked above");
            let mut items: Vec<Item> = Vec::new();
            loop {
                let Some(edge) = self.peek() else {
                    return Err(Self::unclosed(open));
                };
                if let EdgeKind::Close(closing) = edge.kind {
                    self.pos += 1;
                    end = edge.span;
                    if open.kind != EdgeKind::Open(closing) {
                        self.warnings
                            .push((Code::MixedBrackets, open.span.cover(edge.span)));
                    }
                    break;
                }
                let bond = match edge.kind {
                    EdgeKind::Bond(order) if !items.is_empty() => {
                        self.pos += 1;
                        Some((order, edge.span))
                    }
                    _ => None,
                };
                let Some(unit) = self.unit()? else {
                    return Err(Self::unexpected(edge));
                };
                items.push(Item { bond, unit });
            }
            if items.is_empty() {
                return Err(mismatch(span, "this locant has no substituent after it"));
            }
            written.push((span, locant, items));
        }
        let mut arene = number(written, None, ring, self.dictionary, true)?;
        arene.open = Vec::new();
        Ok(Unit {
            kind: UnitKind::Arene(arene),
            count: None,
            span: ring.span.cover(end),
        })
    }
}
