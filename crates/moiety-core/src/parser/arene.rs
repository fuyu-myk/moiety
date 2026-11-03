//! Benzene rings numbered by locants and the substituents written before them,
//! with the position rules every numbered ring shares (grammar, section 6.11).

use super::{Item, Parser, Seq, SyntaxError, Unit, UnitKind};
use crate::diagnostic::Code;
use crate::dictionary::{Dictionary, Entry};
use crate::element::Element;
use crate::lexer::{Edge, EdgeKind, Locant};
use crate::molecule::{AtomKind, BondOrder};
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
    /// The locant written with the substituent.
    pub locant: Span,
    pub span: Span,
}

/// A substituent as written: its locant and its units.
pub(super) type Written = (Span, Locant, Vec<Item>);

impl Parser<'_> {
    pub(super) fn take_locant(&mut self) -> Option<(Span, Locant)> {
        let edge = self.peek()?;
        let EdgeKind::Locant(locant) = edge.kind else {
            return None;
        };
        self.pos += 1;
        Some((edge.span, locant))
    }

    /// Takes a benzene ring written `C6H4`, `C6H3`, `C6H2` or `C6H`, returning its
    /// span and hydrogen count.
    pub(super) fn ring(&mut self) -> Option<(Span, u8)> {
        let (span, hydrogens, taken) = ring_at(&self.tokens[self.pos..])?;
        self.pos += taken;
        Some((span, hydrogens))
    }

    /// Takes a dictionary spelling of phenyl, such as `Ph`, returning its span.
    fn phenyl(&mut self) -> Option<Span> {
        let edge = self.peek()?;
        is_phenyl_at(edge, self.dictionary, self.text).then(|| {
            self.pos += 1;
            edge.span
        })
    }

    /// Whether a benzene ring such as `C6H4` comes later in the current sequence,
    /// in which case a phenyl spelling before it is part of a substituent.
    fn ring_ahead(&self) -> bool {
        let mut depth = 0usize;
        for (index, edge) in self.tokens.iter().enumerate().skip(self.pos) {
            match edge.kind {
                EdgeKind::Open(_) => depth += 1,
                EdgeKind::Close(_) if depth == 0 => return false,
                EdgeKind::Close(_) => depth -= 1,
                _ if depth == 0 && ring_at(&self.tokens[index..]).is_some() => return true,
                _ => {}
            }
        }
        false
    }

    /// Whether the next tokens are a locant directly before a benzene ring, which
    /// may follow a bond inside a sequence, as the `p-C6H4` of `HO-p-C6H4-CO2H`.
    pub(super) fn bare_locant_follows(&self) -> bool {
        matches!(self.peek().map(|e| e.kind), Some(EdgeKind::Locant(_)))
            && ring_at(&self.tokens[self.pos + 1..]).is_some()
    }

    /// Parses a benzene ring numbered by locants, with the substituents written before it.
    pub(super) fn arene(&mut self) -> Result<Option<Unit>, SyntaxError> {
        let Some(first) = self.take_locant() else {
            return Ok(None);
        };
        let mut written: Vec<Written> = vec![(first.0, first.1, Vec::new())];
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
            let ring = match self.ring() {
                Some((span, hydrogens)) => Some((span, Some(hydrogens))),
                None if started && !self.ring_ahead() => self.phenyl().map(|span| (span, None)),
                None => None,
            };
            if let Some((ring, hydrogens)) = ring {
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

/// The span, hydrogen count and token count of a benzene ring written `C6H4`,
/// `C6H3`, `C6H2` or `C6H` at the start of `tokens`.
pub(super) fn ring_at(tokens: &[Edge]) -> Option<(Span, u8, usize)> {
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

/// Whether an edge is a forward spelling of a dictionary entry for phenyl that
/// states no hydrogen count, as `Ph` does and `C6H5` does not.
fn is_phenyl_at(edge: Edge, dictionary: &Dictionary, text: &str) -> bool {
    let EdgeKind::Abbreviation {
        entry,
        reversed: false,
    } = edge.kind
    else {
        return false;
    };
    let counted = text[edge.span.start..edge.span.end]
        .bytes()
        .any(|b| b.is_ascii_digit());
    !counted && dictionary.entries().get(entry).is_some_and(is_phenyl)
}

/// Whether an entry's fragment is phenyl, `*c1ccccc1`: six aromatic carbons in
/// one ring, each with a hydrogen except the attached one.
fn is_phenyl(entry: &Entry) -> bool {
    let fragment = entry.fragment();
    let [attachment] = entry.attachments() else {
        return false;
    };
    let atoms = fragment.atoms();
    let carbons = atoms.len() == 6
        && atoms.iter().enumerate().all(|(index, atom)| {
            atom.kind == AtomKind::Element(Element::C)
                && atom.aromatic
                && atom.charge == 0
                && atom.hydrogens == u8::from(index != attachment.atom)
        });
    if !carbons
        || fragment.bonds().len() != 6
        || fragment
            .bonds()
            .iter()
            .any(|b| b.order != BondOrder::Aromatic)
    {
        return false;
    }
    let adjacency = fragment.adjacency();
    let mut seen = vec![false; atoms.len()];
    let mut stack = vec![0];
    while let Some(atom) = stack.pop() {
        if !std::mem::replace(&mut seen[atom], true) {
            stack.extend(adjacency[atom].iter().map(|&(next, _)| next));
        }
    }
    adjacency.iter().all(|next| next.len() == 2) && seen.iter().all(|&s| s)
}

pub(super) fn mismatch(span: Span, message: &str) -> SyntaxError {
    Parser::error(Code::LocantMismatch, span, message)
}

/// Gives ring positions to a ring's substituents and to the units around it,
/// checking them against the ring's hydrogens, which are `None` for a phenyl
/// spelling. `after` is whether the substituents are written after the ring,
/// where position 1 is never listed.
pub(super) fn number(
    written: Vec<Written>,
    bond: Option<(BondOrder, Span)>,
    ring: Span,
    hydrogens: Option<u8>,
    after: bool,
) -> Result<Arene, SyntaxError> {
    let name = match hydrogens {
        Some(1) => "C6H".to_string(),
        Some(h) => format!("C6H{h}"),
        None => "the phenyl".to_string(),
    };
    let several = written.len() > 1;
    let relative = written
        .iter()
        .find(|(_, locant, _)| matches!(locant, Locant::Relative(_)))
        .map(|(span, _, _)| *span);
    let ortho_meta_para = |substituted: usize| match relative {
        Some(span) if several || substituted != 2 => Err(mismatch(
            span,
            "o-, m- and p- number the two substituents of a C6H4 ring, with no other locant",
        )),
        _ => Ok(()),
    };
    if let Some(h) = hydrogens {
        ortho_meta_para(6 - usize::from(h))?;
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
        let Some(h) = hydrogens else {
            return Err(mismatch(span, "this locant has no substituent after it"));
        };
        let substituted = 6 - usize::from(h);
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
        if positions.contains(&1) {
            if after {
                return Err(mismatch(
                    span,
                    "position 1 is where the ring bonds to what comes before it",
                ));
            }
            if hydrogens.is_none() {
                return Err(mismatch(span, "position 1 is the phenyl's own attachment"));
            }
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
    let substituted = hydrogens.map_or(numbered, |h| 6 - usize::from(h));
    if hydrogens.is_none() {
        ortho_meta_para(substituted)?;
    }
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
