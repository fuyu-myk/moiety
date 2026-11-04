//! Rings numbered by locants and the substituents written before them, with the
//! position rules every numbered ring shares (grammar, sections 6.11 and 6.12).

use super::{Item, Parser, Seq, SyntaxError, Unit, UnitKind};
use crate::diagnostic::Code;
use crate::dictionary::{Dictionary, Entry};
use crate::element::Element;
use crate::lexer::{Edge, EdgeKind, Locant};
use crate::molecule::{AtomKind, BondOrder};
use crate::span::Span;

/// A ring numbered by locants, with its substituents (grammar, sections 6.11 and 6.12).
#[derive(Clone, Debug)]
pub(crate) struct Arene {
    /// The `C6H4` of the ring.
    pub ring: Span,
    /// The ring's template: 0 for benzene, `i + 1` for [`Dictionary::rings`]`[i]`.
    pub template: usize,
    /// How many atoms the ring has.
    pub size: usize,
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

/// A ring as written: its span, its template, and the hydrogens it has left,
/// which are `None` for a phenyl spelling.
#[derive(Clone, Copy, Debug)]
pub(super) struct Ring {
    pub span: Span,
    pub template: usize,
    pub hydrogens: Option<u8>,
}

/// A ring at the start of a token run, and how many tokens it takes.
pub(super) struct RingAt {
    pub ring: Ring,
    pub taken: usize,
}

impl Parser<'_> {
    pub(super) fn take_locant(&mut self) -> Option<(Span, Locant)> {
        let edge = self.peek()?;
        let EdgeKind::Locant(locant) = edge.kind else {
            return None;
        };
        self.pos += 1;
        Some((edge.span, locant))
    }

    /// Takes a ring: benzene written `C6H4` to `C6H`, or a template ring written
    /// as its formula, such as `C5H4N`.
    pub(super) fn ring(&mut self) -> Option<Ring> {
        let found = ring_at(&self.tokens[self.pos..], self.dictionary)?;
        self.pos += found.taken;
        Some(found.ring)
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
                _ if depth == 0 && ring_at(&self.tokens[index..], self.dictionary).is_some() => {
                    return true;
                }
                _ => {}
            }
        }
        false
    }

    /// Whether the next tokens are a locant directly before a ring, which may
    /// follow a bond inside a sequence, as the `p-C6H4` of `HO-p-C6H4-CO2H`.
    pub(super) fn bare_locant_follows(&self) -> bool {
        matches!(self.peek().map(|e| e.kind), Some(EdgeKind::Locant(_)))
            && ring_at(&self.tokens[self.pos + 1..], self.dictionary).is_some()
    }

    /// Parses a ring numbered by locants, with the substituents written before it.
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
                Some(ring) => Some(ring),
                None if started && !self.ring_ahead() => self.phenyl().map(|span| Ring {
                    span,
                    template: 0,
                    hydrogens: None,
                }),
                None => None,
            };
            if let Some(ring) = ring {
                let arene = number(
                    written,
                    bond.map(|(order, edge)| (order, edge.span)),
                    ring,
                    self.dictionary,
                    false,
                )?;
                return Ok(Some(Unit {
                    kind: UnitKind::Arene(arene),
                    count: None,
                    span: first.0.cover(ring.span),
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
                        mismatch(*span, "this locant has no ring after it")
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

/// The ring at the start of `tokens`: benzene written `C6H4`, `C6H3`, `C6H2` or
/// `C6H`, or a dictionary template written as its formula with one hydrogen or
/// more left and at least one position substituted. When several match, the
/// longest formula wins, as `C6H2S2` over the `C6H2` of benzene.
pub(super) fn ring_at(tokens: &[Edge], dictionary: &Dictionary) -> Option<RingAt> {
    let found = (0..dictionary.templates()).rev().filter_map(|template| {
        let parent = dictionary.template(template);
        let (span, hydrogens, taken) = formula_at(tokens, parent.formula())?;
        let most = if template == 0 {
            4
        } else {
            parent.hydrogens() - 1
        };
        (1..=most)
            .contains(&usize::from(hydrogens))
            .then_some(RingAt {
                ring: Ring {
                    span,
                    template,
                    hydrogens: Some(hydrogens),
                },
                taken,
            })
    });
    found.max_by_key(|at| at.taken)
}

/// The span, hydrogen count and token count of `formula` written in Hill order
/// with a hydrogen count after the carbon, at the start of `tokens`.
fn formula_at(tokens: &[Edge], formula: &[(Element, usize)]) -> Option<(Span, u8, usize)> {
    let element = |at: usize, wanted: Element| {
        matches!(
            tokens.get(at).map(|e| e.kind),
            Some(EdgeKind::Element { element, folded: false }) if element == wanted
        )
    };
    let count = |at: usize| match tokens.get(at).map(|e| e.kind) {
        Some(EdgeKind::Count(n)) => Some(n),
        _ => None,
    };
    let mut at = 0;
    let mut hydrogens = None;
    for &(wanted, n) in formula {
        if !element(at, wanted) {
            return None;
        }
        at += 1;
        if n > 1 {
            if count(at) != u32::try_from(n).ok() {
                return None;
            }
            at += 1;
        }
        if hydrogens.is_none() {
            if !element(at, Element::H) {
                return None;
            }
            at += 1;
            hydrogens = Some(match count(at) {
                Some(k) => {
                    at += 1;
                    k
                }
                None => 1,
            });
        }
    }
    if count(at).is_some() {
        return None;
    }
    let span = tokens[0].span.cover(tokens[at - 1].span);
    Some((span, u8::try_from(hydrogens?).ok()?, at))
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
/// checking them against the ring's template and hydrogens. `after` is whether
/// the substituents are written after the ring. A locant written directly before
/// the ring, the last of `written` with no units, lists its attachments.
pub(super) fn number(
    mut written: Vec<Written>,
    bond: Option<(BondOrder, Span)>,
    ring: Ring,
    dictionary: &Dictionary,
    after: bool,
) -> Result<Arene, SyntaxError> {
    let template = dictionary.template(ring.template);
    let implicit = ring.template == 0;
    let name = match ring.hydrogens {
        Some(h) => template.written(usize::from(h)),
        None => "the phenyl".to_string(),
    };
    let total = |h: u8| template.hydrogens() - usize::from(h);
    let several = written.len() > 1;
    let relative = written
        .iter()
        .find(|(_, locant, _)| matches!(locant, Locant::Relative(_)))
        .map(|(span, _, _)| *span);
    let ortho_meta_para = |substituted: usize| match relative {
        Some(span) if !implicit => Err(mismatch(span, "o-, m- and p- number benzene rings only")),
        Some(span) if several || substituted != 2 => Err(mismatch(
            span,
            "o-, m- and p- number the two substituents of a C6H4 ring, with no other locant",
        )),
        _ => Ok(()),
    };
    if let Some(h) = ring.hydrogens {
        ortho_meta_para(total(h))?;
    }
    let attachment = match written.last() {
        Some((_, _, items)) if items.is_empty() => written.pop(),
        _ => None,
    };
    let mut taken = 0u8;
    let kind = |position: u8| {
        template
            .atom(position)
            .map(|atom| template.parent().atoms()[atom].kind)
    };
    let alternative = |position: u8, taken: u8| {
        let wanted = kind(position)?;
        (1..=8u8)
            .filter(|&p| template.substitutable(p) && taken & (1 << (p - 1)) == 0)
            .filter(|&p| kind(p) == Some(wanted))
            .min_by_key(|&p| (p.abs_diff(position), p))
    };
    let take = |taken: &mut u8, span: Span, positions: &[u8]| -> Result<(), SyntaxError> {
        for &position in positions {
            if !template.substitutable(position) {
                let mut error = mismatch(
                    span,
                    &format!("{name} has no position {position} that can take a substituent"),
                );
                let digit = |p: u8| char::from(b'0' + p);
                error.fix = alternative(position, *taken).map(|p| (digit(position), digit(p)));
                return Err(error);
            }
            let bit = 1 << (position - 1);
            if *taken & bit != 0 {
                return Err(mismatch(
                    span,
                    &format!("position {position} is listed twice"),
                ));
            }
            *taken |= bit;
        }
        Ok(())
    };
    let last = written.len().checked_sub(1);
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
        if positions.contains(&1) && implicit {
            if after {
                return Err(mismatch(
                    span,
                    "position 1 is where the ring bonds to what comes before it",
                ));
            }
            if ring.hydrogens.is_none() {
                return Err(mismatch(span, "position 1 is the phenyl's own attachment"));
            }
        }
        take(&mut taken, span, &positions)?;
        substituents.push(Substituent {
            seq,
            positions,
            bond: if Some(index) == last {
                bond.or(own)
            } else {
                own
            },
            locant: span,
            span: span.cover(end),
        });
    }
    let listed = taken.count_ones() as usize;
    let open = match attachment {
        Some((span, locant, _)) => {
            let h = ring
                .hydrogens
                .expect("a phenyl spelling follows a substituent, never a bare locant");
            let open = positions(span, locant, total(h).saturating_sub(listed))?;
            take(&mut taken, span, &open)?;
            if listed == 0 && open.len() != total(h) {
                return Err(mismatch(
                    span,
                    &format!(
                        "{name} has {} substituted positions, and this locant lists {}",
                        total(h),
                        open.len()
                    ),
                ));
            }
            open
        }
        None if implicit && taken & 1 == 0 => vec![1],
        None => Vec::new(),
    };
    let numbered = listed + open.len();
    let substituted = ring.hydrogens.map_or(numbered, total);
    if ring.hydrogens.is_none() {
        ortho_meta_para(substituted)?;
    }
    if numbered != substituted {
        return Err(mismatch(
            ring.span,
            &format!(
                "{name} has {substituted} substituted positions, and its locants number {numbered}"
            ),
        ));
    }
    Ok(Arene {
        ring: ring.span,
        template: ring.template,
        size: template.parent().atoms().len(),
        substituents,
        open,
        after,
    })
}

/// The positions a locant lists, ascending. `o-`, `m-` and `p-` give one position
/// for a single copy and two otherwise.
fn positions(span: Span, locant: Locant, copies: usize) -> Result<Vec<u8>, SyntaxError> {
    match locant {
        Locant::Positions(0) => Err(mismatch(span, "a locant lists positions 1 to 8, each once")),
        Locant::Positions(bits) => Ok((1..=8).filter(|p| bits & (1 << (p - 1)) != 0).collect()),
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
