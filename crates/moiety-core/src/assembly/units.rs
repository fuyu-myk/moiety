use super::state::{AAtom, State};
use super::{Ctx, Failure, assemble_seq, finish_pending};
use crate::diagnostic::Code;
use crate::element::Element;
use crate::molecule::BondOrder;
use crate::parser::{Unit, UnitKind};
use crate::span::Span;

const COMBINATION_LIMIT: usize = 256;

/// How a unit's atoms join the structure built so far.
pub(crate) enum Piece {
    /// `head` bonds to the host, or starts the structure; `push` joins the host stack.
    Link {
        head: usize,
        push: Vec<usize>,
        bond: Option<(BondOrder, Span)>,
    },
    /// Every head bonds to the same host, or waits for an anchor; `push` joins the host stack.
    Pendants {
        heads: Vec<usize>,
        push: Vec<usize>,
        bond: Option<(BondOrder, Span)>,
    },
    /// A complete molecule that must stand alone.
    Whole,
}

/// Every way to create a unit's atoms in `state`, with how each joins the structure.
pub(crate) fn expand(state: &State, unit: &Unit, ctx: &mut Ctx<'_>) -> Vec<(State, Piece)> {
    match &unit.kind {
        UnitKind::Atom { element, symbol } => atom(state, unit, *element, *symbol, ctx),
        UnitKind::Abbreviation {
            entry,
            reversed,
            symbol,
        } => abbreviation(state, unit, *entry, *reversed, *symbol, ctx),
        UnitKind::Group { bond, seq, .. } => {
            let whole_allowed = std::mem::replace(&mut ctx.whole_allowed, false);
            let inner = assemble_seq(seq, ctx, State::default(), None);
            ctx.whole_allowed = whole_allowed;
            let mut resolved: Vec<State> = Vec::new();
            for s in inner {
                resolved.extend(finish_pending(s, ctx));
            }
            let least = resolved
                .iter()
                .map(|s| s.live_free())
                .filter(|&free| free > 0)
                .min();
            let inner: Vec<State> = resolved
                .into_iter()
                .filter(|s| least.is_none_or(|least| s.live_free() == least))
                .collect();
            inner
                .iter()
                .flat_map(|group| repeat_group(state, group, unit, *bond, ctx))
                .collect()
        }
    }
}

fn options(element: Element, charge_left: i32, single: bool) -> Vec<(i8, u8)> {
    let mut out: Vec<(i8, u8)> = element.valences().iter().map(|&v| (0, v)).collect();
    if charge_left != 0 && (single || element != Element::H) {
        let charge = if single {
            i8::try_from(charge_left).unwrap_or(0)
        } else if charge_left > 0 {
            1
        } else {
            -1
        };
        if charge != 0 {
            out.extend(
                element
                    .charged_valences(charge)
                    .into_iter()
                    .map(|v| (charge, v)),
            );
        }
    }
    out
}

fn add_atom(state: &mut State, element: Element, (charge, valence): (i8, u8), span: Span) -> usize {
    let mut atom = AAtom::element(element, charge, valence, span);
    if charge != 0 {
        state.charge_left -= i32::from(charge);
        atom.budget_charge = true;
    }
    state.atoms.push(atom);
    state.atoms.len() - 1
}

/// Creates `n` atoms of `element` in every combination of valence and charge.
fn copies(
    state: &State,
    element: Element,
    n: u32,
    span: Span,
    ctx: &mut Ctx<'_>,
) -> Vec<(State, Vec<usize>)> {
    let mut partial = vec![(state.clone(), Vec::new())];
    for _ in 0..n {
        let mut next = Vec::new();
        for (s, atoms) in partial {
            for option in options(element, s.charge_left, ctx.single_atom) {
                if !ctx.spend(1) {
                    ctx.exhausted = true;
                    return Vec::new();
                }
                let mut s2 = s.clone();
                let index = add_atom(&mut s2, element, option, span);
                let mut atoms2 = atoms.clone();
                atoms2.push(index);
                next.push((s2, atoms2));
            }
        }
        if next.len() > COMBINATION_LIMIT {
            ctx.truncated = true;
            next.truncate(COMBINATION_LIMIT);
        }
        partial = next;
    }
    partial
}

fn atom(
    state: &State,
    unit: &Unit,
    element: Element,
    symbol: Span,
    ctx: &mut Ctx<'_>,
) -> Vec<(State, Piece)> {
    let n = unit.count();
    if element != Element::H && state.would_pass_limit(usize::try_from(n).unwrap_or(usize::MAX)) {
        ctx.oversized.get_or_insert(unit.span);
        return Vec::new();
    }
    if n == 1 {
        return copies(state, element, 1, symbol, ctx)
            .into_iter()
            .map(|(s, atoms)| {
                let head = atoms[0];
                (
                    s,
                    Piece::Link {
                        head,
                        push: atoms,
                        bond: None,
                    },
                )
            })
            .collect();
    }
    copies(state, element, n, unit.span, ctx)
        .into_iter()
        .map(|(s, atoms)| {
            (
                s,
                Piece::Pendants {
                    heads: atoms.clone(),
                    push: atoms,
                    bond: None,
                },
            )
        })
        .collect()
}

fn abbreviation(
    state: &State,
    unit: &Unit,
    entry: usize,
    reversed: bool,
    symbol: Span,
    ctx: &mut Ctx<'_>,
) -> Vec<(State, Piece)> {
    let definition = &ctx.dictionary.entries()[entry];
    let n = unit.count();
    let size = definition.fragment().atoms().len();
    if state.would_pass_limit(
        usize::try_from(n)
            .unwrap_or(usize::MAX)
            .saturating_mul(size),
    ) {
        ctx.oversized.get_or_insert(unit.span);
        return Vec::new();
    }
    let attachments = definition.attachments().len();
    if attachments == 0 {
        if !ctx.whole_allowed || n != 1 || !state.atoms.is_empty() {
            ctx.fail(Failure::new(
                Code::MoleculeAsGroup,
                "a complete molecule cannot be bonded to anything",
                symbol,
                format!("{} has no attachment point", definition.name()),
                symbol.start,
            ));
            return Vec::new();
        }
        let mut s = state.clone();
        instantiate(&mut s, ctx, entry, reversed, symbol);
        return vec![(s, Piece::Whole)];
    }
    let mut s = state.clone();
    let mut sites = Vec::new();
    for _ in 0..n {
        sites.push(instantiate(&mut s, ctx, entry, reversed, symbol));
    }
    if attachments == 1 {
        let heads: Vec<usize> = sites.iter().map(|site| site[0]).collect();
        let piece = if n == 1 {
            Piece::Link {
                head: heads[0],
                push: heads,
                bond: None,
            }
        } else {
            Piece::Pendants {
                heads: heads.clone(),
                push: heads,
                bond: None,
            }
        };
        return vec![(s, piece)];
    }
    let push: Vec<usize> = sites.iter().flatten().copied().collect();
    let mut states = vec![s];
    for pair in sites.windows(2) {
        let (tail, head) = (pair[0][1], pair[1][0]);
        states = states
            .into_iter()
            .flat_map(|st| {
                st.orders(tail, head, None)
                    .into_iter()
                    .take(1)
                    .map(move |k| {
                        let mut s2 = st.clone();
                        s2.link(tail, head, k, None);
                        s2
                    })
            })
            .collect();
    }
    states
        .into_iter()
        .map(|s| {
            (
                s,
                Piece::Link {
                    head: sites[0][0],
                    push: push.clone(),
                    bond: None,
                },
            )
        })
        .collect()
}

/// Copies a dictionary fragment into the state, returning its attachment atoms in bonding order.
fn instantiate(
    state: &mut State,
    ctx: &Ctx<'_>,
    entry: usize,
    reversed: bool,
    symbol: Span,
) -> Vec<usize> {
    let definition = &ctx.dictionary.entries()[entry];
    let fragment = definition.fragment();
    let offset = state.atoms.len();
    for (index, atom) in fragment.atoms().iter().enumerate() {
        let attach: u8 = definition
            .attachments()
            .iter()
            .filter(|a| a.atom == index)
            .map(|a| a.order.valence())
            .sum();
        let used = fragment.bond_valence(index) + atom.hydrogens + u8::from(atom.aromatic);
        state.atoms.push(AAtom {
            kind: atom.kind,
            valence: used + attach,
            used,
            charge: atom.charge,
            aromatic: atom.aromatic,
            hydrogens: atom.hydrogens,
            span: symbol,
            group: symbol,
            reversed_site: reversed && attach > 0,
            budget_charge: false,
            hydrogenated: false,
        });
    }
    for bond in fragment.bonds() {
        state.bonds.push(super::state::ABond {
            a: bond.a + offset,
            b: bond.b + offset,
            order: bond.order,
            token: Some(symbol),
        });
    }
    let mut sites: Vec<usize> = definition
        .attachments()
        .iter()
        .map(|a| a.atom + offset)
        .collect();
    if reversed {
        sites.reverse();
    }
    sites
}

/// Applies rule R1 to an assembled group interior.
fn repeat_group(
    state: &State,
    group: &State,
    unit: &Unit,
    bond: Option<(BondOrder, Span)>,
    ctx: &mut Ctx<'_>,
) -> Vec<(State, Piece)> {
    let free: u32 = group.atoms.iter().map(|a| u32::from(a.free())).sum();
    let Some(head) = group.atoms.iter().position(|a| a.free() > 0) else {
        ctx.fail(Failure::new(
            Code::GroupWithoutValence,
            "this group has no free valence to bond with",
            unit.span,
            "every atom inside is already saturated".into(),
            unit.span.start,
        ));
        return Vec::new();
    };
    let tail = group.host().unwrap_or(head);
    let n = unit.count();
    let copies = usize::try_from(n).unwrap_or(usize::MAX);
    let heavy = group
        .atoms
        .iter()
        .filter(|a| a.element_kind() != Some(Element::H))
        .count();
    if state.would_pass_limit(copies.saturating_mul(heavy)) {
        ctx.oversized.get_or_insert(unit.span);
        return Vec::new();
    }
    // Nested counts multiply, so every copied atom is charged.
    if !ctx.spend(copies.saturating_mul(group.atoms.len())) {
        ctx.exhausted = true;
        return Vec::new();
    }
    let mut out = Vec::new();
    if n == 1 || free == 1 || free >= 3 {
        let mut s = state.clone();
        let heads: Vec<usize> = (0..n).map(|_| embed(&mut s, group) + head).collect();
        if n > 1 && free >= 3 {
            s.forks.push((unit.span, "pendant"));
            s.score -= 1;
        }
        out.push((
            s,
            Piece::Pendants {
                heads,
                push: Vec::new(),
                bond,
            },
        ));
    }
    if n > 1 && free >= 2 {
        let mut s = state.clone();
        let offsets: Vec<usize> = (0..n).map(|_| embed(&mut s, group)).collect();
        if free >= 3 {
            s.forks.push((unit.span, "chain"));
        }
        let push: Vec<usize> = offsets
            .iter()
            .flat_map(|&offset| group.stack.iter().map(move |&i| i + offset))
            .collect();
        let width = group.atoms.len();
        let mut states = vec![s];
        for (i, pair) in offsets.windows(2).enumerate() {
            let (from, to) = (pair[0] + tail, pair[1] + head);
            let mut next = Vec::new();
            for st in states {
                for k in st.orders(from, to, None) {
                    if !ctx.spend(1) {
                        ctx.exhausted = true;
                        return out;
                    }
                    let mut s2 = st.clone();
                    s2.link(from, to, k, None);
                    let saturated = i == 0
                        || s2.atoms[pair[0]..pair[0] + width]
                            .iter()
                            .all(|a| a.free() == 0);
                    if saturated {
                        next.push(s2);
                    }
                }
            }
            if next.len() > COMBINATION_LIMIT {
                ctx.truncated = true;
                next.truncate(COMBINATION_LIMIT);
            }
            states = next;
        }
        for s in states {
            out.push((
                s,
                Piece::Link {
                    head: offsets[0] + head,
                    push: push.clone(),
                    bond,
                },
            ));
        }
    }
    out
}

/// Appends a group interior to the state and returns the index offset of its atoms.
fn embed(state: &mut State, group: &State) -> usize {
    let offset = state.atoms.len();
    state.atoms.extend(group.atoms.iter().cloned());
    state
        .bonds
        .extend(group.bonds.iter().map(|b| super::state::ABond {
            a: b.a + offset,
            b: b.b + offset,
            ..b.clone()
        }));
    state.score += group.score;
    state.forks.extend(group.forks.iter().copied());
    state.forward += group.forward;
    state.backward += group.backward;
    offset
}
