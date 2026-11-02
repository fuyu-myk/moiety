//! Valence resolution: turns a parse tree into every molecule it can denote.

mod state;
mod units;

use state::{AAtom, Pending, State, dedupe};
use units::{Piece, expand};

use crate::diagnostic::Code;
use crate::dictionary::Dictionary;
use crate::element::Element;
use crate::molecule::{Atom, AtomKind, Bond, BondOrder, Molecule};
use crate::parser::{Item, Seq, Tail, Tree, UnitKind};
use crate::reading::Underdetermination;
use crate::span::Span;

const STATE_LIMIT: usize = 20_000;
const BOND_COMBINATION_LIMIT: usize = 256;

/// A molecule derived from one parse tree, with what is needed to rank and explain it.
#[derive(Clone, Debug)]
pub(crate) struct Derivation {
    pub molecule: Molecule,
    pub atom_spans: Vec<Span>,
    pub bond_spans: Vec<Span>,
    pub score: i32,
    pub forks: Vec<(Span, &'static str)>,
    pub forward: u32,
    pub backward: u32,
    pub charge_sites: Vec<Span>,
}

/// Why a branch of the search died, located in normalized coordinates.
#[derive(Clone, Debug)]
pub(crate) struct Failure {
    pub code: Code,
    pub message: String,
    pub primary: (Span, String),
    pub secondary: Vec<(Span, String)>,
    pub progress: usize,
}

impl Failure {
    pub fn new(code: Code, message: &str, span: Span, label: String, progress: usize) -> Self {
        Self {
            code,
            message: message.to_string(),
            primary: (span, label),
            secondary: Vec::new(),
            progress,
        }
    }

    fn with(mut self, span: Span, label: String) -> Self {
        self.secondary.push((span, label));
        self
    }
}

/// Everything one parse tree produced.
#[derive(Debug, Default)]
pub(crate) struct Assembly {
    pub derivations: Vec<Derivation>,
    pub failure: Option<Failure>,
    pub underdetermined: Option<(Span, Underdetermination)>,
    pub truncated: bool,
}

pub(crate) struct Ctx<'a> {
    pub dictionary: &'a Dictionary,
    pub failure: Option<Failure>,
    pub underdetermined: Option<(Span, Underdetermination)>,
    pub truncated: bool,
    pub single_atom: bool,
    pub whole_allowed: bool,
    pub component_end: usize,
}

impl<'a> Ctx<'a> {
    fn new(dictionary: &'a Dictionary, tree: &Tree) -> Self {
        let lone = tree.seq.items.len() == 1 && tree.lead.is_none();
        let single_atom = lone
            && matches!(tree.seq.items[0].unit.kind, UnitKind::Atom { .. })
            && tree.seq.items[0].unit.count() == 1
            && !matches!(tree.tail, Tail::Marker(..));
        Self {
            dictionary,
            failure: None,
            underdetermined: None,
            truncated: false,
            single_atom,
            whole_allowed: lone && tree.tail == Tail::None,
            component_end: tree.span.end,
        }
    }

    pub fn fail(&mut self, failure: Failure) {
        if self
            .failure
            .as_ref()
            .is_none_or(|current| failure.progress > current.progress)
        {
            self.failure = Some(failure);
        }
    }

    fn finish(self, derivations: Vec<Derivation>) -> Assembly {
        Assembly {
            derivations,
            failure: self.failure,
            underdetermined: self.underdetermined,
            truncated: self.truncated,
        }
    }
}

/// Assembles one parse tree under rules R1 to R4.
pub(crate) fn assemble(tree: &Tree, dictionary: &Dictionary) -> Assembly {
    let mut ctx = Ctx::new(dictionary, tree);
    let mut initial = State::default();
    if let Tail::Charge(charge, _) = tree.tail {
        initial.charge_left = charge;
    }
    if let Some((order, span)) = tree.lead {
        initial.atoms.push(attachment(order, span));
        initial.stack.push(0);
        initial.started = true;
    }
    let derivations = run(tree, &tree.seq, &mut ctx, initial)
        .into_iter()
        .map(derive)
        .collect();
    ctx.finish(derivations)
}

/// Reads a component as a salt under rule R5: leading metals become cations and the rest an anion.
pub(crate) fn assemble_salt(tree: &Tree, dictionary: &Dictionary) -> Assembly {
    let metal = |item: &Item| {
        matches!(item.unit.kind, UnitKind::Atom { element, .. } if element.is_metal())
    };
    let split = tree
        .seq
        .items
        .iter()
        .position(|item| !metal(item) || item.bond.is_some());
    let eligible = tree.lead.is_none() && tree.tail == Tail::None;
    let Some(split) = split.filter(|&s| s > 0 && eligible) else {
        return Assembly::default();
    };
    let rest = Seq {
        items: tree.seq.items[split..].to_vec(),
    };
    if rest.items.iter().any(metal) {
        return Assembly::default();
    }
    let cations: Vec<(Element, u32, Span)> = tree.seq.items[..split]
        .iter()
        .filter_map(|item| match item.unit.kind {
            UnitKind::Atom { element, .. } => Some((element, item.unit.count(), item.unit.span)),
            _ => None,
        })
        .collect();
    let mut combos: Vec<Vec<u8>> = vec![Vec::new()];
    for (element, _, _) in &cations {
        combos = combos
            .into_iter()
            .flat_map(|combo| {
                element.valences().iter().filter(|&&v| v >= 1).map(move |&v| {
                    let mut next = combo.clone();
                    next.push(v);
                    next
                })
            })
            .collect();
    }
    let rest_span = Span::new(rest.items[0].unit.span.start, tree.span.end);
    let mut result = Assembly::default();
    for combo in combos {
        let total: i32 = cations
            .iter()
            .zip(&combo)
            .map(|((_, count, _), &v)| i32::try_from(*count).unwrap_or(0) * i32::from(v))
            .sum();
        let anion = Tree {
            seq: rest.clone(),
            tail: Tail::Charge(-total, rest_span),
            ..tree.clone()
        };
        let mut ctx = Ctx::new(dictionary, &anion);
        let mut initial = State {
            charge_left: -total,
            score: -1,
            ..State::default()
        };
        for ((element, count, span), &v) in cations.iter().zip(&combo) {
            for _ in 0..*count {
                let charge = i8::try_from(v).unwrap_or(0);
                initial.atoms.push(AAtom::element(*element, charge, 0, *span));
            }
        }
        result
            .derivations
            .extend(run(&anion, &anion.seq, &mut ctx, initial).into_iter().map(derive));
        result.truncated |= ctx.truncated;
    }
    result
}

fn run(tree: &Tree, seq: &Seq, ctx: &mut Ctx<'_>, initial: State) -> Vec<State> {
    let allowance = match tree.tail {
        Tail::Marker(order, _) => u32::from(order.valence()),
        _ => 0,
    };
    let mut resolved = Vec::new();
    for state in assemble_seq(seq, ctx, initial, Some(allowance)) {
        resolved.extend(finish_pending(state, ctx));
    }
    resolved
        .into_iter()
        .filter_map(|s| complete(s, tree, ctx))
        .collect()
}

fn attachment(order: BondOrder, span: Span) -> AAtom {
    AAtom {
        kind: AtomKind::Attachment,
        valence: order.valence(),
        used: 0,
        charge: 0,
        aromatic: false,
        hydrogens: 0,
        span,
        group: span,
        reversed_site: false,
        budget_charge: false,
        hydrogenated: false,
    }
}

/// Assembles a unit sequence. `allowance` is `Some` at the top level, where it
/// bounds the open valence later units could still absorb, and `None` inside groups.
pub(crate) fn assemble_seq(
    seq: &Seq,
    ctx: &mut Ctx<'_>,
    initial: State,
    allowance: Option<u32>,
) -> Vec<State> {
    let future = future_atoms(seq);
    let mut states = vec![initial];
    for (i, item) in seq.items.iter().enumerate() {
        if let Some(found) = carbon_count(seq, i, allowance.is_none()) {
            ctx.underdetermined.get_or_insert(found);
            return Vec::new();
        }
        let mut next = Vec::new();
        for state in &states {
            for (expanded, piece) in expand(state, &item.unit, ctx) {
                next.extend(attach(expanded, piece, item, ctx));
            }
        }
        let unreachable = |s: &State| match allowance {
            Some(extra) => s.started && s.live_free() > 3 * future[i + 1] + extra,
            None => false,
        };
        if let Some(dead) = next.iter().find(|s| s.has_dead_atom() || unreachable(s)) {
            let progress = if unreachable(dead) { ctx.component_end } else { item.unit.span.end };
            ctx.fail(open_valence(dead, progress));
        }
        next.retain(|s| !s.has_dead_atom() && !unreachable(s));
        let (kept, truncated) = dedupe(next, STATE_LIMIT);
        ctx.truncated |= truncated;
        if kept.is_empty() {
            return kept;
        }
        states = kept;
    }
    states
}

fn future_atoms(seq: &Seq) -> Vec<u32> {
    let mut sums = vec![0; seq.items.len() + 1];
    for (i, item) in seq.items.iter().enumerate().rev() {
        let per_copy = match &item.unit.kind {
            UnitKind::Atom { .. } | UnitKind::Abbreviation { .. } => 1,
            UnitKind::Group { seq, .. } => future_atoms(seq)[0],
        };
        sums[i] = sums[i + 1] + per_copy * item.unit.count();
    }
    sums
}

fn carbon_count(seq: &Seq, i: usize, in_group: bool) -> Option<(Span, Underdetermination)> {
    let unit = &seq.items[i].unit;
    let n = unit.count();
    if !matches!(unit.kind, UnitKind::Atom { element: Element::C, .. }) || n < 2 {
        return None;
    }
    let hydrogens = seq.items.get(i + 1).and_then(|item| match item.unit.kind {
        UnitKind::Atom {
            element: Element::H,
            ..
        } if item.bond.is_none() => Some((item.unit.count(), item.unit.span.end)),
        _ => None,
    });
    let (m, end, after) = hydrogens.map_or((0, unit.span.end, i + 1), |(m, end)| (m, end, i + 2));
    let neighbours = u32::from(i > 0) + u32::from(after < seq.items.len()) + u32::from(in_group);
    let whole = !in_group && i == 0 && after == seq.items.len();
    let unsaturation = i64::from(2 * n + 2) - i64::from(m) - i64::from(neighbours);
    let reason = if !whole && unsaturation > 0 {
        Underdetermination::UnknownRing
    } else {
        Underdetermination::MolecularFormula
    };
    Some((Span::new(unit.span.start, end), reason))
}

fn attach(mut state: State, piece: Piece, item: &Item, ctx: &mut Ctx<'_>) -> Vec<State> {
    let unit = &item.unit;
    match piece {
        Piece::Whole => {
            state.started = true;
            vec![state]
        }
        Piece::Pendants { heads, push, bond } => {
            let bond = bond.or(item.bond);
            if !state.started {
                state.pending.push(Pending {
                    heads,
                    push,
                    bond,
                    span: unit.span,
                });
                return vec![state];
            }
            let Some(host) = state.host() else {
                ctx.fail(no_host(&state, unit.span, ctx.component_end));
                return Vec::new();
            };
            let room = state.atoms[host].free();
            if usize::from(room) < heads.len() {
                ctx.fail(
                    Failure::new(
                        Code::HostFull,
                        "too many groups for one atom",
                        unit.span,
                        format!("these need {} bonds", heads.len()),
                        unit.span.start,
                    )
                    .with(state.atoms[host].group, format!("this atom has room for {room}")),
                );
                return Vec::new();
            }
            let bonded = bond_all(state, host, &heads, bond);
            if bonded.is_empty() {
                ctx.fail(bond_failure(unit.span, bond.is_some()));
            }
            bonded
                .into_iter()
                .map(|mut s| {
                    s.forward += count(heads.len());
                    s.stack.extend(&push);
                    s
                })
                .collect()
        }
        Piece::Link { head, push, bond } => {
            let bond = bond.or(item.bond);
            if !state.started {
                state.started = true;
                let pending = std::mem::take(&mut state.pending);
                if pending.is_empty() {
                    state.stack.extend(push);
                    return vec![state];
                }
                let heads: Vec<usize> = pending.iter().flat_map(|p| p.heads.iter().copied()).collect();
                let room = state.atoms[head].free();
                if usize::from(room) < heads.len() {
                    let mut failure = Failure::new(
                        Code::NoAnchor,
                        "the leading groups have no anchor",
                        unit.span,
                        format!("this atom has room for {room}"),
                        unit.span.start,
                    );
                    for p in &pending {
                        failure = failure.with(p.span, "this group needs an atom to hang from".into());
                    }
                    ctx.fail(failure);
                    return Vec::new();
                }
                let explicit = pending.iter().find_map(|p| p.bond).or(bond);
                let bonded = bond_all(state, head, &heads, explicit);
                if bonded.is_empty() {
                    ctx.fail(bond_failure(unit.span, explicit.is_some()));
                }
                return bonded
                    .into_iter()
                    .map(|mut s| {
                        s.backward += count(heads.len());
                        for p in &pending {
                            s.stack.extend(&p.push);
                        }
                        s.stack.extend(&push);
                        s
                    })
                    .collect();
            }
            let Some(host) = state.host() else {
                ctx.fail(no_host(&state, unit.span, ctx.component_end));
                return Vec::new();
            };
            let orders = state.orders(host, head, bond.map(|b| b.0));
            if orders.is_empty() {
                ctx.fail(
                    bond_failure(unit.span, bond.is_some())
                        .with(state.atoms[host].group, "the atom it would bond to".into()),
                );
                return Vec::new();
            }
            orders
                .into_iter()
                .map(|k| {
                    let mut s = state.clone();
                    s.link(host, head, k, bond.map(|b| b.1));
                    absorb_hydrogen(&mut s, host, head);
                    if s.atoms[host].reversed_site {
                        s.backward += 1;
                    } else {
                        s.forward += 1;
                    }
                    s.stack.extend(&push);
                    s
                })
                .collect()
        }
    }
}

fn count(n: usize) -> u32 {
    u32::try_from(n).unwrap_or(u32::MAX)
}

fn bond_all(state: State, center: usize, heads: &[usize], bond: Option<(BondOrder, Span)>) -> Vec<State> {
    let mut partial = vec![state];
    for &head in heads {
        let mut next = Vec::new();
        for s in partial {
            for k in s.orders(center, head, bond.map(|b| b.0)) {
                let mut s2 = s.clone();
                s2.link(center, head, k, bond.map(|b| b.1));
                absorb_hydrogen(&mut s2, center, head);
                next.push(s2);
            }
        }
        next.truncate(BOND_COMBINATION_LIMIT);
        partial = next;
    }
    partial
}

/// Widens a heavy atom's group span over a hydrogen written right next to it.
fn absorb_hydrogen(state: &mut State, a: usize, b: usize) {
    let is_h = |atom: &AAtom| atom.element_kind() == Some(Element::H);
    let (heavy, hydrogen) = match (is_h(&state.atoms[a]), is_h(&state.atoms[b])) {
        (false, true) => (a, b),
        (true, false) => (b, a),
        _ => return,
    };
    let h = state.atoms[hydrogen].span;
    state.atoms[heavy].hydrogenated = true;
    let group = &mut state.atoms[heavy].group;
    if h.start == group.end || h.end == group.start {
        *group = group.cover(h);
    }
}

fn open_valence(state: &State, progress: usize) -> Failure {
    let reachable = |i: usize| {
        state.stack.contains(&i) || state.pending.iter().any(|p| p.heads.contains(&i))
    };
    let open = (0..state.atoms.len())
        .filter(|&i| state.atoms[i].free() > 0)
        .min_by_key(|&i| reachable(i))
        .map(|i| &state.atoms[i])
        .expect("called only for states with open valence");
    let free = open.free();
    let plural = if free == 1 { "" } else { "s" };
    Failure::new(
        Code::OpenValence,
        "valence left open at the end of the formula",
        open.group,
        format!("this atom still has {free} open valence{plural}"),
        progress,
    )
}

fn no_host(state: &State, unit: Span, end: usize) -> Failure {
    let heavy = |&&atom: &&usize| state.atoms[atom].element_kind() != Some(Element::H);
    let terminal = state
        .stack
        .iter()
        .rev()
        .find(heavy)
        .or(state.stack.last())
        .map_or(unit, |&atom| state.atoms[atom].group);
    Failure::new(
        Code::NoHost,
        "chain cannot continue past a terminal group",
        terminal,
        "this group has no free valence left".into(),
        unit.start,
    )
    .with(Span::new(unit.start, end), "but more units follow it".into())
}

fn bond_failure(unit: Span, explicit: bool) -> Failure {
    let label = if explicit {
        "this bond is stronger than the atoms can take"
    } else {
        "this unit cannot bond to the atom before it"
    };
    Failure::new(
        Code::BondNotPossible,
        "these atoms cannot form this bond",
        unit,
        label.into(),
        unit.start,
    )
}

/// Resolves leading pendants that never met an anchor: the first copy anchors the rest.
pub(crate) fn finish_pending(mut state: State, ctx: &mut Ctx<'_>) -> Vec<State> {
    if state.started || state.pending.is_empty() {
        return vec![state];
    }
    let pending = std::mem::take(&mut state.pending);
    state.started = true;
    let anchor = pending[0].heads[0];
    let heads: Vec<usize> = pending
        .iter()
        .flat_map(|p| p.heads.iter().copied())
        .filter(|&h| h != anchor)
        .collect();
    let explicit = pending.iter().find_map(|p| p.bond);
    let room = state.atoms[anchor].free();
    let bonded = bond_all(state, anchor, &heads, explicit);
    if bonded.is_empty() {
        let mut failure = Failure::new(
            Code::NoAnchor,
            "the leading groups have no anchor",
            pending[0].span,
            format!("the first of these has room for {room}"),
            ctx.component_end,
        );
        for p in &pending[1..] {
            failure = failure.with(p.span, "this group needs an atom to hang from".into());
        }
        ctx.fail(failure);
    }
    bonded
        .into_iter()
        .map(|mut s| {
            s.forward += count(heads.len());
            for p in &pending {
                s.stack.extend(&p.push);
            }
            if !s.stack.contains(&anchor) {
                s.stack.push(anchor);
            }
            s
        })
        .collect()
}

/// Applies the trailing marker and checks rule R4.
fn complete(mut state: State, tree: &Tree, ctx: &mut Ctx<'_>) -> Option<State> {
    if let Tail::Marker(order, span) = tree.tail {
        let star = state.atoms.len();
        state.atoms.push(attachment(order, span));
        let Some(host) = state.host() else {
            ctx.fail(no_host(&state, span, ctx.component_end));
            return None;
        };
        if state.orders(host, star, Some(order)).is_empty() {
            ctx.fail(bond_failure(span, true));
            return None;
        }
        state.bond(host, star, order.valence(), Some(span));
        state.forward += 1;
    }
    if state.atoms.iter().any(|a| a.free() > 0) {
        ctx.fail(open_valence(&state, ctx.component_end));
        return None;
    }
    if state.charge_left != 0 {
        let span = match tree.tail {
            Tail::Charge(_, span) => span,
            _ => tree.span,
        };
        ctx.fail(Failure::new(
            Code::ChargeNotPlaced,
            "the written charge cannot be placed",
            span,
            "no atom can carry this charge".into(),
            ctx.component_end,
        ));
        return None;
    }
    Some(state)
}

/// Builds the final molecule, folding hydrogens into implicit counts.
fn derive(state: State) -> Derivation {
    let mut neighbours = vec![Vec::new(); state.atoms.len()];
    for (index, bond) in state.bonds.iter().enumerate() {
        neighbours[bond.a].push((bond.b, index));
        neighbours[bond.b].push((bond.a, index));
    }
    let collapsible = |i: usize| {
        let atom = &state.atoms[i];
        atom.element_kind() == Some(Element::H)
            && atom.charge == 0
            && neighbours[i].len() == 1
            && state.atoms[neighbours[i][0].0]
                .element_kind()
                .is_some_and(|e| e != Element::H)
    };
    let mut index = vec![None; state.atoms.len()];
    let mut atoms = Vec::new();
    let mut atom_spans = Vec::new();
    for (i, atom) in state.atoms.iter().enumerate() {
        if collapsible(i) {
            continue;
        }
        index[i] = Some(atoms.len());
        atoms.push(Atom {
            kind: atom.kind,
            charge: atom.charge,
            hydrogens: atom.hydrogens,
            aromatic: atom.aromatic,
        });
        atom_spans.push(atom.group);
    }
    for i in 0..state.atoms.len() {
        if collapsible(i) {
            let heavy = index[neighbours[i][0].0].expect("heavy atoms are kept");
            atoms[heavy].hydrogens += 1;
        }
    }
    let mut bonds = Vec::new();
    let mut bond_spans = Vec::new();
    for bond in &state.bonds {
        if let (Some(a), Some(b)) = (index[bond.a], index[bond.b]) {
            bonds.push(Bond { a, b, order: bond.order });
            bond_spans.push(
                bond.token
                    .unwrap_or_else(|| state.atoms[bond.a].group.cover(state.atoms[bond.b].group)),
            );
        }
    }
    let charge_sites = state
        .atoms
        .iter()
        .filter(|a| a.budget_charge)
        .map(|a| a.span)
        .collect();
    Derivation {
        molecule: Molecule::from_parts(atoms, bonds),
        atom_spans,
        bond_spans,
        score: state.score,
        forks: state.forks,
        forward: state.forward,
        backward: state.backward,
        charge_sites,
    }
}
