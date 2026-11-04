//! Valence resolution: turns a parse tree into every molecule it can denote.

mod derive;
mod failure;
mod state;
mod units;

use std::cell::Cell;

use derive::{derive, pentavalent_nitrogen_without_oxo};
pub(crate) use failure::Failure;
use failure::{bond_failure, no_host, open_valence};
use state::{AAtom, Pending, State, dedupe};
use units::{Piece, expand};

use crate::diagnostic::Code;
use crate::dictionary::Dictionary;
use crate::element::Element;
use crate::molecule::{AtomKind, BondOrder, Molecule};
use crate::parser::{Item, Seq, Tail, Tree, UnitKind};
use crate::reading::Underdetermination;
use crate::span::Span;

const STATE_LIMIT: usize = 20_000;
/// Most atoms, summed over its states, that a frontier keeps. Every state copies
/// its molecule, so this bounds memory where the state count alone does not.
const STATE_ATOMS: usize = 400_000;
const BOND_COMBINATION_LIMIT: usize = 256;

/// States that cannot complete, discarded as soon as they are made, allowed per search step.
const DISCARDS_PER_STEP: usize = 10;

/// What one component may spend on its search, or a suggestion search on all of its edits.
pub(crate) struct Budget {
    steps: Cell<usize>,
    discards: Cell<usize>,
}

impl Budget {
    /// `steps` search steps, and [`DISCARDS_PER_STEP`] discarded states for each.
    pub fn new(steps: usize) -> Self {
        Self {
            steps: Cell::new(steps),
            discards: Cell::new(steps.saturating_mul(DISCARDS_PER_STEP)),
        }
    }

    /// Whether every search step is spent.
    pub fn is_spent(&self) -> bool {
        self.steps.get() == 0
    }
}

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

/// Everything one parse tree produced.
#[derive(Debug, Default)]
pub(crate) struct Assembly {
    pub derivations: Vec<Derivation>,
    pub failure: Option<Failure>,
    pub underdetermined: Option<(Span, Underdetermination)>,
    pub truncated: bool,
    pub exhausted: bool,
    /// Where a reading would pass [`MAX_ATOMS`](crate::interpret::MAX_ATOMS).
    pub oversized: Option<Span>,
}

pub(crate) struct Ctx<'a> {
    pub dictionary: &'a Dictionary,
    pub budget: &'a Budget,
    pub exhausted: bool,
    pub failure: Option<Failure>,
    pub underdetermined: Option<(Span, Underdetermination)>,
    pub truncated: bool,
    pub oversized: Option<Span>,
    pub single_atom: bool,
    pub whole_allowed: bool,
    /// Whether the sequence being assembled is a group that starts with a bond,
    /// whose first unit bonds to the group's host by its head.
    pub anchored: bool,
    pub component_end: usize,
}

impl<'a> Ctx<'a> {
    fn new(dictionary: &'a Dictionary, budget: &'a Budget, tree: &Tree) -> Self {
        let lone = tree.seq.items.len() == 1 && tree.lead.is_none();
        let single_atom = lone
            && matches!(tree.seq.items[0].unit.kind, UnitKind::Atom { .. })
            && tree.seq.items[0].unit.count() == 1
            && !matches!(tree.tail, Tail::Marker(..));
        Self {
            dictionary,
            budget,
            exhausted: false,
            failure: None,
            underdetermined: None,
            truncated: false,
            oversized: None,
            single_atom,
            whole_allowed: lone && tree.tail == Tail::None,
            anchored: false,
            component_end: tree.span.end,
        }
    }

    /// Charges `n` search steps to the budget. Returns `false` once the budget is spent.
    pub fn spend(&self, n: usize) -> bool {
        let left = self.budget.steps.get().saturating_sub(n);
        self.budget.steps.set(left);
        left > 0
    }

    /// Charges one state that is discarded as soon as it is made. Returns `false`
    /// once the allowance for discarded states is spent.
    pub fn discard(&self) -> bool {
        let left = self.budget.discards.get().saturating_sub(1);
        self.budget.discards.set(left);
        left > 0
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
            exhausted: self.exhausted || self.budget.is_spent(),
            oversized: self.oversized,
        }
    }
}

/// Assembles one parse tree under rules R1 to R4, spending from `budget`.
pub(crate) fn assemble(tree: &Tree, dictionary: &Dictionary, budget: &Budget) -> Assembly {
    let mut ctx = Ctx::new(dictionary, budget, tree);
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
pub(crate) fn assemble_salt(tree: &Tree, dictionary: &Dictionary, budget: &Budget) -> Assembly {
    let metal = |item: &Item| matches!(item.unit.kind, UnitKind::Atom { element, .. } if element.is_metal());
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
    let choices: Vec<Vec<u8>> = cations
        .iter()
        .map(|(element, _, _)| {
            element
                .valences()
                .iter()
                .copied()
                .filter(|&v| v >= 1)
                .collect()
        })
        .collect();
    if choices.iter().any(Vec::is_empty) {
        return Assembly::default();
    }
    let rest_span = Span::new(rest.items[0].unit.span.start, tree.span.end);
    let mut result = Assembly::default();
    let mut odometer = vec![0usize; choices.len()];
    loop {
        let combo: Vec<u8> = odometer.iter().zip(&choices).map(|(&i, c)| c[i]).collect();
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
        let mut ctx = Ctx::new(dictionary, budget, &anion);
        if !ctx.spend(1) {
            result.exhausted = true;
            return result;
        }
        let mut initial = State {
            charge_left: -total,
            score: -1,
            ..State::default()
        };
        for ((element, count, span), &v) in cations.iter().zip(&combo) {
            for _ in 0..*count {
                let charge = i8::try_from(v).unwrap_or(0);
                initial
                    .atoms
                    .push(AAtom::element(*element, charge, 0, *span));
            }
        }
        result.derivations.extend(
            run(&anion, &anion.seq, &mut ctx, initial)
                .into_iter()
                .map(derive),
        );
        result.truncated |= ctx.truncated;
        result.exhausted |= ctx.exhausted || budget.is_spent();
        result.oversized = result.oversized.or(ctx.oversized);
        if result.exhausted || !advance(&mut odometer, &choices) {
            return result;
        }
    }
}

/// Steps a mixed-radix counter. Returns `false` after the last combination.
fn advance(odometer: &mut [usize], choices: &[Vec<u8>]) -> bool {
    for (digit, options) in odometer.iter_mut().zip(choices).rev() {
        *digit += 1;
        if *digit < options.len() {
            return true;
        }
        *digit = 0;
    }
    false
}

fn run(tree: &Tree, seq: &Seq, ctx: &mut Ctx<'_>, initial: State) -> Vec<State> {
    let allowance = match tree.tail {
        Tail::Marker(order, _) => u32::from(order.valence()),
        Tail::Charge(charge, _) => charge.unsigned_abs(),
        Tail::None => 0,
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
        numbered: None,
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
        let unreachable = |s: &State| match allowance {
            Some(extra) => {
                s.started && s.live_free() > future[i + 1].saturating_mul(3).saturating_add(extra)
            }
            None => false,
        };
        let mut next = Vec::new();
        let mut atoms = 0;
        let mut dead_seen = false;
        for state in &states {
            for (expanded, piece) in expand(state, &item.unit, ctx) {
                let before = expanded.stack.len();
                for mut attached in attach(expanded, piece, item, ctx) {
                    if attached.has_dead_atom() || unreachable(&attached) {
                        if !ctx.discard() {
                            ctx.exhausted = true;
                            return Vec::new();
                        }
                        if !dead_seen {
                            dead_seen = true;
                            let at_end = future[i + 1] == 0;
                            let progress = if unreachable(&attached) && at_end {
                                ctx.component_end
                            } else {
                                item.unit.span.end
                            };
                            ctx.fail(open_valence(&attached, progress));
                        }
                        continue;
                    }
                    if !ctx.spend(1) {
                        ctx.exhausted = true;
                        return Vec::new();
                    }
                    attached.previous = attached.stack[before..].to_vec();
                    atoms += attached.atoms.len();
                    next.push(attached);
                    if atoms >= 2 * STATE_ATOMS {
                        let taken = std::mem::take(&mut next);
                        let (kept, truncated) = dedupe(taken, STATE_LIMIT, STATE_ATOMS);
                        ctx.truncated |= truncated;
                        atoms = kept.iter().map(|s| s.atoms.len()).sum();
                        next = kept;
                    }
                }
            }
        }
        if next.is_empty() && !ctx.spend(1) {
            ctx.exhausted = true;
            return Vec::new();
        }
        let (kept, truncated) = dedupe(next, STATE_LIMIT, STATE_ATOMS);
        ctx.truncated |= truncated;
        if kept.is_empty() {
            return kept;
        }
        states = kept;
    }
    states
}

fn future_atoms(seq: &Seq) -> Vec<u32> {
    let mut sums = vec![0u32; seq.items.len() + 1];
    for (i, item) in seq.items.iter().enumerate().rev() {
        let per_copy = match &item.unit.kind {
            UnitKind::Atom { .. } | UnitKind::Abbreviation { .. } => 1,
            UnitKind::Group { seq, .. } => future_atoms(seq)[0],
            UnitKind::Arene(arene) => arene.substituents.iter().fold(
                u32::try_from(arene.size).unwrap_or(u32::MAX),
                |sum: u32, s| {
                    let copies = u32::try_from(s.positions.len()).unwrap_or(u32::MAX);
                    sum.saturating_add(future_atoms(&s.seq)[0].saturating_mul(copies))
                },
            ),
        };
        sums[i] = sums[i + 1].saturating_add(per_copy.saturating_mul(item.unit.count()));
    }
    sums
}

fn carbon_count(seq: &Seq, i: usize, in_group: bool) -> Option<(Span, Underdetermination)> {
    let unit = &seq.items[i].unit;
    let n = unit.count();
    if !matches!(
        unit.kind,
        UnitKind::Atom {
            element: Element::C,
            ..
        }
    ) || n < 2
    {
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
    let hydrogen = matches!(
        unit.kind,
        UnitKind::Atom {
            element: Element::H,
            ..
        }
    );
    let host_of = |state: &State| {
        if hydrogen {
            state.hydrogen_host()
        } else {
            state.host()
        }
    };
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
            let Some(host) = host_of(&state) else {
                ctx.fail(no_host(&state, unit.span, ctx.component_end));
                return Vec::new();
            };
            let hosts = state
                .numbered_hosts(host, heads.len())
                .unwrap_or_else(|| vec![host; heads.len()]);
            let room = state.atoms[host].free();
            if hosts[0] == hosts[hosts.len() - 1] && usize::from(room) < heads.len() {
                ctx.fail(
                    Failure::new(
                        Code::HostFull,
                        "too many groups for one atom",
                        unit.span,
                        format!("these need {} bonds", heads.len()),
                        unit.span.start,
                    )
                    .with(
                        state.atoms[host].group,
                        format!("this atom has room for {room}"),
                    ),
                );
                return Vec::new();
            }
            let bonded = bond_pairs(state, &hosts, &heads, bond, ctx);
            if bonded.is_empty() {
                ctx.fail(bond_failure(unit.span, bond.map(|b| b.1)));
            }
            let events = if hydrogen { 0 } else { count(heads.len()) };
            bonded
                .into_iter()
                .map(|mut s| {
                    s.forward += events;
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
                let heads: Vec<usize> = pending
                    .iter()
                    .flat_map(|p| p.heads.iter().copied())
                    .collect();
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
                        failure =
                            failure.with(p.span, "this group needs an atom to hang from".into());
                    }
                    ctx.fail(failure);
                    return Vec::new();
                }
                let explicit = pending.iter().find_map(|p| p.bond).or(bond);
                let events = count(
                    heads
                        .iter()
                        .filter(|&&h| state.atoms[h].element_kind() != Some(Element::H))
                        .count(),
                );
                let bonded = bond_all(state, head, &heads, explicit, ctx);
                if bonded.is_empty() {
                    ctx.fail(bond_failure(unit.span, explicit.map(|b| b.1)));
                }
                return bonded
                    .into_iter()
                    .map(|mut s| {
                        s.backward += events;
                        for p in &pending {
                            s.stack.extend(&p.push);
                        }
                        s.stack.extend(&push);
                        s
                    })
                    .collect();
            }
            let Some(host) = host_of(&state) else {
                ctx.fail(no_host(&state, unit.span, ctx.component_end));
                return Vec::new();
            };
            let orders = state.orders(host, head, bond.map(|b| b.0));
            if orders.is_empty() {
                ctx.fail(
                    bond_failure(unit.span, bond.map(|b| b.1))
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
                    if !hydrogen {
                        if s.atoms[host].reversed_site {
                            s.backward += 1;
                        } else {
                            s.forward += 1;
                        }
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

fn bond_all(
    state: State,
    center: usize,
    heads: &[usize],
    bond: Option<(BondOrder, Span)>,
    ctx: &mut Ctx<'_>,
) -> Vec<State> {
    bond_pairs(state, &vec![center; heads.len()], heads, bond, ctx)
}

/// Bonds each head to the center at the same index, in every combination of orders.
fn bond_pairs(
    state: State,
    centers: &[usize],
    heads: &[usize],
    bond: Option<(BondOrder, Span)>,
    ctx: &mut Ctx<'_>,
) -> Vec<State> {
    let mut partial = vec![state];
    for (&center, &head) in centers.iter().zip(heads) {
        let mut next = Vec::new();
        for s in partial {
            for k in s.orders(center, head, bond.map(|b| b.0)) {
                if !ctx.spend(1) {
                    ctx.exhausted = true;
                    return Vec::new();
                }
                let mut s2 = s.clone();
                s2.link(center, head, k, bond.map(|b| b.1));
                absorb_hydrogen(&mut s2, center, head);
                next.push(s2);
            }
        }
        if next.len() > BOND_COMBINATION_LIMIT {
            ctx.truncated = true;
            next.truncate(BOND_COMBINATION_LIMIT);
        }
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
    let bonded = bond_all(state, anchor, &heads, explicit, ctx);
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

/// Places charge still unplaced at the end on the atoms left open, when charging
/// each of them closes its valence exactly. This reaches atoms inside groups,
/// which are assembled before the written charge is known to apply to them.
fn place_remaining_charge(state: &mut State) {
    if state.charge_left == 0 {
        return;
    }
    let sign: i8 = if state.charge_left > 0 { 1 } else { -1 };
    let open: Vec<usize> = (0..state.atoms.len())
        .filter(|&i| state.atoms[i].free() > 0)
        .collect();
    let fits = |i: usize| {
        let atom = &state.atoms[i];
        atom.charge == 0
            && atom.free() == 1
            && atom
                .element_kind()
                .is_some_and(|e| e != Element::H && e.charged_valences(sign).contains(&atom.used))
    };
    if open.len() != state.charge_left.unsigned_abs() as usize || !open.iter().all(|&i| fits(i)) {
        return;
    }
    for i in open {
        let atom = &mut state.atoms[i];
        atom.charge = sign;
        atom.valence = atom.used;
        atom.budget_charge = true;
    }
    state.charge_left = 0;
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
            ctx.fail(bond_failure(span, Some(span)));
            return None;
        }
        state.bond(host, star, order.valence(), Some(span));
        state.forward += 1;
    }
    place_remaining_charge(&mut state);
    if state.atoms.iter().any(|a| a.free() > 0) {
        ctx.fail(open_valence(&state, ctx.component_end));
        return None;
    }
    if let Some(nitrogen) = pentavalent_nitrogen_without_oxo(&state) {
        ctx.fail(Failure::new(
            Code::BondNotPossible,
            "nitrogen takes five bonds only with a double bond to oxygen",
            state.atoms[nitrogen].group,
            "this nitrogen would have five bonds".into(),
            ctx.component_end,
        ));
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
