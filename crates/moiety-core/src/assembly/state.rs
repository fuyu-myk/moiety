use std::collections::HashMap;

use crate::element::Element;
use crate::molecule::{AtomKind, BondOrder};
use crate::span::Span;

/// An atom under assembly.
#[derive(Clone, Debug)]
pub(crate) struct AAtom {
    pub kind: AtomKind,
    pub valence: u8,
    pub used: u8,
    pub charge: i8,
    pub aromatic: bool,
    pub hydrogens: u8,
    pub span: Span,
    pub group: Span,
    pub reversed_site: bool,
    pub budget_charge: bool,
    pub hydrogenated: bool,
}

impl AAtom {
    pub fn element(kind: Element, charge: i8, valence: u8, span: Span) -> Self {
        Self {
            kind: AtomKind::Element(kind),
            valence,
            used: 0,
            charge,
            aromatic: false,
            hydrogens: 0,
            span,
            group: span,
            reversed_site: false,
            budget_charge: false,
            hydrogenated: false,
        }
    }

    pub fn free(&self) -> u8 {
        self.valence.saturating_sub(self.used)
    }

    pub fn element_kind(&self) -> Option<Element> {
        match self.kind {
            AtomKind::Element(element) => Some(element),
            AtomKind::Attachment => None,
        }
    }
}

/// A bond under assembly.
#[derive(Clone, Debug)]
pub(crate) struct ABond {
    pub a: usize,
    pub b: usize,
    pub order: BondOrder,
    pub token: Option<Span>,
}

impl ABond {
    pub fn other(&self, atom: usize) -> Option<usize> {
        if self.a == atom {
            Some(self.b)
        } else if self.b == atom {
            Some(self.a)
        } else {
            None
        }
    }
}

/// Leading pendants waiting for an anchor (rule R2).
#[derive(Clone, Debug)]
pub(crate) struct Pending {
    pub heads: Vec<usize>,
    pub push: Vec<usize>,
    pub bond: Option<(BondOrder, Span)>,
    pub span: Span,
}

/// One branch of the assembly search.
#[derive(Clone, Debug, Default)]
pub(crate) struct State {
    pub atoms: Vec<AAtom>,
    pub bonds: Vec<ABond>,
    pub stack: Vec<usize>,
    pub previous: Vec<usize>,
    pub pending: Vec<Pending>,
    pub started: bool,
    pub charge_left: i32,
    pub score: i32,
    pub forks: Vec<(Span, &'static str)>,
    pub forward: u32,
    pub backward: u32,
}

impl State {
    /// The most recently created atom on the host stack that still has free valence.
    ///
    /// The search steps back over saturated atoms such as the `=O` of `COOH`, but
    /// stops at a saturated atom with written hydrogens: a `CH2` or `NH` written
    /// mid-chain is a chain link, so nothing after it bonds to an atom before it.
    pub fn host(&self) -> Option<usize> {
        for &index in self.stack.iter().rev() {
            let atom = &self.atoms[index];
            if atom.free() > 0 {
                return Some(index);
            }
            if atom.hydrogenated {
                return None;
            }
        }
        None
    }

    /// The host for a hydrogen: the most recent open atom of the unit written just before it.
    pub fn hydrogen_host(&self) -> Option<usize> {
        self.previous
            .iter()
            .rev()
            .copied()
            .find(|&atom| self.atoms[atom].free() > 0)
    }

    /// Bond orders that atoms `a` and `b` could form now.
    pub fn orders(&self, a: usize, b: usize, explicit: Option<BondOrder>) -> Vec<u8> {
        let (x, y) = (&self.atoms[a], &self.atoms[b]);
        if !can_bond(x, y) {
            return Vec::new();
        }
        let mut max = x.free().min(y.free()).min(3);
        if !multiple_allowed(x) || !multiple_allowed(y) {
            max = max.min(1);
        }
        for (center, partner) in [(x, y), (y, x)] {
            if raised(center) {
                let limit = if partner.element_kind() == Some(Element::O) {
                    2
                } else {
                    1
                };
                max = max.min(limit);
            }
        }
        match explicit {
            Some(order) => {
                let k = order.valence();
                if k <= max { vec![k] } else { Vec::new() }
            }
            None => (1..=max).collect(),
        }
    }

    /// Adds a bond between two units, applying the heteroatom-link penalty of ranking policy 1.
    pub fn link(&mut self, a: usize, b: usize, order: u8, token: Option<Span>) {
        if order == 1 && hetero_link(&self.atoms[a]) && hetero_link(&self.atoms[b]) {
            self.score -= 1;
        }
        self.bond(a, b, order, token);
    }

    /// Adds a bond with no ranking effect.
    pub fn bond(&mut self, a: usize, b: usize, order: u8, token: Option<Span>) {
        self.atoms[a].used += order;
        self.atoms[b].used += order;
        self.bonds.push(ABond {
            a,
            b,
            order: BondOrder::from_count(order).expect("orders are 1 to 3"),
            token,
        });
    }

    /// Whether an atom with free valence can no longer be reached by any later unit.
    pub fn has_dead_atom(&self) -> bool {
        self.atoms.iter().enumerate().any(|(index, atom)| {
            atom.free() > 0
                && !self.stack.contains(&index)
                && !self
                    .pending
                    .iter()
                    .any(|p| p.heads.contains(&index) || p.push.contains(&index))
        })
    }

    /// Total free valence of atoms later units can still reach.
    pub fn live_free(&self) -> u32 {
        self.atoms.iter().map(|atom| u32::from(atom.free())).sum()
    }

    fn key(&self) -> StateKey {
        StateKey {
            atoms: self
                .atoms
                .iter()
                .map(|a| (a.kind, a.valence, a.used, a.charge))
                .collect(),
            bonds: self.bonds.iter().map(|b| (b.a, b.b, b.order)).collect(),
            stack: self.stack.clone(),
            pending: self.pending.iter().map(|p| p.heads.clone()).collect(),
            started: self.started,
            charge_left: self.charge_left,
        }
    }
}

#[derive(PartialEq, Eq, Hash)]
struct StateKey {
    atoms: Vec<(AtomKind, u8, u8, i8)>,
    bonds: Vec<(usize, usize, BondOrder)>,
    stack: Vec<usize>,
    pending: Vec<Vec<usize>>,
    started: bool,
    charge_left: i32,
}

/// Removes duplicate states, keeping the best score, then keeps at most `limit` of them.
/// Returns whether states were dropped by the limit.
pub(crate) fn dedupe(states: Vec<State>, limit: usize) -> (Vec<State>, bool) {
    let mut index: HashMap<StateKey, usize> = HashMap::new();
    let mut kept: Vec<State> = Vec::new();
    for state in states {
        match index.entry(state.key()) {
            std::collections::hash_map::Entry::Occupied(slot) => {
                let i = *slot.get();
                if state.score > kept[i].score {
                    kept[i] = state;
                }
            }
            std::collections::hash_map::Entry::Vacant(slot) => {
                slot.insert(kept.len());
                kept.push(state);
            }
        }
    }
    let truncated = kept.len() > limit;
    if truncated {
        kept.sort_by_key(|s| std::cmp::Reverse(s.score));
        kept.truncate(limit);
    }
    (kept, truncated)
}

fn can_bond(x: &AAtom, y: &AAtom) -> bool {
    match (x.element_kind(), y.element_kind()) {
        (Some(a), Some(b)) => !(a.is_metal() && b.is_metal()) && carbon_ok(a, b) && carbon_ok(b, a),
        _ => true,
    }
}

fn carbon_ok(metal: Element, other: Element) -> bool {
    metal.bonds_to_carbon() || !matches!(other, Element::C | Element::H)
}

fn multiple_allowed(atom: &AAtom) -> bool {
    match atom.element_kind() {
        Some(element) => element.allows_multiple_bonds() && !atom.aromatic,
        None => true,
    }
}

/// Whether a neutral atom uses more than its lowest valence, as sulfur does in a sulfone.
fn raised(atom: &AAtom) -> bool {
    atom.charge == 0
        && atom
            .element_kind()
            .is_some_and(|element| atom.valence > element.valences()[0])
}

fn hetero_link(atom: &AAtom) -> bool {
    matches!(
        atom.element_kind(),
        Some(Element::N | Element::O | Element::S | Element::SE)
    )
}
