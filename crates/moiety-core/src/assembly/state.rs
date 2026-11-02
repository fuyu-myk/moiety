use std::collections::HashMap;

use crate::element::Element;
use crate::interpret::MAX_ATOMS;
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
    ///
    /// When that unit is a run of halogens, as the `F2` of `CF2H`, the hydrogen
    /// bonds to the atom they hang from instead.
    pub fn hydrogen_host(&self) -> Option<usize> {
        if let Some(atom) = self
            .previous
            .iter()
            .rev()
            .copied()
            .find(|&atom| self.atoms[atom].free() > 0)
        {
            return Some(atom);
        }
        let halogens = !self.previous.is_empty()
            && self.previous.iter().all(|&atom| {
                matches!(
                    self.atoms[atom].element_kind(),
                    Some(Element::F | Element::CL | Element::BR | Element::I)
                )
            });
        if !halogens {
            return None;
        }
        let host = self
            .bonds
            .iter()
            .find_map(|bond| bond.other(self.previous[0]))?;
        (self.atoms[host].free() > 0).then_some(host)
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

    /// Whether adding `atoms` atoms other than hydrogen would take the reading past
    /// [`MAX_ATOMS`]. Written hydrogens become implicit, so they do not count.
    pub fn would_pass_limit(&self, atoms: usize) -> bool {
        let heavy = self
            .atoms
            .iter()
            .filter(|a| a.element_kind() != Some(Element::H))
            .count();
        heavy.saturating_add(atoms) > MAX_ATOMS
    }

    /// Total free valence of atoms later units can still reach.
    pub fn live_free(&self) -> u32 {
        self.atoms.iter().map(|atom| u32::from(atom.free())).sum()
    }

    /// A hash of what makes two states interchangeable: see [`State::same_key`].
    fn key_hash(&self) -> u64 {
        let mut hash = mix(0, u64::from(self.started));
        hash = mix(
            hash,
            u64::from(self.charge_left.unsigned_abs()) << 1 | u64::from(self.charge_left < 0),
        );
        for atom in &self.atoms {
            let kind = match atom.kind {
                AtomKind::Element(element) => u64::from(element.atomic_number()),
                AtomKind::Attachment => 0,
            };
            let charge = u64::from(atom.charge.unsigned_abs()) << 1 | u64::from(atom.charge < 0);
            hash = mix(
                hash,
                kind << 24 | u64::from(atom.valence) << 16 | u64::from(atom.used) << 8 | charge,
            );
        }
        for bond in &self.bonds {
            hash = mix(hash, bond.a as u64);
            hash = mix(hash, bond.b as u64 ^ (bond.order as u64) << 56);
        }
        for &index in &self.stack {
            hash = mix(hash, index as u64);
        }
        for pending in &self.pending {
            hash = mix(hash, pending.heads.len() as u64);
            for &head in &pending.heads {
                hash = mix(hash, head as u64);
            }
        }
        hash
    }

    /// Whether the two states can only lead to the same molecules, so one may stand for both.
    fn same_key(&self, other: &State) -> bool {
        let atom = |a: &AAtom| (a.kind, a.valence, a.used, a.charge);
        let bond = |b: &ABond| (b.a, b.b, b.order);
        self.started == other.started
            && self.charge_left == other.charge_left
            && self.stack == other.stack
            && self.atoms.len() == other.atoms.len()
            && self.atoms.iter().map(atom).eq(other.atoms.iter().map(atom))
            && self.bonds.len() == other.bonds.len()
            && self.bonds.iter().map(bond).eq(other.bonds.iter().map(bond))
            && self.pending.len() == other.pending.len()
            && self
                .pending
                .iter()
                .zip(&other.pending)
                .all(|(a, b)| a.heads == b.heads)
    }
}

fn mix(hash: u64, value: u64) -> u64 {
    (hash.rotate_left(5) ^ value).wrapping_mul(0x517c_c1b7_2722_0a95)
}

/// Removes duplicate states, keeping the best score, then keeps the best-scored
/// states up to `limit` states and `atom_limit` atoms in all. Returns whether
/// states were dropped by the limits.
pub(crate) fn dedupe(states: Vec<State>, limit: usize, atom_limit: usize) -> (Vec<State>, bool) {
    let mut index: HashMap<u64, Vec<usize>> = HashMap::new();
    let mut kept: Vec<State> = Vec::new();
    for state in states {
        let slots = index.entry(state.key_hash()).or_default();
        match slots.iter().copied().find(|&i| kept[i].same_key(&state)) {
            Some(i) => {
                if state.score > kept[i].score {
                    kept[i] = state;
                }
            }
            None => {
                slots.push(kept.len());
                kept.push(state);
            }
        }
    }
    let atoms: usize = kept.iter().map(|s| s.atoms.len()).sum();
    if kept.len() <= limit && atoms <= atom_limit {
        return (kept, false);
    }
    kept.sort_by_key(|s| std::cmp::Reverse(s.score));
    let mut atoms = 0;
    let fits = kept
        .iter()
        .take(limit)
        .take_while(|s| {
            atoms += s.atoms.len();
            atoms <= atom_limit
        })
        .count();
    kept.truncate(fits);
    (kept, true)
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
