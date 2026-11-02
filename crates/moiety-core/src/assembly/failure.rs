use super::state::State;
use crate::diagnostic::Code;
use crate::element::Element;
use crate::span::Span;

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

    pub fn with(mut self, span: Span, label: String) -> Self {
        self.secondary.push((span, label));
        self
    }
}

pub(crate) fn open_valence(state: &State, progress: usize) -> Failure {
    let reachable =
        |i: usize| state.stack.contains(&i) || state.pending.iter().any(|p| p.heads.contains(&i));
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

pub(crate) fn no_host(state: &State, unit: Span, end: usize) -> Failure {
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
    .with(
        Span::new(unit.start, end),
        "but more units follow it".into(),
    )
}

pub(crate) fn bond_failure(unit: Span, explicit: Option<Span>) -> Failure {
    match explicit {
        Some(token) => Failure::new(
            Code::BondNotPossible,
            "these atoms cannot form this bond",
            token,
            "this bond is stronger than the atoms can take".into(),
            unit.start,
        )
        .with(unit, "the unit it bonds".into()),
        None => Failure::new(
            Code::BondNotPossible,
            "these atoms cannot form this bond",
            unit,
            "this unit cannot bond to the atom before it".into(),
            unit.start,
        ),
    }
}
