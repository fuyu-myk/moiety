//! The inverse writer: condensed formulas from molecular graphs.
//!
//! Rings and charged groups are written with dictionary spellings, everything
//! else atom by atom along the longest chain, with side chains in parentheses.

use std::collections::BTreeMap;

use moiety_core::{AtomKind, BondOrder, Dictionary, Element, Entry, Molecule};

/// A molecule the writer cannot express as a condensed formula.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum WriteError {
    /// A ring that no dictionary entry matches.
    #[error("atom {0} is in a ring that no dictionary entry matches")]
    UnknownRing(usize),
    /// An attachment point that is not at either end of a chain.
    #[error("attachment point {0} is not at the end of a chain")]
    InnerAttachment(usize),
    /// The molecule has no atoms.
    #[error("the molecule is empty")]
    Empty,
}

/// Writes `molecule` as a condensed formula that moiety reads back to the same molecule.
///
/// # Errors
///
/// Fails for rings the dictionary does not cover and for attachment points that are
/// not at the end of a chain.
pub fn condense(molecule: &Molecule, dictionary: &Dictionary) -> Result<String, WriteError> {
    if molecule.is_empty() {
        return Err(WriteError::Empty);
    }
    let groups = match_groups(molecule, dictionary)?;
    let mut parts: Vec<String> = Vec::new();
    for component in molecule.components() {
        parts.push(Writer::new(molecule, &groups, &component).write()?);
    }
    let mut merged: Vec<(usize, String)> = Vec::new();
    for part in parts {
        match merged.last_mut() {
            Some((count, last)) if *last == part => *count += 1,
            _ => merged.push((1, part)),
        }
    }
    Ok(merged
        .into_iter()
        .map(|(count, part)| {
            if count > 1 {
                format!("{count}{part}")
            } else {
                part
            }
        })
        .collect::<Vec<_>>()
        .join("."))
}

/// A dictionary fragment found in the molecule.
struct Group {
    atoms: Vec<usize>,
    spelling: String,
    reversed: Option<String>,
}

fn match_groups(molecule: &Molecule, dictionary: &Dictionary) -> Result<Vec<Group>, WriteError> {
    let adjacency = molecule.adjacency();
    let in_ring = ring_atoms(molecule, &adjacency);
    let mut entries: Vec<&Entry> = dictionary
        .entries()
        .iter()
        .filter(|e| e.attachments().len() == 1)
        .filter(|e| has_ring(e.fragment()) || e.fragment().atoms().iter().any(|a| a.charge != 0))
        .collect();
    entries.sort_by_key(|e| std::cmp::Reverse(e.fragment().atoms().len()));
    let mut covered = vec![false; molecule.atoms().len()];
    let mut groups = Vec::new();
    for entry in &entries {
        for anchor in 0..molecule.atoms().len() {
            if covered[anchor] {
                continue;
            }
            if let Some(atoms) = embed(entry, molecule, &adjacency, anchor, &covered) {
                for &a in &atoms {
                    covered[a] = true;
                }
                groups.push(Group {
                    atoms,
                    spelling: entry.spellings().first().cloned().unwrap_or_default(),
                    reversed: entry.reversed().first().cloned(),
                });
            }
        }
    }
    match (0..molecule.atoms().len()).find(|&a| in_ring[a] && !covered[a]) {
        Some(atom) => Err(WriteError::UnknownRing(atom)),
        None => Ok(groups),
    }
}

/// Maps the fragment of `entry` onto the molecule with its attachment atom at `anchor`.
fn embed(
    entry: &Entry,
    molecule: &Molecule,
    adjacency: &[Vec<(usize, usize)>],
    anchor: usize,
    covered: &[bool],
) -> Option<Vec<usize>> {
    let fragment = entry.fragment();
    let site = entry.attachments()[0];
    let inner = fragment.adjacency();
    let mut matcher = Matcher {
        fragment,
        inner: &inner,
        molecule,
        adjacency,
        covered,
        map: vec![None; fragment.atoms().len()],
        used: vec![false; molecule.atoms().len()],
    };
    if covered[anchor] || !matcher.atom_matches(site.atom, anchor, 1) {
        return None;
    }
    matcher.map[site.atom] = Some(anchor);
    matcher.used[anchor] = true;
    let order = bfs_order(&inner, site.atom);
    matcher.extend(1, &order).then(|| {
        matcher
            .map
            .into_iter()
            .map(|m| m.expect("every fragment atom is mapped"))
            .collect()
    })
}

fn bfs_order(adjacency: &[Vec<(usize, usize)>], start: usize) -> Vec<usize> {
    let mut order = vec![start];
    let mut seen = vec![false; adjacency.len()];
    seen[start] = true;
    let mut i = 0;
    while i < order.len() {
        for &(next, _) in &adjacency[order[i]] {
            if !seen[next] {
                seen[next] = true;
                order.push(next);
            }
        }
        i += 1;
    }
    order
}

/// Backtracking subgraph match of a dictionary fragment onto a molecule.
struct Matcher<'a> {
    fragment: &'a Molecule,
    inner: &'a [Vec<(usize, usize)>],
    molecule: &'a Molecule,
    adjacency: &'a [Vec<(usize, usize)>],
    covered: &'a [bool],
    map: Vec<Option<usize>>,
    used: Vec<bool>,
}

impl Matcher<'_> {
    fn extend(&mut self, depth: usize, order: &[usize]) -> bool {
        let Some(&atom) = order.get(depth) else {
            return true;
        };
        let (parent, bond) = self.inner[atom]
            .iter()
            .copied()
            .find(|&(p, _)| self.map[p].is_some())
            .expect("breadth-first order reaches every atom through a mapped neighbour");
        let target_parent = self.map[parent].expect("checked above");
        for &(candidate, molecule_bond) in &self.adjacency[target_parent] {
            if self.used[candidate]
                || self.covered[candidate]
                || self.molecule.bonds()[molecule_bond].order != self.fragment.bonds()[bond].order
                || !self.atom_matches(atom, candidate, 0)
                || !self.consistent(atom, candidate)
            {
                continue;
            }
            self.map[atom] = Some(candidate);
            self.used[candidate] = true;
            if self.extend(depth + 1, order) {
                return true;
            }
            self.map[atom] = None;
            self.used[candidate] = false;
        }
        false
    }

    fn consistent(&self, atom: usize, candidate: usize) -> bool {
        self.inner[atom]
            .iter()
            .all(|&(other, b)| match self.map[other] {
                Some(target) => self.adjacency[candidate].iter().any(|&(t, mb)| {
                    t == target && self.molecule.bonds()[mb].order == self.fragment.bonds()[b].order
                }),
                None => true,
            })
    }

    fn atom_matches(&self, atom: usize, candidate: usize, extra: usize) -> bool {
        let (f, m) = (
            self.fragment.atoms()[atom],
            self.molecule.atoms()[candidate],
        );
        f.kind == m.kind
            && f.charge == m.charge
            && f.hydrogens == m.hydrogens
            && f.aromatic == m.aromatic
            && self.adjacency[candidate].len() == self.inner[atom].len() + extra
    }
}

fn has_ring(molecule: &Molecule) -> bool {
    let components = molecule.components().len();
    molecule.bonds().len() + components > molecule.atoms().len()
}

/// Atoms on at least one cycle: those with a bond whose removal leaves its ends connected.
fn ring_atoms(molecule: &Molecule, adjacency: &[Vec<(usize, usize)>]) -> Vec<bool> {
    let mut in_ring = vec![false; molecule.atoms().len()];
    for (index, bond) in molecule.bonds().iter().enumerate() {
        let mut seen = vec![false; in_ring.len()];
        let mut stack = vec![bond.a];
        seen[bond.a] = true;
        while let Some(atom) = stack.pop() {
            for &(next, via) in &adjacency[atom] {
                if via != index && !seen[next] {
                    seen[next] = true;
                    stack.push(next);
                }
            }
        }
        if seen[bond.b] {
            in_ring[bond.a] = true;
            in_ring[bond.b] = true;
        }
    }
    in_ring
}

/// A node of the contracted tree: one atom, or one dictionary group.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Node {
    Atom(usize),
    Group(usize),
}

struct Writer<'a> {
    molecule: &'a Molecule,
    groups: &'a [Group],
    nodes: Vec<Node>,
    edges: Vec<Vec<(usize, BondOrder)>>,
    charge: i32,
}

impl<'a> Writer<'a> {
    fn new(molecule: &'a Molecule, groups: &'a [Group], component: &[usize]) -> Self {
        let mut node_of = BTreeMap::new();
        let mut nodes = Vec::new();
        for (g, group) in groups.iter().enumerate() {
            if component.contains(&group.atoms[0]) {
                for &atom in &group.atoms {
                    node_of.insert(atom, nodes.len());
                }
                nodes.push(Node::Group(g));
            }
        }
        for &atom in component {
            node_of.entry(atom).or_insert_with(|| {
                nodes.push(Node::Atom(atom));
                nodes.len() - 1
            });
        }
        let mut edges = vec![Vec::new(); nodes.len()];
        for bond in molecule.bonds() {
            let (Some(&a), Some(&b)) = (node_of.get(&bond.a), node_of.get(&bond.b)) else {
                continue;
            };
            if a != b {
                edges[a].push((b, bond.order));
                edges[b].push((a, bond.order));
            }
        }
        let charge = nodes
            .iter()
            .filter_map(|n| match n {
                Node::Atom(a) => Some(i32::from(molecule.atoms()[*a].charge)),
                Node::Group(_) => None,
            })
            .sum();
        Self {
            molecule,
            groups,
            nodes,
            edges,
            charge,
        }
    }

    fn attachments(&self) -> Vec<usize> {
        (0..self.nodes.len())
            .filter(|&n| {
                matches!(self.nodes[n], Node::Atom(a) if self.molecule.atoms()[a].kind == AtomKind::Attachment)
            })
            .collect()
    }

    fn write(&self) -> Result<String, WriteError> {
        let attachments = self.attachments();
        if let Some(&extra) = attachments.get(2) {
            let Node::Atom(atom) = self.nodes[extra] else {
                unreachable!("attachments are atoms")
            };
            return Err(WriteError::InnerAttachment(atom));
        }
        if self.nodes.len() == 1 {
            if let Node::Atom(a) = self.nodes[0] {
                return Ok(format!(
                    "{}{}",
                    lone_atom(self.molecule, a),
                    charge_suffix(self.charge)
                ));
            }
        }
        let path = self.main_chain();
        let mut out = String::new();
        let first = path[0];
        if let Node::Atom(a) = self.nodes[first] {
            if self.molecule.atoms()[a].kind == AtomKind::Attachment {
                out.push_str(bond_text(self.edges[first][0].1, true));
                self.chain(&path[1..], Some(first), &mut out)?;
                return Ok(self.finish(out));
            }
        }
        self.chain(&path, None, &mut out)?;
        Ok(self.finish(out))
    }

    fn finish(&self, mut out: String) -> String {
        if !out.ends_with('-') {
            out.push_str(&charge_suffix(self.charge));
        }
        out
    }

    fn is_inline(&self, node: usize) -> bool {
        match self.nodes[node] {
            Node::Atom(a) => {
                let atom = self.molecule.atoms()[a];
                let terminal = self.edges[node].len() == 1;
                let halogen = matches!(
                    atom.element(),
                    Some(Element::F | Element::CL | Element::BR | Element::I)
                );
                let oxo = matches!(atom.element(), Some(Element::O | Element::S))
                    && atom.hydrogens == 0
                    && self.edges[node]
                        .first()
                        .is_some_and(|e| e.1 == BondOrder::Double);
                terminal && (halogen || oxo)
            }
            Node::Group(_) => false,
        }
    }

    fn main_chain(&self) -> Vec<usize> {
        let attachments = self.attachments();
        if let [from, to] = attachments[..] {
            return self
                .path_between(from, to, None)
                .expect("a component is connected");
        }
        let skeleton_degree = |n: usize| {
            self.edges[n]
                .iter()
                .filter(|&&(o, _)| !self.is_inline(o))
                .count()
        };
        let terminal: Vec<usize> = (0..self.nodes.len())
            .filter(|&n| !self.is_inline(n) && self.edges[n].len() <= 1)
            .collect();
        let ends: Vec<usize> = if terminal.is_empty() {
            (0..self.nodes.len())
                .filter(|&n| !self.is_inline(n) && skeleton_degree(n) <= 1)
                .collect()
        } else {
            terminal
        };
        let starts = match attachments.first() {
            Some(&n) => vec![n],
            None if ends.is_empty() => vec![0],
            None => ends,
        };
        let mut best: Vec<usize> = Vec::new();
        for &start in &starts {
            let path = self.longest_from(start, None);
            if path.len() > best.len() {
                best = path;
            }
        }
        best
    }

    fn path_between(&self, from: usize, to: usize, parent: Option<usize>) -> Option<Vec<usize>> {
        if from == to {
            return Some(vec![to]);
        }
        self.edges[from]
            .iter()
            .filter(|(next, _)| Some(*next) != parent)
            .find_map(|&(next, _)| self.path_between(next, to, Some(from)))
            .map(|rest| [vec![from], rest].concat())
    }

    fn longest_from(&self, node: usize, parent: Option<usize>) -> Vec<usize> {
        let mut best: Vec<usize> = Vec::new();
        for &(next, _) in &self.edges[node] {
            if Some(next) == parent || self.is_inline(next) {
                continue;
            }
            let path = self.longest_from(next, Some(node));
            if path.len() > best.len() {
                best = path;
            }
        }
        let mut path = vec![node];
        path.extend(best);
        path
    }

    fn chain(
        &self,
        path: &[usize],
        parent: Option<usize>,
        out: &mut String,
    ) -> Result<(), WriteError> {
        for (i, &node) in path.iter().enumerate() {
            let previous = if i == 0 { parent } else { Some(path[i - 1]) };
            let next = path.get(i + 1).copied();
            out.push_str(&self.label(node, previous.is_none() && i == 0, next.is_none()));
            let mut inline: BTreeMap<String, usize> = BTreeMap::new();
            let mut branches: Vec<String> = Vec::new();
            for &(other, order) in &self.edges[node] {
                if Some(other) == previous || Some(other) == next {
                    continue;
                }
                if self.is_inline(other)
                    && (self.is_halogen(other) || !self.has_branches(node, previous, next))
                {
                    *inline.entry(self.label(other, false, true)).or_default() += 1;
                    continue;
                }
                let mut branch = String::from("(");
                branch.push_str(bond_text(order, false));
                let sub = self.longest_from(other, Some(node));
                self.chain(&sub, Some(node), &mut branch)?;
                branch.push(')');
                branches.push(branch);
            }
            for (label, count) in inline {
                out.push_str(&label);
                if count > 1 {
                    out.push_str(&count.to_string());
                }
            }
            push_merged(out, branches);
            if let Some(next) = next {
                let order = self.edges[node]
                    .iter()
                    .find(|(other, _)| *other == next)
                    .map(|(_, order)| *order)
                    .expect("consecutive chain nodes are bonded");
                let terminal = path.len() == i + 2 && self.edges[next].len() == 1;
                if !(terminal && self.is_heteroatom(next)) {
                    out.push_str(bond_text(order, false));
                }
            }
        }
        Ok(())
    }

    /// Whether `node` has side chains that must be written in parentheses.
    fn has_branches(&self, node: usize, previous: Option<usize>, next: Option<usize>) -> bool {
        self.edges[node].iter().any(|&(other, _)| {
            Some(other) != previous && Some(other) != next && !self.is_inline(other)
        })
    }

    fn is_halogen(&self, node: usize) -> bool {
        matches!(self.nodes[node], Node::Atom(a) if matches!(self.molecule.atoms()[a].element(), Some(Element::F | Element::CL | Element::BR | Element::I)))
    }

    fn is_heteroatom(&self, node: usize) -> bool {
        matches!(self.nodes[node], Node::Atom(a) if matches!(self.molecule.atoms()[a].element(), Some(Element::O | Element::S | Element::N)))
    }

    fn label(&self, node: usize, first: bool, last: bool) -> String {
        match self.nodes[node] {
            Node::Group(g) => {
                let group = &self.groups[g];
                match (&group.reversed, first && !last) {
                    (Some(reversed), true) => reversed.clone(),
                    _ => group.spelling.clone(),
                }
            }
            Node::Atom(a) => {
                let atom = self.molecule.atoms()[a];
                let Some(element) = atom.element() else {
                    return "-".to_string();
                };
                let hydrogens = match atom.hydrogens {
                    0 => String::new(),
                    1 => "H".to_string(),
                    n => format!("H{n}"),
                };
                if first && element != Element::C && !hydrogens.is_empty() {
                    format!("{hydrogens}{}", element.symbol())
                } else {
                    format!("{}{hydrogens}", element.symbol())
                }
            }
        }
    }
}

fn push_merged(out: &mut String, branches: Vec<String>) {
    let mut iter = branches.into_iter().peekable();
    while let Some(branch) = iter.next() {
        let mut count = 1;
        while iter.peek() == Some(&branch) {
            iter.next();
            count += 1;
        }
        out.push_str(&branch);
        if count > 1 {
            out.push_str(&count.to_string());
        }
    }
}

fn bond_text(order: BondOrder, marker: bool) -> &'static str {
    match order {
        BondOrder::Double => "=",
        BondOrder::Triple => "#",
        _ if marker => "-",
        _ => "",
    }
}

fn lone_atom(molecule: &Molecule, atom: usize) -> String {
    let a = molecule.atoms()[atom];
    let Some(element) = a.element() else {
        return "*".to_string();
    };
    let hydrogens = match a.hydrogens {
        0 => return element.symbol().to_string(),
        1 => "H".to_string(),
        n => format!("H{n}"),
    };
    let hydrogen_first = matches!(element.symbol(), "O" | "S" | "Se" | "F" | "Cl" | "Br" | "I");
    if hydrogen_first {
        format!("{hydrogens}{}", element.symbol())
    } else {
        format!("{}{hydrogens}", element.symbol())
    }
}

fn charge_suffix(charge: i32) -> String {
    match charge {
        0 => String::new(),
        1 => "+".to_string(),
        -1 => "-".to_string(),
        c if c > 0 => format!("^{c}+"),
        c => format!("^{}-", -c),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{default_dictionary, interpret};

    fn write(smiles: &str) -> String {
        let molecule = moiety_core::smiles::parse(smiles).unwrap();
        condense(&molecule, default_dictionary()).unwrap()
    }

    fn round_trips(smiles: &str) {
        let molecule = moiety_core::smiles::parse(smiles).unwrap();
        let written = condense(&molecule, default_dictionary()).unwrap();
        let top = interpret(&written).top().map(|r| r.smiles.clone());
        assert_eq!(
            top,
            Some(molecule.to_smiles()),
            "{smiles} written as {written}"
        );
    }

    #[test]
    fn writes_natural_condensed_formulas() {
        assert_eq!(write("CCO"), "CH3CH2OH");
        assert_eq!(write("CC(=O)O"), "CH3COOH");
        assert_eq!(write("CC(C)C"), "CH3CH(CH3)CH3");
        assert_eq!(write("OCc1ccccc1"), "BnOH");
        assert_eq!(write("C=CC"), "CH2=CHCH3");
        assert_eq!(write("O.O.O"), "3H2O");
        assert_eq!(write("[Na+].[Cl-]"), "Na+.Cl-");
    }

    #[test]
    fn round_trips_common_molecules() {
        for smiles in [
            "CCCCCC",
            "CC(C)(C)O",
            "OC(=O)CCC(=O)O",
            "CC(=O)OCC",
            "C[N+](=O)[O-]",
            "Cc1ccc(S(=O)(=O)Cl)cc1",
            "FC(F)(F)C(=O)O",
            "CN(C)C",
            "C#CC",
            "*CC(=O)O",
            "OCC(O)CO",
            "CC(C)CC(C)(C)C",
            "O=CCOC(C=O)=O",
        ] {
            round_trips(smiles);
        }
    }

    #[test]
    fn rejects_unknown_rings_and_inner_attachments() {
        let ring = moiety_core::smiles::parse("C1CCCC1").unwrap();
        assert!(matches!(
            condense(&ring, default_dictionary()),
            Err(WriteError::UnknownRing(_))
        ));
        assert_eq!(write("*CC(*)C"), "-CH2CH(CH3)-");
        let inner = moiety_core::smiles::parse("*CC(*)C*").unwrap();
        assert!(matches!(
            condense(&inner, default_dictionary()),
            Err(WriteError::InnerAttachment(_))
        ));
    }
}
