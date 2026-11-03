//! Rings the dictionary's ring templates number, such as pyridine, written with
//! locants (grammar, section 6.12): locants first as `4-N(CH3)2C5H4N` or
//! `2-C5H4NCH2OH` when the ring starts the formula, a trailing attachment as
//! `C5H4N-2` at the end of a chain, and a bare locant as `HOCH2-2,6-C5H3N-CH2OH`
//! inside one.

use std::collections::BTreeMap;

use moiety_core::{BondOrder, Dictionary, Molecule, RingTemplate};

use super::{Node, WriteError, Writer, push_label};

/// Most numberings kept for one ring, more than any bundled template's symmetry gives.
const NUMBERINGS: usize = 24;

/// A ring that a ring template numbers.
pub(super) struct TemplateRing<'a> {
    pub template: &'a RingTemplate,
    /// Every numbering the ring's symmetry allows: for each parent atom, the
    /// molecule atom it is.
    pub numberings: Vec<Vec<usize>>,
}

/// Finds the rings that the dictionary's templates number among atoms nothing
/// covers yet, and marks them covered. Each numbered atom holds its hydrogen or
/// one single bond to a substituent in no ring with it.
pub(super) fn template_rings<'a>(
    molecule: &Molecule,
    adjacency: &[Vec<(usize, usize)>],
    dictionary: &'a Dictionary,
    covered: &mut [bool],
) -> Vec<TemplateRing<'a>> {
    let mut found = Vec::new();
    for template in dictionary.rings() {
        let mut matcher = Matcher {
            template,
            inner: template.parent().adjacency(),
            molecule,
            adjacency,
            covered,
            map: Vec::new(),
            used: vec![false; molecule.atoms().len()],
            found: Vec::new(),
            before: 0,
        };
        let order = matcher.order();
        for start in 0..molecule.atoms().len() {
            matcher.map = vec![None; template.parent().atoms().len()];
            matcher.before = matcher.found.len();
            if matcher.fits(order[0], start) {
                matcher.map[order[0]] = Some(start);
                matcher.used[start] = true;
                matcher.extend(1, &order);
                matcher.used[start] = false;
            }
        }
        let mut groups: BTreeMap<Vec<usize>, Vec<Vec<usize>>> = BTreeMap::new();
        for numbering in matcher.found {
            let mut atoms = numbering.clone();
            atoms.sort_unstable();
            let group = groups.entry(atoms).or_default();
            if group.len() < NUMBERINGS {
                group.push(numbering);
            }
        }
        for (atoms, numberings) in groups {
            if atoms.iter().any(|&a| covered[a])
                || !substituents_are_bridges(molecule, adjacency, &atoms)
            {
                continue;
            }
            for &atom in &atoms {
                covered[atom] = true;
            }
            found.push(TemplateRing {
                template,
                numberings,
            });
        }
    }
    found
}

/// Whether every bond from the ring `atoms` to an atom outside it is single and
/// in no ring, so each substituent hangs from one position.
fn substituents_are_bridges(
    molecule: &Molecule,
    adjacency: &[Vec<(usize, usize)>],
    atoms: &[usize],
) -> bool {
    atoms.iter().all(|&atom| {
        adjacency[atom]
            .iter()
            .filter(|(next, _)| !atoms.contains(next))
            .all(|&(next, bond)| {
                molecule.bonds()[bond].order == BondOrder::Single
                    && !connected_without(adjacency, atom, next, bond)
            })
    })
}

/// Whether `a` reaches `b` without crossing `bond`.
fn connected_without(adjacency: &[Vec<(usize, usize)>], a: usize, b: usize, bond: usize) -> bool {
    let mut seen = vec![false; adjacency.len()];
    let mut stack = vec![a];
    seen[a] = true;
    while let Some(atom) = stack.pop() {
        for &(next, via) in &adjacency[atom] {
            if via != bond && !seen[next] {
                seen[next] = true;
                stack.push(next);
            }
        }
    }
    seen[b]
}

/// Backtracking match of a template's parent ring onto the molecule, collecting
/// every complete numbering.
struct Matcher<'a> {
    template: &'a RingTemplate,
    inner: Vec<Vec<(usize, usize)>>,
    molecule: &'a Molecule,
    adjacency: &'a [Vec<(usize, usize)>],
    covered: &'a [bool],
    map: Vec<Option<usize>>,
    used: Vec<bool>,
    found: Vec<Vec<usize>>,
    /// How many numberings were found before the current start atom.
    before: usize,
}

impl Matcher<'_> {
    /// The parent atoms in breadth-first order, so each after the first has a
    /// neighbour placed before it.
    fn order(&self) -> Vec<usize> {
        let mut order = vec![0];
        let mut seen = vec![false; self.inner.len()];
        seen[0] = true;
        let mut i = 0;
        while i < order.len() {
            for &(next, _) in &self.inner[order[i]] {
                if !seen[next] {
                    seen[next] = true;
                    order.push(next);
                }
            }
            i += 1;
        }
        order
    }

    /// Whether molecule atom `candidate` can be parent atom `atom`.
    fn fits(&self, atom: usize, candidate: usize) -> bool {
        let parent = self.template.parent().atoms()[atom];
        let found = self.molecule.atoms()[candidate];
        if self.covered[candidate]
            || parent.kind != found.kind
            || parent.aromatic != found.aromatic
            || parent.charge != found.charge
        {
            return false;
        }
        let degree = self.adjacency[candidate].len();
        let ring_degree = self.inner[atom].len();
        let substitutable = self
            .template
            .position(atom)
            .is_some_and(|p| self.template.substitutable(p));
        if substitutable {
            (found.hydrogens == 1 && degree == ring_degree)
                || (found.hydrogens == 0 && degree == ring_degree + 1)
        } else {
            found.hydrogens == parent.hydrogens && degree == ring_degree
        }
    }

    fn extend(&mut self, depth: usize, order: &[usize]) {
        if self.found.len() - self.before >= NUMBERINGS {
            return;
        }
        let Some(&atom) = order.get(depth) else {
            self.found.push(
                self.map
                    .iter()
                    .map(|m| m.expect("every atom is placed"))
                    .collect(),
            );
            return;
        };
        let parent_bonds = self.template.parent().bonds();
        let (placed, bond) = self.inner[atom]
            .iter()
            .copied()
            .find(|&(p, _)| self.map[p].is_some())
            .expect("breadth-first order reaches every atom through a placed neighbour");
        let anchor = self.map[placed].expect("checked above");
        for &(candidate, molecule_bond) in &self.adjacency[anchor] {
            let consistent = self.inner[atom]
                .iter()
                .all(|&(other, b)| match self.map[other] {
                    Some(target) => self.adjacency[candidate].iter().any(|&(t, mb)| {
                        t == target && self.molecule.bonds()[mb].order == parent_bonds[b].order
                    }),
                    None => true,
                });
            if self.used[candidate]
                || self.molecule.bonds()[molecule_bond].order != parent_bonds[bond].order
                || !self.fits(atom, candidate)
                || !consistent
            {
                continue;
            }
            self.map[atom] = Some(candidate);
            self.used[candidate] = true;
            self.extend(depth + 1, order);
            self.map[atom] = None;
            self.used[candidate] = false;
        }
    }
}

impl Writer<'_> {
    fn template_ring(&self, node: usize) -> &TemplateRing<'_> {
        let Node::Template(t) = self.nodes[node] else {
            unreachable!("only template nodes are numbered by templates")
        };
        &self.templates[t]
    }

    /// The position of the ring atom bonded to `neighbour` under `numbering`.
    fn template_site(&self, node: usize, numbering: &[usize], neighbour: usize) -> u8 {
        let template = self.template_ring(node).template;
        self.molecule
            .bonds()
            .iter()
            .find_map(|bond| {
                [(bond.a, bond.b), (bond.b, bond.a)]
                    .into_iter()
                    .find_map(|(inside, outside)| {
                        let index = numbering.iter().position(|&a| a == inside)?;
                        (self.node_of.get(&outside) == Some(&neighbour))
                            .then(|| template.position(index))
                            .flatten()
                    })
            })
            .expect("a ring node is bonded to each of its neighbours at a position")
    }

    /// The ring's formula with the hydrogens it has left, as `C5H3N`. A ring with
    /// none left has no formula the reader numbers.
    fn template_formula(&self, node: usize) -> Result<String, WriteError> {
        let template = self.template_ring(node).template;
        match template.hydrogens() - self.edges[node].len() {
            0 => Err(self.template_error(node)),
            left => Ok(template.written(left)),
        }
    }

    /// Writes the ring as a bracketed group on the unit before it, every other
    /// neighbour a substituent with its locant, as `(5-Br-2-C5H3N)`.
    fn template_bracketed(
        &self,
        node: usize,
        previous: usize,
        others: &[usize],
        out: &mut String,
    ) -> Result<(), WriteError> {
        let mut neighbours = others.to_vec();
        neighbours.push(previous);
        let positions = self
            .template_numbering(node, &neighbours, |_| true)
            .ok_or_else(|| self.template_error(node))?;
        let attachment = positions[others.len()];
        let formula = self.template_formula(node)?;
        if others.is_empty() {
            push_label(out, &format!("{formula}-{attachment}"));
            return Ok(());
        }
        let entries: Vec<(u8, usize)> = positions
            .iter()
            .copied()
            .zip(others.iter().copied())
            .collect();
        let prefix = self.template_prefix(node, &entries)?;
        out.push_str(&format!("({prefix}-{attachment}-{formula})"));
        Ok(())
    }

    fn template_error(&self, node: usize) -> WriteError {
        WriteError::UnknownRing(self.template_ring(node).numberings[0][0])
    }

    /// The positions of `neighbours` under the numbering that gives the lowest
    /// locants among those `allowed` accepts.
    fn template_numbering(
        &self,
        node: usize,
        neighbours: &[usize],
        allowed: impl Fn(&[u8]) -> bool,
    ) -> Option<Vec<u8>> {
        self.template_ring(node)
            .numberings
            .iter()
            .map(|numbering| {
                neighbours
                    .iter()
                    .map(|&n| self.template_site(node, numbering, n))
                    .collect::<Vec<u8>>()
            })
            .filter(|positions| allowed(positions))
            .min_by_key(|positions| {
                let mut sorted = positions.clone();
                sorted.sort_unstable();
                sorted
            })
    }

    /// The substituents with their locants, as `2-Cl-4-F`, joined for a prefix.
    fn template_prefix(&self, node: usize, entries: &[(u8, usize)]) -> Result<String, WriteError> {
        Ok(self
            .chunks(node, entries)?
            .iter()
            .map(|(locants, text)| format!("{locants}-{text}"))
            .collect::<Vec<_>>()
            .join("-"))
    }

    /// Writes a template ring that starts the formula: its substituents with
    /// their locants and, when `next` follows, the locant of its attachment, as
    /// `4-N(CH3)2C5H4N` or `2-C5H4NCH2OH`.
    pub(super) fn template_first(
        &self,
        node: usize,
        next: Option<usize>,
        out: &mut String,
    ) -> Result<(), WriteError> {
        let mut neighbours: Vec<usize> = self.edges[node]
            .iter()
            .map(|&(other, _)| other)
            .filter(|&o| Some(o) != next)
            .collect();
        let others = neighbours.len();
        neighbours.extend(next);
        if neighbours.is_empty() {
            return Err(self.template_error(node));
        }
        let positions = self
            .template_numbering(node, &neighbours, |_| true)
            .ok_or_else(|| self.template_error(node))?;
        let entries: Vec<(u8, usize)> = positions[..others]
            .iter()
            .copied()
            .zip(neighbours[..others].iter().copied())
            .collect();
        let prefix = self.template_prefix(node, &entries)?;
        let formula = self.template_formula(node)?;
        let text = match (next.is_some(), prefix.is_empty()) {
            (true, true) => format!("{}-{formula}", positions[others]),
            (true, false) => format!("{prefix}-{}-{formula}", positions[others]),
            (false, _) => format!("{prefix}{formula}"),
        };
        push_label(out, &text);
        Ok(())
    }

    /// Writes a template ring with a unit before it. At the end of a chain it
    /// takes a trailing attachment, as `C5H4N-2`, or brackets its substituents
    /// with locants first; inside a chain it takes a bare locant, as the
    /// `-2,6-C5H3N` of `HOCH2-2,6-C5H3N-CH2OH`, with the units of `rest` after it.
    pub(super) fn template_after(
        &self,
        node: usize,
        previous: usize,
        rest: &[usize],
        out: &mut String,
    ) -> Result<(), WriteError> {
        let next = rest.first().copied();
        let others: Vec<usize> = self.edges[node]
            .iter()
            .map(|&(other, _)| other)
            .filter(|&o| o != previous && Some(o) != next)
            .collect();
        let Some(next) = next else {
            return self.template_bracketed(node, previous, &others, out);
        };
        let formula = self.template_formula(node)?;
        let mut neighbours = others.clone();
        neighbours.push(previous);
        neighbours.push(next);
        let first = others.len();
        let Some(positions) = self.template_numbering(node, &neighbours, |p| {
            let (before, after) = (p[first], p[first + 1]);
            p[..first].iter().all(|&o| before < o && o < after) && before < after
        }) else {
            let mut all = others;
            all.push(next);
            return self.template_bracketed(node, previous, &all, out);
        };
        let mut listed = positions.clone();
        listed.sort_unstable();
        let locants = listed
            .iter()
            .map(u8::to_string)
            .collect::<Vec<_>>()
            .join(",");
        if !out.ends_with('-') {
            out.push('-');
        }
        out.push_str(&format!("{locants}-{formula}"));
        let mut ordered: Vec<(u8, usize)> = positions[..first]
            .iter()
            .copied()
            .zip(others.iter().copied())
            .collect();
        ordered.sort_unstable();
        for (_, other) in ordered {
            let text = self.substituent(node, other)?;
            if text.starts_with('(') && text.ends_with(')') {
                out.push_str(&text);
            } else {
                out.push_str(&format!("({text})"));
            }
        }
        self.chain(rest, Some(node), out)
    }
}
