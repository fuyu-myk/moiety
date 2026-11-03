//! Benzene rings the dictionary does not cover, written with locants (grammar,
//! section 6.11): after the ring as `CH3C6H4-4-CHO` when a unit comes before it,
//! with locants first as `1,4-Cl2C6H4` when it starts the formula, and with a
//! bare locant as `-1,4-C6H4-` between two attachment points.

use moiety_core::{AtomKind, BondOrder, Element, Molecule};

use super::{Node, WriteError, Writer, push_label};

/// A benzene ring the dictionary does not cover.
pub(super) struct Ring {
    /// The ring atoms in cyclic order.
    pub atoms: [usize; 6],
}

/// Finds the benzene rings among atoms no dictionary group covers, and marks them
/// covered: six aromatic carbons in one ring, each with a hydrogen or a single
/// bond to one substituent that is in no other ring with it.
pub(super) fn benzene_rings(
    molecule: &Molecule,
    adjacency: &[Vec<(usize, usize)>],
    covered: &mut [bool],
) -> Vec<Ring> {
    let aromatic = |bond: usize| molecule.bonds()[bond].order == BondOrder::Aromatic;
    let bridge = |bond: usize| {
        let (a, b) = (molecule.bonds()[bond].a, molecule.bonds()[bond].b);
        let mut seen = vec![false; molecule.atoms().len()];
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
        !seen[b]
    };
    let ring_neighbours = |atom: usize| -> Vec<usize> {
        adjacency[atom]
            .iter()
            .filter(|&&(_, bond)| aromatic(bond))
            .map(|&(next, _)| next)
            .collect()
    };
    let candidate: Vec<bool> = (0..molecule.atoms().len())
        .map(|atom| {
            let a = molecule.atoms()[atom];
            let external: Vec<usize> = adjacency[atom]
                .iter()
                .filter(|&&(_, bond)| !aromatic(bond))
                .map(|&(_, bond)| bond)
                .collect();
            !covered[atom]
                && a.kind == AtomKind::Element(Element::C)
                && a.aromatic
                && a.charge == 0
                && ring_neighbours(atom).len() == 2
                && usize::from(a.hydrogens) + external.len() == 1
                && external
                    .iter()
                    .all(|&bond| molecule.bonds()[bond].order == BondOrder::Single && bridge(bond))
        })
        .collect();
    let mut rings = Vec::new();
    for start in 0..molecule.atoms().len() {
        if !candidate[start] || covered[start] {
            continue;
        }
        let mut cycle = vec![start];
        let (mut previous, mut current) = (start, ring_neighbours(start)[0]);
        while current != start && cycle.len() < 6 && candidate[current] {
            cycle.push(current);
            let Some(next) = ring_neighbours(current)
                .into_iter()
                .find(|&n| n != previous)
            else {
                break;
            };
            (previous, current) = (current, next);
        }
        if current == start {
            if let Ok(atoms) = <[usize; 6]>::try_from(cycle) {
                for &atom in &atoms {
                    covered[atom] = true;
                }
                rings.push(Ring { atoms });
            }
        }
    }
    rings
}

/// The position, 1 to 6, of each ring site when the site at cycle index `first` is
/// position 1, numbered in the direction that gives the lowest locants or, with
/// `last`, in a direction that leaves that site the highest position.
fn positions(first: usize, sites: &[usize], last: Option<usize>) -> Option<Vec<u8>> {
    let number = |forward: bool, site: usize| {
        let step = if forward {
            (site + 6 - first) % 6
        } else {
            (first + 6 - site) % 6
        };
        u8::try_from(step).expect("below 6") + 1
    };
    let mut options: Vec<Vec<u8>> = [true, false]
        .into_iter()
        .map(|forward| sites.iter().map(|&s| number(forward, s)).collect())
        .collect();
    if let Some(last) = last {
        options.retain(|numbers: &Vec<u8>| {
            let tail = numbers[last];
            numbers.iter().all(|&n| n <= tail)
        });
    }
    options.into_iter().min_by_key(|numbers| {
        let mut sorted = numbers.clone();
        sorted.sort_unstable();
        sorted
    })
}

fn locants(positions: &[u8]) -> String {
    positions
        .iter()
        .map(u8::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

/// A substituent written `count` times, as `Cl2` or `(OCH3)2`.
fn counted(text: &str, count: usize) -> String {
    if count == 1 {
        return text.to_string();
    }
    let single = matches!(text, "F" | "Cl" | "Br" | "I") || wrapped(text);
    if single {
        format!("{text}{count}")
    } else {
        format!("({text}){count}")
    }
}

/// Whether `text` is one bracketed group from its first character to its last.
fn wrapped(text: &str) -> bool {
    if !text.starts_with('(') {
        return false;
    }
    let mut depth = 0usize;
    for (index, c) in text.char_indices() {
        match c {
            '(' => depth += 1,
            ')' => {
                depth -= 1;
                if depth == 0 {
                    return index == text.len() - 1;
                }
            }
            _ => {}
        }
    }
    false
}

impl Writer<'_> {
    fn ring(&self, node: usize) -> &super::rings::Ring {
        let Node::Ring(r) = self.nodes[node] else {
            unreachable!("only ring nodes are numbered")
        };
        &self.rings[r]
    }

    /// The cycle index of the ring atom bonded to `neighbour`.
    fn site(&self, node: usize, neighbour: usize) -> usize {
        let atoms = self.ring(node).atoms;
        self.molecule
            .bonds()
            .iter()
            .find_map(|bond| {
                let ends = [(bond.a, bond.b), (bond.b, bond.a)];
                ends.into_iter().find_map(|(inside, outside)| {
                    let index = atoms.iter().position(|&a| a == inside)?;
                    (self.node_of.get(&outside) == Some(&neighbour)).then_some(index)
                })
            })
            .expect("a ring node is bonded to each of its neighbours")
    }

    /// The ring's formula, such as `C6H4`, when locants can number it.
    fn ring_formula(&self, node: usize) -> Result<String, WriteError> {
        match 6 - self.edges[node].len() {
            1 => Ok("C6H".to_string()),
            hydrogens @ 2..=4 => Ok(format!("C6H{hydrogens}")),
            _ => Err(WriteError::UnknownRing(self.ring(node).atoms[0])),
        }
    }

    /// Whether any node of the substituent rooted at `root`, away from `node`,
    /// satisfies `test`.
    fn substituent_has(&self, node: usize, root: usize, test: impl Fn(usize) -> bool) -> bool {
        let mut stack = vec![(root, node)];
        while let Some((current, parent)) = stack.pop() {
            if test(current) {
                return true;
            }
            stack.extend(
                self.edges[current]
                    .iter()
                    .filter(|&&(next, _)| next != parent)
                    .map(|&(next, _)| (next, current)),
            );
        }
        false
    }

    /// A substituent of `node` rooted at `root`, bracketed when it holds a ring or
    /// a dictionary group, which could otherwise be read as the ring numbered.
    pub(super) fn substituent(&self, node: usize, root: usize) -> Result<String, WriteError> {
        let path = self.longest_from(root, Some(node));
        let mut text = String::new();
        self.chain(&path, Some(node), &mut text)?;
        let nested = self.substituent_has(node, root, |n| {
            matches!(self.nodes[n], Node::Ring(_) | Node::Group(_))
        });
        Ok(if nested { format!("({text})") } else { text })
    }

    /// The substituents at `entries`, each a position and a neighbour, with equal
    /// uncharged ones written once and counted, ordered by their first position
    /// except that one with a charged atom comes last, where the written charge
    /// can reach it.
    pub(super) fn chunks(
        &self,
        node: usize,
        entries: &[(u8, usize)],
    ) -> Result<Vec<(String, String)>, WriteError> {
        let mut sorted = entries.to_vec();
        sorted.sort_unstable();
        let mut grouped: Vec<(Vec<u8>, String, bool)> = Vec::new();
        for (position, neighbour) in sorted {
            let text = self.substituent(node, neighbour)?;
            let charged = self.substituent_has(
                node,
                neighbour,
                |n| matches!(self.nodes[n], Node::Atom(a) if self.molecule.atoms()[a].charge != 0),
            );
            match grouped
                .iter_mut()
                .find(|(_, t, c)| *t == text && !c && !charged)
            {
                Some((positions, _, _)) => positions.push(position),
                None => grouped.push((vec![position], text, charged)),
            }
        }
        grouped.sort_by_key(|(_, _, charged)| *charged);
        Ok(grouped
            .into_iter()
            .map(|(positions, text, _)| (locants(&positions), counted(&text, positions.len())))
            .collect())
    }

    /// Writes a ring with a unit before it, which bonds at position 1. Its other
    /// substituents follow it as `C6H4-4-OMe`. The units of `rest` continue as the
    /// last when they lead to an attachment point, and are one more substituent
    /// otherwise.
    pub(super) fn ring_after(
        &self,
        node: usize,
        previous: usize,
        rest: &[usize],
        out: &mut String,
    ) -> Result<(), WriteError> {
        let next = rest
            .first()
            .copied()
            .filter(|&n| self.substituent_has(node, n, |m| self.is_attachment(m)));
        let others: Vec<usize> = self.edges[node]
            .iter()
            .map(|&(other, _)| other)
            .filter(|&o| o != previous && Some(o) != next)
            .collect();
        let first = self.site(node, previous);
        let mut sites: Vec<usize> = others.iter().map(|&o| self.site(node, o)).collect();
        let ring_ahead = rest.iter().any(|&n| matches!(self.nodes[n], Node::Ring(_)));
        if let Some(next) = next.filter(|&n| self.is_attachment(n) || ring_ahead) {
            sites.push(self.site(node, next));
            return self.ring_between(node, &others, &sites, first, rest, out);
        }
        if let Some(next) = next {
            sites.push(self.site(node, next));
        }
        let numbers = positions(first, &sites, None).expect("both directions are allowed");
        let entries: Vec<(u8, usize)> = numbers
            .iter()
            .copied()
            .zip(others.iter().copied())
            .collect();
        let chunks = self.chunks(node, &entries)?;
        if chunks.is_empty() && next.is_none() {
            return Err(WriteError::UnknownRing(self.ring(node).atoms[first]));
        }
        push_label(out, &self.ring_formula(node)?);
        for (locants, text) in chunks {
            out.push_str(&format!("-{locants}-{text}"));
        }
        if next.is_some() {
            let position = numbers.last().expect("the next unit has a position");
            out.push_str(&format!("-{position}-"));
            self.chain(rest, Some(node), out)?;
        }
        Ok(())
    }

    /// Writes a ring on the path between two attachment points with a bare locant,
    /// as `-1,4-C6H4-`. Other substituents are bracketed after it and must take
    /// lower positions than the next unit on the path, which bonds last.
    fn ring_between(
        &self,
        node: usize,
        others: &[usize],
        sites: &[usize],
        first: usize,
        rest: &[usize],
        out: &mut String,
    ) -> Result<(), WriteError> {
        let numbers = positions(first, sites, Some(sites.len() - 1))
            .ok_or(WriteError::UnknownRing(self.ring(node).atoms[first]))?;
        let mut listed: Vec<u8> = numbers.clone();
        listed.push(1);
        listed.sort_unstable();
        if !out.ends_with('-') {
            out.push('-');
        }
        out.push_str(&format!(
            "{}-{}",
            locants(&listed),
            self.ring_formula(node)?
        ));
        let mut ordered: Vec<(u8, usize)> = numbers
            .iter()
            .copied()
            .zip(others.iter().copied())
            .collect();
        ordered.sort_unstable();
        for (_, other) in ordered {
            let text = self.substituent(node, other)?;
            if wrapped(&text) {
                out.push_str(&text);
            } else {
                out.push_str(&format!("({text})"));
            }
        }
        self.chain(rest, Some(node), out)
    }

    /// Writes a ring that starts the formula. Its substituents come first with
    /// their locants, as `2,4-Cl2C6H3`, and `next` bonds at position 1; with no
    /// `next`, one substituent takes position 1, as in `1,4-Cl2C6H4`.
    pub(super) fn ring_first(
        &self,
        node: usize,
        next: Option<usize>,
        out: &mut String,
    ) -> Result<(), WriteError> {
        let others: Vec<usize> = self.edges[node]
            .iter()
            .map(|&(other, _)| other)
            .filter(|&o| Some(o) != next)
            .collect();
        if others.is_empty() && next.is_none() {
            push_label(out, "C6H6");
            return Ok(());
        }
        let sites: Vec<usize> = others.iter().map(|&o| self.site(node, o)).collect();
        let numbers = match next {
            Some(next) => positions(self.site(node, next), &sites, None),
            None => sites
                .iter()
                .filter_map(|&first| positions(first, &sites, None))
                .min_by_key(|numbers| {
                    let mut sorted = numbers.clone();
                    sorted.sort_unstable();
                    sorted
                }),
        }
        .expect("both directions are allowed");
        let entries: Vec<(u8, usize)> = numbers
            .iter()
            .copied()
            .zip(others.iter().copied())
            .collect();
        let chunks = self.chunks(node, &entries)?;
        if chunks.is_empty() {
            return Err(WriteError::UnknownRing(self.ring(node).atoms[0]));
        }
        let prefix = chunks
            .iter()
            .map(|(locants, text)| format!("{locants}-{text}"))
            .collect::<Vec<_>>()
            .join("-");
        push_label(out, &prefix);
        out.push_str(&self.ring_formula(node)?);
        Ok(())
    }
}
