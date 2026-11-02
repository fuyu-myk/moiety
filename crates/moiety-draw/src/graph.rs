//! Breadth-first searches over a molecule's adjacency lists.

use std::collections::VecDeque;

use moiety_core::{BondOrder, Molecule};

/// Adjacency lists as returned by [`moiety_core::Molecule::adjacency`]: `(neighbour, bond)` pairs.
pub(crate) type Adjacency = [Vec<(usize, usize)>];

/// The shortest path from `from` to `to` that does not use bond `skip`, endpoints included.
pub(crate) fn shortest_path(
    adjacency: &Adjacency,
    from: usize,
    to: usize,
    skip: Option<usize>,
) -> Option<Vec<usize>> {
    let mut previous = vec![None; adjacency.len()];
    let mut queue = VecDeque::from([from]);
    previous[from] = Some(from);
    while let Some(atom) = queue.pop_front() {
        if atom == to {
            let mut path = vec![to];
            let mut current = to;
            while current != from {
                current = previous[current].expect("every reached atom has a predecessor");
                path.push(current);
            }
            path.reverse();
            return Some(path);
        }
        for &(next, bond) in &adjacency[atom] {
            if Some(bond) != skip && previous[next].is_none() {
                previous[next] = Some(atom);
                queue.push_back(next);
            }
        }
    }
    None
}

/// Every atom reachable from `start` without passing through `blocked`, in breadth-first order.
pub(crate) fn side(adjacency: &Adjacency, start: usize, blocked: usize) -> Vec<usize> {
    let mut seen = vec![false; adjacency.len()];
    seen[start] = true;
    seen[blocked] = true;
    let mut members = vec![start];
    let mut cursor = 0;
    while let Some(&atom) = members.get(cursor) {
        cursor += 1;
        for &(next, _) in &adjacency[atom] {
            if !seen[next] {
                seen[next] = true;
                members.push(next);
            }
        }
    }
    members
}

/// The atom farthest from `start` by bond count, the lowest index on ties.
pub(crate) fn farthest(adjacency: &Adjacency, start: usize) -> usize {
    let mut distance = vec![usize::MAX; adjacency.len()];
    distance[start] = 0;
    let mut queue = VecDeque::from([start]);
    let mut best = start;
    while let Some(atom) = queue.pop_front() {
        let here = distance[atom];
        if here > distance[best] || (here == distance[best] && atom < best) {
            best = atom;
        }
        for &(next, _) in &adjacency[atom] {
            if distance[next] == usize::MAX {
                distance[next] = here + 1;
                queue.push_back(next);
            }
        }
    }
    best
}

/// Whether `atom` has two neighbours in a straight line: a triple bond, or two double bonds.
pub(crate) fn is_linear(molecule: &Molecule, adjacency: &Adjacency, atom: usize) -> bool {
    let orders: Vec<BondOrder> = adjacency[atom]
        .iter()
        .map(|&(_, bond)| molecule.bonds()[bond].order)
        .collect();
    orders.len() == 2
        && (orders.contains(&BondOrder::Triple)
            || orders.iter().all(|&order| order == BondOrder::Double))
}

/// The bond joining `a` and `b`, if they are bonded.
pub(crate) fn bond_between(adjacency: &Adjacency, a: usize, b: usize) -> Option<usize> {
    adjacency[a]
        .iter()
        .find(|&&(other, _)| other == b)
        .map(|&(_, bond)| bond)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chain(n: usize) -> Vec<Vec<(usize, usize)>> {
        let mut adjacency = vec![Vec::new(); n];
        for i in 1..n {
            adjacency[i - 1].push((i, i - 1));
            adjacency[i].push((i - 1, i - 1));
        }
        adjacency
    }

    #[test]
    fn finds_paths_and_sides() {
        let adjacency = chain(5);
        assert_eq!(
            shortest_path(&adjacency, 0, 3, None),
            Some(vec![0, 1, 2, 3])
        );
        assert_eq!(shortest_path(&adjacency, 0, 3, Some(1)), None);
        assert_eq!(side(&adjacency, 2, 1), vec![2, 3, 4]);
        assert_eq!(farthest(&adjacency, 2), 0);
        assert_eq!(bond_between(&adjacency, 3, 2), Some(2));
        assert_eq!(bond_between(&adjacency, 0, 2), None);
    }
}
