//! Stress majorization for ring systems that no template fits, such as cages
//! and bridged polycycles.
//!
//! Every pair of atoms in the system is drawn toward a target distance that
//! grows with the number of bonds between them, bonded pairs most strongly.
//! Each sweep moves one atom at a time to the weighted mean of where its pairs
//! would put it (the localized SMACOF update), which never increases the stress.

use std::collections::VecDeque;

use super::geometry::{EPSILON, Point};
use crate::graph::Adjacency;

/// Number of passes over all atoms.
const SWEEPS: usize = 300;

/// Target distance added per bond beyond the first, in bond lengths.
const STEP: f64 = 0.7;

/// Turn between successive fallback directions, so no two coincide.
const GOLDEN_ANGLE: f64 = 2.399_963_229_728_653;

/// Size of the deterministic nudge that separates atoms placed on the same spot.
const NUDGE: f64 = 0.05;

/// Moves the atoms of `system` (all placed in `local`) to a layout with bonds near unit length and no atoms piled up.
///
/// Two starts are relaxed, the given positions and a circle in depth-first
/// order, and the one that ends with less stress is kept.
pub(super) fn relax(system: &[usize], adjacency: &Adjacency, local: &mut [Option<Point>]) {
    let targets = targets(system, adjacency);
    let given: Vec<Point> = system
        .iter()
        .enumerate()
        .map(|(k, &atom)| {
            local[atom].unwrap_or_default() + Point::from_angle(GOLDEN_ANGLE * k as f64) * NUDGE
        })
        .collect();
    let points = [given, circle(system, adjacency)]
        .into_iter()
        .map(|start| settle(start, &targets))
        .min_by(|a, b| stress(a, &targets).total_cmp(&stress(b, &targets)))
        .expect("relaxation has two starts");
    let bonds: Vec<f64> = (0..points.len())
        .flat_map(|i| (i + 1..points.len()).map(move |j| (i, j)))
        .filter(|&(i, j)| targets[i][j] == 1.0)
        .map(|(i, j)| points[i].distance(points[j]))
        .collect();
    let mean = bonds.iter().sum::<f64>() / bonds.len().max(1) as f64;
    let scale = if mean > EPSILON { mean.recip() } else { 1.0 };
    for (&atom, point) in system.iter().zip(points) {
        local[atom] = Some(point * scale);
    }
}

fn settle(mut points: Vec<Point>, targets: &[Vec<f64>]) -> Vec<Point> {
    for _ in 0..SWEEPS {
        for i in 0..points.len() {
            points[i] = update(i, &points, &targets[i]);
        }
    }
    points
}

/// Weighted squared misfit of every pair against its target distance.
fn stress(points: &[Point], targets: &[Vec<f64>]) -> f64 {
    (0..points.len())
        .flat_map(|i| (i + 1..points.len()).map(move |j| (i, j)))
        .map(|(i, j)| ((points[i].distance(points[j]) - targets[i][j]) / targets[i][j]).powi(2))
        .sum()
}

/// The system atoms on a circle with unit spacing, in depth-first order from the first atom.
fn circle(system: &[usize], adjacency: &Adjacency) -> Vec<Point> {
    let mut order = vec![None; adjacency.len()];
    let mut rank = 0;
    let mut stack = vec![system[0]];
    while let Some(atom) = stack.pop() {
        if order[atom].is_some() || system.binary_search(&atom).is_err() {
            continue;
        }
        order[atom] = Some(rank);
        rank += 1;
        stack.extend(adjacency[atom].iter().rev().map(|&(next, _)| next));
    }
    let count = system.len() as f64;
    system
        .iter()
        .enumerate()
        .map(|(k, &atom)| {
            let place = order[atom].unwrap_or(k) as f64;
            Point::from_angle(std::f64::consts::TAU * place / count)
                * (count / std::f64::consts::TAU)
        })
        .collect()
}

/// Where atom `i` minimizes its share of the stress, given every other atom.
fn update(i: usize, points: &[Point], targets: &[f64]) -> Point {
    let (sum, total) = points
        .iter()
        .zip(targets)
        .enumerate()
        .filter(|&(j, _)| j != i)
        .fold(
            (Point::default(), 0.0),
            |(sum, total), (j, (&other, &target))| {
                let weight = target.powi(-2);
                let direction = (points[i] - other).normalized().unwrap_or_else(|| {
                    Point::from_angle(GOLDEN_ANGLE * (i * points.len() + j) as f64)
                });
                (sum + (other + direction * target) * weight, total + weight)
            },
        );
    if total > 0.0 {
        sum * total.recip()
    } else {
        points[i]
    }
}

/// Target distance between every pair of system atoms, from bond counts inside the system.
fn targets(system: &[usize], adjacency: &Adjacency) -> Vec<Vec<f64>> {
    let mut slot = vec![None; adjacency.len()];
    for (k, &atom) in system.iter().enumerate() {
        slot[atom] = Some(k);
    }
    let far = system.len() as f64;
    system
        .iter()
        .map(|&start| {
            let mut bonds = vec![None; system.len()];
            let mut queue = VecDeque::from([(start, 0usize)]);
            if let Some(k) = slot[start] {
                bonds[k] = Some(0);
            }
            while let Some((atom, count)) = queue.pop_front() {
                for &(next, _) in &adjacency[atom] {
                    if let Some(k) = slot[next].filter(|&k| bonds[k].is_none()) {
                        bonds[k] = Some(count + 1);
                        queue.push_back((next, count + 1));
                    }
                }
            }
            bonds
                .into_iter()
                .map(|count| count.map_or(far, |c| 1.0 + STEP * (c.max(1) - 1) as f64))
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spreads_a_square_drawn_on_one_point() {
        let adjacency: Vec<Vec<(usize, usize)>> = vec![
            vec![(1, 0), (3, 3)],
            vec![(0, 0), (2, 1)],
            vec![(1, 1), (3, 2)],
            vec![(2, 2), (0, 3)],
        ];
        let mut local = vec![Some(Point::default()); 4];
        relax(&[0, 1, 2, 3], &adjacency, &mut local);
        let p: Vec<Point> = local.into_iter().flatten().collect();
        for i in 0..4 {
            let bond = p[i].distance(p[(i + 1) % 4]);
            assert!((bond - 1.0).abs() < 0.05, "bond {i}: {bond}");
        }
        assert!(p[0].distance(p[2]) > 1.2);
    }
}
