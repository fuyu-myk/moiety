use crate::dictionary::Dictionary;
use std::cell::Cell;

use crate::interpret::{Options, interpret_within};
use crate::reading::Interpretation;
use crate::span::Span;

const MAX_SUGGESTION_INPUT: usize = 40;
const SUGGESTION_SEARCH_LIMIT: usize = 20_000;
const MAX_CANDIDATES: usize = 120;
const DIGITS: [char; 3] = ['2', '3', '4'];
const LETTERS: [char; 5] = ['C', 'H', 'O', 'N', 'S'];

/// Searches single-character edits of `input`, nearest to `near` first, for one that reads.
pub(crate) fn suggest(
    input: &str,
    near: Span,
    dictionary: &Dictionary,
    options: &Options,
) -> Option<String> {
    let chars: Vec<(usize, char)> = input.char_indices().collect();
    if chars.is_empty() || input.len() > MAX_SUGGESTION_INPUT {
        return None;
    }
    let quiet = Options {
        suggestions: false,
        ..options.clone()
    };
    let budget = Cell::new(SUGGESTION_SEARCH_LIMIT);
    let distance = |offset: usize| {
        if near.start <= offset && offset < near.end.max(near.start + 1) {
            0
        } else {
            offset.abs_diff(near.start).min(offset.abs_diff(near.end))
        }
    };
    let mut positions: Vec<usize> = (0..chars.len()).collect();
    positions.sort_by_key(|&i| (distance(chars[i].0), std::cmp::Reverse(i)));
    let tried = Cell::new(0);
    let reads = |candidate: &str| {
        tried.set(tried.get() + 1);
        tried.get() <= MAX_CANDIDATES
            && budget.get() > 0
            && matches!(
                interpret_within(candidate, dictionary, &quiet, Some(&budget)),
                Interpretation::Unique(_) | Interpretation::Ambiguous { .. }
            )
    };
    let edit = |i: usize, replacement: Option<char>, insert: bool| {
        let mut out: String = chars[..i].iter().map(|&(_, c)| c).collect();
        if let Some(c) = replacement {
            out.push(c);
        }
        let rest = if insert { i } else { i + 1 };
        out.extend(chars[rest..].iter().map(|&(_, c)| c));
        out
    };
    for &i in &positions {
        let original = chars[i].1;
        let mut replacements: Vec<char> = if original.is_ascii_digit() {
            DIGITS.into_iter().chain(LETTERS).collect()
        } else if original.is_ascii_alphabetic() {
            LETTERS.to_vec()
        } else {
            Vec::new()
        };
        replacements.retain(|&c| c != original);
        if original.is_ascii_alphabetic() {
            let flipped = if original.is_ascii_uppercase() {
                original.to_ascii_lowercase()
            } else {
                original.to_ascii_uppercase()
            };
            replacements.insert(0, flipped);
        }
        for c in replacements {
            let candidate = edit(i, Some(c), false);
            if reads(&candidate) {
                return Some(candidate);
            }
        }
        let deleted = edit(i, None, false);
        if !deleted.is_empty() && reads(&deleted) {
            return Some(deleted);
        }
    }
    for &i in &positions {
        for c in DIGITS.into_iter().chain(LETTERS) {
            let candidate = edit(i, Some(c), true);
            if reads(&candidate) {
                return Some(candidate);
            }
        }
    }
    None
}
