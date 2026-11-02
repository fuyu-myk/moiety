//! Folding Unicode formula notation to the ASCII form the lexer reads (grammar, section 2).

use crate::span::Span;

/// One substitution made by [`normalize`], mapping a normalized span back to the original text.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fold {
    /// Span in the normalized string. Empty when original text was removed.
    pub normalized: Span,
    /// Span in the original input.
    pub original: Span,
}

/// A formula after normalization, with the folds needed to map spans back.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Normalized {
    text: String,
    folds: Vec<Fold>,
}

impl Normalized {
    /// The normalized ASCII-folded text.
    #[must_use]
    pub fn text(&self) -> &str {
        &self.text
    }

    /// Every substitution, in order.
    #[must_use]
    pub fn folds(&self) -> &[Fold] {
        &self.folds
    }

    /// Maps a span of the normalized text to the original input.
    #[must_use]
    pub fn to_original(&self, span: Span) -> Span {
        to_original(&self.folds, span)
    }
}

/// Maps a span of normalized text to the original input through a fold table.
#[must_use]
pub fn to_original(folds: &[Fold], span: Span) -> Span {
    Span::new(map(folds, span.start, false), map(folds, span.end, true))
}

fn map(folds: &[Fold], offset: usize, is_end: bool) -> usize {
    let mut anchor = (0, 0);
    for fold in folds {
        let n = fold.normalized;
        if offset < n.start || (offset == n.start && (is_end || !n.is_empty())) {
            break;
        }
        if offset >= n.end {
            anchor = (n.end, fold.original.end);
            continue;
        }
        return if is_end {
            fold.original.end
        } else {
            fold.original.start
        };
    }
    offset - anchor.0 + anchor.1
}

/// Folds a formula to the ASCII form the lexer reads.
///
/// Subscript digits become digits, runs of superscript digits and signs become a
/// `^` charge, the middle dot, bullet and asterisk become `.`, `≡` becomes `#`,
/// dash-like characters become `-`, and whitespace is removed. Every change is
/// recorded so spans can be mapped back with [`Normalized::to_original`].
#[must_use]
pub fn normalize(input: &str) -> Normalized {
    let mut text = String::with_capacity(input.len());
    let mut folds = Vec::new();
    let mut chars = input.char_indices().peekable();
    while let Some((start, c)) = chars.next() {
        let end = start + c.len_utf8();
        if c.is_whitespace() {
            folds.push(Fold {
                normalized: Span::new(text.len(), text.len()),
                original: Span::new(start, end),
            });
            continue;
        }
        if superscript(c).is_some() {
            let mut run = String::new();
            let mut run_end = end;
            run.push(superscript(c).expect("checked above"));
            while let Some(&(next_start, next)) = chars.peek() {
                let Some(folded) = superscript(next) else {
                    break;
                };
                run.push(folded);
                run_end = next_start + next.len_utf8();
                chars.next();
            }
            let from = text.len();
            text.push('^');
            text.push_str(&run);
            folds.push(Fold {
                normalized: Span::new(from, text.len()),
                original: Span::new(start, run_end),
            });
            continue;
        }
        let folded = match c {
            '₀'..='₉' => char::from_digit(u32::from(c) - u32::from('₀'), 10),
            '·' | '•' | '∙' | '*' => Some('.'),
            '≡' => Some('#'),
            '\u{2010}'..='\u{2014}' | '\u{2212}' => Some('-'),
            _ => None,
        };
        match folded {
            Some(replacement) if replacement != c => {
                let from = text.len();
                text.push(replacement);
                folds.push(Fold {
                    normalized: Span::new(from, text.len()),
                    original: Span::new(start, end),
                });
            }
            _ => text.push(c),
        }
    }
    Normalized { text, folds }
}

fn superscript(c: char) -> Option<char> {
    match c {
        '⁰' => Some('0'),
        '¹' => Some('1'),
        '²' => Some('2'),
        '³' => Some('3'),
        '⁴'..='⁹' => char::from_digit(u32::from(c) - u32::from('⁴') + 4, 10),
        '⁺' => Some('+'),
        '⁻' => Some('-'),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn folds_unicode_to_ascii() {
        assert_eq!(normalize("CH₃CH₂OH").text(), "CH3CH2OH");
        assert_eq!(normalize("HC≡CH").text(), "HC#CH");
        assert_eq!(normalize("CuSO4·5H2O").text(), "CuSO4.5H2O");
        assert_eq!(normalize("SO₄²⁻").text(), "SO4^2-");
        assert_eq!(normalize(" CH3 \u{2013} CH3 ").text(), "CH3-CH3");
        assert_eq!(normalize("CH3CH2OH").folds(), &[]);
    }

    #[test]
    fn maps_spans_back_through_folds() {
        let n = normalize("CH₃CH₂OH");
        assert_eq!(n.to_original(Span::new(0, 3)), Span::new(0, 5));
        assert_eq!(n.to_original(Span::new(3, 6)), Span::new(5, 10));
        assert_eq!(n.to_original(Span::new(6, 8)), Span::new(10, 12));
        assert_eq!(n.to_original(Span::new(2, 3)), Span::new(2, 5));
    }

    #[test]
    fn maps_spans_around_removed_whitespace() {
        let n = normalize("  CH3 CH3");
        assert_eq!(n.text(), "CH3CH3");
        assert_eq!(n.to_original(Span::new(0, 3)), Span::new(2, 5));
        assert_eq!(n.to_original(Span::new(3, 6)), Span::new(6, 9));
    }

    #[test]
    fn maps_superscript_charge_runs() {
        let n = normalize("SO₄²⁻");
        assert_eq!(n.to_original(Span::new(3, 6)), Span::new(5, 10));
    }
}
