//! Escaping text for XML and HTML.

/// `text` with the five XML special characters escaped and control characters
/// other than tab, newline and carriage return replaced by U+FFFD, so the result
/// is safe in element content and in quoted attribute values.
pub(crate) fn escape(text: &str) -> String {
    let mut escaped = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '&' => escaped.push_str("&amp;"),
            '<' => escaped.push_str("&lt;"),
            '>' => escaped.push_str("&gt;"),
            '"' => escaped.push_str("&quot;"),
            '\'' => escaped.push_str("&#39;"),
            '\t' | '\n' | '\r' => escaped.push(c),
            c if c.is_control() => escaped.push(char::REPLACEMENT_CHARACTER),
            c => escaped.push(c),
        }
    }
    escaped
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn escapes_markup_and_controls() {
        assert_eq!(
            escape(r#"<a href="x">&'</a>"#),
            "&lt;a href=&quot;x&quot;&gt;&amp;&#39;&lt;/a&gt;"
        );
        assert_eq!(escape("C\u{0}H\n"), "C\u{FFFD}H\n");
        assert_eq!(escape("CuSO4·5H2O"), "CuSO4·5H2O");
    }
}
