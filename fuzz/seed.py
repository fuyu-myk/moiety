"""Write the fuzzing seed corpus and token dictionary from the repository's own data.

Usage, from the repository root:

    python3 -I fuzz/seed.py

Seeds are every conformance input, every differential label with its
attachment marker, and the playground examples. Tokens are every element
symbol, every dictionary spelling, and the characters the lexer and the
normalizer treat specially.
"""

import hashlib
import re
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CORPUS = ROOT / "fuzz" / "corpus" / "playground"
DICTIONARY = ROOT / "fuzz" / "formula.dict"

PUNCTUATION = [
    "(", ")", "[", "]", "{", "}", "-", "+", "^", ".", "*", "=", "#", " ", ",",
    "o-", "m-", "p-", "C6H4", "C6H3",
    "2", "3", "4", "6", "10", "100",
    "·", "≡", "−", "₂", "₃", "₁₀", "²", "³", "⁺", "⁻", "²⁺",
]


def seeds() -> list[str]:
    cases = tomllib.loads((ROOT / "conformance" / "cases.toml").read_text(encoding="utf-8"))
    found = [case["input"] for case in cases["case"]]
    labels = ROOT / "conformance" / "differential" / "labels.txt"
    for line in labels.read_text(encoding="utf-8").splitlines():
        content = line.split("#", 1)[0].split()
        if len(content) == 2:
            text, side = content
            found.append(f"-{text}" if side == "left" else f"{text}-")
    page = (ROOT / "playground" / "index.html").read_text(encoding="utf-8")
    found += re.findall(r'<button type="button">([^<]+)</button>', page)
    return found


def tokens() -> list[str]:
    elements = (ROOT / "crates" / "moiety-core" / "src" / "element.rs").read_text(encoding="utf-8")
    found = re.findall(r'd\("([A-Z][a-z]?)"', elements)
    dictionary = tomllib.loads(
        (ROOT / "crates" / "moiety" / "data" / "dictionary.toml").read_text(encoding="utf-8")
    )
    for entry in dictionary["entry"]:
        found += entry.get("spellings", []) + entry.get("reversed", [])
    return list(dict.fromkeys(found + PUNCTUATION))


def quoted(token: str) -> str:
    escaped = "".join(
        f"\\x{byte:02X}" if byte > 0x7E or byte < 0x20 or byte in b'"\\' else chr(byte)
        for byte in token.encode("utf-8")
    )
    return f'"{escaped}"'


def main() -> None:
    CORPUS.mkdir(parents=True, exist_ok=True)
    for text in seeds():
        data = text.encode("utf-8")
        (CORPUS / hashlib.sha1(data).hexdigest()).write_bytes(data)
    DICTIONARY.write_text("\n".join(quoted(t) for t in tokens()) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
