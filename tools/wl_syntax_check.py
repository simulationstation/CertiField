"""Very light syntax check for Wolfram Language files.

Does not invoke the kernel (keeps CI free).  Catches the common
authoring bugs: unmatched brackets, stray back-ticks, missing
BeginPackage/EndPackage pairs.
"""
from __future__ import annotations
import sys
from pathlib import Path


def check(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    errors: list[str] = []

    # Bracket balance (ignores brackets inside strings and (* *) comments).
    pairs = {"(": ")", "[": "]", "{": "}"}
    stack: list[tuple[str, int]] = []
    i, n = 0, len(text)
    line = 1
    in_string = False
    comment_depth = 0
    while i < n:
        c = text[i]
        if c == "\n":
            line += 1
        if comment_depth > 0:
            if c == "(" and i + 1 < n and text[i + 1] == "*":
                comment_depth += 1
                i += 2
                continue
            if c == "*" and i + 1 < n and text[i + 1] == ")":
                comment_depth -= 1
                i += 2
                continue
            i += 1
            continue
        if in_string:
            if c == "\\":
                i += 2
                continue
            if c == '"':
                in_string = False
            i += 1
            continue
        if c == '"':
            in_string = True
            i += 1
            continue
        if c == "(" and i + 1 < n and text[i + 1] == "*":
            comment_depth = 1
            i += 2
            continue
        if c in pairs:
            stack.append((c, line))
        elif c in pairs.values():
            if not stack:
                errors.append(f"{path}:{line}: stray '{c}'")
            else:
                opener, opened_at = stack.pop()
                if pairs[opener] != c:
                    errors.append(
                        f"{path}:{line}: expected '{pairs[opener]}' "
                        f"to close '{opener}' from line {opened_at}, got '{c}'"
                    )
        i += 1

    if comment_depth > 0:
        errors.append(f"{path}: unterminated comment (comments nest in Wolfram Language: a '(*' inside a comment needs its own '*)')")
    if in_string:
        errors.append(f"{path}: unterminated string")
    for opener, opened_at in stack:
        errors.append(f"{path}:{opened_at}: unclosed '{opener}'")

    # BeginPackage/EndPackage balance for .wl files.
    if path.suffix == ".wl":
        begins = text.count("BeginPackage[")
        ends = text.count("EndPackage[]")
        if begins != ends:
            errors.append(
                f"{path}: BeginPackage count ({begins}) != EndPackage count ({ends})"
            )

    return errors


def main(argv: list[str]) -> int:
    paths = [Path(p) for p in argv[1:]]
    if not paths:
        print("usage: wl_syntax_check.py FILE [FILE ...]", file=sys.stderr)
        return 2
    all_errors: list[str] = []
    for p in paths:
        if not p.exists():
            all_errors.append(f"{p}: not found")
            continue
        all_errors.extend(check(p))
    if all_errors:
        for e in all_errors:
            print(e, file=sys.stderr)
        return 1
    print(f"OK: {len(paths)} file(s) parsed cleanly.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
