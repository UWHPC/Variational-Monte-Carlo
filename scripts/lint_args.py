#!/usr/bin/env python3
"""
lint_args.py - Enforce block-indent style for functions with 2+ arguments.

Flags function declarations, definitions, and calls where 2+ arguments
sit on a single line instead of being wrapped one-per-line.

Usage:
    python scripts/lint_args.py                   # check src/ and tests/
    python scripts/lint_args.py src/pbc/pbc.hpp   # check specific files
    python scripts/lint_args.py --fix src/         # auto-fix in place

Ignores: comments, preprocessor, for/if/while, return, throw, static_assert,
         template<>, braced-init lists, and single-argument functions.
"""

import re
import sys
from pathlib import Path

INDENT = "    "
EXTENSIONS = {".hpp", ".cpp", ".h", ".cc", ".cxx"}
WRAP_THRESHOLD = 100  # only wrap if the line exceeds this many characters

# Lines to never touch
SKIP_RE = re.compile(
    r"^\s*(?:"
    r"//|/?\*"
    r"|#"
    r"|template\s*<"
    r"|for\s*\("
    r"|if\s*\("
    r"|while\s*\("
    r"|return\b"
    r"|throw\b"
    r"|static_assert"
    r"|using\b"
    r")"
)


def count_top_level_commas(s):
    """Count commas not nested inside <>, (), {}, or []."""
    depth = 0
    count = 0
    for ch in s:
        if ch in "(<[{":
            depth += 1
        elif ch in ")>]}":
            depth = max(0, depth - 1)
        elif ch == "," and depth == 0:
            count += 1
    return count


def split_top_level(s):
    """Split string by top-level commas."""
    parts = []
    current = []
    depth = 0
    for ch in s:
        if ch in "(<[{":
            depth += 1
            current.append(ch)
        elif ch in ")>]}":
            depth = max(0, depth - 1)
            current.append(ch)
        elif ch == "," and depth == 0:
            parts.append("".join(current).strip())
            current = []
        else:
            current.append(ch)
    parts.append("".join(current).strip())
    return [p for p in parts if p]


def find_func_paren(line):
    """
    Find the outermost '(' that looks like a function call/decl.
    Returns (open_col, close_col) or None.
    """
    opens = []
    pairs = []
    in_str = False
    in_char = False
    for i, ch in enumerate(line):
        if ch == '"' and not in_char and (i == 0 or line[i - 1] != '\\'):
            in_str = not in_str
        if ch == "'" and not in_str and (i == 0 or line[i - 1] != '\\'):
            in_char = not in_char
        if in_str or in_char:
            continue
        if ch == '(':
            opens.append(i)
        elif ch == ')' and opens:
            o = opens.pop()
            pairs.append((o, i))

    if not pairs:
        return None

    pairs.sort(key=lambda p: p[0])
    open_col, close_col = pairs[0]

    prefix = line[:open_col].rstrip()
    if not prefix:
        return None
    if not re.search(r"[\w>]$", prefix):
        return None
    if prefix.endswith(("=", "{", ",")):
        return None

    return (open_col, close_col)


def is_declaration_or_definition(line, open_col):
    """
    Heuristic: a function declaration or definition has a return type
    followed by a function name before the '('. A call is just a name.

    We also check the suffix after ')' — declarations end with ';' or '{',
    often preceded by qualifiers like const/noexcept/override.
    Calls end with just ';' or ');' without qualifiers.
    """
    prefix = line[:open_col].strip()
    suffix_start = line.find(")", open_col)
    suffix = line[suffix_start + 1:].strip() if suffix_start >= 0 else ""

    # Check suffix: declarations have qualifiers or '{' after ')'
    # e.g. ") const noexcept {", ") noexcept;", ") override;"
    # Calls just have ";" or nothing meaningful
    has_decl_suffix = bool(re.match(
        r"^(?:const|volatile|noexcept|override|final|\s)*[{;]",
        suffix
    )) and bool(re.search(r"(?:const|volatile|noexcept|override|final|\{)", suffix))

    # Check prefix: must have a return type token before the function name
    # Pattern: [qualifiers] <type> [*&] <name>(
    # e.g. "void foo", "double* bar", "[[nodiscard]] int baz"
    #       "explicit Foo", "static void qux"
    decl_prefix = re.match(
        r"^"
        r"(?:\[\[[\w:]+\]\]\s+)?"                       # optional [[attribute]]
        r"(?:(?:static|virtual|inline|explicit|friend|constexpr|consteval)\s+)*"  # optional qualifiers
        r"(?:void|bool|char|int|float|double|auto|unsigned|signed|long|short"
        r"|std::size_t|std::string|std::unique_ptr"       # common std types
        r"|\w+)"                                          # any type name
        r"[\s\*&]+"                                       # space/pointer/ref separator
        r"[\w:~]+$",                                      # function name (possibly Class::method)
        prefix
    )

    # Constructor/destructor pattern: ClassName::ClassName( or ClassName::~ClassName(
    ctor_pattern = re.match(r"^(\w+)::~?\1$", prefix)

    # For class member declarations (inside class body), prefix is just type + name
    # but we need the suffix check too
    if decl_prefix:
        return True
    if ctor_pattern:
        return True
    # Fallback: if suffix has qualifiers, it's likely a declaration
    if has_decl_suffix:
        return True

    return False


def check_line(line, lineno, filepath):
    """Check a single line. Returns a violation dict or None."""
    if SKIP_RE.match(line):
        return None

    result = find_func_paren(line)
    if result is None:
        return None

    open_col, close_col = result

    # Only wrap declarations and definitions, not calls inside bodies
    if not is_declaration_or_definition(line, open_col):
        return None

    content = line[open_col + 1 : close_col]

    if count_top_level_commas(content) < 1:
        return None

    # Only wrap if the line is long enough to justify it
    if len(line.rstrip()) <= WRAP_THRESHOLD:
        return None

    return {
        "file": filepath,
        "line": lineno,
        "num_args": count_top_level_commas(content) + 1,
        "text": line.rstrip(),
        "open_col": open_col,
        "close_col": close_col,
        "content": content,
    }


def fix_line(line, info):
    """Rewrite a single-line function signature into block-indent style."""
    open_col = info["open_col"]
    close_col = info["close_col"]

    before_paren = line[:open_col]
    after_paren = line[close_col + 1 :].rstrip()

    base_indent = ""
    for ch in line:
        if ch in (" ", "\t"):
            base_indent += ch
        else:
            break

    args = split_top_level(info["content"])
    lines = [before_paren.rstrip() + "("]
    for i, arg in enumerate(args):
        suffix = "," if i < len(args) - 1 else ""
        lines.append(f"{base_indent}{INDENT}{arg}{suffix}")
    lines.append(f"{base_indent}){after_paren}")

    return "\n".join(lines)


CTOR_BODY_RE = re.compile(r"^(.*\S)\s*\{\s*\}\s*$")


def fix_constructor_bodies(text):
    """
    Fix lines like:
        , invL_{1.0 / L} { }
    Into:
        , invL_{1.0 / L}
        { }
    where { } is indented to match the constructor declaration line,
    not the initializer list lines.
    """
    lines = text.splitlines()
    new_lines = []
    for idx, line in enumerate(lines):
        m = CTOR_BODY_RE.match(line)
        if m:
            prefix = m.group(1)
            stripped = prefix.lstrip()
            # Check if this looks like a constructor initializer line
            if stripped.startswith(",") or stripped.startswith(":"):
                # Walk backwards to find the constructor declaration line
                # (the line before the first ':' initializer)
                ctor_indent = ""
                for back in range(len(new_lines) - 1, -1, -1):
                    back_stripped = new_lines[back].lstrip()
                    if back_stripped.startswith(":") or back_stripped.startswith(","):
                        continue
                    # Found the constructor line
                    ctor_indent = new_lines[back][: len(new_lines[back]) - len(back_stripped)]
                    break
                new_lines.append(prefix)
                new_lines.append(f"{ctor_indent}{{ }}")
                continue
        new_lines.append(line)
    return "\n".join(new_lines)


def fix_closing_paren_indent(text):
    """
    Fix closing parens on declarations/definitions to align with the
    function declaration line (the line containing the opening '(').

    Transforms:
        void foo(
            int a,
            int b
        ) const noexcept {

    Into:
        void foo(
            int a,
            int b
            ) const noexcept {

    The ')' indentation = declaration line indent + one INDENT level.
    """
    lines = text.splitlines()
    new_lines = []
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.lstrip()
        # Look for a line that starts with ')'
        if stripped.startswith(")") and i > 0:
            # Walk backwards to find the line with the opening '('
            decl_indent = ""
            found_decl = False
            for back in range(i - 1, -1, -1):
                back_line = new_lines[back] if back < len(new_lines) else lines[back]
                back_stripped = back_line.rstrip()
                if back_stripped.endswith("("):
                    # Check if the opening line looks like a declaration/definition
                    prefix = back_stripped[:-1].strip()
                    # Declarations have return type + name, or qualifiers
                    # e.g. "void evaluateLogPsi", "void Foo::bar"
                    # Calls are just "functionName" with body indent
                    looks_like_decl = bool(re.search(
                        r"(?:void|bool|char|int|float|double|auto|explicit|virtual|static|inline"
                        r"|constexpr|unsigned|signed|long|short|\w+::\w+|std::\w+)"
                        r"\s+[\w:~]+$",
                        prefix
                    )) or bool(re.search(r"\[\[[\w:]+\]\]", prefix))
                    if looks_like_decl:
                        decl_indent = back_line[: len(back_line) - len(back_line.lstrip())]
                        found_decl = True
                    break
            if found_decl:
                target_indent = decl_indent
                curr_indent = len(line) - len(stripped)
                if curr_indent != len(target_indent):
                    new_lines.append(target_indent + stripped)
                    i += 1
                    continue
        new_lines.append(line)
        i += 1
    return "\n".join(new_lines)


def process_file(filepath, fix=False):
    """Process one file. Returns list of violation messages."""
    text = filepath.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    violations = []
    new_lines = []

    for i, line in enumerate(lines):
        info = check_line(line, i + 1, filepath)
        if info is None:
            new_lines.append(line)
            continue

        violations.append(
            f"{filepath.as_posix()}:{i + 1}: "
            f"{info['num_args']} args on one line: {line.strip()}"
        )

        if fix:
            new_lines.append(fix_line(line, info))
        else:
            new_lines.append(line)

    if fix and violations:
        filepath.write_text("\n".join(new_lines) + "\n", encoding="utf-8")

    # Also fix constructor empty bodies and closing paren indentation
    if fix:
        text = filepath.read_text(encoding="utf-8", errors="replace")
        text = fix_constructor_bodies(text)
        text = fix_closing_paren_indent(text)
        filepath.write_text(text + "\n" if not text.endswith("\n") else text, encoding="utf-8")

    return violations


def collect_files(paths):
    """Gather all C++ source files."""
    files = []
    for p in paths:
        path = Path(p)
        if path.is_file() and path.suffix in EXTENSIONS:
            files.append(path)
        elif path.is_dir():
            for ext in EXTENSIONS:
                files.extend(path.rglob(f"*{ext}"))
    return sorted(set(files))


def main():
    args = sys.argv[1:]
    fix = "--fix" in args
    if fix:
        args.remove("--fix")
    if not args:
        args = ["src/", "tests/"]

    files = collect_files(args)
    if not files:
        print("No C++ files found.")
        return 0

    total = []
    for f in files:
        total.extend(process_file(f, fix=fix))

    if total:
        label = "Fixed" if fix else "Violation"
        for v in total:
            print(f"  {label}: {v}")
        print(f"\n{len(total)} issue(s) across {len(files)} file(s).")
        if not fix:
            print("Run with --fix to auto-fix.")
        return 1

    print(f"All clean: {len(files)} file(s) checked.")
    return 0


if __name__ == "__main__":
    sys.exit(main())