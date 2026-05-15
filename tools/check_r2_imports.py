#!/usr/bin/env python3
"""R2 invariant scanner: src/md_analysis/utils MUST NOT runtime-import
md_analysis.engines.

Detects all 5 forms (per Phase 3 design doc + Round 11 review):
  1. `import md_analysis.engines[.X]`           (Import + alias.name)
  2. `from md_analysis.engines[.X] import Y`    (ImportFrom mod prefix)
  3a. `from md_analysis import engines`         (ImportFrom mod + alias)
  3b. `from .... import engines`                (ImportFrom level>=1 + alias)
  4. `from ....engines[.X] import Y`            (ImportFrom level>=1 + relative)
  5. `importlib.import_module("md_analysis.engines...")` (dynamic str)

Allows engines imports only when sitting inside an `if TYPE_CHECKING:`
body. Exit 0 if clean; exit 1 with violations; exit 2 if utils dir
not found.
"""
from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

DYNAMIC_IMPORT_RE = re.compile(
    r"""importlib\.import_module\s*\(\s*["']md_analysis\.engines"""
)


def _import_node_touches_engines(node: ast.AST) -> bool:
    """True iff *node* is an Import/ImportFrom touching md_analysis.engines."""
    if isinstance(node, ast.Import):
        return any(
            alias.name == "md_analysis.engines"
            or alias.name.startswith("md_analysis.engines.")
            for alias in node.names
        )
    if isinstance(node, ast.ImportFrom):
        mod = node.module or ""

        # Form 2: from md_analysis.engines[.X] import Y
        if mod == "md_analysis.engines" or mod.startswith("md_analysis.engines."):
            return True

        # Form 4: relative from ....engines[.X] import Y
        if node.level >= 1:
            if mod == "engines" or mod.startswith("engines."):
                return True

        # Form 3a (absolute): from md_analysis import engines
        if mod == "md_analysis":
            for alias in node.names:
                if alias.name == "engines":
                    return True

        # Form 3b (relative): from .... import engines  (level>=1, mod="")
        if node.level >= 1 and mod == "":
            for alias in node.names:
                if alias.name == "engines":
                    return True

        return False
    return False


def _type_checking_member_ids(tree: ast.AST) -> set[int]:
    """id(node) of every AST node inside any `if TYPE_CHECKING:` body."""
    ids: set[int] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.If):
            try:
                if ast.unparse(node.test) == "TYPE_CHECKING":
                    for sub in ast.walk(node):
                        ids.add(id(sub))
            except Exception:
                pass
    return ids


def scan_file(path: Path) -> list[tuple[Path, int, str]]:
    src = path.read_text(encoding="utf-8")
    violations: list[tuple[Path, int, str]] = []

    # AST-based: forms 1-4
    try:
        tree = ast.parse(src)
    except SyntaxError as exc:
        return [(path, exc.lineno or 0, f"SyntaxError: {exc}")]
    tc_ids = _type_checking_member_ids(tree)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            if not _import_node_touches_engines(node):
                continue
            if id(node) in tc_ids:
                continue
            try:
                text = ast.unparse(node)
            except Exception:
                text = "<unparse failed>"
            violations.append((path, node.lineno, text))

    # String-based: form 5 (dynamic import) — does NOT honour TYPE_CHECKING
    # because importlib.import_module IS a runtime call by definition
    for m in DYNAMIC_IMPORT_RE.finditer(src):
        lineno = src[: m.start()].count("\n") + 1
        violations.append((path, lineno, f"dynamic: {m.group()}"))

    return violations


def main() -> int:
    root = Path("src/md_analysis/utils")
    if not root.is_dir():
        print(f"ERROR: {root} not found", file=sys.stderr)
        return 2
    all_v: list[tuple[Path, int, str]] = []
    for py in root.rglob("*.py"):
        all_v.extend(scan_file(py))
    if all_v:
        print("R2 VIOLATIONS (utils must not runtime-import engines):")
        for path, lineno, text in all_v:
            print(f"  {path}:{lineno}  {text}")
        return 1
    print("R2 OK: no runtime engines import found under src/md_analysis/utils")
    return 0


if __name__ == "__main__":
    sys.exit(main())
