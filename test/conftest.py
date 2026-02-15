"""Shared pytest fixtures and helpers for MD Analysis tests."""

from __future__ import annotations

import re
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def parse_abc_from_md_inp(md_inp_path: Path) -> tuple[float, float, float]:
    """Parse `ABC [angstrom] a b c` from CP2K md.inp."""
    text = md_inp_path.read_text(encoding="utf-8")
    match = re.search(r"^\s*ABC\s+\[angstrom\]\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*$", text, re.MULTILINE)
    if not match:
        raise ValueError(f"Cannot find `ABC [angstrom] ...` in {md_inp_path}")
    a, b, c = (float(match.group(i)) for i in (1, 2, 3))
    return a, b, c
