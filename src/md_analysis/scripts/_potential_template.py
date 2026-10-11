"""Built-in default CP2K potential-SP inp template and template inspection.

Design choice — *template file + element-driven ``&KIND`` block generation*
rather than pure placeholder substitution inside a static template:

- The repo paradigm (``BaderGen``) ships static template files under
  ``md_analysis/scripts/template/`` accessed via ``importlib.resources``, and
  ``_inp_utils.py`` patches inp text with plain string operations.  A template
  file keeps the physics reviewable/diffable against the validated original.
- The element composition of a frame is only known at call time, and the
  element set is open-ended (Pt/O/H today, Cu/Ag/K/Cl tomorrow).  Per-element
  placeholders (``{{KIND_Pt}}`` …) cannot express an open element set; a
  generator function can.  Everything else in the inp stays static text, so
  the whole pipeline remains plain-string based and unit-testable.

Physics baseline for the shipped template
``template/sp_potential.inp``: the validated sp.inp from
``data_example/potential/distributed/`` (PBE + DFTD3, GTH-PBE pseudopotentials,
DZVP-MOLOPT-SR-GTH basis, metallic smearing, ``&V_HARTREE_CUBE``).  The only
deliberate divergence: ``&KIND`` blocks are generated with the bare names
``POTENTIAL GTH-PBE`` / ``BASIS_SET DZVP-MOLOPT-SR-GTH`` instead of the
q-suffixed ones used in the validated file (``GTH-PBE-q6`` …).  CP2K resolves
the bare potential name through the alias/"default GTH PP" mechanism in
``GTH_POTENTIALS`` (H→q1, O→q6, C→q4, Cu→q11, Ag→q11, Pt→q18, …), and
``BASIS_MOLOPT`` ships a valence-matching DZVP-MOLOPT-SR-GTH for each of
those.  This avoids hard-coding an element→valence table in md-analysis; the
q-suffix mapping stays CP2K's responsibility.  (The validated file also set
``MASS 2.0`` on H for a D2O system — system-specific and therefore not
reproduced by the generic renderer.)
"""

from __future__ import annotations

import re
from importlib.resources import files
from typing import Iterable

__all__ = [
    "DEFAULT_BASIS_SET",
    "DEFAULT_POTENTIAL",
    "DEFAULT_TEMPLATE_RESOURCE",
    "KIND_BLOCKS_PLACEHOLDER",
    "extract_kind_elements",
    "has_smear",
    "has_v_hartree_cube",
    "load_default_template_text",
    "render_default_potential_inp",
    "render_kind_blocks",
]

#: File name of the built-in template inside ``md_analysis.scripts.template``.
DEFAULT_TEMPLATE_RESOURCE = "sp_potential.inp"

#: Line placeholder inside the built-in template where &KIND blocks are inserted.
KIND_BLOCKS_PLACEHOLDER = "{{KIND_BLOCKS}}"

#: Basis / potential assigned to every element by :func:`render_kind_blocks`.
DEFAULT_BASIS_SET = "DZVP-MOLOPT-SR-GTH"
DEFAULT_POTENTIAL = "GTH-PBE"

# CP2K keywords are case-insensitive; "&" may be separated from the name by
# whitespace.  Matching is deliberately kept simple (no full inp parser).
_KIND_BLOCK_RE = re.compile(r"&KIND\b.*?&END\s+KIND", re.DOTALL | re.IGNORECASE)
_KIND_NAME_RE = re.compile(r"&KIND\s+(\S+)", re.IGNORECASE)
_KIND_ELEMENT_RE = re.compile(r"^\s*ELEMENT\s+(\S+)", re.MULTILINE | re.IGNORECASE)
_V_HARTREE_CUBE_RE = re.compile(r"&\s*V_HARTREE_CUBE\b", re.IGNORECASE)
_SMEAR_RE = re.compile(r"&\s*SMEAR\b", re.IGNORECASE)


def load_default_template_text() -> str:
    """Return the raw built-in template text (``{{KIND_BLOCKS}}`` unresolved)."""
    resource = files("md_analysis.scripts.template").joinpath(DEFAULT_TEMPLATE_RESOURCE)
    return resource.read_text(encoding="utf-8")


def render_kind_blocks(elements: Iterable[str]) -> str:
    """Generate one ``&KIND`` block per element (first-appearance order, deduped).

    Each block assigns the project defaults ``BASIS_SET DZVP-MOLOPT-SR-GTH`` +
    ``POTENTIAL GTH-PBE``; see the module docstring for how CP2K resolves the
    bare names per element.
    """
    lines: list[str] = []
    for symbol in dict.fromkeys(elements):
        lines.append(f"    &KIND {symbol}")
        lines.append(f"      BASIS_SET {DEFAULT_BASIS_SET}")
        lines.append(f"      POTENTIAL {DEFAULT_POTENTIAL}")
        lines.append("    &END KIND")
    return "\n".join(lines)


def render_default_potential_inp(elements: Iterable[str]) -> str:
    """Return the built-in template with ``&KIND`` blocks rendered for *elements*."""
    template = load_default_template_text()
    if KIND_BLOCKS_PLACEHOLDER not in template:
        raise ValueError(
            f"Built-in template '{DEFAULT_TEMPLATE_RESOURCE}' is missing the "
            f"'{KIND_BLOCKS_PLACEHOLDER}' placeholder; the package template "
            "file is corrupted."
        )
    return template.replace(KIND_BLOCKS_PLACEHOLDER, render_kind_blocks(elements))


def extract_kind_elements(inp_text: str) -> set[str]:
    """Return element symbols (upper-case) covered by ``&KIND`` blocks.

    A kind covers an element when the kind *name* is the element symbol or
    when the block carries an explicit ``ELEMENT <symbol>`` keyword (kinds may
    be named freely, e.g. ``&KIND O_wat ELEMENT O``).
    """
    covered: set[str] = set()
    for match in _KIND_BLOCK_RE.finditer(inp_text):
        block = match.group(0)
        name = _KIND_NAME_RE.search(block)
        if name is not None:
            covered.add(name.group(1).upper())
        element = _KIND_ELEMENT_RE.search(block)
        if element is not None:
            covered.add(element.group(1).upper())
    return covered


def has_v_hartree_cube(inp_text: str) -> bool:
    """True if the inp contains a ``&V_HARTREE_CUBE`` print section (case-insensitive)."""
    return _V_HARTREE_CUBE_RE.search(inp_text) is not None


def has_smear(inp_text: str) -> bool:
    """True if the inp contains a ``&SMEAR`` section (case-insensitive)."""
    return _SMEAR_RE.search(inp_text) is not None
