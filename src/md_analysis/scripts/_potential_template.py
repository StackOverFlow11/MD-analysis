"""Built-in default CP2K potential-SP template: constants and tree checks.

Design — *legal template file + tree-driven ``&KIND`` insertion*: the
shipped ``template/sp_potential.inp`` is plain, parseable CP2K input
(no placeholders).  At generation time it is parsed with the cit
wrapper (``utils.formats.cp2k.inp``), ``&KIND`` blocks are upserted
from the frame's element composition
(``scripts/_inp_tree_edit.py::ensure_kind_blocks``) and the edited tree
is generated back to text.  The element composition of a frame is only
known at call time and the element set is open-ended (Pt/O/H today,
Cu/Ag/K/Cl tomorrow); per-element placeholders cannot express an open
element set, a tree edit can.

Physics baseline for the shipped template
``template/sp_potential.inp``: the validated sp.inp from
``data_example/potential/distributed/`` (PBE + DFTD3, GTH-PBE pseudopotentials,
DZVP-MOLOPT-SR-GTH basis, metallic smearing, ``&V_HARTREE_CUBE``).  The only
deliberate divergence: generated ``&KIND`` blocks use the bare names
``POTENTIAL GTH-PBE`` / ``BASIS_SET DZVP-MOLOPT-SR-GTH`` instead of the
q-suffixed ones used in the validated file (``GTH-PBE-q6`` …).  CP2K resolves
the bare potential name through the alias/"default GTH PP" mechanism in
``GTH_POTENTIALS`` (H→q1, O→q6, C→q4, Cu→q11, Ag→q11, Pt→q18, …), and
``BASIS_MOLOPT`` ships a valence-matching DZVP-MOLOPT-SR-GTH for each of
those.  This avoids hard-coding an element→valence table in md-analysis; the
q-suffix mapping stays CP2K's responsibility.  (The validated file also set
``MASS 2.0`` on H for a D2O system — system-specific and therefore not
reproduced by the generated defaults.)

The checks in this module operate on **cit trees** (see
``utils/formats/cp2k/inp.py`` for the tree conventions), not on inp
text.
"""

from __future__ import annotations

from ..utils.formats.cp2k.inp import get_section

__all__ = [
    "DEFAULT_BASIS_SET",
    "DEFAULT_POTENTIAL",
    "DEFAULT_TEMPLATE_RESOURCE",
    "extract_kind_elements",
    "has_smear",
    "has_v_hartree_cube",
]

#: File name of the built-in template inside ``md_analysis.scripts.template``.
DEFAULT_TEMPLATE_RESOURCE = "sp_potential.inp"

#: Basis / potential assigned to every element by ``ensure_kind_blocks``.
DEFAULT_BASIS_SET = "DZVP-MOLOPT-SR-GTH"
DEFAULT_POTENTIAL = "GTH-PBE"


def _section_exists(tree: dict, path: str) -> bool:
    """True when the section at dotted *path* exists in the tree."""
    try:
        get_section(tree, path)
    except KeyError:
        return False
    return True


def extract_kind_elements(tree: dict) -> set[str]:
    """Return element symbols (upper-case) covered by ``&KIND`` blocks.

    A kind covers an element when the kind *name* is the element symbol
    or when the block carries an explicit ``ELEMENT <symbol>`` keyword
    (kinds may be named freely, e.g. ``&KIND O_wat ELEMENT O``).
    Handles both tree forms of ``force_eval.subsys.kind``: dict-of-dicts
    (parsed form, kind name = dict key) and list-of-dicts (``_``
    section parameter).
    """
    try:
        kinds = get_section(tree, "force_eval.subsys.kind")
    except KeyError:
        return set()

    named: list[tuple[str | None, dict]] = []
    if isinstance(kinds, dict):
        named = [
            (name, block)
            for name, block in kinds.items()
            if isinstance(block, dict)
        ]
    elif isinstance(kinds, list):
        named = [
            (block.get("_"), block)
            for block in kinds
            if isinstance(block, dict)
        ]

    covered: set[str] = set()
    for name, block in named:
        if isinstance(name, str) and name:
            covered.add(name.upper())
        element = block.get("element")
        if isinstance(element, str) and element:
            covered.add(element.upper())
    return covered


def has_smear(tree: dict) -> bool:
    """True if the tree contains a ``&SMEAR`` section (``&SCF`` level)."""
    return _section_exists(tree, "force_eval.dft.scf.smear")


def has_v_hartree_cube(tree: dict) -> bool:
    """True if the tree contains a ``&V_HARTREE_CUBE`` print section."""
    return _section_exists(tree, "force_eval.dft.print.v_hartree_cube")
