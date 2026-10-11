"""Workflow-level CP2K input tree editing helpers (refactor plan §3.3).

Second layer on top of the format primitives in
``md_analysis.utils.formats.cp2k.inp`` (parse / generate / path
walking): everything here knows about *workflow* semantics — which
sections a single-point or constrained-MD template must carry, which
keywords the generators own, and the default basis/potential assigned
to generated ``&KIND`` blocks.

Tree conventions (canonical lowercase keyword names, unit policy,
dict/list section forms) are those of ``inp.py``.  Callers follow the
double-track discipline: TIGen deep-copies the base tree per target;
SpGen / PotentialGen edit once and reuse the generated text (their
edits do not depend on the frame).
"""

from __future__ import annotations

import logging
from typing import Iterable

from ..utils.formats.cp2k.inp import (
    get_section,
    iter_collectives,
    set_keyword,
)
from ._potential_template import DEFAULT_BASIS_SET, DEFAULT_POTENTIAL

__all__ = [
    "ensure_topology_init_xyz",
    "set_project",
    "set_md_steps",
    "set_collective_target",
    "zero_target_growth",
    "ensure_kind_blocks",
]

logger = logging.getLogger(__name__)


def ensure_topology_init_xyz(tree: dict) -> None:
    """Point ``&TOPOLOGY`` at ``init.xyz`` (XYZ format), in place.

    Creates ``force_eval.subsys.topology`` when missing (``setdefault``
    single-instance section, logged); existing ``COORD_FILE_NAME`` /
    ``COORD_FILE_FORMAT`` values are overwritten either way.
    """
    subsys = get_section(tree, "force_eval.subsys")
    topology = subsys.setdefault("topology", {})
    if not topology:
        logger.info(
            "No &TOPOLOGY section in template; created one pointing at "
            "init.xyz"
        )
    topology["coord_file_name"] = "init.xyz"
    topology["coord_file_format"] = "xyz"


def set_project(tree: dict, project_name: str) -> None:
    """Write ``GLOBAL.PROJECT_NAME`` using the canonical key ``project_name``.

    A missing ``&GLOBAL`` section is created (single-instance section).
    """
    tree.setdefault("global", {})["project_name"] = project_name


def set_md_steps(tree: dict, steps: int) -> None:
    """Write ``MOTION.MD.STEPS`` (parent section must already exist)."""
    set_keyword(tree, "motion.md.steps", int(steps))


def set_collective_target(
    tree: dict,
    target_au: float,
    colvar_id: int | None = None,
) -> None:
    """Write a bare (internal-unit) ``TARGET`` into a ``&COLLECTIVE`` block.

    ``colvar_id=None`` means the **primary** (first) COLLECTIVE block —
    aligned with the reader-side ``colvars.primary`` semantics; an
    explicit id selects the matching ``COLVAR`` block.  ``KeyError`` is
    raised when there is no block at all or no block matches.
    """
    blocks = list(iter_collectives(tree))
    if not blocks:
        raise KeyError(
            "motion.constraint.collective: no &COLLECTIVE block found"
        )
    if colvar_id is None:
        blocks[0][1]["target"] = float(target_au)
        return
    for cid, block in blocks:
        if cid == colvar_id:
            block["target"] = float(target_au)
            return
    raise KeyError(
        f"motion.constraint.collective: no &COLLECTIVE block with "
        f"COLVAR {colvar_id}"
    )


def zero_target_growth(tree: dict) -> None:
    """Set ``TARGET_GROWTH`` to zero in **every** ``&COLLECTIVE`` block."""
    for _, block in iter_collectives(tree):
        block["target_growth"] = 0.0


def ensure_kind_blocks(
    tree: dict,
    elements: Iterable[str],
    *,
    basis_set: str = DEFAULT_BASIS_SET,
    potential: str = DEFAULT_POTENTIAL,
) -> None:
    """Upsert one ``&KIND`` block per element (first-appearance, deduped).

    Contract (plan §3.3) on ``force_eval.subsys.kind``:

    - ``kind`` key missing → created as a dict-of-dicts;
    - dict form → per-element upsert: an existing same-name block is
      replaced by the default ``{basis_set, potential}`` (idempotent);
    - list form (only manually built trees, defensive branch) → same
      upsert keyed on the ``_`` section parameter; missing elements are
      appended as ``{"_": symbol, ...}`` blocks.

    Default values are the project constants kept in
    ``_potential_template``; the parent path is located via
    ``inp.get_section``.
    """
    subsys = get_section(tree, "force_eval.subsys")
    kind = subsys.get("kind")
    if kind is None:
        kind = {}
        subsys["kind"] = kind
    for symbol in dict.fromkeys(elements):
        block = {"basis_set": basis_set, "potential": potential}
        if isinstance(kind, list):
            for i, existing in enumerate(kind):
                if isinstance(existing, dict) and existing.get("_") == symbol:
                    kind[i] = {"_": symbol, **block}
                    break
            else:
                kind.append({"_": symbol, **block})
        else:
            kind[symbol] = block
