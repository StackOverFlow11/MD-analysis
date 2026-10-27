"""Parse CP2K COLVAR restart files and LagrangeMultLog files.

The ``.restart`` path parses the file into a cit tree via
:func:`.inp.parse_inp` and reads ``global.project_name``,
``motion.md.{step_start_val, time_start_val, timestep, steps}`` (time
keywords land on their XML default unit, fs),
``motion.constraint.collective.{colvar, target, target_growth}``
(internal CP2K units, direct read -- the same convention as the restart
text itself), ``motion.constraint.lagrange_multipliers.filename`` and
``motion.constraint.fixed_atoms.list`` (integer ranges are expanded by
the inp.py preprocessing; the ascending-tuple contract is enforced
here).  Cell parameters come from :func:`.cell._abc_ang_from_cell_tree`
on the same tree, so the orthogonality policy is shared with
``cell.py``.  Legacy error messages are preserved verbatim for the
"parsed but the tree lacks the block" scenarios; the validation order
is pinned: PROJECT_NAME -> MD -> CONSTRAINT -> COLLECTIVE -> cell (the
lagrange-multipliers and fixed-atoms blocks are optional).  The
LagrangeMultLog parser is plain-text and independent of cit.

Output contract
---------------
Parsers return **CP2K-specific raw dataclasses** (``Cp2kConstraintInfoRaw`` /
``Cp2kColvarInfoRaw`` / ``Cp2kConstraintMetadataRaw`` / ``Cp2kLambdaSeriesRaw``).
The engine-neutral canonical types (``ConstraintInfo`` / ``ColvarInfo`` /
``ConstraintMetadata`` / ``LambdaSeries`` / ``ConstraintRun``) live in
:mod:`md_analysis.engines.models`.  Conversion from raw to canonical is
performed by ``engines.cp2k._cp2k_raw_to_*`` helpers.

Module-level invariant: this file does NOT runtime-import
``md_analysis.engines`` (R2).  Verified by ``tools/check_r2_imports.py``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .cell import _abc_ang_from_cell_tree
from .inp import get_keyword, get_section, iter_collectives, parse_inp

from ....exceptions import MDAnalysisError


class ColvarParseError(MDAnalysisError):
    """Raised when parsing a COLVAR restart or LagrangeMultLog file fails."""


# ---------------------------------------------------------------------------
# CP2K-specific raw dataclasses
#
# Field set is 1:1 with the canonical ConstraintInfo / ColvarInfo /
# ConstraintMetadata / LambdaSeries types in engines.models.  The raw types
# are intentionally accessor-free (no .primary / __getitem__ / accessor
# methods): they are the *parser output*, and engines.cp2k._cp2k_raw_to_*
# conversion helpers build the accessor-bearing canonical models.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Cp2kConstraintInfoRaw:
    """CP2K-specific raw view of one COLLECTIVE constraint block.

    ``target_growth_au`` is the rate of change per atomic unit of time
    (as stored in CP2K restart files).
    """

    colvar_id: int
    target_au: float
    target_growth_au: float
    intermolecular: bool


@dataclass(frozen=True)
class Cp2kColvarInfoRaw:
    """CP2K-specific raw view of the collection of COLLECTIVE constraints."""

    constraints: tuple[Cp2kConstraintInfoRaw, ...]


@dataclass(frozen=True)
class Cp2kConstraintMetadataRaw:
    """CP2K-specific raw view of constraint-MD restart metadata."""

    project_name: str
    step_start: int
    time_start_fs: float
    timestep_fs: float
    total_steps: int
    colvars: Cp2kColvarInfoRaw
    lagrange_filename: str | None
    cell_abc_ang: tuple[float, float, float]
    fixed_atom_indices: tuple[int, ...] | None


@dataclass(frozen=True)
class Cp2kLambdaSeriesRaw:
    """CP2K-specific raw view of a LagrangeMultLog file."""

    shake: np.ndarray
    rattle: np.ndarray
    n_steps: int
    n_constraints: int


# ---------------------------------------------------------------------------
# Private helpers — restart parsing (cit tree path)
# ---------------------------------------------------------------------------


def _require_tree_value(section: dict, tree_key: str, keyword: str, context: str):
    """Fetch *tree_key* from *section* or raise the legacy message."""
    try:
        return section[tree_key]
    except KeyError:
        raise ColvarParseError(f"{keyword} not found in {context}") from None


def _parse_md_from_tree(tree: dict) -> dict:
    """Read ``motion.md.{step_start_val, time_start_val, timestep, steps}``."""
    try:
        md = get_section(tree, "motion.md")
    except KeyError:
        raise ColvarParseError("No &MD block found") from None
    return {
        "step_start": int(
            _require_tree_value(md, "step_start_val", "STEP_START_VAL", "&MD")
        ),
        "time_start_fs": float(
            _require_tree_value(md, "time_start_val", "TIME_START_VAL", "&MD")
        ),
        "timestep_fs": float(
            _require_tree_value(md, "timestep", "TIMESTEP", "&MD")
        ),
        "total_steps": int(_require_tree_value(md, "steps", "STEPS", "&MD")),
    }


def _constraint_from_block(block: dict) -> Cp2kConstraintInfoRaw:
    """Build one raw constraint from a ``&COLLECTIVE`` block dict.

    ``target`` / ``target_growth`` are internal_cp2k keywords: the tree
    holds the internal-unit values (bohr / radian / aut) directly, like
    the ``.restart`` text itself -- no conversion happens here.
    """
    colvar = _require_tree_value(block, "colvar", "COLVAR", "&COLLECTIVE")
    if isinstance(colvar, (tuple, list)):  # multi-value keyword: first id,
        colvar = colvar[0]                 # as the legacy first-token read
    return Cp2kConstraintInfoRaw(
        colvar_id=int(colvar),
        target_au=float(
            _require_tree_value(block, "target", "TARGET", "&COLLECTIVE")
        ),
        target_growth_au=float(
            _require_tree_value(block, "target_growth", "TARGET_GROWTH", "&COLLECTIVE")
        ),
        intermolecular=block.get("intermolecular", False),
    )


def _parse_all_collective_blocks(tree: dict) -> Cp2kColvarInfoRaw:
    """Enumerate the ``&COLLECTIVE`` blocks of ``motion.constraint``.

    ``iter_collectives`` handles both tree forms (single dict or list of
    dicts for repeated ``&COLLECTIVE`` blocks).
    """
    try:
        get_section(tree, "motion.constraint")
    except KeyError:
        raise ColvarParseError("No &CONSTRAINT block found") from None
    try:
        get_section(tree, "motion.constraint.collective")
    except KeyError:
        raise ColvarParseError("No &COLLECTIVE block found in &CONSTRAINT") from None

    constraints = []
    try:
        for _, block in iter_collectives(tree):
            constraints.append(_constraint_from_block(block))
    except KeyError:
        # only remaining failure: a &COLLECTIVE block without a COLVAR id
        raise ColvarParseError("COLVAR not found in &COLLECTIVE") from None
    return Cp2kColvarInfoRaw(constraints=tuple(constraints))


def _parse_lagrange_filename(tree: dict) -> str | None:
    """Read ``motion.constraint.lagrange_multipliers.filename`` (or None)."""
    try:
        lm = get_section(tree, "motion.constraint.lagrange_multipliers")
    except KeyError:
        return None
    if isinstance(lm, list):  # param-less repeats: first block, like the regex
        lm = lm[0]
    filename = lm.get("filename")
    return filename if filename else None


def _parse_fixed_atoms_from_tree(tree: dict) -> tuple[int, ...] | None:
    """Read ``motion.constraint.fixed_atoms.list`` as a sorted int tuple.

    Integer ranges are already expanded by the inp.py preprocessing
    (W1); the explicit ``sorted`` keeps the ascending-tuple contract of
    the legacy helper.  A ``&FIXED_ATOMS`` block without a ``LIST``
    keyword yields ``()`` (legacy semantics), a missing block ``None``.
    """
    try:
        fa = get_section(tree, "motion.constraint.fixed_atoms")
    except KeyError:
        return None
    if isinstance(fa, list):  # param-less repeats: first block, like the regex
        fa = fa[0]
    return tuple(sorted(int(i) for i in fa.get("list", ())))


# ---------------------------------------------------------------------------
# Private helpers — LagrangeMultLog parsing
# ---------------------------------------------------------------------------

_LABEL_RE = re.compile(r"^(Shake|Rattle)\s+Lagrangian\s+Multipliers:", re.IGNORECASE)
_OVERFLOW_RE = re.compile(r"^\*+$")


def _safe_float(token: str) -> float:
    """Convert a token to float, returning nan for CP2K overflow (``***``)."""
    if _OVERFLOW_RE.match(token.strip()):
        return float("nan")
    return float(token)


def _detect_log_format(lines: list[str]) -> str:
    if len(lines) < 2:
        raise ColvarParseError("LagrangeMultLog file is too short")
    if _LABEL_RE.match(lines[1].strip()):
        return "single"
    return "multi"


def _parse_single_constraint_log(
    lines: list[str],
) -> tuple[np.ndarray, np.ndarray, int]:
    shake_vals: list[float] = []
    rattle_vals: list[float] = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("Shake"):
            shake_vals.append(_safe_float(stripped.split(":")[1].strip()))
        elif stripped.startswith("Rattle"):
            rattle_vals.append(_safe_float(stripped.split(":")[1].strip()))
    n = min(len(shake_vals), len(rattle_vals))
    return np.array(shake_vals[:n]), np.array(rattle_vals[:n]), n


def _parse_multi_constraint_log(
    lines: list[str],
) -> tuple[np.ndarray, np.ndarray, int, int]:
    blocks: list[tuple[str, list[float]]] = []
    current_type: str | None = None
    current_values: list[float] = []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        label_m = _LABEL_RE.match(stripped)
        if label_m:
            if current_type is not None:
                blocks.append((current_type, current_values))
            current_type = label_m.group(1).lower()
            after_colon = stripped.split(":", 1)[1]
            current_values = [_safe_float(x) for x in after_colon.split()]
        else:
            current_values.extend(_safe_float(x) for x in stripped.split())

    if current_type is not None:
        blocks.append((current_type, current_values))

    # Pair consecutive Shake/Rattle blocks
    shake_list: list[list[float]] = []
    rattle_list: list[list[float]] = []
    i = 0
    while i < len(blocks) - 1:
        if blocks[i][0] == "shake" and blocks[i + 1][0] == "rattle":
            shake_list.append(blocks[i][1])
            rattle_list.append(blocks[i + 1][1])
            i += 2
        else:
            i += 1  # skip unpaired

    n_steps = len(shake_list)
    n_constraints = len(shake_list[0]) if n_steps > 0 else 0
    return (
        np.array(shake_list),
        np.array(rattle_list),
        n_steps,
        n_constraints,
    )


# ---------------------------------------------------------------------------
# Public API — parsers return RAW types (Phase 4 Commit 1)
# ---------------------------------------------------------------------------


def parse_colvar_restart(restart_path: str | Path) -> Cp2kConstraintMetadataRaw:
    """Parse COLVAR metadata from a CP2K restart file into a raw type.

    Returns a CP2K-specific raw dataclass; conversion to the canonical
    :class:`md_analysis.engines.models.ConstraintMetadata` is performed by
    ``md_analysis.engines.cp2k._cp2k_raw_to_constraint_metadata``.

    The file is parsed once into a cit tree (see :mod:`.inp`); the cell
    parameters are extracted from the same tree via
    :func:`.cell._abc_ang_from_cell_tree`.  Missing blocks raise
    :class:`ColvarParseError` with the legacy messages, in the pinned
    order PROJECT_NAME -> MD -> CONSTRAINT -> COLLECTIVE -> cell (the
    ``&LAGRANGE_MULTIPLIERS`` and ``&FIXED_ATOMS`` blocks are optional;
    cell errors surface as :class:`.cell.CellParseError`).
    """
    tree = parse_inp(restart_path)

    try:
        project_name = get_keyword(tree, "global.project_name")
    except KeyError:
        raise ColvarParseError(f"PROJECT_NAME not found in {restart_path}") from None

    md = _parse_md_from_tree(tree)
    colvars = _parse_all_collective_blocks(tree)
    cell_abc_ang = _abc_ang_from_cell_tree(tree, restart_path)
    lagrange_filename = _parse_lagrange_filename(tree)
    fixed_atoms = _parse_fixed_atoms_from_tree(tree)

    return Cp2kConstraintMetadataRaw(
        project_name=project_name,
        step_start=md["step_start"],
        time_start_fs=md["time_start_fs"],
        timestep_fs=md["timestep_fs"],
        total_steps=md["total_steps"],
        colvars=colvars,
        lagrange_filename=lagrange_filename,
        cell_abc_ang=cell_abc_ang,
        fixed_atom_indices=fixed_atoms,
    )


def parse_lagrange_mult_log(log_path: str | Path) -> Cp2kLambdaSeriesRaw:
    """Parse a LagrangeMultLog file into a raw type.

    Returns a CP2K-specific raw dataclass; conversion to the canonical
    :class:`md_analysis.engines.models.LambdaSeries` is performed by
    ``md_analysis.engines.cp2k._cp2k_raw_to_lambda_series``.

    Auto-detects single/multi constraint log format.
    """
    path = Path(log_path)
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ColvarParseError(f"LagrangeMultLog file is empty: {log_path}")

    fmt = _detect_log_format(lines)

    if fmt == "single":
        shake, rattle, n_steps = _parse_single_constraint_log(lines)
        return Cp2kLambdaSeriesRaw(
            shake=shake, rattle=rattle, n_steps=n_steps, n_constraints=1,
        )

    shake, rattle, n_steps, n_constraints = _parse_multi_constraint_log(lines)
    return Cp2kLambdaSeriesRaw(
        shake=shake, rattle=rattle, n_steps=n_steps, n_constraints=n_constraints,
    )
