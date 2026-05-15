"""Parse CP2K COLVAR restart files and LagrangeMultLog files.

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

from .cell import parse_abc_from_restart

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
# Private helpers — restart parsing
# ---------------------------------------------------------------------------

_CONSTRAINT_BLOCK_RE = re.compile(
    r"&CONSTRAINT\b(.*?)&END\s+CONSTRAINT", re.DOTALL | re.IGNORECASE,
)
_COLLECTIVE_BLOCK_RE = re.compile(
    r"&COLLECTIVE\s*\n(.*?)&END\s+COLLECTIVE", re.DOTALL | re.IGNORECASE,
)
_LAGRANGE_BLOCK_RE = re.compile(
    r"&LAGRANGE_MULTIPLIERS\b[^\n]*\n(.*?)&END\s+LAGRANGE_MULTIPLIERS",
    re.DOTALL | re.IGNORECASE,
)
_FIXED_ATOMS_BLOCK_RE = re.compile(
    r"&FIXED_ATOMS\s*\n(.*?)&END\s+FIXED_ATOMS",
    re.DOTALL | re.IGNORECASE,
)


def _extract_scalar(block: str, key: str) -> str | None:
    m = re.search(rf"^\s*{key}\s+(\S+)", block, re.MULTILINE | re.IGNORECASE)
    return m.group(1) if m else None


def _require_scalar(block: str, key: str, context: str) -> str:
    val = _extract_scalar(block, key)
    if val is None:
        raise ColvarParseError(f"{key} not found in {context}")
    return val


def _parse_md_block(text: str) -> dict:
    md_match = re.search(r"&MD\b(.*?)&END\s+MD", text, re.DOTALL | re.IGNORECASE)
    if not md_match:
        raise ColvarParseError("No &MD block found")
    block = md_match.group(1)
    return {
        "step_start": int(_require_scalar(block, "STEP_START_VAL", "&MD")),
        "time_start_fs": float(_require_scalar(block, "TIME_START_VAL", "&MD")),
        "timestep_fs": float(_require_scalar(block, "TIMESTEP", "&MD")),
        "total_steps": int(_require_scalar(block, "STEPS", "&MD")),
    }


def _parse_single_collective_block(block: str) -> Cp2kConstraintInfoRaw:
    inter_val = _extract_scalar(block, "INTERMOLECULAR")
    intermolecular = False
    if inter_val is not None:
        intermolecular = inter_val.upper().strip(".") in ("T", "TRUE", "YES")

    return Cp2kConstraintInfoRaw(
        colvar_id=int(_require_scalar(block, "COLVAR", "&COLLECTIVE")),
        target_au=float(_require_scalar(block, "TARGET", "&COLLECTIVE")),
        target_growth_au=float(
            _require_scalar(block, "TARGET_GROWTH", "&COLLECTIVE")
        ),
        intermolecular=intermolecular,
    )


def _parse_all_collective_blocks(text: str) -> Cp2kColvarInfoRaw:
    constraint_match = _CONSTRAINT_BLOCK_RE.search(text)
    if not constraint_match:
        raise ColvarParseError("No &CONSTRAINT block found")
    constraint_text = constraint_match.group(1)

    matches = list(_COLLECTIVE_BLOCK_RE.finditer(constraint_text))
    if not matches:
        raise ColvarParseError("No &COLLECTIVE block found in &CONSTRAINT")

    constraints = tuple(
        _parse_single_collective_block(m.group(1)) for m in matches
    )
    return Cp2kColvarInfoRaw(constraints=constraints)


def _parse_lagrange_filename(text: str) -> str | None:
    constraint_match = _CONSTRAINT_BLOCK_RE.search(text)
    if not constraint_match:
        return None
    lag_match = _LAGRANGE_BLOCK_RE.search(constraint_match.group(1))
    if not lag_match:
        return None
    fn = _extract_scalar(lag_match.group(1), "FILENAME")
    return fn.strip() if fn else None


def _parse_fixed_atoms_list(text: str) -> tuple[int, ...] | None:
    constraint_match = _CONSTRAINT_BLOCK_RE.search(text)
    if not constraint_match:
        return None
    fa_match = _FIXED_ATOMS_BLOCK_RE.search(constraint_match.group(1))
    if not fa_match:
        return None
    block = fa_match.group(1)

    # Collect LIST line(s), handling \ continuation
    list_text = ""
    collecting = False
    for line in block.split("\n"):
        stripped = line.strip()
        if not collecting:
            if stripped.upper().startswith("LIST"):
                content = stripped[4:].strip()  # remove 'LIST' prefix
                if content.endswith("\\"):
                    list_text += content[:-1] + " "
                    collecting = True
                else:
                    list_text += content
        else:
            if stripped.endswith("\\"):
                list_text += stripped[:-1] + " "
            else:
                list_text += stripped
                collecting = False

    # Parse tokens, expanding N..M ranges
    indices: list[int] = []
    for token in list_text.replace(",", " ").split():
        if ".." in token:
            parts = token.split("..")
            start, end = int(parts[0]), int(parts[1])
            indices.extend(range(start, end + 1))
        else:
            indices.append(int(token))

    return tuple(sorted(indices))


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

    Reuses :func:`~md_analysis.utils.formats.cp2k.cell.parse_abc_from_restart`
    for cell parameters.
    """
    path = Path(restart_path)
    text = path.read_text(encoding="utf-8")

    proj_match = re.search(
        r"^\s*PROJECT_NAME\s+(\S+)", text, re.MULTILINE | re.IGNORECASE,
    )
    if not proj_match:
        raise ColvarParseError(f"PROJECT_NAME not found in {restart_path}")

    md = _parse_md_block(text)
    colvars = _parse_all_collective_blocks(text)
    lagrange_filename = _parse_lagrange_filename(text)
    cell_abc_ang = parse_abc_from_restart(restart_path)
    fixed_atoms = _parse_fixed_atoms_list(text)

    return Cp2kConstraintMetadataRaw(
        project_name=proj_match.group(1),
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
