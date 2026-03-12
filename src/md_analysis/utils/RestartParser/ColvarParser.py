"""Parse CP2K COLVAR restart files and LagrangeMultLog files."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np

from .CellParser import parse_abc_from_restart


from ...exceptions import MDAnalysisError


class ColvarParseError(MDAnalysisError):
    """Raised when parsing a COLVAR restart or LagrangeMultLog file fails."""


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ConstraintInfo:
    """COLLECTIVE constraint parameters."""

    colvar_id: int
    target_au: float
    target_growth_au: float
    intermolecular: bool


@dataclass(frozen=True)
class ColvarInfo:
    """Collection of collective variable constraints."""

    constraints: tuple[ConstraintInfo, ...]

    def __len__(self) -> int:
        return len(self.constraints)

    def __getitem__(self, colvar_id: int) -> ConstraintInfo:
        for c in self.constraints:
            if c.colvar_id == colvar_id:
                return c
        raise KeyError(f"No constraint with colvar_id={colvar_id}")

    def __iter__(self) -> Iterator[ConstraintInfo]:
        return iter(self.constraints)

    @property
    def primary(self) -> ConstraintInfo:
        """Return the first constraint (primary CV)."""
        return self.constraints[0]


@dataclass(frozen=True)
class ColvarRestart:
    """Metadata parsed from a COLVAR restart file."""

    project_name: str
    step_start: int
    time_start_fs: float
    timestep_fs: float
    total_steps: int
    colvars: ColvarInfo
    lagrange_filename: str
    cell_abc_ang: tuple[float, float, float]
    fixed_atom_indices: tuple[int, ...] | None


@dataclass(frozen=True)
class LagrangeMultLog:
    """Lagrange multiplier time series."""

    shake: np.ndarray
    rattle: np.ndarray
    n_steps: int
    n_constraints: int

    @property
    def collective_shake(self) -> np.ndarray:
        """Shake multiplier for the CV constraint, shape ``(n_steps,)``."""
        return self.shake if self.n_constraints == 1 else self.shake[:, 0]

    @property
    def collective_rattle(self) -> np.ndarray:
        """Rattle multiplier for the CV constraint, shape ``(n_steps,)``."""
        return self.rattle if self.n_constraints == 1 else self.rattle[:, 0]


@dataclass(frozen=True)
class ColvarMDInfo:
    """Complete slow-growth MD session: restart config + Lagrange multiplier data.

    Combines :class:`ColvarRestart` (input configuration) with
    :class:`LagrangeMultLog` (output multiplier data) and provides
    correctly aligned step/time/target arrays.
    """

    restart: ColvarRestart
    lagrange: LagrangeMultLog

    @property
    def n_steps(self) -> int:
        """Number of MD steps (from Lagrange multiplier log)."""
        return self.lagrange.n_steps

    @property
    def steps(self) -> np.ndarray:
        """Absolute step numbers, shape ``(n_steps,)``: ``[0, 1, ..., n_steps-1]``."""
        return np.arange(self.n_steps)

    @property
    def times_fs(self) -> np.ndarray:
        """Absolute times in fs, shape ``(n_steps,)``."""
        return self.steps * self.restart.timestep_fs

    def target_series_au(self, colvar_id: int | None = None) -> np.ndarray:
        """Target CV series in atomic units, shape ``(n_steps,)``.

        ``xi(k) = target_au + (k - step_start) * target_growth_au``

        where *k* are absolute step numbers ``[0, 1, ..., n_steps-1]``.
        """
        c = (
            self.restart.colvars[colvar_id]
            if colvar_id is not None
            else self.restart.colvars.primary
        )
        return c.target_au + (self.steps - self.restart.step_start) * c.target_growth_au

    @classmethod
    def from_paths(
        cls,
        restart_path: str | Path,
        log_path: str | Path,
    ) -> ColvarMDInfo:
        """Parse restart and LagrangeMultLog files into a single object."""
        return cls(
            restart=parse_colvar_restart(restart_path),
            lagrange=parse_lagrange_mult_log(log_path),
        )


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


def _parse_single_collective_block(block: str) -> ConstraintInfo:
    inter_val = _extract_scalar(block, "INTERMOLECULAR")
    intermolecular = False
    if inter_val is not None:
        intermolecular = inter_val.upper().strip(".") in ("T", "TRUE", "YES")

    return ConstraintInfo(
        colvar_id=int(_require_scalar(block, "COLVAR", "&COLLECTIVE")),
        target_au=float(_require_scalar(block, "TARGET", "&COLLECTIVE")),
        target_growth_au=float(
            _require_scalar(block, "TARGET_GROWTH", "&COLLECTIVE")
        ),
        intermolecular=intermolecular,
    )


def _parse_all_collective_blocks(text: str) -> ColvarInfo:
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
    return ColvarInfo(constraints=constraints)


def _parse_lagrange_filename(text: str) -> str:
    constraint_match = _CONSTRAINT_BLOCK_RE.search(text)
    if not constraint_match:
        raise ColvarParseError("No &CONSTRAINT block found")
    lag_match = _LAGRANGE_BLOCK_RE.search(constraint_match.group(1))
    if not lag_match:
        raise ColvarParseError("No &LAGRANGE_MULTIPLIERS block found")
    fn = _require_scalar(lag_match.group(1), "FILENAME", "&LAGRANGE_MULTIPLIERS")
    return fn.strip()


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
# Public API
# ---------------------------------------------------------------------------


def parse_colvar_restart(restart_path: str | Path) -> ColvarRestart:
    """Parse COLVAR metadata from a CP2K restart file.

    Reuses :func:`~md_analysis.utils.CellParser.parse_abc_from_restart` for
    cell parameters.
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

    return ColvarRestart(
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


def parse_lagrange_mult_log(log_path: str | Path) -> LagrangeMultLog:
    """Parse a LagrangeMultLog file. Auto-detects single/multi constraint."""
    path = Path(log_path)
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ColvarParseError(f"LagrangeMultLog file is empty: {log_path}")

    fmt = _detect_log_format(lines)

    if fmt == "single":
        shake, rattle, n_steps = _parse_single_constraint_log(lines)
        return LagrangeMultLog(
            shake=shake, rattle=rattle, n_steps=n_steps, n_constraints=1,
        )

    shake, rattle, n_steps, n_constraints = _parse_multi_constraint_log(lines)
    return LagrangeMultLog(
        shake=shake, rattle=rattle, n_steps=n_steps, n_constraints=n_constraints,
    )


def compute_target_series(
    restart: ColvarRestart,
    n_steps: int,
    *,
    colvar_id: int | None = None,
) -> np.ndarray:
    """Reconstruct the target CV series in atomic units.

    ``xi(k) = target_au + (k - step_start) * target_growth_au``
    where *k* = 0, 1, ..., *n_steps* - 1 (absolute step numbers).

    Parameters
    ----------
    restart : ColvarRestart
        Parsed restart metadata.  ``target_au`` is the target value
        **at** ``step_start`` (the restart snapshot), not the initial value.
    n_steps : int
        Number of steps to generate.
    colvar_id : int, optional
        If given, use the constraint with this ``colvar_id``.
        Defaults to the primary (first) constraint.
    """
    if colvar_id is not None:
        constraint = restart.colvars[colvar_id]
    else:
        constraint = restart.colvars.primary
    k = np.arange(n_steps)
    return (
        constraint.target_au
        + (k - restart.step_start) * constraint.target_growth_au
    )
