"""Parse CP2K slow-growth restart files and LagrangeMultLog files."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .CellParser import parse_abc_from_restart


class SlowGrowthParseError(RuntimeError):
    """Raised when parsing a slow-growth file fails."""


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ColvarDef:
    """Collective variable definition."""

    cv_type: str  # "DISTANCE" | "ANGLE" | "COMBINE_COLVAR"
    atoms: tuple[int, ...] | None
    function: str | None
    variables: tuple[str, ...] | None
    components: tuple[ColvarDef, ...] | None


@dataclass(frozen=True)
class ConstraintInfo:
    """COLLECTIVE constraint parameters."""

    colvar_id: int
    target_au: float
    target_growth_au: float
    intermolecular: bool


@dataclass(frozen=True)
class SlowGrowthRestart:
    """Metadata parsed from a slow-growth restart file."""

    project_name: str
    step_start: int
    time_start_fs: float
    timestep_fs: float
    total_steps: int
    constraint: ConstraintInfo
    lagrange_filename: str
    colvar: ColvarDef
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
_COMBINE_COLVAR_RE = re.compile(
    r"&COMBINE_COLVAR\s*\n(.*?)&END\s+COMBINE_COLVAR",
    re.DOTALL | re.IGNORECASE,
)
_ANGLE_RE = re.compile(
    r"&ANGLE\s*\n\s*ATOMS\s+([\d\s]+?)\s*\n\s*&END\s+ANGLE",
    re.IGNORECASE,
)
_DISTANCE_RE = re.compile(
    r"&DISTANCE\s*\n\s*ATOMS\s+([\d\s]+?)\s*\n\s*&END\s+DISTANCE",
    re.IGNORECASE,
)


def _extract_scalar(block: str, key: str) -> str | None:
    m = re.search(rf"^\s*{key}\s+(\S+)", block, re.MULTILINE | re.IGNORECASE)
    return m.group(1) if m else None


def _require_scalar(block: str, key: str, context: str) -> str:
    val = _extract_scalar(block, key)
    if val is None:
        raise SlowGrowthParseError(f"{key} not found in {context}")
    return val


def _parse_md_block(text: str) -> dict:
    md_match = re.search(r"&MD\b(.*?)&END\s+MD", text, re.DOTALL | re.IGNORECASE)
    if not md_match:
        raise SlowGrowthParseError("No &MD block found")
    block = md_match.group(1)
    return {
        "step_start": int(_require_scalar(block, "STEP_START_VAL", "&MD")),
        "time_start_fs": float(_require_scalar(block, "TIME_START_VAL", "&MD")),
        "timestep_fs": float(_require_scalar(block, "TIMESTEP", "&MD")),
        "total_steps": int(_require_scalar(block, "STEPS", "&MD")),
    }


def _parse_collective_block(text: str) -> ConstraintInfo:
    constraint_match = _CONSTRAINT_BLOCK_RE.search(text)
    if not constraint_match:
        raise SlowGrowthParseError("No &CONSTRAINT block found")
    constraint_text = constraint_match.group(1)

    coll_match = _COLLECTIVE_BLOCK_RE.search(constraint_text)
    if not coll_match:
        raise SlowGrowthParseError("No &COLLECTIVE block found in &CONSTRAINT")
    block = coll_match.group(1)

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


def _parse_lagrange_filename(text: str) -> str:
    constraint_match = _CONSTRAINT_BLOCK_RE.search(text)
    if not constraint_match:
        raise SlowGrowthParseError("No &CONSTRAINT block found")
    lag_match = _LAGRANGE_BLOCK_RE.search(constraint_match.group(1))
    if not lag_match:
        raise SlowGrowthParseError("No &LAGRANGE_MULTIPLIERS block found")
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


def _parse_simple_colvar(text: str) -> ColvarDef:
    angle_m = _ANGLE_RE.search(text)
    if angle_m:
        atoms = tuple(int(x) for x in angle_m.group(1).split())
        return ColvarDef("ANGLE", atoms=atoms, function=None, variables=None, components=None)

    dist_m = _DISTANCE_RE.search(text)
    if dist_m:
        atoms = tuple(int(x) for x in dist_m.group(1).split())
        return ColvarDef("DISTANCE", atoms=atoms, function=None, variables=None, components=None)

    raise SlowGrowthParseError(
        "No supported COLVAR type found (expected ANGLE or DISTANCE)"
    )


def _parse_colvar_block(text: str) -> ColvarDef:
    combine_m = _COMBINE_COLVAR_RE.search(text)
    if combine_m:
        block = combine_m.group(1)
        func_m = re.search(
            r"^\s*FUNCTION\s+(.+?)\s*$", block, re.MULTILINE | re.IGNORECASE,
        )
        vars_m = re.search(
            r"^\s*VARIABLES\s+(.+?)\s*$", block, re.MULTILINE | re.IGNORECASE,
        )
        function = func_m.group(1).strip() if func_m else None
        variables = tuple(vars_m.group(1).split()) if vars_m else None

        # Find nested &COLVAR ... &END COLVAR blocks inside COMBINE_COLVAR
        components: list[ColvarDef] = []
        for nested in re.finditer(
            r"&COLVAR\s*\n(.*?)&END\s+COLVAR", block, re.DOTALL | re.IGNORECASE,
        ):
            components.append(_parse_simple_colvar(nested.group(1)))

        return ColvarDef(
            cv_type="COMBINE_COLVAR",
            atoms=None,
            function=function,
            variables=variables,
            components=tuple(components),
        )

    return _parse_simple_colvar(text)


# ---------------------------------------------------------------------------
# Private helpers — LagrangeMultLog parsing
# ---------------------------------------------------------------------------

_LABEL_RE = re.compile(r"^(Shake|Rattle)\s+Lagrangian\s+Multipliers:", re.IGNORECASE)


def _detect_log_format(lines: list[str]) -> str:
    if len(lines) < 2:
        raise SlowGrowthParseError("LagrangeMultLog file is too short")
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
            shake_vals.append(float(stripped.split(":")[1].strip()))
        elif stripped.startswith("Rattle"):
            rattle_vals.append(float(stripped.split(":")[1].strip()))
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
            current_values = [float(x) for x in after_colon.split()]
        else:
            current_values.extend(float(x) for x in stripped.split())

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


def parse_slowgrowth_restart(restart_path: str | Path) -> SlowGrowthRestart:
    """Parse slow-growth metadata from a CP2K restart file.

    Reuses :func:`~md_analysis.utils.CellParser.parse_abc_from_restart` for
    cell parameters.
    """
    path = Path(restart_path)
    text = path.read_text(encoding="utf-8")

    proj_match = re.search(
        r"^\s*PROJECT_NAME\s+(\S+)", text, re.MULTILINE | re.IGNORECASE,
    )
    if not proj_match:
        raise SlowGrowthParseError(f"PROJECT_NAME not found in {restart_path}")

    md = _parse_md_block(text)
    constraint = _parse_collective_block(text)
    lagrange_filename = _parse_lagrange_filename(text)
    colvar = _parse_colvar_block(text)
    cell_abc_ang = parse_abc_from_restart(restart_path)
    fixed_atoms = _parse_fixed_atoms_list(text)

    return SlowGrowthRestart(
        project_name=proj_match.group(1),
        step_start=md["step_start"],
        time_start_fs=md["time_start_fs"],
        timestep_fs=md["timestep_fs"],
        total_steps=md["total_steps"],
        constraint=constraint,
        lagrange_filename=lagrange_filename,
        colvar=colvar,
        cell_abc_ang=cell_abc_ang,
        fixed_atom_indices=fixed_atoms,
    )


def parse_lagrange_mult_log(log_path: str | Path) -> LagrangeMultLog:
    """Parse a LagrangeMultLog file. Auto-detects single/multi constraint."""
    path = Path(log_path)
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise SlowGrowthParseError(f"LagrangeMultLog file is empty: {log_path}")

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


def compute_target_series(restart: SlowGrowthRestart, n_steps: int) -> np.ndarray:
    """Reconstruct the target CV series in atomic units.

    ``xi(k) = target_au + (k - step_start) * target_growth_au``
    where *k* = 1, 2, ..., *n_steps*.
    """
    k = np.arange(1, n_steps + 1)
    return (
        restart.constraint.target_au
        + (k - restart.step_start) * restart.constraint.target_growth_au
    )
