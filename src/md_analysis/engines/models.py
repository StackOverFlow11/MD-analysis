"""Engine-neutral data models.

Frozen dataclasses returned by engine adapter modules
(``engines.cp2k``, future ``engines.vasp``). Upper analysis layers
(``electrochemical``, ``enhanced_sampling``, ``water``) depend on these
neutral types rather than engine-specific parser output.

Constraint-MD types
-------------------
- ``ConstraintInfo`` / ``ColvarInfo``       — neutral nested types
  describing a single CV constraint and a collection thereof.
- ``ConstraintMetadata`` / ``LambdaSeries`` — canonical top-level
  payloads of a constraint-MD point (restart-side inputs and the
  Lagrange-multiplier time series).
- ``ConstraintRun``                         — composite view of a
  constraint-MD run (``metadata`` + ``lambda_series`` + derived
  target/time series).

Potential / Fermi types
-----------------------
- ``PotentialFrame``                — one frame of potential analysis
  data (cube header + values + Fermi-raw payload).
- ``CenterPotentialScalarFrame``    — scalar-level reduction of a
  ``PotentialFrame`` plus slab geometry.
- ``FermiRecord``                   — typed record mirroring the legacy
  ``parse_md_out_fermi`` dict shape.

Cell type
---------
- ``CellSpec`` — engine-neutral cell descriptor (3x3 matrix + pbc
  + derived ``abc_ang`` / ``is_orthorhombic`` properties).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np

from ..utils.constants import AU_TIME_TO_FS
from ..utils.formats.common.cube import CubeHeader

try:
    from ase import Atoms
except ImportError:  # pragma: no cover
    Atoms = object  # type: ignore[misc]


# ---------------------------------------------------------------------------
# Neutral nested types (Phase 4 Commit 1 — D10 nested-type migration)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ConstraintInfo:
    """COLLECTIVE constraint parameters (engine-neutral).

    ``target_growth_au`` is the rate of change per atomic unit of time
    (CP2K convention).  Multiply by the timestep in a.u. to obtain the
    per-step increment.
    """

    colvar_id: int
    target_au: float
    target_growth_au: float
    intermolecular: bool


@dataclass(frozen=True)
class ColvarInfo:
    """Collection of collective variable constraints (engine-neutral)."""

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


# ---------------------------------------------------------------------------
# Canonical top-level dataclasses (Phase 4 Commit 1 — D10 physical migration)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ConstraintMetadata:
    """Metadata parsed from a constraint-MD restart file (engine-neutral).

    Engine-neutral payload returned by ``ConstraintMDParser.parse_metadata``.
    """

    project_name: str
    step_start: int
    time_start_fs: float
    timestep_fs: float
    total_steps: int
    colvars: ColvarInfo
    lagrange_filename: str | None
    cell_abc_ang: tuple[float, float, float]
    fixed_atom_indices: tuple[int, ...] | None


@dataclass(frozen=True)
class LambdaSeries:
    """Lagrange multiplier (constraint force) time series (engine-neutral).

    Engine-neutral payload returned by ``ConstraintMDParser.parse_lambda_series``.
    The ``LagrangeMultLog`` file suffix (``*.LagrangeMultLog``) is a CP2K
    output filename and is NOT a type name.
    """

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
# Composite view (D8 ConstraintRun)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ConstraintRun:
    """Engine-neutral composite view of a constraint-MD run.

    Combines metadata (restart-time inputs) with the resulting
    Lagrange-multiplier time series.  Derived series are computed via
    ``@property`` (not stored as fields) to keep this object an
    "inputs-only" snapshot.

    Construction: the standard entry point is
    :func:`md_analysis.engines.cp2k.read_constraint_run_from_files`
    (file-level facade over a ``(restart_path, log_path)`` pair).  A
    ``from_paths`` classmethod is intentionally NOT provided on this
    class because it would force ``engines.models`` to import
    ``engines.cp2k`` and create a circular import.

    Derived ``target_series_au`` follows the CP2K SHAKE/RATTLE
    convention (formula:
    ``xi(k) = target_au + (k - step_start) * target_growth_au * dt_au``,
    where ``dt_au = timestep_fs / AU_TIME_TO_FS``).  Other engines must
    only construct this object when their semantics map losslessly onto
    the same formula; otherwise the adapter MUST raise
    :class:`NotImplementedError` rather than silently producing a
    half-valid object.
    """

    metadata: ConstraintMetadata
    lambda_series: LambdaSeries

    # ------------------------------------------------------------------
    # Derived properties (engine-neutral)
    # ------------------------------------------------------------------

    @property
    def n_steps(self) -> int:
        """Number of MD steps (from Lagrange multiplier log)."""
        return self.lambda_series.n_steps

    @property
    def steps(self) -> np.ndarray:
        """Absolute step numbers, shape ``(n_steps,)``: ``[0, 1, ..., n_steps-1]``."""
        return np.arange(self.n_steps)

    @property
    def times_fs(self) -> np.ndarray:
        """Absolute times in fs, shape ``(n_steps,)``."""
        return self.steps * self.metadata.timestep_fs

    def target_series_au(self, colvar_id: int | None = None) -> np.ndarray:
        """Target CV series in atomic units, shape ``(n_steps,)``.

        ``xi(k) = target_au + (k - step_start) * target_growth_au * dt_au``

        where *k* are absolute step numbers ``[0, 1, ..., n_steps-1]``
        and *dt_au* is the MD timestep in atomic time units.
        """
        c = (
            self.metadata.colvars[colvar_id]
            if colvar_id is not None
            else self.metadata.colvars.primary
        )
        dt_au = self.metadata.timestep_fs / AU_TIME_TO_FS
        return c.target_au + (self.steps - self.metadata.step_start) * c.target_growth_au * dt_au


# ---------------------------------------------------------------------------
# Cell descriptor (D7 Layer 1)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CellSpec:
    """Engine-neutral cell descriptor.

    Holds a 3x3 lattice matrix in row-vector convention (one lattice
    vector per row, matching ``ase.Atoms.cell.array`` and CP2K's
    ``&CELL A/B/C`` row order).

    Phase 4 status: the only engine adapter (CP2K's
    :func:`md_analysis.engines.cp2k.read_cell`) currently produces
    orthorhombic matrices only -- both ``parse_abc_from_md_inp`` and
    ``parse_abc_from_restart`` reject non-orthogonal cells.  The
    dataclass shape is intentionally general so future engine adapters
    (or a future non-orthogonal restart parser) can populate
    ``cell_matrix_ang`` without a model rev.
    """

    cell_matrix_ang: np.ndarray  # shape (3, 3), row-vector convention
    pbc: tuple[bool, bool, bool] = (True, True, True)

    @property
    def abc_ang(self) -> tuple[float, float, float]:
        """Norms of the three lattice vectors, in Angstrom.

        WARNING: For non-orthorhombic cells this is NOT the box-edge
        length.  Callers that rely on ``(a, b, c)`` as orthogonal-box
        semantics MUST first check :attr:`is_orthorhombic`; otherwise
        they must consume the full ``cell_matrix_ang`` (e.g. via
        ``ase.Atoms.cell``).
        """
        m = self.cell_matrix_ang
        return (
            float(np.linalg.norm(m[0])),
            float(np.linalg.norm(m[1])),
            float(np.linalg.norm(m[2])),
        )

    @property
    def is_orthorhombic(self) -> bool:
        """True iff off-diagonal elements of ``cell_matrix_ang`` are
        zero within a fixed tolerance of 1e-6 Angstrom (matching
        :func:`md_analysis.utils.formats.cp2k.cell.parse_abc_from_restart`'s
        orthogonality check)."""
        m = self.cell_matrix_ang
        tol = 1e-6
        off = (
            abs(m[0, 1]), abs(m[0, 2]),
            abs(m[1, 0]), abs(m[1, 2]),
            abs(m[2, 0]), abs(m[2, 1]),
        )
        return all(v <= tol for v in off)


# ---------------------------------------------------------------------------
# Potential frame + Fermi record (pre-existing, unchanged)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PotentialFrame:
    """One frame of potential analysis data (immutable, engine-neutral)."""

    step: int
    time_fs: float | None
    cube_path: Path
    header: CubeHeader
    values: np.ndarray
    fermi_raw: float | None  # Hartree; None if unavailable
    atoms: Atoms | None  # for interface detection (with cell set)


@dataclass(frozen=True)
class CenterPotentialScalarFrame:
    """Engine-neutral scalar-level frame for slab-averaged potential analysis.

    Distinct from :class:`PotentialFrame`: this carries scalar aggregates
    (slab-centered Hartree potential, Fermi level, slab geometry)
    instead of the heavy ``(cube_path, values, header)`` raw payload.
    Built by
    :func:`md_analysis.engines.cp2k.read_center_potential_scalar_frame`
    from an already-parsed :class:`PotentialFrame` plus slab geometry.

    Unit suffix convention (codex Round 3 §2.5): both ``phi_center_ev``
    and ``fermi_level_ev`` are in eV so business-layer cSHE formulas
    ``U = -E_Fermi + phi_center + ...`` consume them directly without
    re-scaling.

    Phase 4 R3: this frame intentionally does NOT carry cSHE reference
    quantities (DELTA_PSI_A_H3O_W_EV / MU_HPLUS_G0_EV / DELTA_E_ZP_EV)
    nor any U_vs_<reference> value -- those remain in
    ``electrochemical.potential``.
    """

    step: int
    time_fs: float | None
    center_source: str  # "interface" / "cell" / "manual"
    center_z_ang: float | None  # slab center z (A); None reserved for non-facade producers
    slab_thickness_ang: float
    phi_center_ev: float
    fermi_level_ev: float | None
    phi_z_std_ev: float | None = None
    n_slices: int | None = None


@dataclass(frozen=True)
class FermiRecord:
    """One ``(step, time_fs, fermi_raw)`` row from a CP2K stdout log.

    Engine-neutral typed counterpart to the legacy dict-of-records shape
    returned by :func:`md_analysis.utils.formats.cp2k.stdout.parse_md_out_fermi`.

    ``fermi_raw`` is in Hartree.
    """

    step: int
    time_fs: float | None
    fermi_raw: float

    @classmethod
    def from_legacy_dict(cls, d: dict) -> "FermiRecord":
        """Bridge a legacy ``parse_md_out_fermi`` dict into a typed record."""
        return cls(
            step=int(d["step"]),
            time_fs=d["time_fs"],
            fermi_raw=float(d["fermi_raw"]),
        )


__all__ = [
    # Neutral nested types
    "ConstraintInfo",
    "ColvarInfo",
    # Canonical top-level types
    "ConstraintMetadata",
    "LambdaSeries",
    # Composite
    "ConstraintRun",
    # Cell layer
    "CellSpec",
    # Potential layer
    "PotentialFrame",
    "CenterPotentialScalarFrame",
    "FermiRecord",
]
