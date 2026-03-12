"""Slow-growth thermodynamic integration data structures.

Public API
----------
- ``Slowgrowth``       — base class (analysis-ready arrays)
- ``SlowgrowthFull``   — constructed from parsed restart + Lagrange log
- ``SlowgrowthSegment`` — slice of a Full with re-zeroed free energy
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ...utils.config import AU_TIME_TO_FS
from ...utils.RestartParser.ColvarParser import ColvarMDInfo


# ---------------------------------------------------------------------------
# Base class
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Slowgrowth:
    """Analysis-ready slow-growth data arrays.

    All array fields have shape ``(n_steps,)`` and are aligned by index.

    Attributes
    ----------
    steps : np.ndarray
        Absolute MD step numbers.
    times_fs : np.ndarray
        Simulation times in femtoseconds.
    target_au : np.ndarray
        Target collective-variable values in atomic units.
    lagrange_shake : np.ndarray
        Shake Lagrange multipliers (constraint force).
    free_energy_au : np.ndarray
        Cumulative free-energy change in atomic units (Hartree).
    timestep_fs : float
        MD timestep in femtoseconds.
    target_growth_au : float
        Target growth per step in atomic units (dξ/step).
    """

    steps: np.ndarray
    times_fs: np.ndarray
    target_au: np.ndarray
    lagrange_shake: np.ndarray
    free_energy_au: np.ndarray
    timestep_fs: float
    target_growth_au: float

    @property
    def n_steps(self) -> int:
        return len(self.steps)

    def reversed(self) -> Slowgrowth:
        """Return a new object with initial and final states swapped.

        - ``target_au``, ``lagrange_shake``: flipped
        - ``free_energy_au``: negated, flipped, then re-zeroed
        - ``target_growth_au``: negated
        - ``times_fs``: reset to start from 0
        - ``steps``: reset to ``[0, 1, ..., n-1]``
        """
        fe_reversed = -self.free_energy_au[::-1]
        fe_reversed = fe_reversed - fe_reversed[0]

        return Slowgrowth(
            steps=np.arange(self.n_steps),
            times_fs=np.arange(self.n_steps) * self.timestep_fs,
            target_au=self.target_au[::-1].copy(),
            lagrange_shake=self.lagrange_shake[::-1].copy(),
            free_energy_au=fe_reversed,
            timestep_fs=self.timestep_fs,
            target_growth_au=-self.target_growth_au,
        )


# ---------------------------------------------------------------------------
# Full trajectory
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SlowgrowthFull(Slowgrowth):
    """Complete slow-growth trajectory parsed from restart + Lagrange log.

    Retains a reference to the original :class:`ColvarMDInfo` for access
    to full metadata (cell parameters, fixed atoms, etc.).
    """

    md_info: ColvarMDInfo = None  # type: ignore[assignment]

    @classmethod
    def from_md_info(
        cls,
        md_info: ColvarMDInfo,
        *,
        colvar_id: int | None = None,
    ) -> SlowgrowthFull:
        """Build from a :class:`ColvarMDInfo` instance.

        Parameters
        ----------
        md_info : ColvarMDInfo
            Parsed restart + Lagrange multiplier data.
        colvar_id : int, optional
            Which collective variable to use.  Defaults to the primary CV.
        """
        constraint = (
            md_info.restart.colvars[colvar_id]
            if colvar_id is not None
            else md_info.restart.colvars.primary
        )
        # Convert growth rate from per-a.u.-time to per-step (dξ/step)
        dt_au = md_info.restart.timestep_fs / AU_TIME_TO_FS
        target_growth_per_step = constraint.target_growth_au * dt_au

        steps = md_info.steps
        times_fs = md_info.times_fs
        target_au = md_info.target_series_au(colvar_id)
        lagrange_shake = md_info.lagrange.collective_shake

        free_energy_au = _integrate_midpoint(lagrange_shake, target_growth_per_step)

        return cls(
            steps=steps,
            times_fs=times_fs,
            target_au=target_au,
            lagrange_shake=lagrange_shake,
            free_energy_au=free_energy_au,
            timestep_fs=md_info.restart.timestep_fs,
            target_growth_au=target_growth_per_step,
            md_info=md_info,
        )

    @classmethod
    def from_paths(
        cls,
        restart_path: str,
        log_path: str,
        *,
        colvar_id: int | None = None,
    ) -> SlowgrowthFull:
        """Parse files and build in one step."""
        md_info = ColvarMDInfo.from_paths(restart_path, log_path)
        return cls.from_md_info(md_info, colvar_id=colvar_id)

    def segment(
        self,
        start: int,
        end: int,
        *,
        ref_step: int | None = None,
    ) -> SlowgrowthSegment:
        """Extract a segment with re-zeroed free energy.

        Parameters
        ----------
        start, end : int
            Array indices (not step numbers) defining the half-open slice
            ``[start, end)``.
        ref_step : int, optional
            Array index whose ``free_energy_au`` becomes the new zero.
            Defaults to *start*.
        """
        if ref_step is None:
            ref_step = start

        ref_value = self.free_energy_au[ref_step]

        return SlowgrowthSegment(
            steps=self.steps[start:end].copy(),
            times_fs=self.times_fs[start:end].copy(),
            target_au=self.target_au[start:end].copy(),
            lagrange_shake=self.lagrange_shake[start:end].copy(),
            free_energy_au=self.free_energy_au[start:end] - ref_value,
            timestep_fs=self.timestep_fs,
            target_growth_au=self.target_growth_au,
            parent=self,
        )


# ---------------------------------------------------------------------------
# Segment (slice of a Full)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SlowgrowthSegment(Slowgrowth):
    """A slice of a :class:`SlowgrowthFull` with re-zeroed free energy.

    Created via :meth:`SlowgrowthFull.segment`.
    """

    parent: SlowgrowthFull = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Integration
# ---------------------------------------------------------------------------

def _integrate_midpoint(
    lagrange_shake: np.ndarray,
    target_growth_au: float,
) -> np.ndarray:
    r"""Cumulative free-energy change via midpoint rectangle rule.

    CP2K convention: the free-energy change is the *negative* integral
    of the constraint force (Lagrange multiplier) over the collective
    variable:

    .. math::
        \Delta A_k = -\sum_{i=0}^{k-1}
            \frac{\lambda_i + \lambda_{i+1}}{2} \, \Delta\xi

    where :math:`\Delta\xi` = *target_growth_au* (constant step size).

    Returns an array of length ``len(lagrange_shake)`` with
    ``result[0] = 0``.
    """
    n = len(lagrange_shake)
    free_energy = np.empty(n, dtype=np.float64)
    free_energy[0] = 0.0
    if n > 1:
        midpoints = 0.5 * (lagrange_shake[:-1] + lagrange_shake[1:])
        np.cumsum(-midpoints * target_growth_au, out=free_energy[1:])
    return free_energy
