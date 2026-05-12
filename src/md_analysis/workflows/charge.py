"""Charge analysis workflow facade.

Exposes three programmatic entry points covering the Bader-based
surface-charge analysis surface:

- :func:`run_surface_charge`     — sigma(t) time series via counterion
  or layer method (writes into ``output_dir/<method>/``).
- :func:`run_tracked_charge`     — per-atom Bader-net-charge tracking
  by XYZ index (writes into ``output_dir/tracked/``).
- :func:`run_counterion_charge`  — per-frame counterion detection +
  charge time series (writes into ``output_dir/counterion_tracking/``).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Iterable

from ..utils.constants import CHARGE_METHOD_COUNTERION, DEFAULT_LAYER_TOL_A
from .models import WorkflowResult

logger = logging.getLogger(__name__)


def run_surface_charge(
    *,
    output_dir: Path | str,
    root_dir: str | Path = ".",
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
    method: str = CHARGE_METHOD_COUNTERION,
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    n_surface_layers: int = 1,
    dir_pattern: str = "bader_t*_i*",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> WorkflowResult:
    """Run surface charge density analysis (sigma vs t, CSV + PNG).

    Outputs are written into ``output_dir/<method>/``. The method
    sub-directory is mandatory so ``counterion`` and ``layer`` runs
    do not collide. Sigma statistics and any phi-axis metadata are
    surfaced through ``WorkflowResult.metadata``.
    """
    from ..electrochemical.charge import surface_charge_analysis
    from ..electrochemical.charge.config import DEFAULT_SURFACE_CHARGE_PNG_NAME

    root_p = Path(root_dir)
    base_dir = Path(output_dir)
    charge_dir = base_dir / method
    charge_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "Starting surface charge analysis: method=%s, output_dir=%s",
        method,
        charge_dir,
    )

    result = surface_charge_analysis(
        root_p,
        metal_symbols=metal_symbols,
        normal=normal,
        method=method,
        layer_tol_A=layer_tol_A,
        n_surface_layers=n_surface_layers,
        dir_pattern=dir_pattern,
        output_dir=charge_dir,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )

    csv_path = Path(result.csv_path)
    artifacts: dict[str, Path] = {
        "charge_csv": csv_path,
        "charge_png": csv_path.parent / DEFAULT_SURFACE_CHARGE_PNG_NAME,
    }
    metadata: dict[str, Any] = {
        "method": method,
        "normal": normal,
        "n_frames": result.n_frames,
        "sigma_aligned_mean": result.sigma_aligned_mean,
        "sigma_aligned_std": result.sigma_aligned_std,
        "sigma_opposed_mean": result.sigma_opposed_mean,
        "sigma_opposed_std": result.sigma_opposed_std,
        "phi_cumavg_last": result.phi_cumavg_last,
        "phi_reference": result.phi_reference,
    }
    return WorkflowResult(
        name="surface_charge",
        output_dir=charge_dir,
        artifacts=artifacts,
        metadata=metadata,
    )


def run_tracked_charge(
    *,
    output_dir: Path | str,
    root_dir: str | Path = ".",
    atom_indices_xyz: Iterable[int],
    dir_pattern: str = "bader_t*_i*",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> WorkflowResult:
    """Track Bader net charges for specified XYZ atom indices (CSV + PNG).

    Outputs are written into ``output_dir/tracked/``. Callers that want
    the canonical ``<root>/electrochemical/charge/tracked/`` layout
    should pass ``output_dir=<root>/electrochemical/charge``.
    """
    from ..electrochemical.charge import tracked_atom_charge_analysis
    from ..electrochemical.charge.Bader.AtomCharges import DEFAULT_TRACKED_CHARGE_PNG

    root_p = Path(root_dir)
    base_dir = Path(output_dir)
    tracked_dir = base_dir / "tracked"
    tracked_dir.mkdir(parents=True, exist_ok=True)

    # Materialise atom indices once so we can both forward them and
    # report the count in metadata without iterating a generator twice.
    indices_list = list(atom_indices_xyz)

    logger.info("Starting tracked charge analysis: output_dir=%s", tracked_dir)

    csv_path_raw = tracked_atom_charge_analysis(
        root_p,
        atom_indices_xyz=indices_list,
        dir_pattern=dir_pattern,
        output_dir=tracked_dir,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )

    csv_path = Path(csv_path_raw)
    artifacts: dict[str, Path] = {
        "tracked_charge_csv": csv_path,
        "tracked_charge_png": csv_path.parent / DEFAULT_TRACKED_CHARGE_PNG,
    }
    metadata: dict[str, Any] = {
        "n_atoms_tracked": len(indices_list),
        "atom_indices_xyz": indices_list,
    }
    return WorkflowResult(
        name="tracked_charge",
        output_dir=tracked_dir,
        artifacts=artifacts,
        metadata=metadata,
    )


def run_counterion_charge(
    *,
    output_dir: Path | str,
    root_dir: str | Path = ".",
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    dir_pattern: str = "bader_t*_i*",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> WorkflowResult:
    """Detect counterions per-frame and track their Bader charges (CSV + PNG).

    Outputs are written into ``output_dir/counterion_tracking/``.
    Callers that want the canonical
    ``<root>/electrochemical/charge/counterion_tracking/`` layout
    should pass ``output_dir=<root>/electrochemical/charge``.

    The PNG artifact is only listed when the underlying analysis
    actually produced one. If no counterion was detected in any
    frame, ``plot_counterion_charges`` short-circuits without
    writing a file and ``artifacts`` correspondingly omits
    ``counterion_charge_png`` so :func:`require_artifacts_exist`
    does not flag an empty-counterion run as broken.
    """
    from ..electrochemical.charge.Bader.AtomCharges import (
        counterion_charge_analysis_with_report,
    )

    root_p = Path(root_dir)
    base_dir = Path(output_dir)
    ci_dir = base_dir / "counterion_tracking"
    ci_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Starting counterion charge analysis: output_dir=%s", ci_dir)

    report = counterion_charge_analysis_with_report(
        root_p,
        metal_symbols=metal_symbols,
        normal=normal,
        layer_tol_A=layer_tol_A,
        dir_pattern=dir_pattern,
        output_dir=ci_dir,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )

    artifacts: dict[str, Path] = {
        "counterion_charge_csv": Path(report.csv_path),
        "counterion_summary_csv": Path(report.summary_path),
    }
    if report.png_path is not None:
        artifacts["counterion_charge_png"] = Path(report.png_path)
    metadata: dict[str, Any] = {
        "normal": normal,
        "layer_tol_A": layer_tol_A,
        "n_frames": report.n_frames,
        "n_unique_counterions": report.n_unique_counterions,
    }
    return WorkflowResult(
        name="counterion_charge",
        output_dir=ci_dir,
        artifacts=artifacts,
        metadata=metadata,
    )
