"""Water analysis workflow facade.

Exposes the canonical programmatic entry points for the water
analysis suite:

- :func:`run_water_three_panel`     -- composite (density + orientation +
  adsorbed-layer) producing the three-panel plot + supporting CSVs.
- :func:`run_water_density`         -- single-step mass density.
- :func:`run_water_orientation`     -- single-step orientation-weighted density.
- :func:`run_ad_water_orientation`  -- single-step adsorbed-layer profile.
- :func:`run_ad_water_theta`        -- single-step adsorbed-layer theta distribution.

Business logic lives under :mod:`md_analysis.water`; this module
organises parameters, output paths, and the :class:`WorkflowResult`
contract.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from ..utils.constants import DEFAULT_LAYER_TOL_A
from .models import WorkflowResult

logger = logging.getLogger(__name__)


def run_water_three_panel(
    xyz_path: Path | str,
    md_inp_path: Path | str | None = None,
    *,
    cell_abc: tuple[float, float, float] | None = None,
    output_dir: Path | str,
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    **kwargs: Any,
) -> WorkflowResult:
    """Run the water three-panel analysis (density + orientation + adsorbed).

    Outputs (density CSV, orientation CSV, adsorbed profile/range/theta
    CSV+TXT, three-panel PNG) are written directly into ``output_dir``.
    Callers that want the canonical ``<root>/water/`` layout should pass
    ``output_dir=<root>/water`` explicitly.
    """
    from ..water import plot_water_three_panel_analysis
    from ..water.config import (
        DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME,
        DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME,
        DEFAULT_ADSORBED_WATER_THETA_DISTRIBUTION_CSV_NAME,
        DEFAULT_WATER_MASS_DENSITY_CSV_NAME,
        DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME,
    )

    xyz_path_p = Path(xyz_path)
    md_inp_path_p = Path(md_inp_path) if md_inp_path is not None else None
    water_dir = Path(output_dir)
    water_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "Starting water three-panel analysis: xyz=%s, output_dir=%s",
        xyz_path_p,
        water_dir,
    )

    png_path = plot_water_three_panel_analysis(
        xyz_path=xyz_path_p,
        md_inp_path=md_inp_path_p,
        cell_abc=cell_abc,
        output_dir=water_dir,
        layer_tol_A=layer_tol_A,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
        **kwargs,
    )

    artifacts: dict[str, Path] = {
        "density_csv": water_dir / DEFAULT_WATER_MASS_DENSITY_CSV_NAME,
        "orientation_csv": (
            water_dir / DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME
        ),
        "adsorbed_profile_csv": water_dir / DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME,
        "adsorbed_range_txt": water_dir / DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME,
        "adsorbed_theta_csv": (
            water_dir / DEFAULT_ADSORBED_WATER_THETA_DISTRIBUTION_CSV_NAME
        ),
        "plot_png": Path(png_path),
    }
    metadata: dict[str, Any] = {
        "cell_abc_provided": cell_abc is not None,
        "layer_tol_A": layer_tol_A,
        "frame_start": frame_start,
        "frame_end": frame_end,
        "frame_step": frame_step,
    }
    return WorkflowResult(
        name="water_three_panel",
        output_dir=water_dir,
        artifacts=artifacts,
        metadata=metadata,
    )


# ---------------------------------------------------------------------------
# Single-step facades (Phase 6.1)
# ---------------------------------------------------------------------------


def _common_metadata(
    *,
    cell_abc: tuple[float, float, float] | None,
    layer_tol_A: float,
    frame_start: int | None,
    frame_end: int | None,
    frame_step: int | None,
    **extra_metadata: Any,
) -> dict[str, Any]:
    """Build the common metadata block shared by all water workflows."""
    md: dict[str, Any] = {
        "cell_abc_provided": cell_abc is not None,
        "layer_tol_A": layer_tol_A,
        "frame_start": frame_start,
        "frame_end": frame_end,
        "frame_step": frame_step,
    }
    md.update(extra_metadata)
    return md


def run_water_density(
    xyz_path: Path | str,
    md_inp_path: Path | str | None = None,
    *,
    cell_abc: tuple[float, float, float] | None = None,
    output_dir: Path | str,
    dz_A: float | None = None,
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    **kwargs: Any,
) -> WorkflowResult:
    """Run water mass-density z-distribution analysis (single-step).

    Wraps :func:`md_analysis.water.water_mass_density_z_distribution_analysis`
    and packs the resulting CSV path into a :class:`WorkflowResult` with
    a stable artifact key ``density_csv``.

    Output (one CSV) is written directly into ``output_dir``.
    Callers that want the canonical ``<root>/water/`` layout should
    pass ``output_dir=<root>/water`` explicitly.
    """
    from ..water import water_mass_density_z_distribution_analysis

    xyz_path_p = Path(xyz_path)
    md_inp_path_p = Path(md_inp_path) if md_inp_path is not None else None
    water_dir = Path(output_dir)
    water_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "Starting water density analysis: xyz=%s, output_dir=%s",
        xyz_path_p, water_dir,
    )

    call_kwargs: dict[str, Any] = {
        "xyz_path": xyz_path_p,
        "md_inp_path": md_inp_path_p,
        "cell_abc": cell_abc,
        "output_dir": water_dir,
        "layer_tol_A": layer_tol_A,
        "frame_start": frame_start,
        "frame_end": frame_end,
        "frame_step": frame_step,
    }
    if dz_A is not None:
        call_kwargs["dz_A"] = dz_A
    call_kwargs.update(kwargs)

    csv_path = water_mass_density_z_distribution_analysis(**call_kwargs)

    return WorkflowResult(
        name="water_density",
        output_dir=water_dir,
        artifacts={"density_csv": Path(csv_path)},
        metadata=_common_metadata(
            cell_abc=cell_abc,
            layer_tol_A=layer_tol_A,
            frame_start=frame_start,
            frame_end=frame_end,
            frame_step=frame_step,
            dz_A=dz_A,
        ),
    )


def run_water_orientation(
    xyz_path: Path | str,
    md_inp_path: Path | str | None = None,
    *,
    cell_abc: tuple[float, float, float] | None = None,
    output_dir: Path | str,
    dz_A: float | None = None,
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    **kwargs: Any,
) -> WorkflowResult:
    """Run water orientation-weighted density z-distribution analysis (single-step).

    Wraps :func:`md_analysis.water.water_orientation_weighted_density_z_distribution_analysis`
    and exposes the resulting CSV under the stable artifact key
    ``orientation_csv``.
    """
    from ..water import water_orientation_weighted_density_z_distribution_analysis

    xyz_path_p = Path(xyz_path)
    md_inp_path_p = Path(md_inp_path) if md_inp_path is not None else None
    water_dir = Path(output_dir)
    water_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "Starting water orientation analysis: xyz=%s, output_dir=%s",
        xyz_path_p, water_dir,
    )

    call_kwargs: dict[str, Any] = {
        "xyz_path": xyz_path_p,
        "md_inp_path": md_inp_path_p,
        "cell_abc": cell_abc,
        "output_dir": water_dir,
        "layer_tol_A": layer_tol_A,
        "frame_start": frame_start,
        "frame_end": frame_end,
        "frame_step": frame_step,
    }
    if dz_A is not None:
        call_kwargs["dz_A"] = dz_A
    call_kwargs.update(kwargs)

    csv_path = water_orientation_weighted_density_z_distribution_analysis(**call_kwargs)

    return WorkflowResult(
        name="water_orientation",
        output_dir=water_dir,
        artifacts={"orientation_csv": Path(csv_path)},
        metadata=_common_metadata(
            cell_abc=cell_abc,
            layer_tol_A=layer_tol_A,
            frame_start=frame_start,
            frame_end=frame_end,
            frame_step=frame_step,
            dz_A=dz_A,
        ),
    )


def run_ad_water_orientation(
    xyz_path: Path | str,
    md_inp_path: Path | str | None = None,
    *,
    cell_abc: tuple[float, float, float] | None = None,
    output_dir: Path | str,
    dz_A: float | None = None,
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    **kwargs: Any,
) -> WorkflowResult:
    """Run adsorbed-layer water orientation analysis (single-step).

    Wraps :func:`md_analysis.water.ad_water_orientation_analysis`,
    which returns ``(profile_csv, range_txt)``, and exposes both
    artifacts under the stable keys ``adsorbed_profile_csv`` and
    ``adsorbed_range_txt``.
    """
    from ..water import ad_water_orientation_analysis

    xyz_path_p = Path(xyz_path)
    md_inp_path_p = Path(md_inp_path) if md_inp_path is not None else None
    water_dir = Path(output_dir)
    water_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "Starting adsorbed-water orientation analysis: xyz=%s, output_dir=%s",
        xyz_path_p, water_dir,
    )

    call_kwargs: dict[str, Any] = {
        "xyz_path": xyz_path_p,
        "md_inp_path": md_inp_path_p,
        "cell_abc": cell_abc,
        "output_dir": water_dir,
        "layer_tol_A": layer_tol_A,
        "frame_start": frame_start,
        "frame_end": frame_end,
        "frame_step": frame_step,
    }
    if dz_A is not None:
        call_kwargs["dz_A"] = dz_A
    call_kwargs.update(kwargs)

    profile_csv, range_txt = ad_water_orientation_analysis(**call_kwargs)

    return WorkflowResult(
        name="ad_water_orientation",
        output_dir=water_dir,
        artifacts={
            "adsorbed_profile_csv": Path(profile_csv),
            "adsorbed_range_txt": Path(range_txt),
        },
        metadata=_common_metadata(
            cell_abc=cell_abc,
            layer_tol_A=layer_tol_A,
            frame_start=frame_start,
            frame_end=frame_end,
            frame_step=frame_step,
            dz_A=dz_A,
        ),
    )


def run_ad_water_theta(
    xyz_path: Path | str,
    md_inp_path: Path | str | None = None,
    *,
    cell_abc: tuple[float, float, float] | None = None,
    output_dir: Path | str,
    dz_A: float | None = None,
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    **kwargs: Any,
) -> WorkflowResult:
    """Run adsorbed-layer water theta distribution analysis (single-step).

    Wraps :func:`md_analysis.water.compute_adsorbed_water_theta_distribution`,
    which returns ``(theta_centers, theta_pdf, csv)``.  Only the
    persisted CSV path is exposed on the result (under the stable key
    ``theta_csv``); the in-memory ``(centers, pdf)`` arrays are
    intentionally not surfaced -- callers who need the arrays should
    re-read the CSV.
    """
    from ..water import compute_adsorbed_water_theta_distribution

    xyz_path_p = Path(xyz_path)
    md_inp_path_p = Path(md_inp_path) if md_inp_path is not None else None
    water_dir = Path(output_dir)
    water_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "Starting adsorbed-water theta analysis: xyz=%s, output_dir=%s",
        xyz_path_p, water_dir,
    )

    call_kwargs: dict[str, Any] = {
        "xyz_path": xyz_path_p,
        "md_inp_path": md_inp_path_p,
        "cell_abc": cell_abc,
        "output_dir": water_dir,
        "layer_tol_A": layer_tol_A,
        "frame_start": frame_start,
        "frame_end": frame_end,
        "frame_step": frame_step,
        "verbose": verbose,
    }
    if dz_A is not None:
        call_kwargs["dz_A"] = dz_A
    call_kwargs.update(kwargs)

    _, _, csv_path = compute_adsorbed_water_theta_distribution(**call_kwargs)

    return WorkflowResult(
        name="ad_water_theta",
        output_dir=water_dir,
        artifacts={"theta_csv": Path(csv_path)},
        metadata=_common_metadata(
            cell_abc=cell_abc,
            layer_tol_A=layer_tol_A,
            frame_start=frame_start,
            frame_end=frame_end,
            frame_step=frame_step,
            dz_A=dz_A,
        ),
    )
