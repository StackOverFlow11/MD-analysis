"""Water analysis workflow facade.

Exposes :func:`run_water_three_panel`, the stable programmatic entry
point that combines density + orientation + adsorbed-layer analyses
into the canonical three-panel plot + supporting CSV/TXT artifacts.

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
