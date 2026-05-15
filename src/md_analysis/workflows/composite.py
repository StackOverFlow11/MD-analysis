"""Composite workflows that wire together leaf workflow facades.

A *composite* combines two or more leaf workflows under a single
canonical output layout. Composites only orchestrate parameter
forwarding and artifact aggregation; they never re-implement
science-level logic.

Currently exposes one composite:

- :func:`run_interface_analysis` — runs the water three-panel
  analysis plus the full Hartree-potential pipeline under the
  conventional ``<output_dir>/water/`` and
  ``<output_dir>/electrochemical/potential/`` sub-tree. This
  replaces the legacy ``md_analysis.main.run_all`` (which
  misleadingly suggested "all" workflows but in fact only ran
  water + potential).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from ..electrochemical.potential.config import DEFAULT_THICKNESS_ANG
from ..utils.constants import DEFAULT_LAYER_TOL_A
from .models import WorkflowResult
from .potential import run_potential_full
from .water import run_water_three_panel

logger = logging.getLogger(__name__)


def run_interface_analysis(
    xyz_path: Path | str,
    md_inp_path: Path | str | None = None,
    *,
    cell_abc: tuple[float, float, float] | None = None,
    output_dir: Path | str,
    # --- Potential-side inputs ---
    cube_pattern: str = "md-POTENTIAL-v_hartree-1_*.cube",
    md_out_path: Path | str | None = None,
    thickness_ang: float = DEFAULT_THICKNESS_ANG,
    center_mode: str = "interface",
    metal_elements: set[str] | None = None,
    layer_tol_ang: float = DEFAULT_LAYER_TOL_A,
    fermi_unit: str = "au",
    compute_u: bool = True,
    compute_phi_z: bool = True,
    max_curves: int = 0,
    thickness_end: float = 15.0,
    # --- Shared time-slicing controls ---
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> WorkflowResult:
    """Run water + potential under the canonical interface output layout.

    Pipeline:

    1. :func:`run_water_three_panel` writes into ``<output_dir>/water/``
    2. :func:`run_potential_full` writes into
       ``<output_dir>/electrochemical/potential/``

    Scope: water + potential **only**. This composite intentionally
    does NOT run charge / calibration / constrained TI / Bader
    work-directory generation; those workflows have separate
    facades in ``md_analysis.workflows.*`` and should be composed
    explicitly when needed.

    Artifact merge: water and potential leaf workflows have disjoint
    artifact keys today, so the merged ``WorkflowResult.artifacts``
    preserves each leaf's key verbatim. If a future leaf change
    introduces a collision the composite raises :class:`RuntimeError`
    rather than silently dropping one leaf's file — this mirrors the
    legacy ``run_all_with_report`` invariant.

    Returns
    -------
    WorkflowResult
        ``name="interface_analysis"``;
        ``output_dir`` is the composite root (the caller's
        ``output_dir``); ``artifacts`` merges the two leaf maps;
        ``metadata`` summarises each sub-workflow plus shared
        slicing parameters.

    Raises
    ------
    RuntimeError
        If water and potential leaf workflows produce conflicting
        artifact keys.
    """
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    water_dir = root / "water"
    potential_dir = root / "electrochemical" / "potential"

    xyz_p = Path(xyz_path)
    md_inp_p = Path(md_inp_path) if md_inp_path is not None else None
    md_out_p = Path(md_out_path) if md_out_path is not None else None

    logger.info("Starting interface analysis: output_dir=%s", root)

    water_result = run_water_three_panel(
        xyz_path=xyz_p,
        md_inp_path=md_inp_p,
        cell_abc=cell_abc,
        output_dir=water_dir,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )

    potential_result = run_potential_full(
        output_dir=potential_dir,
        cube_pattern=cube_pattern,
        md_out_path=md_out_p,
        xyz_path=xyz_p,
        thickness_ang=thickness_ang,
        center_mode=center_mode,
        metal_elements=metal_elements,
        layer_tol_ang=layer_tol_ang,
        fermi_unit=fermi_unit,
        compute_u=compute_u,
        compute_phi_z=compute_phi_z,
        max_curves=max_curves,
        thickness_end=thickness_end,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )

    collisions = sorted(
        set(water_result.artifacts) & set(potential_result.artifacts)
    )
    if collisions:
        raise RuntimeError(
            "run_interface_analysis refuses to merge: water + potential "
            f"leaves share artifact keys {collisions!r}. Fix upstream "
            "filename convention instead of silently overwriting one leaf's file."
        )

    merged_artifacts: dict[str, Path] = {}
    merged_artifacts.update(water_result.artifacts)
    merged_artifacts.update(potential_result.artifacts)

    metadata: dict[str, Any] = {
        "sub_workflows": ["water_three_panel", "potential_full"],
        "water_output_dir": str(water_result.output_dir),
        "potential_output_dir": str(potential_result.output_dir),
        "water_n_artifacts": len(water_result.artifacts),
        "potential_n_artifacts": len(potential_result.artifacts),
        "n_artifacts": len(merged_artifacts),
        # Forward selected sub-workflow flags verbatim so callers can
        # tell which potential sub-analyses actually ran without
        # rummaging through the merged artifact map.
        "potential_input_mode": potential_result.metadata.get("input_mode"),
        "potential_sub_analyses": potential_result.metadata.get(
            "sub_analyses", []
        ),
        "potential_has_fermi": potential_result.metadata.get("has_fermi"),
        "frame_start": frame_start,
        "frame_end": frame_end,
        "frame_step": frame_step,
    }
    return WorkflowResult(
        name="interface_analysis",
        output_dir=root,
        artifacts=merged_artifacts,
        metadata=metadata,
    )
