"""Potential analysis workflow facade.

Exposes :func:`run_potential_full`, the stable programmatic entry point
that orchestrates center-slab Hartree potential, Fermi-energy time
series, electrode potential (cSHE), planar-averaged phi(z) overlay,
and the thickness-sensitivity sweep.

Sub-analyses are routed into ``output_dir/<sub>/`` where ``<sub>`` is
one of ``center``, ``fermi``, ``electrode``, ``phi_z``,
``thickness_sensitivity``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from ..electrochemical.potential.config import DEFAULT_THICKNESS_ANG
from ..utils.constants import DEFAULT_LAYER_TOL_A
from .models import WorkflowResult

logger = logging.getLogger(__name__)


def run_potential_full(
    *,
    output_dir: Path | str,
    cube_pattern: str = "md-POTENTIAL-v_hartree-1_*.cube",
    md_out_path: Path | str | None = None,
    xyz_path: Path | str | None = None,
    thickness_ang: float = DEFAULT_THICKNESS_ANG,
    center_mode: str = "interface",
    metal_elements: set[str] | None = None,
    layer_tol_ang: float = DEFAULT_LAYER_TOL_A,
    fermi_unit: str = "au",
    compute_u: bool = True,
    compute_phi_z: bool = True,
    max_curves: int = 0,
    thickness_end: float = 15.0,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    # --- distributed mode params ---
    input_mode: str = "continuous",
    sp_root_dir: Path | str | None = None,
    sp_dir_pattern: str = "potential_t*_i*",
    sp_cube_filename: str = "sp_potential-v_hartree-1_0.cube",
    sp_out_filename: str = "sp.out",
) -> WorkflowResult:
    """Run the full Hartree-potential pipeline.

    The exact set of artifacts depends on the inputs:

    - When Fermi data is available (``input_mode='distributed'`` or
      ``md_out_path`` provided) and ``compute_u`` is True, the
      ``electrode`` sub-analysis is run; otherwise ``center`` runs
      standalone, and ``fermi`` runs separately when ``md_out_path``
      is provided.
    - ``phi_z`` runs when ``compute_phi_z`` is True.
    - ``thickness_sensitivity`` runs whenever Fermi data is available.

    All produced files are reported in ``WorkflowResult.artifacts``;
    sub-analysis flags are reported in ``metadata``.
    """
    from ..electrochemical.potential import (
        center_slab_potential_analysis,
        electrode_potential_analysis,
        fermi_energy_analysis,
        phi_z_planeavg_analysis,
        thickness_sensitivity_analysis,
    )

    pot_dir = Path(output_dir)
    pot_dir.mkdir(parents=True, exist_ok=True)
    md_out_p = Path(md_out_path) if md_out_path is not None else None
    xyz_p = Path(xyz_path) if xyz_path is not None else None
    sp_root_p = Path(sp_root_dir) if sp_root_dir is not None else None

    logger.info("Starting potential analysis: output_dir=%s", pot_dir)

    is_distributed = input_mode == "distributed"
    has_fermi = is_distributed or md_out_p is not None

    _dist = {
        "input_mode": input_mode,
        "sp_root_dir": sp_root_p,
        "sp_dir_pattern": sp_dir_pattern,
        "sp_cube_filename": sp_cube_filename,
        "sp_out_filename": sp_out_filename,
    }

    artifacts: dict[str, Path] = {}
    sub_analyses: list[str] = []

    if compute_u and has_fermi:
        electrode_dir = pot_dir / "electrode"
        electrode_dir.mkdir(parents=True, exist_ok=True)
        u_csv = electrode_potential_analysis(
            cube_pattern,
            md_out_p,
            output_dir=electrode_dir,
            thickness_ang=thickness_ang,
            center_mode=center_mode,
            xyz_path=xyz_p,
            metal_elements=metal_elements,
            layer_tol_ang=layer_tol_ang,
            fermi_unit=fermi_unit,
            frame_start=frame_start,
            frame_end=frame_end,
            frame_step=frame_step,
            verbose=verbose,
            **_dist,
        )
        artifacts["electrode_csv"] = Path(u_csv)
        sub_analyses.append("electrode")
    else:
        center_dir = pot_dir / "center"
        center_dir.mkdir(parents=True, exist_ok=True)
        center_csv = center_slab_potential_analysis(
            cube_pattern,
            output_dir=center_dir,
            thickness_ang=thickness_ang,
            center_mode=center_mode,
            xyz_path=xyz_p,
            metal_elements=metal_elements,
            layer_tol_ang=layer_tol_ang,
            frame_start=frame_start,
            frame_end=frame_end,
            frame_step=frame_step,
            verbose=verbose,
            **_dist,
        )
        artifacts["center_csv"] = Path(center_csv)
        sub_analyses.append("center")

        if md_out_p is not None:
            fermi_dir = pot_dir / "fermi"
            fermi_dir.mkdir(parents=True, exist_ok=True)
            fermi_csv = fermi_energy_analysis(
                md_out_p,
                output_dir=fermi_dir,
                fermi_unit=fermi_unit,
                frame_start=frame_start,
                frame_end=frame_end,
                frame_step=frame_step,
                **_dist,
            )
            artifacts["fermi_csv"] = Path(fermi_csv)
            sub_analyses.append("fermi")

    if compute_phi_z:
        phi_z_dir = pot_dir / "phi_z"
        phi_z_dir.mkdir(parents=True, exist_ok=True)
        phi_z_png = phi_z_planeavg_analysis(
            cube_pattern,
            output_dir=phi_z_dir,
            max_curves=max_curves,
            frame_start=frame_start,
            frame_end=frame_end,
            frame_step=frame_step,
            verbose=verbose,
            **_dist,
        )
        artifacts["phi_z_png"] = Path(phi_z_png)
        sub_analyses.append("phi_z")

    if has_fermi:
        ts_dir = pot_dir / "thickness_sensitivity"
        ts_dir.mkdir(parents=True, exist_ok=True)
        ts_csv = thickness_sensitivity_analysis(
            cube_pattern,
            md_out_p,
            output_dir=ts_dir,
            thickness_end=thickness_end,
            center_mode=center_mode,
            xyz_path=xyz_p,
            metal_elements=metal_elements,
            layer_tol_ang=layer_tol_ang,
            fermi_unit=fermi_unit,
            frame_start=frame_start,
            frame_end=frame_end,
            frame_step=frame_step,
            verbose=verbose,
            **_dist,
        )
        artifacts["thickness_sensitivity_csv"] = Path(ts_csv)
        sub_analyses.append("thickness_sensitivity")

    metadata: dict[str, Any] = {
        "input_mode": input_mode,
        "center_mode": center_mode,
        "has_fermi": has_fermi,
        "compute_u": compute_u,
        "compute_phi_z": compute_phi_z,
        "thickness_ang": thickness_ang,
        "thickness_end": thickness_end,
        "fermi_unit": fermi_unit,
        "sub_analyses": sub_analyses,
    }
    return WorkflowResult(
        name="potential_full",
        output_dir=pot_dir,
        artifacts=artifacts,
        metadata=metadata,
    )
