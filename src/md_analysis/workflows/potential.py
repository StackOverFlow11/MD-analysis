"""Potential analysis workflow facade.

Exposes the canonical programmatic entry points for the potential
analysis suite:

- :func:`run_potential_full`         -- composite (orchestrates all
  sub-analyses; sub-dirs are appended under ``output_dir``).
- :func:`run_center_potential`       -- single-step slab-averaged
  Hartree potential (CSV).
- :func:`run_fermi_energy`           -- single-step Fermi-energy time
  series (CSV).
- :func:`run_electrode_potential`    -- single-step electrode potential
  vs SHE (CSV).
- :func:`run_phi_z_profile`          -- single-step planar-averaged
  phi(z) overlay (PNG).
- :func:`run_thickness_sensitivity`  -- single-step thickness sweep
  (CSV).

Single-step facades write outputs directly into the caller-provided
``output_dir`` (no sub-dir append); composite ``run_potential_full``
appends ``<sub>/`` for each sub-analysis.

All facades support both continuous (cube glob + optional md.out +
xyz) and distributed (``input_mode="distributed"`` + ``sp_*``)
modes; the underlying business functions dispatch on ``input_mode``.
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


# ---------------------------------------------------------------------------
# Single-step facades (Phase 6.2)
# ---------------------------------------------------------------------------


def _potential_distributed_kwargs(
    *,
    input_mode: str,
    sp_root_dir: Path | str | None,
    sp_dir_pattern: str,
    sp_cube_filename: str,
    sp_out_filename: str,
) -> dict[str, Any]:
    """Pack the 5 ``input_mode`` / sp_* keyword arguments shared by every
    potential business function so the single-step facades can forward
    them verbatim without per-facade ``Path(...)`` boilerplate."""
    return {
        "input_mode": input_mode,
        "sp_root_dir": Path(sp_root_dir) if sp_root_dir is not None else None,
        "sp_dir_pattern": sp_dir_pattern,
        "sp_cube_filename": sp_cube_filename,
        "sp_out_filename": sp_out_filename,
    }


def run_center_potential(
    *,
    output_dir: Path | str,
    cube_pattern: str = "md-POTENTIAL-v_hartree-1_*.cube",
    xyz_path: Path | str | None = None,
    thickness_ang: float = DEFAULT_THICKNESS_ANG,
    center_mode: str = "interface",
    metal_elements: set[str] | None = None,
    layer_tol_ang: float = DEFAULT_LAYER_TOL_A,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    input_mode: str = "continuous",
    sp_root_dir: Path | str | None = None,
    sp_dir_pattern: str = "potential_t*_i*",
    sp_cube_filename: str = "sp_potential-v_hartree-1_0.cube",
    sp_out_filename: str = "sp.out",
    **kwargs: Any,
) -> WorkflowResult:
    """Run center-slab Hartree potential analysis (single-step).

    Wraps :func:`md_analysis.electrochemical.potential.center_slab_potential_analysis`
    and packs the resulting CSV path into a :class:`WorkflowResult`
    with the stable artifact key ``center_csv``.

    Supports both continuous (cube_pattern + xyz_path) and distributed
    (input_mode="distributed" + sp_*) modes; the underlying business
    function dispatches on input_mode.
    """
    from ..electrochemical.potential import center_slab_potential_analysis

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    xyz_p = Path(xyz_path) if xyz_path is not None else None

    logger.info("Starting center potential analysis: output_dir=%s", out_dir)

    csv_path = center_slab_potential_analysis(
        cube_pattern,
        output_dir=out_dir,
        thickness_ang=thickness_ang,
        center_mode=center_mode,
        xyz_path=xyz_p,
        metal_elements=metal_elements,
        layer_tol_ang=layer_tol_ang,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
        **_potential_distributed_kwargs(
            input_mode=input_mode,
            sp_root_dir=sp_root_dir,
            sp_dir_pattern=sp_dir_pattern,
            sp_cube_filename=sp_cube_filename,
            sp_out_filename=sp_out_filename,
        ),
        **kwargs,
    )

    return WorkflowResult(
        name="center_potential",
        output_dir=out_dir,
        artifacts={"center_csv": Path(csv_path)},
        metadata={
            "input_mode": input_mode,
            "center_mode": center_mode,
            "thickness_ang": thickness_ang,
            "layer_tol_ang": layer_tol_ang,
            "frame_start": frame_start,
            "frame_end": frame_end,
            "frame_step": frame_step,
        },
    )


def run_fermi_energy(
    *,
    output_dir: Path | str,
    md_out_path: Path | str | None = None,
    fermi_unit: str = "au",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    input_mode: str = "continuous",
    sp_root_dir: Path | str | None = None,
    sp_dir_pattern: str = "potential_t*_i*",
    sp_cube_filename: str = "sp_potential-v_hartree-1_0.cube",
    sp_out_filename: str = "sp.out",
    **kwargs: Any,
) -> WorkflowResult:
    """Run Fermi-energy time-series analysis (single-step).

    Wraps :func:`md_analysis.electrochemical.potential.fermi_energy_analysis`
    and packs the resulting CSV path under the stable artifact key
    ``fermi_csv``.  Continuous mode needs ``md_out_path``; distributed
    mode reads Fermi from each ``sp.out`` in the SP subdirs.
    """
    from ..electrochemical.potential import fermi_energy_analysis

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    md_out_p = Path(md_out_path) if md_out_path is not None else None

    logger.info("Starting fermi energy analysis: output_dir=%s", out_dir)

    csv_path = fermi_energy_analysis(
        md_out_p,
        output_dir=out_dir,
        fermi_unit=fermi_unit,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        **_potential_distributed_kwargs(
            input_mode=input_mode,
            sp_root_dir=sp_root_dir,
            sp_dir_pattern=sp_dir_pattern,
            sp_cube_filename=sp_cube_filename,
            sp_out_filename=sp_out_filename,
        ),
        **kwargs,
    )

    return WorkflowResult(
        name="fermi_energy",
        output_dir=out_dir,
        artifacts={"fermi_csv": Path(csv_path)},
        metadata={
            "input_mode": input_mode,
            "fermi_unit": fermi_unit,
            "frame_start": frame_start,
            "frame_end": frame_end,
            "frame_step": frame_step,
        },
    )


def run_electrode_potential(
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
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    input_mode: str = "continuous",
    sp_root_dir: Path | str | None = None,
    sp_dir_pattern: str = "potential_t*_i*",
    sp_cube_filename: str = "sp_potential-v_hartree-1_0.cube",
    sp_out_filename: str = "sp.out",
    **kwargs: Any,
) -> WorkflowResult:
    """Run electrode potential (U vs SHE) analysis (single-step).

    Wraps :func:`md_analysis.electrochemical.potential.electrode_potential_analysis`
    and packs the resulting CSV path under the stable artifact key
    ``electrode_csv``.  Combines slab-centered Hartree potential with
    Fermi level to compute U vs SHE per the cSHE convention.
    """
    from ..electrochemical.potential import electrode_potential_analysis

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    md_out_p = Path(md_out_path) if md_out_path is not None else None
    xyz_p = Path(xyz_path) if xyz_path is not None else None

    logger.info("Starting electrode potential analysis: output_dir=%s", out_dir)

    csv_path = electrode_potential_analysis(
        cube_pattern,
        md_out_p,
        output_dir=out_dir,
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
        **_potential_distributed_kwargs(
            input_mode=input_mode,
            sp_root_dir=sp_root_dir,
            sp_dir_pattern=sp_dir_pattern,
            sp_cube_filename=sp_cube_filename,
            sp_out_filename=sp_out_filename,
        ),
        **kwargs,
    )

    return WorkflowResult(
        name="electrode_potential",
        output_dir=out_dir,
        artifacts={"electrode_csv": Path(csv_path)},
        metadata={
            "input_mode": input_mode,
            "center_mode": center_mode,
            "thickness_ang": thickness_ang,
            "layer_tol_ang": layer_tol_ang,
            "fermi_unit": fermi_unit,
            "frame_start": frame_start,
            "frame_end": frame_end,
            "frame_step": frame_step,
        },
    )


def run_phi_z_profile(
    *,
    output_dir: Path | str,
    cube_pattern: str = "md-POTENTIAL-v_hartree-1_*.cube",
    max_curves: int = 0,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    input_mode: str = "continuous",
    sp_root_dir: Path | str | None = None,
    sp_dir_pattern: str = "potential_t*_i*",
    sp_cube_filename: str = "sp_potential-v_hartree-1_0.cube",
    sp_out_filename: str = "sp.out",
    **kwargs: Any,
) -> WorkflowResult:
    """Run planar-averaged phi(z) overlay analysis (single-step).

    Wraps :func:`md_analysis.electrochemical.potential.phi_z_planeavg_analysis`
    and packs the resulting PNG path under the stable artifact key
    ``phi_z_png`` (the only potential single-step facade returning an
    image rather than a CSV).
    """
    from ..electrochemical.potential import phi_z_planeavg_analysis

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Starting phi(z) profile analysis: output_dir=%s", out_dir)

    png_path = phi_z_planeavg_analysis(
        cube_pattern,
        output_dir=out_dir,
        max_curves=max_curves,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
        **_potential_distributed_kwargs(
            input_mode=input_mode,
            sp_root_dir=sp_root_dir,
            sp_dir_pattern=sp_dir_pattern,
            sp_cube_filename=sp_cube_filename,
            sp_out_filename=sp_out_filename,
        ),
        **kwargs,
    )

    return WorkflowResult(
        name="phi_z_profile",
        output_dir=out_dir,
        artifacts={"phi_z_png": Path(png_path)},
        metadata={
            "input_mode": input_mode,
            "max_curves": max_curves,
            "frame_start": frame_start,
            "frame_end": frame_end,
            "frame_step": frame_step,
        },
    )


def run_thickness_sensitivity(
    *,
    output_dir: Path | str,
    cube_pattern: str = "md-POTENTIAL-v_hartree-1_*.cube",
    md_out_path: Path | str | None = None,
    xyz_path: Path | str | None = None,
    thickness_end: float = 15.0,
    center_mode: str = "interface",
    metal_elements: set[str] | None = None,
    layer_tol_ang: float = DEFAULT_LAYER_TOL_A,
    fermi_unit: str = "au",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    input_mode: str = "continuous",
    sp_root_dir: Path | str | None = None,
    sp_dir_pattern: str = "potential_t*_i*",
    sp_cube_filename: str = "sp_potential-v_hartree-1_0.cube",
    sp_out_filename: str = "sp.out",
    **kwargs: Any,
) -> WorkflowResult:
    """Run thickness-sensitivity sweep analysis (single-step).

    Wraps :func:`md_analysis.electrochemical.potential.thickness_sensitivity_analysis`
    and packs the resulting CSV path under the stable artifact key
    ``thickness_sensitivity_csv``.  Sweeps slab thickness up to
    ``thickness_end`` (default 15.0 A) and records mean U vs SHE +
    spatial std of phi(z) inside the slab at each thickness.
    """
    from ..electrochemical.potential import thickness_sensitivity_analysis

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    md_out_p = Path(md_out_path) if md_out_path is not None else None
    xyz_p = Path(xyz_path) if xyz_path is not None else None

    logger.info("Starting thickness sensitivity analysis: output_dir=%s", out_dir)

    csv_path = thickness_sensitivity_analysis(
        cube_pattern,
        md_out_p,
        output_dir=out_dir,
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
        **_potential_distributed_kwargs(
            input_mode=input_mode,
            sp_root_dir=sp_root_dir,
            sp_dir_pattern=sp_dir_pattern,
            sp_cube_filename=sp_cube_filename,
            sp_out_filename=sp_out_filename,
        ),
        **kwargs,
    )

    return WorkflowResult(
        name="thickness_sensitivity",
        output_dir=out_dir,
        artifacts={"thickness_sensitivity_csv": Path(csv_path)},
        metadata={
            "input_mode": input_mode,
            "center_mode": center_mode,
            "thickness_end": thickness_end,
            "layer_tol_ang": layer_tol_ang,
            "fermi_unit": fermi_unit,
            "frame_start": frame_start,
            "frame_end": frame_end,
            "frame_step": frame_step,
        },
    )
