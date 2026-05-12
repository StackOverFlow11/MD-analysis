"""Programmatic entry points for running analysis workflows.

These ``run_*_analysis`` names are kept as backwards-compatibility
shims around the canonical workflow facades in
:mod:`md_analysis.workflows`. They return ``dict[str, Path]`` rather
than :class:`~md_analysis.workflows.WorkflowResult`. New code should
prefer the explicit ``md_analysis.workflows.*`` imports.

Public API
----------
- ``run_water_analysis``        — shim → ``workflows.water.run_water_three_panel``
- ``run_potential_analysis``    — shim → ``workflows.potential.run_potential_full``
- ``run_charge_analysis``       — shim → ``workflows.charge.run_surface_charge``
- ``run_tracked_charge_analysis``   — shim → ``workflows.charge.run_tracked_charge``
- ``run_counterion_charge_analysis``— shim → ``workflows.charge.run_counterion_charge``
- ``run_all``                   — deprecated shim → ``workflows.composite.run_interface_analysis``
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Iterable

from .electrochemical.potential.config import DEFAULT_THICKNESS_ANG
from .utils.constants import CHARGE_METHOD_COUNTERION, DEFAULT_LAYER_TOL_A

logger = logging.getLogger(__name__)


def run_water_analysis(
    xyz_path: Path,
    md_inp_path: Path | None = None,
    *,
    cell_abc: tuple[float, float, float] | None = None,
    output_dir: Path,
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    **kwargs: Any,
) -> dict[str, Path]:
    """Run water analysis (three-panel plot + CSVs).

    Thin shim around
    :func:`md_analysis.workflows.water.run_water_three_panel`; returns
    ``dict[str, Path]`` for backwards compatibility. New code should
    prefer the workflow function, which returns ``WorkflowResult``.
    """
    from .workflows.water import run_water_three_panel

    result = run_water_three_panel(
        xyz_path,
        md_inp_path,
        cell_abc=cell_abc,
        output_dir=output_dir,
        layer_tol_A=layer_tol_A,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
        **kwargs,
    )
    return dict(result.artifacts)


def run_potential_analysis(
    *,
    output_dir: Path,
    cube_pattern: str = "md-POTENTIAL-v_hartree-1_*.cube",
    md_out_path: Path | None = None,
    xyz_path: Path | None = None,
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
) -> dict[str, Path]:
    """Run all potential analysis workflows.

    Thin shim around
    :func:`md_analysis.workflows.potential.run_potential_full`; returns
    ``dict[str, Path]`` for backwards compatibility. New code should
    prefer the workflow function, which returns ``WorkflowResult``.
    """
    from .workflows.potential import run_potential_full

    result = run_potential_full(
        output_dir=output_dir,
        cube_pattern=cube_pattern,
        md_out_path=md_out_path,
        xyz_path=xyz_path,
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
        input_mode=input_mode,
        sp_root_dir=sp_root_dir,
        sp_dir_pattern=sp_dir_pattern,
        sp_cube_filename=sp_cube_filename,
        sp_out_filename=sp_out_filename,
    )
    return dict(result.artifacts)


def run_charge_analysis(
    *,
    output_dir: Path,
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
) -> dict[str, Path]:
    """Run surface charge density analysis (CSV + PNG).

    Thin shim around
    :func:`md_analysis.workflows.charge.run_surface_charge`; returns
    ``dict[str, Path]`` for backwards compatibility. New code should
    prefer the workflow function, which returns ``WorkflowResult``.
    """
    from .workflows.charge import run_surface_charge

    result = run_surface_charge(
        output_dir=output_dir,
        root_dir=root_dir,
        metal_symbols=metal_symbols,
        normal=normal,
        method=method,
        layer_tol_A=layer_tol_A,
        n_surface_layers=n_surface_layers,
        dir_pattern=dir_pattern,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )
    return dict(result.artifacts)


def run_tracked_charge_analysis(
    *,
    output_dir: Path,
    root_dir: str | Path = ".",
    atom_indices_xyz: Iterable[int],
    dir_pattern: str = "bader_t*_i*",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> dict[str, Path]:
    """Track Bader net charges for specified XYZ atoms (CSV + PNG).

    Thin shim around
    :func:`md_analysis.workflows.charge.run_tracked_charge`; returns
    ``dict[str, Path]`` for backwards compatibility. New code should
    prefer the workflow function, which returns ``WorkflowResult``.
    """
    from .workflows.charge import run_tracked_charge

    result = run_tracked_charge(
        output_dir=output_dir,
        root_dir=root_dir,
        atom_indices_xyz=atom_indices_xyz,
        dir_pattern=dir_pattern,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )
    return dict(result.artifacts)


def run_counterion_charge_analysis(
    *,
    output_dir: Path,
    root_dir: str | Path = ".",
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
    layer_tol_A: float = DEFAULT_LAYER_TOL_A,
    dir_pattern: str = "bader_t*_i*",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> dict[str, Path]:
    """Detect counterions per-frame and track their Bader charges (CSV + PNG).

    Thin shim around
    :func:`md_analysis.workflows.charge.run_counterion_charge`; returns
    ``dict[str, Path]`` for backwards compatibility. New code should
    prefer the workflow function, which returns ``WorkflowResult``.
    """
    from .workflows.charge import run_counterion_charge

    result = run_counterion_charge(
        output_dir=output_dir,
        root_dir=root_dir,
        metal_symbols=metal_symbols,
        normal=normal,
        layer_tol_A=layer_tol_A,
        dir_pattern=dir_pattern,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )
    return dict(result.artifacts)



_RUN_ALL_POTENTIAL_KWARGS = frozenset({
    "thickness_ang", "center_mode", "metal_elements",
    "layer_tol_ang", "fermi_unit", "compute_u",
    "compute_phi_z", "max_curves", "thickness_end",
})


def run_all(
    xyz_path: Path,
    md_inp_path: Path | None = None,
    *,
    cell_abc: tuple[float, float, float] | None = None,
    output_dir: Path,
    cube_pattern: str = "md-POTENTIAL-v_hartree-1_*.cube",
    md_out_path: Path | None = None,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    **kwargs: Any,
) -> dict[str, Path]:
    """Run water + potential under the canonical interface layout.

    Deprecated thin shim around
    :func:`md_analysis.workflows.composite.run_interface_analysis`;
    returns ``dict[str, Path]`` for backwards compatibility. The name
    ``run_all`` is retained only to keep existing callers (CLI menus,
    integration tests, and documentation references) working until
    Phase 6 (CLI rewire) and Phase 8 (documentation + integration-test
    sweep) finish. New code should prefer
    ``md_analysis.workflows.run_interface_analysis``, which returns
    a :class:`WorkflowResult`.

    The ``**kwargs`` whitelist mirrors the pre-refactor behaviour:
    only the keys in :data:`_RUN_ALL_POTENTIAL_KWARGS` are forwarded
    to the potential leaf; anything else is silently dropped, as it
    was before.
    """
    from .workflows.composite import run_interface_analysis

    pot_kwargs = {
        k: v for k, v in kwargs.items() if k in _RUN_ALL_POTENTIAL_KWARGS
    }
    result = run_interface_analysis(
        xyz_path,
        md_inp_path,
        cell_abc=cell_abc,
        output_dir=output_dir,
        cube_pattern=cube_pattern,
        md_out_path=md_out_path,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
        **pot_kwargs,
    )
    return dict(result.artifacts)
