"""Programmatic entry points for running analysis workflows.

Public API
----------
- ``run_water_analysis``  — water density/orientation/adsorbed-layer analysis
- ``run_potential_analysis`` — Hartree potential / Fermi / electrode potential analysis
- ``run_charge_analysis`` — Bader surface charge density time series
- ``run_all`` — both water + potential
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

from .potential.config import DEFAULT_THICKNESS_ANG


def run_water_analysis(
    xyz_path: Path,
    md_inp_path: Path,
    *,
    output_dir: Path,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    **kwargs: Any,
) -> dict[str, Path]:
    """Run water analysis (three-panel plot + CSVs).

    Outputs are written under ``output_dir/water/`` with sub-directories
    for density, orientation, and adsorbed-layer results.

    Returns a dict mapping output names to file paths.
    """
    from .water import plot_water_three_panel_analysis
    from .water.config import (
        DEFAULT_WATER_MASS_DENSITY_CSV_NAME,
        DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME,
        DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME,
        DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME,
        DEFAULT_ADSORBED_WATER_THETA_DISTRIBUTION_CSV_NAME,
        DEFAULT_WATER_THREE_PANEL_PLOT_PNG_NAME,
    )

    water_dir = Path(output_dir) / "water"
    water_dir.mkdir(parents=True, exist_ok=True)

    png_path = plot_water_three_panel_analysis(
        xyz_path=xyz_path,
        md_inp_path=md_inp_path,
        output_dir=water_dir,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
        **kwargs,
    )

    return {
        "density_csv": water_dir / DEFAULT_WATER_MASS_DENSITY_CSV_NAME,
        "orientation_csv": water_dir / DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME,
        "adsorbed_profile_csv": water_dir / DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME,
        "adsorbed_range_txt": water_dir / DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME,
        "adsorbed_theta_csv": water_dir / DEFAULT_ADSORBED_WATER_THETA_DISTRIBUTION_CSV_NAME,
        "plot_png": png_path,
    }


def run_potential_analysis(
    *,
    output_dir: Path,
    cube_pattern: str = "md-POTENTIAL-v_hartree-1_*.cube",
    md_out_path: Path | None = None,
    xyz_path: Path | None = None,
    thickness_ang: float = DEFAULT_THICKNESS_ANG,
    center_mode: str = "interface",
    metal_elements: set[str] | None = None,
    layer_tol_ang: float = 0.6,
    fermi_unit: str = "au",
    compute_u: bool = True,
    compute_phi_z: bool = True,
    max_curves: int = 0,
    thickness_end: float = 15.0,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> dict[str, Path]:
    """Run all potential analysis workflows.

    Outputs are written under ``output_dir/potential/`` with sub-directories
    for center, fermi, electrode, and phi_z results.

    Returns a dict mapping output names to file paths.
    """
    from .potential import (
        center_slab_potential_analysis,
        fermi_energy_analysis,
        electrode_potential_analysis,
        phi_z_planeavg_analysis,
        thickness_sensitivity_analysis,
    )

    pot_dir = Path(output_dir) / "potential"
    results: dict[str, Path] = {}

    if compute_u and md_out_path is not None:
        # Full electrode potential analysis (includes center + fermi)
        electrode_dir = pot_dir / "electrode"
        electrode_dir.mkdir(parents=True, exist_ok=True)
        u_csv = electrode_potential_analysis(
            cube_pattern,
            md_out_path,
            output_dir=electrode_dir,
            thickness_ang=thickness_ang,
            center_mode=center_mode,
            xyz_path=xyz_path,
            metal_elements=metal_elements,
            layer_tol_ang=layer_tol_ang,
            fermi_unit=fermi_unit,
            frame_start=frame_start,
            frame_end=frame_end,
            frame_step=frame_step,
            verbose=verbose,
        )
        results["electrode_csv"] = u_csv
    else:
        # Run individually
        center_dir = pot_dir / "center"
        center_dir.mkdir(parents=True, exist_ok=True)
        center_csv = center_slab_potential_analysis(
            cube_pattern,
            output_dir=center_dir,
            thickness_ang=thickness_ang,
            center_mode=center_mode,
            xyz_path=xyz_path,
            metal_elements=metal_elements,
            layer_tol_ang=layer_tol_ang,
            frame_start=frame_start,
            frame_end=frame_end,
            frame_step=frame_step,
            verbose=verbose,
        )
        results["center_csv"] = center_csv

        if md_out_path is not None:
            fermi_dir = pot_dir / "fermi"
            fermi_dir.mkdir(parents=True, exist_ok=True)
            fermi_csv = fermi_energy_analysis(
                md_out_path,
                output_dir=fermi_dir,
                fermi_unit=fermi_unit,
                frame_start=frame_start,
                frame_end=frame_end,
                frame_step=frame_step,
            )
            results["fermi_csv"] = fermi_csv

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
        )
        results["phi_z_png"] = phi_z_png

    # Thickness sensitivity sweep (requires Fermi energies from md.out)
    if md_out_path is not None:
        ts_dir = pot_dir / "thickness_sensitivity"
        ts_dir.mkdir(parents=True, exist_ok=True)
        ts_csv = thickness_sensitivity_analysis(
            cube_pattern,
            md_out_path,
            output_dir=ts_dir,
            thickness_end=thickness_end,
            center_mode=center_mode,
            xyz_path=xyz_path,
            metal_elements=metal_elements,
            layer_tol_ang=layer_tol_ang,
            fermi_unit=fermi_unit,
            frame_start=frame_start,
            frame_end=frame_end,
            frame_step=frame_step,
            verbose=verbose,
        )
        results["thickness_sensitivity_csv"] = ts_csv

    return results


def run_charge_analysis(
    *,
    output_dir: Path,
    root_dir: str | Path = ".",
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
    dir_pattern: str = "calc_t*_i*",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> dict[str, Path]:
    """Run surface charge density analysis (CSV + PNG).

    Outputs are written under ``output_dir/charge/``.
    Returns a dict mapping output names to file paths.
    """
    from .charge import surface_charge_analysis
    from .charge.config import (
        DEFAULT_SURFACE_CHARGE_CSV_NAME,
        DEFAULT_SURFACE_CHARGE_PNG_NAME,
    )

    charge_dir = Path(output_dir) / "charge"
    charge_dir.mkdir(parents=True, exist_ok=True)

    csv_path = surface_charge_analysis(
        root_dir,
        metal_symbols=metal_symbols,
        normal=normal,
        dir_pattern=dir_pattern,
        output_dir=charge_dir,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )

    return {
        "charge_csv": csv_path,
        "charge_png": csv_path.parent / DEFAULT_SURFACE_CHARGE_PNG_NAME,
    }


def run_all(
    xyz_path: Path,
    md_inp_path: Path,
    *,
    output_dir: Path,
    cube_pattern: str = "md-POTENTIAL-v_hartree-1_*.cube",
    md_out_path: Path | None = None,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    **kwargs: Any,
) -> dict[str, Path]:
    """Run all analysis workflows (water + potential).

    Returns a merged dict of all output file paths.
    """
    results: dict[str, Path] = {}

    results.update(run_water_analysis(
        xyz_path, md_inp_path, output_dir=output_dir,
        frame_start=frame_start, frame_end=frame_end, frame_step=frame_step,
        verbose=verbose,
    ))

    pot_kwargs = {
        k: v for k, v in kwargs.items()
        if k in {
            "thickness_ang", "center_mode", "metal_elements",
            "layer_tol_ang", "fermi_unit", "compute_u",
            "compute_phi_z", "max_curves", "thickness_end",
        }
    }
    results.update(run_potential_analysis(
        output_dir=output_dir,
        cube_pattern=cube_pattern,
        md_out_path=md_out_path,
        xyz_path=xyz_path,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
        **pot_kwargs,
    ))

    return results
