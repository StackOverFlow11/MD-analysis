"""Water analysis command classes (101-105)."""

from __future__ import annotations

from pathlib import Path

from ._framework import MenuCommand, lazy_import
from ._params import (
    K,
    cell_abc,
    dz_bin,
    frame_slice,
    layer_tol,
    outdir,
    xyz_path,
)

_WATER_PARAMS = (xyz_path, cell_abc, dz_bin)
_WATER_ADVANCED = (layer_tol, outdir, frame_slice)


class WaterDensityCmd(MenuCommand):
    params = _WATER_PARAMS
    advanced_params = _WATER_ADVANCED
    output_subdir = "water"

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import("md_analysis.water",
                              "water_mass_density_z_distribution_analysis")
        csv = analyze(
            xyz_path=ctx[K.XYZ],
            cell_abc=ctx[K.CELL_ABC],
            output_dir=ctx[K.OUTDIR_RESOLVED],
            dz_A=ctx[K.DZ_A],
            layer_tol_A=ctx[K.LAYER_TOL],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
        )
        print(f"\n Analysis complete. Output:\n   density_csv: {csv}")


class WaterOrientationCmd(MenuCommand):
    params = _WATER_PARAMS
    advanced_params = _WATER_ADVANCED
    output_subdir = "water"

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import(
            "md_analysis.water",
            "water_orientation_weighted_density_z_distribution_analysis",
        )
        csv = analyze(
            xyz_path=ctx[K.XYZ],
            cell_abc=ctx[K.CELL_ABC],
            output_dir=ctx[K.OUTDIR_RESOLVED],
            dz_A=ctx[K.DZ_A],
            layer_tol_A=ctx[K.LAYER_TOL],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
        )
        print(f"\n Analysis complete. Output:\n   orientation_csv: {csv}")


class AdWaterOrientationCmd(MenuCommand):
    params = _WATER_PARAMS
    advanced_params = _WATER_ADVANCED
    output_subdir = "water"

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import("md_analysis.water",
                              "ad_water_orientation_analysis")
        profile_csv, range_txt = analyze(
            xyz_path=ctx[K.XYZ],
            cell_abc=ctx[K.CELL_ABC],
            output_dir=ctx[K.OUTDIR_RESOLVED],
            dz_A=ctx[K.DZ_A],
            layer_tol_A=ctx[K.LAYER_TOL],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
        )
        print("\n Analysis complete. Outputs:")
        print(f"   adsorbed_profile_csv: {profile_csv}")
        print(f"   adsorbed_range_txt:   {range_txt}")


class AdWaterThetaCmd(MenuCommand):
    params = _WATER_PARAMS
    advanced_params = _WATER_ADVANCED
    output_subdir = "water"

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import("md_analysis.water",
                              "compute_adsorbed_water_theta_distribution")
        _, _, csv = analyze(
            xyz_path=ctx[K.XYZ],
            cell_abc=ctx[K.CELL_ABC],
            output_dir=ctx[K.OUTDIR_RESOLVED],
            dz_A=ctx[K.DZ_A],
            layer_tol_A=ctx[K.LAYER_TOL],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
            verbose=True,
        )
        print(f"\n Analysis complete. Output:\n   theta_csv: {csv}")


class WaterThreePanelCmd(MenuCommand):
    params = _WATER_PARAMS
    advanced_params = _WATER_ADVANCED
    output_subdir = ""

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import("md_analysis.main", "run_water_analysis")
        results = analyze(
            xyz_path=Path(ctx[K.XYZ]),
            cell_abc=ctx[K.CELL_ABC],
            output_dir=Path(ctx[K.OUTDIR]),
            dz_A=ctx[K.DZ_A],
            layer_tol_A=ctx[K.LAYER_TOL],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
            verbose=True,
        )
        print("\n Analysis complete. Outputs:")
        for name, path in results.items():
            print(f"   {name}: {path}")
