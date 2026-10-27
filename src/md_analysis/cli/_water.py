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
    # output_name inherited from parent MenuGroup("1", output_name="water")

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import(
            "md_analysis.workflows.water", "run_water_density",
        )
        result = analyze(
            xyz_path=Path(ctx[K.XYZ]),
            cell_abc=ctx[K.CELL_ABC],
            output_dir=ctx[K.OUTDIR_RESOLVED],
            dz_A=ctx[K.DZ_A],
            layer_tol_A=ctx[K.LAYER_TOL],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
        )
        print("\n Analysis complete. Outputs:")
        for name, path in result.artifacts.items():
            print(f"   {name}: {path}")


class WaterOrientationCmd(MenuCommand):
    params = _WATER_PARAMS
    advanced_params = _WATER_ADVANCED
    # output_name inherited from parent MenuGroup("1", output_name="water")

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import(
            "md_analysis.workflows.water", "run_water_orientation",
        )
        result = analyze(
            xyz_path=Path(ctx[K.XYZ]),
            cell_abc=ctx[K.CELL_ABC],
            output_dir=ctx[K.OUTDIR_RESOLVED],
            dz_A=ctx[K.DZ_A],
            layer_tol_A=ctx[K.LAYER_TOL],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
        )
        print("\n Analysis complete. Outputs:")
        for name, path in result.artifacts.items():
            print(f"   {name}: {path}")


class AdWaterOrientationCmd(MenuCommand):
    params = _WATER_PARAMS
    advanced_params = _WATER_ADVANCED
    # output_name inherited from parent MenuGroup("1", output_name="water")

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import(
            "md_analysis.workflows.water", "run_ad_water_orientation",
        )
        result = analyze(
            xyz_path=Path(ctx[K.XYZ]),
            cell_abc=ctx[K.CELL_ABC],
            output_dir=ctx[K.OUTDIR_RESOLVED],
            dz_A=ctx[K.DZ_A],
            layer_tol_A=ctx[K.LAYER_TOL],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
        )
        print("\n Analysis complete. Outputs:")
        for name, path in result.artifacts.items():
            print(f"   {name}: {path}")


class AdWaterThetaCmd(MenuCommand):
    params = _WATER_PARAMS
    advanced_params = _WATER_ADVANCED
    # output_name inherited from parent MenuGroup("1", output_name="water")

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import(
            "md_analysis.workflows.water", "run_ad_water_theta",
        )
        result = analyze(
            xyz_path=Path(ctx[K.XYZ]),
            cell_abc=ctx[K.CELL_ABC],
            output_dir=ctx[K.OUTDIR_RESOLVED],
            dz_A=ctx[K.DZ_A],
            layer_tol_A=ctx[K.LAYER_TOL],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
            verbose=True,
        )
        print("\n Analysis complete. Outputs:")
        for name, path in result.artifacts.items():
            print(f"   {name}: {path}")


class WaterThreePanelCmd(MenuCommand):
    params = _WATER_PARAMS
    advanced_params = _WATER_ADVANCED
    # output_name inherited from parent MenuGroup("1", output_name="water")

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import(
            "md_analysis.workflows.water", "run_water_three_panel",
        )
        result = analyze(
            xyz_path=Path(ctx[K.XYZ]),
            cell_abc=ctx[K.CELL_ABC],
            output_dir=ctx[K.OUTDIR_RESOLVED],
            dz_A=ctx[K.DZ_A],
            layer_tol_A=ctx[K.LAYER_TOL],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
            verbose=True,
        )
        print("\n Analysis complete. Outputs:")
        for name, path in result.artifacts.items():
            print(f"   {name}: {path}")
