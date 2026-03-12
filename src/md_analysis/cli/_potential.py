"""Potential analysis command classes (211-216)."""

from __future__ import annotations

from pathlib import Path

from ._framework import MenuCommand, lazy_import
from ._params import (
    K,
    center_mode,
    cube_pattern,
    fermi_unit,
    frame_slice,
    layer_tol,
    max_curves,
    md_out_path,
    metal_elements,
    outdir,
    thickness,
    thickness_end,
    xyz_path_opt,
)


class CenterPotentialCmd(MenuCommand):
    params = (cube_pattern, xyz_path_opt, thickness, center_mode,
              metal_elements, layer_tol)
    advanced_params = (outdir, frame_slice)
    output_subdir = "potential/center"

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import("md_analysis.electrochemical.potential",
                              "center_slab_potential_analysis")
        csv = analyze(
            ctx[K.CUBE_PATTERN],
            output_dir=ctx[K.OUTDIR_RESOLVED],
            thickness_ang=ctx[K.THICKNESS],
            center_mode=ctx[K.CENTER_MODE],
            xyz_path=Path(ctx[K.XYZ]),
            metal_elements=ctx[K.METAL_ELEMENTS],
            layer_tol_ang=ctx[K.LAYER_TOL],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
            verbose=True,
        )
        print(f"\n Analysis complete. Output:\n   center_csv: {csv}")


class FermiEnergyCmd(MenuCommand):
    params = (md_out_path, fermi_unit)
    advanced_params = (outdir, frame_slice)
    output_subdir = "potential/fermi"

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import("md_analysis.electrochemical.potential",
                              "fermi_energy_analysis")
        csv = analyze(
            Path(ctx[K.MD_OUT]),
            output_dir=ctx[K.OUTDIR_RESOLVED],
            fermi_unit=ctx[K.FERMI_UNIT],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
        )
        print(f"\n Analysis complete. Output:\n   fermi_csv: {csv}")


class ElectrodePotentialCmd(MenuCommand):
    params = (cube_pattern, xyz_path_opt, thickness, center_mode,
              metal_elements, layer_tol, md_out_path, fermi_unit)
    advanced_params = (outdir, frame_slice)
    output_subdir = "potential/electrode"

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import("md_analysis.electrochemical.potential",
                              "electrode_potential_analysis")
        csv = analyze(
            ctx[K.CUBE_PATTERN],
            Path(ctx[K.MD_OUT]),
            output_dir=ctx[K.OUTDIR_RESOLVED],
            thickness_ang=ctx[K.THICKNESS],
            center_mode=ctx[K.CENTER_MODE],
            xyz_path=Path(ctx[K.XYZ]),
            metal_elements=ctx[K.METAL_ELEMENTS],
            layer_tol_ang=ctx[K.LAYER_TOL],
            fermi_unit=ctx[K.FERMI_UNIT],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
            verbose=True,
        )
        print(f"\n Analysis complete. Output:\n   electrode_csv: {csv}")


class PhiZProfileCmd(MenuCommand):
    params = (cube_pattern, max_curves)
    advanced_params = (outdir, frame_slice)
    output_subdir = "potential/phi_z"

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import("md_analysis.electrochemical.potential",
                              "phi_z_planeavg_analysis")
        png = analyze(
            ctx[K.CUBE_PATTERN],
            output_dir=ctx[K.OUTDIR_RESOLVED],
            max_curves=ctx[K.MAX_CURVES],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
            verbose=True,
        )
        print(f"\n Analysis complete. Output:\n   phi_z_png: {png}")


class ThicknessSensitivityCmd(MenuCommand):
    params = (cube_pattern, xyz_path_opt, thickness, center_mode,
              metal_elements, layer_tol, md_out_path, fermi_unit,
              thickness_end)
    advanced_params = (outdir, frame_slice)
    output_subdir = "potential/thickness_sensitivity"

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import("md_analysis.electrochemical.potential",
                              "thickness_sensitivity_analysis")
        csv = analyze(
            ctx[K.CUBE_PATTERN],
            Path(ctx[K.MD_OUT]),
            output_dir=ctx[K.OUTDIR_RESOLVED],
            thickness_end=ctx[K.THICKNESS_END],
            center_mode=ctx[K.CENTER_MODE],
            xyz_path=Path(ctx[K.XYZ]),
            metal_elements=ctx[K.METAL_ELEMENTS],
            layer_tol_ang=ctx[K.LAYER_TOL],
            fermi_unit=ctx[K.FERMI_UNIT],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
            verbose=True,
        )
        print(f"\n Analysis complete. Output:\n   thickness_sensitivity_csv: {csv}")


class FullPotentialCmd(MenuCommand):
    params = (cube_pattern, xyz_path_opt, thickness, center_mode,
              metal_elements, layer_tol, md_out_path, fermi_unit,
              max_curves, thickness_end)
    advanced_params = (outdir, frame_slice)
    output_subdir = ""

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import("md_analysis.main", "run_potential_analysis")
        md_out = Path(ctx[K.MD_OUT])
        if not md_out.exists():
            print(f"  WARNING: {md_out} not found, "
                  "skipping Fermi/electrode/sensitivity analyses.")
            md_out = None

        results = analyze(
            output_dir=Path(ctx[K.OUTDIR]),
            cube_pattern=ctx[K.CUBE_PATTERN],
            md_out_path=md_out,
            xyz_path=Path(ctx[K.XYZ]),
            thickness_ang=ctx[K.THICKNESS],
            center_mode=ctx[K.CENTER_MODE],
            metal_elements=ctx[K.METAL_ELEMENTS],
            layer_tol_ang=ctx[K.LAYER_TOL],
            fermi_unit=ctx[K.FERMI_UNIT],
            max_curves=ctx[K.MAX_CURVES],
            thickness_end=ctx[K.THICKNESS_END],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
            verbose=True,
        )
        print("\n Analysis complete. Outputs:")
        for name, path in results.items():
            print(f"   {name}: {path}")
