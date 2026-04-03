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
    input_mode,
    layer_tol,
    max_curves,
    md_out_path,
    metal_elements,
    outdir,
    sp_cube_filename,
    sp_dir_pattern,
    sp_out_filename,
    sp_root_dir,
    thickness,
    thickness_end,
    xyz_path_opt,
)


def _distributed_kwargs(ctx: dict) -> dict:
    """Extract distributed-mode keyword arguments from context."""
    return {
        "input_mode": ctx[K.INPUT_MODE],
        "sp_root_dir": ctx[K.SP_ROOT_DIR],
        "sp_dir_pattern": ctx[K.SP_DIR_PATTERN],
        "sp_cube_filename": ctx[K.SP_CUBE_FILENAME],
        "sp_out_filename": ctx[K.SP_OUT_FILENAME],
    }


def _is_distributed(ctx: dict) -> bool:
    return ctx.get(K.INPUT_MODE) == "distributed"


class CenterPotentialCmd(MenuCommand):
    output_name = "center"
    params = (input_mode, cube_pattern, xyz_path_opt, thickness, center_mode,
              metal_elements, layer_tol,
              sp_root_dir, sp_dir_pattern, sp_cube_filename, sp_out_filename)
    advanced_params = (outdir, frame_slice)

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import("md_analysis.electrochemical.potential",
                              "center_slab_potential_analysis")
        kwargs: dict = {
            "output_dir": ctx[K.OUTDIR_RESOLVED],
            "thickness_ang": ctx[K.THICKNESS],
            "center_mode": ctx[K.CENTER_MODE],
            "metal_elements": ctx[K.METAL_ELEMENTS],
            "layer_tol_ang": ctx[K.LAYER_TOL],
            "frame_start": ctx[K.FRAME_START],
            "frame_end": ctx[K.FRAME_END],
            "frame_step": ctx[K.FRAME_STEP],
            "verbose": True,
            **_distributed_kwargs(ctx),
        }
        if not _is_distributed(ctx):
            kwargs["cube_pattern"] = ctx[K.CUBE_PATTERN]
            kwargs["xyz_path"] = Path(ctx[K.XYZ])
        csv = analyze(**kwargs)
        print(f"\n Analysis complete. Output:\n   center_csv: {csv}")


class FermiEnergyCmd(MenuCommand):
    output_name = "fermi"
    params = (input_mode, md_out_path, fermi_unit,
              sp_root_dir, sp_dir_pattern, sp_cube_filename, sp_out_filename)
    advanced_params = (outdir, frame_slice)

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import("md_analysis.electrochemical.potential",
                              "fermi_energy_analysis")
        kwargs: dict = {
            "output_dir": ctx[K.OUTDIR_RESOLVED],
            "fermi_unit": ctx[K.FERMI_UNIT],
            "frame_start": ctx[K.FRAME_START],
            "frame_end": ctx[K.FRAME_END],
            "frame_step": ctx[K.FRAME_STEP],
            **_distributed_kwargs(ctx),
        }
        if not _is_distributed(ctx):
            kwargs["md_out_path"] = Path(ctx[K.MD_OUT])
        csv = analyze(**kwargs)
        print(f"\n Analysis complete. Output:\n   fermi_csv: {csv}")


class ElectrodePotentialCmd(MenuCommand):
    output_name = "electrode"
    params = (input_mode, cube_pattern, xyz_path_opt, thickness, center_mode,
              metal_elements, layer_tol, md_out_path, fermi_unit,
              sp_root_dir, sp_dir_pattern, sp_cube_filename, sp_out_filename)
    advanced_params = (outdir, frame_slice)

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import("md_analysis.electrochemical.potential",
                              "electrode_potential_analysis")
        kwargs: dict = {
            "output_dir": ctx[K.OUTDIR_RESOLVED],
            "thickness_ang": ctx[K.THICKNESS],
            "center_mode": ctx[K.CENTER_MODE],
            "metal_elements": ctx[K.METAL_ELEMENTS],
            "layer_tol_ang": ctx[K.LAYER_TOL],
            "fermi_unit": ctx[K.FERMI_UNIT],
            "frame_start": ctx[K.FRAME_START],
            "frame_end": ctx[K.FRAME_END],
            "frame_step": ctx[K.FRAME_STEP],
            "verbose": True,
            **_distributed_kwargs(ctx),
        }
        if not _is_distributed(ctx):
            kwargs["cube_pattern"] = ctx[K.CUBE_PATTERN]
            kwargs["md_out_path"] = Path(ctx[K.MD_OUT])
            kwargs["xyz_path"] = Path(ctx[K.XYZ])
        csv = analyze(**kwargs)
        print(f"\n Analysis complete. Output:\n   electrode_csv: {csv}")


class PhiZProfileCmd(MenuCommand):
    output_name = "phi_z"
    params = (input_mode, cube_pattern, max_curves,
              sp_root_dir, sp_dir_pattern, sp_cube_filename, sp_out_filename)
    advanced_params = (outdir, frame_slice)

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import("md_analysis.electrochemical.potential",
                              "phi_z_planeavg_analysis")
        kwargs: dict = {
            "output_dir": ctx[K.OUTDIR_RESOLVED],
            "max_curves": ctx[K.MAX_CURVES],
            "frame_start": ctx[K.FRAME_START],
            "frame_end": ctx[K.FRAME_END],
            "frame_step": ctx[K.FRAME_STEP],
            "verbose": True,
            **_distributed_kwargs(ctx),
        }
        if not _is_distributed(ctx):
            kwargs["cube_pattern"] = ctx[K.CUBE_PATTERN]
        png = analyze(**kwargs)
        print(f"\n Analysis complete. Output:\n   phi_z_png: {png}")


class ThicknessSensitivityCmd(MenuCommand):
    output_name = "thickness_sensitivity"
    params = (input_mode, cube_pattern, xyz_path_opt, thickness, center_mode,
              metal_elements, layer_tol, md_out_path, fermi_unit,
              thickness_end,
              sp_root_dir, sp_dir_pattern, sp_cube_filename, sp_out_filename)
    advanced_params = (outdir, frame_slice)

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import("md_analysis.electrochemical.potential",
                              "thickness_sensitivity_analysis")
        kwargs: dict = {
            "output_dir": ctx[K.OUTDIR_RESOLVED],
            "thickness_end": ctx[K.THICKNESS_END],
            "center_mode": ctx[K.CENTER_MODE],
            "metal_elements": ctx[K.METAL_ELEMENTS],
            "layer_tol_ang": ctx[K.LAYER_TOL],
            "fermi_unit": ctx[K.FERMI_UNIT],
            "frame_start": ctx[K.FRAME_START],
            "frame_end": ctx[K.FRAME_END],
            "frame_step": ctx[K.FRAME_STEP],
            "verbose": True,
            **_distributed_kwargs(ctx),
        }
        if not _is_distributed(ctx):
            kwargs["cube_pattern"] = ctx[K.CUBE_PATTERN]
            kwargs["md_out_path"] = Path(ctx[K.MD_OUT])
            kwargs["xyz_path"] = Path(ctx[K.XYZ])
        csv = analyze(**kwargs)
        print(f"\n Analysis complete. Output:\n   thickness_sensitivity_csv: {csv}")


class FullPotentialCmd(MenuCommand):
    # output_name not set → inherits "electrochemical/potential" from parent chain
    params = (input_mode, cube_pattern, xyz_path_opt, thickness, center_mode,
              metal_elements, layer_tol, md_out_path, fermi_unit,
              max_curves, thickness_end,
              sp_root_dir, sp_dir_pattern, sp_cube_filename, sp_out_filename)
    advanced_params = (outdir, frame_slice)

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import("md_analysis.main", "run_potential_analysis")

        kwargs: dict = {
            "output_dir": ctx[K.OUTDIR_RESOLVED],
            "_nest": False,
            "thickness_ang": ctx[K.THICKNESS],
            "center_mode": ctx[K.CENTER_MODE],
            "metal_elements": ctx[K.METAL_ELEMENTS],
            "layer_tol_ang": ctx[K.LAYER_TOL],
            "fermi_unit": ctx[K.FERMI_UNIT],
            "max_curves": ctx[K.MAX_CURVES],
            "thickness_end": ctx[K.THICKNESS_END],
            "frame_start": ctx[K.FRAME_START],
            "frame_end": ctx[K.FRAME_END],
            "frame_step": ctx[K.FRAME_STEP],
            "verbose": True,
            **_distributed_kwargs(ctx),
        }

        if _is_distributed(ctx):
            kwargs["md_out_path"] = None
            kwargs["cube_pattern"] = ""
            kwargs["xyz_path"] = None
        else:
            kwargs["cube_pattern"] = ctx[K.CUBE_PATTERN]
            kwargs["xyz_path"] = Path(ctx[K.XYZ])
            md_out = Path(ctx[K.MD_OUT])
            if not md_out.exists():
                print(f"  WARNING: {md_out} not found, "
                      "skipping Fermi/electrode/sensitivity analyses.")
                md_out = None
            kwargs["md_out_path"] = md_out

        results = analyze(**kwargs)
        print("\n Analysis complete. Outputs:")
        for name, path in results.items():
            print(f"   {name}: {path}")
