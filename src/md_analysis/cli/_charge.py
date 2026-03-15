"""Charge analysis command classes (221-225)."""

from __future__ import annotations

import csv as _csv
from pathlib import Path

from ._framework import MenuCommand, lazy_import
from ._params import (
    K,
    FixedParam,
    atom_indices_xyz,
    charge_method,
    dir_pattern,
    frame_slice,
    layer_tol,
    metal_elements,
    n_surface_layers,
    normal_axis,
    outdir,
    root_dir,
)


def _print_ensemble_summary(csv_path: Path) -> None:
    """Print ensemble average summary from the charge CSV."""
    if not csv_path.exists():
        return

    with csv_path.open(encoding="utf-8") as f:
        reader = _csv.DictReader(f)
        rows = list(reader)

    if not rows:
        return

    import numpy as np

    aligned = np.array([float(r["sigma_aligned_uC_cm2"]) for r in rows])
    opposed = np.array([float(r["sigma_opposed_uC_cm2"]) for r in rows])
    print(f"\n Ensemble average ({len(rows)} frames):")
    print(f"   sigma_aligned: {aligned.mean():8.4f} +/- {aligned.std():.4f} uC/cm^2")
    print(f"   sigma_opposed: {opposed.mean():8.4f} +/- {opposed.std():.4f} uC/cm^2")


class SurfaceChargeCmd(MenuCommand):
    advanced_params = (root_dir, dir_pattern, normal_axis, metal_elements,
                       layer_tol, n_surface_layers, outdir, frame_slice)
    output_subdir = ""

    def __init__(self, code: str, label: str, *, method: str | None = None):
        super().__init__(code, label)
        if method is not None:
            self.params = (FixedParam(K.METHOD, method),)
        else:
            self.params = (charge_method,)

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import("md_analysis.main", "run_charge_analysis")
        results = analyze(
            output_dir=Path(ctx[K.OUTDIR]),
            root_dir=ctx[K.ROOT_DIR],
            metal_symbols=ctx[K.METAL_ELEMENTS],
            normal=ctx[K.NORMAL],
            method=ctx[K.METHOD],
            layer_tol_A=ctx[K.LAYER_TOL],
            n_surface_layers=ctx[K.N_SURFACE_LAYERS],
            dir_pattern=ctx[K.DIR_PATTERN],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
            verbose=True,
        )
        print("\n Analysis complete. Outputs:")
        for name, path in results.items():
            print(f"   {name}: {path}")

        _print_ensemble_summary(results["charge_csv"])


class TrackedChargeCmd(MenuCommand):
    """Track Bader net charges for specified XYZ atoms."""

    params = (atom_indices_xyz,)
    advanced_params = (root_dir, dir_pattern, outdir, frame_slice)
    output_subdir = ""

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import("md_analysis.main", "run_tracked_charge_analysis")
        results = analyze(
            output_dir=Path(ctx[K.OUTDIR]),
            root_dir=ctx[K.ROOT_DIR],
            atom_indices_xyz=ctx[K.ATOM_INDICES_XYZ],
            dir_pattern=ctx[K.DIR_PATTERN],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
            verbose=True,
        )
        print("\n Analysis complete. Outputs:")
        for name, path in results.items():
            print(f"   {name}: {path}")


class CounterionChargeCmd(MenuCommand):
    """Detect counterions per-frame and track their Bader charges."""

    advanced_params = (root_dir, dir_pattern, normal_axis, metal_elements,
                       layer_tol, outdir, frame_slice)
    output_subdir = ""

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import("md_analysis.main", "run_counterion_charge_analysis")
        results = analyze(
            output_dir=Path(ctx[K.OUTDIR]),
            root_dir=ctx[K.ROOT_DIR],
            metal_symbols=ctx[K.METAL_ELEMENTS],
            normal=ctx[K.NORMAL],
            layer_tol_A=ctx[K.LAYER_TOL],
            dir_pattern=ctx[K.DIR_PATTERN],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
            verbose=True,
        )
        print("\n Analysis complete. Outputs:")
        for name, path in results.items():
            print(f"   {name}: {path}")
