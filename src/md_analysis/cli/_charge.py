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
    target_side,
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

    def __init__(self, code: str, label: str, *, method: str | None = None):
        super().__init__(code, label)
        if method is not None:
            self.params = (FixedParam(K.METHOD, method),)
            self.output_name = method          # e.g. "counterion" or "layer"
        else:
            self.params = (charge_method,)
            # output_name stays "" — dynamic method, handled in execute()

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import(
            "md_analysis.electrochemical.charge",
            "surface_charge_analysis",
        )
        from ..config import (
            KEY_POTENTIAL_PH,
            KEY_POTENTIAL_PHI_PZC,
            KEY_POTENTIAL_REFERENCE,
            KEY_POTENTIAL_TEMPERATURE_K,
            get_config,
        )

        # Static method (221/222): framework resolved via output_name
        # Dynamic method (223): append method to avoid overwriting
        if K.OUTDIR_RESOLVED in ctx:
            out = Path(ctx[K.OUTDIR_RESOLVED])
        else:
            out = (Path(ctx.get(K.OUTDIR, "analysis"))
                   / self.output_subdir)

        if not self.output_name:
            out = out / ctx[K.METHOD]

        out.mkdir(parents=True, exist_ok=True)

        csv_path = analyze(
            ctx[K.ROOT_DIR],
            metal_symbols=ctx[K.METAL_ELEMENTS],
            normal=ctx[K.NORMAL],
            method=ctx[K.METHOD],
            layer_tol_A=ctx[K.LAYER_TOL],
            n_surface_layers=ctx[K.N_SURFACE_LAYERS],
            dir_pattern=ctx[K.DIR_PATTERN],
            output_dir=out,
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
            verbose=True,
            potential_reference=get_config(KEY_POTENTIAL_REFERENCE, "SHE"),
            potential_pH=get_config(KEY_POTENTIAL_PH, 0.0),
            potential_temperature_K=get_config(KEY_POTENTIAL_TEMPERATURE_K, 298.15),
            potential_phi_pzc=get_config(KEY_POTENTIAL_PHI_PZC),
        )
        print(f"\n Analysis complete. Output:\n   charge_csv: {csv_path}")
        _print_ensemble_summary(csv_path)


class SingleSideChargeCmd(MenuCommand):
    """Single-side surface charge density with potential extrapolation."""

    params = (charge_method, target_side)
    advanced_params = (root_dir, dir_pattern, normal_axis, metal_elements,
                       layer_tol, n_surface_layers, outdir, frame_slice)

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import(
            "md_analysis.electrochemical.charge",
            "surface_charge_analysis",
        )
        from ..config import (
            KEY_POTENTIAL_PH,
            KEY_POTENTIAL_PHI_PZC,
            KEY_POTENTIAL_REFERENCE,
            KEY_POTENTIAL_TEMPERATURE_K,
            get_config,
        )

        method = ctx[K.METHOD]
        side = ctx[K.TARGET_SIDE]
        out = Path(ctx.get(K.OUTDIR, "analysis")) / self.output_subdir / f"{method}_{side}"
        if K.OUTDIR_RESOLVED in ctx:
            out = Path(ctx[K.OUTDIR_RESOLVED]) / f"{method}_{side}"
        out.mkdir(parents=True, exist_ok=True)

        csv_path = analyze(
            ctx[K.ROOT_DIR],
            metal_symbols=ctx[K.METAL_ELEMENTS],
            normal=ctx[K.NORMAL],
            method=method,
            layer_tol_A=ctx[K.LAYER_TOL],
            n_surface_layers=ctx[K.N_SURFACE_LAYERS],
            dir_pattern=ctx[K.DIR_PATTERN],
            output_dir=out,
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
            verbose=True,
            potential_reference=get_config(KEY_POTENTIAL_REFERENCE, "SHE"),
            potential_pH=get_config(KEY_POTENTIAL_PH, 0.0),
            potential_temperature_K=get_config(KEY_POTENTIAL_TEMPERATURE_K, 298.15),
            potential_phi_pzc=get_config(KEY_POTENTIAL_PHI_PZC),
            target_side=side,
        )
        print(f"\n Analysis complete ({side} side). Output:\n   charge_csv: {csv_path}")

        # Ensemble summary for single side
        if csv_path.exists():
            import numpy as np
            with csv_path.open(encoding="utf-8") as f:
                reader = _csv.DictReader(f)
                rows = list(reader)
            if rows:
                sigma = np.array([float(r["sigma_uC_cm2"]) for r in rows])
                print(f"\n Ensemble average ({len(rows)} frames, {side} side):")
                print(f"   sigma: {sigma.mean():8.4f} +/- {sigma.std():.4f} uC/cm^2")
                phi_col = [c for c in rows[0] if c.startswith("phi_cumavg")]
                if phi_col:
                    phi = np.array([float(r[phi_col[0]]) for r in rows])
                    ref = phi_col[0].split("_vs_")[-1]
                    print(f"   phi:   {phi[-1]:8.4f} V vs {ref} (cum. avg)")


class TrackedChargeCmd(MenuCommand):
    """Track Bader net charges for specified XYZ atoms."""

    output_name = "tracked"
    params = (atom_indices_xyz,)
    advanced_params = (root_dir, dir_pattern, outdir, frame_slice)

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import(
            "md_analysis.electrochemical.charge",
            "tracked_atom_charge_analysis",
        )
        csv_path = analyze(
            ctx[K.ROOT_DIR],
            atom_indices_xyz=ctx[K.ATOM_INDICES_XYZ],
            dir_pattern=ctx[K.DIR_PATTERN],
            output_dir=ctx[K.OUTDIR_RESOLVED],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
            verbose=True,
        )
        print(f"\n Analysis complete. Output:\n   tracked_csv: {csv_path}")


class CounterionChargeCmd(MenuCommand):
    """Detect counterions per-frame and track their Bader charges."""

    output_name = "counterion_tracking"
    advanced_params = (root_dir, dir_pattern, normal_axis, metal_elements,
                       layer_tol, outdir, frame_slice)

    def execute(self, ctx: dict) -> None:
        analyze = lazy_import(
            "md_analysis.electrochemical.charge",
            "counterion_charge_analysis",
        )
        csv_path = analyze(
            ctx[K.ROOT_DIR],
            metal_symbols=ctx[K.METAL_ELEMENTS],
            normal=ctx[K.NORMAL],
            layer_tol_A=ctx[K.LAYER_TOL],
            dir_pattern=ctx[K.DIR_PATTERN],
            output_dir=ctx[K.OUTDIR_RESOLVED],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
            verbose=True,
        )
        print(f"\n Analysis complete. Output:\n   counterion_csv: {csv_path}")
