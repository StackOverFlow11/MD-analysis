"""Calibration command classes (231-233).

Phase 6.4: routed through ``workflows.calibration.run_*`` facades.
The workflow keeps ``calibration_json_path`` as a required argument
(agent-safe semantics); the CLI resolves the global default
``~/.config/md_analysis/calibration.json`` here so the pre-migration
UX is preserved.
"""

from __future__ import annotations

from pathlib import Path

from ._framework import MenuCommand, lazy_import
from ._params import (
    K,
    calibration_csv_path,
    calibration_json,
    fitting_method,
    outdir,
    ph_value,
    phi_pzc,
    poly_degree,
    potential_reference,
    sigma_value,
    temperature_k,
)
from ._prompt import _read, prompt_float


class CalibrateFromCSVCmd(MenuCommand):
    """Calibrate σ→φ mapping from a CSV file."""

    output_name = "fit"
    params = (calibration_csv_path, fitting_method)
    advanced_params = (poly_degree, outdir)

    def execute(self, ctx: dict) -> None:
        run_fit = lazy_import(
            "md_analysis.workflows.calibration", "run_calibration_fit",
        )
        default_json = lazy_import(
            "md_analysis.electrochemical.calibration.config",
            "DEFAULT_CALIBRATION_FILE",
        )
        print("\n  CSV format: column 1 = potential φ (V vs SHE)")
        print("              column 2 = surface charge density σ (μC/cm²)")
        print("  Header row is optional (auto-detected).\n")

        result = run_fit(
            csv_path=ctx[K.CALIBRATION_CSV],
            method=ctx[K.FITTING_METHOD],
            poly_degree=ctx.get(K.POLY_DEGREE, 2),
            output_dir=ctx.get(K.OUTDIR_RESOLVED),
            calibration_json_path=default_json,
        )
        print(f"\n Calibration complete.")
        print(f"   JSON saved: {result.artifacts['calibration_json']}")


class CalibrateManualCmd(MenuCommand):
    """Calibrate σ→φ mapping from manually entered data points."""

    output_name = "fit"
    params = (fitting_method,)
    advanced_params = (poly_degree, outdir)

    def execute(self, ctx: dict) -> None:
        run_fit = lazy_import(
            "md_analysis.workflows.calibration", "run_calibration_fit",
        )
        default_json = lazy_import(
            "md_analysis.electrochemical.calibration.config",
            "DEFAULT_CALIBRATION_FILE",
        )
        print("\n  Enter calibration data points.")
        print("  φ in V vs SHE, σ in μC/cm².")
        print("  Type 'done' when finished (minimum 2 points).\n")

        points: list[tuple[float, float]] = []
        i = 1
        while True:
            raw = _read(f"  Point {i} — φ (V vs SHE) [or 'done']: ").strip()
            if raw.lower() == "done":
                if len(points) < 2:
                    print("  Need at least 2 points.")
                    continue
                break
            try:
                phi = float(raw)
            except ValueError:
                print(f"  Invalid number: {raw!r}")
                continue
            sigma = prompt_float(f"  Point {i} — σ (μC/cm²)", default=0.0)
            points.append((phi, sigma))
            print(f"    → recorded ({phi}, {sigma})")
            i += 1

        result = run_fit(
            data_points=points,
            method=ctx[K.FITTING_METHOD],
            poly_degree=ctx.get(K.POLY_DEGREE, 2),
            output_dir=ctx.get(K.OUTDIR_RESOLVED),
            calibration_json_path=default_json,
        )
        print(f"\n Calibration complete ({len(points)} points).")
        print(f"   JSON saved: {result.artifacts['calibration_json']}")


class PredictPotentialCmd(MenuCommand):
    """Predict electrode potential from surface charge density."""

    output_name = "predict"
    params = (sigma_value,)
    advanced_params = (calibration_json, potential_reference,
                       temperature_k, ph_value, phi_pzc)

    def execute(self, ctx: dict) -> None:
        run_predict = lazy_import(
            "md_analysis.workflows.calibration", "run_calibration_predict",
        )
        default_json = lazy_import(
            "md_analysis.electrochemical.calibration.config",
            "DEFAULT_CALIBRATION_FILE",
        )
        json_path = ctx.get(K.CALIBRATION_JSON) or default_json
        result = run_predict(
            ctx[K.SIGMA_VALUE],
            calibration_json_path=Path(json_path),
            target_reference=ctx.get(K.POTENTIAL_REFERENCE, "SHE"),
            temperature_K=ctx.get(K.TEMPERATURE_K, 298.15),
            pH=ctx.get(K.PH, 0.0),
            phi_pzc=ctx.get(K.PHI_PZC),
        )
        phi_val = float(result.extra.potential_V[0])
        ref = result.metadata["target_reference"]
        print(f"\n  σ = {ctx[K.SIGMA_VALUE]:.4f} μC/cm²")
        print(f"  φ = {phi_val:.6f} V vs {ref}")
