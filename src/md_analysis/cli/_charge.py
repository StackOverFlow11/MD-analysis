"""Charge analysis sub-menu and handlers."""

from __future__ import annotations

from pathlib import Path

from ..config import KEY_LAYER_TOL_A
from ..utils.config import CHARGE_METHOD_COUNTERION, CHARGE_METHOD_LAYER
from ..charge.config import DEFAULT_N_SURFACE_LAYERS
from ._prompt import (
    _get_effective_default,
    _handle_cmd_error,
    _parse_metal_elements,
    _prompt_bool,
    _prompt_choice,
    _prompt_float,
    _prompt_global_params,
    _prompt_int,
    _prompt_str,
)

_MENU = """\

 ---------- Charge Analysis ----------

 301) Surface Charge (Counterion Method)
 302) Surface Charge (Layer Method)
 303) Full Charge Analysis with Plots

   0) Back / Exit

"""


def _collect_params(*, method: str | None = None) -> dict:
    """Prompt for charge-specific parameters."""
    print()
    params: dict = {}

    if method is not None:
        params["method"] = method
    else:
        params["method"] = _prompt_choice(
            "Charge method",
            [CHARGE_METHOD_COUNTERION, CHARGE_METHOD_LAYER],
            default=CHARGE_METHOD_COUNTERION,
        )

    if _prompt_bool("Modify advanced parameters?", default=False):
        params["root_dir"] = _prompt_str("Root directory", default=".") or "."
        params["dir_pattern"] = (
            _prompt_str("Frame subdirectory pattern", default="bader_t*_i*")
            or "bader_t*_i*"
        )
        params["normal"] = _prompt_choice(
            "Surface normal axis", ["a", "b", "c"], default="c",
        )
        params["metal_elements"] = _parse_metal_elements(
            _prompt_str("Metal elements (comma-separated)", default=None)
        )
        params["layer_tol"] = _prompt_float(
            "Layer clustering tolerance (A)",
            default=_get_effective_default(KEY_LAYER_TOL_A),
        )
        if params["method"] == CHARGE_METHOD_LAYER:
            params["n_surface_layers"] = _prompt_int(
                "Number of surface layers per interface",
                default=DEFAULT_N_SURFACE_LAYERS,
            )
        params.update(_prompt_global_params())
    else:
        params["root_dir"] = "."
        params["dir_pattern"] = "bader_t*_i*"
        params["normal"] = "c"
        params["metal_elements"] = None
        params["layer_tol"] = _get_effective_default(KEY_LAYER_TOL_A)
        params["n_surface_layers"] = DEFAULT_N_SURFACE_LAYERS
        params["outdir"] = "analysis"
        params["frame_start"] = None
        params["frame_end"] = None
        params["frame_step"] = None

    return params


def _print_ensemble_summary(csv_path: Path) -> None:
    """Print ensemble average summary from the charge CSV."""
    import csv as _csv

    if not csv_path.exists():
        return

    with csv_path.open(encoding="utf-8") as f:
        reader = _csv.DictReader(f)
        rows = list(reader)

    if not rows:
        return

    import numpy as _np

    aligned = _np.array([float(r["sigma_aligned_uC_cm2"]) for r in rows])
    opposed = _np.array([float(r["sigma_opposed_uC_cm2"]) for r in rows])
    print(f"\n Ensemble average ({len(rows)} frames):")
    print(f"   sigma_aligned: {aligned.mean():8.4f} +/- {aligned.std():.4f} uC/cm^2")
    print(f"   sigma_opposed: {opposed.mean():8.4f} +/- {opposed.std():.4f} uC/cm^2")


@_handle_cmd_error
def _run_charge(params: dict) -> int:
    """Shared handler for 301/302/303."""
    from ..main import run_charge_analysis

    results = run_charge_analysis(
        output_dir=Path(params["outdir"]),
        root_dir=params["root_dir"],
        metal_symbols=params["metal_elements"],
        normal=params["normal"],
        method=params["method"],
        layer_tol_A=params["layer_tol"],
        n_surface_layers=params.get("n_surface_layers", DEFAULT_N_SURFACE_LAYERS),
        dir_pattern=params["dir_pattern"],
        frame_start=params["frame_start"],
        frame_end=params["frame_end"],
        frame_step=params["frame_step"],
        verbose=True,
    )

    print("\n Analysis complete. Outputs:")
    for name, path in results.items():
        print(f"   {name}: {path}")

    _print_ensemble_summary(results["charge_csv"])
    return 0


def charge_menu() -> int:
    """Display the charge sub-menu and dispatch."""
    print(_MENU)
    choice = input(" Input: ").strip()

    if choice == "0":
        print("\n Bye!")
        return 0

    if choice == "301":
        params = _collect_params(method=CHARGE_METHOD_COUNTERION)
        return _run_charge(params)
    elif choice == "302":
        params = _collect_params(method=CHARGE_METHOD_LAYER)
        return _run_charge(params)
    elif choice == "303":
        params = _collect_params()
        return _run_charge(params)
    else:
        print(f"\n Invalid choice: {choice!r}")
        return 1
