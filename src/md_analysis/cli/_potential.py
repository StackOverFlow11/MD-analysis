"""Potential analysis sub-menu and handlers."""

from __future__ import annotations

from pathlib import Path

from ..potential.config import DEFAULT_THICKNESS_ANG
from ._prompt import (
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

 ---------- Potential Analysis ----------

 201) Center Slab Potential (phi_center)
 202) Fermi Energy Time Series
 203) Electrode Potential (U vs SHE)    (includes 201+202)
 204) Phi(z) Plane-Averaged Profile
 205) Thickness Sensitivity Sweep
 206) Full Potential Analysis           (includes 201-205)

   0) Back / Exit

"""


def _collect_center_params(params: dict) -> None:
    """Prompt for center-slab-potential parameters."""
    params["cube_pattern"] = (
        _prompt_str("Cube file glob pattern", default="md-POTENTIAL-v_hartree-1_*.cube")
        or "md-POTENTIAL-v_hartree-1_*.cube"
    )
    params["xyz"] = _prompt_str("XYZ trajectory file", default="md-pos-1.xyz") or "md-pos-1.xyz"
    params["thickness"] = _prompt_float("Slab averaging thickness (A)", default=DEFAULT_THICKNESS_ANG)
    params["center_mode"] = _prompt_choice("Slab center mode", ["interface", "cell"], default="interface")
    params["metal_elements"] = _parse_metal_elements(
        _prompt_str("Metal elements (comma-separated, e.g. Cu,Ag)", default=None)
    )
    params["layer_tol"] = _prompt_float("Layer clustering tolerance (A)", default=0.6)


def _collect_fermi_params(params: dict) -> None:
    """Prompt for Fermi-energy parameters."""
    params["md_out"] = _prompt_str("CP2K md.out file", default="md.out") or "md.out"
    params["fermi_unit"] = _prompt_choice("Fermi energy unit in md.out", ["au", "ev"], default="au")


def _collect_phi_z_params(params: dict) -> None:
    """Prompt for phi(z) parameters."""
    params["cube_pattern"] = (
        _prompt_str("Cube file glob pattern", default="md-POTENTIAL-v_hartree-1_*.cube")
        or "md-POTENTIAL-v_hartree-1_*.cube"
    )
    params["max_curves"] = _prompt_int("Max curves on overlay (0=all)", default=0) or 0


def _collect_sensitivity_params(params: dict) -> None:
    """Prompt for thickness sensitivity parameters."""
    params["thickness_end"] = _prompt_float("Thickness sweep upper limit (A)", default=15.0)


def _collect_params_for(code: str) -> dict:
    """Prompt for parameters relevant to the chosen analysis code."""
    print()
    params: dict = {}

    if code in ("201", "203", "205", "206"):
        _collect_center_params(params)
    if code in ("202", "203", "205", "206"):
        _collect_fermi_params(params)
    if code in ("204",):
        _collect_phi_z_params(params)
    if code == "204" and "cube_pattern" not in params:
        params["cube_pattern"] = (
            _prompt_str("Cube file glob pattern", default="md-POTENTIAL-v_hartree-1_*.cube")
            or "md-POTENTIAL-v_hartree-1_*.cube"
        )
    if code in ("205", "206"):
        _collect_sensitivity_params(params)
    if code == "206":
        # Full analysis also needs phi_z params
        if "max_curves" not in params:
            params["max_curves"] = _prompt_int("Max curves on phi(z) overlay (0=all)", default=0) or 0

    if _prompt_bool("Modify advanced parameters?", default=False):
        params.update(_prompt_global_params())
    else:
        params["outdir"] = "analysis"
        params["frame_start"] = None
        params["frame_end"] = None
        params["frame_step"] = None

    return params


@_handle_cmd_error
def _cmd_201(params: dict) -> int:
    from ..potential import center_slab_potential_analysis

    outdir = Path(params["outdir"]) / "potential" / "center"
    outdir.mkdir(parents=True, exist_ok=True)

    csv_path = center_slab_potential_analysis(
        params["cube_pattern"],
        output_dir=outdir,
        thickness_ang=params["thickness"],
        center_mode=params["center_mode"],
        xyz_path=Path(params["xyz"]),
        metal_elements=params["metal_elements"],
        layer_tol_ang=params["layer_tol"],
        frame_start=params["frame_start"],
        frame_end=params["frame_end"],
        frame_step=params["frame_step"],
        verbose=True,
    )
    print(f"\n Analysis complete. Output:\n   center_csv: {csv_path}")
    return 0


@_handle_cmd_error
def _cmd_202(params: dict) -> int:
    from ..potential import fermi_energy_analysis

    outdir = Path(params["outdir"]) / "potential" / "fermi"
    outdir.mkdir(parents=True, exist_ok=True)

    csv_path = fermi_energy_analysis(
        Path(params["md_out"]),
        output_dir=outdir,
        fermi_unit=params["fermi_unit"],
        frame_start=params["frame_start"],
        frame_end=params["frame_end"],
        frame_step=params["frame_step"],
    )
    print(f"\n Analysis complete. Output:\n   fermi_csv: {csv_path}")
    return 0


@_handle_cmd_error
def _cmd_203(params: dict) -> int:
    from ..potential import electrode_potential_analysis

    outdir = Path(params["outdir"]) / "potential" / "electrode"
    outdir.mkdir(parents=True, exist_ok=True)

    csv_path = electrode_potential_analysis(
        params["cube_pattern"],
        Path(params["md_out"]),
        output_dir=outdir,
        thickness_ang=params["thickness"],
        center_mode=params["center_mode"],
        xyz_path=Path(params["xyz"]),
        metal_elements=params["metal_elements"],
        layer_tol_ang=params["layer_tol"],
        fermi_unit=params["fermi_unit"],
        frame_start=params["frame_start"],
        frame_end=params["frame_end"],
        frame_step=params["frame_step"],
        verbose=True,
    )
    print(f"\n Analysis complete. Output:\n   electrode_csv: {csv_path}")
    return 0


@_handle_cmd_error
def _cmd_204(params: dict) -> int:
    from ..potential import phi_z_planeavg_analysis

    outdir = Path(params["outdir"]) / "potential" / "phi_z"
    outdir.mkdir(parents=True, exist_ok=True)

    png_path = phi_z_planeavg_analysis(
        params["cube_pattern"],
        output_dir=outdir,
        max_curves=params["max_curves"],
        frame_start=params["frame_start"],
        frame_end=params["frame_end"],
        frame_step=params["frame_step"],
        verbose=True,
    )
    print(f"\n Analysis complete. Output:\n   phi_z_png: {png_path}")
    return 0


@_handle_cmd_error
def _cmd_205(params: dict) -> int:
    from ..potential import thickness_sensitivity_analysis

    outdir = Path(params["outdir"]) / "potential" / "thickness_sensitivity"
    outdir.mkdir(parents=True, exist_ok=True)

    csv_path = thickness_sensitivity_analysis(
        params["cube_pattern"],
        Path(params["md_out"]),
        output_dir=outdir,
        thickness_end=params["thickness_end"],
        center_mode=params["center_mode"],
        xyz_path=Path(params["xyz"]),
        metal_elements=params["metal_elements"],
        layer_tol_ang=params["layer_tol"],
        fermi_unit=params["fermi_unit"],
        frame_start=params["frame_start"],
        frame_end=params["frame_end"],
        frame_step=params["frame_step"],
        verbose=True,
    )
    print(f"\n Analysis complete. Output:\n   thickness_sensitivity_csv: {csv_path}")
    return 0


@_handle_cmd_error
def _cmd_206(params: dict) -> int:
    from ..main import run_potential_analysis

    md_out_path = Path(params["md_out"])
    if not md_out_path.exists():
        print(f"  WARNING: {md_out_path} not found, skipping Fermi/electrode/sensitivity analyses.")
        md_out_path = None

    results = run_potential_analysis(
        output_dir=Path(params["outdir"]),
        cube_pattern=params["cube_pattern"],
        md_out_path=md_out_path,
        xyz_path=Path(params["xyz"]),
        thickness_ang=params["thickness"],
        center_mode=params["center_mode"],
        metal_elements=params["metal_elements"],
        layer_tol_ang=params["layer_tol"],
        fermi_unit=params["fermi_unit"],
        max_curves=params.get("max_curves", 0),
        thickness_end=params.get("thickness_end", 15.0),
        frame_start=params["frame_start"],
        frame_end=params["frame_end"],
        frame_step=params["frame_step"],
        verbose=True,
    )
    print("\n Analysis complete. Outputs:")
    for name, path in results.items():
        print(f"   {name}: {path}")
    return 0


_DISPATCH = {
    "201": _cmd_201,
    "202": _cmd_202,
    "203": _cmd_203,
    "204": _cmd_204,
    "205": _cmd_205,
    "206": _cmd_206,
}


def potential_menu() -> int:
    """Display the potential sub-menu and dispatch."""
    print(_MENU)
    choice = input(" Input: ").strip()

    if choice == "0":
        print("\n Bye!")
        return 0

    handler = _DISPATCH.get(choice)
    if handler is None:
        print(f"\n Invalid choice: {choice!r}")
        return 1

    params = _collect_params_for(choice)
    return handler(params)
