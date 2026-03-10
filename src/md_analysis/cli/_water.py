"""Water analysis sub-menu and handlers."""

from __future__ import annotations

from pathlib import Path

from ..config import KEY_Z_BIN_WIDTH_A
from ._prompt import (
    _get_effective_default,
    _handle_cmd_error,
    _prompt_bool,
    _prompt_cell_abc,
    _prompt_global_params,
    _prompt_str_required,
)

_MENU = """\

 ---------- Water Analysis ----------

 101) Water Mass Density Profile
 102) Water Orientation-Weighted Density Profile
 103) Adsorbed-Water Orientation Profile
 104) Adsorbed-Water Theta Distribution
 105) Full Water Three-Panel Analysis  (includes 101-104)

   0) Back / Exit

"""


def _collect_params() -> dict:
    """Prompt for water-specific required + optional parameters."""
    print()
    params: dict = {}
    params["xyz"] = _prompt_str_required("XYZ trajectory file (e.g. md-pos-1.xyz)")
    params["cell_abc"] = _prompt_cell_abc()
    params["dz_A"] = _get_effective_default(KEY_Z_BIN_WIDTH_A)

    if _prompt_bool("Modify advanced parameters?", default=False):
        params.update(_prompt_global_params())
    else:
        params["outdir"] = "analysis"
        params["frame_start"] = None
        params["frame_end"] = None
        params["frame_step"] = None

    return params


@_handle_cmd_error
def _cmd_101(params: dict) -> int:
    """Run water mass density analysis (menu 101)."""
    from ..water import water_mass_density_z_distribution_analysis

    outdir = Path(params["outdir"]) / "water"
    outdir.mkdir(parents=True, exist_ok=True)

    csv_path = water_mass_density_z_distribution_analysis(
        xyz_path=params["xyz"],
        cell_abc=params["cell_abc"],
        output_dir=outdir,
        dz_A=params["dz_A"],
        frame_start=params["frame_start"],
        frame_end=params["frame_end"],
        frame_step=params["frame_step"],
    )
    print(f"\n Analysis complete. Output:\n   density_csv: {csv_path}")
    return 0


@_handle_cmd_error
def _cmd_102(params: dict) -> int:
    """Run water orientation-weighted density analysis (menu 102)."""
    from ..water import water_orientation_weighted_density_z_distribution_analysis

    outdir = Path(params["outdir"]) / "water"
    outdir.mkdir(parents=True, exist_ok=True)

    csv_path = water_orientation_weighted_density_z_distribution_analysis(
        xyz_path=params["xyz"],
        cell_abc=params["cell_abc"],
        output_dir=outdir,
        dz_A=params["dz_A"],
        frame_start=params["frame_start"],
        frame_end=params["frame_end"],
        frame_step=params["frame_step"],
    )
    print(f"\n Analysis complete. Output:\n   orientation_csv: {csv_path}")
    return 0


@_handle_cmd_error
def _cmd_103(params: dict) -> int:
    """Run adsorbed-water orientation analysis (menu 103)."""
    from ..water import ad_water_orientation_analysis

    outdir = Path(params["outdir"]) / "water"
    outdir.mkdir(parents=True, exist_ok=True)

    profile_csv, range_txt = ad_water_orientation_analysis(
        xyz_path=params["xyz"],
        cell_abc=params["cell_abc"],
        output_dir=outdir,
        dz_A=params["dz_A"],
        frame_start=params["frame_start"],
        frame_end=params["frame_end"],
        frame_step=params["frame_step"],
    )
    print("\n Analysis complete. Outputs:")
    print(f"   adsorbed_profile_csv: {profile_csv}")
    print(f"   adsorbed_range_txt:   {range_txt}")
    return 0


@_handle_cmd_error
def _cmd_104(params: dict) -> int:
    """Run adsorbed-water theta distribution (menu 104)."""
    from ..water import compute_adsorbed_water_theta_distribution

    outdir = Path(params["outdir"]) / "water"
    outdir.mkdir(parents=True, exist_ok=True)

    _, _, csv_path = compute_adsorbed_water_theta_distribution(
        xyz_path=params["xyz"],
        cell_abc=params["cell_abc"],
        output_dir=outdir,
        dz_A=params["dz_A"],
        frame_start=params["frame_start"],
        frame_end=params["frame_end"],
        frame_step=params["frame_step"],
        verbose=True,
    )
    print(f"\n Analysis complete. Output:\n   theta_csv: {csv_path}")
    return 0


@_handle_cmd_error
def _cmd_105(params: dict) -> int:
    """Run full three-panel water analysis (menu 105)."""
    from ..main import run_water_analysis

    results = run_water_analysis(
        xyz_path=Path(params["xyz"]),
        cell_abc=params["cell_abc"],
        output_dir=Path(params["outdir"]),
        dz_A=params["dz_A"],
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
    "101": _cmd_101,
    "102": _cmd_102,
    "103": _cmd_103,
    "104": _cmd_104,
    "105": _cmd_105,
}


def water_menu() -> int:
    """Display the water sub-menu and dispatch."""
    print(_MENU)
    choice = input(" Input: ").strip()

    if choice == "0":
        print("\n Bye!")
        return 0

    handler = _DISPATCH.get(choice)
    if handler is None:
        print(f"\n Invalid choice: {choice!r}")
        return 1

    params = _collect_params()
    return handler(params)
