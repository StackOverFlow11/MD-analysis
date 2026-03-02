"""Command-line interface for md-analysis.

Usage::

    md-analysis water --xyz md-pos-1.xyz --md-inp md.inp
    md-analysis potential --cube-pattern "md-POTENTIAL-v_hartree-1_*.cube"
    md-analysis all --xyz md-pos-1.xyz --md-inp md.inp
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _add_water_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--xyz", type=str, required=True, help="XYZ trajectory file.")
    parser.add_argument("--md-inp", type=str, required=True, help="CP2K md.inp file.")


def _add_potential_args(parser: argparse.ArgumentParser, *, skip_xyz: bool = False) -> None:
    parser.add_argument(
        "--cube-pattern", type=str, default="md-POTENTIAL-v_hartree-1_*.cube",
        help="Glob pattern for cube files (default: md-POTENTIAL-v_hartree-1_*.cube).",
    )
    parser.add_argument("--md-out", type=str, default="md.out", help="CP2K md.out file (default: md.out).")
    if not skip_xyz:
        parser.add_argument("--xyz", type=str, default="md-pos-1.xyz", help="XYZ trajectory for interface detection (default: md-pos-1.xyz).")
    parser.add_argument("--thickness", type=float, default=7.0, help="Slab averaging thickness in Å (default: 7).")
    parser.add_argument("--center-mode", choices=["interface", "cell"], default="interface", help="Slab center mode (default: interface).")
    parser.add_argument("--fermi-unit", choices=["au", "ev"], default="au", help="Fermi energy unit in md.out (default: au).")
    parser.add_argument("--metal-elements", type=str, default=None, help="Comma-separated metal element symbols (e.g., Cu,Ag).")
    parser.add_argument("--layer-tol", type=float, default=0.6, help="Metal layer clustering tolerance in Å (default: 0.6).")
    parser.add_argument("--no-compute-u", action="store_true", help="Skip U vs SHE computation.")
    parser.add_argument("--no-phi-z", action="store_true", help="Skip φ(z) visualization.")
    parser.add_argument("--z-mavg-window", type=float, default=7.0, help="Moving average window for φ(z) in Å (default: 7).")
    parser.add_argument("--max-curves", type=int, default=0, help="Max curves on φ(z) overlay (0=all).")


def _parse_metal_elements(s: str | None) -> set[str] | None:
    if s is None:
        return None
    elements = {e.strip() for e in s.split(",") if e.strip()}
    return elements or None


def _cmd_water(args: argparse.Namespace) -> int:
    from .main import run_water_analysis

    results = run_water_analysis(
        xyz_path=Path(args.xyz),
        md_inp_path=Path(args.md_inp),
        output_dir=Path(args.outdir),
        frame_start=args.frame_start,
        frame_end=args.frame_end,
        frame_step=args.frame_step,
        verbose=True,
    )

    print("Water analysis complete. Outputs:")
    for name, path in results.items():
        print(f"  {name}: {path}")
    return 0


def _cmd_potential(args: argparse.Namespace) -> int:
    from .main import run_potential_analysis

    md_out_path = Path(args.md_out) if hasattr(args, "md_out") else None
    if md_out_path is not None and not md_out_path.exists():
        md_out_path = None

    results = run_potential_analysis(
        output_dir=Path(args.outdir),
        cube_pattern=args.cube_pattern,
        md_out_path=md_out_path,
        xyz_path=Path(args.xyz) if args.xyz else None,
        thickness_ang=args.thickness,
        center_mode=args.center_mode,
        metal_elements=_parse_metal_elements(args.metal_elements),
        layer_tol_ang=args.layer_tol,
        fermi_unit=args.fermi_unit,
        compute_u=not args.no_compute_u,
        compute_phi_z=not args.no_phi_z,
        z_mavg_window_ang=args.z_mavg_window,
        max_curves=args.max_curves,
        frame_start=args.frame_start,
        frame_end=args.frame_end,
        frame_step=args.frame_step,
        verbose=True,
    )

    print("Potential analysis complete. Outputs:")
    for name, path in results.items():
        print(f"  {name}: {path}")
    return 0


def _cmd_all(args: argparse.Namespace) -> int:
    confirm = input("Will run all analyses. Enter 'yes' to confirm: ")
    if confirm.strip().lower() != "yes":
        print("Aborted.")
        return 1

    from .main import run_all

    md_out_path = Path(args.md_out) if hasattr(args, "md_out") else None
    if md_out_path is not None and not md_out_path.exists():
        md_out_path = None

    results = run_all(
        xyz_path=Path(args.xyz),
        md_inp_path=Path(args.md_inp),
        output_dir=Path(args.outdir),
        cube_pattern=args.cube_pattern,
        md_out_path=md_out_path,
        frame_start=args.frame_start,
        frame_end=args.frame_end,
        frame_step=args.frame_step,
        verbose=True,
        thickness_ang=args.thickness,
        center_mode=args.center_mode,
        metal_elements=_parse_metal_elements(args.metal_elements),
        layer_tol_ang=args.layer_tol,
        fermi_unit=args.fermi_unit,
        compute_u=not args.no_compute_u,
        compute_phi_z=not args.no_phi_z,
        z_mavg_window_ang=args.z_mavg_window,
        max_curves=args.max_curves,
    )

    print("All analyses complete. Outputs:")
    for name, path in results.items():
        print(f"  {name}: {path}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="md-analysis",
        description="MD analysis utilities for periodic metal-water interfaces from CP2K simulations.",
    )
    parser.add_argument(
        "--outdir", type=str, default="analysis",
        help="Output root directory (default: ./analysis).",
    )
    parser.add_argument(
        "--frame-start", type=int, default=None,
        help="Start frame index, 0-based (default: 0, i.e. first frame).",
    )
    parser.add_argument(
        "--frame-end", type=int, default=None,
        help="End frame index, exclusive (default: no limit, i.e. all frames).",
    )
    parser.add_argument(
        "--frame-step", type=int, default=None,
        help="Frame step (default: 1).",
    )
    subparsers = parser.add_subparsers(dest="command", help="Analysis sub-commands")

    # water subcommand
    water_parser = subparsers.add_parser("water", help="Water analysis (density + orientation + adsorbed + three-panel plot).")
    _add_water_args(water_parser)

    # potential subcommand
    potential_parser = subparsers.add_parser("potential", help="Potential analysis (center + fermi + electrode + phi_z).")
    _add_potential_args(potential_parser)

    # all subcommand (--xyz is shared between water and potential, so add it via water only)
    all_parser = subparsers.add_parser("all", help="Run all analyses (water + potential). Requires confirmation.")
    _add_water_args(all_parser)
    _add_potential_args(all_parser, skip_xyz=True)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return 0

    if args.command == "water":
        return _cmd_water(args)
    elif args.command == "potential":
        return _cmd_potential(args)
    elif args.command == "all":
        return _cmd_all(args)
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
