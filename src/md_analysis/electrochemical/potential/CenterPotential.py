"""Multi-frame center-slab potential, Fermi energy, and electrode potential analysis.

Public API
----------
- ``center_slab_potential_analysis``
- ``fermi_energy_analysis``
- ``electrode_potential_analysis``
- ``thickness_sensitivity_analysis``
"""

from __future__ import annotations

import csv
import logging
import math
import re
from pathlib import Path
from typing import Iterable, Optional

import numpy as np

logger = logging.getLogger(__name__)

from ...utils.config import (
    BOHR_TO_ANG,
    DEFAULT_LAYER_TOL_A,
    DP_A_H3O_W_EV,
    DELTA_E_ZP_EV,
    HA_TO_EV,
    MU_HPLUS_G0_EV,
    TRANSITION_METAL_SYMBOLS,
)
from ...utils._io_helpers import _cumulative_average, _write_csv
from ...utils.CubeParser import (
    _float,
    discover_cube_files,
    extract_step_from_cube_filename,
    read_cube_header_and_values,
    slab_average_potential_ev,
)
from ...utils.StructureParser.ClusterUtils import gap_midpoint_periodic
from ...utils.StructureParser.LayerParser import detect_interface_layers

try:
    from ase import Atoms
except ImportError:  # pragma: no cover
    Atoms = object  # type: ignore[misc]
from ._frame_source import (
    PotentialFrame,
    discover_continuous_frames,
    discover_distributed_frames,
)
from .config import (
    DEFAULT_CENTER_POTENTIAL_CSV_NAME,
    DEFAULT_CENTER_POTENTIAL_PNG_NAME,
    DEFAULT_ELECTRODE_POTENTIAL_CSV_NAME,
    DEFAULT_ELECTRODE_POTENTIAL_PNG_NAME,
    DEFAULT_FERMI_ENERGY_CSV_NAME,
    DEFAULT_FERMI_ENERGY_PNG_NAME,
    DEFAULT_SLAB_CENTER_CSV_NAME,
    DEFAULT_THICKNESS_SENSITIVITY_CSV_NAME,
    DEFAULT_THICKNESS_SENSITIVITY_PNG_NAME,
    DEFAULT_THICKNESS_ANG,
)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _extract_interface_geometry(detection, axis_length_ang: float) -> dict:
    """Extract interface geometry from a ``SurfaceDetectionResult``."""
    aligned = detection.interface_normal_aligned()
    opposed = detection.interface_normal_opposed()
    z_lower_ang = aligned.center_frac * axis_length_ang
    z_upper_ang = opposed.center_frac * axis_length_ang
    z_mid_ang = gap_midpoint_periodic(
        aligned.center_frac, opposed.center_frac, 1.0,
    ) * axis_length_ang
    water_gap_ang = (
        (opposed.center_frac - aligned.center_frac) % 1.0
    ) * axis_length_ang
    return {
        "z_lower_ang": z_lower_ang,
        "z_upper_ang": z_upper_ang,
        "z_mid_ang": z_mid_ang,
        "water_gap_ang": water_gap_ang,
        "n_layers": len(detection.metal_layers_sorted),
    }


def _plot_series_with_cumavg(
    png_path: Path,
    x: np.ndarray,
    y: np.ndarray,
    y_cum: np.ndarray,
    xlabel: str,
    ylabel: str,
    title: str,
) -> None:
    """Plot instantaneous values with cumulative average overlay."""
    png_path.parent.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 4.8), dpi=160)
    ax.plot(x, y, lw=1.0, alpha=0.65, label="instantaneous")
    ax.plot(x, y_cum, lw=2.0, label="cumulative average")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(png_path)
    plt.close(fig)


def _parse_csv_symbols(s: Optional[str]) -> Optional[set[str]]:
    """Parse comma-separated element symbols into a set, or None."""
    if s is None:
        return None
    raw = [p.strip() for p in s.split(",")]
    out = {p for p in raw if p}
    return out or None


def _resolve_frames(
    *,
    input_mode: str,
    # Mode A params
    cube_pattern: str = "",
    workdir: Path | None = None,
    md_out_path: Path | None = None,
    xyz_path: Path | None = None,
    # Mode B params
    sp_root_dir: Path | str | None = None,
    sp_dir_pattern: str = "potential_t*_i*",
    sp_cube_filename: str = "sp_potential-v_hartree-1_0.cube",
    sp_out_filename: str = "sp.out",
    # Shared params
    center_mode: str = "interface",
    metal_elements: set[str] | None = None,
    layer_tol_ang: float = DEFAULT_LAYER_TOL_A,
    fermi_unit: str = "au",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> list[PotentialFrame]:
    """Route to the correct frame discovery function based on *input_mode*."""
    if input_mode == "distributed":
        if sp_root_dir is None:
            raise ValueError("sp_root_dir is required for distributed mode")
        return discover_distributed_frames(
            sp_root_dir,
            dir_pattern=sp_dir_pattern,
            cube_filename=sp_cube_filename,
            sp_out_filename=sp_out_filename,
            center_mode=center_mode,
            metal_elements=metal_elements,
            layer_tol_ang=layer_tol_ang,
            frame_start=frame_start,
            frame_end=frame_end,
            frame_step=frame_step,
            verbose=verbose,
        )
    # Default: continuous
    return discover_continuous_frames(
        cube_pattern,
        workdir=workdir,
        md_out_path=md_out_path,
        xyz_path=xyz_path,
        center_mode=center_mode,
        metal_elements=metal_elements,
        fermi_unit=fermi_unit,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
    )


def _analyze_center_from_frames(
    frames: list[PotentialFrame],
    *,
    thickness_ang: float,
    center_mode: str,
    metal_elements: set[str] | None,
    layer_tol_ang: float,
    output_dir: Path,
    verbose: bool = False,
) -> Path:
    """Core center-slab potential analysis on pre-loaded frames.

    Writes center_potential.csv + .png + slab_center_and_interfaces.csv.
    Returns the CSV path.
    """
    outdir = output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    use_interface_center = center_mode == "interface"
    metal_elements_used = metal_elements

    cube_rows: list[dict] = []
    cube_steps: list[int] = []
    cube_vals_ev: list[float] = []
    center_rows: list[dict] = []

    frame_iter: Iterable[PotentialFrame] = frames
    if verbose:
        from tqdm import tqdm
        frame_iter = tqdm(frames, desc="Center slab potential", unit="cube", ascii=" =")

    for frame in frame_iter:
        header = frame.header
        values = frame.values
        step = frame.step

        dz_bohr = float(np.linalg.norm(header.vz_bohr))
        lz_ang = header.nz * dz_bohr * BOHR_TO_ANG
        origin_z_ang = float(header.origin_bohr[2]) * BOHR_TO_ANG

        z_center_ang: Optional[float] = None
        iface: Optional[dict] = None
        center_source = "cell"

        if use_interface_center and frame.atoms is not None:
            if metal_elements_used is None:
                metal_elements_used = (
                    set(frame.atoms.get_chemical_symbols())
                    & set(TRANSITION_METAL_SYMBOLS)
                )
            if metal_elements_used:
                detection = detect_interface_layers(
                    frame.atoms, metal_symbols=metal_elements_used,
                    normal="c", layer_tol_A=layer_tol_ang,
                )
                iface = _extract_interface_geometry(detection, lz_ang)
                z_center_ang = origin_z_ang + float(iface["z_mid_ang"])
                center_source = "interface"

        phi_ev, info = slab_average_potential_ev(
            header, values, thickness_ang=thickness_ang, z_center_ang=z_center_ang,
        )
        cube_steps.append(int(step))
        cube_vals_ev.append(float(phi_ev))

        center_rows.append({
            "step": int(step),
            "center_source": center_source,
            "z_center_ang": None if z_center_ang is None else float(z_center_ang),
            "z_iface_lower_ang": None if iface is None else float(iface["z_lower_ang"] + origin_z_ang),
            "z_iface_upper_ang": None if iface is None else float(iface["z_upper_ang"] + origin_z_ang),
            "z_iface_mid_ang": None if iface is None else float(iface["z_mid_ang"] + origin_z_ang),
            "water_gap_ang": None if iface is None else float(iface["water_gap_ang"]),
            "n_metal_layers": None if iface is None else int(iface["n_layers"]),
            "Lz_ang": float(lz_ang),
        })

    # Sort by step
    order = np.argsort(np.array(cube_steps, dtype=int))
    cube_steps_arr = np.array(cube_steps, dtype=int)[order]
    cube_vals_arr = np.array(cube_vals_ev, dtype=float)[order]
    cube_cum = _cumulative_average(cube_vals_arr)

    rows_out: list[dict] = []
    for s, v, c in zip(cube_steps_arr, cube_vals_arr, cube_cum, strict=True):
        rows_out.append({"step": int(s), "phi_center_ev": float(v), "phi_center_cumavg_ev": float(c)})

    csv_path = outdir / DEFAULT_CENTER_POTENTIAL_CSV_NAME
    _write_csv(csv_path, rows_out, ["step", "phi_center_ev", "phi_center_cumavg_ev"])

    if center_rows:
        center_rows_sorted = sorted(center_rows, key=lambda r: int(r["step"]))
        _write_csv(
            outdir / DEFAULT_SLAB_CENTER_CSV_NAME,
            center_rows_sorted,
            ["step", "center_source", "z_center_ang", "z_iface_lower_ang",
             "z_iface_upper_ang", "z_iface_mid_ang", "water_gap_ang",
             "n_metal_layers", "Lz_ang"],
        )

    _plot_series_with_cumavg(
        outdir / DEFAULT_CENTER_POTENTIAL_PNG_NAME,
        x=cube_steps_arr.astype(float),
        y=cube_vals_arr,
        y_cum=cube_cum,
        xlabel="MD step",
        ylabel="Center slab potential (eV)",
        title=f"Center water slab potential (thickness={thickness_ang:g} A)",
    )

    return csv_path


def _analyze_fermi_from_frames(
    frames: list[PotentialFrame],
    *,
    fermi_unit: str = "au",
    output_dir: Path,
) -> Path:
    """Core Fermi energy analysis on pre-loaded frames.

    Writes fermi_energy.csv + .png. Returns the CSV path.
    """
    outdir = output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Filter frames with fermi data
    fermi_frames = [f for f in frames if f.fermi_raw is not None]
    if not fermi_frames:
        raise RuntimeError("No frames with Fermi energy data found.")
    logger.info("Fermi energy: %d records", len(fermi_frames))

    f_steps = np.array([f.step for f in fermi_frames], dtype=int)
    f_time = np.array(
        [f.time_fs if f.time_fs is not None else math.nan for f in fermi_frames],
        dtype=float,
    )
    f_raw = np.array([f.fermi_raw for f in fermi_frames], dtype=float)
    f_vals = f_raw * HA_TO_EV if fermi_unit == "au" else f_raw
    f_cum = _cumulative_average(f_vals)

    f_rows: list[dict] = []
    for s, t, raw, v, c in zip(f_steps, f_time, f_raw, f_vals, f_cum, strict=True):
        f_rows.append({
            "step": int(s),
            "time_fs": None if math.isnan(float(t)) else float(t),
            "fermi_raw": float(raw),
            "fermi_ev": float(v),
            "fermi_cumavg_ev": float(c),
        })

    csv_path = outdir / DEFAULT_FERMI_ENERGY_CSV_NAME
    _write_csv(csv_path, f_rows, ["step", "time_fs", "fermi_raw", "fermi_ev", "fermi_cumavg_ev"])

    _plot_series_with_cumavg(
        outdir / DEFAULT_FERMI_ENERGY_PNG_NAME,
        x=f_steps.astype(float),
        y=f_vals,
        y_cum=f_cum,
        xlabel="MD step",
        ylabel="Fermi energy (eV)",
        title="Fermi energy from md.out",
    )

    return csv_path


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def center_slab_potential_analysis(
    cube_pattern: str = "",
    *,
    output_dir: Path | None = None,
    thickness_ang: float = DEFAULT_THICKNESS_ANG,
    center_mode: str = "interface",
    xyz_path: Path | None = None,
    metal_elements: set[str] | None = None,
    layer_tol_ang: float = DEFAULT_LAYER_TOL_A,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    # --- distributed mode params ---
    input_mode: str = "continuous",
    sp_root_dir: Path | str | None = None,
    sp_dir_pattern: str = "potential_t*_i*",
    sp_cube_filename: str = "sp_potential-v_hartree-1_0.cube",
    sp_out_filename: str = "sp.out",
) -> Path:
    """Multi-frame center-slab potential analysis.

    Reads all cube files matching *cube_pattern* (continuous mode) or
    all subdirectories matching *sp_dir_pattern* (distributed mode),
    computes the plane-averaged slab potential for each, and writes:

    - ``center_potential.csv`` + ``center_potential.png``
    - ``slab_center_and_interfaces.csv``

    Returns the CSV path.
    """
    workdir = Path(".").resolve()
    outdir = (output_dir or workdir).resolve()

    logger.info(
        "Center slab potential: mode=%s, thickness=%.1f A",
        input_mode, thickness_ang,
    )

    frames = _resolve_frames(
        input_mode=input_mode,
        cube_pattern=cube_pattern,
        workdir=workdir,
        xyz_path=xyz_path,
        center_mode=center_mode,
        metal_elements=metal_elements,
        layer_tol_ang=layer_tol_ang,
        sp_root_dir=sp_root_dir,
        sp_dir_pattern=sp_dir_pattern,
        sp_cube_filename=sp_cube_filename,
        sp_out_filename=sp_out_filename,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )

    return _analyze_center_from_frames(
        frames,
        thickness_ang=thickness_ang,
        center_mode=center_mode,
        metal_elements=metal_elements,
        layer_tol_ang=layer_tol_ang,
        output_dir=outdir,
        verbose=verbose,
    )


def fermi_energy_analysis(
    md_out_path: Path | None = None,
    *,
    output_dir: Path | None = None,
    fermi_unit: str = "au",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    # --- distributed mode params ---
    input_mode: str = "continuous",
    sp_root_dir: Path | str | None = None,
    sp_dir_pattern: str = "potential_t*_i*",
    sp_cube_filename: str = "sp_potential-v_hartree-1_0.cube",
    sp_out_filename: str = "sp.out",
) -> Path:
    """Fermi energy time-series analysis.

    In continuous mode, reads from a single CP2K md.out.
    In distributed mode, reads from sp.out in each subdirectory.

    Outputs ``fermi_energy.csv`` + ``fermi_energy.png``.
    Returns the CSV path.
    """
    outdir = (output_dir or Path(".")).resolve()

    if input_mode == "distributed":
        frames = _resolve_frames(
            input_mode="distributed",
            sp_root_dir=sp_root_dir,
            sp_dir_pattern=sp_dir_pattern,
            sp_cube_filename=sp_cube_filename,
            sp_out_filename=sp_out_filename,
            center_mode="cell",  # no interface detection needed for fermi-only
            frame_start=frame_start,
            frame_end=frame_end,
            frame_step=frame_step,
        )
        return _analyze_fermi_from_frames(
            frames, fermi_unit=fermi_unit, output_dir=outdir,
        )

    # Continuous mode: original logic (parse md.out directly)
    if md_out_path is None:
        raise ValueError("md_out_path is required for continuous mode")
    md_out_path = Path(md_out_path).resolve()
    if not md_out_path.exists():
        raise FileNotFoundError(f"md.out not found: {md_out_path}")

    from ._frame_source import _parse_md_out_fermi

    fermi_records = _parse_md_out_fermi(md_out_path)
    if not fermi_records:
        raise RuntimeError(f"No (step, Fermi energy) records parsed from: {md_out_path}")
    fermi_records = fermi_records[frame_start:frame_end:frame_step]
    logger.info("Fermi energy: %d records from %s", len(fermi_records), md_out_path)

    f_steps = np.array([r["step"] for r in fermi_records], dtype=int)
    f_time = np.array([r["time_fs"] if r["time_fs"] is not None else math.nan for r in fermi_records], dtype=float)
    f_raw = np.array([r["fermi_raw"] for r in fermi_records], dtype=float)
    f_vals = f_raw * HA_TO_EV if fermi_unit == "au" else f_raw
    f_cum = _cumulative_average(f_vals)

    outdir.mkdir(parents=True, exist_ok=True)
    f_rows: list[dict] = []
    for s, t, raw, v, c in zip(f_steps, f_time, f_raw, f_vals, f_cum, strict=True):
        f_rows.append({
            "step": int(s),
            "time_fs": None if math.isnan(float(t)) else float(t),
            "fermi_raw": float(raw),
            "fermi_ev": float(v),
            "fermi_cumavg_ev": float(c),
        })

    csv_path = outdir / DEFAULT_FERMI_ENERGY_CSV_NAME
    _write_csv(csv_path, f_rows, ["step", "time_fs", "fermi_raw", "fermi_ev", "fermi_cumavg_ev"])

    _plot_series_with_cumavg(
        outdir / DEFAULT_FERMI_ENERGY_PNG_NAME,
        x=f_steps.astype(float),
        y=f_vals,
        y_cum=f_cum,
        xlabel="MD step",
        ylabel="Fermi energy (eV)",
        title="Fermi energy from md.out",
    )

    return csv_path


def electrode_potential_analysis(
    cube_pattern: str = "",
    md_out_path: Path | None = None,
    *,
    output_dir: Path | None = None,
    thickness_ang: float = DEFAULT_THICKNESS_ANG,
    center_mode: str = "interface",
    xyz_path: Path | None = None,
    metal_elements: set[str] | None = None,
    layer_tol_ang: float = DEFAULT_LAYER_TOL_A,
    fermi_unit: str = "au",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    # --- distributed mode params ---
    input_mode: str = "continuous",
    sp_root_dir: Path | str | None = None,
    sp_dir_pattern: str = "potential_t*_i*",
    sp_cube_filename: str = "sp_potential-v_hartree-1_0.cube",
    sp_out_filename: str = "sp.out",
) -> Path:
    """Full electrode potential analysis (phi_center + E_Fermi -> U vs SHE).

    In distributed mode, both cube files and Fermi energies are loaded
    from the single-point subdirectories. In continuous mode, cube files
    are matched by glob and Fermi energies are parsed from md.out.

    Returns the U vs SHE CSV path.
    """
    outdir = (output_dir or Path(".")).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    if input_mode == "distributed":
        # Unified path: frames contain both cube data and fermi
        frames = _resolve_frames(
            input_mode="distributed",
            sp_root_dir=sp_root_dir,
            sp_dir_pattern=sp_dir_pattern,
            sp_cube_filename=sp_cube_filename,
            sp_out_filename=sp_out_filename,
            center_mode=center_mode,
            metal_elements=metal_elements,
            layer_tol_ang=layer_tol_ang,
            frame_start=frame_start,
            frame_end=frame_end,
            frame_step=frame_step,
            verbose=verbose,
        )
        # Write center potential sub-analysis
        _analyze_center_from_frames(
            frames,
            thickness_ang=thickness_ang,
            center_mode=center_mode,
            metal_elements=metal_elements,
            layer_tol_ang=layer_tol_ang,
            output_dir=outdir,
            verbose=verbose,
        )
        # Write fermi sub-analysis
        _analyze_fermi_from_frames(
            frames, fermi_unit=fermi_unit, output_dir=outdir,
        )
        # Merge: compute U vs SHE directly from frames
        return _compute_electrode_potential_from_frames(
            frames,
            thickness_ang=thickness_ang,
            center_mode=center_mode,
            metal_elements=metal_elements,
            layer_tol_ang=layer_tol_ang,
            fermi_unit=fermi_unit,
            output_dir=outdir,
        )

    # Continuous mode: delegate to sub-analyses then merge via CSV
    center_csv = center_slab_potential_analysis(
        cube_pattern,
        output_dir=outdir,
        thickness_ang=thickness_ang,
        center_mode=center_mode,
        xyz_path=xyz_path,
        metal_elements=metal_elements,
        layer_tol_ang=layer_tol_ang,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )
    fermi_csv = fermi_energy_analysis(
        md_out_path,
        output_dir=outdir,
        fermi_unit=fermi_unit,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
    )

    # Merge CSVs
    center_data = np.atleast_1d(np.genfromtxt(center_csv, delimiter=",", names=True, dtype=None, encoding="utf-8"))
    fermi_data = np.atleast_1d(np.genfromtxt(fermi_csv, delimiter=",", names=True, dtype=None, encoding="utf-8"))

    fermi_by_step = {int(r["step"]): float(r["fermi_ev"]) for r in fermi_data}

    u_rows: list[dict] = []
    u_vals: list[float] = []
    for row in center_data:
        s = int(row["step"])
        phi_ev = float(row["phi_center_ev"])
        if s not in fermi_by_step:
            continue
        ef_ev = fermi_by_step[s]
        u_v = -ef_ev + phi_ev + DP_A_H3O_W_EV - MU_HPLUS_G0_EV - DELTA_E_ZP_EV
        u_vals.append(u_v)
        u_rows.append({"step": s, "U_vs_SHE_V": float(u_v)})

    logger.info("Electrode potential: %d matched frames", len(u_rows))

    u_csv_path = outdir / DEFAULT_ELECTRODE_POTENTIAL_CSV_NAME
    if u_rows:
        u_vals_arr = np.array(u_vals, dtype=float)
        u_cum = _cumulative_average(u_vals_arr)
        for i, r in enumerate(u_rows):
            r["U_cumavg_V"] = float(u_cum[i])
        _write_csv(u_csv_path, u_rows, ["step", "U_vs_SHE_V", "U_cumavg_V"])
        _plot_series_with_cumavg(
            outdir / DEFAULT_ELECTRODE_POTENTIAL_PNG_NAME,
            x=np.array([r["step"] for r in u_rows], dtype=float),
            y=np.array([r["U_vs_SHE_V"] for r in u_rows], dtype=float),
            y_cum=np.array([r["U_cumavg_V"] for r in u_rows], dtype=float),
            xlabel="MD step",
            ylabel="U (V vs SHE)",
            title="Electrode potential U vs SHE (aligned at cube steps)",
        )

    return u_csv_path


def _compute_electrode_potential_from_frames(
    frames: list[PotentialFrame],
    *,
    thickness_ang: float,
    center_mode: str,
    metal_elements: set[str] | None,
    layer_tol_ang: float,
    fermi_unit: str,
    output_dir: Path,
) -> Path:
    """Compute U vs SHE directly from PotentialFrame list (distributed mode)."""
    use_interface = center_mode == "interface"
    metal_used = metal_elements
    cSHE_offset = DP_A_H3O_W_EV - MU_HPLUS_G0_EV - DELTA_E_ZP_EV

    u_rows: list[dict] = []
    u_vals: list[float] = []

    for frame in frames:
        if frame.fermi_raw is None:
            continue
        fermi_ev = frame.fermi_raw * HA_TO_EV if fermi_unit == "au" else frame.fermi_raw

        header = frame.header
        dz_bohr = float(np.linalg.norm(header.vz_bohr))
        lz_ang = header.nz * dz_bohr * BOHR_TO_ANG
        origin_z_ang = float(header.origin_bohr[2]) * BOHR_TO_ANG

        z_center_ang: Optional[float] = None
        if use_interface and frame.atoms is not None:
            if metal_used is None:
                metal_used = (
                    set(frame.atoms.get_chemical_symbols())
                    & set(TRANSITION_METAL_SYMBOLS)
                )
            if metal_used:
                detection = detect_interface_layers(
                    frame.atoms, metal_symbols=metal_used,
                    normal="c", layer_tol_A=layer_tol_ang,
                )
                iface = _extract_interface_geometry(detection, lz_ang)
                z_center_ang = origin_z_ang + float(iface["z_mid_ang"])

        phi_ev, _ = slab_average_potential_ev(
            header, frame.values, thickness_ang=thickness_ang,
            z_center_ang=z_center_ang,
        )
        u_v = -fermi_ev + phi_ev + cSHE_offset
        u_vals.append(u_v)
        u_rows.append({"step": frame.step, "U_vs_SHE_V": float(u_v)})

    logger.info("Electrode potential (distributed): %d frames", len(u_rows))

    u_csv_path = output_dir / DEFAULT_ELECTRODE_POTENTIAL_CSV_NAME
    if u_rows:
        u_vals_arr = np.array(u_vals, dtype=float)
        u_cum = _cumulative_average(u_vals_arr)
        for i, r in enumerate(u_rows):
            r["U_cumavg_V"] = float(u_cum[i])
        _write_csv(u_csv_path, u_rows, ["step", "U_vs_SHE_V", "U_cumavg_V"])
        _plot_series_with_cumavg(
            output_dir / DEFAULT_ELECTRODE_POTENTIAL_PNG_NAME,
            x=np.array([r["step"] for r in u_rows], dtype=float),
            y=np.array([r["U_vs_SHE_V"] for r in u_rows], dtype=float),
            y_cum=np.array([r["U_cumavg_V"] for r in u_rows], dtype=float),
            xlabel="MD step",
            ylabel="U (V vs SHE)",
            title="Electrode potential U vs SHE (distributed SP)",
        )

    return u_csv_path


def thickness_sensitivity_analysis(
    cube_pattern: str = "",
    md_out_path: Path | None = None,
    *,
    output_dir: Path | None = None,
    thickness_start: float = 3.5,
    thickness_end: float = 15.0,
    thickness_step: float = 0.5,
    center_mode: str = "interface",
    xyz_path: Path | None = None,
    metal_elements: set[str] | None = None,
    layer_tol_ang: float = DEFAULT_LAYER_TOL_A,
    fermi_unit: str = "au",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
    # --- distributed mode params ---
    input_mode: str = "continuous",
    sp_root_dir: Path | str | None = None,
    sp_dir_pattern: str = "potential_t*_i*",
    sp_cube_filename: str = "sp_potential-v_hartree-1_0.cube",
    sp_out_filename: str = "sp.out",
) -> Path:
    """Sweep slab-averaging thickness and report mean/std of U vs SHE.

    Returns the CSV path.
    """
    outdir = (output_dir or Path(".")).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Load all frames
    frames = _resolve_frames(
        input_mode=input_mode,
        cube_pattern=cube_pattern,
        workdir=Path(".").resolve() if input_mode == "continuous" else None,
        md_out_path=md_out_path,
        xyz_path=xyz_path,
        sp_root_dir=sp_root_dir,
        sp_dir_pattern=sp_dir_pattern,
        sp_cube_filename=sp_cube_filename,
        sp_out_filename=sp_out_filename,
        center_mode=center_mode,
        metal_elements=metal_elements,
        layer_tol_ang=layer_tol_ang,
        fermi_unit=fermi_unit,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )

    # Build per-frame cache: (header, values, z_center_ang, fermi_ev)
    use_interface = center_mode == "interface"
    metal_used = metal_elements
    frame_cache: list[tuple] = []

    for frame in frames:
        if frame.fermi_raw is None:
            continue
        fermi_ev = frame.fermi_raw * HA_TO_EV if fermi_unit == "au" else frame.fermi_raw

        header = frame.header
        dz_bohr = float(np.linalg.norm(header.vz_bohr))
        lz_ang = header.nz * dz_bohr * BOHR_TO_ANG
        origin_z_ang = float(header.origin_bohr[2]) * BOHR_TO_ANG

        z_center_ang: Optional[float] = None
        if use_interface and frame.atoms is not None:
            if metal_used is None:
                metal_used = (
                    set(frame.atoms.get_chemical_symbols())
                    & set(TRANSITION_METAL_SYMBOLS)
                )
            if metal_used:
                detection = detect_interface_layers(
                    frame.atoms, metal_symbols=metal_used,
                    normal="c", layer_tol_A=layer_tol_ang,
                )
                iface = _extract_interface_geometry(detection, lz_ang)
                z_center_ang = origin_z_ang + float(iface["z_mid_ang"])

        frame_cache.append((header, frame.values, z_center_ang, fermi_ev))

    if not frame_cache:
        raise RuntimeError("No frames with matching cube + Fermi energy data.")

    # Sweep thickness values
    cSHE_offset = DP_A_H3O_W_EV - MU_HPLUS_G0_EV - DELTA_E_ZP_EV
    thicknesses = np.arange(thickness_start, thickness_end + thickness_step * 0.5, thickness_step)
    logger.info("Thickness sensitivity: sweeping %d values", len(thicknesses))

    rows: list[dict] = []
    means: list[float] = []
    spatial_stds: list[float] = []

    thick_iter: Iterable = thicknesses
    if verbose:
        from tqdm import tqdm
        thick_iter = tqdm(thicknesses, desc="Thickness sweep", unit="t", ascii=" =")

    for thick in thick_iter:
        thick = float(thick)
        frame_u_vals: list[float] = []
        frame_phi_stds: list[float] = []
        for hdr, vals, zc, ef_ev in frame_cache:
            phi_ev, info = slab_average_potential_ev(
                hdr, vals, thickness_ang=thick, z_center_ang=zc,
            )
            u_v = -ef_ev + phi_ev + cSHE_offset
            frame_u_vals.append(u_v)
            frame_phi_stds.append(info["phi_z_std_ev"])
        u_arr = np.array(frame_u_vals, dtype=float)
        m = float(np.mean(u_arr))
        spatial_std = float(np.mean(frame_phi_stds))
        means.append(m)
        spatial_stds.append(spatial_std)
        rows.append({
            "thickness_ang": round(thick, 4),
            "mean_U_vs_SHE_V": m,
            "mean_phi_z_spatial_std_eV": spatial_std,
            "n_frames": len(frame_u_vals),
        })

    # Write CSV
    csv_path = outdir / DEFAULT_THICKNESS_SENSITIVITY_CSV_NAME
    _write_csv(csv_path, rows, [
        "thickness_ang", "mean_U_vs_SHE_V", "mean_phi_z_spatial_std_eV", "n_frames",
    ])

    # Dual-axis plot
    png_path = outdir / DEFAULT_THICKNESS_SENSITIVITY_PNG_NAME
    png_path.parent.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax1 = plt.subplots(figsize=(9, 4.8), dpi=160)
    color_mean = "tab:blue"
    ax1.plot(thicknesses, means, "o-", color=color_mean, lw=1.5, markersize=4)
    ax1.set_xlabel("Slab thickness (Å)")
    ax1.set_ylabel("Mean U vs SHE (V)", color=color_mean)
    ax1.tick_params(axis="y", labelcolor=color_mean)

    ax2 = ax1.twinx()
    color_std = "tab:red"
    ax2.plot(thicknesses, spatial_stds, "s--", color=color_std, lw=1.5, markersize=4)
    ax2.set_ylabel("Spatial std of φ(z) in slab (eV)", color=color_std)
    ax2.tick_params(axis="y", labelcolor=color_std)

    ax1.set_title("Electrode potential U vs SHE — thickness sensitivity")
    ax1.grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig(png_path)
    plt.close(fig)

    return csv_path
