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
import math
import re
from pathlib import Path
from typing import Iterable, Optional

import numpy as np

from ..utils.config import (
    BOHR_TO_ANG,
    DP_A_H3O_W_EV,
    DELTA_E_ZP_EV,
    HA_TO_EV,
    MU_HPLUS_G0_EV,
    TRANSITION_METAL_SYMBOLS,
)
from ..utils.CubeParser import (
    extract_step_from_cube_filename,
    read_cube_header_and_values,
    slab_average_potential_ev,
)
from ..utils.ClusterUtils import cluster_1d_periodic, find_largest_gap_periodic, gap_midpoint_periodic
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
# Private helpers (migrated from temp/analyze_potential.py)
# ---------------------------------------------------------------------------

XYZ_STEP_RE = re.compile(r"\bi\s*=\s*(\d+)\b")
STEP_RE = re.compile(r"STEP NUMBER\s*=\s*(\d+)")
TIME_RE = re.compile(r"TIME\s*\[fs\]\s*=\s*([+-]?\d+(?:\.\d*)?(?:[EeDd][+-]?\d+)?)")
FERMI_RE = re.compile(r"Fermi energy:\s*([+-]?\d+(?:\.\d*)?(?:[EeDd][+-]?\d+)?)")


def _float(s: str) -> float:
    return float(s.replace("D", "E").replace("d", "e"))


def _parse_md_out_fermi(md_out_path: Path) -> list[dict]:
    """Parse ``(step, time_fs, fermi_raw)`` records from CP2K md.out."""
    records: list[dict] = []
    fermi_pending_raw: Optional[float] = None
    last_rec: Optional[dict] = None

    with md_out_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            m = FERMI_RE.search(line)
            if m:
                fermi_pending_raw = _float(m.group(1))
                continue

            m = STEP_RE.search(line)
            if m:
                step = int(m.group(1))
                rec = {"step": step, "time_fs": None, "fermi_raw": None}
                if fermi_pending_raw is not None:
                    rec["fermi_raw"] = fermi_pending_raw
                    fermi_pending_raw = None
                records.append(rec)
                last_rec = rec
                continue

            m = TIME_RE.search(line)
            if m and last_rec is not None and last_rec["time_fs"] is None:
                last_rec["time_fs"] = _float(m.group(1))

    return [r for r in records if r["fermi_raw"] is not None]


def _read_xyz_metal_z_for_steps(
    xyz_path: Path,
    steps: set[int],
    *,
    metal_elements: Optional[set[str]] = None,
) -> tuple[dict[int, np.ndarray], set[str]]:
    """Stream-parse a CP2K xyz trajectory and extract metal-atom z (Å) for given steps."""
    if not xyz_path.exists():
        raise FileNotFoundError(xyz_path)

    steps = {int(s) for s in steps}
    out: dict[int, np.ndarray] = {}
    inferred_metal: Optional[set[str]] = set(metal_elements) if metal_elements is not None else None
    inferred_done = inferred_metal is not None

    with xyz_path.open("r", encoding="utf-8", errors="replace") as f:
        while True:
            natoms_line = f.readline()
            if not natoms_line:
                break
            natoms_line = natoms_line.strip()
            if not natoms_line:
                continue
            try:
                natoms = int(natoms_line.split()[0])
            except ValueError:
                continue

            comment = f.readline()
            if not comment:
                break
            m = XYZ_STEP_RE.search(comment)
            step = int(m.group(1)) if m else None
            need_this = (step is not None) and (step in steps)

            if not inferred_done:
                symbols: set[str] = set()
                z_by_symbol: dict[str, list[float]] = {}
                for _ in range(natoms):
                    line = f.readline()
                    if not line:
                        break
                    parts = line.split()
                    if len(parts) < 4:
                        continue
                    sym = parts[0]
                    symbols.add(sym)
                    if need_this:
                        z_by_symbol.setdefault(sym, []).append(_float(parts[3]))

                inferred_metal = set(symbols) & set(TRANSITION_METAL_SYMBOLS)
                inferred_done = True

                if need_this and step is not None:
                    z_list: list[float] = []
                    for sym in inferred_metal:
                        z_list.extend(z_by_symbol.get(sym, []))
                    out[int(step)] = np.array(z_list, dtype=float)
                continue

            if need_this:
                z_list = []
                assert inferred_metal is not None
                for _ in range(natoms):
                    line = f.readline()
                    if not line:
                        break
                    parts = line.split()
                    if len(parts) < 4:
                        continue
                    if parts[0] in inferred_metal:
                        z_list.append(_float(parts[3]))
                if step is not None:
                    out[int(step)] = np.array(z_list, dtype=float)
            else:
                for _ in range(natoms):
                    if not f.readline():
                        break

    if inferred_metal is None:
        inferred_metal = set()
    return out, inferred_metal


def _detect_metal_water_interfaces_z_ang(
    metal_z_ang: np.ndarray,
    lz_ang: float,
    *,
    layer_tol_ang: float = 0.6,
) -> dict:
    """Identify two metal/water interfaces from metal z-coordinates (Å).

    Uses ``cluster_1d_periodic`` + ``find_largest_gap_periodic`` from ClusterUtils.
    """
    if metal_z_ang.size == 0:
        raise ValueError("Empty metal_z array")
    if lz_ang <= 0:
        raise ValueError(f"Invalid Lz: {lz_ang}")

    clusters = cluster_1d_periodic(metal_z_ang, period=lz_ang, tol=layer_tol_ang)
    centers = np.array([c for c, _ in clusters], dtype=float)

    if centers.size < 2:
        raise ValueError(f"Need >= 2 metal layers to detect two interfaces; got {centers.size}")

    low_idx, high_idx, water_gap = find_largest_gap_periodic(centers, period=lz_ang)
    z_lower = float(centers[low_idx])
    z_upper = float(centers[high_idx])
    z_mid = gap_midpoint_periodic(z_lower, z_upper, period=lz_ang)

    return {
        "z_lower_ang": z_lower,
        "z_upper_ang": z_upper,
        "z_mid_ang": z_mid,
        "water_gap_ang": water_gap,
        "n_layers": int(centers.size),
    }


def _cumulative_average(values: np.ndarray) -> np.ndarray:
    csum = np.cumsum(values, dtype=float)
    return csum / np.arange(1, values.size + 1, dtype=float)


def _write_csv(path: Path, rows: Iterable[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _plot_series_with_cumavg(
    png_path: Path,
    x: np.ndarray,
    y: np.ndarray,
    y_cum: np.ndarray,
    xlabel: str,
    ylabel: str,
    title: str,
) -> None:
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
    if s is None:
        return None
    raw = [p.strip() for p in s.split(",")]
    out = {p for p in raw if p}
    return out or None


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def center_slab_potential_analysis(
    cube_pattern: str,
    *,
    output_dir: Path | None = None,
    thickness_ang: float = DEFAULT_THICKNESS_ANG,
    center_mode: str = "interface",
    xyz_path: Path | None = None,
    metal_elements: set[str] | None = None,
    layer_tol_ang: float = 0.6,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> Path:
    """Multi-frame center-slab potential analysis.

    Reads all cube files matching *cube_pattern* in the current working
    directory, computes the plane-averaged slab potential for each, and
    writes:

    - ``center_potential.csv`` + ``center_potential.png``
    - ``slab_center_and_interfaces.csv``

    Returns the CSV path.
    """
    workdir = Path(".").resolve()
    outdir = (output_dir or workdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    cube_paths = [Path(p) for p in sorted(workdir.glob(cube_pattern))]
    if not cube_paths:
        raise FileNotFoundError(f"No cube files matched pattern: {cube_pattern!r} in {workdir}")
    cube_paths = cube_paths[frame_start:frame_end:frame_step]

    # Optional interface detection from xyz
    use_interface_center = center_mode == "interface"
    metal_z_by_step: dict[int, np.ndarray] = {}
    metal_elements_used: Optional[set[str]] = None

    if use_interface_center:
        if xyz_path is None:
            xyz_path = workdir / "md-pos-1.xyz"
        xyz_path = Path(xyz_path).resolve()
        if not xyz_path.exists():
            raise FileNotFoundError(
                f"Interface center requested but xyz not found: {xyz_path}"
            )
        needed_steps = {
            s for s in (extract_step_from_cube_filename(cp) for cp in cube_paths) if s is not None
        }
        metal_z_by_step, inferred_metal = _read_xyz_metal_z_for_steps(
            xyz_path, needed_steps, metal_elements=metal_elements,
        )
        metal_elements_used = metal_elements or inferred_metal
        if not metal_elements_used:
            raise RuntimeError(
                "Interface detection failed: no transition metals inferred from xyz. "
                "Specify metal_elements explicitly."
            )

    cube_rows: list[dict] = []
    cube_steps: list[int] = []
    cube_vals_ev: list[float] = []
    slab_info_first: Optional[dict] = None
    center_rows: list[dict] = []

    cube_iter: Iterable[Path] = cube_paths
    if verbose:
        from tqdm import tqdm
        cube_iter = tqdm(cube_paths, desc="Center slab potential", unit="cube")

    for cp in cube_iter:
        step = extract_step_from_cube_filename(cp)
        if step is None:
            step = len(cube_steps)

        header, values = read_cube_header_and_values(cp)
        dz_bohr = float(np.linalg.norm(header.vz_bohr))
        lz_ang = header.nz * dz_bohr * BOHR_TO_ANG
        origin_z_ang = float(header.origin_bohr[2]) * BOHR_TO_ANG

        z_center_ang: Optional[float] = None
        iface: Optional[dict] = None
        center_source = "cell"

        if use_interface_center and int(step) in metal_z_by_step:
            iface = _detect_metal_water_interfaces_z_ang(
                metal_z_by_step[int(step)],
                lz_ang,
                layer_tol_ang=layer_tol_ang,
            )
            z_center_ang = origin_z_ang + float(iface["z_mid_ang"])
            center_source = "interface"

        phi_ev, info = slab_average_potential_ev(
            header, values, thickness_ang=thickness_ang, z_center_ang=z_center_ang,
        )
        if slab_info_first is None:
            slab_info_first = info
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

    for s, v, c in zip(cube_steps_arr, cube_vals_arr, cube_cum, strict=True):
        cube_rows.append({"step": int(s), "phi_center_ev": float(v), "phi_center_cumavg_ev": float(c)})

    csv_path = outdir / DEFAULT_CENTER_POTENTIAL_CSV_NAME
    _write_csv(csv_path, cube_rows, ["step", "phi_center_ev", "phi_center_cumavg_ev"])

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


def fermi_energy_analysis(
    md_out_path: Path,
    *,
    output_dir: Path | None = None,
    fermi_unit: str = "au",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
) -> Path:
    """Fermi energy time-series analysis from CP2K md.out.

    Outputs ``fermi_energy.csv`` + ``fermi_energy.png``.
    Returns the CSV path.
    """
    md_out_path = Path(md_out_path).resolve()
    if not md_out_path.exists():
        raise FileNotFoundError(f"md.out not found: {md_out_path}")

    outdir = (output_dir or Path(".")).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    fermi_records = _parse_md_out_fermi(md_out_path)
    if not fermi_records:
        raise RuntimeError(f"No (step, Fermi energy) records parsed from: {md_out_path}")
    fermi_records = fermi_records[frame_start:frame_end:frame_step]

    f_steps = np.array([r["step"] for r in fermi_records], dtype=int)
    f_time = np.array([r["time_fs"] if r["time_fs"] is not None else math.nan for r in fermi_records], dtype=float)
    f_raw = np.array([r["fermi_raw"] for r in fermi_records], dtype=float)
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


def electrode_potential_analysis(
    cube_pattern: str,
    md_out_path: Path,
    *,
    output_dir: Path | None = None,
    thickness_ang: float = DEFAULT_THICKNESS_ANG,
    center_mode: str = "interface",
    xyz_path: Path | None = None,
    metal_elements: set[str] | None = None,
    layer_tol_ang: float = 0.6,
    fermi_unit: str = "au",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> Path:
    """Full electrode potential analysis (φ_center + E_Fermi → U vs SHE).

    Internally calls ``center_slab_potential_analysis`` and
    ``fermi_energy_analysis``, then merges to compute U vs SHE.

    Returns the U vs SHE CSV path.
    """
    outdir = (output_dir or Path(".")).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Run sub-analyses
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

    # Read back CSVs to merge (atleast_1d guards against 0-d arrays from single-row CSV)
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


def thickness_sensitivity_analysis(
    cube_pattern: str,
    md_out_path: Path,
    *,
    output_dir: Path | None = None,
    thickness_start: float = 3.5,
    thickness_end: float = 15.0,
    thickness_step: float = 0.5,
    center_mode: str = "interface",
    xyz_path: Path | None = None,
    metal_elements: set[str] | None = None,
    layer_tol_ang: float = 0.6,
    fermi_unit: str = "au",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> Path:
    """Sweep slab-averaging thickness and report mean/std of U vs SHE.

    For each thickness in ``[thickness_start, thickness_end]`` (step
    ``thickness_step``), computes the electrode potential
    ``U = -E_Fermi + φ_center + ΔΨ - μ(H⁺) - ΔE_ZP`` for every
    cube frame, then reports the ensemble mean and standard deviation.

    Outputs:
    - ``thickness_sensitivity.csv``
    - ``thickness_sensitivity.png`` (dual-axis: left=mean U, right=std U)

    Returns the CSV path.
    """
    workdir = Path(".").resolve()
    outdir = (output_dir or workdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    cube_paths = [Path(p) for p in sorted(workdir.glob(cube_pattern))]
    if not cube_paths:
        raise FileNotFoundError(f"No cube files matched pattern: {cube_pattern!r} in {workdir}")
    cube_paths = cube_paths[frame_start:frame_end:frame_step]

    # --- Parse Fermi energies from md.out ---
    md_out_path = Path(md_out_path).resolve()
    if not md_out_path.exists():
        raise FileNotFoundError(f"md.out not found: {md_out_path}")
    fermi_records = _parse_md_out_fermi(md_out_path)
    if not fermi_records:
        raise RuntimeError(f"No Fermi energy records parsed from: {md_out_path}")
    fermi_records = fermi_records[frame_start:frame_end:frame_step]
    fermi_raw_by_step: dict[int, float] = {
        int(r["step"]): float(r["fermi_raw"]) for r in fermi_records
    }

    # --- Interface detection (same logic as center_slab_potential_analysis) ---
    use_interface_center = center_mode == "interface"
    metal_z_by_step: dict[int, np.ndarray] = {}

    if use_interface_center:
        if xyz_path is None:
            xyz_path = workdir / "md-pos-1.xyz"
        xyz_path = Path(xyz_path).resolve()
        if not xyz_path.exists():
            raise FileNotFoundError(
                f"Interface center requested but xyz not found: {xyz_path}"
            )
        needed_steps = {
            s for s in (extract_step_from_cube_filename(cp) for cp in cube_paths) if s is not None
        }
        metal_z_by_step, inferred_metal = _read_xyz_metal_z_for_steps(
            xyz_path, needed_steps, metal_elements=metal_elements,
        )
        metal_elements_used = metal_elements or inferred_metal
        if not metal_elements_used:
            raise RuntimeError(
                "Interface detection failed: no transition metals inferred from xyz. "
                "Specify metal_elements explicitly."
            )

    # --- Cache per-frame cube data + z_center + fermi_ev ---
    frame_cache: list[tuple] = []  # (header, values, z_center_ang, fermi_ev)
    for cp in cube_paths:
        step = extract_step_from_cube_filename(cp)
        if step is None:
            step = len(frame_cache)

        # Skip frames without matching Fermi energy
        if int(step) not in fermi_raw_by_step:
            continue

        fermi_raw = fermi_raw_by_step[int(step)]
        fermi_ev = fermi_raw * HA_TO_EV if fermi_unit == "au" else fermi_raw

        header, values = read_cube_header_and_values(cp)
        dz_bohr = float(np.linalg.norm(header.vz_bohr))
        lz_ang = header.nz * dz_bohr * BOHR_TO_ANG
        origin_z_ang = float(header.origin_bohr[2]) * BOHR_TO_ANG

        z_center_ang: Optional[float] = None
        if use_interface_center and int(step) in metal_z_by_step:
            iface = _detect_metal_water_interfaces_z_ang(
                metal_z_by_step[int(step)],
                lz_ang,
                layer_tol_ang=layer_tol_ang,
            )
            z_center_ang = origin_z_ang + float(iface["z_mid_ang"])

        frame_cache.append((header, values, z_center_ang, fermi_ev))

    if not frame_cache:
        raise RuntimeError("No frames with matching cube + Fermi energy data.")

    # --- Sweep thickness values ---
    cSHE_offset = DP_A_H3O_W_EV - MU_HPLUS_G0_EV - DELTA_E_ZP_EV
    thicknesses = np.arange(thickness_start, thickness_end + thickness_step * 0.5, thickness_step)

    rows: list[dict] = []
    means: list[float] = []
    spatial_stds: list[float] = []

    thick_iter: Iterable = thicknesses
    if verbose:
        from tqdm import tqdm
        thick_iter = tqdm(thicknesses, desc="Thickness sweep", unit="t")

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

    # --- Write CSV ---
    csv_path = outdir / DEFAULT_THICKNESS_SENSITIVITY_CSV_NAME
    _write_csv(csv_path, rows, [
        "thickness_ang", "mean_U_vs_SHE_V", "mean_phi_z_spatial_std_eV", "n_frames",
    ])

    # --- Dual-axis plot ---
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
