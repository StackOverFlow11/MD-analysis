"""
CP2K DFTMD post-processing:
- Center water slab plane-averaged Hartree potential from *.cube snapshots
- Fermi energy time series from md.out

Default usage (run inside a CP2K output directory):
  python analyze_potential.py

Outputs (by default) in ./analysis_outputs:
  - center_potential.csv / center_potential.png
  - fermi_energy.csv / fermi_energy.png

Optional:
  --compute-u    Compute electrode potential U vs SHE using computational SHE formula
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np


HA_TO_EV = 27.211386245988  # CODATA 2018 (commonly used); enough precision for this workflow
BOHR_TO_ANG = 0.529177210903

# computational SHE constants (eV), from method_intro/*
DP_A_H3O_W_EV = 15.35
MU_HPLUS_G0_EV = 15.81
DELTA_E_ZP_EV = 0.35

# Transition metals (d-block) used for default metal inference
TRANSITION_METALS = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
}

@dataclass(frozen=True)
class CubeHeader:
    natoms: int
    origin_bohr: np.ndarray  # (3,)
    nx: int
    ny: int
    nz: int
    vx_bohr: np.ndarray  # (3,) voxel step vector
    vy_bohr: np.ndarray  # (3,)
    vz_bohr: np.ndarray  # (3,)


def _float(s: str) -> float:
    # CP2K sometimes uses Fortran 'D' exponents in some outputs; be tolerant.
    return float(s.replace("D", "E").replace("d", "e"))


def _extract_step_from_cube_filename(path: Path) -> Optional[int]:
    # Typical: md-POTENTIAL-v_hartree-1_1050.cube  -> 1050
    m = re.search(r"_(\d+)\.cube$", path.name)
    if m:
        return int(m.group(1))
    return None


def read_cube_header_and_values(path: Path) -> tuple[CubeHeader, np.ndarray]:
    """
    Read a Gaussian cube file written by CP2K and return header + scalar field values.
    Assumes data ordering: x fastest, then y, then z (standard cube convention).
    Returns values as a flat array of length nx*ny*nz.
    """
    with path.open("r", encoding="utf-8", errors="replace") as f:
        _ = f.readline()  # comment 1
        _ = f.readline()  # comment 2

        t = f.readline().split()
        if len(t) < 4:
            raise ValueError(f"Bad cube header in {path}: natoms/origin line too short")
        natoms = int(t[0])
        origin_bohr = np.array([_float(t[1]), _float(t[2]), _float(t[3])], dtype=float)

        t = f.readline().split()
        nx = int(t[0])
        vx_bohr = np.array([_float(t[1]), _float(t[2]), _float(t[3])], dtype=float)

        t = f.readline().split()
        ny = int(t[0])
        vy_bohr = np.array([_float(t[1]), _float(t[2]), _float(t[3])], dtype=float)

        t = f.readline().split()
        nz = int(t[0])
        vz_bohr = np.array([_float(t[1]), _float(t[2]), _float(t[3])], dtype=float)

        # Skip atom lines
        for _ in range(natoms):
            f.readline()

        rest = f.read()
        if "D" in rest or "d" in rest:
            rest = rest.replace("D", "E").replace("d", "e")
        values = np.fromstring(rest, sep=" ", dtype=float)

    expected = nx * ny * nz
    if values.size != expected:
        raise ValueError(
            f"Cube data size mismatch in {path}: got {values.size}, expected {expected} (nx,ny,nz={nx},{ny},{nz})"
        )

    header = CubeHeader(
        natoms=natoms,
        origin_bohr=origin_bohr,
        nx=nx,
        ny=ny,
        nz=nz,
        vx_bohr=vx_bohr,
        vy_bohr=vy_bohr,
        vz_bohr=vz_bohr,
    )
    return header, values


def slab_average_potential_ev_from_cube(
    header: CubeHeader,
    values: np.ndarray,
    thickness_ang: float,
    *,
    z_center_ang: Optional[float] = None,
) -> tuple[float, dict]:
    """
    Compute slab-average potential from a cube snapshot.

    The cube is assumed to be aligned such that the 3rd axis corresponds to the surface normal (z).
    Steps:
      - reshape to (nz, ny, nx)
      - plane-average in xy for each z slice
      - average within a slab of thickness `thickness_ang` centered at `z_center_ang`
        (if None: use cell geometric center)

    Returns:
      - potential (eV) (i.e., e*phi) averaged in the center slab
      - small info dict for logging/debugging
    """
    dz_bohr = float(np.linalg.norm(header.vz_bohr))
    if dz_bohr <= 0:
        raise ValueError(f"Invalid dz in cube header: |vz|={dz_bohr} bohr")

    # NOTE (CP2K cube ordering):
    # CP2K's V_HARTREE_CUBE is written with z as the fastest-running index (then y, then x).
    # Therefore, the flat value array should be reshaped as (nx, ny, nz) in C-order.
    field = values.reshape((header.nx, header.ny, header.nz))
    phi_z_ha = field.mean(axis=(0, 1))  # (nz,)

    thickness_bohr = thickness_ang / BOHR_TO_ANG
    half_thickness_bohr = 0.5 * thickness_bohr

    z0_bohr = float(header.origin_bohr[2])
    lz_bohr = header.nz * dz_bohr
    z_coords_rel = (np.arange(header.nz, dtype=float) + 0.5) * dz_bohr  # in [0, Lz)

    if z_center_ang is None:
        z_center_rel = 0.5 * lz_bohr
        z_center_source = "cell"
    else:
        z_center_bohr = float(z_center_ang) / BOHR_TO_ANG
        z_center_rel = float((z_center_bohr - z0_bohr) % lz_bohr)
        z_center_source = "interface"

    # Periodic distance to center (robust if center lies close to the box boundary)
    dist = np.abs(z_coords_rel - z_center_rel)
    dist = np.minimum(dist, lz_bohr - dist)
    mask = dist <= half_thickness_bohr
    if not np.any(mask):
        raise ValueError(
            f"Thickness {thickness_ang} A selects 0 z-slices (dz={dz_bohr * BOHR_TO_ANG:.4f} A, nz={header.nz})"
        )

    phi_center_ha = float(phi_z_ha[mask].mean())
    # Convert to eV: in atomic units, a potential of 1 Ha/e corresponds numerically to 1 Ha.
    # Multiplying by HA_TO_EV yields eV (and also V numerically for a +e test charge).
    phi_center_ev = phi_center_ha * HA_TO_EV

    info = {
        "nx": header.nx,
        "ny": header.ny,
        "nz": header.nz,
        "dz_ang": dz_bohr * BOHR_TO_ANG,
        "n_slices": int(mask.sum()),
        "z_center_source": z_center_source,
        "z_center_rel_bohr": z_center_rel,
        "thickness_ang": thickness_ang,
    }
    return phi_center_ev, info


def cube_center_slab_potential_ev(path: Path, thickness_ang: float) -> tuple[float, dict]:
    """
    Backward-compatible wrapper: slab average centered at the cell geometric center.
    """
    header, values = read_cube_header_and_values(path)
    return slab_average_potential_ev_from_cube(header, values, thickness_ang=thickness_ang, z_center_ang=None)


XYZ_STEP_RE = re.compile(r"\bi\s*=\s*(\d+)\b")


def _parse_csv_symbols(s: Optional[str]) -> Optional[set[str]]:
    if s is None:
        return None
    raw = [p.strip() for p in s.split(",")]
    out = {p for p in raw if p}
    return out or None


def read_xyz_metal_z_for_steps(
    xyz_path: Path,
    steps: set[int],
    *,
    metal_elements: Optional[set[str]] = None,
) -> tuple[dict[int, np.ndarray], set[str]]:
    """
    Stream-parse a CP2K xyz trajectory (e.g., md-pos-1.xyz) and extract metal-atom z (Angstrom)
    for frames whose MD step is in `steps`.

    Notes:
    - Step id is parsed from the 2nd line of each frame via pattern: 'i = <step>'.
    - If `metal_elements` is None, we infer it from the *first* frame as (elements) ∩ TRANSITION_METALS.
    """
    if not xyz_path.exists():
        raise FileNotFoundError(xyz_path)

    steps = {int(s) for s in steps}
    out: dict[int, np.ndarray] = {}

    inferred_metal: Optional[set[str]] = (set(metal_elements) if metal_elements is not None else None)
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
                # tolerate occasional junk lines
                continue

            comment = f.readline()
            if not comment:
                break
            m = XYZ_STEP_RE.search(comment)
            step = int(m.group(1)) if m else None
            need_this = (step is not None) and (step in steps)

            if not inferred_done:
                # For the first frame, infer element types and (if needed) extract z by symbol.
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

                inferred_metal = set(symbols) & set(TRANSITION_METALS)
                inferred_done = True

                if need_this and step is not None:
                    z_list: list[float] = []
                    for sym in inferred_metal:
                        z_list.extend(z_by_symbol.get(sym, []))
                    out[int(step)] = np.array(z_list, dtype=float)
                continue

            # inferred already (or user-provided): only parse z for frames we need
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
                    sym = parts[0]
                    if sym in inferred_metal:
                        z_list.append(_float(parts[3]))
                if step is not None:
                    out[int(step)] = np.array(z_list, dtype=float)
            else:
                # fast-skip atom lines
                for _ in range(natoms):
                    if not f.readline():
                        break

    if inferred_metal is None:
        inferred_metal = set()
    return out, inferred_metal


def metal_layers_z_ang(
    metal_z_ang: np.ndarray,
    lz_ang: float,
    *,
    layer_tol_ang: float = 0.6,
) -> np.ndarray:
    """
    Cluster metal atoms into layers along z (periodic in [0, Lz)) using a simple 1D tolerance.

    Returns sorted layer centers in [0, Lz).
    """
    if metal_z_ang.size == 0:
        raise ValueError("Empty metal_z array")
    if lz_ang <= 0:
        raise ValueError(f"Invalid Lz: {lz_ang}")

    z = np.mod(metal_z_ang.astype(float), lz_ang)
    z_sorted = np.sort(z)

    layers: list[list[float]] = []
    cur: list[float] = [float(z_sorted[0])]
    for v in z_sorted[1:]:
        v = float(v)
        if v - cur[-1] <= layer_tol_ang:
            cur.append(v)
        else:
            layers.append(cur)
            cur = [v]
    layers.append(cur)

    # Merge first/last if the layer crosses the periodic boundary
    if len(layers) > 1:
        wrap_gap = (layers[0][0] + lz_ang) - layers[-1][-1]
        if wrap_gap <= layer_tol_ang:
            merged = layers[-1] + [v + lz_ang for v in layers[0]]
            merged_center = float(np.mean(merged) % lz_ang)
            # replace last with merged center, drop first
            layers = layers[1:-1] + [[merged_center]]

    centers = np.array([float(np.mean(layer)) for layer in layers], dtype=float)
    centers = np.mod(centers, lz_ang)
    centers.sort()
    return centers


def detect_metal_water_interfaces_z_ang(
    metal_z_ang: np.ndarray,
    lz_ang: float,
    *,
    layer_tol_ang: float = 0.6,
) -> dict:
    """
    Identify the two metal/water interfaces from metal atom z coordinates (Angstrom).

    Procedure:
    - cluster metal atoms into z-layers with tolerance `layer_tol_ang`
    - compute circular gaps between adjacent layers (PBC)
    - choose the *largest* gap as the water region
    - the two layers bounding this gap are the two metal/water interfaces
    - midpoint along that water segment is used as the water-region center
    """
    centers = metal_layers_z_ang(metal_z_ang, lz_ang, layer_tol_ang=layer_tol_ang)
    if centers.size < 2:
        raise ValueError(f"Need >= 2 metal layers to detect two interfaces; got {centers.size}")

    gaps = np.diff(np.concatenate([centers, [centers[0] + lz_ang]]))
    i = int(np.argmax(gaps))
    z_lower = float(centers[i])
    z_upper_unwrapped = float(centers[i + 1] if i + 1 < centers.size else centers[0] + lz_ang)
    water_gap = float(z_upper_unwrapped - z_lower)
    z_mid = float((z_lower + 0.5 * water_gap) % lz_ang)

    return {
        "z_lower_ang": z_lower,
        "z_upper_ang": float(z_upper_unwrapped % lz_ang),
        "z_mid_ang": z_mid,
        "water_gap_ang": water_gap,
        "n_layers": int(centers.size),
    }


STEP_RE = re.compile(r"STEP NUMBER\s*=\s*(\d+)")
TIME_RE = re.compile(r"TIME\s*\[fs\]\s*=\s*([+-]?\d+(?:\.\d*)?(?:[EeDd][+-]?\d+)?)")
FERMI_RE = re.compile(r"Fermi energy:\s*([+-]?\d+(?:\.\d*)?(?:[EeDd][+-]?\d+)?)")


def parse_md_out_fermi(md_out_path: Path) -> list[dict]:
    """
    Parse CP2K md.out and return a list of records:
      {step, time_fs, fermi_raw}

    NOTE: In many CP2K outputs, 'Fermi energy' appears in the SCF block *before*
    the MD summary that contains 'STEP NUMBER'. We assign the most recently seen
    Fermi energy to the next encountered step number (robust for this common layout).
    """
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

    # Drop any steps where we failed to assign a fermi value (should be rare)
    records = [r for r in records if r["fermi_raw"] is not None]
    return records


def cumulative_average(values: np.ndarray) -> np.ndarray:
    csum = np.cumsum(values, dtype=float)
    return csum / np.arange(1, values.size + 1, dtype=float)


def write_csv(path: Path, rows: Iterable[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def plot_series_with_cumavg(
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


def main() -> int:
    p = argparse.ArgumentParser(
        description="CP2K DFTMD post-processing: center water potential (cube) + Fermi energy (md.out)."
    )
    p.add_argument(
        "--thickness",
        type=float,
        default=7.0,
        help="Center slab thickness in Angstrom for averaging the plane-averaged potential (default: 7).",
    )
    p.add_argument(
        "--cube-pattern",
        type=str,
        default="md-POTENTIAL-v_hartree-1_*.cube",
        help="Glob pattern for Hartree potential cube snapshots (default: md-POTENTIAL-v_hartree-1_*.cube).",
    )
    p.add_argument(
        "--md-out",
        type=str,
        default="md.out",
        help="CP2K output file containing 'Fermi energy:' and 'STEP NUMBER' (default: md.out).",
    )
    p.add_argument(
        "--outdir",
        type=str,
        default="analysis_outputs",
        help="Output directory for CSV/PNG files (default: analysis_outputs).",
    )
    p.add_argument(
        "--compute-u",
        action="store_true",
        default=True,
        help="Compute electrode potential U vs SHE at cube steps (needs both cube potentials and md.out Fermi energies).",
    )
    p.add_argument(
        "--no-compute-u",
        action="store_false",
        dest="compute_u",
        help="Disable U vs SHE computation.",
    )
    p.add_argument(
        "--fermi-unit",
        choices=["au", "ev"],
        default="au",
        help=(
            "Unit of 'Fermi energy' values in md.out. "
            "Use 'au' (default) for typical CP2K outputs (atomic units / Hartree); "
            "use 'ev' if your output is already in eV. "
            "A wrong choice typically shifts U by ~27.2×."
        ),
    )
    p.add_argument(
        "--center-mode",
        choices=["interface", "cell"],
        default="interface",
        help=(
            "How to define the center z for the water-slab averaging. "
            "'interface' uses midpoint between the two metal/water interfaces (needs xyz); "
            "if xyz is missing or interface detection fails, it will fall back to 'cell'. "
            "'cell' uses the geometric cell center (default: interface)."
        ),
    )
    p.add_argument(
        "--xyz",
        type=str,
        default="md-pos-1.xyz",
        help="XYZ trajectory file for interface detection (default: md-pos-1.xyz).",
    )
    p.add_argument(
        "--layer-tol",
        type=float,
        default=0.6,
        help="Metal layer clustering tolerance along z in Angstrom (default: 0.6).",
    )
    p.add_argument(
        "--metal-elements",
        type=str,
        default=None,
        help="Comma-separated metal element symbols for interface detection (e.g., 'Cu,Ag'). If omitted, infer as (all elements) - {H,O} from the first xyz frame.",
    )

    args = p.parse_args()

    workdir = Path(".").resolve()
    outdir = (workdir / args.outdir).resolve()

    # ---- Cube snapshots: center slab potential ----
    cube_paths = [Path(pth) for pth in sorted(workdir.glob(args.cube_pattern))]
    if not cube_paths:
        raise SystemExit(f"No cube files matched pattern: {args.cube_pattern!r} in {workdir}")

    # ---- Optional: interface centers from xyz ----
    use_interface_center = args.center_mode == "interface"
    metal_z_by_step: dict[int, np.ndarray] = {}
    metal_elements_used: Optional[set[str]] = None
    if use_interface_center:
        xyz_path = (workdir / args.xyz).resolve()
        if not xyz_path.exists():
            raise SystemExit(
                f"Interface center requested but xyz not found: {xyz_path}. "
                "Please provide --xyz and --metal-elements."
            )
        else:
            # only steps that we can map (encoded in the filename)
            needed_steps = {
                s for s in (_extract_step_from_cube_filename(cp) for cp in cube_paths) if s is not None
            }
            user_metal = _parse_csv_symbols(args.metal_elements)
            metal_z_by_step, inferred_metal = read_xyz_metal_z_for_steps(
                xyz_path,
                needed_steps,
                metal_elements=user_metal,
            )
            metal_elements_used = user_metal or inferred_metal
            if not metal_elements_used:
                raise SystemExit(
                    "Interface detection failed: no transition metals inferred from xyz. "
                    "Please specify the interface metal elements explicitly with --metal-elements (e.g., 'Cu,Ag')."
                )
            missing_steps = sorted(s for s in needed_steps if s not in metal_z_by_step)
            if missing_steps:
                raise SystemExit(
                    "Interface detection failed: missing metal coordinates for some cube steps. "
                    f"Missing steps (first 10): {missing_steps[:10]}. "
                    "Please verify xyz step indexing or specify --metal-elements."
                )

    cube_rows: list[dict] = []
    cube_steps: list[int] = []
    cube_vals_ev: list[float] = []
    slab_info_first: Optional[dict] = None
    center_rows: list[dict] = []

    for cp in cube_paths:
        step = _extract_step_from_cube_filename(cp)
        if step is None:
            # fallback to stable ordering index (1-based) if filename doesn't encode a step
            step = len(cube_steps)

        header, values = read_cube_header_and_values(cp)
        dz_bohr = float(np.linalg.norm(header.vz_bohr))
        lz_ang = header.nz * dz_bohr * BOHR_TO_ANG
        origin_z_ang = float(header.origin_bohr[2]) * BOHR_TO_ANG

        z_center_ang: Optional[float] = None
        iface: Optional[dict] = None
        center_source = "cell"
        if use_interface_center and int(step) in metal_z_by_step:
            try:
                iface = detect_metal_water_interfaces_z_ang(
                    metal_z_by_step[int(step)],
                    lz_ang,
                    layer_tol_ang=float(args.layer_tol),
                )
                # cube/xyz coordinates are assumed to share the same origin; include origin_z_ang for completeness
                z_center_ang = origin_z_ang + float(iface["z_mid_ang"])
                center_source = "interface"
            except Exception as e:
                inferred = ",".join(sorted(metal_elements_used or [])) or "(none)"
                raise SystemExit(
                    "Interface detection failed for step "
                    f"{step} ({cp.name}): {e}. "
                    f"Inferred metals: {inferred}. "
                    "Please specify interface metals explicitly with --metal-elements (e.g., 'Cu,Ag')."
                ) from e

        phi_ev, info = slab_average_potential_ev_from_cube(
            header, values, thickness_ang=float(args.thickness), z_center_ang=z_center_ang
        )
        if slab_info_first is None:
            slab_info_first = info
        cube_steps.append(int(step))
        cube_vals_ev.append(float(phi_ev))

        center_rows.append(
            {
                "step": int(step),
                "center_source": center_source,
                "z_center_ang": (None if z_center_ang is None else float(z_center_ang)),
                "z_iface_lower_ang": (None if iface is None else float(iface["z_lower_ang"] + origin_z_ang)),
                "z_iface_upper_ang": (None if iface is None else float(iface["z_upper_ang"] + origin_z_ang)),
                "z_iface_mid_ang": (None if iface is None else float(iface["z_mid_ang"] + origin_z_ang)),
                "water_gap_ang": (None if iface is None else float(iface["water_gap_ang"])),
                "n_metal_layers": (None if iface is None else int(iface["n_layers"])),
                "Lz_ang": float(lz_ang),
            }
        )

    # sort by step
    order = np.argsort(np.array(cube_steps, dtype=int))
    cube_steps_arr = np.array(cube_steps, dtype=int)[order]
    cube_vals_arr = np.array(cube_vals_ev, dtype=float)[order]
    cube_cum = cumulative_average(cube_vals_arr)

    for s, v, c in zip(cube_steps_arr, cube_vals_arr, cube_cum, strict=True):
        cube_rows.append(
            {
                "step": int(s),
                "phi_center_ev": float(v),
                "phi_center_cumavg_ev": float(c),
            }
        )

    write_csv(outdir / "center_potential.csv", cube_rows, ["step", "phi_center_ev", "phi_center_cumavg_ev"])
    # keep per-frame center/interface bookkeeping separate to avoid breaking existing CSV consumers
    if center_rows:
        # sort center rows by step
        center_rows_sorted = sorted(center_rows, key=lambda r: int(r["step"]))
        write_csv(
            outdir / "slab_center_and_interfaces.csv",
            center_rows_sorted,
            [
                "step",
                "center_source",
                "z_center_ang",
                "z_iface_lower_ang",
                "z_iface_upper_ang",
                "z_iface_mid_ang",
                "water_gap_ang",
                "n_metal_layers",
                "Lz_ang",
            ],
        )
    plot_series_with_cumavg(
        outdir / "center_potential.png",
        x=cube_steps_arr.astype(float),
        y=cube_vals_arr,
        y_cum=cube_cum,
        xlabel="MD step",
        ylabel="Center slab potential (eV)",
        title=f"Center water slab potential (thickness={args.thickness:g} A)",
    )

    phi_center_mean_ev = float(cube_vals_arr.mean())

    # ---- md.out: Fermi energy ----
    md_out_path = (workdir / args.md_out).resolve()
    if not md_out_path.exists():
        raise SystemExit(f"md.out not found: {md_out_path}")

    fermi_records = parse_md_out_fermi(md_out_path)
    if not fermi_records:
        raise SystemExit(f"No (step, Fermi energy) records parsed from: {md_out_path}")

    f_steps = np.array([r["step"] for r in fermi_records], dtype=int)
    f_time = np.array([r["time_fs"] if r["time_fs"] is not None else math.nan for r in fermi_records], dtype=float)
    f_raw = np.array([r["fermi_raw"] for r in fermi_records], dtype=float)
    if args.fermi_unit == "au":
        f_vals = f_raw * HA_TO_EV
    else:
        f_vals = f_raw
    f_cum = cumulative_average(f_vals)

    f_rows: list[dict] = []
    for s, t, raw, v, c in zip(f_steps, f_time, f_raw, f_vals, f_cum, strict=True):
        f_rows.append(
            {
                "step": int(s),
                "time_fs": (None if math.isnan(float(t)) else float(t)),
                "fermi_raw": float(raw),
                "fermi_ev": float(v),
                "fermi_cumavg_ev": float(c),
            }
        )

    write_csv(outdir / "fermi_energy.csv", f_rows, ["step", "time_fs", "fermi_raw", "fermi_ev", "fermi_cumavg_ev"])
    plot_series_with_cumavg(
        outdir / "fermi_energy.png",
        x=f_steps.astype(float),
        y=f_vals,
        y_cum=f_cum,
        xlabel="MD step",
        ylabel="Fermi energy (eV)",
        title="Fermi energy from md.out",
    )

    fermi_mean_ev = float(f_vals.mean())

    # ---- Optional: U vs SHE ----
    u_rows: list[dict] = []
    u_mean_v: Optional[float] = None
    if args.compute_u:
        # build step->fermi mapping
        fermi_by_step_ev = {int(s): float(v) for s, v in zip(f_steps, f_vals, strict=True)}

        u_vals: list[float] = []
        for row in cube_rows:
            s = int(row["step"])
            phi_ev = float(row["phi_center_ev"])
            if s not in fermi_by_step_ev:
                # Skip if we cannot align this cube snapshot to a MD step in md.out
                continue
            ef_ev = fermi_by_step_ev[s]
            # IMPORTANT (sign convention):
            # For CP2K V_HARTREE_CUBE, the reported "Hartree potential" is an energy-like
            # potential entering the Kohn–Sham equations (i.e., it corresponds to -e*phi
            # rather than +e*phi). Therefore the cSHE formula term "-e*phi_wat" is obtained
            # by *adding* the averaged cube value.
            e_u_ev = -ef_ev + phi_ev + DP_A_H3O_W_EV - MU_HPLUS_G0_EV - DELTA_E_ZP_EV
            u_v = e_u_ev  # numerically identical (eU in eV -> U in V)
            u_vals.append(u_v)
            u_rows.append({"step": s, "U_vs_SHE_V": float(u_v)})

        if u_rows:
            u_vals_arr = np.array(u_vals, dtype=float)
            u_cum = cumulative_average(u_vals_arr)
            for i, r in enumerate(u_rows):
                r["U_cumavg_V"] = float(u_cum[i])
            write_csv(outdir / "electrode_potential_U_vs_SHE.csv", u_rows, ["step", "U_vs_SHE_V", "U_cumavg_V"])
            plot_series_with_cumavg(
                outdir / "electrode_potential_U_vs_SHE.png",
                x=np.array([r["step"] for r in u_rows], dtype=float),
                y=np.array([r["U_vs_SHE_V"] for r in u_rows], dtype=float),
                y_cum=np.array([r["U_cumavg_V"] for r in u_rows], dtype=float),
                xlabel="MD step",
                ylabel="U (V vs SHE)",
                title="Electrode potential U vs SHE (aligned at cube steps)",
            )
            u_mean_v = float(u_vals_arr.mean())

    # ---- Terminal summary ----
    print("=== CP2K post-processing summary ===")
    if slab_info_first:
        print(
            f"Cube grid: nx,ny,nz={slab_info_first['nx']},{slab_info_first['ny']},{slab_info_first['nz']} | "
            f"dz={slab_info_first['dz_ang']:.4f} A | slab slices={slab_info_first['n_slices']} | thickness={args.thickness:g} A"
        )
    if args.center_mode == "interface":
        me = ",".join(sorted(metal_elements_used)) if metal_elements_used else "(none)"
        print(f"Slab center:            interface-midpoint (layer_tol={args.layer_tol:g} A, metal={me})")
    else:
        print("Slab center:            cell center")
    print(f"Center slab potential: mean = {phi_center_mean_ev:.6f} eV  (from {len(cube_rows)} cube snapshots)")
    print(
        f"Fermi energy:          mean = {fermi_mean_ev:.6f} eV  (from {len(f_vals)} MD steps, parsed as {args.fermi_unit})"
    )
    if args.compute_u:
        if u_mean_v is None:
            print("U vs SHE:              could not compute (no cube steps matched md.out steps)")
        else:
            print(f"U vs SHE:              mean = {u_mean_v:.6f} V  (aligned at {len(u_rows)} cube steps)")
            # Very large |U| is often a sign of a wrong unit assumption for the Fermi energy.
            if abs(float(u_mean_v)) > 2.0:
                other = "ev" if args.fermi_unit == "au" else "au"
                print(
                    "WARNING: |U| > 2 V (unusual for metal PZC vs SHE). "
                    f"Please double-check --fermi-unit (currently {args.fermi_unit}); "
                    f"you may want to try --fermi-unit {other}."
                )
    print(f"Outputs written to: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

