#!/share/home/chem-wangyg/shaofl/local/src/installs/miniconda/envs/md_env/bin/python3
"""
Extract a specified frame (by MD step) from a CP2K XYZ trajectory and write it to CIF.

Inputs (default: current working directory):
  - CP2K XYZ trajectory (*-pos-1.xyz by default) for atomic positions
  - CP2K input file (*.inp, auto-detected from the trajectory directory) for the simulation cell

CLI arguments:
  -traj   (optional): XYZ trajectory path. If omitted, auto-detect one file in CWD that ends with "-pos-1.xyz".
  -opath  (optional): output CIF path. If omitted, write to "./traj_<fs>.cif" in the current working directory.
  -fs     (required): requested MD step number (i=... in the XYZ comment line). If not divisible by 5, use the nearest step.

Notes:
  - This script expects XYZ comment lines like: "i =  10, time =  5.000, E = ..."
    (standard CP2K format). If "i=" is missing, the script cannot locate the step reliably.
  - The output CIF uses fractional coordinates wrapped into [0, 1).
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


RE_TIME = re.compile(r"time\s*=\s*([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)", re.I)
RE_STEP = re.compile(r"\bi\s*=\s*([+-]?\d+)\b", re.I)
RE_ENERGY = re.compile(r"\bE\s*=\s*([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)", re.I)


@dataclass(frozen=True)
class Frame:
    n: int
    comment: str
    symbols: List[str]
    coords: np.ndarray  # (n,3) float, Angstrom
    step: Optional[int] = None
    time_fs: Optional[float] = None
    energy: Optional[float] = None


@dataclass(frozen=True)
class FrameMeta:
    # Byte offset (in binary mode) of the beginning of the frame (the atom-count line)
    offset: int
    step: Optional[int]
    time_fs: Optional[float]
    comment: str


def parse_xyz_comment(comment: str) -> Tuple[Optional[int], Optional[float], Optional[float]]:
    step = None
    time_fs = None
    energy = None

    m = RE_STEP.search(comment)
    if m:
        try:
            step = int(m.group(1))
        except ValueError:
            step = None
    m = RE_TIME.search(comment)
    if m:
        try:
            time_fs = float(m.group(1))
        except ValueError:
            time_fs = None
    m = RE_ENERGY.search(comment)
    if m:
        try:
            energy = float(m.group(1))
        except ValueError:
            energy = None
    return step, time_fs, energy


def _readline_text_from_binary(f) -> Optional[str]:
    line = f.readline()
    if not line:
        return None
    # Keep '\r' out if present
    return line.decode("utf-8", errors="ignore").rstrip("\n").rstrip("\r")


def build_xyz_frame_index(xyz_path: Path) -> List[FrameMeta]:
    """
    Build an index of frame byte offsets and meta (step/time/comment).
    """
    metas: List[FrameMeta] = []
    with xyz_path.open("rb") as f:
        while True:
            offset = f.tell()
            line = _readline_text_from_binary(f)
            if line is None:
                break
            if not line.strip():
                continue

            try:
                n = int(line.strip())
            except ValueError:
                raise ValueError(f"Invalid XYZ atom count line at offset {offset}: {line!r}")

            comment = _readline_text_from_binary(f)
            if comment is None:
                raise ValueError("Unexpected EOF while reading XYZ comment line.")
            step, time_fs, _ = parse_xyz_comment(comment)

            # Skip atom coordinate lines
            for _ in range(n):
                atom_line = f.readline()
                if not atom_line:
                    raise ValueError("Unexpected EOF while skipping XYZ atoms.")

            metas.append(FrameMeta(offset=offset, step=step, time_fs=time_fs, comment=comment))

    return metas


def read_xyz_frame_at_offset(xyz_path: Path, offset: int) -> Frame:
    """
    Random-access read of a single XYZ frame.
    """
    with xyz_path.open("rb") as f:
        f.seek(offset)
        line = _readline_text_from_binary(f)
        if line is None:
            raise ValueError(f"Offset out of range: {offset}")
        while line is not None and not line.strip():
            offset = f.tell()
            line = _readline_text_from_binary(f)
        if line is None:
            raise ValueError(f"Offset out of range: {offset}")

        try:
            n = int(line.strip())
        except ValueError:
            raise ValueError(f"Invalid XYZ atom count line at offset {offset}: {line!r}")

        comment = _readline_text_from_binary(f)
        if comment is None:
            raise ValueError("Unexpected EOF while reading XYZ comment line.")
        step, time_fs, energy = parse_xyz_comment(comment)

        symbols: List[str] = []
        coords = np.zeros((n, 3), dtype=float)
        for i in range(n):
            atom_line = _readline_text_from_binary(f)
            if atom_line is None:
                raise ValueError("Unexpected EOF while reading XYZ atoms.")
            parts = atom_line.split()
            if len(parts) < 4:
                raise ValueError(f"Invalid XYZ atom line: {atom_line!r}")
            symbols.append(parts[0].strip())
            coords[i, 0] = float(parts[1])
            coords[i, 1] = float(parts[2])
            coords[i, 2] = float(parts[3])

        return Frame(n=n, comment=comment, symbols=symbols, coords=coords, step=step, time_fs=time_fs, energy=energy)


def autodetect_one(path: Path, patterns: Sequence[str]) -> Optional[Path]:
    matches: List[Path] = []
    for pat in patterns:
        matches.extend(path.glob(pat))
    matches = sorted(set(matches))
    if len(matches) == 1:
        return matches[0]
    return None


def _strip_cp2k_comments(line: str) -> str:
    # CP2K supports '!' and '#'. We keep it simple: drop everything after the first marker.
    # This is good enough for CELL parsing.
    for marker in ("!", "#"):
        if marker in line:
            line = line.split(marker, 1)[0]
    return line.strip()


def _parse_float_tokens(tokens: Sequence[str]) -> List[float]:
    # Drop unit tokens like [angstrom]
    vals: List[float] = []
    for t in tokens:
        if t.startswith("[") and t.endswith("]"):
            continue
        try:
            vals.append(float(t))
        except ValueError:
            # ignore non-float tokens
            pass
    return vals


def parse_cp2k_cell(inp_path: Path) -> np.ndarray:
    """
    Return cell matrix M (3x3) with column vectors A, B, C in Angstrom.
    Supports &CELL with either (A,B,C vectors) or (ABC lengths + angles).
    """
    in_cell = False
    a_vec = b_vec = c_vec = None
    abc = None
    angles = None

    with inp_path.open("r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = _strip_cp2k_comments(raw)
            if not line:
                continue

            u = line.upper()
            if u.startswith("&CELL"):
                in_cell = True
                continue
            if in_cell and u.startswith("&END") and "CELL" in u:
                in_cell = False
                continue

            if not in_cell:
                continue

            parts = line.split()
            if not parts:
                continue
            key = parts[0].upper()
            rest = parts[1:]

            if key in {"A", "B", "C"}:
                vals = _parse_float_tokens(rest)
                if len(vals) >= 3:
                    vec = np.array(vals[:3], dtype=float)
                    if key == "A":
                        a_vec = vec
                    elif key == "B":
                        b_vec = vec
                    else:
                        c_vec = vec
            elif key == "ABC":
                vals = _parse_float_tokens(rest)
                if len(vals) >= 3:
                    abc = tuple(float(x) for x in vals[:3])
            elif key == "ALPHA_BETA_GAMMA":
                vals = _parse_float_tokens(rest)
                if len(vals) >= 3:
                    angles = tuple(float(x) for x in vals[:3])

    if a_vec is not None and b_vec is not None and c_vec is not None:
        return np.column_stack([a_vec, b_vec, c_vec])

    if abc is None:
        raise ValueError(
            f"Cannot find cell in {inp_path}. Expected &CELL with either A/B/C vectors or ABC lengths."
        )

    a, b, c = abc
    if angles is None:
        # Default orthorhombic if angles are not provided
        alpha = beta = gamma = 90.0
    else:
        alpha, beta, gamma = angles

    # Build a right-handed cell with A along x and B in xy plane
    deg = math.pi / 180.0
    ca, cb, cg = math.cos(alpha * deg), math.cos(beta * deg), math.cos(gamma * deg)
    sg = math.sin(gamma * deg)
    if abs(sg) < 1e-12:
        raise ValueError(f"Invalid CELL angles in {inp_path}: gamma too close to 0/180.")

    a_vec = np.array([a, 0.0, 0.0], dtype=float)
    b_vec = np.array([b * cg, b * sg, 0.0], dtype=float)
    c_x = c * cb
    c_y = c * (ca - cb * cg) / sg
    c_z2 = c * c - c_x * c_x - c_y * c_y
    c_z = math.sqrt(max(c_z2, 0.0))
    c_vec = np.array([c_x, c_y, c_z], dtype=float)
    return np.column_stack([a_vec, b_vec, c_vec])


def cell_to_lengths_angles(M: np.ndarray) -> Tuple[float, float, float, float, float, float]:
    """
    Convert cell matrix (columns are A,B,C) to (a,b,c, alpha, beta, gamma).
    Angles are in degrees.
    """
    A = np.asarray(M[:, 0], dtype=float)
    B = np.asarray(M[:, 1], dtype=float)
    C = np.asarray(M[:, 2], dtype=float)

    a = float(np.linalg.norm(A))
    b = float(np.linalg.norm(B))
    c = float(np.linalg.norm(C))

    def angle(u: np.ndarray, v: np.ndarray) -> float:
        nu = float(np.linalg.norm(u))
        nv = float(np.linalg.norm(v))
        if nu <= 1e-15 or nv <= 1e-15:
            return 90.0
        cosx = float(np.dot(u, v) / (nu * nv))
        if cosx > 1.0:
            cosx = 1.0
        elif cosx < -1.0:
            cosx = -1.0
        return float(np.degrees(np.arccos(cosx)))

    alpha = angle(B, C)
    beta = angle(A, C)
    gamma = angle(A, B)
    return a, b, c, alpha, beta, gamma


def snap_step_to_nearest_5(req_step: int) -> int:
    # "If not divisible by 5, pick the nearest frame" (CP2K traj written every 5 steps).
    return int(round(float(req_step) / 5.0)) * 5


def choose_meta_by_step(metas: Sequence[FrameMeta], target_step: int) -> FrameMeta:
    metas_step = [m for m in metas if m.step is not None]
    if not metas_step:
        raise ValueError("Cannot find 'i=' (MD step) in XYZ comment lines; cannot select by -fs.")

    # Prefer exact match
    for m in metas_step:
        if int(m.step) == int(target_step):
            return m

    # Else choose nearest; tie-breaker: smaller step
    return min(metas_step, key=lambda m: (abs(int(m.step) - int(target_step)), int(m.step)))


def write_cif(
    out_path: Path,
    M: np.ndarray,
    symbols: Sequence[str],
    coords: np.ndarray,
    data_name: str,
) -> None:
    invM = np.linalg.inv(M)
    frac = coords @ invM.T
    frac = frac - np.floor(frac)  # wrap into [0,1)

    a, b, c, alpha, beta, gamma = cell_to_lengths_angles(M)

    # Build per-element labels like Cu1, Cu2, O1, H1...
    counts: Dict[str, int] = defaultdict(int)
    labels: List[str] = []
    for sym in symbols:
        s = sym.strip()
        counts[s] += 1
        labels.append(f"{s}{counts[s]}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="\n") as f:
        f.write(f"data_{data_name}\n")
        f.write(f"_audit_creation_date              {date.today().isoformat()}\n")
        f.write(f"_audit_creation_method            'traj_frame_to_cif.py'\n")
        f.write(f"_symmetry_space_group_name_H-M    'P1'\n")
        f.write(f"_symmetry_Int_Tables_number       1\n")
        f.write(f"_symmetry_cell_setting            triclinic\n")
        f.write("loop_\n")
        f.write("_symmetry_equiv_pos_as_xyz\n")
        f.write("  x,y,z\n")
        f.write(f"_cell_length_a                    {a:.4f}\n")
        f.write(f"_cell_length_b                    {b:.4f}\n")
        f.write(f"_cell_length_c                    {c:.4f}\n")
        f.write(f"_cell_angle_alpha                 {alpha:.4f}\n")
        f.write(f"_cell_angle_beta                  {beta:.4f}\n")
        f.write(f"_cell_angle_gamma                 {gamma:.4f}\n")

        f.write("loop_\n")
        f.write("_atom_site_label\n")
        f.write("_atom_site_type_symbol\n")
        f.write("_atom_site_fract_x\n")
        f.write("_atom_site_fract_y\n")
        f.write("_atom_site_fract_z\n")
        f.write("_atom_site_U_iso_or_equiv\n")
        f.write("_atom_site_adp_type\n")
        f.write("_atom_site_occupancy\n")

        for label, sym, (fx, fy, fz) in zip(labels, symbols, frac.tolist()):
            # Uiso is not known; keep 0.00000 (good enough for geometry / visualization).
            f.write(f"{label:<6s} {sym:<3s} {fx: .5f} {fy: .5f} {fz: .5f}  0.00000  Uiso   1.00\n")


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(
        prog="traj_frame_to_cif.py",
        description="Extract a specified MD step frame from a CP2K XYZ trajectory and write to CIF.",
    )
    p.add_argument("-traj", type=str, default="", help='Trajectory XYZ path (default: auto-detect "*-pos-1.xyz" in CWD).')
    p.add_argument("-opath", type=str, default="", help="Output CIF path (default: ./traj_<fs>.cif).")
    p.add_argument("-fs", type=int, required=True, help="Requested MD step (i=...) to extract. If not divisible by 5, use nearest.")

    args = p.parse_args(argv)

    cwd = Path(".").resolve()
    if args.traj:
        traj_path = Path(args.traj).expanduser().resolve()
    else:
        traj_path = autodetect_one(cwd, patterns=["*-pos-1.xyz"])
        if traj_path is None:
            print('[ERROR] Cannot auto-detect trajectory in current directory. Expected exactly one "*-pos-1.xyz".', file=sys.stderr)
            return 2

    if not traj_path.exists():
        print(f"[ERROR] Trajectory not found: {traj_path}", file=sys.stderr)
        return 2

    req_step = int(args.fs)
    if req_step < 0:
        print("[ERROR] -fs must be a non-negative integer MD step.", file=sys.stderr)
        return 2

    target_step = snap_step_to_nearest_5(req_step)
    if target_step != req_step:
        print(f"[WARN] Requested fs={req_step} is not divisible by 5. Using nearest fs={target_step}.")

    metas = build_xyz_frame_index(traj_path)
    if not metas:
        print("[ERROR] No frames found in trajectory.", file=sys.stderr)
        return 2

    try:
        meta = choose_meta_by_step(metas, target_step=target_step)
    except Exception as e:
        print(f"[ERROR] {type(e).__name__}: {e}", file=sys.stderr)
        return 2

    chosen_step = meta.step
    # Choose output path (in current working directory by default)
    if args.opath:
        out_path = Path(args.opath).expanduser().resolve()
    else:
        # Use the actually selected step if available, otherwise fall back to the requested target_step.
        step_tag = target_step if chosen_step is None else int(chosen_step)
        out_path = (cwd / f"traj_{step_tag}.cif").resolve()

    # Parse cell from a CP2K input file in the trajectory directory (same as water_profile.py style)
    traj_dir = traj_path.parent
    inp_path = autodetect_one(traj_dir, patterns=["*.inp"])
    if inp_path is None or not inp_path.exists():
        print(f"[ERROR] Cannot auto-detect a unique *.inp in trajectory directory: {traj_dir}", file=sys.stderr)
        print("[ERROR] CIF requires cell parameters. Please ensure the directory contains exactly one CP2K input (*.inp).", file=sys.stderr)
        return 2

    try:
        M = parse_cp2k_cell(inp_path)
    except Exception as e:
        print(f"[ERROR] Failed to parse cell from input: {inp_path}", file=sys.stderr)
        print(f"[ERROR] {type(e).__name__}: {e}", file=sys.stderr)
        return 2

    fr = read_xyz_frame_at_offset(traj_path, meta.offset)
    if fr.step is None:
        print("[ERROR] The selected frame has no 'i=' step in the comment line; cannot report selected step.", file=sys.stderr)
    else:
        if int(fr.step) != int(target_step) and req_step == target_step:
            # Rare case: exact multiple-of-5 requested but not found; meta may be nearest.
            print(f"[WARN] Exact step {target_step} not found. Using nearest available step {int(fr.step)}.")

    # Build a simple data_ name (CIF identifier): avoid spaces/special chars
    data_name = f"traj_{int(fr.step) if fr.step is not None else target_step}"

    write_cif(out_path, M=M, symbols=fr.symbols, coords=fr.coords, data_name=data_name)

    print(f"[INFO] Trajectory : {traj_path}")
    print(f"[INFO] Input cell : {inp_path}")
    print(f"[INFO] Requested : fs={req_step} (snapped={target_step})")
    if fr.step is not None:
        print(f"[INFO] Selected  : step={int(fr.step)}")
    if fr.time_fs is not None:
        print(f"[INFO] Time      : {fr.time_fs:.6f} fs")
    print(f"[INFO] Output CIF : {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


