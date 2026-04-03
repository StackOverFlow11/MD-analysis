"""Potential frame data source abstraction.

Provides a unified interface for discovering and loading potential
analysis frames from two input modes:

- **Continuous MD** (mode A): cube files + md.out in a single directory.
- **Distributed single-point** (mode B): ``potential_t{time}_i{step}/``
  subdirectories, each containing one cube file and one sp.out.

Public API
----------
- ``PotentialFrame`` — frozen dataclass holding one frame's data
- ``discover_continuous_frames`` — mode A discovery
- ``discover_distributed_frames`` — mode B discovery
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from ...utils.CubeParser import (
    CubeHeader,
    _float,
    discover_cube_files,
    extract_step_from_cube_filename,
    read_cube_atoms,
    read_cube_header_and_values,
)
from ...utils.config import (
    BOHR_TO_ANG,
    HA_TO_EV,
    TRANSITION_METAL_SYMBOLS,
)

try:
    from ase import Atoms
except ImportError:  # pragma: no cover
    Atoms = object  # type: ignore[misc]

logger = logging.getLogger(__name__)

# Regex shared with CenterPotential for md.out parsing
FERMI_RE = re.compile(r"Fermi energy:\s*([+-]?\d+(?:\.\d*)?(?:[EeDd][+-]?\d+)?)")
STEP_RE = re.compile(r"STEP NUMBER\s*=\s*(\d+)")
TIME_RE = re.compile(r"TIME\s*\[fs\]\s*=\s*([+-]?\d+(?:\.\d*)?(?:[EeDd][+-]?\d+)?)")
XYZ_STEP_RE = re.compile(r"\bi\s*=\s*(\d+)\b")

# Directory name pattern: potential_t{time}_i{step}
_SP_DIR_RE = re.compile(r"_t(\d+)_i(\d+)")


@dataclass(frozen=True)
class PotentialFrame:
    """One frame of potential analysis data (immutable)."""

    step: int
    time_fs: float | None
    cube_path: Path
    header: CubeHeader
    values: np.ndarray
    fermi_raw: float | None  # Hartree; None if unavailable
    atoms: Atoms | None  # for interface detection (with cell set)


# ---------------------------------------------------------------------------
# Mode A: Continuous MD
# ---------------------------------------------------------------------------

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


def _read_xyz_atoms_for_steps(
    xyz_path: Path,
    steps: set[int],
    *,
    metal_elements: Optional[set[str]] = None,
) -> tuple[dict[int, Atoms], set[str]]:
    """Stream-parse a CP2K xyz trajectory and build Atoms for given steps."""
    if not xyz_path.exists():
        raise FileNotFoundError(xyz_path)

    steps = {int(s) for s in steps}
    out: dict[int, Atoms] = {}
    inferred_metal: Optional[set[str]] = (
        set(metal_elements) if metal_elements is not None else None
    )
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

            need_parse = need_this or not inferred_done
            if need_parse:
                symbols: list[str] = []
                positions: list[list[float]] = []
                for _ in range(natoms):
                    line = f.readline()
                    if not line:
                        break
                    parts = line.split()
                    if len(parts) < 4:
                        continue
                    symbols.append(parts[0])
                    if need_this:
                        positions.append(
                            [_float(parts[1]), _float(parts[2]), _float(parts[3])]
                        )

                if not inferred_done:
                    inferred_metal = set(symbols) & set(TRANSITION_METAL_SYMBOLS)
                    inferred_done = True

                if need_this and step is not None:
                    out[int(step)] = Atoms(symbols=symbols, positions=positions)
            else:
                for _ in range(natoms):
                    if not f.readline():
                        break

    if inferred_metal is None:
        inferred_metal = set()
    return out, inferred_metal


def discover_continuous_frames(
    cube_pattern: str,
    *,
    workdir: Path | None = None,
    md_out_path: Path | None = None,
    xyz_path: Path | None = None,
    center_mode: str = "interface",
    metal_elements: set[str] | None = None,
    fermi_unit: str = "au",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
) -> list[PotentialFrame]:
    """Discover frames from a continuous MD run (mode A).

    Parameters
    ----------
    cube_pattern
        Glob pattern for cube files (relative to *workdir*).
    workdir
        Base directory. Defaults to cwd.
    md_out_path
        Path to CP2K md.out for Fermi energy extraction.
        If ``None``, frames will have ``fermi_raw=None``.
    xyz_path
        Path to XYZ trajectory for interface detection.
        Required when ``center_mode="interface"``.
    center_mode
        ``"interface"`` or ``"cell"``.
    metal_elements
        Explicit metal element set; auto-detected if ``None``.
    fermi_unit
        ``"au"`` (Hartree) or ``"ev"`` for raw Fermi values in md.out.

    Returns
    -------
    list[PotentialFrame]
        Sorted by step.
    """
    workdir = (workdir or Path(".")).resolve()

    cube_paths = discover_cube_files(
        cube_pattern,
        workdir=workdir,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
    )

    # Parse Fermi energies
    fermi_by_step: dict[int, float] = {}
    if md_out_path is not None:
        md_out_path = Path(md_out_path).resolve()
        if md_out_path.exists():
            records = _parse_md_out_fermi(md_out_path)
            records = records[frame_start:frame_end:frame_step]
            for r in records:
                fermi_by_step[int(r["step"])] = float(r["fermi_raw"])

    # Parse xyz for interface detection
    use_interface = center_mode == "interface"
    atoms_by_step: dict[int, Atoms] = {}
    metal_used: set[str] | None = None

    if use_interface:
        if xyz_path is None:
            xyz_path = workdir / "md-pos-1.xyz"
        xyz_path = Path(xyz_path).resolve()
        if xyz_path.exists():
            needed_steps = {
                s
                for s in (extract_step_from_cube_filename(cp) for cp in cube_paths)
                if s is not None
            }
            atoms_by_step, inferred = _read_xyz_atoms_for_steps(
                xyz_path, needed_steps, metal_elements=metal_elements
            )
            metal_used = metal_elements or inferred

    # Build frames
    frames: list[PotentialFrame] = []
    for cp in cube_paths:
        step = extract_step_from_cube_filename(cp)
        if step is None:
            step = len(frames)

        header, values = read_cube_header_and_values(cp)

        # Fermi energy
        fermi_raw = fermi_by_step.get(int(step))

        # Atoms with cell from cube header
        frame_atoms: Atoms | None = None
        if use_interface and int(step) in atoms_by_step:
            frame_atoms = atoms_by_step[int(step)]
            frame_atoms.set_cell([
                header.vx_bohr * header.nx * BOHR_TO_ANG,
                header.vy_bohr * header.ny * BOHR_TO_ANG,
                header.vz_bohr * header.nz * BOHR_TO_ANG,
            ])
            frame_atoms.set_pbc(True)

        frames.append(PotentialFrame(
            step=int(step),
            time_fs=None,
            cube_path=cp,
            header=header,
            values=values,
            fermi_raw=fermi_raw,
            atoms=frame_atoms,
        ))

    frames.sort(key=lambda f: f.step)
    return frames


# ---------------------------------------------------------------------------
# Mode B: Distributed single-point
# ---------------------------------------------------------------------------

def _parse_sp_out_fermi(sp_out_path: Path) -> float | None:
    """Extract the last Fermi energy (Hartree) from a single-point sp.out.

    Returns ``None`` if no ``Fermi energy:`` line is found.
    """
    fermi_raw: float | None = None
    with sp_out_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            m = FERMI_RE.search(line)
            if m:
                fermi_raw = _float(m.group(1))
    return fermi_raw


def discover_distributed_frames(
    root_dir: Path | str,
    *,
    dir_pattern: str = "potential_t*_i*",
    cube_filename: str = "sp_potential-v_hartree-1_0.cube",
    sp_out_filename: str = "sp.out",
    center_mode: str = "interface",
    metal_elements: set[str] | None = None,
    layer_tol_ang: float = 0.6,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> list[PotentialFrame]:
    """Discover frames from distributed single-point directories (mode B).

    Each subdirectory matching *dir_pattern* under *root_dir* should contain
    a cube file (*cube_filename*) and optionally an output file
    (*sp_out_filename*) with Fermi energy.

    Parameters
    ----------
    root_dir
        Parent directory containing ``potential_t*_i*`` subdirectories.
    dir_pattern
        Glob pattern for subdirectory names.
    cube_filename
        Name of the cube file inside each subdirectory.
    sp_out_filename
        Name of the CP2K output file (for Fermi energy).
    center_mode
        ``"interface"`` or ``"cell"``.
    metal_elements
        Explicit metal element set; auto-detected from first frame if ``None``.
    layer_tol_ang
        Layer clustering tolerance for interface detection.

    Returns
    -------
    list[PotentialFrame]
        Sorted by step.
    """
    root = Path(root_dir).resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Root directory not found: {root}")

    # Discover and sort subdirectories
    subdirs: list[tuple[Path, int, int]] = []  # (path, time_fs, step)
    for d in sorted(root.glob(dir_pattern)):
        if not d.is_dir():
            continue
        m = _SP_DIR_RE.search(d.name)
        if m is None:
            logger.warning("Cannot parse step/time from directory name: %s", d.name)
            continue
        time_fs = int(m.group(1))
        step = int(m.group(2))
        cube_path = d / cube_filename
        if not cube_path.exists():
            logger.debug("Cube file missing in %s, skipping", d.name)
            continue
        subdirs.append((d, time_fs, step))

    # Sort by step
    subdirs.sort(key=lambda x: x[2])

    if not subdirs:
        raise FileNotFoundError(
            f"No valid single-point directories found matching {dir_pattern!r} "
            f"in {root} (with {cube_filename})"
        )

    # Apply frame slice
    subdirs = subdirs[frame_start:frame_end:frame_step]
    logger.info(
        "Distributed SP: %d directories in %s", len(subdirs), root
    )

    # Build frames
    use_interface = center_mode == "interface"
    metal_detected: set[str] | None = (
        set(metal_elements) if metal_elements else None
    )
    frames: list[PotentialFrame] = []

    sp_iter = subdirs
    if verbose:
        from tqdm import tqdm

        sp_iter = tqdm(subdirs, desc="Loading SP frames", unit="dir", ascii=" =")

    for d, time_fs, step in sp_iter:
        cube_path = d / cube_filename
        header, values = read_cube_header_and_values(cube_path)

        # Fermi energy from sp.out
        sp_out = d / sp_out_filename
        fermi_raw: float | None = None
        if sp_out.exists():
            fermi_raw = _parse_sp_out_fermi(sp_out)

        # Atoms from cube file for interface detection
        frame_atoms: Atoms | None = None
        if use_interface:
            frame_atoms = read_cube_atoms(cube_path, header)
            if metal_detected is None:
                metal_detected = (
                    set(frame_atoms.get_chemical_symbols())
                    & set(TRANSITION_METAL_SYMBOLS)
                )
                if metal_detected:
                    logger.info(
                        "Auto-detected metal elements: %s", metal_detected
                    )

        frames.append(PotentialFrame(
            step=step,
            time_fs=float(time_fs),
            cube_path=cube_path,
            header=header,
            values=values,
            fermi_raw=fermi_raw,
            atoms=frame_atoms,
        ))

    return frames
