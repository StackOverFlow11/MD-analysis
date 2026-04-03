"""Gaussian cube file I/O and single-frame potential utilities.

Provides reading of CP2K Hartree-potential cube files, plane-averaged
potential extraction φ(z), and slab-averaged potential computation.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from .config import BOHR_TO_ANG, HA_TO_EV

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class CubeHeader:
    """Gaussian cube file header information."""

    natoms: int
    origin_bohr: np.ndarray  # (3,)
    nx: int
    ny: int
    nz: int
    vx_bohr: np.ndarray  # (3,) voxel step vector
    vy_bohr: np.ndarray  # (3,)
    vz_bohr: np.ndarray  # (3,)


def _float(s: str) -> float:
    """Parse a float, tolerating Fortran 'D' exponents."""
    return float(s.replace("D", "E").replace("d", "e"))


def read_cube_header_and_values(path: Path) -> tuple[CubeHeader, np.ndarray]:
    """Read a Gaussian cube file and return ``(header, flat_values)``.

    The flat values array has length ``nx * ny * nz``.  CP2K writes
    V_HARTREE_CUBE with z as the fastest-running index, so the array
    should be reshaped as ``(nx, ny, nz)`` in C-order.
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

        for _ in range(natoms):
            f.readline()

        rest = f.read()
        if "D" in rest or "d" in rest:
            rest = rest.replace("D", "E").replace("d", "e")
        values = np.fromstring(rest, sep=" ", dtype=float)

    expected = nx * ny * nz
    if values.size != expected:
        raise ValueError(
            f"Cube data size mismatch in {path}: got {values.size}, "
            f"expected {expected} (nx,ny,nz={nx},{ny},{nz})"
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


def slab_average_potential_ev(
    header: CubeHeader,
    values: np.ndarray,
    thickness_ang: float,
    *,
    z_center_ang: Optional[float] = None,
) -> tuple[float, dict]:
    """Compute the slab-averaged Hartree potential (eV) from a cube snapshot.

    The cube is assumed to be aligned such that the 3rd axis is the
    surface normal (z).

    Parameters
    ----------
    header
        Cube file header.
    values
        Flat scalar-field array of length ``nx * ny * nz``.
    thickness_ang
        Slab averaging thickness in Ångström.
    z_center_ang
        Absolute z-center for the slab (Å).  If ``None``, uses the
        geometric cell center.

    Returns
    -------
    (phi_center_ev, info_dict)
    """
    dz_bohr = float(np.linalg.norm(header.vz_bohr))
    if dz_bohr <= 0:
        raise ValueError(f"Invalid dz in cube header: |vz|={dz_bohr} bohr")

    field = values.reshape((header.nx, header.ny, header.nz))
    phi_z_ha = field.mean(axis=(0, 1))  # (nz,)

    thickness_bohr = thickness_ang / BOHR_TO_ANG
    half_thickness_bohr = 0.5 * thickness_bohr

    z0_bohr = float(header.origin_bohr[2])
    lz_bohr = header.nz * dz_bohr
    z_coords_rel = (np.arange(header.nz, dtype=float) + 0.5) * dz_bohr

    if z_center_ang is None:
        z_center_rel = 0.5 * lz_bohr
        z_center_source = "cell"
    else:
        z_center_bohr = float(z_center_ang) / BOHR_TO_ANG
        z_center_rel = float((z_center_bohr - z0_bohr) % lz_bohr)
        z_center_source = "interface"

    dist = np.abs(z_coords_rel - z_center_rel)
    dist = np.minimum(dist, lz_bohr - dist)
    mask = dist <= half_thickness_bohr
    if not np.any(mask):
        raise ValueError(
            f"Thickness {thickness_ang} A selects 0 z-slices "
            f"(dz={dz_bohr * BOHR_TO_ANG:.4f} A, nz={header.nz})"
        )

    phi_center_ha = float(phi_z_ha[mask].mean())
    phi_center_ev = phi_center_ha * HA_TO_EV
    phi_std_ev = float(phi_z_ha[mask].std()) * HA_TO_EV

    info = {
        "nx": header.nx,
        "ny": header.ny,
        "nz": header.nz,
        "dz_ang": dz_bohr * BOHR_TO_ANG,
        "n_slices": int(mask.sum()),
        "z_center_source": z_center_source,
        "z_center_rel_bohr": z_center_rel,
        "thickness_ang": thickness_ang,
        "phi_z_std_ev": phi_std_ev,
    }
    return phi_center_ev, info


def plane_avg_phi_z_ev(header: CubeHeader, values: np.ndarray) -> np.ndarray:
    """Return the xy plane-averaged potential φ(z) in eV, shape ``(nz,)``."""
    field = values.reshape((header.nx, header.ny, header.nz))
    phi_z_ha = field.mean(axis=(0, 1))
    return phi_z_ha * HA_TO_EV


def z_coords_ang(header: CubeHeader) -> np.ndarray:
    """Return the z-coordinate array (Å) from a cube header, shape ``(nz,)``."""
    dz_bohr = float(np.linalg.norm(header.vz_bohr))
    dz_ang = dz_bohr * BOHR_TO_ANG
    origin_z_ang = float(header.origin_bohr[2]) * BOHR_TO_ANG
    return origin_z_ang + (np.arange(header.nz, dtype=float) + 0.5) * dz_ang


def discover_cube_files(
    cube_pattern: str,
    *,
    workdir: Path | None = None,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
) -> list[Path]:
    """Discover and slice cube files matching *cube_pattern*.

    Parameters
    ----------
    cube_pattern
        Glob pattern relative to *workdir*.
    workdir
        Base directory. Defaults to ``Path(".").resolve()``.
    frame_start, frame_end, frame_step
        Optional slice parameters applied to the sorted file list.

    Returns
    -------
    list[Path]
        Sorted (lexicographic) list of matched cube file paths,
        after slicing.

    Raises
    ------
    FileNotFoundError
        If no files match the pattern.
    """
    workdir = (workdir or Path(".")).resolve()
    cube_paths = [Path(p) for p in sorted(workdir.glob(cube_pattern))]
    logger.debug("Discovered %d cube files matching %r", len(cube_paths), cube_pattern)
    if not cube_paths:
        raise FileNotFoundError(
            f"No cube files matched pattern: {cube_pattern!r} in {workdir}"
        )
    return cube_paths[frame_start:frame_end:frame_step]


_CUBE_STEP_RE = re.compile(r"_(\d+)\.cube$")


def read_cube_atoms(path: Path, header: CubeHeader):
    """Parse atomic positions from a cube file header and return ``ase.Atoms``.

    Reuses the already-parsed *header* for cell vectors; only re-reads the
    atom lines from *path*.

    Parameters
    ----------
    path
        Path to the cube file.
    header
        Previously parsed ``CubeHeader`` from ``read_cube_header_and_values``.

    Returns
    -------
    ase.Atoms
        Atoms object with positions (Å) and cell set from cube header.
    """
    from ase import Atoms as _Atoms
    from ase.data import chemical_symbols

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for _ in range(6):  # 2 comments + natoms/origin + 3 voxel lines
            f.readline()
        symbols: list[str] = []
        positions_ang: list[list[float]] = []
        for _ in range(header.natoms):
            parts = f.readline().split()
            z_num = int(parts[0])
            x, y, z = _float(parts[2]), _float(parts[3]), _float(parts[4])
            symbols.append(chemical_symbols[z_num])
            positions_ang.append([x * BOHR_TO_ANG, y * BOHR_TO_ANG, z * BOHR_TO_ANG])

    cell_ang = np.array([
        header.vx_bohr * header.nx * BOHR_TO_ANG,
        header.vy_bohr * header.ny * BOHR_TO_ANG,
        header.vz_bohr * header.nz * BOHR_TO_ANG,
    ])
    return _Atoms(symbols=symbols, positions=positions_ang, cell=cell_ang, pbc=True)


def extract_step_from_cube_filename(path: Path) -> Optional[int]:
    """Extract the MD step number from a cube filename.

    Typical: ``md-POTENTIAL-v_hartree-1_1050.cube`` → ``1050``.
    Returns ``None`` if the pattern is not matched.
    """
    m = _CUBE_STEP_RE.search(path.name)
    if m:
        return int(m.group(1))
    return None
