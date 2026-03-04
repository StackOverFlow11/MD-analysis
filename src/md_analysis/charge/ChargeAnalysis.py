"""Surface charge density and per-atom charge analysis from Bader output."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
from ase import Atoms

from ..utils.BaderParser import load_bader_atoms
from ..utils.LayerParser import detect_interface_layers
from .config import (
    DEFAULT_ACF_FILENAME,
    DEFAULT_DIR_PATTERN,
    DEFAULT_POTCAR_FILENAME,
    DEFAULT_STRUCTURE_FILENAME,
    E_PER_A2_TO_UC_PER_CM2,
)

# Map normal axis to the two cell-vector indices spanning the surface plane
_AREA_VECTORS = {"a": (1, 2), "b": (0, 2), "c": (0, 1)}


# ---------------------------------------------------------------------------
# Single-frame analysis
# ---------------------------------------------------------------------------

def compute_frame_surface_charge(
    atoms: Atoms,
    *,
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
) -> Atoms:
    """Compute surface charge density for a single frame.

    The input ``atoms`` must already carry ``"bader_net_charge"`` in
    ``atoms.arrays`` (as produced by :func:`load_bader_atoms`).

    Results are stored in-place on the returned Atoms:
    - ``atoms.info["surface_charge_density_e_A2"]`` — [σ_bottom, σ_top] in e/Å²
    - ``atoms.info["surface_charge_density_uC_cm2"]`` — same in μC/cm²

    Parameters
    ----------
    atoms
        ASE Atoms with ``"bader_net_charge"`` per-atom array.
    metal_symbols
        Override default metal symbols for layer detection.
    normal
        Cell axis perpendicular to the surface (``"a"``, ``"b"``, or ``"c"``).

    Returns
    -------
    Atoms
        The same object, mutated with results in ``atoms.info``.
    """
    if normal not in _AREA_VECTORS:
        raise ValueError(
            f"normal must be one of {set(_AREA_VECTORS)}, got {normal!r}"
        )

    if "bader_net_charge" not in atoms.arrays:
        raise ValueError(
            "atoms.arrays must contain 'bader_net_charge'. "
            "Use load_bader_atoms() first."
        )

    net_charge = atoms.arrays["bader_net_charge"]

    # Detect interface layers
    det = detect_interface_layers(atoms, metal_symbols=metal_symbols, normal=normal)
    iface_layers = det.interface_layers()
    if len(iface_layers) != 2:
        raise ValueError(
            f"Expected exactly 2 interface layers, got {len(iface_layers)}"
        )

    # Surface area from the two cell vectors spanning the surface plane
    i0, i1 = _AREA_VECTORS[normal]
    cell = np.asarray(atoms.cell.array, dtype=float)
    area_A2 = float(np.linalg.norm(np.cross(cell[i0], cell[i1])))

    # Sort interface layers by center_s (bottom first)
    sorted_iface = sorted(iface_layers, key=lambda L: L.center_s)

    sigma_e_A2 = np.empty(2)
    for i, layer in enumerate(sorted_iface):
        idx = np.array(layer.atom_indices)
        sigma_e_A2[i] = net_charge[idx].sum() / area_A2

    sigma_uC_cm2 = sigma_e_A2 * E_PER_A2_TO_UC_PER_CM2

    atoms.info["surface_charge_density_e_A2"] = sigma_e_A2.tolist()
    atoms.info["surface_charge_density_uC_cm2"] = sigma_uC_cm2.tolist()

    return atoms


# ---------------------------------------------------------------------------
# Trajectory indexed atom charges
# ---------------------------------------------------------------------------

def trajectory_indexed_atom_charges(
    root_dir: str | Path,
    atom_index_matrix: np.ndarray,
    *,
    dir_pattern: str = DEFAULT_DIR_PATTERN,
    structure_filename: str = DEFAULT_STRUCTURE_FILENAME,
    acf_filename: str = DEFAULT_ACF_FILENAME,
    potcar_filename: str = DEFAULT_POTCAR_FILENAME,
) -> np.ndarray:
    """Extract net charges for specified atom indices across trajectory frames.

    Parameters
    ----------
    root_dir
        Parent directory containing per-frame subdirectories.
    atom_index_matrix
        2-D integer array of shape ``(t, N)`` with 0-based atom indices.
        Row *i* lists the atom indices to query in frame *i*.
    dir_pattern
        Glob pattern for discovering frame subdirectories.
    structure_filename
        Name of the structure file (POSCAR) in each subdirectory.
    acf_filename
        Name of the ACF.dat file in each subdirectory.
    potcar_filename
        Name of the POTCAR file in each subdirectory.

    Returns
    -------
    np.ndarray
        Array of shape ``(t, N, 2)`` where ``[:, :, 0]`` contains the
        (echoed-back) atom indices and ``[:, :, 1]`` contains the
        corresponding Bader net charges.
    """
    # --- Validate atom_index_matrix ---
    arr = np.asarray(atom_index_matrix)
    if arr.ndim != 2:
        raise ValueError(
            f"atom_index_matrix must be 2-D, got {arr.ndim}-D"
        )
    if not np.issubdtype(arr.dtype, np.integer):
        raise ValueError(
            f"atom_index_matrix must have integer dtype, got {arr.dtype}"
        )
    if np.any(arr < 0):
        raise ValueError("atom_index_matrix contains negative indices")

    root = Path(root_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"root_dir does not exist: {root}")

    frame_dirs = sorted(root.glob(dir_pattern))
    if not frame_dirs:
        raise FileNotFoundError(
            f"No subdirectories matching '{dir_pattern}' in {root}"
        )

    t_expected = arr.shape[0]
    if t_expected != len(frame_dirs):
        raise ValueError(
            f"atom_index_matrix has {t_expected} rows but found "
            f"{len(frame_dirs)} frame directories"
        )

    t, n = arr.shape
    result = np.empty((t, n, 2), dtype=float)

    for i, frame_dir in enumerate(frame_dirs):
        fname = frame_dir.name
        poscar = frame_dir / structure_filename
        acf = frame_dir / acf_filename
        potcar = frame_dir / potcar_filename

        for path, label in [(poscar, structure_filename),
                            (acf, acf_filename),
                            (potcar, potcar_filename)]:
            if not path.exists():
                raise FileNotFoundError(
                    f"{label} not found in frame {fname}: {path}"
                )

        atoms = load_bader_atoms(poscar, acf, potcar)
        net_charge = atoms.arrays["bader_net_charge"]
        n_atoms = len(atoms)

        indices_row = arr[i]
        oob = indices_row[indices_row >= n_atoms]
        if oob.size > 0:
            raise IndexError(
                f"Atom index {int(oob[0])} out of bounds for frame {fname} "
                f"with {n_atoms} atoms"
            )

        result[i, :, 0] = indices_row.astype(float)
        result[i, :, 1] = net_charge[indices_row]

    return result
