"""Surface charge density and per-atom charge analysis from Bader output."""

from __future__ import annotations

import csv
import warnings
from abc import ABC, abstractmethod
from dataclasses import dataclass
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
    DEFAULT_SELECTED_ATOM_CHARGES_CSV_NAME,
    DEFAULT_STRUCTURE_FILENAME,
    DEFAULT_SURFACE_CHARGE_CSV_NAME,
    E_PER_A2_TO_UC_PER_CM2,
)


# ---------------------------------------------------------------------------
# Atom selectors
# ---------------------------------------------------------------------------

class AtomSelector(ABC):
    """Abstract base class for selecting atoms of interest."""

    @abstractmethod
    def select(self, atoms: Atoms) -> np.ndarray:
        """Return 0-based indices of selected atoms. May return an empty array."""
        ...


class ElementSelector(AtomSelector):
    """Select atoms by element symbol(s)."""

    def __init__(self, symbols: Iterable[str]) -> None:
        self._symbols = frozenset(symbols)
        if not self._symbols:
            raise ValueError("symbols must be non-empty")

    def select(self, atoms: Atoms) -> np.ndarray:
        sym = np.asarray(atoms.get_chemical_symbols())
        mask = np.isin(sym, list(self._symbols))
        return np.where(mask)[0]


class IndexSelector(AtomSelector):
    """Select atoms by fixed 0-based indices."""

    def __init__(self, indices: Iterable[int]) -> None:
        self._indices = np.asarray(sorted(set(indices)), dtype=int)

    def select(self, atoms: Atoms) -> np.ndarray:
        n = len(atoms)
        valid = self._indices[(self._indices >= 0) & (self._indices < n)]
        return valid


# ---------------------------------------------------------------------------
# Trajectory result dataclass
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TrajectoryChargeResult:
    """Aggregated results from multi-frame charge analysis."""

    frame_labels: tuple[str, ...]
    surface_charge_density_uC_cm2: np.ndarray       # (n_frames, 2) [bottom, top]
    mean_surface_charge_density_uC_cm2: np.ndarray   # (2,)
    std_surface_charge_density_uC_cm2: np.ndarray    # (2,)
    selected_atom_net_charges: np.ndarray             # (n_frames, n_selected) or empty
    mean_selected_atom_net_charges: np.ndarray        # (n_selected,) or empty
    selected_atom_indices: tuple[int, ...]            # empty tuple if no selector


# ---------------------------------------------------------------------------
# Single-frame analysis
# ---------------------------------------------------------------------------

def compute_frame_surface_charge(
    atoms: Atoms,
    *,
    atom_selector: AtomSelector | None = None,
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
) -> Atoms:
    """Compute surface charge density for a single frame.

    The input ``atoms`` must already carry ``"bader_net_charge"`` in
    ``atoms.arrays`` (as produced by :func:`load_bader_atoms`).

    Results are stored in-place on the returned Atoms:
    - ``atoms.info["surface_charge_density_e_A2"]`` — [σ_bottom, σ_top] in e/Å²
    - ``atoms.info["surface_charge_density_uC_cm2"]`` — same in μC/cm²
    - ``atoms.info["selected_atom_indices"]`` — (if selector provided)
    - ``atoms.info["selected_atom_net_charges"]`` — (if selector provided)

    Parameters
    ----------
    atoms
        ASE Atoms with ``"bader_net_charge"`` per-atom array.
    atom_selector
        Optional selector for atoms of interest.
    metal_symbols
        Override default metal symbols for layer detection.
    normal
        Cell axis for layer detection (default ``"c"``).

    Returns
    -------
    Atoms
        The same object, mutated with results in ``atoms.info``.
    """
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

    # Surface area from cell vectors a × b (orthogonal cell, normal="c")
    cell = np.asarray(atoms.cell.array, dtype=float)
    a_vec = cell[0]
    b_vec = cell[1]
    area_A2 = float(np.linalg.norm(np.cross(a_vec, b_vec)))

    # Sort interface layers by center_s (bottom first)
    sorted_iface = sorted(iface_layers, key=lambda L: L.center_s)

    sigma_e_A2 = np.empty(2)
    for i, layer in enumerate(sorted_iface):
        idx = np.array(layer.atom_indices)
        sigma_e_A2[i] = net_charge[idx].sum() / area_A2

    sigma_uC_cm2 = sigma_e_A2 * E_PER_A2_TO_UC_PER_CM2

    atoms.info["surface_charge_density_e_A2"] = sigma_e_A2.tolist()
    atoms.info["surface_charge_density_uC_cm2"] = sigma_uC_cm2.tolist()

    # Optional atom selection
    if atom_selector is not None:
        sel_idx = atom_selector.select(atoms)
        atoms.info["selected_atom_indices"] = sel_idx.tolist()
        atoms.info["selected_atom_net_charges"] = net_charge[sel_idx].tolist()

    return atoms


# ---------------------------------------------------------------------------
# Trajectory analysis
# ---------------------------------------------------------------------------

def trajectory_charge_analysis(
    root_dir: str | Path,
    *,
    atom_selector: AtomSelector | None = None,
    metal_symbols: Iterable[str] | None = None,
    normal: str = "c",
    dir_pattern: str = DEFAULT_DIR_PATTERN,
    structure_filename: str = DEFAULT_STRUCTURE_FILENAME,
    acf_filename: str = DEFAULT_ACF_FILENAME,
    potcar_filename: str = DEFAULT_POTCAR_FILENAME,
    output_dir: str | Path | None = None,
    surface_charge_csv_name: str = DEFAULT_SURFACE_CHARGE_CSV_NAME,
    selected_atom_charges_csv_name: str = DEFAULT_SELECTED_ATOM_CHARGES_CSV_NAME,
) -> TrajectoryChargeResult:
    """Run charge analysis over a trajectory of Bader calculation directories.

    Parameters
    ----------
    root_dir
        Parent directory containing per-frame subdirectories.
    atom_selector
        Optional selector for tracking specific atom charges across frames.
    metal_symbols
        Override default metal symbols for layer detection.
    normal
        Cell axis for layer detection.
    dir_pattern
        Glob pattern for discovering frame subdirectories.
    structure_filename
        Name of the structure file (POSCAR) in each subdirectory.
    acf_filename
        Name of the ACF.dat file in each subdirectory.
    potcar_filename
        Name of the POTCAR file in each subdirectory.
    output_dir
        Directory for CSV output. Defaults to *root_dir*.
    surface_charge_csv_name
        Filename for the surface charge density CSV.
    selected_atom_charges_csv_name
        Filename for the selected-atom charges CSV.

    Returns
    -------
    TrajectoryChargeResult
    """
    root = Path(root_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"root_dir does not exist: {root}")

    frame_dirs = sorted(root.glob(dir_pattern))
    if not frame_dirs:
        raise FileNotFoundError(
            f"No subdirectories matching '{dir_pattern}' in {root}"
        )

    out_dir = Path(output_dir) if output_dir is not None else root
    out_dir.mkdir(parents=True, exist_ok=True)

    labels: list[str] = []
    sigma_list: list[np.ndarray] = []
    sel_charges_list: list[np.ndarray] = []
    sel_indices: np.ndarray | None = None

    for frame_dir in frame_dirs:
        acf_path = frame_dir / acf_filename
        if not acf_path.exists():
            warnings.warn(
                f"Skipping {frame_dir.name}: {acf_filename} not found",
                stacklevel=2,
            )
            continue

        atoms = load_bader_atoms(
            frame_dir / structure_filename,
            acf_path,
            frame_dir / potcar_filename,
        )
        compute_frame_surface_charge(
            atoms,
            atom_selector=atom_selector,
            metal_symbols=metal_symbols,
            normal=normal,
        )

        labels.append(frame_dir.name)
        sigma_list.append(
            np.array(atoms.info["surface_charge_density_uC_cm2"])
        )

        if atom_selector is not None:
            frame_sel_idx = np.array(atoms.info["selected_atom_indices"])
            frame_sel_q = np.array(atoms.info["selected_atom_net_charges"])

            if sel_indices is None:
                sel_indices = frame_sel_idx
            else:
                if not np.array_equal(sel_indices, frame_sel_idx):
                    raise ValueError(
                        f"Selected atom indices changed between frames: "
                        f"first frame {sel_indices.tolist()}, "
                        f"frame {frame_dir.name} {frame_sel_idx.tolist()}"
                    )
            sel_charges_list.append(frame_sel_q)

    if not labels:
        raise RuntimeError("No valid frames found in trajectory")

    sigma_arr = np.stack(sigma_list)  # (n_frames, 2)
    mean_sigma = sigma_arr.mean(axis=0)
    std_sigma = sigma_arr.std(axis=0)

    if sel_charges_list:
        sel_arr = np.stack(sel_charges_list)  # (n_frames, n_selected)
        mean_sel = sel_arr.mean(axis=0)
        final_sel_indices = tuple(int(i) for i in sel_indices)  # type: ignore[union-attr]
    else:
        sel_arr = np.empty((len(labels), 0))
        mean_sel = np.empty(0)
        final_sel_indices = ()

    # Write surface charge CSV
    _write_surface_charge_csv(
        out_dir / surface_charge_csv_name, labels, sigma_arr, mean_sigma, std_sigma
    )

    # Write selected atom charges CSV (if applicable)
    if final_sel_indices:
        _write_selected_atom_charges_csv(
            out_dir / selected_atom_charges_csv_name,
            labels,
            final_sel_indices,
            sel_arr,
            mean_sel,
        )

    return TrajectoryChargeResult(
        frame_labels=tuple(labels),
        surface_charge_density_uC_cm2=sigma_arr,
        mean_surface_charge_density_uC_cm2=mean_sigma,
        std_surface_charge_density_uC_cm2=std_sigma,
        selected_atom_net_charges=sel_arr,
        mean_selected_atom_net_charges=mean_sel,
        selected_atom_indices=final_sel_indices,
    )


# ---------------------------------------------------------------------------
# CSV writers
# ---------------------------------------------------------------------------

def _write_surface_charge_csv(
    path: Path,
    labels: list[str],
    sigma_arr: np.ndarray,
    mean_sigma: np.ndarray,
    std_sigma: np.ndarray,
) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["frame", "sigma_bottom_uC_cm2", "sigma_top_uC_cm2"])
        for i, label in enumerate(labels):
            writer.writerow([label, f"{sigma_arr[i, 0]:.6f}", f"{sigma_arr[i, 1]:.6f}"])
        writer.writerow(["mean", f"{mean_sigma[0]:.6f}", f"{mean_sigma[1]:.6f}"])
        writer.writerow(["std", f"{std_sigma[0]:.6f}", f"{std_sigma[1]:.6f}"])


def _write_selected_atom_charges_csv(
    path: Path,
    labels: list[str],
    indices: tuple[int, ...],
    sel_arr: np.ndarray,
    mean_sel: np.ndarray,
) -> None:
    header = ["frame"] + [f"atom_{idx}" for idx in indices]
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for i, label in enumerate(labels):
            writer.writerow([label] + [f"{v:.6f}" for v in sel_arr[i]])
        writer.writerow(["mean"] + [f"{v:.6f}" for v in mean_sel])
