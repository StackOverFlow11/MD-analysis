"""Core data structure for trajectory-level Bader charge data in XYZ order."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ....scripts.utils.IndexMapper import read_index_map_from_poscar, remap_array
from ....utils.BaderParser import load_bader_atoms
from ..config import (
    DEFAULT_ACF_FILENAME,
    DEFAULT_DIR_PATTERN,
    DEFAULT_POTCAR_FILENAME,
    DEFAULT_STRUCTURE_FILENAME,
)
from ._frame_utils import _extract_step_and_time, _sorted_frame_dirs

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class BaderTrajectoryData:
    """Trajectory Bader net charges remapped to original XYZ atom ordering.

    Attributes
    ----------
    steps : np.ndarray
        Shape ``(n_frames,)``, MD step numbers (int).
    times : np.ndarray
        Shape ``(n_frames,)``, simulation time in fs (int, from directory name).
    atom_indices_xyz : np.ndarray
        Shape ``(n_atoms,)``, ``[0, 1, ..., N-1]`` in XYZ ordering.
    net_charges : np.ndarray
        Shape ``(n_frames, n_atoms)``, Bader net charges remapped to XYZ order.
    """

    steps: np.ndarray
    times: np.ndarray
    atom_indices_xyz: np.ndarray
    net_charges: np.ndarray


def load_bader_trajectory(
    root_dir: str | Path = ".",
    *,
    dir_pattern: str = DEFAULT_DIR_PATTERN,
    structure_filename: str = DEFAULT_STRUCTURE_FILENAME,
    acf_filename: str = DEFAULT_ACF_FILENAME,
    potcar_filename: str = DEFAULT_POTCAR_FILENAME,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> BaderTrajectoryData:
    """Load Bader data across trajectory frames, remapped to XYZ atom order.

    Parameters
    ----------
    root_dir
        Parent directory containing per-frame subdirectories.
    dir_pattern
        Glob pattern for frame subdirectories (default ``bader_t*_i*``).
    structure_filename, acf_filename, potcar_filename
        Per-frame file names.
    frame_start, frame_end, frame_step
        Slice parameters applied to the sorted frame list.
    verbose
        Show tqdm progress bar.

    Returns
    -------
    BaderTrajectoryData
        Frozen dataclass with steps, times, atom indices, and net charges
        all in original XYZ ordering.
    """
    root = Path(root_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"root_dir does not exist: {root}")

    frame_dirs = _sorted_frame_dirs(root, dir_pattern)
    frame_dirs = frame_dirs[frame_start:frame_end:frame_step]
    if not frame_dirs:
        raise FileNotFoundError("No frame directories after slicing")

    logger.info("Loading Bader trajectory: %d frames", len(frame_dirs))

    steps_list: list[int] = []
    times_list: list[int] = []
    charges_list: list[np.ndarray] = []

    iterator: enumerate = enumerate(frame_dirs)
    if verbose:
        from tqdm import tqdm
        iterator = enumerate(tqdm(frame_dirs, desc="Bader load", unit="frame", ascii=" ="))

    n_atoms: int | None = None

    for _idx, frame_dir in iterator:
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
        imap = read_index_map_from_poscar(poscar)

        # Remap net charges from POSCAR order → XYZ order
        net_charge_poscar = atoms.arrays["bader_net_charge"]
        net_charge_xyz = remap_array(net_charge_poscar, imap, "poscar_to_xyz")

        if n_atoms is None:
            n_atoms = len(atoms)
        elif len(atoms) != n_atoms:
            raise ValueError(
                f"Frame {fname} has {len(atoms)} atoms, "
                f"expected {n_atoms} (from first frame)"
            )

        step, time_fs = _extract_step_and_time(fname)
        steps_list.append(step)
        times_list.append(time_fs)
        charges_list.append(net_charge_xyz)

    return BaderTrajectoryData(
        steps=np.array(steps_list, dtype=int),
        times=np.array(times_list, dtype=int),
        atom_indices_xyz=np.arange(n_atoms, dtype=int),
        net_charges=np.stack(charges_list, axis=0),
    )
