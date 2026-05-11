"""Potential frame data source — thin compatibility wrapper (Phase 7b2).

The canonical implementation now lives in :mod:`md_analysis.engines.cp2k`:

- ``read_continuous_potential_frames``  (was ``discover_continuous_frames``)
- ``read_distributed_potential_frames`` (was ``discover_distributed_frames``)

This module keeps the legacy ``discover_*_frames`` names as thin
forwarding wrappers so existing call sites in
``electrochemical/potential/CenterPotential.py`` /
``electrochemical/potential/PhiZProfile.py`` /
``electrochemical/potential/__init__.py`` continue to work unchanged
during the engines/ migration.  All discovery, slicing and frame
construction logic is delegated; this file no longer carries any CP2K
parsing detail.

``PotentialFrame`` is the engine-neutral dataclass from
``md_analysis.engines.models`` and is re-exported here for backwards
compatibility.

The CP2K stdout / xyz parser primitives (``FERMI_RE`` /
``parse_md_out_fermi`` / ``read_xyz_atoms_for_steps`` etc.) are no
longer re-exported from this module; consumers should import them from
``md_analysis.utils.formats.cp2k_stdout`` / ``cp2k_xyz`` directly.
"""

from __future__ import annotations

from pathlib import Path

from ...engines.cp2k import (
    read_continuous_potential_frames,
    read_distributed_potential_frames,
)
from ...engines.models import PotentialFrame


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
    """Legacy name; delegates to
    :func:`md_analysis.engines.cp2k.read_continuous_potential_frames`."""
    return read_continuous_potential_frames(
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
    """Legacy name; delegates to
    :func:`md_analysis.engines.cp2k.read_distributed_potential_frames`."""
    return read_distributed_potential_frames(
        root_dir,
        dir_pattern=dir_pattern,
        cube_filename=cube_filename,
        sp_out_filename=sp_out_filename,
        center_mode=center_mode,
        metal_elements=metal_elements,
        layer_tol_ang=layer_tol_ang,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        verbose=verbose,
    )


__all__ = [
    "PotentialFrame",
    "discover_continuous_frames",
    "discover_distributed_frames",
]
