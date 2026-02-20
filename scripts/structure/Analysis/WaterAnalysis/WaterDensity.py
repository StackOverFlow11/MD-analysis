"""High-level water density analysis along interface-to-midpoint direction."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np

from ...utils.config import DEFAULT_Z_BIN_WIDTH_A
from ..config import DEFAULT_OUTPUT_DIR, DEFAULT_WATER_MASS_DENSITY_CSV_NAME
from ..config import DEFAULT_START_INTERFACE
from ._common import (
    StartInterface,
    _compute_density_orientation_ensemble,
    _detect_low_high_interface_fractions,
    _parse_abc_from_md_inp,
)

# Re-export StartInterface for backward compatibility with external importers.
__all__ = ["StartInterface", "water_mass_density_z_distribution_analysis"]


def water_mass_density_z_distribution_analysis(
    xyz_path: str | Path,
    md_inp_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    output_csv_name: str = DEFAULT_WATER_MASS_DENSITY_CSV_NAME,
    start_interface: StartInterface = DEFAULT_START_INTERFACE,
    dz_A: float = DEFAULT_Z_BIN_WIDTH_A,
    metal_symbols: Iterable[str] | None = None,
) -> Path:
    """
    Ensemble-average water mass density profile from a selected metal interface
    toward the midpoint between two interfaces.

    Output CSV columns:
    - path_fraction_center: normalized distance in [0, 1] along interface->midpoint
    - distance_A: path_fraction_center * mean_path_length_A
    - rho_ensemble_avg_g_cm3: ensemble-averaged water mass density
    """
    xyz_path = Path(xyz_path)
    md_inp_path = Path(md_inp_path)
    output_dir_path = Path(output_dir) if output_dir is not None else Path.cwd()
    output_dir_path.mkdir(parents=True, exist_ok=True)

    common_centers_u, mean_path_A, rho_ensemble, _ = _compute_density_orientation_ensemble(
        xyz_path,
        md_inp_path,
        start_interface=start_interface,
        dz_A=dz_A,
        metal_symbols=metal_symbols,
    )
    distance_A = common_centers_u * mean_path_A

    out_csv_path = output_dir_path / output_csv_name
    np.savetxt(
        out_csv_path,
        np.column_stack([common_centers_u, distance_A, rho_ensemble]),
        delimiter=",",
        header="path_fraction_center,distance_A,rho_ensemble_avg_g_cm3",
        comments="",
    )
    return out_csv_path
