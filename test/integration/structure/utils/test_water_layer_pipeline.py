"""Integration tests for structure utils pipeline."""

from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pytest

ase_io = pytest.importorskip("ase.io")

from scripts.structure.utils.LayerParser import detect_interface_layers
from scripts.structure.utils.WaterParser import compute_water_mass_density_z_distribution
from scripts.structure.utils.WaterParser import compute_water_orientation_weighted_density_z_distribution
from scripts.structure.utils.WaterParser import detect_water_molecule_indices
from scripts.structure.utils.WaterParser import get_water_oxygen_indices_array


def _parse_abc_from_md_inp(md_inp_path: Path) -> tuple[float, float, float]:
    text = md_inp_path.read_text(encoding="utf-8")
    match = re.search(r"^\s*ABC\s+\[angstrom\]\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*$", text, re.MULTILINE)
    if not match:
        raise ValueError(f"Cannot find `ABC [angstrom] ...` in {md_inp_path}")
    return float(match.group(1)), float(match.group(2)), float(match.group(3))


def test_last_frame_pipeline_outputs_have_consistent_shapes() -> None:
    repo_root = Path(__file__).resolve().parents[4]
    xyz_path = repo_root / "data_example" / "potential" / "md-pos-1.xyz"
    md_inp_path = repo_root / "data_example" / "potential" / "md.inp"

    atoms = ase_io.read(str(xyz_path), index=-1)
    a, b, c = _parse_abc_from_md_inp(md_inp_path)
    atoms.set_cell([a, b, c])
    atoms.set_pbc([True, True, True])

    layer_result = detect_interface_layers(atoms, normal="c")
    assert len(layer_result.metal_layers_sorted) >= 2

    water_idx = detect_water_molecule_indices(atoms)
    oxygen_n1 = get_water_oxygen_indices_array(water_idx)
    rho = compute_water_mass_density_z_distribution(atoms, oxygen_n1)
    orient = compute_water_orientation_weighted_density_z_distribution(
        atoms,
        oxygen_n1,
        water_molecule_indices=water_idx,
    )

    assert rho.ndim == 2 and rho.shape[1] == 1
    assert orient.ndim == 2 and orient.shape[1] == 1
    assert rho.shape == orient.shape
    assert np.all(np.isfinite(rho))
    assert np.all(np.isfinite(orient))
