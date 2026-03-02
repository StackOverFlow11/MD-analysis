"""Integration tests for center_slab_potential_analysis and related functions.

Uses data_example/potential/ as input data.
Can also be run as a standalone script: python test/integration/potential/test_center_potential.py
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import numpy as np
import pytest

from src.potential import (
    center_slab_potential_analysis,
    fermi_energy_analysis,
    electrode_potential_analysis,
)

# Resolve data directory relative to this file
_DATA_DIR = Path(__file__).resolve().parents[3] / "data_example" / "potential"

# Skip all tests if data directory doesn't exist
pytestmark = pytest.mark.skipif(
    not _DATA_DIR.exists(),
    reason=f"data_example/potential/ not found at {_DATA_DIR}",
)


class TestCenterSlabPotentialAnalysis:

    def test_basic_cell_center(self, tmp_path: Path) -> None:
        """Run with center_mode='cell' (no xyz needed)."""
        old_cwd = os.getcwd()
        try:
            os.chdir(_DATA_DIR)
            csv_path = center_slab_potential_analysis(
                "md-POTENTIAL-v_hartree-1_*.cube",
                output_dir=tmp_path,
                thickness_ang=7.0,
                center_mode="cell",
            )
        finally:
            os.chdir(old_cwd)

        assert csv_path.exists()
        data = np.genfromtxt(csv_path, delimiter=",", names=True, dtype=None, encoding="utf-8")
        assert data.size > 0
        assert "step" in data.dtype.names
        assert "phi_center_ev" in data.dtype.names

    def test_interface_center(self, tmp_path: Path) -> None:
        """Run with center_mode='interface' using xyz trajectory."""
        old_cwd = os.getcwd()
        try:
            os.chdir(_DATA_DIR)
            csv_path = center_slab_potential_analysis(
                "md-POTENTIAL-v_hartree-1_*.cube",
                output_dir=tmp_path,
                thickness_ang=7.0,
                center_mode="interface",
                xyz_path=_DATA_DIR / "md-pos-1.xyz",
            )
        finally:
            os.chdir(old_cwd)

        assert csv_path.exists()
        # Should also produce slab_center_and_interfaces.csv
        slab_csv = tmp_path / "slab_center_and_interfaces.csv"
        assert slab_csv.exists()


class TestFermiEnergyAnalysis:

    def test_basic(self, tmp_path: Path) -> None:
        csv_path = fermi_energy_analysis(
            _DATA_DIR / "md.out",
            output_dir=tmp_path,
            fermi_unit="au",
        )
        assert csv_path.exists()
        data = np.genfromtxt(csv_path, delimiter=",", names=True, dtype=None, encoding="utf-8")
        assert data.size > 0
        assert "fermi_ev" in data.dtype.names


class TestElectrodePotentialAnalysis:

    def test_full_pipeline(self, tmp_path: Path) -> None:
        old_cwd = os.getcwd()
        try:
            os.chdir(_DATA_DIR)
            u_csv = electrode_potential_analysis(
                "md-POTENTIAL-v_hartree-1_*.cube",
                _DATA_DIR / "md.out",
                output_dir=tmp_path,
                thickness_ang=7.0,
                center_mode="interface",
                xyz_path=_DATA_DIR / "md-pos-1.xyz",
                fermi_unit="au",
            )
        finally:
            os.chdir(old_cwd)

        assert u_csv.exists()
        # Should also have center_potential.csv and fermi_energy.csv
        assert (tmp_path / "center_potential.csv").exists()
        assert (tmp_path / "fermi_energy.csv").exists()


# ---------------------------------------------------------------------------
# Standalone runner
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    _PREVIEW = Path(__file__).resolve().parents[2] / "_tmp_preview"
    _PREVIEW.mkdir(exist_ok=True)

    old_cwd = os.getcwd()
    try:
        os.chdir(_DATA_DIR)
        csv_path = electrode_potential_analysis(
            "md-POTENTIAL-v_hartree-1_*.cube",
            _DATA_DIR / "md.out",
            output_dir=_PREVIEW / "potential",
            thickness_ang=7.0,
            center_mode="interface",
            xyz_path=_DATA_DIR / "md-pos-1.xyz",
        )
    finally:
        os.chdir(old_cwd)

    print(f"Electrode potential CSV: {csv_path}")
    sys.exit(0)
