"""Integration tests for phi_z_planeavg_analysis.

Uses data_example/potential/ as input data.
Can also be run as a standalone script: python test/integration/potential/test_phi_z_profile.py
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from src.potential import phi_z_planeavg_analysis

# Resolve data directory relative to this file
_DATA_DIR = Path(__file__).resolve().parents[3] / "data_example" / "potential"

pytestmark = pytest.mark.skipif(
    not _DATA_DIR.exists(),
    reason=f"data_example/potential/ not found at {_DATA_DIR}",
)


class TestPhiZPlaneavgAnalysis:

    def test_basic(self, tmp_path: Path) -> None:
        old_cwd = os.getcwd()
        try:
            os.chdir(_DATA_DIR)
            png_path = phi_z_planeavg_analysis(
                "md-POTENTIAL-v_hartree-1_*.cube",
                output_dir=tmp_path,
                z_mavg_window_ang=7.0,
            )
        finally:
            os.chdir(old_cwd)

        assert png_path.exists()
        assert png_path.suffix == ".png"
        # Should also produce stats CSV
        assert (tmp_path / "phi_z_planeavg_stats.csv").exists()
        assert (tmp_path / "phi_z_planeavg_mavg_7A_stats.csv").exists()


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
        png = phi_z_planeavg_analysis(
            "md-POTENTIAL-v_hartree-1_*.cube",
            output_dir=_PREVIEW / "phi_z",
        )
    finally:
        os.chdir(old_cwd)

    print(f"φ(z) plot: {png}")
    sys.exit(0)
