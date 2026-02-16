"""Integration test: generate an integrated three-panel water analysis figure."""

from __future__ import annotations

from pathlib import Path
import sys

import pytest

pytest.importorskip("matplotlib")

REPO_ROOT = Path(__file__).resolve().parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.structure.Analysis import plot_water_three_panel_analysis


def test_water_three_panel_plot_generation() -> None:
    xyz_path = REPO_ROOT / "data_example" / "potential" / "md-pos-1.xyz"
    md_inp_path = REPO_ROOT / "data_example" / "potential" / "md.inp"
    out_dir = REPO_ROOT / "test" / "_tmp_preview"
    out_dir.mkdir(parents=True, exist_ok=True)

    out_png_path = plot_water_three_panel_analysis(
        xyz_path=xyz_path,
        md_inp_path=md_inp_path,
        output_dir=out_dir,
        output_png_name="water_three_panel_analysis.png",
    )

    assert out_png_path.exists()
    assert out_png_path.stat().st_size > 0


if __name__ == "__main__":
    test_water_three_panel_plot_generation()
    print("saved: test/_tmp_preview/water_three_panel_analysis.png")
