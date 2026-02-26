"""Integration test: generate integrated three-panel analysis outputs.

This file targets the revised `scripts.structure.Analysis` API and validates:
1) `plot_water_three_panel_analysis(...)` writes one PNG
2) side-output CSV/TXT files are present and content is sane
"""

from __future__ import annotations

from pathlib import Path
import re
import sys

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
pytest.importorskip("ase")

REPO_ROOT = Path(__file__).resolve().parents[4]

from src.structure import DEFAULT_THETA_BIN_DEG
from src.structure.Analysis import DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME
from src.structure.Analysis import DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME
from src.structure.Analysis import DEFAULT_ADSORBED_WATER_THETA_DISTRIBUTION_CSV_NAME
from src.structure.Analysis import DEFAULT_WATER_THREE_PANEL_PLOT_PNG_NAME
from src.structure.Analysis import DEFAULT_WATER_MASS_DENSITY_CSV_NAME
from src.structure.Analysis import DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME
from src.structure.Analysis import plot_water_three_panel_analysis


def _as_2d(data: np.ndarray) -> np.ndarray:
    data = np.asarray(data)
    if data.ndim == 1:
        return data.reshape(1, -1)
    return data


def _run(out_dir: Path) -> dict[str, Path]:
    out_dir = Path(out_dir)
    xyz_path = REPO_ROOT / "data_example" / "potential" / "md-pos-1.xyz"
    md_inp_path = REPO_ROOT / "data_example" / "potential" / "md.inp"
    out_dir.mkdir(parents=True, exist_ok=True)

    out_png_path = plot_water_three_panel_analysis(
        xyz_path=xyz_path,
        md_inp_path=md_inp_path,
        output_dir=out_dir,
        output_png_name=DEFAULT_WATER_THREE_PANEL_PLOT_PNG_NAME,
    )

    return {
        "png": out_png_path,
        "density_csv": out_dir / DEFAULT_WATER_MASS_DENSITY_CSV_NAME,
        "orientation_csv": out_dir / DEFAULT_WATER_ORIENTATION_WEIGHTED_DENSITY_CSV_NAME,
        "adsorbed_profile_csv": out_dir / DEFAULT_ADSORBED_WATER_PROFILE_CSV_NAME,
        "adsorbed_range_txt": out_dir / DEFAULT_ADSORBED_WATER_RANGE_TXT_NAME,
        "theta_csv": out_dir / DEFAULT_ADSORBED_WATER_THETA_DISTRIBUTION_CSV_NAME,
    }


def test_water_three_panel_plot_generation() -> None:
    out_dir = REPO_ROOT / "test" / "_tmp_preview"
    outputs = _run(out_dir)

    # --- File existence ---
    for name, path in outputs.items():
        assert path.exists(), f"missing output: {name} -> {path}"
        assert path.stat().st_size > 0, f"empty output: {name} -> {path}"

    # --- Basic content sanity checks (shapes / finite / non-negative where applicable) ---
    density = _as_2d(np.loadtxt(outputs["density_csv"], delimiter=",", skiprows=1))
    assert density.shape[1] == 3
    assert np.all(np.isfinite(density))
    assert np.all(density[:, 2] >= 0.0)
    assert np.all(density[:, 0] >= 0.0) and np.all(density[:, 0] <= 1.0)

    orientation = _as_2d(np.loadtxt(outputs["orientation_csv"], delimiter=",", skiprows=1))
    assert orientation.shape[1] == 3
    assert np.all(np.isfinite(orientation))
    assert np.all(orientation[:, 0] >= 0.0) and np.all(orientation[:, 0] <= 1.0)

    ads_profile = _as_2d(np.loadtxt(outputs["adsorbed_profile_csv"], delimiter=",", skiprows=1))
    assert ads_profile.shape[1] == 4
    assert np.all(np.isfinite(ads_profile[:, 0:3]))
    mask = ads_profile[:, 3]
    mask_rounded = np.rint(mask).astype(int)
    assert np.allclose(mask, mask_rounded)
    assert set(np.unique(mask_rounded).tolist()).issubset({0, 1})

    # --- Range TXT sanity ---
    txt = outputs["adsorbed_range_txt"].read_text(encoding="utf-8")
    m_start = re.search(r"adsorbed_layer_start_A=([0-9.eE+-]+)", txt)
    m_end = re.search(r"adsorbed_layer_end_A=([0-9.eE+-]+)", txt)
    m_peak = re.search(r"main_peak_distance_A=([0-9.eE+-]+)", txt)
    assert m_start and m_end and m_peak
    d_start = float(m_start.group(1))
    d_end = float(m_end.group(1))
    d_peak = float(m_peak.group(1))
    assert d_end > d_start
    assert d_start <= d_peak <= d_end
    assert "near_zero_ratio=" in txt
    assert "smoothing_window_bins=" in txt

    # --- Theta distribution sanity ---
    theta = _as_2d(np.loadtxt(outputs["theta_csv"], delimiter=",", skiprows=1))
    assert theta.shape[1] == 2
    assert np.all(np.isfinite(theta))
    assert np.all(theta[:, 1] >= 0.0)
    assert np.any(theta[:, 1] > 0.0)
    assert np.isclose(float(np.sum(theta[:, 1]) * float(DEFAULT_THETA_BIN_DEG)), 1.0, rtol=1.0e-2, atol=1.0e-2)


if __name__ == "__main__":
    out_dir = REPO_ROOT / "test" / "_tmp_preview"
    outputs = _run(out_dir)
    for name, path in outputs.items():
        try:
            rel = path.relative_to(REPO_ROOT)
        except Exception:
            rel = path
        print(f"saved ({name}): {rel}")
