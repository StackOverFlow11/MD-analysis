"""Unit tests for ``md_analysis.workflows.water`` single-step facades.

Phase 6.1 introduces 4 single-step run_* functions:
``run_water_density`` / ``run_water_orientation`` /
``run_ad_water_orientation`` / ``run_ad_water_theta``.

These tests monkeypatch the underlying business functions (mock,
no I/O) and pin the facade contract:

- ``WorkflowResult.name`` is the documented identifier
- ``output_dir`` is the caller-provided path, mkdir'd
- ``artifacts`` use the documented stable keys and point to the
  business-function return path(s)
- ``metadata`` contains ``cell_abc_provided`` plus the kwargs the
  CLI forwards (``dz_A`` / ``layer_tol_A`` / frame slice)

The composite ``run_water_three_panel`` has its own integration test
(`test/integration/water`) and is not re-covered here.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from md_analysis.workflows import (
    WorkflowResult,
    run_ad_water_orientation,
    run_ad_water_theta,
    run_water_density,
    run_water_orientation,
)


def _xyz_path(tmp_path: Path) -> Path:
    """Return a placeholder xyz path (file is never read in mocked tests)."""
    return tmp_path / "md-pos-1.xyz"


# ---------------------------------------------------------------------------
# run_water_density
# ---------------------------------------------------------------------------


class TestRunWaterDensity:
    def test_workflow_result_contract(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        captured: dict[str, Any] = {}

        def fake_business(**kwargs: Any) -> Path:
            captured.update(kwargs)
            return tmp_path / "water" / "water_mass_density.csv"

        monkeypatch.setattr(
            "md_analysis.water.water_mass_density_z_distribution_analysis",
            fake_business,
        )

        result = run_water_density(
            xyz_path=_xyz_path(tmp_path),
            cell_abc=(10.0, 10.0, 30.0),
            output_dir=tmp_path / "water",
            dz_A=0.1,
            layer_tol_A=0.5,
            frame_start=0,
            frame_end=100,
            frame_step=2,
        )

        # WorkflowResult shape
        assert isinstance(result, WorkflowResult)
        assert result.name == "water_density"
        assert result.output_dir == tmp_path / "water"
        assert (tmp_path / "water").exists()  # mkdir(parents=True)

        # artifact contract
        assert set(result.artifacts) == {"density_csv"}
        assert result.artifacts["density_csv"] == (
            tmp_path / "water" / "water_mass_density.csv"
        )

        # metadata contract
        assert result.metadata["cell_abc_provided"] is True
        assert result.metadata["dz_A"] == 0.1
        assert result.metadata["layer_tol_A"] == 0.5
        assert result.metadata["frame_start"] == 0
        assert result.metadata["frame_end"] == 100
        assert result.metadata["frame_step"] == 2

        # underlying business call kwargs (mode-agnostic forwarding)
        assert captured["output_dir"] == tmp_path / "water"
        assert captured["dz_A"] == 0.1


# ---------------------------------------------------------------------------
# run_water_orientation
# ---------------------------------------------------------------------------


class TestRunWaterOrientation:
    def test_workflow_result_contract(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        def fake_business(**kwargs: Any) -> Path:
            return tmp_path / "water" / "water_orientation.csv"

        monkeypatch.setattr(
            "md_analysis.water"
            ".water_orientation_weighted_density_z_distribution_analysis",
            fake_business,
        )

        result = run_water_orientation(
            xyz_path=_xyz_path(tmp_path),
            cell_abc=None,
            output_dir=tmp_path / "water",
        )

        assert isinstance(result, WorkflowResult)
        assert result.name == "water_orientation"
        assert result.output_dir == tmp_path / "water"
        assert set(result.artifacts) == {"orientation_csv"}
        assert result.artifacts["orientation_csv"] == (
            tmp_path / "water" / "water_orientation.csv"
        )
        # cell_abc=None should be reflected in metadata
        assert result.metadata["cell_abc_provided"] is False


# ---------------------------------------------------------------------------
# run_ad_water_orientation
# ---------------------------------------------------------------------------


class TestRunAdWaterOrientation:
    def test_workflow_result_contract_two_artifacts(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        profile_csv = tmp_path / "water" / "ad_water_profile.csv"
        range_txt = tmp_path / "water" / "ad_water_range.txt"

        def fake_business(**kwargs: Any) -> tuple[Path, Path]:
            return profile_csv, range_txt

        monkeypatch.setattr(
            "md_analysis.water.ad_water_orientation_analysis",
            fake_business,
        )

        result = run_ad_water_orientation(
            xyz_path=_xyz_path(tmp_path),
            cell_abc=(10.0, 10.0, 30.0),
            output_dir=tmp_path / "water",
        )

        assert isinstance(result, WorkflowResult)
        assert result.name == "ad_water_orientation"
        # Two artifacts: documented stable keys
        assert set(result.artifacts) == {
            "adsorbed_profile_csv",
            "adsorbed_range_txt",
        }
        assert result.artifacts["adsorbed_profile_csv"] == profile_csv
        assert result.artifacts["adsorbed_range_txt"] == range_txt


# ---------------------------------------------------------------------------
# run_ad_water_theta
# ---------------------------------------------------------------------------


class TestRunAdWaterTheta:
    def test_workflow_result_contract_drops_arrays(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """Business returns ``(centers, pdf, csv_path)``; the facade
        exposes only the CSV under ``theta_csv``.  The in-memory
        arrays are intentionally not surfaced (callers re-read CSV).
        """
        import numpy as np

        centers = np.linspace(0.0, 180.0, 73)
        pdf = np.zeros(73)
        csv_path = tmp_path / "water" / "ad_water_theta.csv"

        def fake_business(**kwargs: Any) -> tuple[np.ndarray, np.ndarray, Path]:
            return centers, pdf, csv_path

        monkeypatch.setattr(
            "md_analysis.water.compute_adsorbed_water_theta_distribution",
            fake_business,
        )

        result = run_ad_water_theta(
            xyz_path=_xyz_path(tmp_path),
            cell_abc=(10.0, 10.0, 30.0),
            output_dir=tmp_path / "water",
            verbose=True,
        )

        assert isinstance(result, WorkflowResult)
        assert result.name == "ad_water_theta"
        # Phase 6.1: artifact key stays ``theta_csv`` (codex-approved
        # to avoid CLI UI string churn).
        assert set(result.artifacts) == {"theta_csv"}
        assert result.artifacts["theta_csv"] == csv_path
        # Arrays are NOT in metadata / extra / artifacts
        assert "centers" not in result.metadata
        assert "pdf" not in result.metadata
        assert result.extra is None
