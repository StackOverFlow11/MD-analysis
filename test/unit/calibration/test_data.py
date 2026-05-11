"""Tests for calibration data loading and JSON I/O."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from md_analysis.electrochemical.calibration._data import (
    CalibrationData,
    _detect_header,
    load_calibration_csv,
    load_calibration_json,
    save_calibration_json,
)

DATA_DIR = Path(__file__).parent / "data"


# ---------------------------------------------------------------------------
# CalibrationData validation
# ---------------------------------------------------------------------------

class TestCalibrationData:

    def test_valid_creation(self):
        cd = CalibrationData(
            potentials=np.array([0.0, 0.5]),
            charge_densities=np.array([-1.0, 1.0]),
        )
        assert cd.n_points == 2

    def test_mismatched_lengths_raises(self):
        with pytest.raises(ValueError, match="equal length"):
            CalibrationData(
                potentials=np.array([0.0, 0.5, 1.0]),
                charge_densities=np.array([-1.0, 1.0]),
            )

    def test_single_point_raises(self):
        with pytest.raises(ValueError, match="At least 2"):
            CalibrationData(
                potentials=np.array([0.5]),
                charge_densities=np.array([1.0]),
            )

    def test_coerces_to_float(self):
        cd = CalibrationData(
            potentials=[1, 2],
            charge_densities=[3, 4],
        )
        assert cd.potentials.dtype == float


# ---------------------------------------------------------------------------
# Header detection
# ---------------------------------------------------------------------------

class TestDetectHeader:

    def test_numeric_row(self):
        assert _detect_header(["0.1", "-3.2"]) is False

    def test_text_row(self):
        assert _detect_header(["potential_V", "sigma_uC_cm2"]) is True

    def test_mixed_row(self):
        assert _detect_header(["potential_V", "0.5"]) is True


# ---------------------------------------------------------------------------
# CSV loading
# ---------------------------------------------------------------------------

class TestLoadCSV:

    def test_with_header(self):
        cd = load_calibration_csv(DATA_DIR / "calibration_with_header.csv")
        assert cd.n_points == 4
        np.testing.assert_allclose(cd.potentials[0], -0.2)
        np.testing.assert_allclose(cd.charge_densities[0], -8.5)

    def test_without_header(self):
        cd = load_calibration_csv(DATA_DIR / "calibration_no_header.csv")
        assert cd.n_points == 4
        np.testing.assert_allclose(cd.potentials[-1], 0.7)

    def test_same_data_regardless_of_header(self):
        cd_h = load_calibration_csv(DATA_DIR / "calibration_with_header.csv")
        cd_n = load_calibration_csv(DATA_DIR / "calibration_no_header.csv")
        np.testing.assert_array_equal(cd_h.potentials, cd_n.potentials)
        np.testing.assert_array_equal(cd_h.charge_densities, cd_n.charge_densities)

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            load_calibration_csv("/nonexistent/path.csv")


# ---------------------------------------------------------------------------
# JSON round-trip
# ---------------------------------------------------------------------------

class TestJsonRoundTrip:

    def test_save_and_load(self, tmp_path: Path):
        cd = CalibrationData(
            potentials=np.array([0.0, 0.5, 1.0]),
            charge_densities=np.array([-5.0, 0.0, 5.0]),
            reference="RHE",
            metadata={"system": "test"},
        )
        fit_params = {"method": "linear", "slope": 0.1, "intercept": 0.5}

        json_path = save_calibration_json(cd, fit_params, tmp_path / "cal.json")
        assert json_path.exists()

        cd2, fp2 = load_calibration_json(json_path)
        assert cd2.n_points == 3
        assert cd2.reference == "RHE"
        assert cd2.metadata["system"] == "test"
        np.testing.assert_allclose(cd2.potentials, cd.potentials)
        np.testing.assert_allclose(cd2.charge_densities, cd.charge_densities)
        assert fp2["method"] == "linear"
        assert fp2["slope"] == pytest.approx(0.1)

    def test_load_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            load_calibration_json("/nonexistent/cal.json")
