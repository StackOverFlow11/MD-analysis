"""End-to-end tests for calibration workflow."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest

from md_analysis.electrochemical.calibration import calibrate, predict_potential

DATA_DIR = Path(__file__).parent / "data"


class TestCalibrateFromCSV:

    def test_creates_json_and_png(self, tmp_path: Path):
        json_path = calibrate(
            csv_path=DATA_DIR / "calibration_with_header.csv",
            method="linear",
            output_dir=tmp_path / "out",
            calibration_json_path=tmp_path / "cal.json",
        )
        assert json_path.exists()
        assert (tmp_path / "out" / "calibration_data.csv").exists()
        assert (tmp_path / "out" / "calibration_fit.png").exists()

    def test_polynomial_method(self, tmp_path: Path):
        json_path = calibrate(
            csv_path=DATA_DIR / "calibration_with_header.csv",
            method="polynomial",
            poly_degree=2,
            output_dir=tmp_path / "out",
            calibration_json_path=tmp_path / "cal.json",
        )
        assert json_path.exists()

    def test_both_inputs_raises(self, tmp_path: Path):
        with pytest.raises(ValueError, match="not both"):
            calibrate(
                csv_path=DATA_DIR / "calibration_with_header.csv",
                data_points=[(0.0, 0.0), (1.0, 1.0)],
                calibration_json_path=tmp_path / "cal.json",
            )

    def test_no_inputs_raises(self, tmp_path: Path):
        with pytest.raises(ValueError, match="either"):
            calibrate(calibration_json_path=tmp_path / "cal.json")


class TestCalibrateFromManual:

    def test_manual_points(self, tmp_path: Path):
        points = [(-0.2, -8.5), (0.1, -3.2), (0.4, 2.1), (0.7, 7.4)]
        json_path = calibrate(
            data_points=points,
            method="linear",
            output_dir=tmp_path / "out",
            calibration_json_path=tmp_path / "cal.json",
        )
        assert json_path.exists()


class TestPredictPotential:

    def test_predict_with_saved_calibration(self, tmp_path: Path):
        # Calibrate first
        json_path = calibrate(
            data_points=[(0.0, 0.0), (1.0, 10.0)],
            method="linear",
            calibration_json_path=tmp_path / "cal.json",
        )

        # Predict
        result = predict_potential(
            5.0,
            calibration_json_path=json_path,
        )
        np.testing.assert_allclose(result, 0.5, atol=1e-10)

    def test_predict_array(self, tmp_path: Path):
        json_path = calibrate(
            data_points=[(0.0, 0.0), (1.0, 10.0)],
            method="linear",
            calibration_json_path=tmp_path / "cal.json",
        )
        result = predict_potential(
            np.array([0.0, 5.0, 10.0]),
            calibration_json_path=json_path,
        )
        np.testing.assert_allclose(result, [0.0, 0.5, 1.0], atol=1e-10)

    def test_predict_with_reference_conversion(self, tmp_path: Path):
        json_path = calibrate(
            data_points=[(0.0, 0.0), (1.0, 10.0)],
            method="linear",
            calibration_json_path=tmp_path / "cal.json",
        )
        result_she = predict_potential(5.0, calibration_json_path=json_path)
        result_rhe = predict_potential(
            5.0,
            calibration_json_path=json_path,
            target_reference="RHE",
            pH=7.0,
        )
        # RHE should be larger than SHE at pH=7
        assert float(result_rhe) > float(result_she)
