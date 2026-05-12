"""Unit tests for ``md_analysis.workflows.calibration``.

Covers the workflow facade for calibration fit / predict:

- ``WorkflowResult`` shape and artifact bookkeeping
- ``metadata`` carries fit metrics and reference info
- ``extra`` carries the strongly-typed report
- ``run_calibration_predict`` produces no artifacts (read-only query)
- Optional / required output_dir handling
- Error propagation when required inputs are missing
- ``require_artifacts_exist`` accepts a successful fit
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from md_analysis.electrochemical.calibration.CalibrationWorkflow import (
    CalibrationFitReport,
    CalibrationPredictionReport,
)
from md_analysis.workflows import (
    MissingArtifactError,
    WorkflowResult,
    require_artifacts_exist,
    run_calibration_fit,
    run_calibration_predict,
)

DATA_DIR = Path(__file__).resolve().parents[1] / "calibration" / "data"
CSV_WITH_HEADER = DATA_DIR / "calibration_with_header.csv"


# ---------------------------------------------------------------------------
# run_calibration_fit — happy paths
# ---------------------------------------------------------------------------


class TestRunCalibrationFitFromCSV:
    def test_returns_workflow_result_with_artifacts(self, tmp_path: Path) -> None:
        out_dir = tmp_path / "fit"
        json_path = tmp_path / "cal.json"
        result = run_calibration_fit(
            calibration_json_path=json_path,
            csv_path=CSV_WITH_HEADER,
            method="linear",
            output_dir=out_dir,
        )
        assert isinstance(result, WorkflowResult)
        assert result.name == "calibration_fit"
        assert result.output_dir == out_dir
        # All three artifacts should be present and on disk
        for key in ("calibration_json", "calibration_csv", "calibration_png"):
            assert key in result.artifacts
            assert result.artifacts[key].is_file()
        require_artifacts_exist(result)  # should not raise

    def test_metadata_carries_fit_metrics(self, tmp_path: Path) -> None:
        result = run_calibration_fit(
            calibration_json_path=tmp_path / "cal.json",
            csv_path=CSV_WITH_HEADER,
            method="linear",
            output_dir=tmp_path / "fit",
        )
        meta = result.metadata
        assert meta["input_source"] == "csv"
        assert meta["method"] == "linear"
        assert meta["n_points"] >= 2
        assert isinstance(meta["r_squared"], float)
        assert isinstance(meta["rmse"], float)
        assert isinstance(meta["equation"], str)
        assert meta["reference"]  # non-empty string

    def test_extra_is_calibration_fit_report(self, tmp_path: Path) -> None:
        result = run_calibration_fit(
            calibration_json_path=tmp_path / "cal.json",
            csv_path=CSV_WITH_HEADER,
            method="linear",
            output_dir=tmp_path / "fit",
        )
        assert isinstance(result.extra, CalibrationFitReport)
        assert Path(result.extra.calibration_json).is_file()


class TestRunCalibrationFitFromDataPoints:
    def test_manual_points_no_output_dir(self, tmp_path: Path) -> None:
        json_path = tmp_path / "cal.json"
        result = run_calibration_fit(
            calibration_json_path=json_path,
            data_points=[(0.0, 0.0), (1.0, 10.0)],
            method="linear",
        )
        # Without output_dir we still get the JSON, but no CSV/PNG
        assert "calibration_json" in result.artifacts
        assert "calibration_csv" not in result.artifacts
        assert "calibration_png" not in result.artifacts
        # Result output_dir falls back to the JSON's parent
        assert result.output_dir == json_path.parent
        assert result.metadata["input_source"] == "data_points"
        require_artifacts_exist(result)

    def test_polynomial_method(self, tmp_path: Path) -> None:
        result = run_calibration_fit(
            calibration_json_path=tmp_path / "cal.json",
            data_points=[(-0.2, -8.5), (0.1, -3.2), (0.4, 2.1), (0.7, 7.4)],
            method="polynomial",
            poly_degree=2,
            output_dir=tmp_path / "fit",
        )
        assert result.metadata["method"] == "polynomial"
        assert result.metadata["poly_degree"] == 2


# ---------------------------------------------------------------------------
# run_calibration_fit — error propagation
# ---------------------------------------------------------------------------


class TestRunCalibrationFitErrors:
    def test_both_csv_and_data_points_raises(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="not both"):
            run_calibration_fit(
                calibration_json_path=tmp_path / "cal.json",
                csv_path=CSV_WITH_HEADER,
                data_points=[(0.0, 0.0), (1.0, 1.0)],
            )

    def test_neither_csv_nor_data_points_raises(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="either"):
            run_calibration_fit(calibration_json_path=tmp_path / "cal.json")


# ---------------------------------------------------------------------------
# run_calibration_predict
# ---------------------------------------------------------------------------


class TestRunCalibrationPredict:
    @pytest.fixture
    def linear_cal_json(self, tmp_path: Path) -> Path:
        json_path = tmp_path / "cal.json"
        run_calibration_fit(
            calibration_json_path=json_path,
            data_points=[(0.0, 0.0), (1.0, 10.0)],
            method="linear",
        )
        return json_path

    def test_scalar_returns_no_artifacts(self, linear_cal_json: Path) -> None:
        result = run_calibration_predict(
            5.0,
            calibration_json_path=linear_cal_json,
        )
        assert isinstance(result, WorkflowResult)
        assert result.name == "calibration_predict"
        assert result.artifacts == {}
        # Predict is read-only; require_artifacts_exist passes trivially
        require_artifacts_exist(result)

    def test_scalar_metadata_and_extra(self, linear_cal_json: Path) -> None:
        result = run_calibration_predict(
            5.0,
            calibration_json_path=linear_cal_json,
        )
        meta = result.metadata
        assert meta["is_scalar_input"] is True
        assert meta["n_values"] == 1
        assert meta["target_reference"] == meta["stored_reference"]
        assert meta["calibration_json"] == str(linear_cal_json)
        assert isinstance(result.extra, CalibrationPredictionReport)
        np.testing.assert_allclose(result.extra.potential_V, (0.5,), atol=1e-10)

    def test_array_input(self, linear_cal_json: Path) -> None:
        result = run_calibration_predict(
            [0.0, 5.0, 10.0],
            calibration_json_path=linear_cal_json,
        )
        assert result.metadata["is_scalar_input"] is False
        assert result.metadata["n_values"] == 3
        report = result.extra
        assert isinstance(report, CalibrationPredictionReport)
        np.testing.assert_allclose(report.potential_V, (0.0, 0.5, 1.0), atol=1e-10)

    def test_target_reference_normalized_upper(
        self, linear_cal_json: Path
    ) -> None:
        result = run_calibration_predict(
            5.0,
            calibration_json_path=linear_cal_json,
            target_reference="rhe",
            pH=7.0,
            temperature_K=298.15,
        )
        assert result.metadata["target_reference"] == "RHE"

    def test_default_output_dir_is_json_parent(
        self, linear_cal_json: Path
    ) -> None:
        result = run_calibration_predict(
            5.0, calibration_json_path=linear_cal_json
        )
        assert result.output_dir == linear_cal_json.parent

    def test_explicit_output_dir_is_preserved(
        self, linear_cal_json: Path, tmp_path: Path
    ) -> None:
        custom = tmp_path / "predict_out"
        # Note: predict does not create the directory; the workflow
        # only records the path for traceability.
        result = run_calibration_predict(
            5.0,
            calibration_json_path=linear_cal_json,
            output_dir=custom,
        )
        assert result.output_dir == custom


# ---------------------------------------------------------------------------
# to_dict serialisation round-trip
# ---------------------------------------------------------------------------


class TestWorkflowResultSerialisation:
    def test_fit_to_dict_is_json_friendly(self, tmp_path: Path) -> None:
        import json

        result = run_calibration_fit(
            calibration_json_path=tmp_path / "cal.json",
            data_points=[(0.0, 0.0), (1.0, 10.0)],
            method="linear",
        )
        # extra has its own to_dict() — exercise the WorkflowResult path
        dumped = result.to_dict()
        json.dumps(dumped)  # must not raise

    def test_predict_to_dict_is_json_friendly(self, tmp_path: Path) -> None:
        import json

        json_path = tmp_path / "cal.json"
        run_calibration_fit(
            calibration_json_path=json_path,
            data_points=[(0.0, 0.0), (1.0, 10.0)],
            method="linear",
        )
        result = run_calibration_predict(
            [0.0, 5.0, 10.0],
            calibration_json_path=json_path,
        )
        json.dumps(result.to_dict())


# ---------------------------------------------------------------------------
# Negative — require_artifacts_exist detects fake file
# ---------------------------------------------------------------------------


def test_require_artifacts_exist_flags_missing_after_unlink(
    tmp_path: Path,
) -> None:
    """If a fit's artifact is deleted out-of-band, the validator catches it."""
    result = run_calibration_fit(
        calibration_json_path=tmp_path / "cal.json",
        data_points=[(0.0, 0.0), (1.0, 10.0)],
        method="linear",
        output_dir=tmp_path / "fit",
    )
    # Remove the PNG (if produced)
    if "calibration_png" in result.artifacts:
        result.artifacts["calibration_png"].unlink()
        with pytest.raises(MissingArtifactError):
            require_artifacts_exist(result)
