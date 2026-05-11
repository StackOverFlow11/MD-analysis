"""Batch 2: contract-backed calibration_fit_csv / calibration_predict.

Covers the new ``calibrate_with_report`` / ``predict_potential_with_report``
wrappers, their :class:`CalibrationFitReport` / :class:`CalibrationPredictionReport`
dataclasses, and the end-to-end dispatch behaviour through
``md_analysis.agent.dispatch``.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pytest

from md_analysis.agent import dispatch, get_task_schema
from md_analysis.electrochemical.calibration import (
    CalibrationFitReport,
    CalibrationPredictionReport,
    calibrate_with_report,
    predict_potential_with_report,
)

DATA_DIR = Path(__file__).parent / "data"
CSV_WITH_HEADER = DATA_DIR / "calibration_with_header.csv"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def fit_env(tmp_path):
    return {
        "csv": CSV_WITH_HEADER,
        "cal_json": tmp_path / "cal.json",
        "outdir": tmp_path / "out",
        "tmp": tmp_path,
    }


@pytest.fixture
def fitted_cal(tmp_path):
    """Run a real linear fit so dispatch tests have a loadable JSON."""
    cal_json = tmp_path / "fitted_cal.json"
    calibrate_with_report(
        csv_path=CSV_WITH_HEADER,
        method="linear",
        calibration_json_path=cal_json,
    )
    return cal_json


# ---------------------------------------------------------------------------
# Wrapper level — calibrate_with_report
# ---------------------------------------------------------------------------


class TestCalibrateWithReport:
    def test_csv_success_writes_all_artifacts(self, fit_env):
        r = calibrate_with_report(
            csv_path=fit_env["csv"],
            method="linear",
            output_dir=fit_env["outdir"],
            calibration_json_path=fit_env["cal_json"],
        )
        assert isinstance(r, CalibrationFitReport)
        assert r.calibration_json.is_file()
        assert r.calibration_csv is not None and r.calibration_csv.is_file()
        assert r.calibration_png is not None and r.calibration_png.is_file()
        assert r.n_points == 4
        assert r.method == "linear"
        assert r.reference == "SHE"
        # r_squared should be a native Python float
        assert isinstance(r.r_squared, float)

    def test_csv_success_no_output_dir_skips_csv_png(self, fit_env):
        r = calibrate_with_report(
            csv_path=fit_env["csv"],
            method="linear",
            calibration_json_path=fit_env["cal_json"],
        )
        assert r.calibration_json.is_file()
        assert r.calibration_csv is None
        assert r.calibration_png is None

    def test_manual_data_points_success(self, fit_env):
        points = [(-0.2, -8.5), (0.1, -3.2), (0.4, 2.1), (0.7, 7.4)]
        r = calibrate_with_report(
            data_points=points,
            method="linear",
            calibration_json_path=fit_env["cal_json"],
        )
        assert r.n_points == 4
        assert r.calibration_json.is_file()

    def test_missing_json_path_raises_value_error(self, fit_env):
        with pytest.raises(ValueError, match="calibration_json_path"):
            calibrate_with_report(csv_path=fit_env["csv"])

    def test_both_inputs_raises_value_error(self, fit_env):
        with pytest.raises(ValueError):
            calibrate_with_report(
                csv_path=fit_env["csv"],
                data_points=[(0.0, 0.0), (1.0, 1.0)],
                calibration_json_path=fit_env["cal_json"],
            )

    def test_neither_input_raises_value_error(self, fit_env):
        with pytest.raises(ValueError):
            calibrate_with_report(
                calibration_json_path=fit_env["cal_json"],
            )

    def test_fit_params_exposes_full_dict(self, fit_env):
        r = calibrate_with_report(
            csv_path=fit_env["csv"],
            method="polynomial",
            poly_degree=2,
            calibration_json_path=fit_env["cal_json"],
        )
        # Polynomial fit stores coefficients — make sure we didn't drop them
        assert isinstance(r.fit_params, dict)
        assert "r_squared" in r.fit_params
        assert "rmse" in r.fit_params
        assert "equation" in r.fit_params

    def test_report_to_dict_is_json_serializable(self, fit_env):
        r = calibrate_with_report(
            csv_path=fit_env["csv"],
            method="linear",
            output_dir=fit_env["outdir"],
            calibration_json_path=fit_env["cal_json"],
        )
        json.dumps(r.to_dict())


# ---------------------------------------------------------------------------
# Wrapper level — predict_potential_with_report
# ---------------------------------------------------------------------------


class TestPredictPotentialWithReport:
    def test_scalar_input(self, fitted_cal):
        r = predict_potential_with_report(
            5.0, calibration_json_path=fitted_cal,
        )
        assert isinstance(r, CalibrationPredictionReport)
        assert r.is_scalar_input is True
        assert r.n_values == 1
        assert len(r.potential_V) == 1
        assert isinstance(r.potential_V[0], float)

    def test_array_input(self, fitted_cal):
        r = predict_potential_with_report(
            [1.0, 3.0, 5.0], calibration_json_path=fitted_cal,
        )
        assert r.is_scalar_input is False
        assert r.n_values == 3
        assert len(r.potential_V) == 3

    def test_default_target_reference_resolves_to_stored(self, fitted_cal):
        r = predict_potential_with_report(
            5.0, calibration_json_path=fitted_cal,
        )
        # Stored reference is SHE in our fixture
        assert r.stored_reference == "SHE"
        assert r.target_reference == "SHE"

    def test_target_reference_rhe(self, fitted_cal):
        r = predict_potential_with_report(
            5.0, calibration_json_path=fitted_cal,
            target_reference="RHE", pH=7.0,
        )
        assert r.target_reference == "RHE"

    def test_pzc_without_phi_pzc_raises(self, fitted_cal):
        with pytest.raises(ValueError):
            predict_potential_with_report(
                5.0, calibration_json_path=fitted_cal,
                target_reference="PZC",
            )

    def test_missing_json_path_raises(self, fitted_cal):
        with pytest.raises(ValueError, match="calibration_json_path"):
            predict_potential_with_report(5.0)

    def test_missing_file_raises_file_not_found(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            predict_potential_with_report(
                5.0, calibration_json_path=tmp_path / "nope.json",
            )

    def test_malformed_json_raises_json_decode_error(self, tmp_path):
        import json as _json

        bad = tmp_path / "bad.json"
        bad.write_text("{not valid json", encoding="utf-8")
        # JSONDecodeError is the exact exception raised by json.loads;
        # narrowing here prevents silent drift if some future refactor
        # starts raising something unrelated.
        with pytest.raises(_json.JSONDecodeError):
            predict_potential_with_report(
                5.0, calibration_json_path=bad,
            )

    def test_lowercase_target_reference_normalized_to_uppercase(
        self, fitted_cal,
    ):
        """Contract enum is uppercase; the wrapper must not echo back
        a lowercase ``target_reference`` even though
        ``convert_reference()`` is case-insensitive internally."""
        r_lower = predict_potential_with_report(
            5.0, calibration_json_path=fitted_cal,
            target_reference="rhe", pH=7.0,
        )
        assert r_lower.target_reference == "RHE"

        # Compare potential_V against the upper-case call to lock the
        # invariant that normalisation does NOT change the predicted value.
        r_upper = predict_potential_with_report(
            5.0, calibration_json_path=fitted_cal,
            target_reference="RHE", pH=7.0,
        )
        assert r_lower.potential_V == r_upper.potential_V

    def test_report_to_dict_is_json_serializable(self, fitted_cal):
        r = predict_potential_with_report(
            [1.0, 2.0], calibration_json_path=fitted_cal,
        )
        json.dumps(r.to_dict())


# ---------------------------------------------------------------------------
# Contract / dispatch — calibration_fit_csv
# ---------------------------------------------------------------------------


class TestCalibrationFitCsvSchema:
    def test_schema_keys(self):
        s = get_task_schema("calibration_fit_csv")
        assert set(s.keys()) == {"name", "description", "parameters"}

    def test_required_field_is_calibration_json_path(self):
        s = get_task_schema("calibration_fit_csv")
        assert s["parameters"]["required"] == ["calibration_json_path"]

    def test_includes_exactly_expected_fields(self):
        s = get_task_schema("calibration_fit_csv")
        assert set(s["parameters"]["properties"].keys()) == {
            "csv_path", "data_points", "method", "poly_degree",
            "output_dir", "calibration_json_path",
        }

    def test_method_enum_has_all_four_methods(self):
        s = get_task_schema("calibration_fit_csv")
        assert set(
            s["parameters"]["properties"]["method"].get("enum", [])
        ) == {"linear", "polynomial", "spline", "differential_capacitance"}

    def test_is_contract_backed(self):
        from md_analysis.agent._core import _TASK_REGISTRY
        assert _TASK_REGISTRY["calibration_fit_csv"].contract is not None


class TestCalibrationFitCsvDispatch:
    def test_csv_success_returns_all_artifacts(self, fit_env):
        r = dispatch("calibration_fit_csv", {
            "csv_path": str(fit_env["csv"]),
            "method": "linear",
            "output_dir": str(fit_env["outdir"]),
            "calibration_json_path": str(fit_env["cal_json"]),
        })
        assert r.success, r.errors
        assert "calibration_json" in r.outputs
        assert "calibration_csv" in r.outputs
        assert "calibration_png" in r.outputs
        for k in (
            "n_points", "method", "r_squared", "rmse", "equation",
            "fit_params", "reference",
        ):
            assert k in r.summary
        # summary must be JSON-serializable
        json.dumps(r.summary)

    def test_missing_json_path_returns_validation(self, fit_env):
        r = dispatch("calibration_fit_csv", {
            "csv_path": str(fit_env["csv"]),
            "method": "linear",
        })
        assert not r.success
        assert r.error_type == "validation"

    def test_both_inputs_returns_validation(self, fit_env):
        r = dispatch("calibration_fit_csv", {
            "csv_path": str(fit_env["csv"]),
            "data_points": [[0.0, 0.0], [1.0, 1.0]],
            "calibration_json_path": str(fit_env["cal_json"]),
        })
        assert not r.success
        assert r.error_type == "validation"

    def test_neither_input_returns_validation(self, fit_env):
        r = dispatch("calibration_fit_csv", {
            "calibration_json_path": str(fit_env["cal_json"]),
        })
        assert not r.success
        assert r.error_type == "validation"

    def test_missing_csv_returns_file_not_found(self, fit_env):
        r = dispatch("calibration_fit_csv", {
            "csv_path": str(fit_env["tmp"] / "missing.csv"),
            "calibration_json_path": str(fit_env["cal_json"]),
        })
        assert not r.success
        assert r.error_type == "file_not_found"

    def test_invalid_method_returns_validation(self, fit_env):
        r = dispatch("calibration_fit_csv", {
            "csv_path": str(fit_env["csv"]),
            "method": "bogus",
            "calibration_json_path": str(fit_env["cal_json"]),
        })
        assert not r.success
        assert r.error_type == "validation"

    def test_malformed_csv_returns_validation(self, fit_env):
        bad_csv = fit_env["tmp"] / "bad.csv"
        bad_csv.write_text("not,numeric,at,all\nalso,bad", encoding="utf-8")
        r = dispatch("calibration_fit_csv", {
            "csv_path": str(bad_csv),
            "method": "linear",
            "calibration_json_path": str(fit_env["cal_json"]),
        })
        assert not r.success
        assert r.error_type == "validation"

    def test_manual_data_points_success(self, fit_env):
        r = dispatch("calibration_fit_csv", {
            "data_points": [[-0.2, -8.5], [0.1, -3.2], [0.4, 2.1], [0.7, 7.4]],
            "method": "linear",
            "calibration_json_path": str(fit_env["cal_json"]),
        })
        assert r.success, r.errors
        assert r.summary["n_points"] == 4


# ---------------------------------------------------------------------------
# Contract / dispatch — calibration_predict
# ---------------------------------------------------------------------------


class TestCalibrationPredictSchema:
    def test_schema_keys(self):
        s = get_task_schema("calibration_predict")
        assert set(s.keys()) == {"name", "description", "parameters"}

    def test_required_fields(self):
        s = get_task_schema("calibration_predict")
        assert set(s["parameters"]["required"]) == {
            "sigma", "calibration_json_path",
        }

    def test_includes_exactly_expected_fields(self):
        s = get_task_schema("calibration_predict")
        assert set(s["parameters"]["properties"].keys()) == {
            "sigma", "calibration_json_path", "target_reference",
            "temperature_K", "pH", "phi_pzc",
        }

    def test_target_reference_enum_includes_uppercase_and_null(self):
        s = get_task_schema("calibration_predict")
        enum_vals = s["parameters"]["properties"]["target_reference"].get("enum")
        assert enum_vals is not None
        assert set(enum_vals) == {"SHE", "RHE", "PZC", None}

    def test_is_contract_backed(self):
        from md_analysis.agent._core import _TASK_REGISTRY
        assert _TASK_REGISTRY["calibration_predict"].contract is not None


class TestCalibrationPredictDispatch:
    def test_scalar_success(self, fitted_cal):
        r = dispatch("calibration_predict", {
            "sigma": 5.0,
            "calibration_json_path": str(fitted_cal),
        })
        assert r.success, r.errors
        assert r.summary["is_scalar_input"] is True
        assert r.summary["n_values"] == 1
        assert len(r.summary["potential_V"]) == 1
        # No artifacts
        assert r.outputs == {}
        # summary is JSON-serializable
        json.dumps(r.summary)

    def test_array_success(self, fitted_cal):
        r = dispatch("calibration_predict", {
            "sigma": [1.0, 3.0, 5.0],
            "calibration_json_path": str(fitted_cal),
        })
        assert r.success, r.errors
        assert r.summary["is_scalar_input"] is False
        assert r.summary["n_values"] == 3
        assert len(r.summary["potential_V"]) == 3

    def test_default_target_reference_resolves(self, fitted_cal):
        r = dispatch("calibration_predict", {
            "sigma": 5.0,
            "calibration_json_path": str(fitted_cal),
        })
        assert r.success
        # Resolved, not raw None — stored reference is 'SHE' in fixture
        assert r.summary["target_reference"] == "SHE"
        assert r.summary["stored_reference"] == "SHE"

    def test_rhe_conversion_reports_rhe(self, fitted_cal):
        r = dispatch("calibration_predict", {
            "sigma": 5.0,
            "calibration_json_path": str(fitted_cal),
            "target_reference": "RHE",
            "pH": 7.0,
        })
        assert r.success
        assert r.summary["target_reference"] == "RHE"

    def test_lowercase_target_reference_reported_as_uppercase(
        self, fitted_cal,
    ):
        """dispatch-level end-to-end guard: even if the caller sends a
        lowercase ``target_reference``, the summary must echo back a
        value inside the contract enum (uppercase)."""
        r = dispatch("calibration_predict", {
            "sigma": 5.0,
            "calibration_json_path": str(fitted_cal),
            "target_reference": "rhe",
            "pH": 7.0,
        })
        assert r.success, r.errors
        assert r.summary["target_reference"] == "RHE"

    def test_pzc_without_phi_pzc_returns_validation(self, fitted_cal):
        r = dispatch("calibration_predict", {
            "sigma": 5.0,
            "calibration_json_path": str(fitted_cal),
            "target_reference": "PZC",
        })
        assert not r.success
        assert r.error_type == "validation"

    def test_missing_json_path_returns_validation(self):
        r = dispatch("calibration_predict", {"sigma": 5.0})
        assert not r.success
        assert r.error_type == "validation"

    def test_missing_calibration_file_returns_file_not_found(self, tmp_path):
        r = dispatch("calibration_predict", {
            "sigma": 5.0,
            "calibration_json_path": str(tmp_path / "nope.json"),
        })
        assert not r.success
        assert r.error_type == "file_not_found"

    def test_malformed_json_returns_validation(self, tmp_path):
        bad = tmp_path / "bad.json"
        bad.write_text("{not valid json", encoding="utf-8")
        r = dispatch("calibration_predict", {
            "sigma": 5.0,
            "calibration_json_path": str(bad),
        })
        assert not r.success
        assert r.error_type == "validation"
