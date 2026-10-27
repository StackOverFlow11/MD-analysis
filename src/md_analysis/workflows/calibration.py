"""Charge ↔ potential calibration workflow facade.

Exposes two programmatic entry points covering the σ → φ calibration
surface:

- :func:`run_calibration_fit`     — fit a charge → potential mapper
  from a CSV or in-memory data points and persist it as a JSON
  calibration file (optionally writing a CSV copy + diagnostic PNG).
- :func:`run_calibration_predict` — query an existing calibration:
  predict potential from a given σ value, optionally converting the
  reference scale (SHE / RHE / PZC). No files are written.

Business logic lives in
:mod:`md_analysis.electrochemical.calibration`; this module organises
parameters, output paths, and the :class:`WorkflowResult` contract.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np

from ..electrochemical.calibration.config import (
    DEFAULT_FITTING_METHOD,
    DEFAULT_PH,
    DEFAULT_POLY_DEGREE,
    DEFAULT_TEMPERATURE_K,
)
from .models import WorkflowResult

logger = logging.getLogger(__name__)


def run_calibration_fit(
    *,
    calibration_json_path: Path | str,
    csv_path: Path | str | None = None,
    data_points: list[tuple[float, float]] | None = None,
    method: str = DEFAULT_FITTING_METHOD,
    poly_degree: int = DEFAULT_POLY_DEGREE,
    output_dir: Path | str | None = None,
) -> WorkflowResult:
    """Fit a σ → φ calibration mapper and persist it as JSON.

    Exactly one of ``csv_path`` / ``data_points`` must be provided.
    Input potentials are interpreted as V vs SHE (see
    :mod:`md_analysis.electrochemical.calibration` for unit conventions).

    Artifacts:

    - ``calibration_json`` — always produced; contains the raw data
      points and fitted mapper parameters.
    - ``calibration_csv``  — produced when ``output_dir`` is given.
    - ``calibration_png``  — produced when ``output_dir`` is given.

    Fit-quality metrics (``r_squared``, ``rmse``, ``equation``) are
    reported facts in ``metadata`` — this workflow does not judge
    whether the fit is scientifically "good".
    """
    from ..electrochemical.calibration.CalibrationWorkflow import (
        calibrate_with_report,
    )

    json_path = Path(calibration_json_path)
    csv_p = Path(csv_path) if csv_path is not None else None
    out_dir = Path(output_dir) if output_dir is not None else None

    logger.info(
        "Starting calibration fit: json=%s, output_dir=%s, method=%s",
        json_path,
        out_dir,
        method,
    )

    report = calibrate_with_report(
        csv_path=csv_p,
        data_points=data_points,
        method=method,
        poly_degree=poly_degree,
        output_dir=out_dir,
        calibration_json_path=json_path,
    )

    artifacts: dict[str, Path] = {
        "calibration_json": Path(report.calibration_json),
    }
    if report.calibration_csv is not None:
        artifacts["calibration_csv"] = Path(report.calibration_csv)
    if report.calibration_png is not None:
        artifacts["calibration_png"] = Path(report.calibration_png)

    # The workflow's output_dir identifies where the produced files
    # were written. If the caller provided one we honour it; otherwise
    # fall back to the directory holding the persisted JSON so the
    # WorkflowResult.output_dir always points to an existing folder.
    result_output_dir = out_dir if out_dir is not None else json_path.parent

    metadata: dict[str, Any] = {
        "n_points": report.n_points,
        "reference": report.reference,
        "method": report.method,
        "r_squared": report.r_squared,
        "rmse": report.rmse,
        "equation": report.equation,
        "poly_degree": poly_degree,
        "input_source": "csv" if csv_p is not None else "data_points",
    }
    return WorkflowResult(
        name="calibration_fit",
        output_dir=result_output_dir,
        artifacts=artifacts,
        metadata=metadata,
        extra=report,
    )


def run_calibration_predict(
    sigma: float | list[float] | np.ndarray,
    *,
    calibration_json_path: Path | str,
    target_reference: str | None = None,
    temperature_K: float = DEFAULT_TEMPERATURE_K,
    pH: float = DEFAULT_PH,
    phi_pzc: float | None = None,
    output_dir: Path | str | None = None,
) -> WorkflowResult:
    """Predict electrode potential from surface charge density.

    This is a **read-only query** against an existing calibration JSON;
    no files are written and ``WorkflowResult.artifacts`` is always
    empty. Predicted values live on ``WorkflowResult.extra`` (a
    :class:`CalibrationPredictionReport`) and are also surfaced in
    ``metadata`` for callers that only need scalar/short-array results.

    ``target_reference`` selects the output reference scale (``"SHE"``,
    ``"RHE"``, or ``"PZC"``). When ``None`` the stored reference is
    reused. ``temperature_K`` / ``pH`` / ``phi_pzc`` participate only
    in reference-conversion arithmetic.

    ``output_dir`` is reported on the result for traceability; when
    unset it defaults to the directory holding the calibration JSON.
    """
    from ..electrochemical.calibration.CalibrationWorkflow import (
        predict_potential_with_report,
    )

    json_path = Path(calibration_json_path)
    out_dir = Path(output_dir) if output_dir is not None else json_path.parent

    logger.info(
        "Starting calibration predict: json=%s, target_reference=%s",
        json_path,
        target_reference,
    )

    report = predict_potential_with_report(
        sigma,
        calibration_json_path=json_path,
        target_reference=target_reference,
        temperature_K=temperature_K,
        pH=pH,
        phi_pzc=phi_pzc,
    )

    metadata: dict[str, Any] = {
        "n_values": report.n_values,
        "is_scalar_input": report.is_scalar_input,
        "stored_reference": report.stored_reference,
        "target_reference": report.target_reference,
        "temperature_K": report.temperature_K,
        "pH": report.pH,
        "phi_pzc": report.phi_pzc,
        "calibration_json": str(json_path),
    }
    return WorkflowResult(
        name="calibration_predict",
        output_dir=out_dir,
        artifacts={},
        metadata=metadata,
        extra=report,
    )
