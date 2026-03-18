"""High-level calibration workflow functions.

Public API
----------
- ``calibrate``           — load data, fit, save JSON, plot
- ``predict_potential``   — given σ, return φ using saved calibration
- ``convert_reference``   — convert between SHE, RHE, PZC
"""

from __future__ import annotations

import csv as _csv
import logging
from pathlib import Path

import numpy as np

from .config import (
    DEFAULT_CALIBRATION_CSV_NAME,
    DEFAULT_CALIBRATION_FILE,
    DEFAULT_FITTING_METHOD,
    DEFAULT_PH,
    DEFAULT_POLY_DEGREE,
    DEFAULT_REFERENCE,
    DEFAULT_TEMPERATURE_K,
    F_C_PER_MOL,
    LN10,
    R_J_PER_MOL_K,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Reference conversion
# ---------------------------------------------------------------------------

def convert_reference(
    phi: float | np.ndarray,
    *,
    from_ref: str,
    to_ref: str,
    temperature_K: float = DEFAULT_TEMPERATURE_K,
    pH: float = DEFAULT_PH,
    phi_pzc: float | None = None,
) -> np.ndarray:
    """Convert potential between reference scales.

    Supported scales: ``"SHE"``, ``"RHE"``, ``"PZC"``.

    Conversions go through SHE as hub:

    - SHE → RHE:  ``φ_RHE = φ_SHE + (RT/F) · ln(10) · pH``
    - SHE → PZC:  ``φ_PZC = φ_SHE − φ_pzc``

    Parameters
    ----------
    phi : float or array
        Input potential value(s).
    from_ref, to_ref : str
        Source / target reference scale.
    temperature_K : float
        Temperature in K (only needed for RHE).
    pH : float
        Solution pH (only needed for RHE).
    phi_pzc : float, optional
        Potential of zero charge in V vs SHE (required for PZC conversions).

    Returns
    -------
    np.ndarray
    """
    phi_arr = np.asarray(phi, dtype=float)
    from_ref = from_ref.upper()
    to_ref = to_ref.upper()

    if from_ref == to_ref:
        return phi_arr.copy()

    _VALID = {"SHE", "RHE", "PZC"}
    if from_ref not in _VALID:
        raise ValueError(f"Unknown reference {from_ref!r}; choose from {_VALID}")
    if to_ref not in _VALID:
        raise ValueError(f"Unknown reference {to_ref!r}; choose from {_VALID}")

    rhe_shift = R_J_PER_MOL_K * temperature_K / F_C_PER_MOL * LN10 * pH

    # Convert to SHE first
    if from_ref == "SHE":
        phi_she = phi_arr
    elif from_ref == "RHE":
        phi_she = phi_arr - rhe_shift
    elif from_ref == "PZC":
        if phi_pzc is None:
            raise ValueError("phi_pzc is required for PZC → SHE conversion")
        phi_she = phi_arr + phi_pzc

    # Convert SHE to target
    if to_ref == "SHE":
        return phi_she
    elif to_ref == "RHE":
        return phi_she + rhe_shift
    elif to_ref == "PZC":
        if phi_pzc is None:
            raise ValueError("phi_pzc is required for SHE → PZC conversion")
        return phi_she - phi_pzc

    return phi_she  # pragma: no cover


# ---------------------------------------------------------------------------
# Calibrate
# ---------------------------------------------------------------------------

def calibrate(
    csv_path: str | Path | None = None,
    *,
    data_points: list[tuple[float, float]] | None = None,
    method: str = DEFAULT_FITTING_METHOD,
    poly_degree: int = DEFAULT_POLY_DEGREE,
    output_dir: Path | None = None,
    calibration_json_path: Path | None = None,
) -> Path:
    """Run the full calibration workflow.

    1. Load data from *csv_path* **or** *data_points* (exactly one).
    2. Fit the σ → φ relationship.
    3. Save calibration JSON (raw data + fit parameters).
    4. Write a CSV copy of the data points.
    5. Plot the fit.

    Input potentials are always in V vs SHE.

    Parameters
    ----------
    csv_path : path, optional
        CSV file (col 1 = φ vs SHE, col 2 = σ). Auto-detects header.
    data_points : list of (φ, σ) tuples, optional
        Manual data input. Mutually exclusive with *csv_path*.
    method : str
        ``"linear"`` | ``"polynomial"`` | ``"spline"``.
    poly_degree : int
        Polynomial degree (only for ``method="polynomial"``).
    output_dir : Path, optional
        Where to write CSV + PNG output.
    calibration_json_path : Path, optional
        Where to save the calibration JSON.

    Returns
    -------
    Path
        Path to the saved calibration JSON file.
    """
    from ._data import CalibrationData, load_calibration_csv, save_calibration_json
    from ._mapper import FittingInfo, create_mapper
    from ._plot import plot_calibration

    # --- 1. Load data ---
    if csv_path is not None and data_points is not None:
        raise ValueError("Provide csv_path or data_points, not both")
    if csv_path is None and data_points is None:
        raise ValueError("Provide either csv_path or data_points")

    if csv_path is not None:
        cal = load_calibration_csv(csv_path)
    else:
        potentials = np.array([p[0] for p in data_points])
        charge_densities = np.array([p[1] for p in data_points])
        cal = CalibrationData(
            potentials=potentials,
            charge_densities=charge_densities,
        )

    logger.info(
        "Loaded %d calibration points (reference=%s)", cal.n_points, cal.reference
    )

    # --- 2. Fit ---
    kwargs = {}
    if method == "polynomial":
        kwargs["degree"] = poly_degree
    mapper = create_mapper(method, **kwargs)
    info = mapper.fit(cal.charge_densities, cal.potentials)
    logger.info("Fit complete: R²=%.6f, RMSE=%.4e", info.r_squared, info.rmse)

    # --- 3. Save JSON ---
    fit_params = mapper.to_dict()
    fit_params["r_squared"] = info.r_squared
    fit_params["rmse"] = info.rmse
    fit_params["equation"] = info.equation_str

    json_path = save_calibration_json(cal, fit_params, calibration_json_path)

    # --- 4. Write CSV copy ---
    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        csv_out = out / DEFAULT_CALIBRATION_CSV_NAME
        with csv_out.open("w", newline="", encoding="utf-8") as fh:
            writer = _csv.writer(fh)
            writer.writerow(["potential_V", "sigma_uC_cm2"])
            for phi_val, sig_val in zip(
                cal.potentials.tolist(), cal.charge_densities.tolist()
            ):
                writer.writerow([phi_val, sig_val])
        logger.info("Calibration CSV saved to %s", csv_out)

    # --- 5. Plot ---
    if output_dir is not None:
        plot_calibration(
            cal.charge_densities,
            cal.potentials,
            mapper,
            info,
            reference=cal.reference,
            output_dir=Path(output_dir),
        )

    return json_path


# ---------------------------------------------------------------------------
# Predict
# ---------------------------------------------------------------------------

def predict_potential(
    sigma: float | np.ndarray,
    *,
    calibration_json_path: Path | None = None,
    target_reference: str | None = None,
    temperature_K: float = DEFAULT_TEMPERATURE_K,
    pH: float = DEFAULT_PH,
    phi_pzc: float | None = None,
) -> np.ndarray:
    """Predict electrode potential from surface charge density.

    Loads the calibration from *calibration_json_path* (or the global
    default), reconstructs the fitted mapper, and predicts φ.

    Parameters
    ----------
    sigma : float or array
        Surface charge density in μC/cm².
    calibration_json_path : Path, optional
        Path to calibration JSON.
    target_reference : str, optional
        Convert result to this reference. If ``None``, uses the
        reference stored in the calibration file.
    temperature_K : float
        Temperature for RHE conversion.
    pH : float
        pH for RHE conversion.
    phi_pzc : float, optional
        Potential of zero charge (V vs SHE); required for PZC output.

    Returns
    -------
    np.ndarray
        Predicted potential in V (vs *target_reference*).
    """
    from ._data import load_calibration_json
    from ._mapper import mapper_from_dict

    cal, fit_params = load_calibration_json(calibration_json_path)
    mapper = mapper_from_dict(fit_params)
    phi = mapper.predict(sigma)

    stored_ref = cal.reference
    out_ref = target_reference or stored_ref
    if out_ref.upper() != stored_ref.upper():
        phi = convert_reference(
            phi,
            from_ref=stored_ref,
            to_ref=out_ref,
            temperature_K=temperature_K,
            pH=pH,
            phi_pzc=phi_pzc,
        )

    return np.asarray(phi, dtype=float)
