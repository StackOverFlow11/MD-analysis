"""Calibration data storage and I/O."""

from __future__ import annotations

import csv
import json
import logging
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from .config import DEFAULT_CALIBRATION_FILE

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data class
# ---------------------------------------------------------------------------

@dataclass
class CalibrationData:
    """Raw calibration data: paired (φ, σ) measurements.

    Attributes
    ----------
    potentials : np.ndarray
        Shape ``(n,)``, electrode potential values in V vs SHE.
    charge_densities : np.ndarray
        Shape ``(n,)``, surface charge densities in μC/cm².
    reference : str
        Always ``"SHE"`` — input potentials must be vs SHE.
    metadata : dict
        Optional metadata (system name, temperature, pH, notes, …).
    """

    potentials: np.ndarray
    charge_densities: np.ndarray
    reference: str = "SHE"  # input always vs SHE
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.potentials = np.asarray(self.potentials, dtype=float)
        self.charge_densities = np.asarray(self.charge_densities, dtype=float)
        if len(self.potentials) != len(self.charge_densities):
            raise ValueError(
                f"potentials ({len(self.potentials)}) and charge_densities "
                f"({len(self.charge_densities)}) must have equal length"
            )
        if len(self.potentials) < 2:
            raise ValueError("At least 2 calibration points are required")

    @property
    def n_points(self) -> int:
        return len(self.potentials)


# ---------------------------------------------------------------------------
# CSV I/O
# ---------------------------------------------------------------------------

def _detect_header(first_row: list[str]) -> bool:
    """Return True if *first_row* looks like a header (non-numeric)."""
    try:
        for cell in first_row[:2]:
            float(cell.strip())
        return False
    except (ValueError, IndexError):
        return True


def load_calibration_csv(csv_path: str | Path) -> CalibrationData:
    """Load calibration data from a two-column CSV file.

    Expected format: column 1 = potential (V), column 2 = σ (μC/cm²).
    The first row is auto-detected as header or data.

    Parameters
    ----------
    csv_path : path
        Path to the CSV file.

    Returns
    -------
    CalibrationData
    """
    path = Path(csv_path)
    if not path.is_file():
        raise FileNotFoundError(f"Calibration CSV not found: {path}")

    with path.open(encoding="utf-8") as fh:
        reader = csv.reader(fh)
        rows = list(reader)

    if not rows:
        raise ValueError(f"Empty CSV file: {path}")

    start = 1 if _detect_header(rows[0]) else 0
    potentials: list[float] = []
    charge_densities: list[float] = []

    for i, row in enumerate(rows[start:], start=start):
        if len(row) < 2:
            raise ValueError(f"Row {i + 1} has fewer than 2 columns: {row}")
        potentials.append(float(row[0].strip()))
        charge_densities.append(float(row[1].strip()))

    return CalibrationData(
        potentials=np.array(potentials),
        charge_densities=np.array(charge_densities),
    )


# ---------------------------------------------------------------------------
# JSON I/O
# ---------------------------------------------------------------------------

def save_calibration_json(
    data: CalibrationData,
    fit_params: dict[str, Any],
    json_path: str | Path | None = None,
) -> Path:
    """Persist calibration data **and** fitted parameters to JSON.

    Both raw data points and fit results are saved so that re-fitting
    with a different method does not require re-entering the data.

    Parameters
    ----------
    data : CalibrationData
        Raw calibration measurements.
    fit_params : dict
        Serialised mapper state (from ``mapper.to_dict()``) plus fit
        quality metrics (``r_squared``, ``rmse``, ``equation``).
    json_path : path, optional
        Defaults to :data:`DEFAULT_CALIBRATION_FILE`.

    Returns
    -------
    Path
        The path where the JSON was written.
    """
    path = Path(json_path) if json_path is not None else DEFAULT_CALIBRATION_FILE
    payload: dict[str, Any] = {
        "version": 1,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "reference": data.reference,
        "metadata": data.metadata,
        "data": {
            "potentials_V": data.potentials.tolist(),
            "charge_densities_uC_cm2": data.charge_densities.tolist(),
        },
        "fit": fit_params,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    logger.info("Calibration saved to %s", path)
    return path


def load_calibration_json(
    json_path: str | Path | None = None,
) -> tuple[CalibrationData, dict[str, Any]]:
    """Load calibration data and fitted parameters from JSON.

    Parameters
    ----------
    json_path : path, optional
        Defaults to :data:`DEFAULT_CALIBRATION_FILE`.

    Returns
    -------
    (CalibrationData, fit_params)
    """
    path = Path(json_path) if json_path is not None else DEFAULT_CALIBRATION_FILE
    if not path.is_file():
        raise FileNotFoundError(f"Calibration JSON not found: {path}")

    payload = json.loads(path.read_text(encoding="utf-8"))

    data_section = payload["data"]
    cal = CalibrationData(
        potentials=np.array(data_section["potentials_V"]),
        charge_densities=np.array(data_section["charge_densities_uC_cm2"]),
        reference=payload.get("reference", "SHE"),
        metadata=payload.get("metadata", {}),
    )

    return cal, payload.get("fit", {})
