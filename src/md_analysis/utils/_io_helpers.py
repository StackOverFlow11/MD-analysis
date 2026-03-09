"""Private shared I/O and numerical helpers."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

import numpy as np


def _cumulative_average(values: np.ndarray) -> np.ndarray:
    """Element-wise cumulative average of a 1-D array."""
    csum = np.cumsum(values, dtype=float)
    return csum / np.arange(1, values.size + 1, dtype=float)


def _write_csv(path: Path, rows: Iterable[dict], fieldnames: list[str]) -> None:
    """Write *rows* (list of dicts) to a CSV file, creating parent dirs."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
