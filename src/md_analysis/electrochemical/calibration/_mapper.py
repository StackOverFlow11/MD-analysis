"""Charge-to-potential mapping: ABC and concrete implementations."""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Fitting quality container
# ---------------------------------------------------------------------------

@dataclass
class FittingInfo:
    """Mutable container for fit quality metrics.

    Attributes
    ----------
    method : str
        Fitting method name (``"linear"``, ``"polynomial"``, ``"spline"``).
    r_squared : float
        Coefficient of determination R².
    residuals : np.ndarray
        Per-point residuals (φ_predicted − φ_actual).
    rmse : float
        Root mean squared error.
    equation_str : str
        Human-readable equation for plot annotation.
    params : dict
        Method-specific fitted parameters.
    """

    method: str = ""
    r_squared: float = 0.0
    residuals: np.ndarray = field(default_factory=lambda: np.array([]))
    rmse: float = 0.0
    equation_str: str = ""
    params: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _compute_r_squared(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    ss_res = float(np.sum((y_true - y_pred) ** 2))
    ss_tot = float(np.sum((y_true - y_true.mean()) ** 2))
    if ss_tot == 0.0:
        return 1.0 if ss_res == 0.0 else 0.0
    return 1.0 - ss_res / ss_tot


def _compute_rmse(residuals: np.ndarray) -> float:
    return float(np.sqrt(np.mean(residuals ** 2)))


# ---------------------------------------------------------------------------
# Abstract base class
# ---------------------------------------------------------------------------

class ChargePotentialMapper(ABC):
    """Abstract interface for σ → φ mapping."""

    @abstractmethod
    def fit(self, sigma: np.ndarray, phi: np.ndarray) -> FittingInfo:
        """Fit the σ→φ relationship from calibration data.

        Parameters
        ----------
        sigma : (n,) array
            Surface charge densities in μC/cm².
        phi : (n,) array
            Electrode potentials in V.

        Returns
        -------
        FittingInfo
        """

    @abstractmethod
    def predict(self, sigma: np.ndarray | float) -> np.ndarray:
        """Predict φ from σ using the fitted model.

        Parameters
        ----------
        sigma : float or array
            Surface charge density in μC/cm².

        Returns
        -------
        np.ndarray
            Predicted potential in V.
        """

    @abstractmethod
    def to_dict(self) -> dict[str, Any]:
        """Serialise fitted state to a JSON-compatible dict."""

    @classmethod
    @abstractmethod
    def from_dict(cls, d: dict[str, Any]) -> ChargePotentialMapper:
        """Reconstruct a fitted mapper from a serialised dict."""

    def _check_fitted(self) -> None:
        if not getattr(self, "_fitted", False):
            raise RuntimeError(
                f"{type(self).__name__} has not been fitted — call fit() first"
            )


# ---------------------------------------------------------------------------
# Linear mapper
# ---------------------------------------------------------------------------

class LinearMapper(ChargePotentialMapper):
    """Linear regression: φ = slope · σ + intercept."""

    def __init__(self) -> None:
        self._slope: float = 0.0
        self._intercept: float = 0.0
        self._fitted: bool = False

    def fit(self, sigma: np.ndarray, phi: np.ndarray) -> FittingInfo:
        coeffs = np.polyfit(sigma, phi, 1)
        self._slope = float(coeffs[0])
        self._intercept = float(coeffs[1])
        self._fitted = True

        phi_pred = np.polyval(coeffs, sigma)
        residuals = phi_pred - phi
        r2 = _compute_r_squared(phi, phi_pred)
        rmse = _compute_rmse(residuals)

        sign = "+" if self._intercept >= 0 else "−"
        eq = f"φ = {self._slope:.4f} · σ {sign} {abs(self._intercept):.4f}"

        return FittingInfo(
            method="linear",
            r_squared=r2,
            residuals=residuals,
            rmse=rmse,
            equation_str=eq,
            params={"slope": self._slope, "intercept": self._intercept},
        )

    def predict(self, sigma: np.ndarray | float) -> np.ndarray:
        self._check_fitted()
        s = np.asarray(sigma, dtype=float)
        return self._slope * s + self._intercept

    def to_dict(self) -> dict[str, Any]:
        return {
            "method": "linear",
            "slope": self._slope,
            "intercept": self._intercept,
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> LinearMapper:
        m = cls()
        m._slope = float(d["slope"])
        m._intercept = float(d["intercept"])
        m._fitted = True
        return m


# ---------------------------------------------------------------------------
# Polynomial mapper
# ---------------------------------------------------------------------------

class PolynomialMapper(ChargePotentialMapper):
    """Polynomial regression: φ = Σ c_i · σ^i."""

    def __init__(self, degree: int = 2) -> None:
        self._degree = degree
        self._coeffs: np.ndarray | None = None
        self._fitted: bool = False

    def fit(self, sigma: np.ndarray, phi: np.ndarray) -> FittingInfo:
        self._coeffs = np.polyfit(sigma, phi, self._degree)
        self._fitted = True

        phi_pred = np.polyval(self._coeffs, sigma)
        residuals = phi_pred - phi
        r2 = _compute_r_squared(phi, phi_pred)
        rmse = _compute_rmse(residuals)

        # Build equation string (highest degree first)
        terms: list[str] = []
        for i, c in enumerate(self._coeffs):
            power = self._degree - i
            if power == 0:
                terms.append(f"{c:+.4f}")
            elif power == 1:
                terms.append(f"{c:+.4f}·σ")
            else:
                terms.append(f"{c:+.4f}·σ^{power}")
        eq = "φ = " + " ".join(terms)

        return FittingInfo(
            method="polynomial",
            r_squared=r2,
            residuals=residuals,
            rmse=rmse,
            equation_str=eq,
            params={
                "degree": self._degree,
                "coefficients": self._coeffs.tolist(),
            },
        )

    def predict(self, sigma: np.ndarray | float) -> np.ndarray:
        self._check_fitted()
        s = np.asarray(sigma, dtype=float)
        return np.polyval(self._coeffs, s)

    def to_dict(self) -> dict[str, Any]:
        return {
            "method": "polynomial",
            "degree": self._degree,
            "coefficients": self._coeffs.tolist() if self._coeffs is not None else [],
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> PolynomialMapper:
        m = cls(degree=int(d["degree"]))
        m._coeffs = np.array(d["coefficients"])
        m._fitted = True
        return m


# ---------------------------------------------------------------------------
# Spline mapper
# ---------------------------------------------------------------------------

class SplineMapper(ChargePotentialMapper):
    """Cubic spline interpolation (requires scipy)."""

    def __init__(self) -> None:
        self._spline: Any = None
        self._sigma_min: float = 0.0
        self._sigma_max: float = 0.0
        # Store raw data for serialisation (spline objects are not trivially
        # JSON-serialisable).
        self._sigma_data: np.ndarray | None = None
        self._phi_data: np.ndarray | None = None
        self._fitted: bool = False

    @staticmethod
    def _import_cubic_spline():
        try:
            from scipy.interpolate import CubicSpline
        except ImportError:
            raise ImportError(
                "SplineMapper requires scipy. "
                "Install it with: pip install scipy"
            ) from None
        return CubicSpline

    def fit(self, sigma: np.ndarray, phi: np.ndarray) -> FittingInfo:
        CubicSpline = self._import_cubic_spline()

        # Sort by sigma for spline construction
        order = np.argsort(sigma)
        sigma_sorted = sigma[order]
        phi_sorted = phi[order]

        self._spline = CubicSpline(sigma_sorted, phi_sorted)
        self._sigma_min = float(sigma_sorted[0])
        self._sigma_max = float(sigma_sorted[-1])
        self._sigma_data = sigma_sorted.copy()
        self._phi_data = phi_sorted.copy()
        self._fitted = True

        phi_pred = self._spline(sigma)
        residuals = phi_pred - phi
        r2 = _compute_r_squared(phi, phi_pred)
        rmse = _compute_rmse(residuals)

        return FittingInfo(
            method="spline",
            r_squared=r2,
            residuals=residuals,
            rmse=rmse,
            equation_str="cubic spline interpolation",
            params={
                "sigma_data": sigma_sorted.tolist(),
                "phi_data": phi_sorted.tolist(),
            },
        )

    def predict(self, sigma: np.ndarray | float) -> np.ndarray:
        self._check_fitted()
        s = np.asarray(sigma, dtype=float)
        return np.asarray(self._spline(s), dtype=float)

    def to_dict(self) -> dict[str, Any]:
        return {
            "method": "spline",
            "sigma_data": (self._sigma_data.tolist()
                           if self._sigma_data is not None else []),
            "phi_data": (self._phi_data.tolist()
                         if self._phi_data is not None else []),
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> SplineMapper:
        CubicSpline = cls._import_cubic_spline()
        m = cls()
        m._sigma_data = np.array(d["sigma_data"])
        m._phi_data = np.array(d["phi_data"])
        m._spline = CubicSpline(m._sigma_data, m._phi_data)
        m._sigma_min = float(m._sigma_data[0])
        m._sigma_max = float(m._sigma_data[-1])
        m._fitted = True
        return m


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------

_MAPPER_REGISTRY: dict[str, type[ChargePotentialMapper]] = {
    "linear": LinearMapper,
    "polynomial": PolynomialMapper,
    "spline": SplineMapper,
}


def create_mapper(method: str, **kwargs: Any) -> ChargePotentialMapper:
    """Create a mapper by method name.

    Parameters
    ----------
    method : str
        One of ``"linear"``, ``"polynomial"``, ``"spline"``.
    **kwargs
        Forwarded to the constructor (e.g. ``degree=3`` for polynomial).
    """
    cls = _MAPPER_REGISTRY.get(method)
    if cls is None:
        raise ValueError(
            f"Unknown fitting method {method!r}. "
            f"Choose from: {', '.join(_MAPPER_REGISTRY)}"
        )
    return cls(**kwargs)


def mapper_from_dict(d: dict[str, Any]) -> ChargePotentialMapper:
    """Reconstruct a fitted mapper from a serialised dict."""
    method = d.get("method", "")
    cls = _MAPPER_REGISTRY.get(method)
    if cls is None:
        raise ValueError(f"Unknown method in serialised data: {method!r}")
    return cls.from_dict(d)
