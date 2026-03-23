"""Trapezoid quadrature, SEM targets, error propagation, optimal allocation."""

from __future__ import annotations

import numpy as np


def compute_trapezoid_weights(xi: np.ndarray) -> np.ndarray:
    """Non-uniform trapezoid weights for TI quadrature.

    Parameters
    ----------
    xi : np.ndarray, shape (K,)
        Constraint point CV values, must be strictly monotonic.

    Returns
    -------
    np.ndarray, shape (K,)
        Trapezoid weights.

    Raises
    ------
    ValueError
        If K < 2 or xi is not strictly monotonic.
    """
    k = len(xi)
    if k < 2:
        raise ValueError(
            f"Trapezoid rule requires at least 2 points, got {k}."
        )

    diffs = np.diff(xi)
    if not (np.all(diffs > 0) or np.all(diffs < 0)):
        raise ValueError(
            "xi must be strictly monotonic (all increasing or all decreasing)."
        )

    weights = np.empty(k)
    weights[0] = diffs[0] / 2.0
    weights[-1] = diffs[-1] / 2.0
    if k > 2:
        weights[1:-1] = (diffs[:-1] + diffs[1:]) / 2.0

    return weights


def _compute_sem_targets(
    weights: np.ndarray,
    epsilon_tol_au: float,
) -> np.ndarray:
    """Per-point SEM upper bound (Scheme A).

    SEM_{max,k} = epsilon_tol_au / (|w_k| * sqrt(K))

    A safety floor prevents degenerate tolerance when w_k -> 0.

    Parameters
    ----------
    weights : np.ndarray, shape (K,)
    epsilon_tol_au : float
        Tolerance in the same unit as lambda_series (Hartree).

    Returns
    -------
    np.ndarray, shape (K,)
        Per-point SEM targets.
    """
    k = len(weights)
    abs_w = np.abs(weights)
    sqrt_k = np.sqrt(k)

    sem_targets = np.where(
        abs_w > 0,
        epsilon_tol_au / (abs_w * sqrt_k),
        np.inf,
    )

    # Safety floor: cap at a reasonable maximum
    safety_cap = 10.0 * epsilon_tol_au
    sem_targets = np.minimum(sem_targets, safety_cap)

    return sem_targets


def _integrate_free_energy(
    forces: np.ndarray,
    weights: np.ndarray,
    sems: np.ndarray,
) -> tuple[float, float]:
    """Integrate free energy with error propagation.

    Parameters
    ----------
    forces : np.ndarray, shape (K,)
        Mean constraint force (lambda_mean) at each point.
    weights : np.ndarray, shape (K,)
        Trapezoid quadrature weights.
    sems : np.ndarray, shape (K,)
        SEM_final at each point.

    Returns
    -------
    tuple[float, float]
        (delta_A, sigma_A) — integrated free energy and propagated error.
    """
    delta_A = float(np.sum(weights * forces))
    sigma_A = float(np.sqrt(np.sum(weights**2 * sems**2)))
    return delta_A, sigma_A


def _suggest_time_allocation(
    weights: np.ndarray,
    sigmas: np.ndarray,
    tau_corrs: np.ndarray,
) -> np.ndarray:
    """Optimal simulation time allocation.

    Time ratios proportional to |w_k| * sigma_k * sqrt(tau_k).

    Parameters
    ----------
    weights : np.ndarray, shape (K,)
    sigmas : np.ndarray, shape (K,)
        Sample standard deviation at each point.
    tau_corrs : np.ndarray, shape (K,)
        Autocorrelation time at each point.

    Returns
    -------
    np.ndarray, shape (K,)
        Normalized allocation ratios (sum to 1).
    """
    raw = np.abs(weights) * sigmas * np.sqrt(np.maximum(tau_corrs, 0.5))
    total = np.sum(raw)
    if total > 0:
        return raw / total
    return np.ones_like(raw) / len(raw)
