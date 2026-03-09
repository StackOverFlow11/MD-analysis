"""1D periodic clustering and gap-detection utilities.

Provides generic tools for grouping values along a periodic axis
and identifying the largest inter-cluster gap. Used by
``LayerParser`` for metal-layer detection and by ``CubeParser`` /
``potential`` for interface identification from Cartesian coordinates.
"""

from __future__ import annotations

import numpy as np


def cluster_1d_periodic(
    values: np.ndarray,
    period: float,
    tol: float,
) -> list[tuple[float, np.ndarray]]:
    """Cluster 1D values on a periodic interval ``[0, period)``.

    Parameters
    ----------
    values
        Raw 1D coordinates (arbitrary range; will be folded into ``[0, period)``).
    period
        Period length (Å or 1.0 for fractional coordinates).
    tol
        Clustering tolerance.  A point joins the current cluster when its
        distance to the running cluster mean is ``<= tol``.

    Returns
    -------
    list of (center, member_indices)
        ``center`` is the circular-mean position in ``[0, period)``.
        ``member_indices`` is an ``int64`` array of indices into the
        *original* ``values`` array.  The list is sorted by center in
        ascending order.

    Raises
    ------
    ValueError
        If ``values`` is empty, ``period <= 0``, or ``tol <= 0``.
    """
    values = np.asarray(values, dtype=float).ravel()
    if values.size == 0:
        raise ValueError("values must be non-empty")
    if period <= 0:
        raise ValueError(f"period must be > 0, got {period}")
    if tol <= 0:
        raise ValueError(f"tol must be > 0, got {tol}")

    # Fold into [0, period)
    folded = np.mod(values, period)

    # Sort and cluster
    order = np.argsort(folded)
    sorted_vals = folded[order]

    # Build clusters on the sorted array
    clusters_raw: list[tuple[list[float], list[int]]] = []
    cur_vals: list[float] = [float(sorted_vals[0])]
    cur_idxs: list[int] = [int(order[0])]
    for i in range(1, sorted_vals.size):
        v = float(sorted_vals[i])
        center = float(np.mean(cur_vals))
        if v - center <= tol:
            cur_vals.append(v)
            cur_idxs.append(int(order[i]))
        else:
            clusters_raw.append((cur_vals, cur_idxs))
            cur_vals = [v]
            cur_idxs = [int(order[i])]
    clusters_raw.append((cur_vals, cur_idxs))

    # Merge first/last cluster if they wrap across the periodic boundary
    if len(clusters_raw) > 1:
        first_vals, first_idxs = clusters_raw[0]
        last_vals, last_idxs = clusters_raw[-1]
        # Wrap distance: first cluster mean is near 0, last cluster mean is near period
        wrap_gap = (first_vals[0] + period) - last_vals[-1]
        if wrap_gap <= tol:
            merged_vals = last_vals + [v + period for v in first_vals]
            merged_idxs = last_idxs + first_idxs
            clusters_raw = clusters_raw[1:-1] + [(merged_vals, merged_idxs)]

    # Compute circular-mean centers and package results
    result: list[tuple[float, np.ndarray]] = []
    for vals, idxs in clusters_raw:
        center = _circular_mean(np.array(vals, dtype=float), period)
        result.append((center, np.array(idxs, dtype=int)))

    # Sort by center
    result.sort(key=lambda x: x[0])
    return result


def find_largest_gap_periodic(
    centers_sorted: np.ndarray,
    period: float,
) -> tuple[int, int, float]:
    """Find the largest gap between adjacent cluster centers on a periodic axis.

    Parameters
    ----------
    centers_sorted
        Cluster centers sorted in ascending order within ``[0, period)``.
    period
        Period length.

    Returns
    -------
    (low_idx, high_idx, gap_size)
        ``low_idx`` and ``high_idx`` are indices into ``centers_sorted``
        for the two centers bounding the largest gap.
        ``gap_size`` is the gap width (in the same units as ``period``).

    Raises
    ------
    ValueError
        If fewer than 2 centers are provided.
    """
    centers_sorted = np.asarray(centers_sorted, dtype=float).ravel()
    n = centers_sorted.size
    if n < 2:
        raise ValueError(f"Need >= 2 centers to find a gap, got {n}")

    # Compute all gaps (including the wrap-around gap)
    gaps = np.empty(n, dtype=float)
    gaps[:-1] = np.diff(centers_sorted)
    gaps[-1] = (centers_sorted[0] + period) - centers_sorted[-1]

    k = int(np.argmax(gaps))
    low_idx = k
    high_idx = (k + 1) % n
    return low_idx, high_idx, float(gaps[k])


def gap_midpoint_periodic(
    center_low: float,
    center_high: float,
    period: float,
) -> float:
    """Return the midpoint of a gap on a periodic axis, in ``[0, period)``.

    Parameters
    ----------
    center_low
        Lower bound of the gap.
    center_high
        Upper bound of the gap (may wrap around).
    period
        Period length.
    """
    if period <= 0:
        raise ValueError(f"period must be > 0, got {period}")
    gap = (center_high - center_low) % period
    mid = (center_low + 0.5 * gap) % period
    return float(mid)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _circular_mean(values: np.ndarray, period: float) -> float:
    """Circular mean on ``[0, period)``.

    Raises ``ValueError`` if *values* is empty.
    """
    values = np.asarray(values, dtype=float).ravel()
    if values.size == 0:
        raise ValueError("cannot compute circular mean of empty array")
    angles = 2.0 * np.pi * values / period
    z = np.mean(np.cos(angles)) + 1j * np.mean(np.sin(angles))
    if z == 0:
        return float(np.mod(np.mean(values), period))
    mean_angle = np.angle(z)
    if mean_angle < 0:
        mean_angle += 2.0 * np.pi
    return float(mean_angle * period / (2.0 * np.pi))
