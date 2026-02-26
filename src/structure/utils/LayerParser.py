"""
Surface / interface layer detection for a *single frame* (ASE Atoms input).

Goal
----
Given an `ase.Atoms` for one MD frame, identify metal layers and mark which
layers are at the metal/water interfaces. Under periodic boundary conditions
for a water–metal–water slice, there are typically **two** interfaces.

Design principles
-----------------
- Default metal symbols are centralized in `config.py` and can be overridden
  explicitly by caller input.
- Works with non-orthogonal cells by projecting positions onto a chosen normal.
- Layering is done by 1D clustering on the projected coordinate.

This module is intentionally minimal and can be extended later.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Literal, Sequence

import numpy as np

try:
    from ase import Atoms
except ImportError:  # pragma: no cover
    Atoms = object  # type: ignore[misc]

from .config import DEFAULT_METAL_SYMBOLS


NormalSpec = Literal["a", "b", "c"] | Sequence[float]


@dataclass(frozen=True)
class Layer:
    """One atomic layer grouped along a surface normal."""

    atom_indices: tuple[int, ...]
    center_s: float  # layer center along the chosen normal-projection coordinate s (Å)
    # Under PBC for water–metal–water, the slab has two interface sides.
    # Mark interface layers and attach the outward normal (metal -> water).
    is_interface: bool = False
    normal_unit: tuple[float, float, float] | None = None  # only set if is_interface is True


@dataclass(frozen=True)
class SurfaceDetectionResult:
    """Result of metal layer detection with interface-layer annotations (both sides under PBC)."""

    # Use an immutable container for deep immutability (np.ndarray is mutable even in frozen dataclass).
    axis_unit: tuple[float, float, float]  # (x, y, z), unit length, used to sort layers by projection coordinate s
    metal_indices: tuple[int, ...]
    metal_layers_sorted: tuple[Layer, ...]  # sorted from low-s to high-s

    def axis_unit_vec(self) -> np.ndarray:
        """Return a numpy view/copy for computations (always safe to mutate externally)."""
        return np.asarray(self.axis_unit, dtype=float)

    def interface_layers(self) -> tuple[Layer, ...]:
        """All layers marked as interface layers (both sides)."""
        return tuple(layer for layer in self.metal_layers_sorted if layer.is_interface)


class SurfaceGeometryError(RuntimeError):
    pass


def _normal_unit_from_atoms(atoms: Atoms, normal: NormalSpec) -> np.ndarray:
    if isinstance(normal, str):
        axis_map = {"a": 0, "b": 1, "c": 2}
        if normal not in axis_map:
            raise ValueError(f"normal must be one of 'a'/'b'/'c' or a vector, got: {normal!r}")
        axis = axis_map[normal]
        cell = np.asarray(atoms.cell.array, dtype=float)
        v = cell[axis]
    else:
        v = np.asarray(normal, dtype=float).reshape(3)

    norm = float(np.linalg.norm(v))
    if norm == 0.0:
        raise ValueError("normal vector must be non-zero")
    return v / norm


def _project_s(positions: np.ndarray, normal_unit: np.ndarray) -> np.ndarray:
    # positions: (N, 3), normal_unit: (3,)
    return positions @ normal_unit


def _cluster_1d(values_sorted: np.ndarray, tol: float) -> list[tuple[float, int, int]]:
    """
    Cluster sorted 1D values into contiguous groups.

    Returns list of (center, start_idx, end_idx_exclusive) in *sorted-array* indexing.
    """
    if values_sorted.size == 0:
        return []
    if tol <= 0:
        raise ValueError("tol must be > 0")

    clusters: list[tuple[float, int, int]] = []
    start = 0
    acc = [float(values_sorted[0])]
    for i in range(1, values_sorted.size):
        v = float(values_sorted[i])
        center = float(np.mean(acc))
        if abs(v - center) <= tol:
            acc.append(v)
        else:
            clusters.append((float(np.mean(acc)), start, i))
            start = i
            acc = [v]
    clusters.append((float(np.mean(acc)), start, values_sorted.size))
    return clusters


def _circular_mean_fractional(f: np.ndarray) -> float:
    """
    Mean of fractional coordinates on a circle (robust near 0/1 wrap).

    Returns a value in [0, 1).
    """
    f = np.asarray(f, dtype=float).ravel()
    if f.size == 0:
        raise ValueError("cannot compute circular mean of empty array")
    angles = 2.0 * np.pi * f
    z = np.mean(np.cos(angles)) + 1j * np.mean(np.sin(angles))
    if z == 0:
        return float(np.mod(np.mean(f), 1.0))
    mean_angle = np.angle(z)
    if mean_angle < 0:
        mean_angle += 2.0 * np.pi
    return float(mean_angle / (2.0 * np.pi))


def _mic_delta_fractional(df: np.ndarray) -> np.ndarray:
    """
    Minimum-image delta for fractional coordinates in 1D.

    Returns values in [-0.5, 0.5).
    """
    df = np.asarray(df, dtype=float)
    return (df + 0.5) % 1.0 - 0.5


def detect_interface_layers(
    atoms: Atoms,
    *,
    metal_symbols: Iterable[str] | None = None,
    normal: NormalSpec = "c",
    layer_tol_A: float = 0.6,
    nonmetal_symbols_hint: Iterable[str] | None = None,
) -> SurfaceDetectionResult:
    """
    Detect metal layers from a single-frame ASE Atoms, and mark interface layers.

    Notes
    -----
    In the periodic water–metal–water setting, there are exactly **two** interfaces
    per frame. This function always marks exactly one directly water-facing layer per
    side (two interface layers total), and assigns `Layer.normal_unit` for those
    layers pointing metal -> water. The interface selection strategy is fixed and
    not configurable.

    Parameters
    ----------
    atoms
        One frame as `ase.Atoms`.
    metal_symbols
        Symbols considered as *metal slab* atoms, e.g. {"Cu", "Ag"}.
        If omitted, uses `DEFAULT_METAL_SYMBOLS` from `config.py`.
    normal
        Surface normal spec:
        - "a"/"b"/"c" -> use the corresponding cell vector
        - a length-3 vector -> use it directly
    layer_tol_A
        1D clustering tolerance (Å) when grouping metal atoms into layers
        along the normal projection coordinate.
    nonmetal_symbols_hint
        Optional hint for which symbols represent "environment" (water/ions/etc.).
        If omitted, we use all atoms not in `metal_symbols` as non-metal for deciding
        which side is the interface.

    Returns
    -------
    SurfaceDetectionResult
    """
    metal_symbols_iter = DEFAULT_METAL_SYMBOLS if metal_symbols is None else metal_symbols
    metal_symbols_set = {str(s) for s in metal_symbols_iter}
    if not metal_symbols_set:
        raise ValueError("metal_symbols must be non-empty")

    symbols = np.asarray(atoms.get_chemical_symbols())
    metal_mask = np.isin(symbols, list(metal_symbols_set))
    metal_idx = np.where(metal_mask)[0]
    if metal_idx.size == 0:
        raise SurfaceGeometryError(f"No metal atoms found for symbols={sorted(metal_symbols_set)}")

    if nonmetal_symbols_hint is None:
        nonmetal_mask = ~metal_mask
    else:
        nonmetal_set = {str(s) for s in nonmetal_symbols_hint}
        nonmetal_mask = np.isin(symbols, list(nonmetal_set))

    nonmetal_idx = np.where(nonmetal_mask)[0]
    if nonmetal_idx.size == 0:
        raise SurfaceGeometryError(
            "No non-metal atoms found to determine interface normal. "
            "Provide nonmetal_symbols_hint or ensure the frame contains environment atoms."
        )

    axis_unit_vec = _normal_unit_from_atoms(atoms, normal)
    pos = np.asarray(atoms.get_positions(), dtype=float)

    s_all = _project_s(pos, axis_unit_vec)
    s_m = s_all[metal_idx]

    # Build layers from metal atoms
    order = np.argsort(s_m)
    s_m_sorted = s_m[order]
    metal_idx_sorted = metal_idx[order]

    clusters = _cluster_1d(s_m_sorted, tol=layer_tol_A)
    if len(clusters) < 2:
        # Still return what we have, but caller should know it looks odd.
        # (Small slabs or tol too large.)
        pass

    layers: list[Layer] = []
    for center, start, end in clusters:
        idxs = tuple(int(i) for i in metal_idx_sorted[start:end])
        layers.append(Layer(atom_indices=idxs, center_s=float(center)))

    layers_sorted = tuple(sorted(layers, key=lambda L: L.center_s))

    # Mark only the two layers that directly face the non-metal region.
    n_layers = len(layers_sorted)
    if n_layers == 0:
        raise SurfaceGeometryError("No metal layers detected (unexpected).")

    # Assign per-layer interface normal using fractional-coordinate MIC along the chosen axis, if possible.
    axis_map = {"a": 0, "b": 1, "c": 2}
    axis_idx: int | None = axis_map.get(normal) if isinstance(normal, str) else None

    if axis_idx is None:
        scaled = None
        env_frac_axis = None
    else:
        scaled = np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float)
        env_frac_axis = scaled[nonmetal_idx, axis_idx]

    interface_layer_ids: set[int]
    normal_sign_by_layer: dict[int, float] = {}

    if axis_idx is not None and scaled is not None:
        # Use circular ordering along fractional axis and pick the pair across the
        # largest gap, which corresponds to the non-metal region between two surfaces.
        layer_frac: list[tuple[int, float]] = []
        for i, layer in enumerate(layers_sorted):
            fvals = scaled[list(layer.atom_indices), axis_idx]
            layer_frac.append((i, _circular_mean_fractional(fvals)))

        if len(layer_frac) == 1:
            interface_layer_ids = {0}
            normal_sign_by_layer[0] = -1.0
        else:
            layer_frac_sorted = sorted(layer_frac, key=lambda x: x[1])
            ids = [item[0] for item in layer_frac_sorted]
            fracs = [item[1] for item in layer_frac_sorted]

            gaps = []
            for k in range(len(fracs)):
                f0 = fracs[k]
                f1 = fracs[(k + 1) % len(fracs)]
                gap = (f1 - f0) % 1.0
                gaps.append(gap)

            k_max = int(np.argmax(np.asarray(gaps, dtype=float)))
            low_side_id = ids[k_max]
            high_side_id = ids[(k_max + 1) % len(ids)]
            interface_layer_ids = {low_side_id, high_side_id}
            # Moving forward (increasing fractional coordinate) through the largest gap
            # points from low-side interface toward the non-metal region.
            normal_sign_by_layer[low_side_id] = 1.0
            normal_sign_by_layer[high_side_id] = -1.0
    else:
        # Fallback without fractional-axis context: keep two geometric end layers.
        if n_layers == 1:
            interface_layer_ids = {0}
            normal_sign_by_layer[0] = -1.0
        else:
            interface_layer_ids = {0, n_layers - 1}
            normal_sign_by_layer[0] = -1.0
            normal_sign_by_layer[n_layers - 1] = 1.0

    marked_layers: list[Layer] = []
    for i, layer in enumerate(layers_sorted):
        if i not in interface_layer_ids:
            marked_layers.append(layer)
            continue

        # Decide normal sign (+axis or -axis) for the selected interface layer.
        sign = normal_sign_by_layer.get(i)
        if sign is None:
            if scaled is None or env_frac_axis is None or env_frac_axis.size == 0:
                sign = -1.0 if i < n_layers / 2 else 1.0
            else:
                m_frac = _circular_mean_fractional(scaled[list(layer.atom_indices), axis_idx])
                df_mic = _mic_delta_fractional(env_frac_axis - m_frac)
                j = int(np.argmin(np.abs(df_mic)))
                sign = float(np.sign(df_mic[j]))
                if sign == 0.0:
                    sign = -1.0 if i < n_layers / 2 else 1.0

        nvec = axis_unit_vec * sign
        marked_layers.append(
            Layer(
                atom_indices=layer.atom_indices,
                center_s=layer.center_s,
                is_interface=True,
                normal_unit=(float(nvec[0]), float(nvec[1]), float(nvec[2])),
            )
        )

    return SurfaceDetectionResult(
        axis_unit=(float(axis_unit_vec[0]), float(axis_unit_vec[1]), float(axis_unit_vec[2])),
        metal_indices=tuple(int(i) for i in metal_idx.tolist()),
        metal_layers_sorted=tuple(marked_layers),
    )


def format_detection_summary(result: SurfaceDetectionResult) -> str:
    """Human-readable summary for quick debugging/printing."""
    lines: list[str] = []
    lines.append(f"Axis unit (sorting): [{result.axis_unit[0]:.4f}, {result.axis_unit[1]:.4f}, {result.axis_unit[2]:.4f}]")
    lines.append(f"Metal atoms: {len(result.metal_indices)}")
    lines.append(f"Metal layers detected: {len(result.metal_layers_sorted)}")
    for i, layer in enumerate(result.metal_layers_sorted, start=1):
        if layer.is_interface and layer.normal_unit is not None:
            lines.append(
                f"  Layer {i:02d}: center_s={layer.center_s:.3f} Å, n={len(layer.atom_indices)} "
                f"<== interface, normal_unit=[{layer.normal_unit[0]:.0f},{layer.normal_unit[1]:.0f},{layer.normal_unit[2]:.0f}]"
            )
        else:
            lines.append(f"  Layer {i:02d}: center_s={layer.center_s:.3f} Å, n={len(layer.atom_indices)}")
    return "\n".join(lines)
