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
- All internal computation uses fractional coordinates along the chosen cell axis.
- Layering is done by periodic 1D clustering on fractional coordinates
  (``cluster_1d_periodic`` from ``ClusterUtils``).
- Layer ordering in ``metal_layers_sorted``:
  ``[interface_normal_aligned, slab_interior…, interface_normal_opposed]``,
  traversing through the slab from the +axis-facing interface to the
  −axis-facing interface.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Literal, Sequence

import numpy as np

try:
    from ase import Atoms
except ImportError:  # pragma: no cover
    Atoms = object  # type: ignore[misc]

from ..config import (
    AXIS_MAP,
    DEFAULT_METAL_SYMBOLS,
    INTERFACE_NORMAL_ALIGNED,
    INTERFACE_NORMAL_OPPOSED,
)
from .ClusterUtils import _circular_mean, cluster_1d_periodic, find_largest_gap_periodic


NormalSpec = Literal["a", "b", "c"] | Sequence[float]


@dataclass(frozen=True)
class Layer:
    """One atomic layer grouped along a surface normal."""

    atom_indices: tuple[int, ...]
    center_frac: float  # layer center in fractional coordinate along the normal axis, [0, 1)
    is_interface: bool = False
    interface_label: str | None = None  # "normal_aligned" | "normal_opposed" | None
    # Outward normal (metal → water), only set if is_interface is True.
    normal_unit: tuple[float, float, float] | None = None


@dataclass(frozen=True)
class SurfaceDetectionResult:
    """Result of metal layer detection with interface-layer annotations (both sides under PBC).

    ``metal_layers_sorted`` is ordered as
    ``[interface_normal_aligned, slab_interior…, interface_normal_opposed]``.
    """

    # Use an immutable container for deep immutability (np.ndarray is mutable even in frozen dataclass).
    axis_unit: tuple[float, float, float]  # (x, y, z), unit length
    metal_indices: tuple[int, ...]
    metal_layers_sorted: tuple[Layer, ...]

    def axis_unit_vec(self) -> np.ndarray:
        """Return a numpy view/copy for computations (always safe to mutate externally)."""
        return np.asarray(self.axis_unit, dtype=float)

    def interface_layers(self) -> tuple[Layer, ...]:
        """All layers marked as interface layers, in order [normal_aligned, normal_opposed]."""
        return tuple(layer for layer in self.metal_layers_sorted if layer.is_interface)

    def interface_normal_aligned(self) -> Layer:
        """The interface layer with outward normal aligned with +axis."""
        for layer in self.metal_layers_sorted:
            if layer.interface_label == INTERFACE_NORMAL_ALIGNED:
                return layer
        raise SurfaceGeometryError("No normal_aligned interface layer found")

    def interface_normal_opposed(self) -> Layer:
        """The interface layer with outward normal opposed to +axis."""
        for layer in self.metal_layers_sorted:
            if layer.interface_label == INTERFACE_NORMAL_OPPOSED:
                return layer
        raise SurfaceGeometryError("No normal_opposed interface layer found")


class SurfaceGeometryError(RuntimeError):
    pass


def _normal_unit_from_atoms(atoms: Atoms, normal: NormalSpec) -> np.ndarray:
    if isinstance(normal, str):
        if normal not in AXIS_MAP:
            raise ValueError(f"normal must be one of 'a'/'b'/'c' or a vector, got: {normal!r}")
        axis = AXIS_MAP[normal]
        cell = np.asarray(atoms.cell.array, dtype=float)
        v = cell[axis]
    else:
        v = np.asarray(normal, dtype=float).reshape(3)

    norm = float(np.linalg.norm(v))
    if norm == 0.0:
        raise ValueError("normal vector must be non-zero")
    return v / norm


def _circular_mean_fractional(f: np.ndarray) -> float:
    """
    Mean of fractional coordinates on a circle (robust near 0/1 wrap).

    Returns a value in [0, 1).  Delegates to
    ``ClusterUtils._circular_mean(values, period=1.0)``.
    """
    return _circular_mean(np.asarray(f, dtype=float).ravel(), period=1.0)


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

    Uses periodic 1D clustering on fractional coordinates along the chosen cell
    axis and identifies the largest inter-layer gap as the water region.

    Layer ordering
    --------------
    ``metal_layers_sorted`` is ordered as
    ``[interface_normal_aligned, slab_interior…, interface_normal_opposed]``,
    traversing through the slab from the +axis-facing interface to the
    −axis-facing interface.

    Parameters
    ----------
    atoms
        One frame as `ase.Atoms`.
    metal_symbols
        Symbols considered as *metal slab* atoms, e.g. {"Cu", "Ag"}.
        If omitted, uses `DEFAULT_METAL_SYMBOLS` from `config.py`.
    normal
        Surface normal spec: ``"a"``/``"b"``/``"c"`` → corresponding cell axis.
        Custom vector normals are not supported.
    layer_tol_A
        1D clustering tolerance (Å); converted to fractional coordinates
        internally using the cell-axis length.
    nonmetal_symbols_hint
        Preserved for API compatibility; not used for gap detection.

    Returns
    -------
    SurfaceDetectionResult
    """
    # --- Validate normal ---
    if not isinstance(normal, str) or normal not in AXIS_MAP:
        raise ValueError(
            f"normal must be 'a', 'b', or 'c', got {normal!r}. "
            "Custom vector normals are not supported; use a cell-axis label."
        )
    axis_idx = AXIS_MAP[normal]

    # --- Metal atom selection ---
    metal_symbols_iter = DEFAULT_METAL_SYMBOLS if metal_symbols is None else metal_symbols
    metal_symbols_set = {str(s) for s in metal_symbols_iter}
    if not metal_symbols_set:
        raise ValueError("metal_symbols must be non-empty")

    symbols = np.asarray(atoms.get_chemical_symbols())
    metal_mask = np.isin(symbols, list(metal_symbols_set))
    metal_idx = np.where(metal_mask)[0]
    if metal_idx.size == 0:
        raise SurfaceGeometryError(f"No metal atoms found for symbols={sorted(metal_symbols_set)}")

    # --- Compute axis unit vector and convert tolerance to fractional ---
    axis_unit_vec = _normal_unit_from_atoms(atoms, normal)
    cell = np.asarray(atoms.cell.array, dtype=float)
    axis_length_A = float(np.linalg.norm(cell[axis_idx]))
    if axis_length_A <= 0:
        raise SurfaceGeometryError("Cell axis length must be positive")
    tol_frac = layer_tol_A / axis_length_A

    # --- Fractional coordinates of metal atoms ---
    scaled = np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float)
    metal_frac = scaled[metal_idx, axis_idx]

    # --- Periodic 1D clustering ---
    clusters = cluster_1d_periodic(metal_frac, period=1.0, tol=tol_frac)

    # Build Layer objects (sorted by center_frac, matching cluster_1d_periodic output)
    layers_sorted: list[Layer] = []
    for center_frac, member_indices in clusters:
        atom_idxs = tuple(int(metal_idx[j]) for j in member_indices)
        layers_sorted.append(Layer(atom_indices=atom_idxs, center_frac=float(center_frac)))

    N = len(layers_sorted)
    if N == 0:
        raise SurfaceGeometryError("No metal layers detected.")

    axis_unit_tuple = (float(axis_unit_vec[0]), float(axis_unit_vec[1]), float(axis_unit_vec[2]))

    # --- Single layer: mark as normal_aligned ---
    if N == 1:
        layer = layers_sorted[0]
        nvec = axis_unit_vec
        marked = Layer(
            atom_indices=layer.atom_indices,
            center_frac=layer.center_frac,
            is_interface=True,
            interface_label=INTERFACE_NORMAL_ALIGNED,
            normal_unit=(float(nvec[0]), float(nvec[1]), float(nvec[2])),
        )
        return SurfaceDetectionResult(
            axis_unit=axis_unit_tuple,
            metal_indices=tuple(int(i) for i in metal_idx.tolist()),
            metal_layers_sorted=(marked,),
        )

    # --- Find largest gap (water region) ---
    centers_arr = np.array([L.center_frac for L in layers_sorted], dtype=float)
    low_k, high_k, _gap = find_largest_gap_periodic(centers_arr, period=1.0)

    # --- Build ordered layer list ---
    # Traverse through the slab from normal_aligned (low_k) to normal_opposed (high_k),
    # going in the decreasing-frac direction (through slab interior, not through water gap).
    ordered_indices: list[int] = []
    i = low_k
    while True:
        ordered_indices.append(i)
        if i == high_k:
            break
        i = (i - 1) % N

    # --- Mark interfaces and build result ---
    marked_layers: list[Layer] = []
    for idx in ordered_indices:
        layer = layers_sorted[idx]
        if idx == low_k:
            nvec = axis_unit_vec  # outward normal = +axis
            marked_layers.append(Layer(
                atom_indices=layer.atom_indices,
                center_frac=layer.center_frac,
                is_interface=True,
                interface_label=INTERFACE_NORMAL_ALIGNED,
                normal_unit=(float(nvec[0]), float(nvec[1]), float(nvec[2])),
            ))
        elif idx == high_k:
            nvec = -axis_unit_vec  # outward normal = -axis
            marked_layers.append(Layer(
                atom_indices=layer.atom_indices,
                center_frac=layer.center_frac,
                is_interface=True,
                interface_label=INTERFACE_NORMAL_OPPOSED,
                normal_unit=(float(nvec[0]), float(nvec[1]), float(nvec[2])),
            ))
        else:
            marked_layers.append(layer)

    return SurfaceDetectionResult(
        axis_unit=axis_unit_tuple,
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
        label = ""
        if layer.is_interface and layer.normal_unit is not None:
            label_str = layer.interface_label or "interface"
            label = (
                f" <== {label_str}, "
                f"normal_unit=[{layer.normal_unit[0]:.0f},{layer.normal_unit[1]:.0f},{layer.normal_unit[2]:.0f}]"
            )
        lines.append(
            f"  Layer {i:02d}: center_frac={layer.center_frac:.4f}, n={len(layer.atom_indices)}{label}"
        )
    return "\n".join(lines)
