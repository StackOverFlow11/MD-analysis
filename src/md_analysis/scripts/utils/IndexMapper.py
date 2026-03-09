"""Bijective index mapping between CP2K XYZ and VASP POSCAR atom orderings."""

from __future__ import annotations

import base64
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Literal
from urllib.parse import quote, unquote

import numpy as np
from ase import Atoms
from ase.io import write as ase_write


from ...exceptions import MDAnalysisError


class IndexMapError(MDAnalysisError):
    """Base exception for index-mapping operations."""


class IndexMapParseError(IndexMapError):
    """Raised when a POSCAR comment line cannot be decoded into an IndexMap."""


_TAG = "md_analysis::v1"
_FIELD_RE = re.compile(
    r"^md_analysis::v1\s+"
    r"frame=(\d+)\s+"
    r"source=(\S*)\s+"
    r"n=(\d+)\s+"
    r"order=(\S+)\s+"
    r"p2x=(\S+)$"
)


@dataclass(frozen=True)
class IndexMap:
    """Bijective mapping between XYZ and POSCAR atom indices.

    Invariants
    ----------
    - ``poscar_to_xyz[xyz_to_poscar[i]] == i`` for all *i*.
    - ``xyz_to_poscar[poscar_to_xyz[j]] == j`` for all *j*.
    """

    xyz_to_poscar: np.ndarray
    poscar_to_xyz: np.ndarray
    frame: int
    source: str
    n_atoms: int
    element_order: tuple[str, ...]


def _choose_dtype(n: int) -> np.dtype:
    return np.dtype(np.uint16) if n <= 65535 else np.dtype(np.uint32)


def compute_index_map(
    atoms_xyz: Atoms,
    *,
    frame: int = 0,
    source: str = "",
    element_order: tuple[str, ...] | None = None,
) -> IndexMap:
    """Compute the permutation that sorts XYZ atoms into POSCAR element-grouped order.

    Parameters
    ----------
    atoms_xyz : ase.Atoms
        Atoms in CP2K XYZ ordering.
    frame : int
        0-based trajectory frame number.
    source : str
        Source XYZ file path (stored as metadata).
    element_order : tuple[str, ...] or None
        Element grouping order for POSCAR.  ``None`` (default) uses the
        first-occurrence order of symbols in *atoms_xyz*.

    Returns
    -------
    IndexMap
    """
    symbols = atoms_xyz.get_chemical_symbols()
    n = len(symbols)

    if element_order is None:
        seen: list[str] = []
        for s in symbols:
            if s not in seen:
                seen.append(s)
        element_order = tuple(seen)

    elem_rank = {s: i for i, s in enumerate(element_order)}
    for s in symbols:
        if s not in elem_rank:
            raise ValueError(
                f"Element '{s}' in atoms not found in element_order {element_order}"
            )
    if len(elem_rank) != len(set(symbols)):
        missing = set(elem_rank) - set(symbols)
        if missing:
            raise ValueError(
                f"element_order contains elements not present in atoms: {missing}"
            )

    ranks = np.array([elem_rank[s] for s in symbols], dtype=np.int64)
    dt = _choose_dtype(n)
    poscar_to_xyz = np.argsort(ranks, kind="stable").astype(dt)

    xyz_to_poscar = np.empty(n, dtype=dt)
    xyz_to_poscar[poscar_to_xyz] = np.arange(n, dtype=dt)

    return IndexMap(
        xyz_to_poscar=xyz_to_poscar,
        poscar_to_xyz=poscar_to_xyz,
        frame=frame,
        source=source,
        n_atoms=n,
        element_order=element_order,
    )


def encode_comment_line(index_map: IndexMap) -> str:
    """Encode an :class:`IndexMap` into a single POSCAR comment line.

    Format
    ------
    ``md_analysis::v1 frame=<int> source=<urlencoded> n=<int> order=<csv> p2x=<base64>``
    """
    dt = _choose_dtype(index_map.n_atoms)
    raw = index_map.poscar_to_xyz.astype(dt).tobytes()
    b64 = base64.b64encode(raw).decode("ascii")
    source_encoded = quote(index_map.source, safe="./")
    order_csv = ",".join(index_map.element_order)
    return (
        f"{_TAG} frame={index_map.frame} source={source_encoded} "
        f"n={index_map.n_atoms} order={order_csv} p2x={b64}"
    )


def decode_comment_line(comment: str) -> IndexMap:
    """Decode a POSCAR comment line produced by :func:`encode_comment_line`.

    Raises
    ------
    IndexMapParseError
        If the line does not match the expected format or contains an invalid
        permutation.
    """
    comment = comment.strip()
    m = _FIELD_RE.match(comment)
    if m is None:
        if not comment.startswith(_TAG):
            raise IndexMapParseError(
                f"Comment line does not start with '{_TAG}': {comment!r}"
            )
        raise IndexMapParseError(f"Cannot parse comment line fields: {comment!r}")

    frame = int(m.group(1))
    source = unquote(m.group(2))
    n = int(m.group(3))
    element_order = tuple(m.group(4).split(","))
    b64 = m.group(5)

    dt = _choose_dtype(n)
    try:
        raw = base64.b64decode(b64)
    except Exception as exc:
        raise IndexMapParseError(f"Invalid base64 in p2x field: {exc}") from exc

    expected_bytes = n * dt.itemsize
    if len(raw) != expected_bytes:
        raise IndexMapParseError(
            f"p2x byte length {len(raw)} != expected {expected_bytes} "
            f"for n={n}, dtype={dt}"
        )

    poscar_to_xyz = np.frombuffer(raw, dtype=dt).copy()

    # Validate permutation
    if not np.array_equal(np.sort(poscar_to_xyz), np.arange(n, dtype=dt)):
        raise IndexMapParseError(
            "p2x is not a valid permutation of 0..n-1"
        )

    xyz_to_poscar = np.empty(n, dtype=dt)
    xyz_to_poscar[poscar_to_xyz] = np.arange(n, dtype=dt)

    return IndexMap(
        xyz_to_poscar=xyz_to_poscar,
        poscar_to_xyz=poscar_to_xyz,
        frame=frame,
        source=source,
        n_atoms=n,
        element_order=element_order,
    )


def write_poscar_with_map(
    atoms_xyz: Atoms,
    output_path: str | Path,
    index_map: IndexMap,
    *,
    direct: bool = True,
) -> Path:
    """Write a POSCAR file with atoms reordered per *index_map* and the mapping
    encoded in the comment line.

    Parameters
    ----------
    atoms_xyz : ase.Atoms
        Atoms in CP2K XYZ ordering.
    output_path : str or Path
        Destination POSCAR path.
    index_map : IndexMap
        Mapping (typically from :func:`compute_index_map`).
    direct : bool
        If ``True``, write fractional coordinates; otherwise Cartesian.

    Returns
    -------
    Path
        The written file path.
    """
    output_path = Path(output_path)
    atoms_poscar = atoms_xyz[index_map.poscar_to_xyz.tolist()]
    ase_write(str(output_path), atoms_poscar, format="vasp", direct=direct, sort=False)

    # Replace first line with encoded comment
    lines = output_path.read_text().splitlines(keepends=True)
    lines[0] = encode_comment_line(index_map) + "\n"
    output_path.write_text("".join(lines))

    return output_path


def read_index_map_from_poscar(poscar_path: str | Path) -> IndexMap:
    """Read and decode the :class:`IndexMap` from a POSCAR comment line.

    Raises
    ------
    IndexMapParseError
        If the first line is not a valid encoded mapping.
    """
    first_line = Path(poscar_path).read_text().split("\n", 1)[0]
    return decode_comment_line(first_line)


def remap_array(
    data: np.ndarray,
    index_map: IndexMap,
    direction: Literal["xyz_to_poscar", "poscar_to_xyz"],
) -> np.ndarray:
    """Reorder the first axis of *data* according to *index_map*.

    Parameters
    ----------
    data : np.ndarray
        Array whose first axis has length ``index_map.n_atoms``.
    index_map : IndexMap
        The bijective mapping.
    direction : ``"xyz_to_poscar"`` or ``"poscar_to_xyz"``
        Direction of remapping.

    Returns
    -------
    np.ndarray
        Reordered copy of *data*.
    """
    if data.shape[0] != index_map.n_atoms:
        raise IndexMapError(
            f"data length {data.shape[0]} != index_map.n_atoms {index_map.n_atoms}"
        )
    if direction == "xyz_to_poscar":
        return data[index_map.poscar_to_xyz]
    elif direction == "poscar_to_xyz":
        return data[index_map.xyz_to_poscar]
    else:
        raise ValueError(f"Unknown direction: {direction!r}")
