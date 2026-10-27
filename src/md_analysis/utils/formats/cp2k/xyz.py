"""CP2K xyz trajectory streaming helpers.

Parses the ``i = <step>`` comment-line tag that CP2K writes into every
frame's comment line in ``*-pos-1.xyz`` files, and walks the (possibly
large) trajectory file once to materialise ``ase.Atoms`` objects for a
specific set of MD steps without holding the full trajectory in memory.

Formats-extraction status
-------------------------
Plain mechanical move from ``electrochemical.potential._frame_source`` —
no behaviour change. The ``(out, inferred_metal)`` return shape of
``read_xyz_atoms_for_steps`` is preserved exactly.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

from ...constants import TRANSITION_METAL_SYMBOLS
from ..common.cube import _float

try:
    from ase import Atoms
except ImportError:  # pragma: no cover
    Atoms = object  # type: ignore[misc]


XYZ_STEP_RE = re.compile(r"\bi\s*=\s*(\d+)\b")


def read_xyz_atoms_for_steps(
    xyz_path: Path,
    steps: set[int],
    *,
    metal_elements: Optional[set[str]] = None,
) -> tuple[dict[int, Atoms], set[str]]:
    """Stream-parse a CP2K xyz trajectory and build :class:`ase.Atoms`
    objects for the requested ``steps``.

    Parameters
    ----------
    xyz_path
        Path to a CP2K ``*-pos-1.xyz`` file. Each frame is expected to
        embed its MD step in the comment line as ``i = <int>``.
    steps
        Set of absolute step numbers to materialise. Frames whose step
        does not appear in ``steps`` are skipped without parsing the
        coordinate block.
    metal_elements
        If given, used verbatim as the metal-element hint. Otherwise the
        first frame's element list is intersected with
        :data:`TRANSITION_METAL_SYMBOLS` to derive an automatic guess.

    Returns
    -------
    tuple
        ``(atoms_by_step, inferred_metals)`` — a dict from step number
        to ``ase.Atoms`` (cell + pbc are NOT set; caller's job) plus the
        metal-element set (the auto-detected guess, or the caller-supplied
        ``metal_elements`` echoed back).

    Raises
    ------
    FileNotFoundError
        If *xyz_path* does not exist.
    """
    if not xyz_path.exists():
        raise FileNotFoundError(xyz_path)

    steps = {int(s) for s in steps}
    out: dict[int, Atoms] = {}
    inferred_metal: Optional[set[str]] = (
        set(metal_elements) if metal_elements is not None else None
    )
    inferred_done = inferred_metal is not None

    with xyz_path.open("r", encoding="utf-8", errors="replace") as f:
        while True:
            natoms_line = f.readline()
            if not natoms_line:
                break
            natoms_line = natoms_line.strip()
            if not natoms_line:
                continue
            try:
                natoms = int(natoms_line.split()[0])
            except ValueError:
                continue

            comment = f.readline()
            if not comment:
                break
            m = XYZ_STEP_RE.search(comment)
            step = int(m.group(1)) if m else None
            need_this = (step is not None) and (step in steps)

            need_parse = need_this or not inferred_done
            if need_parse:
                symbols: list[str] = []
                positions: list[list[float]] = []
                for _ in range(natoms):
                    line = f.readline()
                    if not line:
                        break
                    parts = line.split()
                    if len(parts) < 4:
                        continue
                    symbols.append(parts[0])
                    if need_this:
                        positions.append(
                            [_float(parts[1]), _float(parts[2]), _float(parts[3])]
                        )

                if not inferred_done:
                    inferred_metal = set(symbols) & set(TRANSITION_METAL_SYMBOLS)
                    inferred_done = True

                if need_this and step is not None:
                    out[int(step)] = Atoms(symbols=symbols, positions=positions)
            else:
                for _ in range(natoms):
                    if not f.readline():
                        break

    if inferred_metal is None:
        inferred_metal = set()
    return out, inferred_metal


__all__ = [
    "XYZ_STEP_RE",
    "read_xyz_atoms_for_steps",
]
