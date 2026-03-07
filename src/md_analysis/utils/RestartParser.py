"""Parse cell parameters from CP2K .restart files."""

from __future__ import annotations

import re
from pathlib import Path

_CELL_BLOCK_RE = re.compile(
    r"&CELL\s*\n(.*?)&END\s+CELL", re.DOTALL | re.IGNORECASE,
)
_VECTOR_RE = re.compile(
    r"^\s*(A|B|C)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)",
    re.MULTILINE | re.IGNORECASE,
)


class RestartParseError(RuntimeError):
    """Raised when parsing a CP2K .restart file fails."""


def parse_abc_from_restart(
    restart_path: str | Path,
) -> tuple[float, float, float]:
    """Parse orthogonal cell lengths (A) from a CP2K .restart file.

    Reads the ``&CELL`` section, extracts diagonal elements of A/B/C vectors.
    Raises :class:`RestartParseError` if the cell block is not found, any
    vector is missing, or the cell is not orthogonal.

    Parameters
    ----------
    restart_path : str or Path
        Path to the ``.restart`` file.

    Returns
    -------
    tuple[float, float, float]
        ``(a, b, c)`` cell lengths in angstrom.
    """
    text = Path(restart_path).read_text(encoding="utf-8")

    match = _CELL_BLOCK_RE.search(text)
    if match is None:
        raise RestartParseError(
            f"No &CELL ... &END CELL block found in {restart_path}"
        )
    cell_block = match.group(1)

    vectors: dict[str, tuple[float, float, float]] = {}
    for m in _VECTOR_RE.finditer(cell_block):
        label = m.group(1).upper()
        vectors[label] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))

    for label in ("A", "B", "C"):
        if label not in vectors:
            raise RestartParseError(
                f"Vector {label} not found in &CELL block of {restart_path}"
            )

    # Validate orthogonality: off-diagonal elements must be negligible.
    tol = 1e-6
    ax, ay, az = vectors["A"]
    bx, by, bz = vectors["B"]
    cx, cy, cz = vectors["C"]
    off_diag = [ay, az, bx, bz, cx, cy]
    if any(abs(v) > tol for v in off_diag):
        raise RestartParseError(
            f"Non-orthogonal cell detected in {restart_path}. "
            f"Off-diagonal elements: A=({ay}, {az}), B=({bx}, {bz}), C=({cx}, {cy})"
        )

    return (ax, by, cz)
