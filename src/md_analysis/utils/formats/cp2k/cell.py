"""Parse cell parameters from CP2K input and restart files.

Both entry points parse the file with :func:`.inp.parse_inp` (cit tree,
see :mod:`.inp` for tree conventions) and read the tree:

- ``.restart`` / vector-form inputs: ``cell.{a,b,c}`` vectors, diagonal
  elements, with the orthogonality validation in
  :func:`_abc_ang_from_cell_tree` (shared with ``colvar.py``).
- ``md.inp`` / ABC-form inputs: ``cell.abc`` scalars.  cit converts any
  explicit unit token to the XML default unit (angstrom for ``ABC``), so
  ``ABC [bohr]``-style lines are handled correctly too (superset of the
  legacy regex, which only matched ``ABC [angstrom]``).

Legacy error messages are preserved verbatim for the "parsed but the tree
lacks the block" scenarios; files that fail to parse raise
:class:`.inp.InpParseError` (with cit line information) instead.
"""

from __future__ import annotations

from pathlib import Path

from .inp import get_keyword, get_section, parse_inp

from ....exceptions import MDAnalysisError


class CellParseError(MDAnalysisError):
    """Raised when parsing cell parameters from a CP2K file fails."""


def _abc_ang_from_cell_tree(
    tree: dict,
    src_path: str | Path,
) -> tuple[float, float, float]:
    """Extract orthogonal cell lengths from the ``&CELL`` vectors of a tree.

    Reads the ``force_eval.subsys.cell.{a,b,c}`` vectors (angstrom
    semantics in cit trees, both for ``.inp`` and ``.restart`` files),
    takes their diagonal elements and validates orthogonality
    (off-diagonal tolerance 1e-6).  Shared by the restart and colvar
    read paths so the orthogonality policy is identical everywhere.

    Parameters
    ----------
    tree : dict
        A cit simplified tree (see :func:`.inp.parse_inp`).
    src_path : str or Path
        Path interpolated into the error messages.

    Returns
    -------
    tuple[float, float, float]
        ``(a, b, c)`` cell lengths in angstrom.

    Raises
    ------
    CellParseError
        If the cell section is not found, any A/B/C vector is missing,
        or the cell is not orthogonal.
    """
    try:
        cell = get_section(tree, "force_eval.subsys.cell")
    except KeyError:
        raise CellParseError(
            f"No &CELL ... &END CELL block found in {src_path}"
        ) from None
    if not isinstance(cell, dict):
        raise CellParseError(
            f"No &CELL ... &END CELL block found in {src_path}"
        )

    vectors: dict[str, tuple[float, float, float]] = {}
    for label in ("a", "b", "c"):
        vec = cell.get(label)
        if not (isinstance(vec, (tuple, list)) and len(vec) == 3):
            raise CellParseError(
                f"Vector {label.upper()} not found in &CELL block of {src_path}"
            )
        vectors[label.upper()] = (float(vec[0]), float(vec[1]), float(vec[2]))

    # Validate orthogonality: off-diagonal elements must be negligible.
    tol = 1e-6
    ax, ay, az = vectors["A"]
    bx, by, bz = vectors["B"]
    cx, cy, cz = vectors["C"]
    off_diag = [ay, az, bx, bz, cx, cy]
    if any(abs(v) > tol for v in off_diag):
        raise CellParseError(
            f"Non-orthogonal cell detected in {src_path}. "
            f"Off-diagonal elements: A=({ay}, {az}), B=({bx}, {bz}), C=({cx}, {cy})"
        )

    return (ax, by, cz)


def parse_abc_from_restart(
    restart_path: str | Path,
) -> tuple[float, float, float]:
    """Parse orthogonal cell lengths (A) from a CP2K .restart file.

    Parses the file into a cit tree, reads the ``&CELL`` section and
    extracts the diagonal elements of the A/B/C vectors (angstrom).
    Raises :class:`CellParseError` if the cell block is not found, any
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
    tree = parse_inp(restart_path)
    return _abc_ang_from_cell_tree(tree, restart_path)


def parse_abc_from_md_inp(
    md_inp_path: str | Path,
) -> tuple[float, float, float]:
    """Parse orthogonal cell lengths (A) from a CP2K md.inp file.

    Parses the file into a cit tree and reads
    ``force_eval.subsys.cell.abc``.  Any explicit unit token on the
    ``ABC`` line is converted by cit to the default unit (angstrom), so
    ``ABC [angstrom]``, ``ABC [bohr]`` and bare values all yield the
    angstrom triple.

    Parameters
    ----------
    md_inp_path : str or Path
        Path to the CP2K input file (e.g. ``md.inp``).

    Returns
    -------
    tuple[float, float, float]
        ``(a, b, c)`` cell lengths in angstrom.

    Raises
    ------
    CellParseError
        If the ``ABC ...`` line is not found in the ``&CELL`` block.
    """
    md_inp_path = Path(md_inp_path)
    tree = parse_inp(md_inp_path)
    try:
        abc = get_keyword(tree, "force_eval.subsys.cell.abc")
    except KeyError:
        raise CellParseError(
            f"Cannot find `ABC [angstrom] ...` in {md_inp_path}"
        ) from None
    a, b, c = abc
    return float(a), float(b), float(c)
