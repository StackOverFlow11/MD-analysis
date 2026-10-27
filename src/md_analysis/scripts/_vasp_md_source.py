"""VASP MD work directory → ase frame source for batch Gen scripts.

Bridges the gap between a plain VASP MD working directory and the
``ase.Atoms`` frame stream consumed by the batch generators: the
existing :func:`._frame_selector.iter_selected_frames` only reads CP2K
XYZ trajectories, while the planned VASP-MD → CP2K SP pipeline
(``context4agent/requirements/vasp_md_to_cp2k_potential_sp.md``,
section A) needs to slice VASP trajectories the same way.

Route policy
------------
- Trajectory file priority: ``vasprun.xml`` > ``XDATCAR``. Single-file
  format parsing stays with ``ase.io``; this module owns work-directory
  resolution, route policy, and frame-selection semantics (that is why
  it lives in ``scripts/`` and not in ``utils/formats/vasp/``, whose
  Phase 8 placeholders are single-file parsers).
- The XDATCAR route requires a sibling ``POSCAR`` as the element-order
  authority: VASP4-style XDATCAR headers carry no element symbols and
  ase cannot parse them. A synthetic species line derived from the
  POSCAR is injected into the header stream handed to ase; VASP5-style
  headers are cross-checked against the POSCAR and must agree.
- Only ``FrameSelection`` index mode is supported. VASP trajectories
  carry no per-frame simulation-time metadata (no ``atoms.info['time']``
  like CP2K XYZ comment lines), so time mode fails fast with
  :class:`VaspMdSourceError` instead of silently falling back to index
  mode.
- Package boundary: orthorhombic cells only. md-analysis writes
  ``&CELL ABC`` triples from trajectory cells; a non-orthorhombic VASP
  cell would silently produce a wrong box in the generated sp.inp, so
  every yielded frame is checked (off-diagonal tolerance
  ``ORTHOGONALITY_TOL_ANG``, matching
  :attr:`engines.models.CellSpec.is_orthorhombic`).

Memory notes
------------
- ``vasprun.xml``: ``ase.io.iread`` streams ``Atoms`` lazily (one frame
  at a time), but the underlying ase reader parses the whole document
  into an ElementTree before yielding the first frame — all
  ``<calculation>`` blocks stay resident in memory even though the
  ``Atoms`` objects are not accumulated. For very large runs prefer the
  XDATCAR route when available.
- ``XDATCAR``: ase's XDATCAR reader materialises every frame (positions
  only, no per-frame metadata) before returning; acceptable because
  XDATCAR files are typically far smaller than vasprun.xml.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import IO, Iterator

from ase import Atoms
from ase.io import iread, read

from ..exceptions import MDAnalysisError
from ._frame_selector import FrameSelection

logger = logging.getLogger(__name__)


class VaspMdSourceError(MDAnalysisError):
    """Raised when a VASP MD work directory cannot serve as a frame source."""


# Trajectory file candidates, in resolution-priority order.
VASPRUN_FILENAME = "vasprun.xml"
XDATCAR_FILENAME = "XDATCAR"
POSCAR_FILENAME = "POSCAR"

# Off-diagonal cell tolerance (Angstrom). Matches
# ``engines.models.CellSpec.is_orthorhombic`` / the orthogonality check
# in ``utils.formats.cp2k.cell.parse_abc_from_restart``.
ORTHOGONALITY_TOL_ANG = 1e-6


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def select_trajectory_file(md_dir: str | Path) -> Path:
    """Pick the trajectory file of a VASP MD work directory.

    Pure filesystem lookup — no file contents are read. Priority:
    ``vasprun.xml`` (rich per-frame metadata) over ``XDATCAR``.

    Parameters
    ----------
    md_dir : str or Path
        VASP MD working directory.

    Returns
    -------
    Path
        Path to the resolved trajectory file.

    Raises
    ------
    VaspMdSourceError
        If *md_dir* is not a directory or contains neither candidate.
    """
    md_dir = Path(md_dir)
    if not md_dir.is_dir():
        raise VaspMdSourceError(f"VASP MD directory not found: {md_dir}")
    for name in (VASPRUN_FILENAME, XDATCAR_FILENAME):
        candidate = md_dir / name
        if candidate.is_file():
            return candidate
    raise VaspMdSourceError(
        f"No VASP trajectory file found in {md_dir}: looked for "
        f"'{VASPRUN_FILENAME}' (preferred) and '{XDATCAR_FILENAME}'."
    )


def read_element_symbols(md_dir: str | Path) -> list[str]:
    """Element symbols of a VASP MD work directory, in file/POSCAR order.

    Route policy mirrors :func:`iter_vasp_md_frames`:

    - ``vasprun.xml`` present → species list from the first frame's
      ``<atominfo>`` block (parsed via ase; note this parses the whole
      XML document — see the module docstring's memory notes).
    - otherwise ``XDATCAR`` → species list from the sibling ``POSCAR``
      (VASP4-style XDATCAR headers carry no species line, so the POSCAR
      is the element-order authority there).

    Intended for the KIND-coverage check of the VASP-MD → CP2K SP
    pipeline (requirement doc, section B): every element returned here
    must have a ``&KIND`` block in the sp.inp template.

    Parameters
    ----------
    md_dir : str or Path
        VASP MD working directory.

    Returns
    -------
    list[str]
        Expanded chemical symbols, one per atom, in trajectory order.

    Raises
    ------
    VaspMdSourceError
        If no trajectory file is found, the XDATCAR route lacks a valid
        POSCAR, or parsing fails.
    """
    md_dir = Path(md_dir)
    traj = select_trajectory_file(md_dir)

    if traj.name == VASPRUN_FILENAME:
        try:
            first = read(str(traj), index=0, format="vasp-xml")
        except StopIteration:
            raise VaspMdSourceError(
                f"No frames found in {traj} — the vasprun.xml contains no "
                "<calculation> blocks (incomplete or crashed run?)."
            ) from None
        except Exception as exc:
            raise VaspMdSourceError(
                f"Failed to read element symbols from {traj}: {exc}"
            ) from exc
        return first.get_chemical_symbols()

    poscar = _require_poscar(md_dir)
    return _read_poscar_atoms(poscar).get_chemical_symbols()


def iter_vasp_md_frames(
    md_dir: str | Path,
    selection: FrameSelection,
) -> Iterator[tuple[int, Atoms]]:
    """Iterate selected ``(frame_idx, atoms)`` pairs from a VASP MD directory.

    Parameters
    ----------
    md_dir : str or Path
        VASP MD working directory. The trajectory file is resolved as
        ``vasprun.xml`` > ``XDATCAR`` (see :func:`select_trajectory_file`).
        The XDATCAR route additionally requires a sibling ``POSCAR``
        for element order.
    selection : FrameSelection
        Index mode only (``frame_start`` / ``frame_end`` / ``frame_step``,
        0-based, matching :func:`._frame_selector.iter_selected_frames`).

    Yields
    ------
    tuple[int, Atoms]
        0-based frame index within the trajectory file and the
        corresponding ase ``Atoms`` object with ``pbc=True`` and the
        cell read by ase from the trajectory (not from the POSCAR).

    Raises
    ------
    VaspMdSourceError
        On time-mode selections, missing/unresolvable inputs (raised
        eagerly at call time), non-orthorhombic cells, or trajectory
        parse failures (raised during iteration).

    FileNotFoundError
        If the trajectory file vanishes between resolution and reading.
    """
    if selection.mode == "time":
        raise VaspMdSourceError(
            "FrameSelection time mode is not supported for VASP MD frame "
            "sources: vasprun.xml / XDATCAR frames carry no simulation-time "
            "metadata (unlike CP2K XYZ, which embeds atoms.info['time'] in "
            "the comment line). Use index mode (frame_start / frame_end / "
            "frame_step) instead."
        )

    md_dir = Path(md_dir)
    traj = select_trajectory_file(md_dir)
    logger.info("VASP MD frame source: using %s", traj)

    if traj.name == VASPRUN_FILENAME:
        route_frames: Iterator[Atoms] = _iter_vasprun_frames(traj)
    else:
        poscar = _require_poscar(md_dir)
        poscar_symbols = _read_poscar_atoms(poscar).get_chemical_symbols()
        route_frames = _iter_xdatcar_frames(traj, poscar, poscar_symbols)

    return _slice_frames(route_frames, selection, traj)


# ---------------------------------------------------------------------------
# Frame slicing + cell guard (route-independent)
# ---------------------------------------------------------------------------

def _slice_frames(
    frames: Iterator[Atoms],
    selection: FrameSelection,
    source: Path,
) -> Iterator[tuple[int, Atoms]]:
    """Apply index-mode selection to a route frame stream.

    Mirrors :func:`._frame_selector._iter_by_index` slicing semantics;
    additionally guards the package-level orthorhombic-cell boundary on
    every yielded frame.
    """
    next_yield = selection.frame_start
    for idx, atoms in enumerate(frames):
        if selection.frame_end is not None and idx >= selection.frame_end:
            break
        if idx == next_yield:
            _require_orthorhombic(atoms, idx, source)
            atoms.set_pbc(True)
            yield idx, atoms
            next_yield += selection.frame_step


def _require_orthorhombic(atoms: Atoms, frame_idx: int, source: Path) -> None:
    """Raise :class:`VaspMdSourceError` if *atoms* has a non-orthogonal cell."""
    cell = atoms.get_cell()
    worst = max(
        ((i, j, abs(cell[i, j])) for i in range(3) for j in range(3) if i != j),
        key=lambda t: t[2],
    )
    if worst[2] > ORTHOGONALITY_TOL_ANG:
        i, j, val = worst
        raise VaspMdSourceError(
            f"Non-orthorhombic VASP cell at frame {frame_idx} of {source}: "
            f"largest off-diagonal element cell[{i}][{j}] = {val:.6f} A "
            f"exceeds the {ORTHOGONALITY_TOL_ANG:g} A tolerance. "
            "md-analysis only supports orthorhombic cells — proceeding "
            "would write wrong ABC lengths into the CP2K sp.inp &CELL "
            "block. Re-orthogonalize the cell or pick another trajectory."
        )


# ---------------------------------------------------------------------------
# Route: vasprun.xml
# ---------------------------------------------------------------------------

def _iter_vasprun_frames(vasprun: Path) -> Iterator[Atoms]:
    """Yield all frames of a vasprun.xml via ``ase.io.iread`` (lazy)."""
    try:
        yield from iread(str(vasprun), index=":", format="vasp-xml")
    except FileNotFoundError:
        raise
    except Exception as exc:
        raise VaspMdSourceError(
            f"Failed to parse vasprun.xml {vasprun}: {exc}"
        ) from exc


# ---------------------------------------------------------------------------
# Route: XDATCAR (+ sibling POSCAR)
# ---------------------------------------------------------------------------

def _iter_xdatcar_frames(
    xdatcar: Path,
    poscar: Path,
    poscar_symbols: list[str],
) -> Iterator[Atoms]:
    """Yield all frames of an XDATCAR, with POSCAR-validated element order.

    VASP4-style XDATCAR headers (no species line) are transparently
    patched with a POSCAR-derived species line before being handed to
    ase (see :class:`_SpeciesLineInjector`); VASP5-style headers are
    parsed as-is. Either way, every frame's expanded symbol list must
    equal *poscar_symbols*.
    """
    def _ase_frames() -> Iterator[Atoms]:
        try:
            with xdatcar.open("r", encoding="utf-8") as fh:
                stream: IO[str] = fh
                if _xdatcar_lacks_symbols_line(fh):
                    logger.debug(
                        "VASP4-style XDATCAR header in %s: injecting "
                        "species line from %s", xdatcar, poscar,
                    )
                    stream = _SpeciesLineInjector(
                        fh, _species_line_from_symbols(poscar_symbols),
                    )
                yield from iread(stream, index=":", format="vasp-xdatcar")
        except FileNotFoundError:
            raise
        except Exception as exc:
            raise VaspMdSourceError(
                f"Failed to parse XDATCAR {xdatcar}: {exc}"
            ) from exc

    for atoms in _ase_frames():
        if atoms.get_chemical_symbols() != poscar_symbols:
            raise VaspMdSourceError(
                f"XDATCAR {xdatcar} element order "
                f"{atoms.get_chemical_symbols()} disagrees with POSCAR "
                f"{poscar} order {poscar_symbols}. Element order must come "
                "from the POSCAR; make sure both files belong to the same "
                "run."
            )
        yield atoms


def _require_poscar(md_dir: Path) -> Path:
    """Return the sibling POSCAR path, raising if absent (XDATCAR route)."""
    poscar = md_dir / POSCAR_FILENAME
    if not poscar.is_file():
        raise VaspMdSourceError(
            f"The XDATCAR route requires a sibling POSCAR to determine "
            f"element order, but {poscar} was not found. VASP4-style "
            "XDATCAR headers carry no element symbols. Provide the POSCAR "
            "of the same run (or a vasprun.xml)."
        )
    return poscar


def _read_poscar_atoms(poscar: Path) -> Atoms:
    """Read a POSCAR via ase, rejecting VASP4-style headers.

    ase silently *guesses* species for VASP4-style POSCARs (no symbols
    line) from the comment line, which would corrupt the element order
    this module is responsible for — so such headers fail loudly here.
    """
    if _poscar_lacks_symbols_line(poscar):
        raise VaspMdSourceError(
            f"POSCAR {poscar} has no element-symbol line (VASP4-style "
            "header): element order cannot be determined reliably. Add "
            "the symbols line, e.g. '   Pt    O    H', above the counts "
            "line."
        )
    try:
        return read(str(poscar), format="vasp")
    except Exception as exc:
        raise VaspMdSourceError(
            f"Failed to read POSCAR {poscar}: {exc}"
        ) from exc


def _poscar_lacks_symbols_line(poscar: Path) -> bool:
    """True if POSCAR line 6 is the counts line (VASP4-style header).

    Line 6 tokens all parsing as integers means the species-symbols
    line is missing. Truncated headers return False and let ase raise
    the actual parse error.
    """
    with poscar.open("r", encoding="utf-8") as fh:
        for _ in range(5):
            if not fh.readline():
                return False
        tokens = fh.readline().split()
    return bool(tokens) and all(_is_int_token(t) for t in tokens)


def _xdatcar_lacks_symbols_line(fh: IO[str]) -> bool:
    """True if the XDATCAR header omits the species-symbols line (VASP4).

    Peeks line 6 (all-integer tokens → counts line → VASP4 style) and
    rewinds the stream to its start.
    """
    for _ in range(5):
        fh.readline()
    line6 = fh.readline()
    fh.seek(0)
    tokens = line6.split()
    return bool(tokens) and all(_is_int_token(t) for t in tokens)


class _SpeciesLineInjector:
    """``readline()``-only file-like for patching VASP4-style XDATCAR headers.

    ase's XDATCAR reader requires the VASP5-style element-symbol line;
    VASP4-style files omit it. This proxy emits the first five header
    lines unchanged (comment, scale, a1, a2, a3), then a POSCAR-derived
    species line, then delegates every remaining line to the wrapped
    file (which continues with its own counts line).
    """

    _HEADER_LINES_BEFORE_SPECIES = 5

    def __init__(self, fd: IO[str], species_line: str) -> None:
        self._fd = fd
        self._species_line = species_line + "\n"
        self._emitted = 0
        self._injected = False

    def readline(self) -> str:
        if self._emitted < self._HEADER_LINES_BEFORE_SPECIES:
            self._emitted += 1
            return self._fd.readline()
        if not self._injected:
            self._injected = True
            return self._species_line
        return self._fd.readline()


def _species_line_from_symbols(symbols: list[str]) -> str:
    """Collapse an expanded symbol list into one token per consecutive run.

    ``['Pt', 'Pt', 'O', 'H', 'H']`` → ``'Pt O H'`` — aligned with the
    species-counts line the XDATCAR already carries.
    """
    runs: list[str] = []
    for sym in symbols:
        if not runs or runs[-1] != sym:
            runs.append(sym)
    return " ".join(runs)


def _is_int_token(token: str) -> bool:
    """True if *token* parses as an integer (species-counts line check)."""
    try:
        int(token)
    except ValueError:
        return False
    return True
