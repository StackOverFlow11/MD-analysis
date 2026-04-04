"""Generate CP2K single-point work directories for Hartree potential analysis."""

from __future__ import annotations

import logging
import re
import shutil
from pathlib import Path

from ase import Atoms
from ase.io import iread, write

from ..config import KEY_CP2K_SCRIPT_PATH, KEY_SP_INP_TEMPLATE_PATH, get_config
from ..exceptions import MDAnalysisError

logger = logging.getLogger(__name__)


class PotentialGenError(MDAnalysisError):
    """Raised when potential work directory generation fails."""


# ---------------------------------------------------------------------------
# Regex for CELL ABC replacement
# ---------------------------------------------------------------------------

_CELL_BLOCK_RE = re.compile(
    r"(&CELL\b)(.*?)(&END\s+CELL)",
    re.DOTALL | re.IGNORECASE,
)
_ABC_LINE_RE = re.compile(
    r"^(\s*ABC\s+)(\[.*?\]\s+)?\S+\s+\S+\s+\S+",
    re.MULTILINE | re.IGNORECASE,
)
_TOPOLOGY_BLOCK_RE = re.compile(
    r"(&TOPOLOGY\b)(.*?)(&END\s+TOPOLOGY)",
    re.DOTALL | re.IGNORECASE,
)
_COORD_FILE_NAME_RE = re.compile(
    r"^(\s*COORD_FILE_NAME)\s+\S+", re.MULTILINE | re.IGNORECASE,
)
_COORD_FILE_FORMAT_RE = re.compile(
    r"^(\s*COORD_FILE_FORMAT)\s+\S+", re.MULTILINE | re.IGNORECASE,
)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _replace_cell_abc(inp_text: str, a: float, b: float, c: float) -> str:
    """Replace the ``ABC`` line inside ``&CELL`` with new values."""

    def _replace_in_cell(match: re.Match) -> str:
        header, body, footer = match.group(1), match.group(2), match.group(3)
        new_body = _ABC_LINE_RE.sub(
            rf"\g<1>[angstrom] {a:.4f}   {b:.4f}   {c:.4f}",
            body,
        )
        return header + new_body + footer

    result = _CELL_BLOCK_RE.sub(_replace_in_cell, inp_text)
    if result == inp_text:
        logger.warning("No &CELL ABC line found in template; cell not updated")
    return result


def _ensure_topology(inp_text: str) -> str:
    """Ensure ``&TOPOLOGY`` has ``COORD_FILE_NAME init.xyz`` and ``COORD_FILE_FORMAT XYZ``."""

    def _replace_topology(match: re.Match) -> str:
        header, body, footer = match.group(1), match.group(2), match.group(3)

        if _COORD_FILE_NAME_RE.search(body):
            body = _COORD_FILE_NAME_RE.sub(r"\1 init.xyz", body)
        else:
            body += "      COORD_FILE_NAME init.xyz\n"

        if _COORD_FILE_FORMAT_RE.search(body):
            body = _COORD_FILE_FORMAT_RE.sub(r"\1 XYZ", body)
        else:
            body += "      COORD_FILE_FORMAT XYZ\n"

        return header + body + footer

    return _TOPOLOGY_BLOCK_RE.sub(_replace_topology, inp_text)


def _modify_inp_for_sp(inp_text: str, cell_abc: tuple[float, float, float]) -> str:
    """Return *inp_text* modified for single-point potential calculation.

    Modifications:
    - ``&CELL ABC`` → updated cell parameters
    - ``&TOPOLOGY``: ensure ``COORD_FILE_NAME init.xyz`` + ``COORD_FILE_FORMAT XYZ``
    """
    text = _replace_cell_abc(inp_text, *cell_abc)
    text = _ensure_topology(text)
    return text


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_potential_workdir(
    atoms: Atoms,
    output_dir: str | Path,
    *,
    inp_template_path: str | Path | None = None,
    cell_abc: tuple[float, float, float] | None = None,
    script_path: str | Path | None = None,
    workdir_name: str = "potential",
    frame: int = 0,
    source: str = "",
) -> Path:
    """Create a CP2K single-point work directory for Hartree potential analysis.

    Parameters
    ----------
    atoms : ase.Atoms
        Single frame (positions only; cell is set from *cell_abc*).
    output_dir : str or Path
        Parent directory under which *workdir_name* will be created.
    inp_template_path : str, Path or None
        Path to the CP2K sp.inp template file.
        If ``None``, falls back to the persisted config value.
    cell_abc : (float, float, float) or None
        Orthogonal cell lengths (A). If ``None``, uses cell from *atoms*.
    script_path : str, Path or None
        Path to a job submission script to copy as ``script.sh``.
        If ``None``, falls back to the persisted config value.
    workdir_name : str
        Name of the work directory (default ``"potential"``).
    frame : int
        0-based trajectory frame number (metadata only).
    source : str
        Source XYZ file path (metadata only).

    Returns
    -------
    Path
        The created work directory.

    Raises
    ------
    PotentialGenError
        If the inp template is not found or invalid.
    """
    logger.info("Generating potential workdir: frame=%d, workdir=%s", frame, workdir_name)

    # Resolve inp template
    if inp_template_path is None:
        cfg_val = get_config(KEY_SP_INP_TEMPLATE_PATH)
        if cfg_val is not None:
            inp_template_path = cfg_val
    if inp_template_path is None:
        raise PotentialGenError(
            "No SP inp template specified. Provide inp_template_path or "
            "set it via Settings → Set SP Inp Template Path."
        )
    inp_template_path = Path(inp_template_path)
    if not inp_template_path.is_file():
        raise FileNotFoundError(f"SP inp template not found: {inp_template_path}")

    # Resolve cell
    if cell_abc is None:
        cell = atoms.get_cell()
        cell_abc = (float(cell[0, 0]), float(cell[1, 1]), float(cell[2, 2]))

    workdir = Path(output_dir) / workdir_name
    workdir.mkdir(parents=True, exist_ok=True)

    # 1. Write init.xyz
    write(str(workdir / "init.xyz"), atoms, format="xyz")

    # 2. Modify and write sp.inp
    inp_text = inp_template_path.read_text(encoding="utf-8")
    modified = _modify_inp_for_sp(inp_text, cell_abc)
    (workdir / "sp.inp").write_text(modified, encoding="utf-8")

    # 3. Submission script
    if script_path is None:
        cfg_val = get_config(KEY_CP2K_SCRIPT_PATH)
        if cfg_val is not None:
            script_path = cfg_val

    if script_path is not None:
        script_path = Path(script_path)
        if not script_path.is_file():
            raise FileNotFoundError(
                f"Submission script not found: {script_path}"
            )
        shutil.copy2(script_path, workdir / "script.sh")

    return workdir


def batch_generate_potential_workdirs(
    xyz_path: str | Path,
    cell_abc: tuple[float, float, float],
    output_dir: str | Path,
    *,
    inp_template_path: str | Path | None = None,
    frame_start: int = 0,
    frame_end: int | None = None,
    frame_step: int = 1,
    script_path: str | Path | None = None,
    verbose: bool = False,
) -> list[Path]:
    """Batch-generate SP potential work directories from a CP2K XYZ trajectory.

    Parameters
    ----------
    xyz_path : str or Path
        CP2K XYZ trajectory file.
    cell_abc : (float, float, float)
        Orthogonal cell lengths (A), e.g. from ``parse_abc_from_restart``.
    output_dir : str or Path
        Parent directory; sub-directories ``potential_t{time}_i{step}`` are created.
    inp_template_path : str, Path or None
        Path to the CP2K sp.inp template file.
        If ``None``, falls back to the persisted config value.
    frame_start : int
        0-based index of first frame (default 0).
    frame_end : int or None
        0-based exclusive upper bound (default: all frames).
    frame_step : int
        Step between frames (default 1).
    script_path : str, Path or None
        Submission script to copy (falls back to config).
    verbose : bool
        If True, show a tqdm progress bar.

    Returns
    -------
    list[Path]
        Created work directory paths.
    """
    xyz_path = Path(xyz_path)
    output_dir = Path(output_dir)
    source = str(xyz_path)

    # Pre-read and modify inp template once
    if inp_template_path is None:
        cfg_val = get_config(KEY_SP_INP_TEMPLATE_PATH)
        if cfg_val is not None:
            inp_template_path = cfg_val
    if inp_template_path is None:
        raise PotentialGenError(
            "No SP inp template specified. Provide inp_template_path or "
            "set it via Settings → Set SP Inp Template Path."
        )
    inp_template_path = Path(inp_template_path)
    if not inp_template_path.is_file():
        raise FileNotFoundError(f"SP inp template not found: {inp_template_path}")

    inp_text = inp_template_path.read_text(encoding="utf-8")
    modified_inp = _modify_inp_for_sp(inp_text, cell_abc)

    # Resolve script path once
    if script_path is None:
        cfg_val = get_config(KEY_CP2K_SCRIPT_PATH)
        if cfg_val is not None:
            script_path = cfg_val
    if script_path is not None:
        script_path = Path(script_path)
        if not script_path.is_file():
            raise FileNotFoundError(
                f"Submission script not found: {script_path}"
            )

    # Collect frames
    next_yield = frame_start
    frames: list[tuple[int, Atoms]] = []
    for idx, atoms in enumerate(iread(str(xyz_path), index=":")):
        if frame_end is not None and idx >= frame_end:
            break
        if idx == next_yield:
            atoms.set_cell(cell_abc)
            atoms.set_pbc(True)
            frames.append((idx, atoms))
            next_yield += frame_step

    logger.info("Batch potential: %d frames from %s", len(frames), xyz_path)

    iterator: list[tuple[int, Atoms]] | object = frames
    if verbose:
        from tqdm import tqdm
        iterator = tqdm(frames, desc="Potential workdirs", unit="frame", ascii=" =")

    result: list[Path] = []
    for frame_idx, atoms in iterator:
        step = int(atoms.info.get("i", frame_idx))
        time_fs = float(atoms.info.get("time", 0.0))
        workdir_name = f"potential_t{int(time_fs)}_i{step}"

        workdir = Path(output_dir) / workdir_name
        workdir.mkdir(parents=True, exist_ok=True)

        # Write init.xyz
        write(str(workdir / "init.xyz"), atoms, format="xyz")

        # Write pre-modified sp.inp (cell already substituted)
        (workdir / "sp.inp").write_text(modified_inp, encoding="utf-8")

        # Copy script
        if script_path is not None:
            shutil.copy2(script_path, workdir / "script.sh")

        result.append(workdir)

    return result
