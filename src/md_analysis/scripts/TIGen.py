"""Generate CP2K constrained-MD work directories for Thermodynamic Integration."""

from __future__ import annotations

import logging
import re
import shutil
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import iread, write

from ..config import KEY_CP2K_SCRIPT_PATH, get_config
from ..exceptions import MDAnalysisError
from ..utils.config import AU_TIME_TO_FS
from ..utils.RestartParser.ColvarParser import (
    ColvarRestart,
    parse_colvar_restart,
)

logger = logging.getLogger(__name__)


class TIGenError(MDAnalysisError):
    """Raised when TI work directory generation fails."""


# ---------------------------------------------------------------------------
# Compiled regex patterns for inp file modification
# ---------------------------------------------------------------------------

_PROJECT_RE = re.compile(
    r"^(\s*PROJECT)\s+\S+", re.MULTILINE | re.IGNORECASE,
)
_STEPS_IN_MD_RE = re.compile(
    r"(&MD\b)(.*?)(&END\s+MD)", re.DOTALL | re.IGNORECASE,
)
_STEPS_RE = re.compile(
    r"^(\s*STEPS)\s+\d+", re.MULTILINE | re.IGNORECASE,
)
# TARGET but NOT TARGET_GROWTH — negative lookahead
_TARGET_VALUE_RE = re.compile(
    r"^(\s*TARGET)(?!_GROWTH)\s+(\[.*?\]\s+)?\S+",
    re.MULTILINE | re.IGNORECASE,
)
_TARGET_GROWTH_RE = re.compile(
    r"^(\s*TARGET_GROWTH)\s+(\[.*?\]\s+)?\S+",
    re.MULTILINE | re.IGNORECASE,
)
_COLLECTIVE_BLOCK_RE = re.compile(
    r"(&COLLECTIVE\s*\n)(.*?)(&END\s+COLLECTIVE)",
    re.DOTALL | re.IGNORECASE,
)
_COLVAR_ID_RE = re.compile(
    r"^\s*COLVAR\s+(\d+)", re.MULTILINE | re.IGNORECASE,
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


def _cv_at_step(restart: ColvarRestart, step: int,
                colvar_id: int | None = None) -> float:
    """Compute the target CV value (a.u.) at a given absolute step."""
    c = (restart.colvars[colvar_id]
         if colvar_id is not None
         else restart.colvars.primary)
    dt_au = restart.timestep_fs / AU_TIME_TO_FS
    return c.target_au + (step - restart.step_start) * c.target_growth_au * dt_au


def _load_trajectory_cv(
    xyz_path: str | Path,
    restart: ColvarRestart,
    colvar_id: int | None = None,
) -> list[tuple[int, float, Atoms]]:
    """Load all trajectory frames, returning ``(step, cv_au, atoms)`` triples.

    Each frame gets cell and PBC set from the restart metadata.
    """
    frames: list[tuple[int, float, Atoms]] = []
    ref_symbols: tuple[str, ...] | None = None
    try:
        for atoms in iread(str(xyz_path), index=":"):
            step = int(atoms.info.get("i", 0))
            symbols = tuple(atoms.get_chemical_symbols())
            if ref_symbols is None:
                ref_symbols = symbols
            elif symbols != ref_symbols:
                raise TIGenError(
                    f"Atom ordering changed at step {step}: "
                    f"expected {len(ref_symbols)} atoms "
                    f"({ref_symbols[:3]}...), "
                    f"got ({symbols[:3]}...)"
                )
            cv = _cv_at_step(restart, step, colvar_id)
            atoms.set_cell(restart.cell_abc_ang)
            atoms.set_pbc(True)
            frames.append((step, cv, atoms))
    except TIGenError:
        raise
    except Exception as exc:
        if not frames:
            raise TIGenError(
                f"No frames found in trajectory: {xyz_path}"
            ) from exc
    if not frames:
        raise TIGenError(f"No frames found in trajectory: {xyz_path}")
    return frames


def _snap_to_nearest_frame(
    frames: list[tuple[int, float, Atoms]],
    target_au: float,
) -> tuple[int, float, Atoms]:
    """Find the frame whose CV is closest to *target_au*.

    Returns ``(step, snapped_cv_au, atoms)``.
    """
    cvs = np.array([cv for _, cv, _ in frames])
    idx = int(np.argmin(np.abs(cvs - target_au)))
    step, snapped_cv, atoms = frames[idx]
    return step, snapped_cv, atoms


def _modify_collective_block(
    block_body: str,
    colvar_id_to_set: int | None,
    target_au: float,
) -> str:
    """Modify a single ``&COLLECTIVE`` block body.

    - If *colvar_id_to_set* is ``None``, modify the first block
      (caller passes the primary colvar_id).
    - Always zero ``TARGET_GROWTH``; set ``TARGET`` only for the
      matching ``COLVAR`` id.
    """
    m = _COLVAR_ID_RE.search(block_body)
    block_cid = int(m.group(1)) if m else 1

    # Always zero TARGET_GROWTH (strip unit)
    block_body = _TARGET_GROWTH_RE.sub(r"\1 0", block_body)

    # Set TARGET only for the matching CV
    if colvar_id_to_set is None or block_cid == colvar_id_to_set:
        block_body = _TARGET_VALUE_RE.sub(
            rf"\1 {target_au:.10E}", block_body,
        )

    return block_body


def _modify_inp_for_ti(
    inp_text: str,
    target_au: float,
    steps: int,
    colvar_id: int | None = None,
) -> str:
    """Return *inp_text* modified for constrained-MD (TI sampling point).

    Modifications
    -------------
    - ``PROJECT`` → ``cMD``
    - ``STEPS`` (inside ``&MD``) → *steps*
    - ``TARGET`` → *target_au* (strip ``[unit]``, bare a.u.)
    - ``TARGET_GROWTH`` → ``0`` (strip ``[unit]``)
    - ``&TOPOLOGY``: ensure ``COORD_FILE_NAME init.xyz`` and
      ``COORD_FILE_FORMAT XYZ``
    """
    # 1. PROJECT → cMD
    text = _PROJECT_RE.sub(r"\1 cMD", inp_text)

    # 2. STEPS inside &MD only
    def _replace_steps_in_md(match: re.Match) -> str:
        md_open, md_body, md_close = match.group(1), match.group(2), match.group(3)
        md_body = _STEPS_RE.sub(rf"\g<1> {steps}", md_body)
        return md_open + md_body + md_close

    text = _STEPS_IN_MD_RE.sub(_replace_steps_in_md, text)

    # 3. TARGET and TARGET_GROWTH inside &COLLECTIVE blocks
    #    Determine which colvar_id to match for TARGET replacement.
    primary_cid = colvar_id  # None means "primary (first)"

    def _replace_collective(match: re.Match) -> str:
        header, body, footer = match.group(1), match.group(2), match.group(3)
        body = _modify_collective_block(body, primary_cid, target_au)
        return header + body + footer

    text = _COLLECTIVE_BLOCK_RE.sub(_replace_collective, text)

    # 4. &TOPOLOGY: ensure COORD_FILE_NAME and COORD_FILE_FORMAT
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

    text = _TOPOLOGY_BLOCK_RE.sub(_replace_topology, text)

    return text


def _format_target_dirname(cv_au: float) -> str:
    """Format a directory name from a CV value in a.u."""
    return f"ti_target_{cv_au:.6f}"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_ti_workdir(
    inp_path: str | Path,
    xyz_path: str | Path,
    restart_path: str | Path,
    target_au: float,
    output_dir: str | Path,
    *,
    steps: int = 10000,
    colvar_id: int | None = None,
    workdir_name: str | None = None,
    script_path: str | Path | None = None,
    _preloaded: list[tuple[int, float, Atoms]] | None = None,
    _restart: ColvarRestart | None = None,
    _inp_text: str | None = None,
) -> Path:
    """Create a CP2K constrained-MD work directory for one TI sampling point.

    Parameters
    ----------
    inp_path : str or Path
        CP2K input file to use as template (SG inp).
    xyz_path : str or Path
        SG trajectory file (e.g. ``slowgrowth-pos-1.xyz``).
    restart_path : str or Path
        CP2K ``.restart`` file from the SG simulation.
    target_au : float
        Desired CV target value in atomic units.  Will be snapped to the
        nearest trajectory frame.
    output_dir : str or Path
        Parent directory under which the work directory is created.
    steps : int
        Number of MD steps for constrained-MD (default 10000).
    colvar_id : int or None
        If given, target the constraint with this ``colvar_id``.
        Defaults to the primary (first) constraint.
    workdir_name : str or None
        Name of the work directory.  Auto-generated from the snapped CV
        value if ``None``.
    script_path : str, Path or None
        Path to a job submission script to copy as ``script.sh``.
        If ``None``, falls back to the persisted config value.

    Returns
    -------
    Path
        The created work directory.
    """
    restart = _restart or parse_colvar_restart(restart_path)
    inp_text = _inp_text or Path(inp_path).read_text(encoding="utf-8")
    frames = _preloaded or _load_trajectory_cv(xyz_path, restart, colvar_id)

    step, snapped_cv, atoms = _snap_to_nearest_frame(frames, target_au)
    logger.info(
        "TI workdir: requested target=%.6e, snapped to %.6e (step %d)",
        target_au, snapped_cv, step,
    )

    if workdir_name is None:
        workdir_name = _format_target_dirname(snapped_cv)

    workdir = Path(output_dir) / workdir_name
    workdir.mkdir(parents=True, exist_ok=True)

    # Write modified inp
    modified = _modify_inp_for_ti(inp_text, snapped_cv, steps, colvar_id)
    (workdir / "cMD.inp").write_text(modified, encoding="utf-8")

    # Write init.xyz
    write(str(workdir / "init.xyz"), atoms, format="xyz")

    # Submission script
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


def batch_generate_ti_workdirs(
    inp_path: str | Path,
    xyz_path: str | Path,
    restart_path: str | Path,
    output_dir: str | Path,
    *,
    targets_au: list[float] | np.ndarray | None = None,
    time_initial_fs: float | None = None,
    time_final_fs: float | None = None,
    n_points: int | None = None,
    steps: int = 10000,
    colvar_id: int | None = None,
    script_path: str | Path | None = None,
    verbose: bool = False,
) -> list[Path]:
    """Batch-generate TI constrained-MD work directories.

    Targets can be specified in two ways (mutually exclusive):

    **Numeric mode** — pass *targets_au* directly.

    **Time mode** — pass *time_initial_fs*, *time_final_fs*, *n_points*;
    CV values are computed from the SG target series at evenly spaced
    times, then snapped to the nearest trajectory frame.

    Parameters
    ----------
    inp_path, xyz_path, restart_path
        See :func:`generate_ti_workdir`.
    output_dir : str or Path
        Parent directory for all work directories.
    targets_au : array-like or None
        Explicit CV target values in atomic units (numeric mode).
    time_initial_fs, time_final_fs : float or None
        Time range in femtoseconds (time mode).
    n_points : int or None
        Number of sampling points (time mode).
    steps : int
        Number of MD steps per TI point (default 10000).
    colvar_id : int or None
        Constraint colvar_id (default: primary).
    script_path : str, Path or None
        Submission script to copy (falls back to config).
    verbose : bool
        If True, show a tqdm progress bar.

    Returns
    -------
    list[Path]
        Created work directory paths.
    """
    time_mode = (time_initial_fs is not None
                 or time_final_fs is not None
                 or n_points is not None)
    numeric_mode = targets_au is not None

    if time_mode and numeric_mode:
        raise TIGenError(
            "Cannot specify both targets_au and time_initial_fs/"
            "time_final_fs/n_points. Use one mode only."
        )
    if not time_mode and not numeric_mode:
        raise TIGenError(
            "Must specify either targets_au or "
            "(time_initial_fs, time_final_fs, n_points)."
        )

    restart = parse_colvar_restart(restart_path)
    inp_text = Path(inp_path).read_text(encoding="utf-8")
    frames = _load_trajectory_cv(xyz_path, restart, colvar_id)

    if time_mode:
        if time_initial_fs is None or time_final_fs is None or n_points is None:
            raise TIGenError(
                "Time mode requires all of time_initial_fs, "
                "time_final_fs, and n_points."
            )
        times = np.linspace(time_initial_fs, time_final_fs, n_points)
        # Convert times → step indices → CV values
        resolved: list[float] = []
        for t in times:
            step_idx = round(t / restart.timestep_fs)
            cv = _cv_at_step(restart, step_idx, colvar_id)
            resolved.append(cv)
        targets = resolved
    else:
        targets = list(targets_au)  # type: ignore[arg-type]

    logger.info("Batch TI: %d target points from %s", len(targets), xyz_path)

    iterator: list[float] | object = targets
    if verbose:
        from tqdm import tqdm
        iterator = tqdm(targets, desc="TI workdirs", unit="point", ascii=" =")

    result: list[Path] = []
    for tgt in iterator:
        workdir = generate_ti_workdir(
            inp_path, xyz_path, restart_path, tgt, output_dir,
            steps=steps,
            colvar_id=colvar_id,
            script_path=script_path,
            _preloaded=frames,
            _restart=restart,
            _inp_text=inp_text,
        )
        result.append(workdir)

    return result
