"""Generate CP2K single-point work directories for Hartree potential analysis."""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

from ase import Atoms
from ase.io import write

from ..config import KEY_CP2K_SCRIPT_PATH, KEY_SP_INP_TEMPLATE_PATH, get_config
from ..exceptions import MDAnalysisError
from ._frame_selector import FrameSelection, iter_selected_frames
from ._inp_utils import modify_inp_for_sp

logger = logging.getLogger(__name__)


class PotentialGenError(MDAnalysisError):
    """Raised when potential work directory generation fails."""


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
    modified = modify_inp_for_sp(inp_text, cell_abc)
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
    mode: str = "index",
    frame_start: int = 0,
    frame_end: int | None = None,
    frame_step: int = 1,
    time_start_fs: float | None = None,
    time_end_fs: float | None = None,
    time_step_fs: float | None = None,
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
    mode : {"index", "time"}
        Frame selection mode. See ``FrameSelection`` for details.
    frame_start, frame_end, frame_step : int
        Index mode: 0-based frame slice (default: all frames, step=1).
    time_start_fs, time_end_fs, time_step_fs : float or None
        Time mode: inclusive [start, end] range with greedy stepping.
        All three must be provided together.
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
    modified_inp = modify_inp_for_sp(inp_text, cell_abc)

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

    # Collect frames via unified selector
    selection = FrameSelection(
        mode=mode,  # type: ignore[arg-type]
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        time_start_fs=time_start_fs,
        time_end_fs=time_end_fs,
        time_step_fs=time_step_fs,
    )
    frames: list[tuple[int, Atoms]] = []
    for idx, atoms in iter_selected_frames(xyz_path, selection):
        atoms.set_cell(cell_abc)
        atoms.set_pbc(True)
        frames.append((idx, atoms))

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
