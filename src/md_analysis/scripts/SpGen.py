"""Generate CP2K single-point work directories for DeePMD training data."""

from __future__ import annotations

import logging
import shutil
from dataclasses import dataclass
from pathlib import Path

from ase import Atoms
from ase.io import write

from ..config import KEY_CP2K_SCRIPT_PATH, KEY_DP_SP_INP_TEMPLATE_PATH, get_config
from ..exceptions import MDAnalysisError
from ._frame_selector import FrameSelection, iter_selected_frames
from ._inp_utils import modify_inp_for_sp

logger = logging.getLogger(__name__)


class SpGenError(MDAnalysisError):
    """Raised when SP work directory generation fails."""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_sp_workdir(
    atoms: Atoms,
    output_dir: str | Path,
    *,
    inp_template_path: str | Path | None = None,
    cell_abc: tuple[float, float, float] | None = None,
    script_path: str | Path | None = None,
    workdir_name: str = "sp",
    frame: int = 0,
    source: str = "",
) -> Path:
    """Create a CP2K single-point work directory for DeePMD training.

    Parameters
    ----------
    atoms : ase.Atoms
        Single frame (positions only; cell is set from *cell_abc*).
    output_dir : str or Path
        Parent directory under which *workdir_name* will be created.
    inp_template_path : str, Path or None
        Path to the CP2K sp.inp template file.
        If ``None``, falls back to ``KEY_DP_SP_INP_TEMPLATE_PATH`` in config.
    cell_abc : (float, float, float) or None
        Orthogonal cell lengths (A). If ``None``, uses cell from *atoms*.
    script_path : str, Path or None
        Path to a job submission script to copy as ``script.sh``.
        If ``None``, falls back to ``KEY_CP2K_SCRIPT_PATH`` in config.
    workdir_name : str
        Name of the work directory (default ``"sp"``).
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
    SpGenError
        If the inp template is not found or invalid.
    """
    logger.info("Generating SP workdir: frame=%d, workdir=%s", frame, workdir_name)

    # Resolve inp template
    if inp_template_path is None:
        cfg_val = get_config(KEY_DP_SP_INP_TEMPLATE_PATH)
        if cfg_val is not None:
            inp_template_path = cfg_val
    if inp_template_path is None:
        raise SpGenError(
            "No DP SP inp template specified. Provide inp_template_path or "
            "set it via Settings → Set DP SP Inp Template Path."
        )
    inp_template_path = Path(inp_template_path)
    if not inp_template_path.is_file():
        raise FileNotFoundError(f"DP SP inp template not found: {inp_template_path}")

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


def batch_generate_sp_workdirs(
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
    """Batch-generate SP work directories from a CP2K XYZ trajectory for DP training.

    Parameters
    ----------
    xyz_path : str or Path
        CP2K XYZ trajectory file.
    cell_abc : (float, float, float)
        Orthogonal cell lengths (A), e.g. from ``parse_abc_from_restart``.
    output_dir : str or Path
        Parent directory; sub-directories ``sp_t{time}_i{step}`` are created.
    inp_template_path : str, Path or None
        Path to the CP2K sp.inp template file.
        If ``None``, falls back to ``KEY_DP_SP_INP_TEMPLATE_PATH`` in config.
    mode : {"index", "time"}
        Frame selection mode. ``"index"`` uses ``frame_start/end/step``;
        ``"time"`` uses ``time_start_fs/end_fs/step_fs`` (all required).
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

    # Pre-read and modify inp template once
    if inp_template_path is None:
        cfg_val = get_config(KEY_DP_SP_INP_TEMPLATE_PATH)
        if cfg_val is not None:
            inp_template_path = cfg_val
    if inp_template_path is None:
        raise SpGenError(
            "No DP SP inp template specified. Provide inp_template_path or "
            "set it via Settings → Set DP SP Inp Template Path."
        )
    inp_template_path = Path(inp_template_path)
    if not inp_template_path.is_file():
        raise FileNotFoundError(f"DP SP inp template not found: {inp_template_path}")

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

    logger.info("Batch SP: %d frames from %s", len(frames), xyz_path)

    iterator: list[tuple[int, Atoms]] | object = frames
    if verbose:
        from tqdm import tqdm
        iterator = tqdm(frames, desc="SP workdirs", unit="frame", ascii=" =")

    result: list[Path] = []
    for frame_idx, atoms in iterator:
        step = int(atoms.info.get("i", frame_idx))
        time_fs = float(atoms.info.get("time", 0.0))
        workdir_name = f"sp_t{int(time_fs)}_i{step}"

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


# ---------------------------------------------------------------------------
# Agent-facing batch wrapper: structured report + preflight
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SpGenBatchReport:
    """Return value of :func:`generate_sp_batch_with_report`.

    Captures the created directories plus frame metadata so the agent
    layer can surface MCP metrics without reparsing directory names.
    """

    workdirs: tuple[Path, ...]
    n_frames: int
    frame_indices: tuple[int, ...]
    steps: tuple[int, ...]
    times_fs: tuple[float, ...]

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-serializable view of the report.

        ``workdirs`` is converted to ``list[str]``; numeric fields are
        coerced to native Python types so ``json.dumps(report.to_dict())``
        never needs a custom encoder.  Direct ``dataclasses.asdict`` is
        intentionally **not** supported as a serialization boundary.
        """
        return {
            "workdirs": [str(p) for p in self.workdirs],
            "n_frames": int(self.n_frames),
            "frame_indices": [int(v) for v in self.frame_indices],
            "steps": [int(v) for v in self.steps],
            "times_fs": [float(v) for v in self.times_fs],
        }


def generate_sp_batch_with_report(
    xyz_path: str | Path,
    cell_abc: tuple[float, float, float] | list[float],
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
) -> SpGenBatchReport:
    """Batch-generate SP work directories with a structured report.

    Agent-safe wrapper around :func:`batch_generate_sp_workdirs`:

    - Resolves ``inp_template_path`` from ``KEY_DP_SP_INP_TEMPLATE_PATH``
      when omitted; raises :class:`SpGenError` (``validation``) if neither
      parameter nor config supplies one.
    - Resolves ``script_path`` from ``KEY_CP2K_SCRIPT_PATH`` when omitted.
    - Preflights ``xyz_path``, the resolved template path, and the
      resolved script path before any trajectory loading or ``mkdir`` so
      missing files surface uniformly as :class:`FileNotFoundError`
      (mapped to ``file_not_found``).
    - Normalises ``cell_abc`` to a length-3 tuple with positive
      components before loading frames so malformed input fails cheaply.
    - Returns a :class:`SpGenBatchReport` carrying workdirs plus
      per-frame metadata (``frame_indices`` / ``steps`` / ``times_fs``)
      derived from ASE's ``atoms.info`` — **not** from the directory
      names (which truncate the time to an integer).

    Existing same-name ``sp_t{time}_i{step}/`` directories are
    **overwritten** for ``init.xyz`` / ``sp.inp`` / ``script.sh``; this
    mirrors :func:`batch_generate_sp_workdirs` current behaviour.
    """
    # 1. Preflight xyz_path before any other work.
    xyz = Path(xyz_path)
    if not xyz.is_file():
        raise FileNotFoundError(f"xyz_path not found or not a file: {xyz}")

    # 2. Resolve inp template (explicit → user config).  Missing → SpGenError.
    if inp_template_path is None:
        cfg_val = get_config(KEY_DP_SP_INP_TEMPLATE_PATH)
        if cfg_val is not None:
            inp_template_path = cfg_val
    if inp_template_path is None:
        raise SpGenError(
            "No DP SP inp template specified. Provide inp_template_path or "
            "set it via Settings → Set DP SP Inp Template Path."
        )
    tmpl = Path(inp_template_path)
    if not tmpl.is_file():
        raise FileNotFoundError(f"DP SP inp template not found: {tmpl}")

    # 3. Resolve script path (explicit → user config).  Missing → FileNotFoundError.
    if script_path is None:
        cfg_val = get_config(KEY_CP2K_SCRIPT_PATH)
        if cfg_val is not None:
            script_path = cfg_val
    if script_path is not None:
        sp = Path(script_path)
        if not sp.is_file():
            raise FileNotFoundError(
                f"Submission script not found: {sp}"
            )

    # 4. Normalise cell_abc before loading the trajectory.
    cell = tuple(float(x) for x in cell_abc)
    if len(cell) != 3:
        raise ValueError(
            f"cell_abc must have length 3, got length {len(cell)}"
        )
    if any(c <= 0.0 for c in cell):
        raise ValueError(
            f"cell_abc components must be positive, got {cell}"
        )

    # 5. Pre-modify sp.inp text once (cell substitution shared across frames).
    inp_text = tmpl.read_text(encoding="utf-8")
    modified_inp = modify_inp_for_sp(inp_text, cell)

    # 6. Collect frames separately to build report metadata without
    #    reparsing directory names.
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
    for idx, atoms in iter_selected_frames(xyz, selection):
        atoms.set_cell(cell)
        atoms.set_pbc(True)
        frames.append((idx, atoms))

    logger.info(
        "SP wrapper: %d frames selected from %s (mode=%s)",
        len(frames), xyz, mode,
    )

    iterator: list[tuple[int, Atoms]] | object = frames
    if verbose:
        from tqdm import tqdm
        iterator = tqdm(frames, desc="SP workdirs", unit="frame", ascii=" =")

    out_dir = Path(output_dir)

    workdirs: list[Path] = []
    frame_indices: list[int] = []
    steps: list[int] = []
    times_fs: list[float] = []

    for frame_idx, atoms in iterator:
        step = int(atoms.info.get("i", frame_idx))
        time_fs = float(atoms.info.get("time", 0.0))
        workdir_name = f"sp_t{int(time_fs)}_i{step}"

        workdir = out_dir / workdir_name
        workdir.mkdir(parents=True, exist_ok=True)

        write(str(workdir / "init.xyz"), atoms, format="xyz")
        (workdir / "sp.inp").write_text(modified_inp, encoding="utf-8")

        if script_path is not None:
            shutil.copy2(script_path, workdir / "script.sh")

        workdirs.append(workdir)
        frame_indices.append(int(frame_idx))
        steps.append(step)
        times_fs.append(time_fs)

    return SpGenBatchReport(
        workdirs=tuple(workdirs),
        n_frames=len(workdirs),
        frame_indices=tuple(frame_indices),
        steps=tuple(steps),
        times_fs=tuple(times_fs),
    )
