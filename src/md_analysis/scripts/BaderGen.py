"""Generate VASP Bader-charge work directories from MD frames."""

from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass
from importlib.resources import as_file, files
from pathlib import Path

from ase import Atoms

from ..config import KEY_VASP_SCRIPT_PATH, get_config
from ._frame_selector import FrameSelection, iter_selected_frames
from .utils.IndexMapper import compute_index_map, write_poscar_with_map

logger = logging.getLogger(__name__)

DEFAULT_WORKDIR_NAME = "bader"


from ..exceptions import MDAnalysisError


class BaderGenError(MDAnalysisError):
    """Raised when Bader work directory generation fails."""


def generate_bader_workdir(
    atoms: Atoms,
    output_dir: str | Path,
    *,
    script_path: str | Path | None = None,
    workdir_name: str = DEFAULT_WORKDIR_NAME,
    frame: int = 0,
    source: str = "",
    element_order: tuple[str, ...] | None = None,
    generate_potcar: bool = True,
    direct: bool = True,
) -> Path:
    """Create a VASP single-point work directory for Bader charge analysis.

    Parameters
    ----------
    atoms : ase.Atoms
        Single frame with cell and PBC set.
    output_dir : str or Path
        Parent directory under which *workdir_name* will be created.
    script_path : str, Path or None
        Path to a job submission script to copy as ``script.sh``.
        If ``None``, falls back to the persisted config value.
    workdir_name : str
        Name of the work directory (default ``"bader"``).
    frame : int
        0-based trajectory frame number (metadata for IndexMap).
    source : str
        Source XYZ file path (metadata for IndexMap).
    element_order : tuple[str, ...] or None
        Element grouping order for POSCAR.
    generate_potcar : bool
        If ``True``, invoke ``vaspkit 103`` to generate POTCAR.
    direct : bool
        If ``True``, write fractional coordinates in POSCAR.

    Returns
    -------
    Path
        The created work directory.

    Raises
    ------
    BaderGenError
        If vaspkit is not found or POTCAR generation fails.
    FileNotFoundError
        If *script_path* does not exist.
    """
    logger.info("Generating Bader workdir: frame=%d", frame)

    workdir = Path(output_dir) / workdir_name
    workdir.mkdir(parents=True, exist_ok=True)

    # 1. POSCAR via IndexMapper
    index_map = compute_index_map(
        atoms, frame=frame, source=source, element_order=element_order,
    )
    write_poscar_with_map(atoms, workdir / "POSCAR", index_map, direct=direct)

    # 2. Copy template INCAR / KPOINTS
    template_pkg = files("md_analysis.scripts.template")
    for name in ("INCAR", "KPOINTS"):
        with as_file(template_pkg / name) as src:
            shutil.copy2(src, workdir / name)

    # 3. POTCAR via vaspkit
    if generate_potcar:
        if shutil.which("vaspkit") is None:
            raise BaderGenError(
                "vaspkit not found in PATH; cannot generate POTCAR. "
                "Install vaspkit or set generate_potcar=False."
            )
        result = subprocess.run(
            ["vaspkit"],
            input="103\n",
            capture_output=True,
            text=True,
            cwd=workdir,
            timeout=60,
        )
        if result.returncode != 0:
            raise BaderGenError(
                f"vaspkit exited with code {result.returncode}:\n{result.stderr}"
            )
        if not (workdir / "POTCAR").exists():
            raise BaderGenError(
                "vaspkit completed but POTCAR was not generated.\n"
                f"stdout: {result.stdout[-500:]}"
            )

    # 4. Submission script
    if script_path is None:
        cfg_val = get_config(KEY_VASP_SCRIPT_PATH)
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


def batch_generate_bader_workdirs(
    xyz_path: str | Path,
    cell_abc: tuple[float, float, float],
    output_dir: str | Path,
    *,
    mode: str = "index",
    frame_start: int = 0,
    frame_end: int | None = None,
    frame_step: int = 1,
    time_start_fs: float | None = None,
    time_end_fs: float | None = None,
    time_step_fs: float | None = None,
    script_path: str | Path | None = None,
    element_order: tuple[str, ...] | None = None,
    generate_potcar: bool = True,
    direct: bool = True,
    verbose: bool = False,
) -> list[Path]:
    """Batch-generate Bader work directories from a CP2K XYZ trajectory.

    Parameters
    ----------
    xyz_path : str or Path
        CP2K XYZ trajectory file.
    cell_abc : (float, float, float)
        Orthogonal cell lengths (A), e.g. from ``parse_abc_from_restart``.
    output_dir : str or Path
        Parent directory; sub-directories ``bader_t{time}_i{step}`` are created.
    mode : {"index", "time"}
        Frame selection mode. See ``FrameSelection`` for details.
    frame_start, frame_end, frame_step : int
        Index mode: 0-based frame slice (default: all frames, step=1).
    time_start_fs, time_end_fs, time_step_fs : float or None
        Time mode: inclusive [start, end] range with greedy stepping.
        All three must be provided together.
    script_path : str, Path or None
        Submission script to copy (falls back to config).
    element_order : tuple[str, ...] or None
        Element grouping order for POSCAR.
    generate_potcar : bool
        If True, invoke ``vaspkit 103`` per directory.
    direct : bool
        If True, write fractional coordinates in POSCAR.
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

    logger.info("Batch Bader: %d frames from %s", len(frames), xyz_path)

    iterator: list[tuple[int, Atoms]] | object = frames
    if verbose:
        from tqdm import tqdm
        iterator = tqdm(frames, desc="Bader workdirs", unit="frame", ascii=" =")

    result: list[Path] = []
    for frame_idx, atoms in iterator:
        step = int(atoms.info.get("i", frame_idx))
        time_fs = float(atoms.info.get("time", 0.0))
        workdir_name = f"bader_t{int(time_fs)}_i{step}"

        workdir = generate_bader_workdir(
            atoms,
            output_dir,
            workdir_name=workdir_name,
            frame=frame_idx,
            source=source,
            element_order=element_order,
            script_path=script_path,
            generate_potcar=generate_potcar,
            direct=direct,
        )
        result.append(workdir)

    return result


# ---------------------------------------------------------------------------
# Agent-facing batch wrapper: structured report + preflight
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BaderGenBatchReport:
    """Return value of :func:`generate_bader_batch_with_report`.

    Captures the created directories plus frame metadata so the agent
    layer can surface MCP metrics without reparsing directory names.
    """

    workdirs: tuple[Path, ...]
    n_frames: int
    frame_indices: tuple[int, ...]
    steps: tuple[int, ...]
    times_fs: tuple[float, ...]
    generate_potcar: bool

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
            "generate_potcar": bool(self.generate_potcar),
        }


def generate_bader_batch_with_report(
    xyz_path: str | Path,
    cell_abc: tuple[float, float, float] | list[float],
    output_dir: str | Path,
    *,
    mode: str = "index",
    frame_start: int = 0,
    frame_end: int | None = None,
    frame_step: int = 1,
    time_start_fs: float | None = None,
    time_end_fs: float | None = None,
    time_step_fs: float | None = None,
    script_path: str | Path | None = None,
    element_order: list[str] | tuple[str, ...] | None = None,
    generate_potcar: bool = True,
    direct: bool = True,
    verbose: bool = False,
) -> BaderGenBatchReport:
    """Batch-generate Bader work directories with a structured report.

    Agent-safe wrapper around :func:`batch_generate_bader_workdirs`:

    - Preflight-validates ``xyz_path`` and ``script_path`` (if given) before
      selecting frames, so a missing file surfaces uniformly as
      :class:`FileNotFoundError` (and is mapped to ``file_not_found`` by
      the agent contract) without any directory being created.
    - Accepts ``element_order`` as a JSON-friendly ``list[str]`` and
      normalises to ``tuple[str, ...]`` for the underlying generator.
    - Returns a :class:`BaderGenBatchReport` carrying workdirs plus
      per-frame metadata (``frame_indices`` / ``steps`` / ``times_fs``)
      derived from ASE's ``atoms.info`` — **not** from the directory
      names (which truncate the time to an integer).

    This wrapper **only prepares** VASP work directories; it does not
    submit jobs and does not parse Bader output.  Setting
    ``generate_potcar=True`` invokes ``vaspkit 103`` once per directory
    as a local side effect.

    Parameters and raises match :func:`batch_generate_bader_workdirs`
    except ``element_order`` which also accepts ``list[str]``.
    """
    xyz = Path(xyz_path)
    if not xyz.is_file():
        raise FileNotFoundError(f"xyz_path not found or not a file: {xyz}")
    if script_path is not None:
        sp = Path(script_path)
        if not sp.is_file():
            raise FileNotFoundError(
                f"script_path not found or not a file: {sp}"
            )

    # Normalise cell_abc before any trajectory loading so malformed
    # inputs fail cheaply.  MCP schema delivers a list; ASE wants a tuple.
    cell = tuple(float(x) for x in cell_abc)
    if len(cell) != 3:
        raise ValueError(
            f"cell_abc must have length 3, got length {len(cell)}"
        )

    # Normalise JSON-friendly element_order → tuple for downstream.
    element_order_tuple: tuple[str, ...] | None
    if element_order is None:
        element_order_tuple = None
    else:
        element_order_tuple = tuple(element_order)

    # Collect frames separately so we can build report metadata without
    # reparsing directory names.
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
        "Bader wrapper: %d frames selected from %s (mode=%s)",
        len(frames), xyz, mode,
    )

    iterator: list[tuple[int, Atoms]] | object = frames
    if verbose:
        from tqdm import tqdm
        iterator = tqdm(frames, desc="Bader workdirs", unit="frame", ascii=" =")

    out_dir = Path(output_dir)
    source = str(xyz)

    workdirs: list[Path] = []
    frame_indices: list[int] = []
    steps: list[int] = []
    times_fs: list[float] = []

    for frame_idx, atoms in iterator:
        step = int(atoms.info.get("i", frame_idx))
        time_fs = float(atoms.info.get("time", 0.0))
        workdir_name = f"bader_t{int(time_fs)}_i{step}"

        wd = generate_bader_workdir(
            atoms,
            out_dir,
            workdir_name=workdir_name,
            frame=frame_idx,
            source=source,
            element_order=element_order_tuple,
            script_path=script_path,
            generate_potcar=generate_potcar,
            direct=direct,
        )
        workdirs.append(wd)
        frame_indices.append(int(frame_idx))
        steps.append(step)
        times_fs.append(time_fs)

    return BaderGenBatchReport(
        workdirs=tuple(workdirs),
        n_frames=len(workdirs),
        frame_indices=tuple(frame_indices),
        steps=tuple(steps),
        times_fs=tuple(times_fs),
        generate_potcar=generate_potcar,
    )
