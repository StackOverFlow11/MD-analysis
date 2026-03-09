"""Generate VASP Bader-charge work directories from MD frames."""

from __future__ import annotations

import logging
import shutil
import subprocess
from importlib.resources import as_file, files
from pathlib import Path

from ase import Atoms
from ase.io import iread

from ..config import KEY_VASP_SCRIPT_PATH, get_config
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
    frame_start: int = 0,
    frame_end: int | None = None,
    frame_step: int = 1,
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
    frame_start : int
        0-based index of first frame (default 0).
    frame_end : int or None
        0-based exclusive upper bound (default: all frames).
    frame_step : int
        Step between frames (default 1).
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

    # Collect frames to process using start/stop/step logic.
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
