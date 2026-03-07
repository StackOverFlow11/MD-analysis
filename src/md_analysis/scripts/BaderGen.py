"""Generate a VASP Bader-charge work directory from a single MD frame."""

from __future__ import annotations

import shutil
import subprocess
from importlib.resources import as_file, files
from pathlib import Path

from ase import Atoms

from ..config import KEY_VASP_SCRIPT_PATH, get_config
from .utils.IndexMapper import compute_index_map, write_poscar_with_map

DEFAULT_WORKDIR_NAME = "bader"


class BaderGenError(Exception):
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
