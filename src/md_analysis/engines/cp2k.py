"""CP2K engine adapter.

Implements :class:`ConstraintMDParser` for CP2K constraint-MD point
directories (``.restart`` + ``.LagrangeMultLog`` file pair) and provides
a small module-level facade for the most frequently used CP2K output
readers.

Phase 7b2 status
----------------
- ``CP2KParser``            — Protocol implementation (Phase 5a).
- ``read_constraint_metadata`` / ``read_lambda_series`` —
  Path-friendly facade over ``CP2KParser`` that accepts ``str | Path``
  and constructs the parser internally.
- ``read_fermi_series``     — Reads ``md.out`` Fermi entries via
  ``utils.formats.cp2k.stdout.parse_md_out_fermi`` and converts each
  legacy dict into a typed ``FermiRecord``. The underlying parser
  function still returns ``list[dict]`` for the existing
  dict-based callers (see ``electrochemical.potential.CenterPotential``).
- ``read_continuous_potential_frames`` /
  ``read_distributed_potential_frames`` — Multi-frame discovery for
  the two CP2K potential workflows (mode A: continuous MD cube series
  + ``md.out``; mode B: per-frame ``potential_t*_i*`` SP subdirs). The
  legacy ``electrochemical.potential._frame_source.discover_*_frames``
  names are now thin wrappers that delegate here.
"""

from __future__ import annotations

import logging
from pathlib import Path

from ..utils.constants import BOHR_TO_ANG, TRANSITION_METAL_SYMBOLS
from ..utils.formats.cp2k.colvar import (
    parse_colvar_restart,
    parse_lagrange_mult_log,
)
from ..utils.formats.cp2k.stdout import parse_md_out_fermi, parse_sp_out_fermi
from ..utils.formats.cp2k.xyz import read_xyz_atoms_for_steps
from ..utils.formats.common.cube import (
    discover_cube_files,
    extract_step_from_cube_filename,
    read_cube_atoms,
    read_cube_header_and_values,
)
from ..utils.io._frame_discovery import extract_step_time_from_dirname
from .models import (
    ConstraintMetadata,
    FermiRecord,
    LambdaSeries,
    PotentialFrame,
)

try:
    from ase import Atoms
except ImportError:  # pragma: no cover
    Atoms = object  # type: ignore[misc]

logger = logging.getLogger(__name__)


class CP2KParser:
    """Parser for CP2K constraint-MD point directories.

    Recognises directories containing ``*.restart`` and
    ``*.LagrangeMultLog`` files (the standard CP2K output pair).
    """

    name = "cp2k"

    def is_constraint_directory(self, directory: Path) -> bool:
        if not directory.is_dir():
            return False
        try:
            self._find_restart(directory)
            self._find_log(directory)
        except FileNotFoundError:
            return False
        return True

    def parse_metadata(self, directory: Path) -> ConstraintMetadata:
        return parse_colvar_restart(self._find_restart(directory))

    def parse_lambda_series(self, directory: Path) -> LambdaSeries:
        return parse_lagrange_mult_log(self._find_log(directory))

    # ------------------------------------------------------------------
    # File discovery (private)
    # ------------------------------------------------------------------

    @staticmethod
    def _find_restart(directory: Path) -> Path:
        """Find the primary .restart file (skips .bak and .RESTART.wfn)."""
        candidates = sorted(directory.glob("*.restart"))
        candidates = [
            p for p in candidates
            if ".bak" not in p.name and "RESTART.wfn" not in p.name
        ]
        if not candidates:
            raise FileNotFoundError(f"No .restart file in {directory}")
        # Prefer the one with highest suffix number (e.g. cMD-1_1500.restart)
        return candidates[-1]

    @staticmethod
    def _find_log(directory: Path) -> Path:
        """Find the .LagrangeMultLog file."""
        candidates = list(directory.glob("*.LagrangeMultLog"))
        if not candidates:
            raise FileNotFoundError(f"No .LagrangeMultLog file in {directory}")
        return candidates[0]


# ---------------------------------------------------------------------------
# Module-level facade (Phase 7b1)
# ---------------------------------------------------------------------------


def read_constraint_metadata(
    directory: str | Path,
) -> ConstraintMetadata:
    """Read constraint-MD metadata from a CP2K point directory.

    Thin path-friendly wrapper around :meth:`CP2KParser.parse_metadata`;
    accepts ``str`` or ``Path`` and uses a fresh ``CP2KParser`` instance.
    """
    return CP2KParser().parse_metadata(Path(directory))


def read_lambda_series(
    directory: str | Path,
) -> LambdaSeries:
    """Read the Lagrange-multiplier (λ(t)) series from a CP2K point directory.

    Thin path-friendly wrapper around
    :meth:`CP2KParser.parse_lambda_series`; accepts ``str`` or ``Path``
    and uses a fresh ``CP2KParser`` instance.
    """
    return CP2KParser().parse_lambda_series(Path(directory))


def read_fermi_series(
    md_out_path: str | Path,
) -> list[FermiRecord]:
    """Read the Fermi-energy series from a CP2K ``md.out`` file.

    Returns a list of engine-neutral :class:`FermiRecord` rows
    (``step``, ``time_fs``, ``fermi_raw`` in Hartree). The underlying
    parser (:func:`utils.formats.cp2k.stdout.parse_md_out_fermi`) still
    returns the legacy ``list[dict]`` shape and is unchanged.
    """
    legacy = parse_md_out_fermi(Path(md_out_path))
    return [FermiRecord.from_legacy_dict(d) for d in legacy]


# ---------------------------------------------------------------------------
# Potential-frame facade (Phase 7b2)
# ---------------------------------------------------------------------------


def read_continuous_potential_frames(
    cube_pattern: str,
    *,
    workdir: Path | None = None,
    md_out_path: Path | None = None,
    xyz_path: Path | None = None,
    center_mode: str = "interface",
    metal_elements: set[str] | None = None,
    fermi_unit: str = "au",
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
) -> list[PotentialFrame]:
    """Discover potential frames from a continuous MD run (mode A).

    Parameters
    ----------
    cube_pattern
        Glob pattern for cube files (relative to *workdir*).
    workdir
        Base directory. Defaults to cwd.
    md_out_path
        Path to CP2K ``md.out`` for Fermi energy extraction. If ``None``,
        frames will have ``fermi_raw=None``.
    xyz_path
        Path to XYZ trajectory for interface detection. Required when
        ``center_mode="interface"``.
    center_mode
        ``"interface"`` or ``"cell"``.
    metal_elements
        Explicit metal element set; auto-detected if ``None``.
    fermi_unit
        ``"au"`` (Hartree) or ``"ev"``; placeholder kept for backwards
        compatibility — currently unused downstream (callers convert to
        eV themselves).

    Returns
    -------
    list[PotentialFrame]
        Sorted by step.
    """
    workdir = (workdir or Path(".")).resolve()

    cube_paths = discover_cube_files(
        cube_pattern,
        workdir=workdir,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
    )

    # Parse Fermi energies
    fermi_by_step: dict[int, float] = {}
    if md_out_path is not None:
        md_out_path = Path(md_out_path).resolve()
        if md_out_path.exists():
            records = parse_md_out_fermi(md_out_path)
            records = records[frame_start:frame_end:frame_step]
            for r in records:
                fermi_by_step[int(r["step"])] = float(r["fermi_raw"])

    # Parse xyz for interface detection
    use_interface = center_mode == "interface"
    atoms_by_step: dict[int, Atoms] = {}
    metal_used: set[str] | None = None  # noqa: F841 (kept for parity)

    if use_interface:
        if xyz_path is None:
            xyz_path = workdir / "md-pos-1.xyz"
        xyz_path = Path(xyz_path).resolve()
        if xyz_path.exists():
            needed_steps = {
                s
                for s in (extract_step_from_cube_filename(cp) for cp in cube_paths)
                if s is not None
            }
            atoms_by_step, inferred = read_xyz_atoms_for_steps(
                xyz_path, needed_steps, metal_elements=metal_elements
            )
            metal_used = metal_elements or inferred  # noqa: F841

    # Build frames
    frames: list[PotentialFrame] = []
    for cp in cube_paths:
        step = extract_step_from_cube_filename(cp)
        if step is None:
            step = len(frames)

        header, values = read_cube_header_and_values(cp)

        # Fermi energy
        fermi_raw = fermi_by_step.get(int(step))

        # Atoms with cell from cube header
        frame_atoms: Atoms | None = None
        if use_interface and int(step) in atoms_by_step:
            frame_atoms = atoms_by_step[int(step)]
            frame_atoms.set_cell([
                header.vx_bohr * header.nx * BOHR_TO_ANG,
                header.vy_bohr * header.ny * BOHR_TO_ANG,
                header.vz_bohr * header.nz * BOHR_TO_ANG,
            ])
            frame_atoms.set_pbc(True)

        frames.append(PotentialFrame(
            step=int(step),
            time_fs=None,
            cube_path=cp,
            header=header,
            values=values,
            fermi_raw=fermi_raw,
            atoms=frame_atoms,
        ))

    frames.sort(key=lambda f: f.step)
    return frames


def read_distributed_potential_frames(
    root_dir: Path | str,
    *,
    dir_pattern: str = "potential_t*_i*",
    cube_filename: str = "sp_potential-v_hartree-1_0.cube",
    sp_out_filename: str = "sp.out",
    center_mode: str = "interface",
    metal_elements: set[str] | None = None,
    layer_tol_ang: float = 0.6,
    frame_start: int | None = None,
    frame_end: int | None = None,
    frame_step: int | None = None,
    verbose: bool = False,
) -> list[PotentialFrame]:
    """Discover potential frames from per-frame SP subdirectories (mode B).

    Each subdirectory matching *dir_pattern* under *root_dir* should
    contain a cube file (*cube_filename*) and optionally an output file
    (*sp_out_filename*) with Fermi energy.

    Parameters
    ----------
    root_dir
        Parent directory containing ``potential_t*_i*`` subdirectories.
    dir_pattern
        Glob pattern for subdirectory names.
    cube_filename
        Name of the cube file inside each subdirectory.
    sp_out_filename
        Name of the CP2K output file (for Fermi energy).
    center_mode
        ``"interface"`` or ``"cell"``.
    metal_elements
        Explicit metal element set; auto-detected from first frame if
        ``None``.
    layer_tol_ang
        Layer clustering tolerance (kept for parity; consumer-side
        knob).

    Returns
    -------
    list[PotentialFrame]
        Sorted by step.
    """
    root = Path(root_dir).resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Root directory not found: {root}")

    # Discover and sort subdirectories
    subdirs: list[tuple[Path, int, int]] = []  # (path, time_fs, step)
    for d in sorted(root.glob(dir_pattern)):
        if not d.is_dir():
            continue
        parsed = extract_step_time_from_dirname(d.name)
        if parsed is None:
            logger.warning(
                "Cannot parse step/time from directory name: %s", d.name,
            )
            continue
        step, time_fs = parsed
        cube_path = d / cube_filename
        if not cube_path.exists():
            logger.debug("Cube file missing in %s, skipping", d.name)
            continue
        subdirs.append((d, time_fs, step))

    # Sort by step
    subdirs.sort(key=lambda x: x[2])

    if not subdirs:
        raise FileNotFoundError(
            f"No valid single-point directories found matching {dir_pattern!r} "
            f"in {root} (with {cube_filename})"
        )

    # Apply frame slice
    subdirs = subdirs[frame_start:frame_end:frame_step]
    logger.info(
        "Distributed SP: %d directories in %s", len(subdirs), root
    )

    # Build frames
    use_interface = center_mode == "interface"
    metal_detected: set[str] | None = (
        set(metal_elements) if metal_elements else None
    )
    frames: list[PotentialFrame] = []

    sp_iter = subdirs
    if verbose:
        from tqdm import tqdm

        sp_iter = tqdm(subdirs, desc="Loading SP frames", unit="dir", ascii=" =")

    for d, time_fs, step in sp_iter:
        cube_path = d / cube_filename
        header, values = read_cube_header_and_values(cube_path)

        # Fermi energy from sp.out
        sp_out = d / sp_out_filename
        fermi_raw: float | None = None
        if sp_out.exists():
            fermi_raw = parse_sp_out_fermi(sp_out)

        # Atoms from cube file for interface detection
        frame_atoms: Atoms | None = None
        if use_interface:
            frame_atoms = read_cube_atoms(cube_path, header)
            if metal_detected is None:
                metal_detected = (
                    set(frame_atoms.get_chemical_symbols())
                    & set(TRANSITION_METAL_SYMBOLS)
                )
                if metal_detected:
                    logger.info(
                        "Auto-detected metal elements: %s", metal_detected
                    )

        frames.append(PotentialFrame(
            step=step,
            time_fs=float(time_fs),
            cube_path=cube_path,
            header=header,
            values=values,
            fermi_raw=fermi_raw,
            atoms=frame_atoms,
        ))

    return frames


__all__ = [
    "CP2KParser",
    "read_constraint_metadata",
    "read_lambda_series",
    "read_fermi_series",
    "read_continuous_potential_frames",
    "read_distributed_potential_frames",
]
