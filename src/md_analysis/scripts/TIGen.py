"""Generate CP2K constrained-MD work directories for Thermodynamic Integration."""

from __future__ import annotations

import copy
import logging
import shutil
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import iread, write

from ..config import KEY_CP2K_SCRIPT_PATH, get_config
from ..engines.models import ConstraintMetadata
from ..exceptions import MDAnalysisError
from ..utils.constants import AU_TIME_TO_FS
from ..utils.formats.cp2k.inp import (
    InpParseError,
    generate_inp,
    parse_inp,
)
from ..engines.cp2k import read_constraint_metadata_from_restart
from ._inp_tree_edit import (
    ensure_topology_init_xyz,
    set_collective_target,
    set_md_steps,
    set_project,
    zero_target_growth,
)

logger = logging.getLogger(__name__)


class TIGenError(MDAnalysisError):
    """Raised when TI work directory generation fails."""


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _cv_at_step(restart: ConstraintMetadata, step: int,
                colvar_id: int | None = None) -> float:
    """Compute the target CV value (a.u.) at a given absolute step."""
    c = (restart.colvars[colvar_id]
         if colvar_id is not None
         else restart.colvars.primary)
    dt_au = restart.timestep_fs / AU_TIME_TO_FS
    return c.target_au + (step - restart.step_start) * c.target_growth_au * dt_au


def _load_trajectory_cv(
    xyz_path: str | Path,
    restart: ConstraintMetadata,
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


def _parse_sg_inp(inp_path: str | Path) -> dict:
    """Parse the SG inp template into a cit tree.

    A template that cannot be parsed raises :class:`TIGenError` carrying
    the template path, with the underlying :class:`InpParseError`
    preserved as ``__cause__`` (plan §3.3); generate-side failures are
    programming errors and propagate untranslated.
    """
    try:
        return parse_inp(inp_path)
    except InpParseError as e:
        raise TIGenError(
            f"Cannot parse SG inp template {inp_path}: {e}"
        ) from e


def _render_ti_inp(
    base_tree: dict,
    target_au: float,
    steps: int,
    colvar_id: int | None = None,
) -> str:
    """Return the ``cMD.inp`` text for one TI sampling point.

    *base_tree* is deep-copied first (per-target discipline, plan §3.3:
    batch entry points parse the template once and reuse the tree for
    every target, so the caller's tree must never be mutated).  Edits:

    - ``GLOBAL.PROJECT_NAME`` → ``cMD`` (canonical key ``project_name``)
    - ``MOTION.MD.STEPS`` → *steps*
    - ``TARGET`` → *target_au* as a bare internal-unit value, in the
      ``&COLLECTIVE`` block selected by *colvar_id* (``None`` = the
      primary / first block, matching the reader-side ``colvars.primary``
      semantics)
    - ``TARGET_GROWTH`` → ``0`` in **every** ``&COLLECTIVE`` block
    - ``&TOPOLOGY`` → ``COORD_FILE_NAME init.xyz`` / ``COORD_FILE_FORMAT
      XYZ`` (the section is created when missing)
    """
    tree = copy.deepcopy(base_tree)
    set_project(tree, "cMD")
    set_md_steps(tree, steps)
    set_collective_target(tree, target_au, colvar_id)
    zero_target_growth(tree)
    ensure_topology_init_xyz(tree)
    return generate_inp(tree)


def _format_target_dirname(cv_au: float) -> str:
    """Format a directory name from a CV value in a.u."""
    return f"ti_target_{cv_au:.6f}"


# ---------------------------------------------------------------------------
# Agent-facing planning helper (MVP: primary CV only)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _PlannedTarget:
    """One planned TI target after snap-to-nearest-frame.

    Used internally by :func:`generate_ti_batch_with_report` so the agent
    layer can surface ``requested_cv`` / ``snapped_cv`` / ``snap_delta``
    as metrics without reparsing directory names (which lose precision
    at 6 decimal places).
    """

    requested_cv: float
    snapped_cv: float
    snap_delta: float
    frame_step: int
    dirname: str


def _plan_ti_targets(
    inp_path: str | Path,  # noqa: ARG001 — kept in signature for parity w/ batch
    xyz_path: str | Path,
    restart_path: str | Path,
    *,
    targets_au: list[float] | np.ndarray | None = None,
    time_range: dict[str, float | int] | None = None,
    colvar_id: int | None = None,
) -> list[_PlannedTarget]:
    """Plan a list of TI targets.

    Accepts exactly one target specification:

    - ``targets_au``: explicit list of CV values in atomic units
    - ``time_range``: dict ``{time_initial_fs, time_final_fs, n_points}``

    ``colvar_id`` selects which constraint the CV is read from
    (``None`` = primary / first colvar).  Each planned target carries
    the snapped CV (from the nearest SG trajectory frame) plus metadata
    for agent-layer metrics and collision checks.  ``inp_path`` is
    accepted for signature symmetry with
    :func:`batch_generate_ti_workdirs` but not read here (the inp file
    is consumed only at generation time).
    """
    numeric_mode = targets_au is not None
    time_mode = time_range is not None
    if numeric_mode and time_mode:
        raise TIGenError(
            "Cannot specify both targets_au and time_range. Use one mode only."
        )
    if not numeric_mode and not time_mode:
        raise TIGenError(
            "Must specify either targets_au or time_range."
        )

    restart = read_constraint_metadata_from_restart(restart_path)
    frames = _load_trajectory_cv(xyz_path, restart, colvar_id)

    if time_mode:
        required_keys = {"time_initial_fs", "time_final_fs", "n_points"}
        missing = required_keys - set(time_range.keys())
        if missing:
            raise TIGenError(
                f"time_range missing required keys: {sorted(missing)}"
            )
        t_initial = float(time_range["time_initial_fs"])
        t_final = float(time_range["time_final_fs"])
        n_points = int(time_range["n_points"])
        if n_points < 2:
            raise TIGenError(
                f"time_range.n_points must be >= 2, got {n_points}"
            )
        if t_initial > t_final:
            raise TIGenError(
                f"time_range.time_initial_fs ({t_initial}) must be <= "
                f"time_final_fs ({t_final})"
            )
        times = np.linspace(t_initial, t_final, n_points)
        requested: list[float] = []
        for t in times:
            step_idx = round(t / restart.timestep_fs)
            requested.append(_cv_at_step(restart, step_idx, colvar_id))
    else:
        requested = [float(x) for x in targets_au]  # type: ignore[union-attr]

    planned: list[_PlannedTarget] = []
    for rq in requested:
        step, snapped, _atoms = _snap_to_nearest_frame(frames, rq)
        planned.append(_PlannedTarget(
            requested_cv=float(rq),
            snapped_cv=float(snapped),
            snap_delta=abs(float(rq) - float(snapped)),
            frame_step=int(step),
            dirname=_format_target_dirname(float(snapped)),
        ))
    return planned


# ---------------------------------------------------------------------------
# Agent-facing batch wrapper: report + collision check
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TIGenBatchReport:
    """Return value of :func:`generate_ti_batch_with_report`.

    All numeric fields preserve full floating-point precision (unlike
    ``ti_target_<cv>`` directory names, which round to 6 decimals).
    """

    workdirs: tuple[Path, ...]
    requested_targets_au: tuple[float, ...]
    snapped_targets_au: tuple[float, ...]
    snap_deltas_au: tuple[float, ...]
    steps: int

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-serializable view of the report.

        ``workdirs`` is converted to ``list[str]``; numeric fields are
        coerced to native Python types so ``json.dumps(report.to_dict())``
        never needs a custom encoder.  Direct ``dataclasses.asdict`` is
        intentionally **not** supported as a serialization boundary.
        """
        return {
            "workdirs": [str(p) for p in self.workdirs],
            "requested_targets_au": [float(v) for v in self.requested_targets_au],
            "snapped_targets_au": [float(v) for v in self.snapped_targets_au],
            "snap_deltas_au": [float(v) for v in self.snap_deltas_au],
            "steps": int(self.steps),
        }


def _preflight_file_inputs(
    *,
    inp_path: str | Path,
    xyz_path: str | Path,
    restart_path: str | Path,
    script_path: str | Path | None,
) -> None:
    """Raise ``FileNotFoundError`` for any missing file input.

    Runs before any planning / mkdir / filesystem write so callers see a
    uniform file-missing signal.  ``script_path=None`` is allowed (falls
    back to user config at generation time).
    """
    for label, path in (
        ("inp_path", inp_path),
        ("xyz_path", xyz_path),
        ("restart_path", restart_path),
    ):
        p = Path(path)
        if not p.is_file():
            raise FileNotFoundError(f"{label} not found or not a file: {p}")
    if script_path is not None:
        sp = Path(script_path)
        if not sp.is_file():
            raise FileNotFoundError(
                f"script_path not found or not a file: {sp}"
            )


def generate_ti_batch_with_report(
    inp_path: str | Path,
    xyz_path: str | Path,
    restart_path: str | Path,
    output_dir: str | Path,
    *,
    targets_au: list[float] | np.ndarray | None = None,
    time_range: dict[str, float | int] | None = None,
    steps: int = 10000,
    script_path: str | Path | None = None,
    colvar_id: int | None = None,
    overwrite: bool = False,
    verbose: bool = False,
) -> TIGenBatchReport:
    """Agent-safe batch generator with collision check and structured report.

    This wrapper is the backend for the ``ti_gen_batch`` agent task and
    (via ``workflows.scripts.run_ti_batch``) for CLI 422.  Differences
    from :func:`batch_generate_ti_workdirs`:

    - Uses an object-form ``time_range`` dict instead of three separate
      positional-ish args (more JSON Schema friendly).
    - Returns a :class:`TIGenBatchReport` so metrics can flow into the
      agent layer without reparsing directory names.
    - **Collision check**: if any planned target directory name already
      exists under ``output_dir`` *and* ``overwrite`` is False, raises
      :class:`ValueError` *before* writing anything.

    ``colvar_id`` selects which constraint the CV target is read from
    (``None`` = primary / first colvar; default).  This is pure
    plumbing — snapping, ``ti_target_<cv>`` dirname formatting, and the
    inp-modification rules are unchanged.

    ``overwrite`` (default False = agent-safe collision guard) only
    skips the pre-write collision check.  Per-target generation still
    goes through :func:`generate_ti_workdir`, which uses
    ``mkdir(exist_ok=True)`` and rewrites ``cMD.inp`` / ``init.xyz`` /
    ``script.sh`` in place.  It is a per-file overwrite — **no
    ``rmtree``, no directory wipe**; any other file already present in a
    colliding ``ti_target_*`` directory is left untouched.  CLI 422
    passes ``overwrite=True`` to preserve its historical
    overwrite-on-regeneration behavior.

    ``verbose`` shows a tqdm progress bar (CLI 422 passes True; the
    agent task does not expose it).

    Parameters
    ----------
    inp_path, xyz_path, restart_path
        Same as :func:`batch_generate_ti_workdirs`.
    output_dir : str or Path
        Parent directory for all work directories.
    targets_au : list of float or None
        Explicit CV targets in atomic units.  Mutually exclusive with
        ``time_range``.
    time_range : dict or None
        ``{"time_initial_fs": float, "time_final_fs": float,
        "n_points": int}``.  Mutually exclusive with ``targets_au``.
    steps : int
        MD steps per TI point (default 10000).
    script_path : str, Path or None
        Submission script to copy; falls back to ``KEY_CP2K_SCRIPT_PATH``
        in user config.

    Returns
    -------
    TIGenBatchReport

    Raises
    ------
    TIGenError
        Invalid target specification (both modes / neither mode / bad
        ``time_range``).
    ValueError
        A planned target directory already exists under ``output_dir``
        and ``overwrite`` is False (collision check).  No files are
        written in this case.
    FileNotFoundError
        ``inp_path`` / ``xyz_path`` / ``restart_path`` / ``script_path``
        missing on disk.
    """
    output_dir = Path(output_dir)

    # Preflight: verify all file inputs exist before any planning / filesystem
    # side effect.  This guarantees that missing-file errors surface uniformly
    # as FileNotFoundError (mapped to ``file_not_found`` by the contract),
    # regardless of which later step would have raised them.
    _preflight_file_inputs(
        inp_path=inp_path,
        xyz_path=xyz_path,
        restart_path=restart_path,
        script_path=script_path,
    )

    planned = _plan_ti_targets(
        inp_path=inp_path,
        xyz_path=xyz_path,
        restart_path=restart_path,
        targets_au=targets_au,
        time_range=time_range,
        colvar_id=colvar_id,
    )

    # Collision check — compare planned dirnames against existing ti_target_*
    # directory names.  Use exact string match (6-decimal dirname) as key.
    # overwrite=True only skips this pre-write guard; per-target generation
    # below still does a per-file overwrite via mkdir(exist_ok=True) — it
    # never removes or wipes a colliding directory.
    if not overwrite:
        existing_names = {p.name for p in output_dir.glob("ti_target_*") if p.is_dir()} \
            if output_dir.is_dir() else set()
        collisions = [pl.dirname for pl in planned if pl.dirname in existing_names]
        if collisions:
            raise ValueError(
                "TI workdir collision — the following target directories already "
                f"exist under {output_dir}: {sorted(collisions)}. "
                "Remove them, choose different targets, or pass overwrite=True "
                "to proceed."
            )

    # Pre-load restart / template tree / frames once for efficiency.
    # The template is parsed once into base_tree; each target deep-copies
    # it inside _render_ti_inp before editing (per-target discipline).
    restart = read_constraint_metadata_from_restart(restart_path)
    base_tree = _parse_sg_inp(inp_path)
    frames = _load_trajectory_cv(xyz_path, restart, colvar_id)

    iterator: list[_PlannedTarget] | object = planned
    if verbose:
        from tqdm import tqdm
        iterator = tqdm(planned, desc="TI workdirs", unit="point", ascii=" =")

    workdirs: list[Path] = []
    for pl in iterator:
        wd = generate_ti_workdir(
            inp_path, xyz_path, restart_path,
            target_au=pl.requested_cv,
            output_dir=output_dir,
            steps=steps,
            colvar_id=colvar_id,
            workdir_name=pl.dirname,
            script_path=script_path,
            _preloaded=frames,
            _restart=restart,
            _base_tree=base_tree,
        )
        workdirs.append(wd)

    return TIGenBatchReport(
        workdirs=tuple(workdirs),
        requested_targets_au=tuple(pl.requested_cv for pl in planned),
        snapped_targets_au=tuple(pl.snapped_cv for pl in planned),
        snap_deltas_au=tuple(pl.snap_delta for pl in planned),
        steps=steps,
    )


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
    _restart: ConstraintMetadata | None = None,
    _base_tree: dict | None = None,
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
    # Parse the template first (fail fast on a bad template, before any
    # trajectory loading or mkdir); batch callers pass the shared,
    # never-mutated base tree via _base_tree instead.
    base_tree = _base_tree if _base_tree is not None else _parse_sg_inp(inp_path)
    restart = _restart or read_constraint_metadata_from_restart(restart_path)
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

    # Write modified inp (per-target deep copy of the base tree)
    modified = _render_ti_inp(base_tree, snapped_cv, steps, colvar_id)
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

    restart = read_constraint_metadata_from_restart(restart_path)
    base_tree = _parse_sg_inp(inp_path)
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
            _base_tree=base_tree,
        )
        result.append(workdir)

    return result
