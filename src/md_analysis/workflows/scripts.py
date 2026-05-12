"""Work-directory generation workflow facade.

Exposes eight programmatic entry points covering the four
work-directory generators (Bader / TI / Potential / SP) in both
single-frame and batch modes:

- :func:`run_bader_single`     / :func:`run_bader_batch`
- :func:`run_ti_single`        / :func:`run_ti_batch`
- :func:`run_potential_single` / :func:`run_potential_batch`
- :func:`run_sp_single`        / :func:`run_sp_batch`

Each function **only prepares input files / work directories**. None
of these workflows submit PBS / SLURM / VASP / CP2K jobs; job
submission belongs in the operator's queue tooling. The
``BaderGen.generate_potcar=True`` path invokes ``vaspkit 103``
locally to write a POTCAR file but does not enqueue any job.

Conventions:

- Path arguments are normalised to :class:`pathlib.Path` at the
  workflow boundary.
- Single workflows expose ``artifacts = {"workdir": <dir>}``; batch
  workflows expose ``artifacts = {"workdir_<i>": <dir_i>, ...}``.
  ``require_artifacts_exist`` accepts directories.
- Batch ``metadata`` always carries ``n_successful`` /
  ``n_skipped`` / ``n_failed`` counts plus a JSON-serialisable
  ``workdir_paths`` list, even when the underlying generator does
  not currently support skipping or partial failure (today they
  succeed-all-or-raise, so ``n_skipped == n_failed == 0``).
- Strongly-typed batch reports are surfaced on ``extra`` when the
  underlying API provides one (Bader / TI / SP); for Potential no
  ``*BatchReport`` exists upstream so ``extra=None`` and all
  metadata is in the dict.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Iterable

from .models import WorkflowResult

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _resolve_single_frame(
    xyz_path: Path,
    *,
    mode: str,
    frame: int,
    time_fs: float | None,
    time_tol_fs: float,
) -> tuple[int, Any, list[str]]:
    """Wrap :func:`scripts._frame_selector.resolve_single_frame`.

    Returns ``(frame_index, atoms, warnings)``.
    """
    from ..scripts._frame_selector import resolve_single_frame

    return resolve_single_frame(
        xyz_path,
        mode=mode,
        frame=frame,
        time_fs=time_fs,
        time_tol_fs=time_tol_fs,
    )


def _build_batch_metadata(
    *,
    workdirs: list[Path] | tuple[Path, ...],
    extra_fields: dict[str, Any],
) -> dict[str, Any]:
    """Assemble the standard batch metadata block.

    Every batch workflow surfaces ``n_successful`` / ``n_skipped`` /
    ``n_failed`` counts plus a JSON-friendly ``workdir_paths`` list
    so downstream callers can inspect outcomes without re-walking
    ``artifacts``. ``extra_fields`` carries per-workflow specifics
    (mode, frame ranges, target counts, ...).
    """
    return {
        "n_successful": len(workdirs),
        "n_skipped": 0,
        "n_failed": 0,
        "workdir_paths": [str(w) for w in workdirs],
        **extra_fields,
    }


def _flatten_workdirs_to_artifacts(
    workdirs: Iterable[Path],
) -> dict[str, Path]:
    """Pack a sequence of generated workdirs into the flat artifact map.

    Keys are stable, zero-padded ``workdir_<i>`` so ``require_artifacts_exist``
    failures can be cross-referenced with ``metadata["workdir_paths"]``
    by index.
    """
    return {f"workdir_{i}": Path(w) for i, w in enumerate(workdirs)}


# ===========================================================================
# Bader
# ===========================================================================


def run_bader_single(
    *,
    xyz_path: Path | str,
    cell_abc: tuple[float, float, float] | list[float],
    output_dir: Path | str,
    frame: int = 0,
    mode: str = "index",
    time_fs: float | None = None,
    time_tol_fs: float = 1e-6,
    workdir_name: str = "bader",
    script_path: Path | str | None = None,
    element_order: tuple[str, ...] | list[str] | None = None,
    generate_potcar: bool = True,
    direct: bool = True,
) -> WorkflowResult:
    """Generate one VASP Bader single-point work directory.

    Resolves a single frame from ``xyz_path`` via
    :func:`scripts._frame_selector.resolve_single_frame`, applies the
    user-supplied ``cell_abc``, and delegates to
    :func:`scripts.generate_bader_workdir`. Only directory creation
    and (optionally) ``vaspkit 103`` POTCAR generation are performed;
    no job submission.
    """
    from ..scripts import generate_bader_workdir

    xyz_p = Path(xyz_path)
    out_dir = Path(output_dir)
    script_p = Path(script_path) if script_path is not None else None

    frame_idx, atoms, warnings = _resolve_single_frame(
        xyz_p,
        mode=mode,
        frame=frame,
        time_fs=time_fs,
        time_tol_fs=time_tol_fs,
    )
    for msg in warnings:
        logger.warning(msg)

    atoms.set_cell(tuple(cell_abc))
    atoms.set_pbc(True)

    workdir = generate_bader_workdir(
        atoms,
        out_dir,
        script_path=script_p,
        workdir_name=workdir_name,
        frame=frame_idx,
        source=str(xyz_p),
        element_order=tuple(element_order) if element_order is not None else None,
        generate_potcar=generate_potcar,
        direct=direct,
    )

    return WorkflowResult(
        name="bader_single",
        output_dir=out_dir,
        artifacts={"workdir": Path(workdir)},
        metadata={
            "frame_index": frame_idx,
            "mode": mode,
            "time_fs": time_fs,
            "workdir_name": workdir_name,
            "generate_potcar": generate_potcar,
            "source_xyz": str(xyz_p),
            "frame_warnings": warnings,
        },
    )


def run_bader_batch(
    *,
    xyz_path: Path | str,
    cell_abc: tuple[float, float, float] | list[float],
    output_dir: Path | str,
    mode: str = "index",
    frame_start: int = 0,
    frame_end: int | None = None,
    frame_step: int = 1,
    time_start_fs: float | None = None,
    time_end_fs: float | None = None,
    time_step_fs: float | None = None,
    script_path: Path | str | None = None,
    element_order: tuple[str, ...] | list[str] | None = None,
    generate_potcar: bool = True,
    direct: bool = True,
    verbose: bool = False,
) -> WorkflowResult:
    """Batch-generate VASP Bader work directories from a trajectory.

    Delegates to :func:`scripts.generate_bader_batch_with_report` so
    the ``BaderGenBatchReport`` (workdirs + per-frame step/time
    metadata) is preserved on ``extra``. ``n_skipped`` / ``n_failed``
    are always zero with the current generator (it succeeds-all or
    raises) but are surfaced for contract stability.
    """
    from ..scripts.BaderGen import generate_bader_batch_with_report

    xyz_p = Path(xyz_path)
    out_dir = Path(output_dir)
    script_p = Path(script_path) if script_path is not None else None

    report = generate_bader_batch_with_report(
        xyz_p,
        tuple(cell_abc),
        out_dir,
        mode=mode,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        time_start_fs=time_start_fs,
        time_end_fs=time_end_fs,
        time_step_fs=time_step_fs,
        script_path=script_p,
        element_order=list(element_order) if element_order is not None else None,
        generate_potcar=generate_potcar,
        direct=direct,
        verbose=verbose,
    )

    workdirs = [Path(w) for w in report.workdirs]
    metadata = _build_batch_metadata(
        workdirs=workdirs,
        extra_fields={
            "mode": mode,
            "frame_indices": list(report.frame_indices),
            "steps": list(report.steps),
            "times_fs": list(report.times_fs),
            "n_frames": report.n_frames,
            "generate_potcar": report.generate_potcar,
            "source_xyz": str(xyz_p),
        },
    )
    return WorkflowResult(
        name="bader_batch",
        output_dir=out_dir,
        artifacts=_flatten_workdirs_to_artifacts(workdirs),
        metadata=metadata,
        extra=report,
    )


# ===========================================================================
# TI
# ===========================================================================


def run_ti_single(
    *,
    inp_path: Path | str,
    xyz_path: Path | str,
    restart_path: Path | str,
    target_au: float,
    output_dir: Path | str,
    steps: int = 10000,
    colvar_id: int | None = None,
    workdir_name: str | None = None,
    script_path: Path | str | None = None,
) -> WorkflowResult:
    """Generate one CP2K constrained-MD work directory for TI sampling.

    Snaps ``target_au`` to the nearest trajectory frame and writes
    a ``cMD.inp`` + ``init.xyz`` (+ ``script.sh``) under
    ``output_dir/<workdir_name>``. No job submission.
    """
    from ..scripts import generate_ti_workdir

    inp_p = Path(inp_path)
    xyz_p = Path(xyz_path)
    restart_p = Path(restart_path)
    out_dir = Path(output_dir)
    script_p = Path(script_path) if script_path is not None else None

    workdir = generate_ti_workdir(
        inp_p,
        xyz_p,
        restart_p,
        target_au,
        out_dir,
        steps=steps,
        colvar_id=colvar_id,
        workdir_name=workdir_name,
        script_path=script_p,
    )

    return WorkflowResult(
        name="ti_single",
        output_dir=out_dir,
        artifacts={"workdir": Path(workdir)},
        metadata={
            "requested_target_au": float(target_au),
            "steps": steps,
            "colvar_id": colvar_id,
            "workdir_name": workdir.name,
            "source_inp": str(inp_p),
            "source_xyz": str(xyz_p),
            "source_restart": str(restart_p),
        },
    )


def run_ti_batch(
    *,
    inp_path: Path | str,
    xyz_path: Path | str,
    restart_path: Path | str,
    output_dir: Path | str,
    targets_au: list[float] | None = None,
    time_range: dict[str, float | int] | None = None,
    steps: int = 10000,
    script_path: Path | str | None = None,
) -> WorkflowResult:
    """Batch-generate CP2K constrained-MD work directories for TI.

    Exactly one of ``targets_au`` / ``time_range`` must be provided
    (the underlying generator enforces this). The wrapper performs
    a pre-write collision check and refuses to overwrite existing
    target directories under ``output_dir``.
    """
    from ..scripts.TIGen import generate_ti_batch_with_report

    inp_p = Path(inp_path)
    xyz_p = Path(xyz_path)
    restart_p = Path(restart_path)
    out_dir = Path(output_dir)
    script_p = Path(script_path) if script_path is not None else None

    report = generate_ti_batch_with_report(
        inp_p,
        xyz_p,
        restart_p,
        out_dir,
        targets_au=targets_au,
        time_range=time_range,
        steps=steps,
        script_path=script_p,
    )

    workdirs = [Path(w) for w in report.workdirs]
    metadata = _build_batch_metadata(
        workdirs=workdirs,
        extra_fields={
            "requested_targets_au": list(report.requested_targets_au),
            "snapped_targets_au": list(report.snapped_targets_au),
            "snap_deltas_au": list(report.snap_deltas_au),
            "steps": report.steps,
            "input_source": "targets_au" if targets_au is not None else "time_range",
            "source_inp": str(inp_p),
            "source_xyz": str(xyz_p),
            "source_restart": str(restart_p),
        },
    )
    return WorkflowResult(
        name="ti_batch",
        output_dir=out_dir,
        artifacts=_flatten_workdirs_to_artifacts(workdirs),
        metadata=metadata,
        extra=report,
    )


# ===========================================================================
# Potential (Hartree SP)
# ===========================================================================


def run_potential_single(
    *,
    xyz_path: Path | str,
    output_dir: Path | str,
    inp_template_path: Path | str | None = None,
    cell_abc: tuple[float, float, float] | list[float] | None = None,
    frame: int = 0,
    mode: str = "index",
    time_fs: float | None = None,
    time_tol_fs: float = 1e-6,
    workdir_name: str = "potential",
    script_path: Path | str | None = None,
) -> WorkflowResult:
    """Generate one CP2K Hartree-potential SP work directory.

    When ``cell_abc`` is ``None`` the cell from the resolved trajectory
    frame is used. When ``inp_template_path`` is ``None`` the workflow
    delegates to ``scripts.generate_potential_workdir`` which in turn
    falls back to the persisted ``KEY_SP_INP_TEMPLATE_PATH`` config
    entry; this is a documented persistent-config dependency.
    """
    from ..scripts import generate_potential_workdir

    xyz_p = Path(xyz_path)
    out_dir = Path(output_dir)
    inp_p = Path(inp_template_path) if inp_template_path is not None else None
    script_p = Path(script_path) if script_path is not None else None

    frame_idx, atoms, warnings = _resolve_single_frame(
        xyz_p,
        mode=mode,
        frame=frame,
        time_fs=time_fs,
        time_tol_fs=time_tol_fs,
    )
    for msg in warnings:
        logger.warning(msg)

    workdir = generate_potential_workdir(
        atoms,
        out_dir,
        inp_template_path=inp_p,
        cell_abc=tuple(cell_abc) if cell_abc is not None else None,
        script_path=script_p,
        workdir_name=workdir_name,
        frame=frame_idx,
        source=str(xyz_p),
    )

    return WorkflowResult(
        name="potential_single",
        output_dir=out_dir,
        artifacts={"workdir": Path(workdir)},
        metadata={
            "frame_index": frame_idx,
            "mode": mode,
            "time_fs": time_fs,
            "workdir_name": workdir_name,
            "source_xyz": str(xyz_p),
            "inp_template_path_resolved_from": (
                "argument" if inp_p is not None else "config_or_default"
            ),
            "frame_warnings": warnings,
        },
    )


def run_potential_batch(
    *,
    xyz_path: Path | str,
    cell_abc: tuple[float, float, float] | list[float],
    output_dir: Path | str,
    inp_template_path: Path | str | None = None,
    mode: str = "index",
    frame_start: int = 0,
    frame_end: int | None = None,
    frame_step: int = 1,
    time_start_fs: float | None = None,
    time_end_fs: float | None = None,
    time_step_fs: float | None = None,
    script_path: Path | str | None = None,
    verbose: bool = False,
) -> WorkflowResult:
    """Batch-generate CP2K Hartree-potential SP work directories.

    The underlying :func:`scripts.batch_generate_potential_workdirs`
    does not yet have a ``*_with_report`` companion (PotentialGen is
    CLI-only at the agent layer), so the workflow synthesises the
    standard batch metadata from the returned ``list[Path]`` itself.
    ``extra`` is therefore ``None`` for this workflow.
    """
    from ..scripts import batch_generate_potential_workdirs

    xyz_p = Path(xyz_path)
    out_dir = Path(output_dir)
    inp_p = Path(inp_template_path) if inp_template_path is not None else None
    script_p = Path(script_path) if script_path is not None else None

    workdir_list = batch_generate_potential_workdirs(
        xyz_p,
        tuple(cell_abc),
        out_dir,
        inp_template_path=inp_p,
        mode=mode,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        time_start_fs=time_start_fs,
        time_end_fs=time_end_fs,
        time_step_fs=time_step_fs,
        script_path=script_p,
        verbose=verbose,
    )

    workdirs = [Path(w) for w in workdir_list]
    metadata = _build_batch_metadata(
        workdirs=workdirs,
        extra_fields={
            "mode": mode,
            "source_xyz": str(xyz_p),
            "inp_template_path_resolved_from": (
                "argument" if inp_p is not None else "config_or_default"
            ),
        },
    )
    return WorkflowResult(
        name="potential_batch",
        output_dir=out_dir,
        artifacts=_flatten_workdirs_to_artifacts(workdirs),
        metadata=metadata,
        extra=None,
    )


# ===========================================================================
# SP (DeePMD training-data SP)
# ===========================================================================


def run_sp_single(
    *,
    xyz_path: Path | str,
    cell_abc: tuple[float, float, float] | list[float],
    output_dir: Path | str,
    inp_template_path: Path | str | None = None,
    frame: int = 0,
    mode: str = "index",
    time_fs: float | None = None,
    time_tol_fs: float = 1e-6,
    workdir_name: str = "sp",
    script_path: Path | str | None = None,
) -> WorkflowResult:
    """Generate one CP2K SP work directory for DeePMD training data.

    Distinct from :func:`run_potential_single` even though both write
    a CP2K SP input — see ``scripts/CLAUDE.md`` for the contract
    differences (config key ``KEY_DP_SP_INP_TEMPLATE_PATH`` and
    ``sp_t*_i*`` directory prefix for batches).
    """
    from ..scripts import generate_sp_workdir

    xyz_p = Path(xyz_path)
    out_dir = Path(output_dir)
    inp_p = Path(inp_template_path) if inp_template_path is not None else None
    script_p = Path(script_path) if script_path is not None else None

    frame_idx, atoms, warnings = _resolve_single_frame(
        xyz_p,
        mode=mode,
        frame=frame,
        time_fs=time_fs,
        time_tol_fs=time_tol_fs,
    )
    for msg in warnings:
        logger.warning(msg)

    workdir = generate_sp_workdir(
        atoms,
        out_dir,
        cell_abc=tuple(cell_abc),
        inp_template_path=inp_p,
        script_path=script_p,
        workdir_name=workdir_name,
        frame=frame_idx,
        source=str(xyz_p),
    )

    return WorkflowResult(
        name="sp_single",
        output_dir=out_dir,
        artifacts={"workdir": Path(workdir)},
        metadata={
            "frame_index": frame_idx,
            "mode": mode,
            "time_fs": time_fs,
            "workdir_name": workdir_name,
            "source_xyz": str(xyz_p),
            "inp_template_path_resolved_from": (
                "argument" if inp_p is not None else "config_or_default"
            ),
            "frame_warnings": warnings,
        },
    )


def run_sp_batch(
    *,
    xyz_path: Path | str,
    cell_abc: tuple[float, float, float] | list[float],
    output_dir: Path | str,
    inp_template_path: Path | str | None = None,
    mode: str = "index",
    frame_start: int = 0,
    frame_end: int | None = None,
    frame_step: int = 1,
    time_start_fs: float | None = None,
    time_end_fs: float | None = None,
    time_step_fs: float | None = None,
    script_path: Path | str | None = None,
    verbose: bool = False,
) -> WorkflowResult:
    """Batch-generate CP2K DeePMD-SP work directories from a trajectory."""
    from ..scripts.SpGen import generate_sp_batch_with_report

    xyz_p = Path(xyz_path)
    out_dir = Path(output_dir)
    inp_p = Path(inp_template_path) if inp_template_path is not None else None
    script_p = Path(script_path) if script_path is not None else None

    report = generate_sp_batch_with_report(
        xyz_p,
        tuple(cell_abc),
        out_dir,
        inp_template_path=inp_p,
        mode=mode,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        time_start_fs=time_start_fs,
        time_end_fs=time_end_fs,
        time_step_fs=time_step_fs,
        script_path=script_p,
        verbose=verbose,
    )

    workdirs = [Path(w) for w in report.workdirs]
    metadata = _build_batch_metadata(
        workdirs=workdirs,
        extra_fields={
            "mode": mode,
            "frame_indices": list(report.frame_indices),
            "steps": list(report.steps),
            "times_fs": list(report.times_fs),
            "n_frames": report.n_frames,
            "source_xyz": str(xyz_p),
            "inp_template_path_resolved_from": (
                "argument" if inp_p is not None else "config_or_default"
            ),
        },
    )
    return WorkflowResult(
        name="sp_batch",
        output_dir=out_dir,
        artifacts=_flatten_workdirs_to_artifacts(workdirs),
        metadata=metadata,
        extra=report,
    )
