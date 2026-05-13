"""Scripts / Tools command classes (Bader 411-412, TI 421-422, Potential 431-432, SpGen 441-442)."""

from __future__ import annotations

from pathlib import Path

from ._framework import MenuCommand, lazy_import
from ._params import (
    BoolParam,
    ConditionalParam,
    DisplayAction,
    FloatParam,
    IntParam,
    K,
    StrParam,
    cell_abc,
    cp2k_script,
    dp_sp_inp_template,
    frame_mode,
    gen_potcar,
    single_time_fs,
    sp_inp_template,
    time_end_fs,
    time_start_fs,
    time_step_fs,
    vasp_script,
)


# ---------------------------------------------------------------------------
# Helpers for frame selection mode (shared by Bader/Potential/SpGen commands)
# ---------------------------------------------------------------------------

_IS_INDEX_MODE = lambda ctx: ctx[K.FRAME_MODE] == "index"
_IS_TIME_MODE = lambda ctx: ctx[K.FRAME_MODE] == "time"


def _resolve_single_frame_from_ctx(ctx: dict):
    """Resolve (frame_idx, atoms) from CLI ctx using frame_mode.

    Logs any warnings returned by ``resolve_single_frame`` at WARNING level.
    """
    import logging as _logging
    resolve_single_frame = lazy_import(
        "md_analysis.scripts._frame_selector", "resolve_single_frame",
    )

    mode = ctx[K.FRAME_MODE]
    if mode == "index":
        idx, atoms, warnings = resolve_single_frame(
            ctx[K.XYZ], mode="index", frame=ctx[K.FRAME],
        )
    else:
        idx, atoms, warnings = resolve_single_frame(
            ctx[K.XYZ], mode="time", time_fs=ctx[K.SINGLE_TIME_FS],
        )

    _logger = _logging.getLogger("md_analysis.cli")
    for w in warnings:
        _logger.warning(w)
        print(f"  WARNING: {w}")

    return idx, atoms


def _frame_selection_params() -> tuple:
    """Return the common frame-selection param block for Batch commands.

    Order: mode prompt → conditional (index triplet) → conditional (time triplet).
    """
    return (
        frame_mode,
        ConditionalParam(
            IntParam(K.FRAME_START, "Frame start (0-based)", default=0),
            _IS_INDEX_MODE,
        ),
        ConditionalParam(
            IntParam(K.FRAME_END, "Frame end (exclusive, empty=all)", default=None),
            _IS_INDEX_MODE,
        ),
        ConditionalParam(
            IntParam(K.FRAME_STEP, "Frame step", default=1),
            _IS_INDEX_MODE,
        ),
        ConditionalParam(time_start_fs, _IS_TIME_MODE),
        ConditionalParam(time_end_fs, _IS_TIME_MODE),
        ConditionalParam(time_step_fs, _IS_TIME_MODE),
    )


def _single_frame_params() -> tuple:
    """Return the common frame-selection param block for Single commands."""
    return (
        frame_mode,
        ConditionalParam(
            IntParam(K.FRAME, "Frame number (0-based)", default=0),
            _IS_INDEX_MODE,
        ),
        ConditionalParam(single_time_fs, _IS_TIME_MODE),
    )


def _batch_frame_kwargs_from_ctx(ctx: dict) -> dict:
    """Build kwargs for batch_generate_*_workdirs from CLI ctx.

    Always passes mode + both parameter sets; FrameSelection internally
    ignores the unused set based on mode.
    """
    return {
        "mode": ctx[K.FRAME_MODE],
        "frame_start": ctx[K.FRAME_START],
        "frame_end": ctx[K.FRAME_END],
        "frame_step": ctx[K.FRAME_STEP],
        "time_start_fs": ctx[K.TIME_START_FS]
            if ctx[K.FRAME_MODE] == "time" else None,
        "time_end_fs": ctx[K.TIME_END_FS]
            if ctx[K.FRAME_MODE] == "time" else None,
        "time_step_fs": ctx[K.TIME_STEP_FS]
            if ctx[K.FRAME_MODE] == "time" else None,
    }


def _single_frame_workflow_kwargs(ctx: dict) -> dict:
    """Build kwargs for run_*_single workflows from CLI ctx.

    The workflow facades take ``mode`` plus ``frame`` / ``time_fs``
    rather than the pre-resolved ``atoms`` object the legacy
    ``generate_*_workdir`` functions used. Frame resolution lives
    inside the workflow.
    """
    mode = ctx[K.FRAME_MODE]
    kwargs: dict = {"mode": mode}
    if mode == "index":
        kwargs["frame"] = ctx[K.FRAME]
    else:
        kwargs["time_fs"] = ctx[K.SINGLE_TIME_FS]
    return kwargs
from ._prompt import (
    prompt_choice,
    prompt_float,
    prompt_int,
    prompt_str,
    prompt_str_required,
)


def _print_trajectory_info(xyz_path: str) -> None:
    """Peek at XYZ trajectory and print frame/step/time metadata."""
    p = Path(xyz_path)
    if not p.is_file():
        return

    with open(p) as fh:
        first_line = fh.readline().strip()
        try:
            natoms = int(first_line)
        except ValueError:
            return
        total_lines = sum(1 for _ in fh) + 1
    block_size = natoms + 2
    total_frames = total_lines // block_size

    from ase.io import iread

    frames_meta: list[tuple[int, float]] = []
    for idx, atoms in enumerate(iread(str(p), index=":")):
        if idx >= 2:
            break
        frames_meta.append((
            int(atoms.info.get("i", idx)),
            float(atoms.info.get("time", 0.0)),
        ))

    print(f"\n  Trajectory: {total_frames} frames, {natoms} atoms/frame")
    if len(frames_meta) >= 2:
        step_interval = frames_meta[1][0] - frames_meta[0][0]
        time_interval = frames_meta[1][1] - frames_meta[0][1]
        if step_interval > 0:
            dt = time_interval / step_interval
            print(f"  Frame interval: {step_interval} MD steps, "
                  f"{time_interval:.1f} fs/frame (dt = {dt:.1f} fs/step)")
    print()


def _resolve_cp2k_script_path() -> str | None:
    """Prompt for CP2K submission script path with config default."""
    from ..config import KEY_CP2K_SCRIPT_PATH, get_config

    default_script = get_config(KEY_CP2K_SCRIPT_PATH)
    return prompt_str("Submission script path", default=default_script)


class BaderSingleCmd(MenuCommand):

    params = (
        StrParam(K.XYZ, "XYZ trajectory file (e.g. md-pos-1.xyz)", required=True),
        DisplayAction(lambda ctx: _print_trajectory_info(ctx[K.XYZ])),
        cell_abc,
        *_single_frame_params(),
        StrParam(K.OUTDIR, "Output directory", default="."),
        StrParam(K.WORKDIR_NAME, "Work directory name", default="bader"),
        vasp_script,
        gen_potcar,
    )

    def execute(self, ctx: dict) -> None:
        run_bader_single = lazy_import(
            "md_analysis.workflows.scripts", "run_bader_single",
        )
        single_kwargs = _single_frame_workflow_kwargs(ctx)
        result = run_bader_single(
            xyz_path=ctx[K.XYZ],
            cell_abc=ctx[K.CELL_ABC],
            output_dir=ctx[K.OUTDIR],
            workdir_name=ctx[K.WORKDIR_NAME],
            script_path=ctx[K.SCRIPT_PATH],
            generate_potcar=ctx[K.GEN_POTCAR],
            **single_kwargs,
        )
        workdir = result.artifacts["workdir"]
        print(f"\n Bader work directory created: {workdir}")
        contents = sorted(p.name for p in workdir.iterdir())
        print(f"  Contents: {', '.join(contents)}")


class BaderBatchCmd(MenuCommand):

    params = (
        StrParam(K.XYZ, "XYZ trajectory file (e.g. md-pos-1.xyz)", required=True),
        DisplayAction(lambda ctx: _print_trajectory_info(ctx[K.XYZ])),
        cell_abc,
        *_frame_selection_params(),
        StrParam(K.OUTDIR, "Output directory", default="."),
        vasp_script,
        gen_potcar,
    )

    def execute(self, ctx: dict) -> None:
        run_bader_batch = lazy_import(
            "md_analysis.workflows.scripts", "run_bader_batch",
        )
        kwargs = _batch_frame_kwargs_from_ctx(ctx)
        result = run_bader_batch(
            xyz_path=ctx[K.XYZ],
            cell_abc=ctx[K.CELL_ABC],
            output_dir=ctx[K.OUTDIR],
            script_path=ctx[K.SCRIPT_PATH],
            generate_potcar=ctx[K.GEN_POTCAR],
            verbose=True,
            **kwargs,
        )
        workdirs = [
            v for k, v in result.artifacts.items() if k.startswith("workdir_")
        ]
        print(f"\n Created {len(workdirs)} Bader work directories:")
        for d in workdirs:
            print(f"  {d}")


# ---------------------------------------------------------------------------
# TI commands (42x)
# ---------------------------------------------------------------------------

def _print_sg_cv_info(restart_path: str, xyz_path: str) -> None:
    """Display SG trajectory CV range and frame info for TI target selection."""
    parse_colvar_restart = lazy_import(
        "md_analysis.utils.formats.cp2k.colvar", "parse_colvar_restart",
    )
    try:
        restart = parse_colvar_restart(restart_path)
    except Exception as exc:
        print(f"  (Could not parse restart: {exc})")
        return

    cv = restart.colvars.primary
    au_time_to_fs = lazy_import("md_analysis.utils.constants", "AU_TIME_TO_FS")
    dt_au = restart.timestep_fs / au_time_to_fs
    growth_per_step = cv.target_growth_au * dt_au

    print(f"\n  SG info:")
    print(f"    Timestep:     {restart.timestep_fs} fs")
    print(f"    CV target:    {cv.target_au:.6f} a.u. (at step {restart.step_start})")
    print(f"    CV growth:    {growth_per_step:.6e} a.u./step")
    print(f"    Cell:         {restart.cell_abc_ang[0]:.4f} x "
          f"{restart.cell_abc_ang[1]:.4f} x {restart.cell_abc_ang[2]:.4f} A")
    print()


class TISingleCmd(MenuCommand):
    """Generate one TI constrained-MD work directory."""


    def _collect_all_params(self) -> dict:
        print()
        ctx: dict = {}
        ctx[K.RESTART_PATH] = prompt_str_required(
            "SG restart file (e.g. slowgrowth-1.restart)"
        )
        ctx[K.INP_PATH] = prompt_str_required(
            "SG input file (e.g. sg.inp)"
        )
        ctx[K.XYZ] = prompt_str_required(
            "SG trajectory file (e.g. slowgrowth-pos-1.xyz)"
        )
        _print_trajectory_info(ctx[K.XYZ])
        _print_sg_cv_info(ctx[K.RESTART_PATH], ctx[K.XYZ])

        ctx[K.TARGET_AU] = prompt_float("Target CV value (a.u.)", default=0.0)
        ctx[K.STEPS] = prompt_int("MD steps for constrained-MD", default=10000) or 10000
        ctx[K.COLVAR_ID] = prompt_int("Colvar ID (empty=primary)", default=None)
        ctx[K.OUTDIR] = prompt_str("Output directory", default=".") or "."
        ctx[K.WORKDIR_NAME] = prompt_str(
            "Work directory name (empty=auto)", default=None,
        )
        ctx[K.SCRIPT_PATH] = _resolve_cp2k_script_path()
        return ctx

    def execute(self, ctx: dict) -> None:
        run_ti_single = lazy_import(
            "md_analysis.workflows.scripts", "run_ti_single",
        )
        result = run_ti_single(
            inp_path=ctx[K.INP_PATH],
            xyz_path=ctx[K.XYZ],
            restart_path=ctx[K.RESTART_PATH],
            target_au=ctx[K.TARGET_AU],
            output_dir=ctx[K.OUTDIR],
            steps=ctx[K.STEPS],
            colvar_id=ctx[K.COLVAR_ID],
            workdir_name=ctx[K.WORKDIR_NAME],
            script_path=ctx[K.SCRIPT_PATH],
        )
        workdir = result.artifacts["workdir"]
        print(f"\n TI work directory created: {workdir}")
        contents = sorted(p.name for p in workdir.iterdir())
        print(f"  Contents: {', '.join(contents)}")


class TIBatchCmd(MenuCommand):
    """Batch-generate TI constrained-MD work directories."""


    def _collect_all_params(self) -> dict:
        print()
        ctx: dict = {}
        ctx[K.RESTART_PATH] = prompt_str_required(
            "SG restart file (e.g. slowgrowth-1.restart)"
        )
        ctx[K.INP_PATH] = prompt_str_required(
            "SG input file (e.g. sg.inp)"
        )
        ctx[K.XYZ] = prompt_str_required(
            "SG trajectory file (e.g. slowgrowth-pos-1.xyz)"
        )
        _print_trajectory_info(ctx[K.XYZ])
        _print_sg_cv_info(ctx[K.RESTART_PATH], ctx[K.XYZ])

        mode = prompt_choice(
            "Target specification mode",
            ["time", "values"],
            default="time",
        )

        if mode == "time":
            ctx[K.TIME_INITIAL_FS] = prompt_float(
                "Initial time (fs)", default=0.0,
            )
            ctx[K.TIME_FINAL_FS] = prompt_float(
                "Final time (fs)", default=0.0,
            )
            ctx[K.N_POINTS] = prompt_int(
                "Number of TI points", default=10,
            ) or 10
            ctx[K.TARGETS_AU] = None
        else:
            raw = prompt_str_required(
                "Target CV values in a.u. (space-separated)"
            )
            ctx[K.TARGETS_AU] = [float(x) for x in raw.split()]
            ctx[K.TIME_INITIAL_FS] = None
            ctx[K.TIME_FINAL_FS] = None
            ctx[K.N_POINTS] = None

        ctx[K.STEPS] = prompt_int("MD steps for constrained-MD", default=10000) or 10000
        ctx[K.COLVAR_ID] = prompt_int("Colvar ID (empty=primary)", default=None)
        ctx[K.OUTDIR] = prompt_str("Output directory", default=".") or "."
        ctx[K.SCRIPT_PATH] = _resolve_cp2k_script_path()
        return ctx

    def execute(self, ctx: dict) -> None:
        batch = lazy_import("md_analysis.scripts", "batch_generate_ti_workdirs")
        dirs = batch(
            ctx[K.INP_PATH],
            ctx[K.XYZ],
            ctx[K.RESTART_PATH],
            ctx[K.OUTDIR],
            targets_au=ctx[K.TARGETS_AU],
            time_initial_fs=ctx[K.TIME_INITIAL_FS],
            time_final_fs=ctx[K.TIME_FINAL_FS],
            n_points=ctx[K.N_POINTS],
            steps=ctx[K.STEPS],
            colvar_id=ctx[K.COLVAR_ID],
            script_path=ctx[K.SCRIPT_PATH],
            verbose=True,
        )
        print(f"\n Created {len(dirs)} TI work directories:")
        for d in dirs:
            print(f"  {d}")


# ---------------------------------------------------------------------------
# Potential SP commands (43x)
# ---------------------------------------------------------------------------

class PotentialSingleCmd(MenuCommand):
    """Generate one SP potential work directory."""

    params = (
        StrParam(K.XYZ, "XYZ trajectory file (e.g. md-pos-1.xyz)", required=True),
        DisplayAction(lambda ctx: _print_trajectory_info(ctx[K.XYZ])),
        cell_abc,
        *_single_frame_params(),
        sp_inp_template,
        StrParam(K.OUTDIR, "Output directory", default="."),
        StrParam(K.WORKDIR_NAME, "Work directory name", default="potential"),
        cp2k_script,
    )

    def execute(self, ctx: dict) -> None:
        run_potential_single = lazy_import(
            "md_analysis.workflows.scripts", "run_potential_single",
        )
        single_kwargs = _single_frame_workflow_kwargs(ctx)
        result = run_potential_single(
            xyz_path=ctx[K.XYZ],
            output_dir=ctx[K.OUTDIR],
            cell_abc=ctx[K.CELL_ABC],
            inp_template_path=ctx[K.INP_TEMPLATE],
            workdir_name=ctx[K.WORKDIR_NAME],
            script_path=ctx[K.SCRIPT_PATH],
            **single_kwargs,
        )
        workdir = result.artifacts["workdir"]
        print(f"\n Potential work directory created: {workdir}")
        contents = sorted(p.name for p in workdir.iterdir())
        print(f"  Contents: {', '.join(contents)}")


class PotentialBatchCmd(MenuCommand):
    """Batch-generate SP potential work directories."""

    params = (
        StrParam(K.XYZ, "XYZ trajectory file (e.g. md-pos-1.xyz)", required=True),
        DisplayAction(lambda ctx: _print_trajectory_info(ctx[K.XYZ])),
        cell_abc,
        *_frame_selection_params(),
        sp_inp_template,
        StrParam(K.OUTDIR, "Output directory", default="."),
        cp2k_script,
    )

    def execute(self, ctx: dict) -> None:
        run_potential_batch = lazy_import(
            "md_analysis.workflows.scripts", "run_potential_batch",
        )
        kwargs = _batch_frame_kwargs_from_ctx(ctx)
        result = run_potential_batch(
            xyz_path=ctx[K.XYZ],
            cell_abc=ctx[K.CELL_ABC],
            output_dir=ctx[K.OUTDIR],
            inp_template_path=ctx[K.INP_TEMPLATE],
            script_path=ctx[K.SCRIPT_PATH],
            verbose=True,
            **kwargs,
        )
        workdirs = [
            v for k, v in result.artifacts.items() if k.startswith("workdir_")
        ]
        print(f"\n Created {len(workdirs)} potential work directories:")
        for d in workdirs:
            print(f"  {d}")


# ---------------------------------------------------------------------------
# SpGen SP commands for DeePMD training (44x)
# ---------------------------------------------------------------------------


class SpGenSingleCmd(MenuCommand):
    """Generate one SP work directory for DeePMD training."""

    params = (
        StrParam(K.XYZ, "XYZ trajectory file (e.g. md-pos-1.xyz)", required=True),
        DisplayAction(lambda ctx: _print_trajectory_info(ctx[K.XYZ])),
        cell_abc,
        *_single_frame_params(),
        dp_sp_inp_template,
        StrParam(K.OUTDIR, "Output directory", default="."),
        StrParam(K.WORKDIR_NAME, "Work directory name", default="sp"),
        cp2k_script,
    )

    def execute(self, ctx: dict) -> None:
        run_sp_single = lazy_import(
            "md_analysis.workflows.scripts", "run_sp_single",
        )
        single_kwargs = _single_frame_workflow_kwargs(ctx)
        result = run_sp_single(
            xyz_path=ctx[K.XYZ],
            cell_abc=ctx[K.CELL_ABC],
            output_dir=ctx[K.OUTDIR],
            inp_template_path=ctx[K.INP_TEMPLATE],
            workdir_name=ctx[K.WORKDIR_NAME],
            script_path=ctx[K.SCRIPT_PATH],
            **single_kwargs,
        )
        workdir = result.artifacts["workdir"]
        print(f"\n SP work directory for DeePMD training: {workdir}")
        contents = sorted(p.name for p in workdir.iterdir())
        print(f"  Contents: {', '.join(contents)}")


class SpGenBatchCmd(MenuCommand):
    """Batch-generate SP work directories for DeePMD training."""

    params = (
        StrParam(K.XYZ, "XYZ trajectory file (e.g. md-pos-1.xyz)", required=True),
        DisplayAction(lambda ctx: _print_trajectory_info(ctx[K.XYZ])),
        cell_abc,
        *_frame_selection_params(),
        dp_sp_inp_template,
        StrParam(K.OUTDIR, "Output directory", default="."),
        cp2k_script,
    )

    def execute(self, ctx: dict) -> None:
        run_sp_batch = lazy_import(
            "md_analysis.workflows.scripts", "run_sp_batch",
        )
        kwargs = _batch_frame_kwargs_from_ctx(ctx)
        result = run_sp_batch(
            xyz_path=ctx[K.XYZ],
            cell_abc=ctx[K.CELL_ABC],
            output_dir=ctx[K.OUTDIR],
            inp_template_path=ctx[K.INP_TEMPLATE],
            script_path=ctx[K.SCRIPT_PATH],
            verbose=True,
            **kwargs,
        )
        workdirs = [
            v for k, v in result.artifacts.items() if k.startswith("workdir_")
        ]
        dirs = workdirs  # backward-compat name for the print loop below
        print(f"\n Created {len(dirs)} SP work directories for DeePMD training:")
        for d in dirs:
            print(f"  {d}")
