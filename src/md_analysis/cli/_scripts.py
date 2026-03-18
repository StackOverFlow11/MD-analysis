"""Scripts / Tools command classes (Bader 411-412, TI 421-422)."""

from __future__ import annotations

from pathlib import Path

from ._framework import MenuCommand, lazy_import
from ._params import K, cell_abc
from ._prompt import (
    prompt_bool,
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


def _resolve_script_path() -> str | None:
    """Prompt for VASP submission script path with config default."""
    from ..config import KEY_VASP_SCRIPT_PATH, get_config

    default_script = get_config(KEY_VASP_SCRIPT_PATH)
    return prompt_str("Submission script path", default=default_script)


def _resolve_cp2k_script_path() -> str | None:
    """Prompt for CP2K submission script path with config default."""
    from ..config import KEY_CP2K_SCRIPT_PATH, get_config

    default_script = get_config(KEY_CP2K_SCRIPT_PATH)
    return prompt_str("Submission script path", default=default_script)


class BaderSingleCmd(MenuCommand):

    def _collect_all_params(self) -> dict:
        """Custom flow: show trajectory info between prompts."""
        print()
        ctx: dict = {}
        ctx[K.XYZ] = prompt_str_required("XYZ trajectory file (e.g. md-pos-1.xyz)")
        _print_trajectory_info(ctx[K.XYZ])
        cell_abc.collect(ctx)
        ctx[K.FRAME] = prompt_int("Frame number (0-based)", default=0) or 0
        ctx[K.OUTDIR] = prompt_str("Output directory", default=".") or "."
        ctx[K.WORKDIR_NAME] = prompt_str("Work directory name", default="bader") or "bader"
        ctx[K.SCRIPT_PATH] = _resolve_script_path()
        ctx[K.GEN_POTCAR] = prompt_bool("Generate POTCAR via vaspkit?", default=True)
        return ctx

    def execute(self, ctx: dict) -> None:
        iread = lazy_import("ase.io", "iread")
        generate = lazy_import("md_analysis.scripts", "generate_bader_workdir")

        print(f"\n Reading frame {ctx[K.FRAME]} from {ctx[K.XYZ]} ...")
        atoms = None
        for i, a in enumerate(iread(ctx[K.XYZ], index=":")):
            if i == ctx[K.FRAME]:
                atoms = a
                break
        if atoms is None:
            print(f"  Error: frame {ctx[K.FRAME]} not found in {ctx[K.XYZ]}")
            return

        atoms.set_cell(ctx[K.CELL_ABC])
        atoms.set_pbc(True)

        workdir = generate(
            atoms,
            ctx[K.OUTDIR],
            script_path=ctx[K.SCRIPT_PATH],
            workdir_name=ctx[K.WORKDIR_NAME],
            frame=ctx[K.FRAME],
            source=ctx[K.XYZ],
            generate_potcar=ctx[K.GEN_POTCAR],
        )

        print(f"\n Bader work directory created: {workdir}")
        contents = sorted(p.name for p in workdir.iterdir())
        print(f"  Contents: {', '.join(contents)}")


class BaderBatchCmd(MenuCommand):

    def _collect_all_params(self) -> dict:
        """Custom flow: show trajectory info between prompts."""
        print()
        ctx: dict = {}
        ctx[K.XYZ] = prompt_str_required("XYZ trajectory file (e.g. md-pos-1.xyz)")
        _print_trajectory_info(ctx[K.XYZ])
        cell_abc.collect(ctx)
        ctx[K.FRAME_START] = prompt_int("Frame start (0-based)", default=0) or 0
        ctx[K.FRAME_END] = prompt_int("Frame end (exclusive, empty=all)", default=None)
        ctx[K.FRAME_STEP] = prompt_int("Frame step", default=1) or 1
        ctx[K.OUTDIR] = prompt_str("Output directory", default=".") or "."
        ctx[K.SCRIPT_PATH] = _resolve_script_path()
        ctx[K.GEN_POTCAR] = prompt_bool("Generate POTCAR via vaspkit?", default=True)
        return ctx

    def execute(self, ctx: dict) -> None:
        batch = lazy_import("md_analysis.scripts", "batch_generate_bader_workdirs")
        dirs = batch(
            ctx[K.XYZ],
            ctx[K.CELL_ABC],
            ctx[K.OUTDIR],
            frame_start=ctx[K.FRAME_START],
            frame_end=ctx[K.FRAME_END],
            frame_step=ctx[K.FRAME_STEP],
            script_path=ctx[K.SCRIPT_PATH],
            generate_potcar=ctx[K.GEN_POTCAR],
            verbose=True,
        )
        print(f"\n Created {len(dirs)} Bader work directories:")
        for d in dirs:
            print(f"  {d}")


# ---------------------------------------------------------------------------
# TI commands (42x)
# ---------------------------------------------------------------------------

def _print_sg_cv_info(restart_path: str, xyz_path: str) -> None:
    """Display SG trajectory CV range and frame info for TI target selection."""
    parse_colvar_restart = lazy_import(
        "md_analysis.utils.RestartParser.ColvarParser", "parse_colvar_restart",
    )
    try:
        restart = parse_colvar_restart(restart_path)
    except Exception as exc:
        print(f"  (Could not parse restart: {exc})")
        return

    cv = restart.colvars.primary
    au_time_to_fs = lazy_import("md_analysis.utils.config", "AU_TIME_TO_FS")
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
        generate = lazy_import("md_analysis.scripts", "generate_ti_workdir")
        workdir = generate(
            ctx[K.INP_PATH],
            ctx[K.XYZ],
            ctx[K.RESTART_PATH],
            ctx[K.TARGET_AU],
            ctx[K.OUTDIR],
            steps=ctx[K.STEPS],
            colvar_id=ctx[K.COLVAR_ID],
            workdir_name=ctx[K.WORKDIR_NAME],
            script_path=ctx[K.SCRIPT_PATH],
        )
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
