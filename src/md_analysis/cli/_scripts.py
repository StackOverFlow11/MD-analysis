"""Scripts / Tools command classes (401-402)."""

from __future__ import annotations

from pathlib import Path

from ._framework import MenuCommand, lazy_import
from ._params import K, cell_abc
from ._prompt import prompt_bool, prompt_int, prompt_str, prompt_str_required


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


class BaderSingleCmd(MenuCommand):
    output_subdir = ""

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
    output_subdir = ""

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
