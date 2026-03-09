"""Scripts / Tools sub-menu."""

from __future__ import annotations

from pathlib import Path

from ._prompt import _handle_cmd_error, _prompt_bool, _prompt_cell_abc, _prompt_int, _prompt_str, _prompt_str_required

_MENU = """\

 ---------- Scripts / Tools ----------

 401) Generate Bader Work Directory (single frame)
 402) Batch Generate Bader Work Directories

   0) Back / Exit

"""


def scripts_menu() -> int:
    """Display the scripts sub-menu and dispatch."""
    print(_MENU)
    choice = input(" Input: ").strip()

    if choice == "0":
        print("\n Bye!")
        return 0

    if choice == "401":
        return _cmd_401()
    elif choice == "402":
        return _cmd_402()
    else:
        print(f"\n Invalid choice: {choice!r}")
        return 1


def _print_trajectory_info(xyz_path: str) -> None:
    """Peek at XYZ trajectory and print frame/step/time metadata."""
    from pathlib import Path

    p = Path(xyz_path)
    if not p.is_file():
        return

    # Count total frames via line count (fast, no parsing).
    with open(p) as fh:
        first_line = fh.readline().strip()
        try:
            natoms = int(first_line)
        except ValueError:
            return
        total_lines = sum(1 for _ in fh) + 1
    block_size = natoms + 2
    total_frames = total_lines // block_size

    # Read first 2 frames for step/time interval.
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
            print(f"  Frame interval: {step_interval} MD steps, {time_interval:.1f} fs/frame (dt = {dt:.1f} fs/step)")
    print()


@_handle_cmd_error
def _cmd_401() -> int:
    """Generate single-frame Bader work directory (menu 401)."""
    print()
    xyz_path = _prompt_str_required("XYZ trajectory file (e.g. md-pos-1.xyz)")
    _print_trajectory_info(xyz_path)
    abc = _prompt_cell_abc()

    frame = _prompt_int("Frame number (0-based)", default=0) or 0
    output_dir = _prompt_str("Output directory", default=".") or "."
    workdir_name = _prompt_str("Work directory name", default="bader") or "bader"

    from ..config import KEY_VASP_SCRIPT_PATH, get_config

    default_script = get_config(KEY_VASP_SCRIPT_PATH)
    script_path = _prompt_str(
        "Submission script path",
        default=default_script,
    )

    gen_potcar = _prompt_bool("Generate POTCAR via vaspkit?", default=True)

    # Read frame from trajectory
    from ase.io import iread

    print(f"\n Reading frame {frame} from {xyz_path} ...")
    atoms = None
    for i, a in enumerate(iread(xyz_path, index=":")):
        if i == frame:
            atoms = a
            break
    if atoms is None:
        print(f"  Error: frame {frame} not found in {xyz_path}")
        return 1

    atoms.set_cell(abc)
    atoms.set_pbc(True)

    from ..scripts import generate_bader_workdir

    workdir = generate_bader_workdir(
        atoms,
        output_dir,
        script_path=script_path,
        workdir_name=workdir_name,
        frame=frame,
        source=xyz_path,
        generate_potcar=gen_potcar,
    )

    print(f"\n Bader work directory created: {workdir}")
    contents = sorted(p.name for p in workdir.iterdir())
    print(f"  Contents: {', '.join(contents)}")
    return 0


@_handle_cmd_error
def _cmd_402() -> int:
    """Batch generate Bader work directories (menu 402)."""
    print()
    xyz_path = _prompt_str_required("XYZ trajectory file (e.g. md-pos-1.xyz)")
    _print_trajectory_info(xyz_path)
    abc = _prompt_cell_abc()

    frame_start = _prompt_int("Frame start (0-based)", default=0) or 0
    frame_end = _prompt_int("Frame end (exclusive, empty=all)", default=None)
    frame_step = _prompt_int("Frame step", default=1) or 1
    output_dir = _prompt_str("Output directory", default=".") or "."

    from ..config import KEY_VASP_SCRIPT_PATH, get_config

    default_script = get_config(KEY_VASP_SCRIPT_PATH)
    script_path = _prompt_str(
        "Submission script path",
        default=default_script,
    )

    gen_potcar = _prompt_bool("Generate POTCAR via vaspkit?", default=True)

    from ..scripts import batch_generate_bader_workdirs

    dirs = batch_generate_bader_workdirs(
        xyz_path,
        abc,
        output_dir,
        frame_start=frame_start,
        frame_end=frame_end,
        frame_step=frame_step,
        script_path=script_path,
        generate_potcar=gen_potcar,
        verbose=True,
    )

    print(f"\n Created {len(dirs)} Bader work directories:")
    for d in dirs:
        print(f"  {d}")
    return 0
