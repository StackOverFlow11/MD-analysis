"""Scripts / Tools sub-menu."""

from __future__ import annotations

from pathlib import Path

from ._prompt import _prompt_bool, _prompt_int, _prompt_str, _prompt_str_required

_MENU = """\

 ---------- Scripts / Tools ----------

 401) Generate Bader Work Directory

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
    else:
        print(f"\n Invalid choice: {choice!r}")
        return 1


def _cmd_401() -> int:
    print()
    xyz_path = _prompt_str_required("XYZ trajectory file (e.g. md-pos-1.xyz)")
    md_inp_path = _prompt_str_required("CP2K input file (e.g. md.inp)")
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
    for i, a in enumerate(iread(xyz_path, format="xyz")):
        if i == frame:
            atoms = a
            break
    if atoms is None:
        print(f"  Error: frame {frame} not found in {xyz_path}")
        return 1

    # Parse cell from md.inp
    from ..water.WaterAnalysis._common import _parse_abc_from_md_inp

    abc = _parse_abc_from_md_inp(md_inp_path)
    atoms.set_cell(abc)
    atoms.set_pbc(True)

    from ..scripts import generate_bader_workdir

    try:
        workdir = generate_bader_workdir(
            atoms,
            output_dir,
            script_path=script_path,
            workdir_name=workdir_name,
            frame=frame,
            source=xyz_path,
            generate_potcar=gen_potcar,
        )
    except Exception as exc:
        print(f"\n  Error: {exc}")
        return 1

    print(f"\n Bader work directory created: {workdir}")
    contents = sorted(p.name for p in workdir.iterdir())
    print(f"  Contents: {', '.join(contents)}")
    return 0
