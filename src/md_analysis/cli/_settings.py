"""Settings sub-menu for persistent user configuration."""

from __future__ import annotations

from pathlib import Path

from ._prompt import _prompt_str

_MENU = """\

 ---------- Settings ----------

 901) Set VASP Submission Script Path
 902) Show Current Configuration

   0) Back / Exit

"""


def settings_menu() -> int:
    """Display the settings sub-menu and dispatch."""
    print(_MENU)
    choice = input(" Input: ").strip()

    if choice == "0":
        print("\n Bye!")
        return 0

    if choice == "901":
        return _cmd_901()
    elif choice == "902":
        return _cmd_902()
    else:
        print(f"\n Invalid choice: {choice!r}")
        return 1


def _cmd_901() -> int:
    from ..config import KEY_VASP_SCRIPT_PATH, get_config, set_config

    current = get_config(KEY_VASP_SCRIPT_PATH)
    print()
    if current:
        print(f"  Current: {current}")

    raw = _prompt_str("VASP submission script path", default=current)
    if raw is None:
        print("  No path provided, skipping.")
        return 0

    p = Path(raw).expanduser().resolve()
    if not p.is_file():
        print(f"  Warning: file does not exist: {p}")
        return 1

    set_config(KEY_VASP_SCRIPT_PATH, str(p))
    print(f"  Saved: {KEY_VASP_SCRIPT_PATH} = {p}")
    return 0


def _cmd_902() -> int:
    from ..config import load_config

    cfg = load_config()
    print()
    if not cfg:
        print("  (no configuration set)")
    else:
        for key, val in cfg.items():
            print(f"  {key} = {val}")
    return 0
