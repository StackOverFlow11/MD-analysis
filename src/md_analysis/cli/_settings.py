"""Settings sub-menu for persistent user configuration."""

from __future__ import annotations

from pathlib import Path

from ._prompt import _handle_cmd_error, _prompt_float, _prompt_str

_MENU = """\

 ---------- Settings ----------

 901) Set VASP Submission Script Path
 902) Show Current Configuration

 --- Analysis Defaults ---
 903) Layer Clustering Tolerance (A)
 904) Z-axis Bin Width (A)
 905) Theta Bin Width (deg)
 906) Water O-H Cutoff (A)
 907) Reset All Defaults

   0) Back / Exit

"""

_DISPATCH: dict[str, object] = {}  # populated after function definitions


def settings_menu() -> int:
    """Display the settings sub-menu and dispatch."""
    print(_MENU)
    choice = input(" Input: ").strip()

    if choice == "0":
        print("\n Bye!")
        return 0

    handler = _DISPATCH.get(choice)
    if handler is None:
        print(f"\n Invalid choice: {choice!r}")
        return 1

    return handler()


@_handle_cmd_error
def _cmd_901() -> int:
    """Set VASP submission script path (menu 901)."""
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


@_handle_cmd_error
def _cmd_902() -> int:
    """Show current configuration (menu 902)."""
    from ..config import CONFIGURABLE_DEFAULTS, get_config, load_config

    cfg = load_config()
    print()
    if not cfg:
        print("  (no user configuration set)")
    else:
        for key, val in cfg.items():
            print(f"  {key} = {val}")

    # Show effective analysis defaults
    print("\n  --- Analysis Defaults ---")
    for key, entry in CONFIGURABLE_DEFAULTS.items():
        user_val = get_config(key)
        if user_val is not None:
            print(f"  {entry['label']}: {user_val}  (custom)")
        else:
            print(f"  {entry['label']}: {entry['default']}  (default)")
    return 0


def _cmd_set_default(key: str) -> int:
    """Generic handler for setting a single configurable analysis default."""
    from ..config import CONFIGURABLE_DEFAULTS, get_config, set_config

    entry = CONFIGURABLE_DEFAULTS[key]
    current = get_config(key, entry["default"])
    print()
    value = _prompt_float(entry["label"], default=current)
    if value <= 0:
        print("  Value must be > 0, skipping.")
        return 1
    set_config(key, value)
    print(f"  Saved: {key} = {value}")
    return 0


@_handle_cmd_error
def _cmd_903() -> int:
    """Set layer clustering tolerance (menu 903)."""
    from ..config import KEY_LAYER_TOL_A
    return _cmd_set_default(KEY_LAYER_TOL_A)


@_handle_cmd_error
def _cmd_904() -> int:
    """Set z-axis bin width (menu 904)."""
    from ..config import KEY_Z_BIN_WIDTH_A
    return _cmd_set_default(KEY_Z_BIN_WIDTH_A)


@_handle_cmd_error
def _cmd_905() -> int:
    """Set theta bin width (menu 905)."""
    from ..config import KEY_THETA_BIN_DEG
    return _cmd_set_default(KEY_THETA_BIN_DEG)


@_handle_cmd_error
def _cmd_906() -> int:
    """Set water O-H cutoff (menu 906)."""
    from ..config import KEY_WATER_OH_CUTOFF_A
    return _cmd_set_default(KEY_WATER_OH_CUTOFF_A)


@_handle_cmd_error
def _cmd_907() -> int:
    """Reset all analysis defaults to hardcoded values (menu 907)."""
    from ..config import CONFIGURABLE_DEFAULTS, delete_config

    for key in CONFIGURABLE_DEFAULTS:
        delete_config(key)
    print("\n  All analysis defaults reset to hardcoded values.")
    return 0


_DISPATCH.update({
    "901": _cmd_901,
    "902": _cmd_902,
    "903": _cmd_903,
    "904": _cmd_904,
    "905": _cmd_905,
    "906": _cmd_906,
    "907": _cmd_907,
})
