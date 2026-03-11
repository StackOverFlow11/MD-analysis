"""Settings command classes (901-907)."""

from __future__ import annotations

from pathlib import Path

from ._framework import MenuCommand, lazy_import
from ._prompt import prompt_float, prompt_str


class SetVaspScriptCmd(MenuCommand):
    def execute(self, ctx: dict) -> None:
        from ..config import KEY_VASP_SCRIPT_PATH, get_config, set_config

        current = get_config(KEY_VASP_SCRIPT_PATH)
        if current:
            print(f"  Current: {current}")

        raw = prompt_str("VASP submission script path", default=current)
        if raw is None:
            print("  No path provided, skipping.")
            return

        p = Path(raw).expanduser().resolve()
        if not p.is_file():
            print(f"  Warning: file does not exist: {p}")
            return

        set_config(KEY_VASP_SCRIPT_PATH, str(p))
        print(f"  Saved: {KEY_VASP_SCRIPT_PATH} = {p}")


class ShowConfigCmd(MenuCommand):
    def execute(self, ctx: dict) -> None:
        from ..config import CONFIGURABLE_DEFAULTS, get_config, load_config

        cfg = load_config()
        if not cfg:
            print("  (no user configuration set)")
        else:
            for key, val in cfg.items():
                print(f"  {key} = {val}")

        print("\n  --- Analysis Defaults ---")
        for key, entry in CONFIGURABLE_DEFAULTS.items():
            user_val = get_config(key)
            if user_val is not None:
                print(f"  {entry['label']}: {user_val}  (custom)")
            else:
                print(f"  {entry['label']}: {entry['default']}  (default)")


class SetAnalysisDefaultCmd(MenuCommand):
    """Generic command for setting one configurable analysis default."""

    def __init__(self, code: str, label: str, *, config_key: str):
        super().__init__(code, label)
        self.config_key = config_key

    def execute(self, ctx: dict) -> None:
        from ..config import CONFIGURABLE_DEFAULTS, get_config, set_config

        entry = CONFIGURABLE_DEFAULTS[self.config_key]
        current = get_config(self.config_key, entry["default"])
        value = prompt_float(entry["label"], default=current)
        if value <= 0:
            print("  Value must be > 0, skipping.")
            return
        set_config(self.config_key, value)
        print(f"  Saved: {self.config_key} = {value}")


class ResetDefaultsCmd(MenuCommand):
    def execute(self, ctx: dict) -> None:
        from ..config import CONFIGURABLE_DEFAULTS, delete_config

        for key in CONFIGURABLE_DEFAULTS:
            delete_config(key)
        print("\n  All analysis defaults reset to hardcoded values.")
