"""Settings command classes (901-907)."""

from __future__ import annotations

from pathlib import Path

from ._framework import MenuCommand, lazy_import
from ._prompt import prompt_choice, prompt_float, prompt_str


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


class SetCp2kScriptCmd(MenuCommand):
    def execute(self, ctx: dict) -> None:
        from ..config import KEY_CP2K_SCRIPT_PATH, get_config, set_config

        current = get_config(KEY_CP2K_SCRIPT_PATH)
        if current:
            print(f"  Current: {current}")

        raw = prompt_str("CP2K submission script path", default=current)
        if raw is None:
            print("  No path provided, skipping.")
            return

        p = Path(raw).expanduser().resolve()
        if not p.is_file():
            print(f"  Warning: file does not exist: {p}")
            return

        set_config(KEY_CP2K_SCRIPT_PATH, str(p))
        print(f"  Saved: {KEY_CP2K_SCRIPT_PATH} = {p}")


class SetSpInpTemplateCmd(MenuCommand):
    def execute(self, ctx: dict) -> None:
        from ..config import KEY_SP_INP_TEMPLATE_PATH, get_config, set_config

        current = get_config(KEY_SP_INP_TEMPLATE_PATH)
        if current:
            print(f"  Current: {current}")

        raw = prompt_str("SP inp template path (e.g. sp.inp)", default=current)
        if raw is None:
            print("  No path provided, skipping.")
            return

        p = Path(raw).expanduser().resolve()
        if not p.is_file():
            print(f"  Warning: file does not exist: {p}")
            return

        set_config(KEY_SP_INP_TEMPLATE_PATH, str(p))
        print(f"  Saved: {KEY_SP_INP_TEMPLATE_PATH} = {p}")


class ShowConfigCmd(MenuCommand):
    def execute(self, ctx: dict) -> None:
        from ..config import (
            CONFIGURABLE_DEFAULTS,
            KEY_POTENTIAL_PH,
            KEY_POTENTIAL_PHI_PZC,
            KEY_POTENTIAL_REFERENCE,
            KEY_POTENTIAL_TEMPERATURE_K,
            get_config,
            load_config,
        )

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

        print("\n  --- Potential Output ---")
        ref = get_config(KEY_POTENTIAL_REFERENCE, "SHE")
        is_custom = get_config(KEY_POTENTIAL_REFERENCE) is not None
        print(f"  Reference: {ref}  {'(custom)' if is_custom else '(default)'}")
        if ref == "RHE":
            print(f"  pH: {get_config(KEY_POTENTIAL_PH, 0.0)}")
            print(f"  Temperature: {get_config(KEY_POTENTIAL_TEMPERATURE_K, 298.15)} K")
        elif ref == "PZC":
            print(f"  φ_PZC: {get_config(KEY_POTENTIAL_PHI_PZC, 0.0)} V vs SHE")


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


class SetPotentialReferenceCmd(MenuCommand):
    """Configure default potential output reference (SHE/RHE/PZC)."""

    def execute(self, ctx: dict) -> None:
        from ..config import (
            KEY_POTENTIAL_PH,
            KEY_POTENTIAL_PHI_PZC,
            KEY_POTENTIAL_REFERENCE,
            KEY_POTENTIAL_TEMPERATURE_K,
            get_config,
            set_config,
        )

        current_ref = get_config(KEY_POTENTIAL_REFERENCE, "SHE")
        print(f"  Current reference: {current_ref}")

        ref = prompt_choice("Potential reference", ["SHE", "RHE", "PZC"],
                            default=current_ref)
        set_config(KEY_POTENTIAL_REFERENCE, ref)

        if ref == "RHE":
            current_pH = get_config(KEY_POTENTIAL_PH, 0.0)
            current_T = get_config(KEY_POTENTIAL_TEMPERATURE_K, 298.15)
            pH = prompt_float("pH", default=current_pH)
            T = prompt_float("Temperature (K)", default=current_T)
            set_config(KEY_POTENTIAL_PH, pH)
            set_config(KEY_POTENTIAL_TEMPERATURE_K, T)
            print(f"\n  Saved: reference={ref}, pH={pH}, T={T} K")
        elif ref == "PZC":
            current_pzc = get_config(KEY_POTENTIAL_PHI_PZC, 0.0)
            pzc = prompt_float("Potential of zero charge (V vs SHE)",
                               default=current_pzc)
            set_config(KEY_POTENTIAL_PHI_PZC, pzc)
            print(f"\n  Saved: reference={ref}, φ_PZC={pzc} V vs SHE")
        else:
            print(f"\n  Saved: reference={ref}")


class ResetDefaultsCmd(MenuCommand):
    def execute(self, ctx: dict) -> None:
        from ..config import (
            CONFIGURABLE_DEFAULTS,
            KEY_POTENTIAL_PH,
            KEY_POTENTIAL_PHI_PZC,
            KEY_POTENTIAL_REFERENCE,
            KEY_POTENTIAL_TEMPERATURE_K,
            delete_config,
        )

        for key in CONFIGURABLE_DEFAULTS:
            delete_config(key)
        for key in (KEY_POTENTIAL_REFERENCE, KEY_POTENTIAL_PH,
                    KEY_POTENTIAL_TEMPERATURE_K, KEY_POTENTIAL_PHI_PZC):
            delete_config(key)
        print("\n  All analysis defaults reset to hardcoded values.")
