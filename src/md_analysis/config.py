"""Persistent user configuration stored at ~/.config/md_analysis/config.json."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

CONFIG_DIR = Path.home() / ".config" / "md_analysis"
CONFIG_FILE = CONFIG_DIR / "config.json"

KEY_VASP_SCRIPT_PATH = "vasp_script_path"
KEY_CP2K_SCRIPT_PATH = "cp2k_script_path"

# Keys for configurable analysis defaults (persisted in user config)
KEY_LAYER_TOL_A = "layer_tol_A"
KEY_Z_BIN_WIDTH_A = "z_bin_width_A"
KEY_THETA_BIN_DEG = "theta_bin_deg"
KEY_WATER_OH_CUTOFF_A = "water_oh_cutoff_A"

# Registry of configurable defaults: maps config key → metadata.
# The "default" values mirror the hardcoded constants in utils/config.py.
CONFIGURABLE_DEFAULTS: dict[str, dict] = {
    KEY_LAYER_TOL_A: {"default": 0.6, "label": "Layer clustering tolerance (A)"},
    KEY_Z_BIN_WIDTH_A: {"default": 0.1, "label": "Z-axis bin width (A)"},
    KEY_THETA_BIN_DEG: {"default": 5.0, "label": "Theta bin width (deg)"},
    KEY_WATER_OH_CUTOFF_A: {"default": 1.25, "label": "Water O-H cutoff (A)"},
}


from .exceptions import MDAnalysisError


class ConfigError(MDAnalysisError):
    """Raised on configuration read/write failures."""


def load_config(config_path: Path | None = None) -> dict[str, Any]:
    """Load configuration from disk. Returns ``{}`` if file does not exist."""
    path = config_path or CONFIG_FILE
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        raise ConfigError(f"Failed to load config from {path}: {exc}") from exc


def save_config(config: dict[str, Any], config_path: Path | None = None) -> Path:
    """Write *config* dict to disk, creating parent directories as needed."""
    path = config_path or CONFIG_FILE
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(config, indent=2, ensure_ascii=False) + "\n",
                        encoding="utf-8")
    except OSError as exc:
        raise ConfigError(f"Failed to save config to {path}: {exc}") from exc
    return path


def get_config(key: str, default: Any = None, *, config_path: Path | None = None) -> Any:
    """Read a single key from the persisted config."""
    return load_config(config_path).get(key, default)


def set_config(key: str, value: Any, *, config_path: Path | None = None) -> None:
    """Atomically set a single key (load -> merge -> save)."""
    cfg = load_config(config_path)
    cfg[key] = value
    save_config(cfg, config_path)


def delete_config(key: str, *, config_path: Path | None = None) -> None:
    """Remove *key* from the persisted config. No-op if absent."""
    cfg = load_config(config_path)
    if key in cfg:
        del cfg[key]
        save_config(cfg, config_path)
