"""Persistent user configuration stored at ~/.config/md_analysis/config.json."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

CONFIG_DIR = Path.home() / ".config" / "md_analysis"
CONFIG_FILE = CONFIG_DIR / "config.json"

KEY_VASP_SCRIPT_PATH = "vasp_script_path"


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
