"""Tests for md_analysis.config — persistent user configuration."""

from __future__ import annotations

import json

import pytest

from md_analysis.config import (
    CONFIGURABLE_DEFAULTS,
    ConfigError,
    KEY_LAYER_TOL_A,
    KEY_THETA_BIN_DEG,
    KEY_WATER_OH_CUTOFF_A,
    KEY_Z_BIN_WIDTH_A,
    delete_config,
    get_config,
    load_config,
    save_config,
    set_config,
)


def test_load_missing_file_returns_empty(tmp_path):
    """load_config returns {} when the config file does not exist."""
    assert load_config(tmp_path / "nonexistent.json") == {}


def test_save_and_load_roundtrip(tmp_path):
    """save_config followed by load_config preserves data."""
    cfg_path = tmp_path / "sub" / "config.json"
    data = {"key1": "value1", "key2": 42, "nested": {"a": True}}
    save_config(data, cfg_path)
    loaded = load_config(cfg_path)
    assert loaded == data


def test_get_config_default(tmp_path):
    """get_config returns default when key is absent."""
    cfg_path = tmp_path / "config.json"
    assert get_config("missing_key", "fallback", config_path=cfg_path) == "fallback"
    assert get_config("missing_key", config_path=cfg_path) is None


def test_set_config_creates_file(tmp_path):
    """set_config creates the config file and parent directories."""
    cfg_path = tmp_path / "deep" / "dir" / "config.json"
    set_config("foo", "bar", config_path=cfg_path)
    assert cfg_path.exists()
    data = json.loads(cfg_path.read_text())
    assert data == {"foo": "bar"}


def test_set_config_merges(tmp_path):
    """set_config preserves existing keys when adding a new one."""
    cfg_path = tmp_path / "config.json"
    set_config("key1", "val1", config_path=cfg_path)
    set_config("key2", "val2", config_path=cfg_path)
    data = load_config(cfg_path)
    assert data == {"key1": "val1", "key2": "val2"}


def test_load_corrupt_file_raises(tmp_path):
    """load_config raises ConfigError on invalid JSON."""
    cfg_path = tmp_path / "config.json"
    cfg_path.write_text("not valid json {{{")
    with pytest.raises(ConfigError):
        load_config(cfg_path)


# ---------- delete_config ----------

def test_delete_config_removes_key(tmp_path):
    """delete_config removes an existing key."""
    cfg_path = tmp_path / "config.json"
    set_config("a", 1, config_path=cfg_path)
    set_config("b", 2, config_path=cfg_path)
    delete_config("a", config_path=cfg_path)
    assert load_config(cfg_path) == {"b": 2}


def test_delete_config_noop_missing_key(tmp_path):
    """delete_config is a no-op when the key does not exist."""
    cfg_path = tmp_path / "config.json"
    set_config("x", 10, config_path=cfg_path)
    delete_config("nonexistent", config_path=cfg_path)
    assert load_config(cfg_path) == {"x": 10}


def test_delete_config_noop_no_file(tmp_path):
    """delete_config is a no-op when the config file does not exist."""
    cfg_path = tmp_path / "config.json"
    delete_config("anything", config_path=cfg_path)
    assert not cfg_path.exists()


# ---------- CONFIGURABLE_DEFAULTS registry ----------

def test_configurable_defaults_keys():
    """CONFIGURABLE_DEFAULTS contains exactly the expected keys."""
    expected = {KEY_LAYER_TOL_A, KEY_Z_BIN_WIDTH_A, KEY_THETA_BIN_DEG, KEY_WATER_OH_CUTOFF_A}
    assert set(CONFIGURABLE_DEFAULTS.keys()) == expected


def test_configurable_defaults_have_required_fields():
    """Each entry has 'default' (float) and 'label' (str)."""
    for key, entry in CONFIGURABLE_DEFAULTS.items():
        assert isinstance(entry["default"], float), f"{key}: default is not float"
        assert isinstance(entry["label"], str), f"{key}: label is not str"
        assert entry["default"] > 0, f"{key}: default must be positive"
