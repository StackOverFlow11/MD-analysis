"""Tests for settings menu 903-907 and _get_effective_default."""

from __future__ import annotations

from unittest.mock import patch

import pytest

from md_analysis.config import (
    CONFIGURABLE_DEFAULTS,
    KEY_LAYER_TOL_A,
    KEY_Z_BIN_WIDTH_A,
    delete_config,
    get_config,
    load_config,
    set_config,
)


# ---------- _get_effective_default ----------


class TestGetEffectiveDefault:
    """Tests for _get_effective_default helper."""

    def test_returns_hardcoded_when_no_user_config(self, tmp_path):
        cfg_path = tmp_path / "config.json"
        with patch("md_analysis.config.CONFIG_FILE", cfg_path):
            from md_analysis.cli._prompt import _get_effective_default
            val = _get_effective_default(KEY_LAYER_TOL_A)
        assert val == CONFIGURABLE_DEFAULTS[KEY_LAYER_TOL_A]["default"]

    def test_returns_user_value_when_set(self, tmp_path):
        cfg_path = tmp_path / "config.json"
        set_config(KEY_LAYER_TOL_A, 0.8, config_path=cfg_path)
        with patch("md_analysis.config.CONFIG_FILE", cfg_path):
            from md_analysis.cli._prompt import _get_effective_default
            val = _get_effective_default(KEY_LAYER_TOL_A)
        assert val == 0.8

    def test_raises_on_unknown_key(self, tmp_path):
        cfg_path = tmp_path / "config.json"
        with patch("md_analysis.config.CONFIG_FILE", cfg_path):
            from md_analysis.cli._prompt import _get_effective_default
            with pytest.raises(KeyError):
                _get_effective_default("nonexistent_key")


# ---------- _cmd_903 to _cmd_907 ----------


class TestSettingsHandlers:
    """Tests for individual settings menu handlers."""

    def test_cmd_903_sets_layer_tol(self, tmp_path):
        cfg_path = tmp_path / "config.json"
        with (
            patch("md_analysis.config.CONFIG_FILE", cfg_path),
            patch("builtins.input", return_value="0.8"),
        ):
            from md_analysis.cli._settings import _cmd_903
            result = _cmd_903()
        assert result == 0
        assert get_config(KEY_LAYER_TOL_A, config_path=cfg_path) == 0.8

    def test_cmd_904_sets_z_bin_width(self, tmp_path):
        cfg_path = tmp_path / "config.json"
        with (
            patch("md_analysis.config.CONFIG_FILE", cfg_path),
            patch("builtins.input", return_value="0.2"),
        ):
            from md_analysis.cli._settings import _cmd_904
            result = _cmd_904()
        assert result == 0
        assert get_config(KEY_Z_BIN_WIDTH_A, config_path=cfg_path) == 0.2

    def test_cmd_905_sets_theta_bin(self, tmp_path):
        from md_analysis.config import KEY_THETA_BIN_DEG
        cfg_path = tmp_path / "config.json"
        with (
            patch("md_analysis.config.CONFIG_FILE", cfg_path),
            patch("builtins.input", return_value="10.0"),
        ):
            from md_analysis.cli._settings import _cmd_905
            result = _cmd_905()
        assert result == 0
        assert get_config(KEY_THETA_BIN_DEG, config_path=cfg_path) == 10.0

    def test_cmd_906_sets_oh_cutoff(self, tmp_path):
        from md_analysis.config import KEY_WATER_OH_CUTOFF_A
        cfg_path = tmp_path / "config.json"
        with (
            patch("md_analysis.config.CONFIG_FILE", cfg_path),
            patch("builtins.input", return_value="1.3"),
        ):
            from md_analysis.cli._settings import _cmd_906
            result = _cmd_906()
        assert result == 0
        assert get_config(KEY_WATER_OH_CUTOFF_A, config_path=cfg_path) == 1.3

    def test_cmd_rejects_zero(self, tmp_path):
        cfg_path = tmp_path / "config.json"
        with (
            patch("md_analysis.config.CONFIG_FILE", cfg_path),
            patch("builtins.input", return_value="0"),
        ):
            from md_analysis.cli._settings import _cmd_903
            result = _cmd_903()
        assert result == 1
        assert not cfg_path.exists()

    def test_cmd_rejects_negative(self, tmp_path):
        cfg_path = tmp_path / "config.json"
        with (
            patch("md_analysis.config.CONFIG_FILE", cfg_path),
            patch("builtins.input", return_value="-1"),
        ):
            from md_analysis.cli._settings import _cmd_904
            result = _cmd_904()
        assert result == 1

    def test_cmd_907_resets_all(self, tmp_path):
        cfg_path = tmp_path / "config.json"
        # Set some custom values first
        set_config(KEY_LAYER_TOL_A, 0.9, config_path=cfg_path)
        set_config(KEY_Z_BIN_WIDTH_A, 0.5, config_path=cfg_path)
        set_config("vasp_script_path", "/some/path", config_path=cfg_path)

        with patch("md_analysis.config.CONFIG_FILE", cfg_path):
            from md_analysis.cli._settings import _cmd_907
            result = _cmd_907()
        assert result == 0

        data = load_config(cfg_path)
        # Analysis defaults removed, but vasp_script_path preserved
        assert KEY_LAYER_TOL_A not in data
        assert KEY_Z_BIN_WIDTH_A not in data
        assert data["vasp_script_path"] == "/some/path"

    def test_cmd_902_shows_defaults(self, tmp_path, capsys):
        cfg_path = tmp_path / "config.json"
        set_config(KEY_LAYER_TOL_A, 0.9, config_path=cfg_path)

        with patch("md_analysis.config.CONFIG_FILE", cfg_path):
            from md_analysis.cli._settings import _cmd_902
            result = _cmd_902()
        assert result == 0

        output = capsys.readouterr().out
        assert "(custom)" in output
        assert "(default)" in output
        assert "0.9" in output
