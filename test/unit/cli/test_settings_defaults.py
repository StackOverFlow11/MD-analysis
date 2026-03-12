"""Tests for settings commands and ConfigDefaultParam."""

from __future__ import annotations

from unittest.mock import patch

import pytest

from md_analysis.cli._prompt import set_input_source
from md_analysis.config import (
    CONFIGURABLE_DEFAULTS,
    KEY_LAYER_TOL_A,
    KEY_THETA_BIN_DEG,
    KEY_WATER_OH_CUTOFF_A,
    KEY_Z_BIN_WIDTH_A,
    delete_config,
    get_config,
    load_config,
    set_config,
)


@pytest.fixture(autouse=True)
def _restore_input():
    """Reset input source after each test."""
    yield
    set_input_source(input)


def _scripted(*responses):
    """Return a callable that yields scripted responses in order."""
    it = iter(responses)
    return lambda _: next(it)


# ---------- ConfigDefaultParam._get_effective_default ----------


class TestConfigDefaultParam:
    """Tests for ConfigDefaultParam effective default logic."""

    def test_returns_hardcoded_when_no_user_config(self, tmp_path):
        cfg_path = tmp_path / "config.json"
        with patch("md_analysis.config.CONFIG_FILE", cfg_path):
            from md_analysis.cli._params import ConfigDefaultParam

            p = ConfigDefaultParam("k", "label", config_key=KEY_LAYER_TOL_A)
            val = p._get_effective_default()
        assert val == CONFIGURABLE_DEFAULTS[KEY_LAYER_TOL_A]["default"]

    def test_returns_user_value_when_set(self, tmp_path):
        cfg_path = tmp_path / "config.json"
        set_config(KEY_LAYER_TOL_A, 0.8, config_path=cfg_path)
        with patch("md_analysis.config.CONFIG_FILE", cfg_path):
            from md_analysis.cli._params import ConfigDefaultParam

            p = ConfigDefaultParam("k", "label", config_key=KEY_LAYER_TOL_A)
            val = p._get_effective_default()
        assert val == 0.8

    def test_raises_on_unknown_key(self, tmp_path):
        cfg_path = tmp_path / "config.json"
        with patch("md_analysis.config.CONFIG_FILE", cfg_path):
            from md_analysis.cli._params import ConfigDefaultParam

            p = ConfigDefaultParam("k", "label", config_key="nonexistent_key")
            with pytest.raises(KeyError):
                p._get_effective_default()


# ---------- Settings commands ----------


class TestSettingsCommands:
    """Tests for individual settings menu command classes."""

    def test_set_layer_tol(self, tmp_path):
        from md_analysis.cli._settings import SetAnalysisDefaultCmd

        cfg_path = tmp_path / "config.json"
        set_input_source(_scripted("0.8"))
        cmd = SetAnalysisDefaultCmd("903", "test", config_key=KEY_LAYER_TOL_A)
        with patch("md_analysis.config.CONFIG_FILE", cfg_path):
            cmd.execute({})
        assert get_config(KEY_LAYER_TOL_A, config_path=cfg_path) == 0.8

    def test_set_z_bin_width(self, tmp_path):
        from md_analysis.cli._settings import SetAnalysisDefaultCmd

        cfg_path = tmp_path / "config.json"
        set_input_source(_scripted("0.2"))
        cmd = SetAnalysisDefaultCmd("904", "test", config_key=KEY_Z_BIN_WIDTH_A)
        with patch("md_analysis.config.CONFIG_FILE", cfg_path):
            cmd.execute({})
        assert get_config(KEY_Z_BIN_WIDTH_A, config_path=cfg_path) == 0.2

    def test_set_theta_bin(self, tmp_path):
        from md_analysis.cli._settings import SetAnalysisDefaultCmd

        cfg_path = tmp_path / "config.json"
        set_input_source(_scripted("10.0"))
        cmd = SetAnalysisDefaultCmd("905", "test", config_key=KEY_THETA_BIN_DEG)
        with patch("md_analysis.config.CONFIG_FILE", cfg_path):
            cmd.execute({})
        assert get_config(KEY_THETA_BIN_DEG, config_path=cfg_path) == 10.0

    def test_set_oh_cutoff(self, tmp_path):
        from md_analysis.cli._settings import SetAnalysisDefaultCmd

        cfg_path = tmp_path / "config.json"
        set_input_source(_scripted("1.3"))
        cmd = SetAnalysisDefaultCmd("906", "test", config_key=KEY_WATER_OH_CUTOFF_A)
        with patch("md_analysis.config.CONFIG_FILE", cfg_path):
            cmd.execute({})
        assert get_config(KEY_WATER_OH_CUTOFF_A, config_path=cfg_path) == 1.3

    def test_rejects_zero(self, tmp_path, capsys):
        from md_analysis.cli._settings import SetAnalysisDefaultCmd

        cfg_path = tmp_path / "config.json"
        set_input_source(_scripted("0"))
        cmd = SetAnalysisDefaultCmd("903", "test", config_key=KEY_LAYER_TOL_A)
        with patch("md_analysis.config.CONFIG_FILE", cfg_path):
            cmd.execute({})
        assert not cfg_path.exists()
        assert "must be > 0" in capsys.readouterr().out

    def test_rejects_negative(self, tmp_path, capsys):
        from md_analysis.cli._settings import SetAnalysisDefaultCmd

        cfg_path = tmp_path / "config.json"
        set_input_source(_scripted("-1"))
        cmd = SetAnalysisDefaultCmd("904", "test", config_key=KEY_Z_BIN_WIDTH_A)
        with patch("md_analysis.config.CONFIG_FILE", cfg_path):
            cmd.execute({})
        assert "must be > 0" in capsys.readouterr().out

    def test_reset_all_defaults(self, tmp_path):
        from md_analysis.cli._settings import ResetDefaultsCmd

        cfg_path = tmp_path / "config.json"
        set_config(KEY_LAYER_TOL_A, 0.9, config_path=cfg_path)
        set_config(KEY_Z_BIN_WIDTH_A, 0.5, config_path=cfg_path)
        set_config("vasp_script_path", "/some/path", config_path=cfg_path)

        cmd = ResetDefaultsCmd("907", "test")
        with patch("md_analysis.config.CONFIG_FILE", cfg_path):
            cmd.execute({})

        data = load_config(cfg_path)
        assert KEY_LAYER_TOL_A not in data
        assert KEY_Z_BIN_WIDTH_A not in data
        assert data["vasp_script_path"] == "/some/path"

    def test_show_config(self, tmp_path, capsys):
        from md_analysis.cli._settings import ShowConfigCmd

        cfg_path = tmp_path / "config.json"
        set_config(KEY_LAYER_TOL_A, 0.9, config_path=cfg_path)

        cmd = ShowConfigCmd("902", "test")
        with patch("md_analysis.config.CONFIG_FILE", cfg_path):
            cmd.execute({})

        output = capsys.readouterr().out
        assert "(custom)" in output
        assert "(default)" in output
        assert "0.9" in output
