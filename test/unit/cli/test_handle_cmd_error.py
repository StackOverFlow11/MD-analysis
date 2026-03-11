"""Tests for MenuCommand.run() error handling (replaces _handle_cmd_error)."""

import logging

import pytest

from md_analysis.cli._framework import MenuCommand
from md_analysis.exceptions import MDAnalysisError


def _make_cmd(exc_class=None, msg="test error"):
    """Create a concrete MenuCommand that raises *exc_class* in execute()."""

    class _Cmd(MenuCommand):
        def _collect_all_params(self):
            return {}

        def execute(self, ctx):
            if exc_class is not None:
                raise exc_class(msg)

    return _Cmd("t", "test")


class TestMenuCommandErrorHandling:
    def test_normal_execution(self, capsys):
        cmd = _make_cmd()
        cmd.run()
        # no error output
        assert "Error" not in capsys.readouterr().out

    def test_md_analysis_error(self, capsys):
        cmd = _make_cmd(MDAnalysisError, "bad cube")
        cmd.run()
        assert "Error: bad cube" in capsys.readouterr().out

    def test_file_not_found(self, capsys):
        cmd = _make_cmd(FileNotFoundError, "missing.xyz")
        cmd.run()
        assert "Error: missing.xyz" in capsys.readouterr().out

    def test_value_error(self, capsys):
        cmd = _make_cmd(ValueError, "invalid axis")
        cmd.run()
        assert "Error: invalid axis" in capsys.readouterr().out

    def test_runtime_error(self, capsys):
        cmd = _make_cmd(RuntimeError, "something broke")
        cmd.run()
        assert "Error: something broke" in capsys.readouterr().out

    def test_unexpected_error(self, capsys):
        cmd = _make_cmd(TypeError, "bad type")
        cmd.run()
        out = capsys.readouterr().out
        assert "Unexpected error (TypeError):" in out
        assert "bad type" in out

    def test_unexpected_error_logged_with_traceback(self, caplog):
        cmd = _make_cmd(TypeError, "logged error")
        with caplog.at_level(logging.ERROR, logger="md_analysis.cli._framework"):
            cmd.run()
        assert any(
            r.levelno == logging.ERROR and "logged error" in r.getMessage()
            for r in caplog.records
        )
        assert any(
            r.exc_info is not None and r.exc_info[0] is TypeError
            for r in caplog.records
        )
