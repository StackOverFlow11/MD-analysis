"""Tests for the logging infrastructure (NullHandler, CLI error logging)."""

import logging

import pytest

from md_analysis.cli._framework import MenuCommand
from md_analysis.exceptions import MDAnalysisError


class TestNullHandler:
    def test_md_analysis_logger_has_null_handler(self):
        import md_analysis  # noqa: F401
        lg = logging.getLogger("md_analysis")
        assert any(isinstance(h, logging.NullHandler) for h in lg.handlers)

    def test_null_handler_absorbs_messages(self, capfd):
        """With only NullHandler, log messages produce no output."""
        import md_analysis  # noqa: F401
        lg = logging.getLogger("md_analysis.test_silent")
        lg.info("should be silent")
        captured = capfd.readouterr()
        assert captured.out == ""
        assert captured.err == ""


class TestMenuCommandErrorLogging:
    def test_unexpected_exception_logs_error(self, caplog):

        class _RaiseTypeCmd(MenuCommand):
            def _collect_all_params(self):
                return {}

            def execute(self, ctx):
                raise TypeError("log me")

        cmd = _RaiseTypeCmd("t", "test")
        with caplog.at_level(logging.ERROR, logger="md_analysis.cli._framework"):
            cmd.run()

        assert any("log me" in r.message for r in caplog.records)
        assert any(r.levelno == logging.ERROR for r in caplog.records)
        assert any(r.exc_info is not None and r.exc_info[0] is TypeError for r in caplog.records)
