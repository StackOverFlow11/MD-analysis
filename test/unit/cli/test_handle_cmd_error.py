"""Tests for the _handle_cmd_error decorator."""

import logging

import pytest

from md_analysis.cli._prompt import _handle_cmd_error
from md_analysis.exceptions import MDAnalysisError


@_handle_cmd_error
def _ok():
    return 0


@_handle_cmd_error
def _raise_md(msg):
    raise MDAnalysisError(msg)


@_handle_cmd_error
def _raise_fnf(msg):
    raise FileNotFoundError(msg)


@_handle_cmd_error
def _raise_val(msg):
    raise ValueError(msg)


@_handle_cmd_error
def _raise_rt(msg):
    raise RuntimeError(msg)


@_handle_cmd_error
def _raise_type(msg):
    raise TypeError(msg)


class TestHandleCmdError:
    def test_normal_return(self):
        assert _ok() == 0

    def test_md_analysis_error(self, capsys):
        assert _raise_md("bad cube") == 1
        assert "Error: bad cube" in capsys.readouterr().out

    def test_file_not_found(self, capsys):
        assert _raise_fnf("missing.xyz") == 1
        assert "Error: missing.xyz" in capsys.readouterr().out

    def test_value_error(self, capsys):
        assert _raise_val("invalid axis") == 1
        assert "Error: invalid axis" in capsys.readouterr().out

    def test_runtime_error(self, capsys):
        assert _raise_rt("something broke") == 1
        assert "Error: something broke" in capsys.readouterr().out

    def test_unexpected_error(self, capsys):
        assert _raise_type("bad type") == 1
        out = capsys.readouterr().out
        assert "Unexpected error (TypeError):" in out
        assert "bad type" in out

    def test_unexpected_error_logged_with_traceback(self, caplog):
        with caplog.at_level(logging.ERROR, logger="md_analysis.cli._prompt"):
            _raise_type("logged error")
        assert any(
            r.levelno == logging.ERROR and "logged error" in r.message
            for r in caplog.records
        )
        assert any(
            r.exc_info is not None and r.exc_info[0] is TypeError
            for r in caplog.records
        )
