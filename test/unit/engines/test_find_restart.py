"""Regression tests for :meth:`CP2KParser._find_restart`.

CP2K writes ``<PROJECT>-<RUN>.restart`` (rolling, per-step — the latest
state) and ``<PROJECT>-<RUN>_<STEP>.restart`` (RESTART_HISTORY periodic
snapshots — always older than the rolling file). The picker must prefer
the rolling file and never select a snapshot when a rolling file
exists.
"""

from __future__ import annotations

import pytest

from md_analysis.engines import CP2KParser


def _touch(directory, name: str) -> None:
    (directory / name).write_text("", encoding="utf-8")


class TestFindRestart:
    def test_prefers_bare_over_history_snapshots(self, tmp_path):
        """Bug regression: lex-sort previously picked the snapshot."""
        _touch(tmp_path, "cMD-1.restart")
        _touch(tmp_path, "cMD-1_500.restart")
        _touch(tmp_path, "cMD-1_1500.restart")
        assert CP2KParser._find_restart(tmp_path).name == "cMD-1.restart"

    def test_ignores_bak_backups(self, tmp_path):
        _touch(tmp_path, "cMD-1.restart")
        _touch(tmp_path, "cMD-1.restart.bak-1")
        _touch(tmp_path, "cMD-1.restart.bak-2")
        assert CP2KParser._find_restart(tmp_path).name == "cMD-1.restart"

    def test_only_snapshots_present_raises(self, tmp_path):
        """User workflow keeps only a bare rolling file via mv on restart;
        a snapshot-only directory is an exceptional state and must error
        explicitly rather than silently return a stale snapshot."""
        _touch(tmp_path, "cMD-1_500.restart")
        _touch(tmp_path, "cMD-1_1500.restart")
        with pytest.raises(FileNotFoundError):
            CP2KParser._find_restart(tmp_path)

    def test_empty_directory_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            CP2KParser._find_restart(tmp_path)

    def test_multi_run_picks_highest_run_index(self, tmp_path):
        """If multiple bare-RUN files coexist, the highest RUN wins."""
        _touch(tmp_path, "cMD-1.restart")
        _touch(tmp_path, "cMD-2.restart")
        assert CP2KParser._find_restart(tmp_path).name == "cMD-2.restart"
