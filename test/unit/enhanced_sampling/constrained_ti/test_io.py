"""Unit tests for constrained_ti.io — discover_ti_points and parser-driven dispatch."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from md_analysis.engines import (
    CP2KParser,
    ParserInferenceError,
)
from md_analysis.enhanced_sampling.constrained_ti.io import (
    discover_ti_points,
    load_ti_series,
)


REPO_ROOT = Path(__file__).resolve().parents[4]
EXAMPLE_TI_ROOT = (
    REPO_ROOT / "data_example" / "ti" / "double_cv" / "1k"
)
EXAMPLE_POINT = EXAMPLE_TI_ROOT / "ti_target_0.031369"


_skip_no_example = pytest.mark.skipif(
    not EXAMPLE_TI_ROOT.is_dir(),
    reason=f"example TI data not at {EXAMPLE_TI_ROOT}",
)


def _seed_point(parent: Path, name: str) -> Path:
    """Copy a real CP2K TI-point directory under *parent*/<name>.

    Uses the example fixture so parse_colvar_restart succeeds.  Tests
    that need to check ξ ordering can rename the destination at will.
    """
    src = EXAMPLE_POINT
    dst = parent / name
    shutil.copytree(src, dst)
    return dst


# ---------------------------------------------------------------------------
# Default discovery (parser="auto", dir_filter=None)
# ---------------------------------------------------------------------------


@_skip_no_example
class TestDiscoverDefault:
    """parser='auto' + dir_filter=None — recognises any directory the
    parser can read, regardless of name."""

    def test_finds_real_ti_root(self):
        points = discover_ti_points(EXAMPLE_TI_ROOT)
        assert len(points) >= 2
        assert all(p.parser.name == "cp2k" for p in points)

    def test_default_ascending_order(self):
        points = discover_ti_points(EXAMPLE_TI_ROOT)
        xis = [p.xi for p in points]
        assert xis == sorted(xis)

    def test_reverse_descending_order(self):
        points = discover_ti_points(EXAMPLE_TI_ROOT, reverse=True)
        xis = [p.xi for p in points]
        assert xis == sorted(xis, reverse=True)

    def test_reverse_same_set_as_default(self):
        asc = discover_ti_points(EXAMPLE_TI_ROOT)
        desc = discover_ti_points(EXAMPLE_TI_ROOT, reverse=True)
        assert {p.xi for p in asc} == {p.xi for p in desc}
        assert [p.xi for p in desc] == list(reversed([p.xi for p in asc]))


# ---------------------------------------------------------------------------
# Engine independence: directory name must NOT determine ξ
# ---------------------------------------------------------------------------


@_skip_no_example
class TestDirectoryNamingFreedom:
    """ξ comes from the restart file's TARGET, not the directory name."""

    def test_arbitrary_name_still_works(self, tmp_path):
        """Renaming the directory to something unrelated (no 'ti_target_'
        or 'xi_' prefix) must not change ξ."""
        d = _seed_point(tmp_path, "pt_001")
        points = discover_ti_points(tmp_path)
        assert len(points) == 1
        # Real point's ξ is 0.031369 (from directory name in source) —
        # but we verify against the metadata, not the directory name.
        assert points[0].xi == pytest.approx(0.031369, abs=1e-5)
        assert points[0].directory == d

    def test_misleading_name_uses_metadata(self, tmp_path):
        """If a user names the dir misleadingly, ξ still comes from the
        actual restart content."""
        _seed_point(tmp_path, "ti_target_999.0")  # misleading name
        points = discover_ti_points(tmp_path)
        assert len(points) == 1
        assert points[0].xi == pytest.approx(0.031369, abs=1e-5)
        assert points[0].xi != 999.0


# ---------------------------------------------------------------------------
# dir_filter resolution
# ---------------------------------------------------------------------------


@_skip_no_example
class TestDirFilter:
    def test_glob_string_filter(self, tmp_path):
        _seed_point(tmp_path, "keep_one")
        _seed_point(tmp_path, "skip_two")
        points = discover_ti_points(tmp_path, dir_filter="keep_*")
        assert len(points) == 1
        assert points[0].directory.name == "keep_one"

    def test_callable_filter(self, tmp_path):
        _seed_point(tmp_path, "alpha")
        _seed_point(tmp_path, "beta")
        points = discover_ti_points(
            tmp_path, dir_filter=lambda d: d.name.startswith("a"),
        )
        assert len(points) == 1
        assert points[0].directory.name == "alpha"

    def test_none_means_content_based_auto(self, tmp_path):
        """dir_filter=None uses parser.is_constraint_directory; helper
        directories without the right files are silently ignored."""
        _seed_point(tmp_path, "real_point")
        # Empty helper dir — no restart, no log
        (tmp_path / "scratch").mkdir()
        points = discover_ti_points(tmp_path)
        assert len(points) == 1
        assert points[0].directory.name == "real_point"


# ---------------------------------------------------------------------------
# Strict mode
# ---------------------------------------------------------------------------


def _make_named_dir(parent: Path, name: str, *, with_restart=True, with_log=True):
    """Create a directory with optional CP2K-recognisable empty files.
    Used for strict-mode tests where actual content doesn't matter
    (only file presence is checked when paired with name-based filter)."""
    d = parent / name
    d.mkdir(parents=True)
    if with_restart:
        (d / "cMD-1.restart").write_text("placeholder\n")
    if with_log:
        (d / "cMD-1.LagrangeMultLog").write_text("placeholder\n")
    return d


@_skip_no_example
class TestStrictMode:
    """When users force a name-based filter, strict=True surfaces parse
    failures rather than silently skipping."""

    def test_lenient_skips_unparseable_dir(self, tmp_path):
        # A real point + an incomplete dir with the right name pattern
        _seed_point(tmp_path, "real_0.031369")  # parses fine
        _make_named_dir(tmp_path, "fake_0.5")    # placeholder content fails parse

        points = discover_ti_points(
            tmp_path, dir_filter="*", strict=False,
        )
        # Only real_0.031369 survives; fake_0.5 silently skipped
        assert len(points) == 1
        assert points[0].directory.name == "real_0.031369"

    def test_strict_raises_on_unparseable_dir(self, tmp_path):
        _seed_point(tmp_path, "real_0.031369")
        _make_named_dir(tmp_path, "fake_0.5")

        with pytest.raises(FileNotFoundError, match="fake_0.5"):
            discover_ti_points(tmp_path, dir_filter="*", strict=True)

    def test_strict_raises_on_missing_log(self, tmp_path):
        """Glob filter forces selection; missing log file → strict failure."""
        _make_named_dir(tmp_path, "ti_target_0.1", with_log=False)

        with pytest.raises(FileNotFoundError):
            discover_ti_points(
                tmp_path, dir_filter="ti_target_*", strict=True,
            )


# ---------------------------------------------------------------------------
# parser argument resolution
# ---------------------------------------------------------------------------


@_skip_no_example
class TestParserArg:
    def test_explicit_instance(self, tmp_path):
        """Explicit parser bypasses sniffing."""
        _seed_point(tmp_path, "p")
        points = discover_ti_points(tmp_path, parser=CP2KParser())
        assert len(points) == 1
        assert points[0].parser.name == "cp2k"

    def test_explicit_name_string(self, tmp_path):
        _seed_point(tmp_path, "p")
        points = discover_ti_points(tmp_path, parser="cp2k")
        assert len(points) == 1

    def test_auto_sniffs_first_recognised_dir(self, tmp_path):
        """A non-TI helper dir at the top should not break sniffing —
        the sniffer keeps scanning until it finds a recognisable one."""
        (tmp_path / "aaa_helper").mkdir()  # alphabetically first, not a TI dir
        _seed_point(tmp_path, "zzz_real")
        points = discover_ti_points(tmp_path, parser="auto")
        assert len(points) == 1
        assert points[0].directory.name == "zzz_real"

    def test_auto_raises_when_nothing_recognised(self, tmp_path):
        """A directory of non-TI subdirs → ParserInferenceError."""
        (tmp_path / "junk1").mkdir()
        (tmp_path / "junk2").mkdir()
        with pytest.raises(ParserInferenceError):
            discover_ti_points(tmp_path, parser="auto")


# ---------------------------------------------------------------------------
# load_ti_series
# ---------------------------------------------------------------------------


@_skip_no_example
class TestLoadTiSeries:
    def test_returns_xi_lambda_dt_tuples(self, tmp_path):
        # Use a real 2-point root
        points = discover_ti_points(EXAMPLE_TI_ROOT)
        # Restrict to the first 2 to keep test fast
        results = load_ti_series(points[:2])
        assert len(results) == 2
        for xi, lam, dt in results:
            assert isinstance(xi, float)
            assert lam.ndim == 1
            assert lam.size > 0
            assert dt > 0
