"""Tests for ``run_ti_full_from_root`` agent wrapper and slice helper."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from md_analysis.enhanced_sampling.constrained_ti.workflow import (
    TIFullAnalysisReport,
    _parse_point_slice,
    run_ti_full_from_root,
)


_EXAMPLE_TI_ROOT = Path(__file__).resolve().parents[4] / "data_example" / "ti" / "double_cv" / "1k"
_EXAMPLE_POINT = _EXAMPLE_TI_ROOT / "ti_target_0.031369"


def _seed_real_point(parent: Path, name: str) -> Path:
    """Copy a real CP2K TI-point under *parent*/<name>."""
    dst = parent / name
    shutil.copytree(_EXAMPLE_POINT, dst)
    return dst

_skip_no_example = pytest.mark.skipif(
    not _EXAMPLE_TI_ROOT.is_dir(),
    reason=f"example TI data not available at {_EXAMPLE_TI_ROOT}",
)


# ---------------------------------------------------------------------------
# _parse_point_slice
# ---------------------------------------------------------------------------


class TestParsePointSlice:
    """Accept valid Python slice syntax; reject malformed forms."""

    def test_start_stop(self):
        assert _parse_point_slice("0:2") == slice(0, 2, None)

    def test_stop_only(self):
        assert _parse_point_slice(":2") == slice(None, 2, None)

    def test_step_only(self):
        assert _parse_point_slice("::2") == slice(None, None, 2)

    def test_all_three(self):
        assert _parse_point_slice("1:5:2") == slice(1, 5, 2)

    def test_negative_indices(self):
        assert _parse_point_slice("-3:-1") == slice(-3, -1, None)

    def test_bare_integer_rejected(self):
        with pytest.raises(ValueError, match="must contain at least one ':'"):
            _parse_point_slice("2")

    def test_four_parts_rejected(self):
        with pytest.raises(ValueError, match="expected 2 or 3 colon-separated parts"):
            _parse_point_slice("1:2:3:4")

    def test_non_integer_part_rejected(self):
        with pytest.raises(ValueError, match="not an integer"):
            _parse_point_slice("0:abc")


# ---------------------------------------------------------------------------
# run_ti_full_from_root — input validation (no filesystem required)
# ---------------------------------------------------------------------------


class TestRunTIFullInputValidation:
    def test_missing_root_dir_raises_file_not_found(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="root_dir"):
            run_ti_full_from_root(
                root_dir=tmp_path / "does_not_exist",
                output_dir=tmp_path / "out",
            )

    def test_unknown_parser_raises(self, tmp_path):
        (tmp_path / "root").mkdir()
        from md_analysis.enhanced_sampling._parsers import (
            ParserInferenceError,
        )
        with pytest.raises(ParserInferenceError):
            run_ti_full_from_root(
                root_dir=tmp_path / "root",
                output_dir=tmp_path / "out",
                parser="bogus_engine",
            )

    @_skip_no_example
    def test_invalid_point_slice_raises_value_error(self, tmp_path):
        with pytest.raises(ValueError, match="Invalid point_slice"):
            run_ti_full_from_root(
                root_dir=_EXAMPLE_TI_ROOT,
                output_dir=tmp_path / "out",
                point_slice="2",
            )

    @_skip_no_example
    def test_too_few_points_after_slice_raises_value_error(self, tmp_path):
        with pytest.raises(ValueError, match="at least 2 TI points"):
            run_ti_full_from_root(
                root_dir=_EXAMPLE_TI_ROOT,
                output_dir=tmp_path / "out",
                point_slice="0:1",  # only 1 point remaining
            )


# ---------------------------------------------------------------------------
# Strict discovery — matched-but-incomplete TI directories
# ---------------------------------------------------------------------------


def _make_minimal_ti_dir(
    parent: Path,
    name: str,
    *,
    with_restart: bool = True,
    with_log: bool = True,
) -> Path:
    """Create a bare ti_target_<cv>/ directory with optional placeholder files."""
    d = parent / name
    d.mkdir(parents=True)
    if with_restart:
        (d / "cMD-1.restart").write_text("placeholder restart\n")
    if with_log:
        (d / "cMD-1.LagrangeMultLog").write_text("placeholder log\n")
    return d


class TestStrictDiscovery:
    """discover_ti_points(strict=False) keeps skip-and-warn;
    strict=True raises FileNotFoundError for incomplete matched dirs."""

    def test_lenient_skips_unparseable_dir(self, tmp_path):
        """Default strict=False with a name-based filter: dirs whose
        files exist but contain bad content are silently skipped."""
        from md_analysis.enhanced_sampling.constrained_ti.io import (
            discover_ti_points,
        )

        # Real directory + 2 placeholder dirs. With dir_filter="ti_target_*"
        # all three are name-selected, but only the real one's content
        # parses; placeholders are skipped under lenient mode.
        _seed_real_point(tmp_path, "ti_target_0.031369")
        _make_minimal_ti_dir(tmp_path, "ti_target_0.200000")
        _make_minimal_ti_dir(tmp_path, "ti_target_0.300000")

        points = discover_ti_points(
            tmp_path, dir_filter="ti_target_*", strict=False,
        )
        # Only the real one survives (placeholders fail metadata parse).
        assert len(points) == 1
        assert points[0].directory.name == "ti_target_0.031369"

    def test_strict_raises_for_unparseable_dir(self, tmp_path):
        """strict=True + dir_filter glob: a name-selected but unparseable
        dir surfaces as FileNotFoundError."""
        from md_analysis.enhanced_sampling.constrained_ti.io import (
            discover_ti_points,
        )

        _seed_real_point(tmp_path, "ti_target_0.031369")
        _make_minimal_ti_dir(tmp_path, "ti_target_0.200000")

        with pytest.raises(FileNotFoundError, match="ti_target_0.200000"):
            discover_ti_points(
                tmp_path, dir_filter="ti_target_*", strict=True,
            )

    def test_strict_raises_for_missing_log(self, tmp_path):
        from md_analysis.enhanced_sampling.constrained_ti.io import (
            discover_ti_points,
        )

        _make_minimal_ti_dir(tmp_path, "ti_target_0.100000", with_log=False)

        with pytest.raises(FileNotFoundError):
            discover_ti_points(
                tmp_path, dir_filter="ti_target_*", strict=True,
            )

    def test_run_ti_full_surfaces_missing_file_as_filenotfound(self, tmp_path):
        """run_ti_full_from_root uses strict=True internally; an
        unparseable matched dir must raise FileNotFoundError (not be
        silently skipped). ``output_dir`` must not have been created."""
        root = tmp_path / "root"
        root.mkdir()
        _seed_real_point(root, "ti_target_0.031369")
        _make_minimal_ti_dir(root, "ti_target_0.200000")  # bad content
        _seed_real_point(root, "ti_target_0.302356")

        with pytest.raises(FileNotFoundError, match="ti_target_0.200000"):
            run_ti_full_from_root(
                root_dir=root,
                output_dir=tmp_path / "out",
                dir_filter="ti_target_*",
            )
        # Strict discovery must fail before any output side effects.
        assert not (tmp_path / "out").exists()

    def test_dispatch_maps_missing_file_to_file_not_found(self, tmp_path):
        """Dispatch-level mirror of the above."""
        from md_analysis.agent import dispatch

        root = tmp_path / "root"
        root.mkdir()
        _seed_real_point(root, "ti_target_0.031369")
        _make_minimal_ti_dir(root, "ti_target_0.200000", with_log=False)

        result = dispatch("ti_full_analysis", {
            "root_dir": str(root),
            "output_dir": str(tmp_path / "out"),
            "dir_filter": "ti_target_*",
        })
        assert not result.success
        assert result.error_type == "file_not_found"


# ---------------------------------------------------------------------------
# run_ti_full_from_root — end-to-end on real example data
# ---------------------------------------------------------------------------


@_skip_no_example
class TestRunTIFullEndToEnd:
    def test_returns_report_with_artifacts(self, tmp_path):
        report = run_ti_full_from_root(
            root_dir=_EXAMPLE_TI_ROOT,
            output_dir=tmp_path / "out",
            point_slice="0:2",
        )
        assert isinstance(report, TIFullAnalysisReport)
        assert report.n_points == 2
        assert report.convergence_csv.is_file()
        assert report.free_energy_csv.is_file()
        assert report.free_energy_png.is_file()
        assert len(report.diagnostics_pngs) == 2
        for p in report.diagnostics_pngs:
            assert p.is_file()

    def test_per_point_has_all_required_fields(self, tmp_path):
        report = run_ti_full_from_root(
            root_dir=_EXAMPLE_TI_ROOT,
            output_dir=tmp_path / "out",
            point_slice="0:2",
        )
        required = {
            "point_index", "xi", "n_analyzed",
            "time_start_fs", "time_end_fs", "time_total_fs",
            "tau_corr", "n_eff", "sem_final_au", "sem_max_au",
            "geweke_z", "drift_D", "passed", "failure_reasons",
        }
        for pt in report.per_point:
            assert set(pt.keys()) == required

    def test_per_point_is_json_serializable(self, tmp_path):
        report = run_ti_full_from_root(
            root_dir=_EXAMPLE_TI_ROOT,
            output_dir=tmp_path / "out",
            point_slice="0:2",
        )
        # Should not raise; would catch numpy scalars leaking through.
        json.dumps(list(report.per_point))

    def test_time_total_matches_n_analyzed_times_dt(self, tmp_path):
        """time_total_fs = n_analyzed * dt (within floating tolerance)."""
        report = run_ti_full_from_root(
            root_dir=_EXAMPLE_TI_ROOT,
            output_dir=tmp_path / "out",
            point_slice="0:2",
        )
        for pt in report.per_point:
            # dt ≈ time_total / n_analyzed; check consistency
            if pt["n_analyzed"] > 0:
                assert pt["time_total_fs"] > 0

    def test_csv_filenames(self, tmp_path):
        report = run_ti_full_from_root(
            root_dir=_EXAMPLE_TI_ROOT,
            output_dir=tmp_path / "out",
            point_slice="0:2",
        )
        assert report.convergence_csv.name == "ti_convergence_report.csv"
        assert report.free_energy_csv.name == "ti_free_energy.csv"
        assert report.free_energy_png.name == "ti_free_energy.png"


# ---------------------------------------------------------------------------
# Dispatch-level integration
# ---------------------------------------------------------------------------


@_skip_no_example
class TestDispatchTIFullAnalysis:
    def test_dispatch_success(self, tmp_path):
        from md_analysis.agent import dispatch

        result = dispatch("ti_full_analysis", {
            "root_dir": str(_EXAMPLE_TI_ROOT),
            "output_dir": str(tmp_path / "out"),
            "point_slice": "0:2",
        })
        assert result.success, f"errors: {result.errors}"
        assert "convergence_csv" in result.outputs
        assert "free_energy_csv" in result.outputs
        assert "free_energy_png" in result.outputs
        assert "diagnostics_png_0" in result.outputs

    def test_dispatch_summary_has_all_metrics(self, tmp_path):
        from md_analysis.agent import dispatch

        result = dispatch("ti_full_analysis", {
            "root_dir": str(_EXAMPLE_TI_ROOT),
            "output_dir": str(tmp_path / "out"),
            "point_slice": "0:2",
        })
        assert result.success
        expected = {"n_points", "delta_A_eV", "sigma_A_eV",
                    "all_passed", "failing_indices", "per_point"}
        assert expected <= set(result.summary.keys())

    def test_dispatch_summary_is_json_serializable(self, tmp_path):
        from md_analysis.agent import dispatch

        result = dispatch("ti_full_analysis", {
            "root_dir": str(_EXAMPLE_TI_ROOT),
            "output_dir": str(tmp_path / "out"),
            "point_slice": "0:2",
        })
        assert result.success
        json.dumps(result.summary)   # raises TypeError if not serializable

    def test_dispatch_missing_root_file_not_found(self, tmp_path):
        from md_analysis.agent import dispatch

        result = dispatch("ti_full_analysis", {
            "root_dir": str(tmp_path / "nonexistent"),
            "output_dir": str(tmp_path / "out"),
            "point_slice": "0:2",
        })
        assert not result.success
        assert result.error_type == "file_not_found"

    def test_dispatch_invalid_point_slice_validation(self, tmp_path):
        from md_analysis.agent import dispatch

        result = dispatch("ti_full_analysis", {
            "root_dir": str(_EXAMPLE_TI_ROOT),
            "output_dir": str(tmp_path / "out"),
            "point_slice": "2",   # no colon → invalid
        })
        assert not result.success
        assert result.error_type == "validation"

    def test_dispatch_too_few_points_validation(self, tmp_path):
        from md_analysis.agent import dispatch

        result = dispatch("ti_full_analysis", {
            "root_dir": str(_EXAMPLE_TI_ROOT),
            "output_dir": str(tmp_path / "out"),
            "point_slice": "0:1",   # < 2 points after slice
        })
        assert not result.success
        assert result.error_type == "validation"
