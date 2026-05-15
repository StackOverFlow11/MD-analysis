"""Smoke tests for the ``md_analysis.engines`` package-level facade.

These tests pin the public symbol surface so that future refactors do
not silently break ``from md_analysis.engines import ...`` consumers.
"""

from __future__ import annotations

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]
TI_POINT = REPO_ROOT / "data_example" / "ti" / "double_cv" / "1k" / "ti_target_0.031369"
DENSE_DIR = REPO_ROOT / "data_example" / "potential" / "dense"
MD_OUT = DENSE_DIR / "md.out"
DISTRIBUTED_DIR = REPO_ROOT / "data_example" / "potential" / "distributed"


def test_public_symbols_exported_via_dunder_all() -> None:
    """All names in ``engines.__all__`` resolve at the package root.

    Iterates ``engines.__all__`` and asserts ``hasattr`` + ``is not None``
    for every entry, then spot-checks callability on the known
    function/class surface.  This way the test cannot drift away from
    ``__all__``: any new addition is covered automatically.
    """
    import md_analysis.engines as engines_pkg

    assert engines_pkg.__all__, "engines.__all__ must be non-empty"

    for name in engines_pkg.__all__:
        assert hasattr(engines_pkg, name), (
            f"engines.__all__ lists {name!r} but the package has no "
            f"such attribute (stale re-export?)"
        )
        attr = getattr(engines_pkg, name)
        assert attr is not None, f"engines.{name} resolved to None"

    # Sanity: the well-known function surface is callable. This is a
    # cheap additional check on top of the __all__ scan; new names added
    # in future commits are covered by the loop above without needing
    # to update this list.
    callable_names = (
        "register_parser",
        "get_parser",
        "infer_parser",
        "resolve_parser",
        "read_constraint_metadata",
        "read_lambda_series",
        "read_constraint_run",
        "read_constraint_metadata_from_restart",
        "read_lambda_series_from_log",
        "read_constraint_run_from_files",
        "compute_target_series",
        "read_fermi_series",
        "read_cell",
        "read_continuous_potential_frames",
        "read_distributed_potential_frames",
    )
    for fn_name in callable_names:
        assert callable(getattr(engines_pkg, fn_name)), (
            f"engines.{fn_name} is not callable"
        )


def test_cp2k_submodule_dunder_all_exports_are_live() -> None:
    """All names in ``engines.cp2k.__all__`` resolve at the submodule.

    Prevents future drift where a new facade is added to
    ``engines/cp2k.py`` but not added to its ``__all__`` list
    (Phase 4 Commit 2 added ``read_cell`` without updating
    ``cp2k.__all__`` -- caught by codex Commit 3 v1 MEDIUM and folded
    into Commit 3).
    """
    import md_analysis.engines.cp2k as cp2k_mod

    assert cp2k_mod.__all__, "engines.cp2k.__all__ must be non-empty"

    for name in cp2k_mod.__all__:
        assert hasattr(cp2k_mod, name), (
            f"engines.cp2k.__all__ lists {name!r} but the module "
            f"has no such attribute"
        )
        attr = getattr(cp2k_mod, name)
        assert attr is not None, f"engines.cp2k.{name} resolved to None"


def test_default_registration_makes_cp2k_resolvable() -> None:
    """Importing the package must auto-register the CP2K parser."""
    from md_analysis.engines import CP2KParser, get_parser

    parser = get_parser("cp2k")
    assert isinstance(parser, CP2KParser)
    assert parser.name == "cp2k"


def test_unknown_parser_name_raises_parser_inference_error() -> None:
    """``get_parser`` on an unregistered name raises ``ParserInferenceError``."""
    from md_analysis.engines import ParserInferenceError, get_parser

    with pytest.raises(ParserInferenceError):
        get_parser("definitely-not-a-real-engine")


def test_vasp_placeholder_is_not_auto_registered() -> None:
    """VASP must NOT be in the default registry — calling ``infer_parser``
    on a directory the CP2K parser doesn't recognise should raise rather
    than silently return the VASP placeholder."""
    from md_analysis.engines import ParserInferenceError, get_parser

    with pytest.raises(ParserInferenceError):
        get_parser("vasp")


def test_explicit_vasp_construction_raises_not_implemented(tmp_path) -> None:
    """The VASP placeholder implements the Protocol shape but every
    method must raise ``NotImplementedError`` until Phase 8 fills it in."""
    from md_analysis.engines.vasp import VASPParser

    parser = VASPParser()
    assert parser.name == "vasp"

    with pytest.raises(NotImplementedError):
        parser.is_constraint_directory(tmp_path)

    with pytest.raises(NotImplementedError):
        parser.parse_metadata(tmp_path)

    with pytest.raises(NotImplementedError):
        parser.parse_lambda_series(tmp_path)


# ---------------------------------------------------------------------------
# Phase 7b1 — module-level facade behaviour
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not TI_POINT.exists(), reason="TI fixture missing"
)
def test_read_constraint_metadata_matches_parser() -> None:
    """The facade must return the same result as the underlying parser."""
    from md_analysis.engines import (
        CP2KParser,
        ConstraintMetadata,
        read_constraint_metadata,
    )

    via_facade = read_constraint_metadata(TI_POINT)
    via_parser = CP2KParser().parse_metadata(TI_POINT)

    assert isinstance(via_facade, ConstraintMetadata)
    # Values must match exactly — facade is a thin wrapper, not a new pass.
    assert via_facade.project_name == via_parser.project_name
    assert via_facade.step_start == via_parser.step_start
    assert via_facade.timestep_fs == via_parser.timestep_fs
    assert via_facade.total_steps == via_parser.total_steps
    assert (
        via_facade.colvars.primary.target_au
        == via_parser.colvars.primary.target_au
    )


@pytest.mark.skipif(
    not TI_POINT.exists(), reason="TI fixture missing"
)
def test_read_lambda_series_matches_parser() -> None:
    """The facade must return the same result as the underlying parser."""
    import numpy as np

    from md_analysis.engines import (
        CP2KParser,
        LambdaSeries,
        read_lambda_series,
    )

    via_facade = read_lambda_series(TI_POINT)
    via_parser = CP2KParser().parse_lambda_series(TI_POINT)

    assert isinstance(via_facade, LambdaSeries)
    assert via_facade.n_steps == via_parser.n_steps
    assert via_facade.n_constraints == via_parser.n_constraints
    np.testing.assert_array_equal(via_facade.shake, via_parser.shake)
    np.testing.assert_array_equal(via_facade.rattle, via_parser.rattle)


def test_read_constraint_metadata_accepts_string_path(tmp_path) -> None:
    """``str`` paths should work as well as ``Path`` (path-friendly facade)."""
    from md_analysis.engines import ParserInferenceError, read_constraint_metadata

    # Pointing at an empty directory must raise FileNotFoundError because
    # CP2KParser._find_restart cannot find a *.restart file there. This
    # also implicitly verifies that the facade did accept a str argument
    # without TypeError.
    with pytest.raises(FileNotFoundError):
        read_constraint_metadata(str(tmp_path))


@pytest.mark.skipif(
    not MD_OUT.exists(), reason="md.out fixture missing"
)
def test_read_fermi_series_returns_typed_records() -> None:
    """``read_fermi_series`` must yield ``FermiRecord`` instances with the
    same data as the legacy dict-based parser."""
    from md_analysis.engines import FermiRecord, read_fermi_series
    from md_analysis.utils.formats.cp2k.stdout import parse_md_out_fermi

    typed = read_fermi_series(MD_OUT)
    legacy = parse_md_out_fermi(MD_OUT)

    assert len(typed) == len(legacy)
    assert all(isinstance(r, FermiRecord) for r in typed)

    # Field-by-field agreement
    for record, raw in zip(typed, legacy):
        assert record.step == raw["step"]
        assert record.time_fs == raw["time_fs"]
        assert record.fermi_raw == raw["fermi_raw"]


def test_parse_md_out_fermi_still_returns_dict() -> None:
    """The parser contract remains dict-shaped: ``parse_md_out_fermi``
    is the utils-layer output format that the ``engines.cp2k`` facades
    convert to typed records.  Business now consumes the typed facade
    (``read_fermi_series``); this test pins the parser contract so a
    future utils-layer refactor cannot silently change the shape."""
    from md_analysis.utils.formats.cp2k.stdout import parse_md_out_fermi

    if not MD_OUT.exists():
        pytest.skip("md.out fixture missing")

    records = parse_md_out_fermi(MD_OUT)
    assert isinstance(records, list)
    if records:
        first = records[0]
        assert isinstance(first, dict)
        assert {"step", "time_fs", "fermi_raw"} <= set(first.keys())


def test_fermi_record_from_legacy_dict_round_trip() -> None:
    """``FermiRecord.from_legacy_dict`` bridges a legacy dict 1-to-1."""
    from md_analysis.engines import FermiRecord

    d = {"step": 42, "time_fs": 12.5, "fermi_raw": -0.123456}
    rec = FermiRecord.from_legacy_dict(d)

    assert rec.step == 42
    assert rec.time_fs == 12.5
    assert rec.fermi_raw == pytest.approx(-0.123456)


def test_fermi_record_handles_none_time_fs() -> None:
    """``time_fs`` is allowed to be ``None`` (matches legacy dict semantics)."""
    from md_analysis.engines import FermiRecord

    rec = FermiRecord.from_legacy_dict(
        {"step": 7, "time_fs": None, "fermi_raw": 0.1}
    )
    assert rec.time_fs is None


# ---------------------------------------------------------------------------
# Phase 7b2 — potential-frame facade
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not DENSE_DIR.exists() or not MD_OUT.exists(),
    reason="dense potential fixture missing",
)
def test_read_continuous_potential_frames_matches_legacy_wrapper() -> None:
    """The facade and the legacy thin wrapper must produce the same
    list of ``PotentialFrame`` objects (same steps, same fermi_raw,
    same cube paths). This pins the wrapper as a true pass-through."""
    from md_analysis.electrochemical.potential._frame_source import (
        discover_continuous_frames,
    )
    from md_analysis.engines import read_continuous_potential_frames

    via_facade = read_continuous_potential_frames(
        "md-POTENTIAL-v_hartree-1_*.cube",
        workdir=DENSE_DIR,
        md_out_path=MD_OUT,
        xyz_path=DENSE_DIR / "md-pos-1.xyz",
        center_mode="cell",
    )
    via_wrapper = discover_continuous_frames(
        "md-POTENTIAL-v_hartree-1_*.cube",
        workdir=DENSE_DIR,
        md_out_path=MD_OUT,
        xyz_path=DENSE_DIR / "md-pos-1.xyz",
        center_mode="cell",
    )

    assert len(via_facade) == len(via_wrapper)
    for a, b in zip(via_facade, via_wrapper):
        assert a.step == b.step
        assert a.cube_path == b.cube_path
        assert a.fermi_raw == b.fermi_raw


@pytest.mark.skipif(
    not DISTRIBUTED_DIR.exists(),
    reason="distributed potential fixture missing",
)
def test_read_distributed_potential_frames_matches_legacy_wrapper() -> None:
    """Same parity check for mode B (distributed SP subdirectories)."""
    from md_analysis.electrochemical.potential._frame_source import (
        discover_distributed_frames,
    )
    from md_analysis.engines import read_distributed_potential_frames

    via_facade = read_distributed_potential_frames(
        DISTRIBUTED_DIR,
        center_mode="cell",
    )
    via_wrapper = discover_distributed_frames(
        DISTRIBUTED_DIR,
        center_mode="cell",
    )

    assert len(via_facade) == len(via_wrapper)
    for a, b in zip(via_facade, via_wrapper):
        assert a.step == b.step
        assert a.time_fs == b.time_fs
        assert a.cube_path == b.cube_path
        assert a.fermi_raw == b.fermi_raw


def test_read_distributed_potential_frames_raises_on_missing_root(
    tmp_path,
) -> None:
    """Behaviour preservation: ``FileNotFoundError`` when root_dir doesn't
    exist, matching the legacy wrapper."""
    from md_analysis.engines import read_distributed_potential_frames

    with pytest.raises(FileNotFoundError):
        read_distributed_potential_frames(tmp_path / "missing")


def test_read_distributed_potential_frames_raises_on_empty_root(
    tmp_path,
) -> None:
    """Behaviour preservation: ``FileNotFoundError`` when no subdirs
    match the pattern."""
    from md_analysis.engines import read_distributed_potential_frames

    with pytest.raises(FileNotFoundError):
        read_distributed_potential_frames(tmp_path, center_mode="cell")


# ---------------------------------------------------------------------------
# Phase 8 — VASP placeholders
# ---------------------------------------------------------------------------


def test_vasp_report_module_raises_not_implemented(tmp_path) -> None:
    """``vasp_report`` placeholders must hard-error rather than return a
    silently-empty default."""
    from md_analysis.utils.formats.vasp.report import (
        parse_vasp_report_lambda_series,
        parse_vasp_report_metadata,
    )

    fake_report = tmp_path / "REPORT"

    with pytest.raises(NotImplementedError):
        parse_vasp_report_metadata(fake_report)

    with pytest.raises(NotImplementedError):
        parse_vasp_report_lambda_series(fake_report)


def test_vasp_outcar_module_raises_not_implemented(tmp_path) -> None:
    """``vasp_outcar.parse_outcar_fermi`` must hard-error."""
    from md_analysis.utils.formats.vasp.outcar import parse_outcar_fermi

    fake_outcar = tmp_path / "OUTCAR"

    with pytest.raises(NotImplementedError):
        parse_outcar_fermi(fake_outcar)


def test_vasp_locpot_module_raises_not_implemented(tmp_path) -> None:
    """``vasp_locpot`` placeholders must hard-error."""
    from md_analysis.utils.formats.vasp.locpot import (
        read_locpot,
        read_locpot_plane_avg,
    )

    fake_locpot = tmp_path / "LOCPOT"

    with pytest.raises(NotImplementedError):
        read_locpot(fake_locpot)

    with pytest.raises(NotImplementedError):
        read_locpot_plane_avg(fake_locpot)


def test_vasp_format_modules_not_imported_during_auto_discovery(
    tmp_path,
) -> None:
    """A directory containing VASP marker files must NOT be matched by
    the default ``infer_parser`` — otherwise the placeholder stubs would
    fire when callers pass ``parser='auto'``.

    Builds a directory that looks like a VASP run (REPORT + OUTCAR +
    LOCPOT) and verifies ``infer_parser`` raises ``ParserInferenceError``
    rather than returning a parser that would then crash on use.
    """
    from md_analysis.engines import ParserInferenceError, infer_parser

    (tmp_path / "REPORT").write_text("")
    (tmp_path / "OUTCAR").write_text("")
    (tmp_path / "LOCPOT").write_text("")

    with pytest.raises(ParserInferenceError):
        infer_parser(tmp_path)


# =========================================================================
# Phase 4 Commit 1 — D10 raw -> canonical conversion + D8 ConstraintRun
# =========================================================================


def test_raw_to_constraint_info_field_equality() -> None:
    """_cp2k_raw_to_constraint_info preserves field values 1:1."""
    import numpy as np

    from md_analysis.engines.cp2k import _cp2k_raw_to_constraint_info
    from md_analysis.engines.models import ConstraintInfo
    from md_analysis.utils.formats.cp2k.colvar import Cp2kConstraintInfoRaw

    raw = Cp2kConstraintInfoRaw(
        colvar_id=7,
        target_au=1.5,
        target_growth_au=3.14e-6,
        intermolecular=True,
    )
    canon = _cp2k_raw_to_constraint_info(raw)
    assert isinstance(canon, ConstraintInfo)
    assert canon.colvar_id == 7
    assert canon.target_au == 1.5
    assert canon.target_growth_au == 3.14e-6
    assert canon.intermolecular is True


def test_raw_to_colvar_info_field_equality_and_nested() -> None:
    """_cp2k_raw_to_colvar_info preserves nested ConstraintInfo conversions."""
    from md_analysis.engines.cp2k import _cp2k_raw_to_colvar_info
    from md_analysis.engines.models import ColvarInfo, ConstraintInfo
    from md_analysis.utils.formats.cp2k.colvar import (
        Cp2kColvarInfoRaw,
        Cp2kConstraintInfoRaw,
    )

    raw = Cp2kColvarInfoRaw(
        constraints=(
            Cp2kConstraintInfoRaw(
                colvar_id=1, target_au=2.0,
                target_growth_au=1e-5, intermolecular=False,
            ),
            Cp2kConstraintInfoRaw(
                colvar_id=2, target_au=3.0,
                target_growth_au=-2e-5, intermolecular=True,
            ),
        ),
    )
    canon = _cp2k_raw_to_colvar_info(raw)
    assert isinstance(canon, ColvarInfo)
    assert len(canon) == 2
    assert all(isinstance(c, ConstraintInfo) for c in canon)
    assert canon.primary.colvar_id == 1
    assert canon[2].colvar_id == 2
    assert canon[2].intermolecular is True


def test_raw_to_constraint_metadata_field_equality() -> None:
    """_cp2k_raw_to_constraint_metadata preserves all 9 fields."""
    from md_analysis.engines.cp2k import _cp2k_raw_to_constraint_metadata
    from md_analysis.engines.models import ConstraintMetadata
    from md_analysis.utils.formats.cp2k.colvar import (
        Cp2kColvarInfoRaw,
        Cp2kConstraintInfoRaw,
        Cp2kConstraintMetadataRaw,
    )

    raw = Cp2kConstraintMetadataRaw(
        project_name="phase4_test",
        step_start=100,
        time_start_fs=50.0,
        timestep_fs=0.5,
        total_steps=200,
        colvars=Cp2kColvarInfoRaw(constraints=(
            Cp2kConstraintInfoRaw(
                colvar_id=1, target_au=1.5,
                target_growth_au=1e-5, intermolecular=True,
            ),
        )),
        lagrange_filename="cf.dat",
        cell_abc_ang=(10.0, 10.0, 30.0),
        fixed_atom_indices=(1, 2, 3),
    )
    canon = _cp2k_raw_to_constraint_metadata(raw)
    assert isinstance(canon, ConstraintMetadata)
    assert canon.project_name == "phase4_test"
    assert canon.step_start == 100
    assert canon.time_start_fs == 50.0
    assert canon.timestep_fs == 0.5
    assert canon.total_steps == 200
    assert canon.lagrange_filename == "cf.dat"
    assert canon.cell_abc_ang == (10.0, 10.0, 30.0)
    assert canon.fixed_atom_indices == (1, 2, 3)
    assert canon.colvars.primary.colvar_id == 1
    assert canon.colvars.primary.target_au == 1.5


def test_raw_to_lambda_series_field_equality() -> None:
    """_cp2k_raw_to_lambda_series preserves shake/rattle arrays + counts."""
    import numpy as np

    from md_analysis.engines.cp2k import _cp2k_raw_to_lambda_series
    from md_analysis.engines.models import LambdaSeries
    from md_analysis.utils.formats.cp2k.colvar import Cp2kLambdaSeriesRaw

    shake = np.array([1.0, 2.0, 3.0])
    rattle = np.array([0.1, 0.2, 0.3])
    raw = Cp2kLambdaSeriesRaw(
        shake=shake, rattle=rattle, n_steps=3, n_constraints=1,
    )
    canon = _cp2k_raw_to_lambda_series(raw)
    assert isinstance(canon, LambdaSeries)
    assert canon.n_steps == 3
    assert canon.n_constraints == 1
    np.testing.assert_array_equal(canon.shake, shake)
    np.testing.assert_array_equal(canon.rattle, rattle)
    # canonical accessor (LambdaSeries property) works
    np.testing.assert_array_equal(canon.collective_shake, shake)


def test_constraint_run_target_series_au_literal_formula() -> None:
    """ConstraintRun.target_series_au matches the literal documented formula.

    Self-contained literal-formula reference:
    references the formula
        xi(k) = target_au + (k - step_start) * target_growth_au * dt_au
    with dt_au = timestep_fs / AU_TIME_TO_FS,
    using AU_TIME_TO_FS = 0.02418884326585 (CODATA).
    """
    import numpy as np

    from md_analysis.engines.models import (
        ColvarInfo,
        ConstraintInfo,
        ConstraintMetadata,
        ConstraintRun,
        LambdaSeries,
    )

    meta = ConstraintMetadata(
        project_name="formula_test",
        step_start=100,
        time_start_fs=50.0,
        timestep_fs=0.5,
        total_steps=10,
        colvars=ColvarInfo(constraints=(
            ConstraintInfo(
                colvar_id=1, target_au=2.0,
                target_growth_au=1e-5, intermolecular=True,
            ),
        )),
        lagrange_filename=None,
        cell_abc_ang=(10.0, 10.0, 30.0),
        fixed_atom_indices=None,
    )
    ls = LambdaSeries(
        shake=np.zeros(10), rattle=np.zeros(10),
        n_steps=10, n_constraints=1,
    )
    run = ConstraintRun(metadata=meta, lambda_series=ls)

    # Literal-formula reference (CODATA AU_TIME_TO_FS)
    AU_TIME_TO_FS = 0.02418884326585
    dt_au = 0.5 / AU_TIME_TO_FS
    k = np.arange(10)
    expected = 2.0 + (k - 100) * 1e-5 * dt_au

    np.testing.assert_allclose(
        run.target_series_au(), expected, atol=1e-12,
    )
