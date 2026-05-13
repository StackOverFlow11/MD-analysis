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


def test_public_symbols_importable_from_package_root() -> None:
    """All names listed in ``engines.__all__`` resolve at the package root."""
    from md_analysis.engines import (
        CP2KParser,
        ConstraintMDParser,
        ConstraintMetadata,
        FermiRecord,
        LambdaSeries,
        ParserInferenceError,
        PotentialFrame,
        get_parser,
        infer_parser,
        read_constraint_metadata,
        read_continuous_potential_frames,
        read_distributed_potential_frames,
        read_fermi_series,
        read_lambda_series,
        register_parser,
        resolve_parser,
    )

    # Every symbol is non-None so e.g. typos in the __init__ re-export
    # would surface here even if the import line itself succeeded
    # (which it would for a stale alias).
    assert ConstraintMDParser is not None
    assert ParserInferenceError is not None
    assert CP2KParser is not None
    assert PotentialFrame is not None
    assert ConstraintMetadata is not None
    assert LambdaSeries is not None
    assert FermiRecord is not None
    assert callable(register_parser)
    assert callable(get_parser)
    assert callable(infer_parser)
    assert callable(resolve_parser)
    assert callable(read_constraint_metadata)
    assert callable(read_lambda_series)
    assert callable(read_fermi_series)
    assert callable(read_continuous_potential_frames)
    assert callable(read_distributed_potential_frames)


def test_legacy_dataclass_aliases_resolve_to_renamed_types() -> None:
    """``ColvarRestart`` / ``LagrangeMultLog`` aliases in ``cp2k.colvar``
    must point at the renamed engines-neutral types so existing imports
    transparently see the same class (single source of truth)."""
    from md_analysis.engines import ConstraintMetadata, LambdaSeries
    from md_analysis.utils.formats.cp2k.colvar import (
        ColvarRestart,
        LagrangeMultLog,
    )

    assert ColvarRestart is ConstraintMetadata
    assert LagrangeMultLog is LambdaSeries


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
    """The underlying legacy parser MUST keep returning ``list[dict]``;
    Phase 7b1 intentionally does not migrate ``CenterPotential`` away
    from dict-style access, so the dict shape must be preserved."""
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
