"""Tests for ``cli._params.CellAbcParam.collect()`` after the Phase 5
Commit 1 migration to ``engines.cp2k.read_cell``.

Covers (codex Commit 1 v1 MEDIUM 2 fold):
  - .restart branch: prompt_choice -> ".restart", path -> real
    data_example fixture, ctx[K.CELL_ABC] populated.
  - md.inp branch: prompt_choice -> "md.inp", path -> real
    data_example fixture, ctx[K.CELL_ABC] populated.
  - source-vs-suffix mismatch: prompt_choice -> ".restart" but path
    has no .restart suffix -> CellParseError (the explicit guard that
    keeps read_cell's suffix dispatch from masking the user's source
    selection).

Monkeypatches ``cli._params``-local references to ``prompt_choice`` /
``prompt_str_required`` / ``prompt_bool`` so the function-local
``from ..engines.cp2k import read_cell`` path is actually exercised.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from md_analysis.cli._params import CellAbcParam, K
from md_analysis.utils.formats.cp2k.cell import CellParseError

REPO_ROOT = Path(__file__).resolve().parents[3]
RESTART_FIXTURE = REPO_ROOT / "data_example" / "sg" / "angle" / "slowgrowth-1.restart"
MD_INP_FIXTURE = REPO_ROOT / "data_example" / "potential" / "dense" / "md.inp"


def _patch_prompts(
    monkeypatch: pytest.MonkeyPatch,
    *,
    source: str,
    path: str,
    retry: bool = False,
) -> None:
    """Patch the three prompt helpers imported into cli._params."""
    monkeypatch.setattr(
        "md_analysis.cli._params.prompt_choice",
        lambda label, choices, default=None: source,
    )
    monkeypatch.setattr(
        "md_analysis.cli._params.prompt_str_required",
        lambda label: path,
    )
    monkeypatch.setattr(
        "md_analysis.cli._params.prompt_bool",
        lambda label, default=True: retry,
    )


@pytest.mark.skipif(
    not RESTART_FIXTURE.exists(), reason="restart fixture missing"
)
def test_restart_branch_populates_cell_abc(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """User selects '.restart' + provides a real .restart path."""
    _patch_prompts(monkeypatch, source=".restart", path=str(RESTART_FIXTURE))

    ctx: dict = {}
    CellAbcParam().collect(ctx)

    abc = ctx[K.CELL_ABC]
    assert len(abc) == 3
    # Baseline values pinned by test_slowgrowth_parser.py:88 and
    # test_cell_facade.py::TestReadCellFromRestart::test_angle_restart.
    assert abc[0] == pytest.approx(10.2239, rel=1e-4)
    assert abc[1] == pytest.approx(10.2239, rel=1e-4)
    assert abc[2] == pytest.approx(26.422, rel=1e-4)


@pytest.mark.skipif(
    not MD_INP_FIXTURE.exists(), reason="md.inp fixture missing"
)
def test_md_inp_branch_populates_cell_abc(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """User selects 'md.inp' + provides a real md.inp path."""
    _patch_prompts(monkeypatch, source="md.inp", path=str(MD_INP_FIXTURE))

    ctx: dict = {}
    CellAbcParam().collect(ctx)

    abc = ctx[K.CELL_ABC]
    assert len(abc) == 3
    assert all(isinstance(v, float) for v in abc)
    assert all(v > 0 for v in abc)


@pytest.mark.skipif(
    not MD_INP_FIXTURE.exists(), reason="md.inp fixture missing"
)
def test_restart_source_rejects_non_restart_path(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """User selects '.restart' but provides an md.inp path -> raise.

    Without the explicit suffix gate, read_cell would silently
    dispatch to the md.inp parser based on the suffix, masking the
    user's source selection. The CellAbcParam adds an explicit guard
    so the source prompt stays meaningful (codex Commit 1 v1 MEDIUM 1).
    """
    # retry=False so the second attempt is skipped and the second
    # raise (the final CellParseError) propagates out.
    _patch_prompts(
        monkeypatch,
        source=".restart",
        path=str(MD_INP_FIXTURE),
        retry=False,
    )

    ctx: dict = {}
    with pytest.raises(CellParseError, match=".restart suffix"):
        CellAbcParam().collect(ctx)
