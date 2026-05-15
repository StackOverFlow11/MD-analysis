"""Unit tests for ``md_analysis.workflows.charge.run_surface_charge``.

Phase 6.3 extends the workflow with 4 ``potential_*`` parameters and a
``target_side`` argument that toggles per-side analysis (and the
matching ``<method>_<side>`` sub-directory).  These tests monkeypatch
the underlying ``surface_charge_analysis`` business function and pin:

  - the new parameters are forwarded verbatim
  - ``target_side=None`` produces ``output_dir/<method>/``
  - ``target_side="aligned"`` / ``"opposed"`` produce
    ``output_dir/<method>_<side>/``
  - ``WorkflowResult.metadata["target_side"]`` reflects the request
  - sigma / phi metadata are surfaced from the business return object
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest

from md_analysis.workflows import WorkflowResult, run_surface_charge


@dataclass
class _FakeResult:
    """Minimal stand-in for SurfaceChargeResult."""

    csv_path: Path
    n_frames: int = 17
    sigma_aligned_mean: float = -1.25
    sigma_aligned_std: float = 0.05
    sigma_opposed_mean: float = 1.30
    sigma_opposed_std: float = 0.06
    phi_cumavg_last: float | None = None
    phi_reference: str | None = None


def _patch_business(monkeypatch: pytest.MonkeyPatch, captured: dict[str, Any]) -> None:
    """Patch surface_charge_analysis to capture kwargs + return a fake."""

    def fake(root_dir: Path, **kwargs: Any) -> _FakeResult:
        captured["root_dir"] = root_dir
        captured.update(kwargs)
        out_dir = Path(kwargs["output_dir"])
        out_dir.mkdir(parents=True, exist_ok=True)
        csv_path = out_dir / "surface_charge.csv"
        csv_path.write_text("")
        return _FakeResult(csv_path=csv_path)

    monkeypatch.setattr(
        "md_analysis.electrochemical.charge.surface_charge_analysis",
        fake,
    )


class TestRunSurfaceChargeTwoSided:
    """target_side=None (default) -> output_dir/<method>/."""

    def test_counterion_subdir(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        captured: dict[str, Any] = {}
        _patch_business(monkeypatch, captured)

        result = run_surface_charge(
            output_dir=tmp_path / "charge",
            root_dir=tmp_path / "data",
            method="counterion",
        )

        # Sub-dir composition: <output_dir>/<method>/
        expected_dir = tmp_path / "charge" / "counterion"
        assert expected_dir.exists()
        assert captured["output_dir"] == expected_dir

        # WorkflowResult contract
        assert isinstance(result, WorkflowResult)
        assert result.name == "surface_charge"
        assert result.output_dir == expected_dir
        assert set(result.artifacts) == {"charge_csv", "charge_png"}
        assert result.artifacts["charge_csv"].name == "surface_charge.csv"

        # metadata
        assert result.metadata["method"] == "counterion"
        assert result.metadata["target_side"] is None
        assert result.metadata["n_frames"] == 17
        assert result.metadata["sigma_aligned_mean"] == pytest.approx(-1.25)

    def test_layer_subdir(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        captured: dict[str, Any] = {}
        _patch_business(monkeypatch, captured)

        result = run_surface_charge(
            output_dir=tmp_path / "charge",
            root_dir=tmp_path / "data",
            method="layer",
        )

        expected_dir = tmp_path / "charge" / "layer"
        assert captured["output_dir"] == expected_dir
        assert result.output_dir == expected_dir
        assert result.metadata["target_side"] is None


class TestRunSurfaceChargeSingleSide:
    """target_side="aligned"/"opposed" -> output_dir/<method>_<side>/."""

    @pytest.mark.parametrize("method", ["counterion", "layer"])
    @pytest.mark.parametrize("side", ["aligned", "opposed"])
    def test_method_side_subdir(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
        method: str,
        side: str,
    ) -> None:
        captured: dict[str, Any] = {}
        _patch_business(monkeypatch, captured)

        result = run_surface_charge(
            output_dir=tmp_path / "charge",
            root_dir=tmp_path / "data",
            method=method,
            target_side=side,
        )

        expected_dir = tmp_path / "charge" / f"{method}_{side}"
        assert expected_dir.exists()
        assert captured["output_dir"] == expected_dir
        assert captured["target_side"] == side

        assert result.output_dir == expected_dir
        assert result.metadata["target_side"] == side
        assert result.metadata["method"] == method


class TestRunSurfaceChargeInvalidSide:
    """Invalid target_side rejected BEFORE any filesystem side-effect.

    Phase 6.3 follow-up: the facade composes
    ``<output_dir>/<method>_<side>/`` and calls ``mkdir`` before
    delegating to the business layer, where target_side is actually
    validated.  Without an early check, a bogus value (e.g.
    ``target_side="bad"``) leaves an empty ``<method>_bad/`` directory
    on disk before the ValueError surfaces.  These tests pin that the
    facade rejects invalid sides up-front and emits the same wording
    as the business layer.
    """

    def test_invalid_target_side_raises_before_mkdir(
        self, tmp_path: Path
    ) -> None:
        base = tmp_path / "charge"
        with pytest.raises(ValueError) as exc:
            run_surface_charge(
                output_dir=base,
                root_dir=tmp_path / "data",
                method="counterion",
                target_side="bad",
            )
        # byte-equal error wording — pins the contract that facade and
        # business layer raise identical messages
        assert str(exc.value) == (
            "target_side must be 'aligned', 'opposed', or None, got 'bad'"
        )
        # No <method>_bad/ subdir leaked
        assert not (base / "counterion_bad").exists()
        # base_dir itself shouldn't be created either: mkdir runs only
        # after validation passes
        assert not base.exists()


class TestRunSurfaceChargePotentialKwargs:
    """The 4 potential_* kwargs are forwarded verbatim to the underlying
    surface_charge_analysis call."""

    def test_potential_kwargs_forwarded(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        captured: dict[str, Any] = {}
        _patch_business(monkeypatch, captured)

        run_surface_charge(
            output_dir=tmp_path / "charge",
            root_dir=tmp_path / "data",
            method="counterion",
            potential_reference="RHE",
            potential_pH=2.5,
            potential_temperature_K=310.15,
            potential_phi_pzc=-0.42,
        )

        assert captured["potential_reference"] == "RHE"
        assert captured["potential_pH"] == pytest.approx(2.5)
        assert captured["potential_temperature_K"] == pytest.approx(310.15)
        assert captured["potential_phi_pzc"] == pytest.approx(-0.42)

    def test_potential_phi_pzc_default_is_none(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        captured: dict[str, Any] = {}
        _patch_business(monkeypatch, captured)

        run_surface_charge(
            output_dir=tmp_path / "charge",
            root_dir=tmp_path / "data",
        )

        assert captured["potential_phi_pzc"] is None
        assert captured["potential_reference"] == "SHE"
