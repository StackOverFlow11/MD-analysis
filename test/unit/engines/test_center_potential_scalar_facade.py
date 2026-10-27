"""Tests for ``CenterPotentialScalarFrame`` + ``read_center_potential_scalar_frame``.

Phase 4 Commit 3 (D11).  Covers:

  - dataclass default fields
  - facade numerical equality with ``slab_average_potential_ev``
  - Fermi Hartree -> eV via ``HA_TO_EV``
  - ``fermi_raw=None`` pass-through
  - ``center_source`` metadata pass-through (no runtime validation)
  - ``step`` / ``time_fs`` pass-through
  - ``phi_z_std_ev`` / ``n_slices`` come from the ``slab_average_potential_ev``
    info dict
  - keyword-only signature constraint
  - explicit ``center_z_ang=None`` None guard (Phase 3 §5.4 lock; codex
    Commit 3 v1 MEDIUM 2)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from md_analysis.engines import (
    CenterPotentialScalarFrame,
    PotentialFrame,
    read_center_potential_scalar_frame,
)
from md_analysis.utils.constants import HA_TO_EV
from md_analysis.utils.formats.common.cube import (
    read_cube_header_and_values,
    slab_average_potential_ev,
)

REPO_ROOT = Path(__file__).resolve().parents[3]
DENSE_DIR = REPO_ROOT / "data_example" / "potential" / "dense"
SAMPLE_CUBE = DENSE_DIR / "md-POTENTIAL-v_hartree-1_0.cube"


# =========================================================================
# Dataclass surface
# =========================================================================


class TestDataclassDefaults:
    """``phi_z_std_ev`` and ``n_slices`` are the only optional fields."""

    def test_optional_fields_default_to_none(self) -> None:
        f = CenterPotentialScalarFrame(
            step=0,
            time_fs=None,
            center_source="manual",
            center_z_ang=5.0,
            slab_thickness_ang=2.0,
            phi_center_ev=-3.0,
            fermi_level_ev=None,
        )
        assert f.phi_z_std_ev is None
        assert f.n_slices is None

    def test_all_fields_populated_round_trip(self) -> None:
        f = CenterPotentialScalarFrame(
            step=42,
            time_fs=12.5,
            center_source="interface",
            center_z_ang=15.0,
            slab_thickness_ang=4.0,
            phi_center_ev=-1.25,
            fermi_level_ev=-0.5,
            phi_z_std_ev=0.01,
            n_slices=17,
        )
        assert f.step == 42
        assert f.time_fs == 12.5
        assert f.center_source == "interface"
        assert f.center_z_ang == 15.0
        assert f.slab_thickness_ang == 4.0
        assert f.phi_center_ev == -1.25
        assert f.fermi_level_ev == -0.5
        assert f.phi_z_std_ev == 0.01
        assert f.n_slices == 17


# =========================================================================
# Helpers — build a real PotentialFrame from data_example/potential/dense
# =========================================================================


def _build_frame(*, fermi_raw: float | None) -> PotentialFrame:
    header, values = read_cube_header_and_values(SAMPLE_CUBE)
    return PotentialFrame(
        step=0,
        time_fs=0.0,
        cube_path=SAMPLE_CUBE,
        header=header,
        values=values,
        fermi_raw=fermi_raw,
        atoms=None,
    )


# =========================================================================
# facade — numerical equality with slab_average_potential_ev
# =========================================================================


@pytest.mark.skipif(
    not SAMPLE_CUBE.exists(), reason="dense potential fixture missing"
)
class TestFacadeNumericalEquality:
    """Facade outputs match the underlying utils helper byte-for-byte."""

    def test_phi_center_ev_matches_helper(self) -> None:
        frame = _build_frame(fermi_raw=-0.18)
        scalar = read_center_potential_scalar_frame(
            frame, center_z_ang=15.0, slab_thickness_ang=4.0,
        )
        phi_direct, _ = slab_average_potential_ev(
            frame.header, frame.values, 4.0, z_center_ang=15.0,
        )
        # Bit-for-bit equality (facade is a thin reshape over the helper).
        assert scalar.phi_center_ev == phi_direct

    def test_phi_z_std_ev_and_n_slices_from_info_dict(self) -> None:
        frame = _build_frame(fermi_raw=-0.18)
        scalar = read_center_potential_scalar_frame(
            frame, center_z_ang=15.0, slab_thickness_ang=4.0,
        )
        _, info = slab_average_potential_ev(
            frame.header, frame.values, 4.0, z_center_ang=15.0,
        )
        assert scalar.phi_z_std_ev == info["phi_z_std_ev"]
        assert scalar.n_slices == info["n_slices"]


# =========================================================================
# facade — Fermi Hartree -> eV
# =========================================================================


@pytest.mark.skipif(
    not SAMPLE_CUBE.exists(), reason="dense potential fixture missing"
)
class TestFermiHaToEv:
    """``fermi_level_ev = frame.fermi_raw * HA_TO_EV`` (literal constant ref)."""

    def test_negative_fermi_raw(self) -> None:
        frame = _build_frame(fermi_raw=-0.183456)
        scalar = read_center_potential_scalar_frame(
            frame, center_z_ang=15.0, slab_thickness_ang=4.0,
        )
        assert scalar.fermi_level_ev == pytest.approx(
            -0.183456 * HA_TO_EV, abs=1e-12,
        )

    def test_zero_fermi_raw_not_none(self) -> None:
        # Edge case: 0.0 is a valid Fermi value; must NOT trip the
        # ``frame.fermi_raw is not None`` branch into None.
        frame = _build_frame(fermi_raw=0.0)
        scalar = read_center_potential_scalar_frame(
            frame, center_z_ang=15.0, slab_thickness_ang=4.0,
        )
        assert scalar.fermi_level_ev == 0.0

    def test_none_fermi_raw_passes_through(self) -> None:
        frame = _build_frame(fermi_raw=None)
        scalar = read_center_potential_scalar_frame(
            frame, center_z_ang=15.0, slab_thickness_ang=4.0,
        )
        assert scalar.fermi_level_ev is None


# =========================================================================
# facade — metadata pass-through (no runtime validation)
# =========================================================================


@pytest.mark.skipif(
    not SAMPLE_CUBE.exists(), reason="dense potential fixture missing"
)
class TestMetadataPassThrough:
    @pytest.mark.parametrize("source", ["manual", "interface", "cell"])
    def test_center_source_pass_through(self, source: str) -> None:
        frame = _build_frame(fermi_raw=-0.18)
        scalar = read_center_potential_scalar_frame(
            frame,
            center_z_ang=15.0,
            slab_thickness_ang=4.0,
            center_source=source,
        )
        assert scalar.center_source == source

    def test_center_source_unknown_value_not_validated(self) -> None:
        """Facade does NOT enforce the enum at runtime (Phase 3 §5.2 note):
        any string is accepted as metadata."""
        frame = _build_frame(fermi_raw=-0.18)
        scalar = read_center_potential_scalar_frame(
            frame,
            center_z_ang=15.0,
            slab_thickness_ang=4.0,
            center_source="bogus-not-an-enum",
        )
        assert scalar.center_source == "bogus-not-an-enum"

    def test_step_and_time_fs_pass_through(self) -> None:
        header, values = read_cube_header_and_values(SAMPLE_CUBE)
        frame = PotentialFrame(
            step=42,
            time_fs=12.5,
            cube_path=SAMPLE_CUBE,
            header=header,
            values=values,
            fermi_raw=-0.18,
            atoms=None,
        )
        scalar = read_center_potential_scalar_frame(
            frame, center_z_ang=15.0, slab_thickness_ang=4.0,
        )
        assert scalar.step == 42
        assert scalar.time_fs == 12.5


# =========================================================================
# facade — signature constraints
# =========================================================================


@pytest.mark.skipif(
    not SAMPLE_CUBE.exists(), reason="dense potential fixture missing"
)
class TestSignatureConstraints:
    def test_center_z_ang_is_keyword_only(self) -> None:
        frame = _build_frame(fermi_raw=-0.18)
        with pytest.raises(TypeError):
            # center_z_ang / slab_thickness_ang are keyword-only.
            read_center_potential_scalar_frame(frame, 15.0, 4.0)  # type: ignore[misc]

    def test_explicit_none_center_z_ang_rejected(self) -> None:
        """Phase 3 §5.4 lock (codex Commit 3 v1 MEDIUM 2):
        facade explicitly refuses ``center_z_ang=None`` so the
        underlying ``slab_average_potential_ev`` cell-center fallback
        cannot leak through this entry point."""
        frame = _build_frame(fermi_raw=-0.18)
        with pytest.raises(TypeError, match="center_z_ang"):
            read_center_potential_scalar_frame(
                frame, center_z_ang=None, slab_thickness_ang=4.0,
            )
