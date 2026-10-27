"""Ground-truth regression for CV target reconstruction at k=0.

End-to-end check that ties three things together:

1. ``CP2KParser._find_restart`` picks the rolling
   ``<PROJECT>-<RUN>.restart`` (not a RESTART_HISTORY snapshot).
2. ``parse_colvar_restart`` extracts the right (target_au,
   step_start, target_growth_au, timestep_fs) tuple.
3. ``compute_target_series`` reconstructs ``xi(k=0)`` matching the
   user-specified initial TARGET in the original CP2K input file.

If ``_find_restart`` were silently returning a snapshot instead of
the rolling file, this test would still pass mathematically because
``(target_au, step_start)`` co-vary along a snapshot — but it would
fail loudly if the picker returned a ``.bak`` / ``RESTART.wfn`` or
raised. Hence this is the *absolute-anchor* regression that
complements ``test_find_restart.py``.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from md_analysis.engines.cp2k import (
    compute_target_series,
    read_constraint_metadata,
)
from md_analysis.utils.constants import BOHR_TO_ANG

REPO_ROOT = Path(__file__).resolve().parents[3]

# Each fixture has: bare slowgrowth-1.restart + 11 RESTART_HISTORY
# snapshots + a sg.inp whose &COLLECTIVE block defines the initial
# TARGET. Numbers below are taken verbatim from the inp files.
ANGLE_DIR = REPO_ROOT / "data_example" / "sg" / "angle"
DIST_DIR = REPO_ROOT / "data_example" / "sg" / "distance_combinedCV"


@pytest.mark.parametrize(
    "directory, inp_target_user_unit, unit",
    [
        (ANGLE_DIR, 122.504, "deg"),        # TARGET [deg] 122.504
        (DIST_DIR, 0.58, "angstrom"),       # TARGET [angstrom] 0.58
    ],
    ids=["angle_deg", "distance_combinedCV_angstrom"],
)
class TestComputeTargetSeriesInitialValue:
    def test_xi_at_k_zero_matches_inp_target(
        self, directory, inp_target_user_unit, unit,
    ):
        """ξ(k=0) reconstructed from the bare restart must equal the
        initial TARGET written in the &COLLECTIVE block of sg.inp."""
        if not directory.exists():
            pytest.skip(f"SG fixture missing: {directory}")

        metadata = read_constraint_metadata(directory)
        # n_steps=1 is enough to evaluate xi at k=0; the rest of the
        # series is not under test here.
        xi = compute_target_series(metadata, n_steps=1)
        assert xi.shape == (1,)

        xi_at_zero_au = float(xi[0])

        if unit == "angstrom":
            xi_at_zero_user_unit = xi_at_zero_au * BOHR_TO_ANG
        elif unit == "deg":
            xi_at_zero_user_unit = math.degrees(xi_at_zero_au)
        else:  # pragma: no cover
            raise ValueError(f"Unknown unit {unit!r}")

        # Tolerance: TARGET in the inp is written to 3 decimals
        # (122.504, 0.58) and the growth term has been integrated
        # over thousands of steps; a 1e-3 absolute tolerance in the
        # user unit comfortably brackets round-off without masking
        # real reconstruction errors.
        assert math.isclose(
            xi_at_zero_user_unit,
            inp_target_user_unit,
            abs_tol=1e-3,
        ), (
            f"ξ(k=0) = {xi_at_zero_user_unit:.6f} {unit}, "
            f"expected inp TARGET = {inp_target_user_unit} {unit}"
        )

    def test_picker_returned_bare_restart(
        self, directory, inp_target_user_unit, unit,
    ):
        """Belt-and-suspenders: confirm the picker actually selected
        the rolling ``slowgrowth-1.restart`` file (not a snapshot)."""
        if not directory.exists():
            pytest.skip(f"SG fixture missing: {directory}")

        from md_analysis.engines.cp2k import CP2KParser

        picked = CP2KParser._find_restart(directory)
        assert picked.name == "slowgrowth-1.restart", (
            f"_find_restart returned {picked.name!r}; should be the "
            "bare rolling file."
        )
