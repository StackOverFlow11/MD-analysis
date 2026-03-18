"""Tests for potential reference conversion."""

from __future__ import annotations

import numpy as np
import pytest

from md_analysis.electrochemical.calibration.CalibrationWorkflow import (
    convert_reference,
)
from md_analysis.electrochemical.calibration.config import (
    F_C_PER_MOL,
    LN10,
    R_J_PER_MOL_K,
)


class TestConvertReference:

    def test_same_reference_noop(self):
        phi = np.array([0.5, 1.0])
        result = convert_reference(phi, from_ref="SHE", to_ref="SHE")
        np.testing.assert_array_equal(result, phi)

    def test_she_to_rhe(self):
        """φ_RHE = φ_SHE + (RT/F)·ln(10)·pH."""
        T = 298.15
        pH = 7.0
        expected_shift = R_J_PER_MOL_K * T / F_C_PER_MOL * LN10 * pH
        phi_she = np.array([0.0])
        result = convert_reference(
            phi_she, from_ref="SHE", to_ref="RHE",
            temperature_K=T, pH=pH,
        )
        np.testing.assert_allclose(result, expected_shift, rtol=1e-10)

    def test_rhe_to_she(self):
        T = 298.15
        pH = 7.0
        shift = R_J_PER_MOL_K * T / F_C_PER_MOL * LN10 * pH
        phi_rhe = np.array([shift])
        result = convert_reference(
            phi_rhe, from_ref="RHE", to_ref="SHE",
            temperature_K=T, pH=pH,
        )
        np.testing.assert_allclose(result, 0.0, atol=1e-14)

    def test_she_to_pzc(self):
        phi_she = np.array([0.5])
        phi_pzc_val = 0.3
        result = convert_reference(
            phi_she, from_ref="SHE", to_ref="PZC", phi_pzc=phi_pzc_val,
        )
        np.testing.assert_allclose(result, 0.2)

    def test_pzc_to_she(self):
        phi_pzc_input = np.array([0.2])
        phi_pzc_val = 0.3
        result = convert_reference(
            phi_pzc_input, from_ref="PZC", to_ref="SHE", phi_pzc=phi_pzc_val,
        )
        np.testing.assert_allclose(result, 0.5)

    def test_roundtrip_she_rhe(self):
        phi = np.array([0.123, -0.456])
        T, pH = 310.0, 4.5
        via_rhe = convert_reference(phi, from_ref="SHE", to_ref="RHE",
                                    temperature_K=T, pH=pH)
        back = convert_reference(via_rhe, from_ref="RHE", to_ref="SHE",
                                 temperature_K=T, pH=pH)
        np.testing.assert_allclose(back, phi, atol=1e-14)

    def test_pzc_without_phi_pzc_raises(self):
        with pytest.raises(ValueError, match="phi_pzc"):
            convert_reference(np.array([0.5]), from_ref="SHE", to_ref="PZC")

    def test_unknown_reference_raises(self):
        with pytest.raises(ValueError, match="Unknown reference"):
            convert_reference(np.array([0.5]), from_ref="WEIRD", to_ref="SHE")

    def test_scalar_input(self):
        result = convert_reference(0.5, from_ref="SHE", to_ref="SHE")
        assert isinstance(result, np.ndarray)
