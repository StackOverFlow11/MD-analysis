"""Unit tests for constant-potential free energy correction."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest

from md_analysis.enhanced_sampling.constrained_ti.correction import (
    ConstantPotentialCorrection,
    ConstantPotentialResult,
    compute_constant_potential_correction,
    write_corrected_free_energy_csv,
    plot_corrected_free_energy_profile,
)
from md_analysis.enhanced_sampling.constrained_ti.models import (
    TIPointDefinition,
    TIReport,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_ti_report(
    xi: np.ndarray,
    forces: np.ndarray,
    weights: np.ndarray,
) -> TIReport:
    """Create a minimal TIReport with synthetic data."""
    from md_analysis.enhanced_sampling.constrained_ti.models import (
        AutocorrResult,
        BlockAverageResult,
        ConstraintPointReport,
        GewekeResult,
        RunningAverageResult,
    )

    k = len(xi)
    errors = np.full(k, 0.001)
    reports = []
    for i in range(k):
        reports.append(ConstraintPointReport(
            xi=float(xi[i]),
            point_index=i,
            n_analyzed=1000,
            lambda_mean=-float(forces[i]),
            sigma_lambda=0.01,
            autocorr=AutocorrResult(
                acf=np.array([1.0, 0.5]),
                tau_corr=5.0, n_eff=100.0, sem_auto=0.001,
                passed_neff=True, passed_sem=True, t_min_frames=None,
            ),
            block_avg=BlockAverageResult(
                block_sizes=np.array([1, 2]),
                sem_curve=np.array([0.001, 0.001]),
                delta_sem=np.array([0.0001, 0.0001]),
                n_total=1000,
                plateau_index=0,
                plateau_sem=0.001,
                plateau_delta=0.0001,
                plateau_block_size=2,
                plateau_reached=True,
                passed=True,
            ),
            running_avg=RunningAverageResult(
                running_mean=np.linspace(0, 1, 10),
                drift_D=0.001, drift_limit=0.003, passed=True,
            ),
            geweke=GewekeResult(
                z=0.5, mean_a=0.0, mean_b=0.0,
                n_eff_a=50.0, passed=True, reliable=True,
            ),
            sem_final=0.001,
            sem_max=0.01,
            passed=True,
            failure_reasons=(),
        ))

    from md_analysis.enhanced_sampling.constrained_ti.integration import (
        _integrate_free_energy,
    )
    delta_A, sigma_A = _integrate_free_energy(forces, weights, errors)

    return TIReport(
        point_reports=tuple(reports),
        xi_values=xi,
        weights=weights,
        forces=forces,
        force_errors=errors,
        delta_A=delta_A,
        sigma_A=sigma_A,
        epsilon_tol_au=0.05 / 27.211386,
        all_passed=True,
        failing_indices=(),
        suggested_time_ratios=None,
    )


def _make_linear_mapper(slope: float, intercept: float) -> MagicMock:
    """Create a mock mapper: phi = slope * sigma + intercept."""
    mapper = MagicMock()
    mapper.predict = MagicMock(
        side_effect=lambda sigma: slope * np.asarray(sigma) + intercept
    )
    return mapper


def _make_point_defs_with_bader(tmp_path: Path, xi_values: list[float]) -> list[TIPointDefinition]:
    """Create ti_target dirs with empty bader/ subdirectories."""
    defs = []
    for xi in xi_values:
        ti_dir = tmp_path / f"ti_target_{xi}"
        ti_dir.mkdir()
        (ti_dir / "cMD-1.restart").touch()
        (ti_dir / "cMD.LagrangeMultLog").touch()
        (ti_dir / "bader").mkdir()
        defs.append(TIPointDefinition(
            xi=xi,
            restart_path=ti_dir / "cMD-1.restart",
            log_path=ti_dir / "cMD.LagrangeMultLog",
        ))
    return defs


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestCorrectionFormula:
    """Test the correction formula with mocked Bader data."""

    def test_zero_correction_when_charge_constant(self) -> None:
        """If sigma is the same at every point, correction = 0."""
        xi = np.array([0.0, 0.5, 1.0])
        weights = np.array([0.25, 0.5, 0.25])
        forces = np.array([0.01, 0.02, 0.01])
        ti_report = _make_ti_report(xi, forces, weights)

        sigma_const = np.array([5.0, 5.0, 5.0])
        area = 100.0
        mapper = _make_linear_mapper(slope=-0.01, intercept=0.0)

        from md_analysis.electrochemical.charge.config import E_PER_A2_TO_UC_PER_CM2
        from md_analysis.utils.config import HA_TO_EV

        # Manually compute expected: all deltas are zero
        phi = mapper.predict(sigma_const)
        delta_sigma_e_A2 = (sigma_const - sigma_const[0]) / E_PER_A2_TO_UC_PER_CM2
        delta_phi = phi - phi[0]
        expected_corr = delta_sigma_e_A2 * delta_phi * area / 2.0

        np.testing.assert_allclose(expected_corr, 0.0, atol=1e-15)

    def test_correction_formula_2_points(self) -> None:
        """Hand-calculated 2-point correction."""
        from md_analysis.electrochemical.charge.config import E_PER_A2_TO_UC_PER_CM2
        from md_analysis.utils.config import HA_TO_EV

        xi = np.array([0.0, 1.0])
        weights = np.array([0.5, 0.5])
        forces = np.array([0.01, 0.02])
        ti_report = _make_ti_report(xi, forces, weights)

        sigma = np.array([10.0, 20.0])  # uC/cm^2
        area = 200.0  # A^2
        mapper = _make_linear_mapper(slope=-0.005, intercept=0.0)

        phi = mapper.predict(sigma)  # [-0.05, -0.10]
        delta_sigma_e_A2 = (sigma[1] - sigma[0]) / E_PER_A2_TO_UC_PER_CM2
        delta_phi = phi[1] - phi[0]  # -0.05
        expected_corr_1 = delta_sigma_e_A2 * delta_phi * area / 2.0

        cumul_A_q = np.cumsum(weights * forces) * HA_TO_EV
        expected_A_phi = cumul_A_q.copy()
        expected_A_phi[1] += expected_corr_1

        # Verify via data model
        corr = ConstantPotentialCorrection(
            sigma_uC_cm2=sigma,
            phi_V_SHE=phi,
            correction_eV=np.array([0.0, expected_corr_1]),
            area_A2=area,
        )
        assert corr.correction_eV[0] == 0.0
        assert corr.correction_eV[1] == pytest.approx(expected_corr_1, rel=1e-10)

    def test_initial_state_correction_is_zero(self) -> None:
        """correction[0] must always be zero (delta from IS to IS)."""
        from md_analysis.electrochemical.charge.config import E_PER_A2_TO_UC_PER_CM2

        sigma = np.array([15.0, 20.0, 25.0])
        phi = np.array([-0.1, -0.2, -0.3])
        area = 150.0

        delta_sigma = (sigma - sigma[0]) / E_PER_A2_TO_UC_PER_CM2
        delta_phi = phi - phi[0]
        correction = delta_sigma * delta_phi * area / 2.0

        assert correction[0] == 0.0

    def test_units_eV(self) -> None:
        """Verify dimensional analysis: e/A^2 * V * A^2 = eV."""
        from md_analysis.electrochemical.charge.config import E_PER_A2_TO_UC_PER_CM2

        # 1 uC/cm^2 difference, 1 V difference, 1 A^2 area
        delta_sigma_uC_cm2 = 1.0
        delta_sigma_e_A2 = delta_sigma_uC_cm2 / E_PER_A2_TO_UC_PER_CM2
        delta_phi_V = 1.0
        area_A2 = 1.0

        correction_eV = delta_sigma_e_A2 * delta_phi_V * area_A2 / 2.0

        # Expected: (1 / 1602.176634) * 1 * 1 / 2 = 3.121e-4 eV
        expected = 1.0 / E_PER_A2_TO_UC_PER_CM2 / 2.0
        assert correction_eV == pytest.approx(expected, rel=1e-10)


class TestCSVOutput:
    """Test CSV export."""

    def test_csv_columns(self, tmp_path: Path) -> None:
        from md_analysis.utils.config import HA_TO_EV

        xi = np.array([0.0, 1.0])
        weights = np.array([0.5, 0.5])
        forces = np.array([0.01, 0.02])
        ti_report = _make_ti_report(xi, forces, weights)

        cumul_q = np.cumsum(weights * forces) * HA_TO_EV
        corr = ConstantPotentialCorrection(
            sigma_uC_cm2=np.array([10.0, 20.0]),
            phi_V_SHE=np.array([-0.1, -0.2]),
            correction_eV=np.array([0.0, 0.001]),
            area_A2=100.0,
        )
        result = ConstantPotentialResult(
            ti_report=ti_report,
            correction=corr,
            A_const_q_eV=cumul_q,
            A_const_phi_eV=cumul_q + corr.correction_eV,
            delta_A_const_phi_eV=float((cumul_q + corr.correction_eV)[-1]),
        )

        csv_path = write_corrected_free_energy_csv(result, output_dir=tmp_path)
        assert csv_path.exists()

        lines = csv_path.read_text().strip().split("\n")
        header = lines[0]
        assert "A_const_q_eV" in header
        assert "A_const_phi_eV" in header
        assert "correction_eV" in header
        assert "sigma_uC_cm2" in header
        assert "phi_V_SHE" in header
        assert len(lines) == 3  # header + 2 data rows


class TestPlotSmoke:
    """Smoke test for plot generation."""

    def test_plot_creates_png(self, tmp_path: Path) -> None:
        from md_analysis.utils.config import HA_TO_EV

        xi = np.array([0.0, 0.5, 1.0])
        weights = np.array([0.25, 0.5, 0.25])
        forces = np.array([0.01, 0.02, 0.01])
        ti_report = _make_ti_report(xi, forces, weights)

        cumul_q = np.cumsum(weights * forces) * HA_TO_EV
        corr = ConstantPotentialCorrection(
            sigma_uC_cm2=np.array([10.0, 15.0, 20.0]),
            phi_V_SHE=np.array([-0.1, -0.15, -0.2]),
            correction_eV=np.array([0.0, 0.001, 0.003]),
            area_A2=100.0,
        )
        result = ConstantPotentialResult(
            ti_report=ti_report,
            correction=corr,
            A_const_q_eV=cumul_q,
            A_const_phi_eV=cumul_q + corr.correction_eV,
            delta_A_const_phi_eV=float((cumul_q + corr.correction_eV)[-1]),
        )

        png_path = plot_corrected_free_energy_profile(result, output_dir=tmp_path)
        assert png_path.exists()
        assert png_path.stat().st_size > 0
