"""Tests for integration: trapezoid weights, SEM targets, error propagation."""

import numpy as np
import pytest

from md_analysis.enhanced_sampling.constrained_ti.integration import (
    _compute_sem_targets,
    _integrate_free_energy,
    _suggest_time_allocation,
    compute_trapezoid_weights,
    sem_chi2_upper_bound,
)


class TestTrapezoidWeights:
    def test_uniform_2_points(self):
        xi = np.array([0.0, 1.0])
        w = compute_trapezoid_weights(xi)
        assert w == pytest.approx([0.5, 0.5])

    def test_uniform_3_points(self):
        xi = np.array([0.0, 1.0, 2.0])
        w = compute_trapezoid_weights(xi)
        assert w == pytest.approx([0.5, 1.0, 0.5])

    def test_sum_equals_range(self):
        xi = np.array([1.0, 2.5, 4.0, 6.0])
        w = compute_trapezoid_weights(xi)
        assert np.sum(w) == pytest.approx(xi[-1] - xi[0])

    def test_k_less_than_2_raises(self):
        with pytest.raises(ValueError, match="at least 2"):
            compute_trapezoid_weights(np.array([1.0]))

    def test_non_monotonic_raises(self):
        with pytest.raises(ValueError, match="monotonic"):
            compute_trapezoid_weights(np.array([1.0, 3.0, 2.0]))

    def test_decreasing_ok(self):
        xi = np.array([3.0, 2.0, 1.0])
        w = compute_trapezoid_weights(xi)
        assert np.sum(np.abs(w)) > 0


class TestSEMTargets:
    def test_uniform_weights(self):
        w = np.array([0.5, 1.0, 0.5])
        targets = _compute_sem_targets(w, epsilon_tol_au=0.01)
        # Endpoints have smaller |w|, so larger SEM_max
        assert targets[0] > targets[1]
        assert targets[2] > targets[1]


class TestIntegrateFreeEnergy:
    def test_2_point_exact(self):
        forces = np.array([1.0, 3.0])
        weights = np.array([0.5, 0.5])
        sems = np.array([0.1, 0.2])
        delta_A, sigma_A = _integrate_free_energy(forces, weights, sems)
        assert delta_A == pytest.approx(0.5 * 1.0 + 0.5 * 3.0)
        assert sigma_A == pytest.approx(np.sqrt(0.5**2 * 0.1**2 + 0.5**2 * 0.2**2))


class TestTimeAllocation:
    def test_normalized(self):
        w = np.array([0.5, 1.0, 0.5])
        s = np.array([1.0, 2.0, 1.0])
        t = np.array([10.0, 5.0, 10.0])
        ratios = _suggest_time_allocation(w, s, t)
        assert np.sum(ratios) == pytest.approx(1.0)
        assert len(ratios) == 3


class TestSEMChi2UpperBound:
    """Chi-square one-sided 95% upper bound of the block SEM."""

    def test_nu15_factor_magnitude(self):
        # chi2_0.05(15) = 7.261 → factor = sqrt(15/7.261) ≈ 1.44
        sem_report, factor = sem_chi2_upper_bound(1.0, 16)
        assert factor == pytest.approx(1.437, abs=1e-3)
        assert sem_report == pytest.approx(1.437, abs=1e-3)

    def test_factor_decreases_towards_one(self):
        _, f5 = sem_chi2_upper_bound(1.0, 6)
        _, f15 = sem_chi2_upper_bound(1.0, 16)
        _, f100 = sem_chi2_upper_bound(1.0, 101)
        assert f5 > f15 > f100 > 1.0

    def test_sem_scaled_by_factor(self):
        sem_report, factor = sem_chi2_upper_bound(0.01, 16)
        assert sem_report == pytest.approx(0.01 * factor, rel=1e-12)

    def test_fewer_than_2_blocks_raises(self):
        with pytest.raises(ValueError, match="n_blocks"):
            sem_chi2_upper_bound(1.0, 1)
