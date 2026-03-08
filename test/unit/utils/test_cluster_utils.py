"""Unit tests for ClusterUtils."""

from __future__ import annotations

import numpy as np
import pytest

from md_analysis.utils.StructureParser.ClusterUtils import (
    cluster_1d_periodic,
    find_largest_gap_periodic,
    gap_midpoint_periodic,
)


# ---------------------------------------------------------------------------
# cluster_1d_periodic
# ---------------------------------------------------------------------------

class TestCluster1dPeriodic:

    def test_basic_two_clusters(self) -> None:
        values = np.array([1.0, 1.1, 1.2, 5.0, 5.1, 5.2])
        result = cluster_1d_periodic(values, period=10.0, tol=0.5)
        assert len(result) == 2
        # first cluster near 1.1, second near 5.1
        assert abs(result[0][0] - 1.1) < 0.2
        assert abs(result[1][0] - 5.1) < 0.2
        # check member counts
        assert result[0][1].size == 3
        assert result[1][1].size == 3

    def test_pbc_wrap_merge(self) -> None:
        """Two clusters that straddle the periodic boundary should merge."""
        values = np.array([0.1, 0.2, 9.8, 9.9])
        result = cluster_1d_periodic(values, period=10.0, tol=0.5)
        # Should merge into one cluster
        assert len(result) == 1
        assert result[0][1].size == 4

    def test_single_cluster(self) -> None:
        values = np.array([3.0, 3.1, 3.2])
        result = cluster_1d_periodic(values, period=10.0, tol=0.5)
        assert len(result) == 1
        assert result[0][1].size == 3

    def test_single_value(self) -> None:
        values = np.array([5.0])
        result = cluster_1d_periodic(values, period=10.0, tol=0.5)
        assert len(result) == 1
        assert result[0][1].size == 1

    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError, match="non-empty"):
            cluster_1d_periodic(np.array([]), period=10.0, tol=0.5)

    def test_negative_period_raises(self) -> None:
        with pytest.raises(ValueError, match="period"):
            cluster_1d_periodic(np.array([1.0]), period=-1.0, tol=0.5)

    def test_negative_tol_raises(self) -> None:
        with pytest.raises(ValueError, match="tol"):
            cluster_1d_periodic(np.array([1.0]), period=10.0, tol=-0.1)

    def test_member_indices_map_to_original(self) -> None:
        """Member indices should index back into the original values array."""
        values = np.array([8.0, 2.0, 8.1, 2.1])
        result = cluster_1d_periodic(values, period=10.0, tol=0.5)
        assert len(result) == 2
        all_indices = np.sort(np.concatenate([r[1] for r in result]))
        np.testing.assert_array_equal(all_indices, [0, 1, 2, 3])

    def test_values_outside_period_folded(self) -> None:
        """Values outside [0, period) should be folded."""
        values = np.array([11.0, 11.1, -0.1])  # 11.0 % 10 = 1.0, etc.
        result = cluster_1d_periodic(values, period=10.0, tol=0.5)
        assert len(result) == 2  # cluster near 1.0 and cluster near 9.9


# ---------------------------------------------------------------------------
# find_largest_gap_periodic
# ---------------------------------------------------------------------------

class TestFindLargestGapPeriodic:

    def test_basic(self) -> None:
        centers = np.array([1.0, 3.0, 9.0])
        low_idx, high_idx, gap = find_largest_gap_periodic(centers, period=10.0)
        # Gap from 3.0 to 9.0 is 6.0
        assert low_idx == 1
        assert high_idx == 2
        assert abs(gap - 6.0) < 1e-10

    def test_wrap_around_gap(self) -> None:
        """Largest gap wraps around the boundary."""
        centers = np.array([1.0, 2.0, 3.0])
        low_idx, high_idx, gap = find_largest_gap_periodic(centers, period=10.0)
        # Wrap gap: 1.0 + (10.0 - 3.0) = 8.0
        assert low_idx == 2
        assert high_idx == 0
        assert abs(gap - 8.0) < 1e-10

    def test_two_centers(self) -> None:
        centers = np.array([2.0, 8.0])
        low_idx, high_idx, gap = find_largest_gap_periodic(centers, period=10.0)
        assert gap == 6.0

    def test_too_few_raises(self) -> None:
        with pytest.raises(ValueError, match="Need >= 2"):
            find_largest_gap_periodic(np.array([1.0]), period=10.0)


# ---------------------------------------------------------------------------
# gap_midpoint_periodic
# ---------------------------------------------------------------------------

class TestGapMidpointPeriodic:

    def test_basic(self) -> None:
        mid = gap_midpoint_periodic(2.0, 8.0, period=10.0)
        assert abs(mid - 5.0) < 1e-10

    def test_wrap_around(self) -> None:
        mid = gap_midpoint_periodic(8.0, 2.0, period=10.0)
        # Gap is 4.0 (8→10→2), midpoint at 10.0 (=0.0 mod 10)
        assert abs(mid - 0.0) < 1e-10 or abs(mid - 10.0) < 1e-10

    def test_negative_period_raises(self) -> None:
        with pytest.raises(ValueError, match="period"):
            gap_midpoint_periodic(1.0, 2.0, period=-1.0)
