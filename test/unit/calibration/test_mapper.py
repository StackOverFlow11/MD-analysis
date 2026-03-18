"""Tests for charge-potential mapper ABC and implementations."""

from __future__ import annotations

import numpy as np
import pytest

from md_analysis.electrochemical.calibration._mapper import (
    FittingInfo,
    LinearMapper,
    PolynomialMapper,
    SplineMapper,
    create_mapper,
    mapper_from_dict,
)


# ---------------------------------------------------------------------------
# LinearMapper
# ---------------------------------------------------------------------------

class TestLinearMapper:

    def test_perfect_fit(self):
        """Two colinear points → R²=1."""
        m = LinearMapper()
        sigma = np.array([0.0, 10.0])
        phi = np.array([0.5, 1.5])
        info = m.fit(sigma, phi)
        assert info.r_squared == pytest.approx(1.0)
        assert info.method == "linear"

    def test_slope_intercept(self):
        m = LinearMapper()
        sigma = np.array([0.0, 10.0])
        phi = np.array([0.5, 1.5])
        m.fit(sigma, phi)
        assert m._slope == pytest.approx(0.1)
        assert m._intercept == pytest.approx(0.5)

    def test_predict(self):
        m = LinearMapper()
        m.fit(np.array([0.0, 10.0]), np.array([0.5, 1.5]))
        result = m.predict(5.0)
        np.testing.assert_allclose(result, 1.0)

    def test_predict_array(self):
        m = LinearMapper()
        m.fit(np.array([0.0, 10.0]), np.array([0.5, 1.5]))
        result = m.predict(np.array([0.0, 5.0, 10.0]))
        np.testing.assert_allclose(result, [0.5, 1.0, 1.5])

    def test_predict_before_fit_raises(self):
        m = LinearMapper()
        with pytest.raises(RuntimeError, match="not been fitted"):
            m.predict(1.0)

    def test_serialization_roundtrip(self):
        m = LinearMapper()
        m.fit(np.array([0.0, 10.0]), np.array([0.5, 1.5]))
        d = m.to_dict()
        m2 = LinearMapper.from_dict(d)
        np.testing.assert_allclose(m2.predict(5.0), m.predict(5.0))

    def test_equation_str(self):
        m = LinearMapper()
        info = m.fit(np.array([0.0, 10.0]), np.array([0.5, 1.5]))
        assert "φ" in info.equation_str
        assert "σ" in info.equation_str


# ---------------------------------------------------------------------------
# PolynomialMapper
# ---------------------------------------------------------------------------

class TestPolynomialMapper:

    def test_quadratic_fit(self):
        """3 points on y = x² → R²=1."""
        sigma = np.array([-1.0, 0.0, 1.0])
        phi = np.array([1.0, 0.0, 1.0])
        m = PolynomialMapper(degree=2)
        info = m.fit(sigma, phi)
        assert info.r_squared == pytest.approx(1.0, abs=1e-10)

    def test_predict(self):
        sigma = np.array([-1.0, 0.0, 1.0])
        phi = np.array([1.0, 0.0, 1.0])
        m = PolynomialMapper(degree=2)
        m.fit(sigma, phi)
        np.testing.assert_allclose(m.predict(2.0), 4.0, atol=1e-10)

    def test_predict_before_fit_raises(self):
        m = PolynomialMapper(degree=2)
        with pytest.raises(RuntimeError):
            m.predict(1.0)

    def test_serialization_roundtrip(self):
        sigma = np.array([-1.0, 0.0, 1.0, 2.0])
        phi = np.array([1.0, 0.0, 1.0, 4.0])
        m = PolynomialMapper(degree=2)
        m.fit(sigma, phi)
        d = m.to_dict()
        m2 = PolynomialMapper.from_dict(d)
        np.testing.assert_allclose(m2.predict(1.5), m.predict(1.5))


# ---------------------------------------------------------------------------
# SplineMapper
# ---------------------------------------------------------------------------

class TestSplineMapper:

    @pytest.fixture(autouse=True)
    def _check_scipy(self):
        pytest.importorskip("scipy")

    def test_interpolates_through_points(self):
        sigma = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
        phi = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        m = SplineMapper()
        m.fit(sigma, phi)
        np.testing.assert_allclose(m.predict(sigma), phi, atol=1e-12)

    def test_predict_before_fit_raises(self):
        m = SplineMapper()
        with pytest.raises(RuntimeError):
            m.predict(1.0)

    def test_serialization_roundtrip(self):
        sigma = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
        phi = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        m = SplineMapper()
        m.fit(sigma, phi)
        d = m.to_dict()
        m2 = SplineMapper.from_dict(d)
        test_s = np.array([-1.5, 0.5, 1.5])
        np.testing.assert_allclose(m2.predict(test_s), m.predict(test_s))


# ---------------------------------------------------------------------------
# Factory + mapper_from_dict
# ---------------------------------------------------------------------------

class TestFactory:

    def test_create_linear(self):
        m = create_mapper("linear")
        assert isinstance(m, LinearMapper)

    def test_create_polynomial(self):
        m = create_mapper("polynomial", degree=3)
        assert isinstance(m, PolynomialMapper)
        assert m._degree == 3

    def test_create_unknown_raises(self):
        with pytest.raises(ValueError, match="Unknown fitting method"):
            create_mapper("magic")

    def test_mapper_from_dict(self):
        m = LinearMapper()
        m.fit(np.array([0.0, 1.0]), np.array([0.0, 1.0]))
        d = m.to_dict()
        m2 = mapper_from_dict(d)
        assert isinstance(m2, LinearMapper)
        np.testing.assert_allclose(m2.predict(0.5), 0.5)


# ---------------------------------------------------------------------------
# FittingInfo
# ---------------------------------------------------------------------------

class TestFittingInfo:

    def test_defaults(self):
        fi = FittingInfo()
        assert fi.method == ""
        assert fi.r_squared == 0.0
        assert fi.rmse == 0.0

    def test_populated(self):
        m = LinearMapper()
        info = m.fit(np.array([0.0, 10.0]), np.array([0.5, 1.5]))
        assert info.r_squared == pytest.approx(1.0)
        assert info.rmse == pytest.approx(0.0, abs=1e-14)
        assert len(info.residuals) == 2
