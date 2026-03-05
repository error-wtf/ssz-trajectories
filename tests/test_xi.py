"""Tests for Xi functions and derived metric quantities."""
import math
import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from ssz_trajectories import (
    PHI, Xi_strong, Xi_weak, Xi_hard, Xi_blend, D, s, dD_dr,
)


class TestXiStrong:
    def test_zero_radius(self):
        assert Xi_strong(0.0) == 0.0

    def test_negative_radius(self):
        assert Xi_strong(-1.0) == 0.0

    def test_at_rs(self):
        """Xi(r_s) = 1 - exp(-phi) ~ 0.802."""
        val = Xi_strong(1.0)
        assert abs(val - 0.802) < 0.001

    def test_monotone_increasing(self):
        """Xi_strong must increase with r."""
        prev = 0.0
        for r in [0.1, 0.5, 1, 2, 5, 10, 50, 100]:
            cur = Xi_strong(float(r))
            assert cur >= prev
            prev = cur

    def test_asymptotic_to_one(self):
        """Xi_strong -> 1 for large r."""
        assert Xi_strong(100.0) > 0.999

    def test_positive(self):
        for r in [0.01, 0.1, 1, 10, 100]:
            assert Xi_strong(float(r)) >= 0.0


class TestXiWeak:
    def test_at_large_r(self):
        """Xi_weak(200) = 1/(2*200) = 0.0025."""
        assert abs(Xi_weak(200.0) - 0.0025) < 1e-10

    def test_inversely_proportional(self):
        assert Xi_weak(100.0) == pytest.approx(2 * Xi_weak(200.0), rel=1e-10)

    def test_zero_radius(self):
        assert Xi_weak(0.0) == 0.0


class TestXiHard:
    def test_strong_regime(self):
        """r/r_s < 100 uses strong formula."""
        assert Xi_hard(50.0) == pytest.approx(Xi_strong(50.0))

    def test_weak_regime(self):
        """r/r_s >= 100 uses weak formula."""
        assert Xi_hard(100.0) == pytest.approx(Xi_weak(100.0))
        assert Xi_hard(200.0) == pytest.approx(Xi_weak(200.0))

    def test_discontinuity_at_100(self):
        """Hard switch has a jump at r/r_s = 100."""
        below = Xi_hard(99.99)
        above = Xi_hard(100.0)
        assert abs(below - above) > 0.5  # strong ~1.0, weak ~0.005


class TestXiBlend:
    def test_pure_strong(self):
        """Below blend zone: matches Xi_strong."""
        assert Xi_blend(50.0) == pytest.approx(Xi_strong(50.0))

    def test_pure_weak(self):
        """Above blend zone: matches Xi_weak."""
        assert Xi_blend(200.0) == pytest.approx(Xi_weak(200.0))

    def test_smooth_in_blend(self):
        """In blend zone [80,120]: intermediate value."""
        val = Xi_blend(100.0)
        assert Xi_weak(100.0) < val < Xi_strong(100.0)

    def test_c2_continuity(self):
        """No jump at blend boundaries."""
        lo_below = Xi_blend(79.99)
        lo_above = Xi_blend(80.01)
        assert abs(lo_below - lo_above) < 0.01

        hi_below = Xi_blend(119.99)
        hi_above = Xi_blend(120.01)
        assert abs(hi_below - hi_above) < 0.01

    def test_monotone_decreasing_in_blend(self):
        """Xi must decrease through the blend zone."""
        prev = Xi_blend(80.0)
        for r in range(85, 121, 5):
            cur = Xi_blend(float(r))
            assert cur <= prev + 1e-10
            prev = cur


class TestMetricD:
    def test_D_at_rs(self):
        """D(r_s) ~ 0.555."""
        d = D(1.0, Xi_strong)
        assert abs(d - 0.555) < 0.001

    def test_D_range(self):
        """0 < D <= 1 for all r > 0."""
        for r in [0.1, 1, 10, 100, 1000]:
            d = D(float(r), Xi_blend)
            assert 0 < d <= 1.0

    def test_s_inverse_of_D(self):
        """s(r) = 1/D(r)."""
        for r in [1, 10, 100]:
            assert s(float(r), Xi_blend) == pytest.approx(
                1.0 / D(float(r), Xi_blend)
            )


class TestDDerivative:
    def test_positive_in_strong(self):
        """D increases with r in strong field."""
        val = dD_dr(Xi_blend, 50.0)
        assert isinstance(val, float)

    def test_finite(self):
        for r in [1, 10, 100, 200]:
            val = dD_dr(Xi_blend, float(r))
            assert math.isfinite(val)
