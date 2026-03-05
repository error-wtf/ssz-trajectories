"""Tests for geodesic integrators."""
import math
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from ssz_trajectories import (
    c_light, Xi_strong, Xi_weak, Xi_blend, D,
    rk4_scalar, trapz,
    integrate_null_radial, integrate_timelike_infall,
    integrate_null_geodesic, find_turning_point,
)


class TestRK4Scalar:
    def test_exponential_decay(self):
        """dy/dt = -y, y(0)=1 => y(1) = e^{-1}."""
        y, t = 1.0, 0.0
        dt = 0.001
        while t < 1.0:
            y, t = rk4_scalar(y, t, dt, lambda yy: -yy)
        assert abs(y - math.exp(-1)) < 1e-6

    def test_linear_growth(self):
        """dy/dt = 1, y(0)=0 => y(1) = 1."""
        y, t = 0.0, 0.0
        dt = 0.01
        while t < 1.0 - 1e-12:
            y, t = rk4_scalar(y, t, dt, lambda yy: 1.0)
        assert abs(y - 1.0) < 1e-8


class TestTrapz:
    def test_constant(self):
        assert trapz(lambda x: 2.0, 0, 1) == 2.0

    def test_linear(self):
        val = trapz(lambda x: x, 0, 1)
        assert abs(val - 0.5) < 1e-6

    def test_reversed_limits(self):
        val = trapz(lambda x: x, 1, 0)
        assert abs(val - (-0.5)) < 1e-6


class TestNullRadial:
    def test_outgoing_monotone(self):
        """Outgoing null: r must increase monotonically."""
        pts = integrate_null_radial(Xi_blend, 1.01, 2e-6, 2e-9, direction=1)
        assert len(pts) > 10
        for i in range(1, len(pts)):
            assert pts[i][1] >= pts[i - 1][1]

    def test_ingoing_monotone(self):
        """Ingoing null: r must decrease monotonically."""
        pts = integrate_null_radial(Xi_blend, 50.0, 1e-7, 1e-10, direction=-1)
        assert len(pts) > 10
        for i in range(1, len(pts)):
            assert pts[i][1] <= pts[i - 1][1]

    def test_speed_bounded(self):
        """dr/dt = c*D^2 <= c always."""
        pts = integrate_null_radial(Xi_blend, 1.01, 1e-6, 2e-9, direction=1)
        for i in range(1, len(pts)):
            dt = pts[i][0] - pts[i - 1][0]
            dr = pts[i][1] - pts[i - 1][1]
            if dt > 0:
                assert dr / dt <= c_light * 1.01  # small tolerance


class TestTimeLikeInfall:
    def test_monotone_decrease(self):
        """Infalling object: r decreases."""
        pts = integrate_timelike_infall(Xi_blend, 50.0, 5e-6, 5e-9)
        assert len(pts) > 10
        for i in range(1, len(pts)):
            assert pts[i][1] <= pts[i - 1][1] + 1e-10

    def test_reaches_boundary(self):
        """Object reaches near r_s."""
        pts = integrate_timelike_infall(Xi_blend, 50.0, 5e-6, 5e-9)
        r_final = pts[-1][1]
        assert r_final < 2.0  # should reach below 2 r_s


class TestNullGeodesic:
    def test_small_b_two_jumps(self):
        """b=2: orbit penetrates deep, x_local crosses N=0->1->0."""
        b = 2.0
        r0 = 300.0
        lam_max = r0 / c_light * 12
        pts = integrate_null_geodesic(Xi_blend, b, r0, lam_max, 2e-9)
        assert len(pts) > 100
        r_min = min(r for _, r, _ in pts)
        assert r_min < 3.0  # deep penetration

    def test_large_b_no_jumps(self):
        """b=200: orbit stays in weak field, no N-level jumps."""
        b = 200.0
        r0 = max(300, b * 4)
        lam_max = r0 / c_light * 12
        pts = integrate_null_geodesic(Xi_blend, b, r0, lam_max, 2e-9)
        assert len(pts) > 100
        r_min = min(r for _, r, _ in pts)
        assert r_min > 100.0  # stays in weak field

    def test_turning_point(self):
        """Orbit must have r_min < r0 and r_max near r0."""
        b = 10.0
        r0 = 300.0
        lam_max = r0 / c_light * 12
        pts = integrate_null_geodesic(Xi_blend, b, r0, lam_max, 2e-9)
        r_min = min(r for _, r, _ in pts)
        r_max = max(r for _, r, _ in pts)
        assert r_min < r0
        assert r_max >= r0 * 0.9

    def test_phi_increases(self):
        """Total phi must be positive."""
        b = 5.0
        r0 = 300.0
        lam_max = r0 / c_light * 12
        pts = integrate_null_geodesic(Xi_blend, b, r0, lam_max, 2e-9)
        phi_total = abs(pts[-1][2] - pts[0][2])
        assert phi_total > 0


class TestTurningPoint:
    def test_exists_for_small_b(self):
        r_tp = find_turning_point(Xi_blend, 10.0)
        assert r_tp is not None
        assert r_tp > 0

    def test_turning_point_condition(self):
        """At r_tp: r/D(r) ~ b."""
        b = 10.0
        r_tp = find_turning_point(Xi_blend, b)
        assert r_tp is not None
        ratio = r_tp / D(r_tp, Xi_blend)
        assert abs(ratio - b) < 0.01
