"""Tests for orbit analysis and bridge identity."""
import math
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from ssz_trajectories import (
    c_light, Xi_strong, Xi_weak, Xi_blend, D,
    integrate_null_geodesic, integrate_null_radial,
    analyze_orbit, deflection_angle, bridge_identity,
    proper_radial_length, tortoise_coordinate,
)


def _orbit(b, r0=None):
    if r0 is None:
        r0 = max(500.0, b * 10)
    lam_max = r0 / c_light * 12
    dlam = 2e-9 if b < 20 else 5e-8
    pts = integrate_null_geodesic(Xi_blend, b, r0, lam_max, dlam)
    return pts


class TestAnalyzeOrbit:
    def test_small_b_jumps(self):
        """b=3: must have 2 N-level jumps (entry + exit)."""
        res = analyze_orbit(Xi_blend, _orbit(3.0))
        assert res["xl_jumps"] == 2

    def test_b50_jumps(self):
        """b=50: entry + exit must both be counted."""
        res = analyze_orbit(Xi_blend, _orbit(50.0))
        assert res["xl_jumps"] == 2

    def test_b80_jumps(self):
        """b=80: entry + exit must both be counted."""
        res = analyze_orbit(Xi_blend, _orbit(80.0))
        assert res["xl_jumps"] == 2

    def test_large_b_no_jumps(self):
        """b=200: weak field only, no jumps."""
        res = analyze_orbit(Xi_blend, _orbit(200.0))
        assert res["xl_jumps"] == 0

    def test_phi_total_positive(self):
        res = analyze_orbit(Xi_blend, _orbit(5.0))
        assert res["phi_total"] > 0

    def test_r_range(self):
        res = analyze_orbit(Xi_blend, _orbit(10.0))
        assert res["r_min"] < 20
        assert res["r_max"] > 200


class TestDeflection:
    def test_is_dphi_minus_pi(self):
        """deflection_angle must return phi_total - pi."""
        assert abs(deflection_angle(4.0) - (4.0 - math.pi)) < 1e-15

    def test_positive_deep_orbit(self):
        """Deep orbits (b=5) must have delta > 0."""
        pts = _orbit(5.0, r0=2000.0)
        phi_total = abs(pts[-1][2] - pts[0][2])
        assert deflection_angle(phi_total) > 0

    def test_increases_below_barrier(self):
        """In SSZ strong field, deflection increases with b
        (more strong-field path traversed) up to the blend-zone barrier."""
        defls = []
        for b in [5.0, 50.0, 80.0]:
            pts = _orbit(b, r0=5000.0)
            phi_total = abs(pts[-1][2] - pts[0][2])
            defls.append(deflection_angle(phi_total))
        assert defls[0] < defls[1] < defls[2]


class TestBridgeIdentity:
    def test_exact_at_all_radii(self):
        """dl/dt = c*D must hold exactly (up to floating-point)."""
        for r in [1.01, 2.0, 10.0, 50.0, 100.0, 1000.0]:
            res = bridge_identity(Xi_blend, r)
            assert res["relative_error"] < 1e-12

    def test_components(self):
        """dr/dt = c*D^2, dl/dt = s*dr/dt, and dl/dt = c*D."""
        res = bridge_identity(Xi_blend, 1.0)
        assert res["dr_dt"] == c_light * res["D"] ** 2
        assert abs(res["dl_dt"] - res["cD"]) < 1e-6


class TestProperLength:
    def test_finite_to_boundary(self):
        """Proper length from r_s to 2*r_s must be finite."""
        ell = proper_radial_length(Xi_blend, 1.001, 2.0)
        assert 0 < ell < 1e10
        assert math.isfinite(ell)

    def test_positive(self):
        ell = proper_radial_length(Xi_blend, 1.0, 10.0)
        assert ell > 0


class TestTortoise:
    def test_finite_no_horizon(self):
        """r* from r_s to 2*r_s must be finite (SSZ has no horizon)."""
        rstar = tortoise_coordinate(Xi_blend, 1.001, 2.0)
        assert math.isfinite(rstar)
        assert rstar > 0
