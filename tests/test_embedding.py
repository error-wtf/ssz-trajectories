"""Tests for phi-Euler embedding functions."""
import math
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from ssz_trajectories import (
    PHI, x_local, N_level, epsilon_residual, count_jumps,
    Xi_strong, Xi_weak, Xi_blend,
)


class TestXLocal:
    def test_at_rs_strong(self):
        """x_local(r_s) with strong Xi ~ 1.44 (N=1)."""
        x = x_local(1.0, Xi_strong)
        assert 1.0 < x < 2.0

    def test_far_field_weak(self):
        """x_local -> 0 at large r in weak field."""
        x = x_local(1000.0, Xi_weak)
        assert abs(x) < 0.01

    def test_monotone_decreasing_outward_blend(self):
        """x_local decreases as r increases past blend zone."""
        prev = x_local(100.0, Xi_blend)
        for r in [150, 200, 500, 1000]:
            cur = x_local(float(r), Xi_blend)
            assert cur <= prev
            prev = cur


class TestNLevel:
    def test_N1_at_rs(self):
        assert N_level(1.0, Xi_strong) == 1

    def test_N0_far(self):
        assert N_level(1000.0, Xi_weak) == 0


class TestEpsilon:
    def test_range(self):
        """Epsilon must be in [-0.5, 0.5]."""
        for r in [0.5, 1, 5, 50, 200]:
            eps = epsilon_residual(float(r), Xi_blend)
            assert -0.5 <= eps <= 0.5


class TestCountJumps:
    def test_no_jumps(self):
        assert count_jumps([0.1, 0.2, 0.3, 0.4]) == 0

    def test_one_jump(self):
        assert count_jumps([0.1, 0.4, 0.6, 0.9]) == 1

    def test_two_jumps(self):
        assert count_jumps([0.1, 0.6, 1.2, 0.6, 0.1]) == 2

    def test_empty(self):
        assert count_jumps([]) == 0

    def test_single(self):
        assert count_jumps([1.0]) == 0
