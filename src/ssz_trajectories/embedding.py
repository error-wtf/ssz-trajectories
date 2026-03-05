"""Phi-Euler embedding: x_local, x_cumulative, N-level counting.

The local embedding maps the segment factor s(r) to a continuous
real-valued 'level':
    x_local(r) = ln(s(r)) / ln(phi)

The integer part N = round(x_local) counts how many full phi-levels
the local geometry occupies.  Jumps in N along a trajectory indicate
transitions between phi-level domains.

Authors: Carmen N. Wrede, Lino P. Casu
"""
import math

from .constants import LN_PHI
from .xi import D, s as s_fn


def x_local(r, xi_fn=None, r_s=1.0):
    """Local embedding: x = ln(s(r)) / ln(phi)."""
    sv = s_fn(r, xi_fn, r_s)
    if sv <= 1e-30:
        return 0.0
    return math.log(sv) / LN_PHI


def N_level(r, xi_fn=None, r_s=1.0):
    """Integer phi-level: N = round(x_local)."""
    return round(x_local(r, xi_fn, r_s))


def epsilon_residual(r, xi_fn=None, r_s=1.0):
    """Fractional part: eps = x - N."""
    x = x_local(r, xi_fn, r_s)
    return x - round(x)


def x_cumulative(pts, xi_fn):
    """Cumulative embedding along a trajectory.

    Parameters
    ----------
    pts : list of (lam, r, phi)
        Orbit points from the integrator.
    xi_fn : callable
        Xi function.

    Returns
    -------
    list of float
        Cumulative x values.
    """
    xc = [0.0]
    for i in range(1, len(pts)):
        dr = abs(pts[i][1] - pts[i - 1][1])
        xi_mid = 0.5 * (xi_fn(pts[i][1]) + xi_fn(pts[i - 1][1]))
        xc.append(xc[-1] + xi_mid * dr / LN_PHI)
    return xc


def count_jumps(arr):
    """Count transitions where round(arr[i]) changes."""
    return sum(
        1 for i in range(1, len(arr))
        if round(arr[i]) != round(arr[i - 1])
    )
