"""SSZ Xi(r) segment-density functions and derived metric quantities.

Authors: Carmen N. Wrede, Lino P. Casu
"""
import math
from .constants import PHI, BLEND_LO, BLEND_HI, RS_DEFAULT


def Xi_strong(r, r_s=RS_DEFAULT):
    """Strong-field: Xi = 1 - exp(-phi * r / r_s)."""
    if r <= 0:
        return 0.0
    return 1.0 - math.exp(-PHI * r / r_s)


def Xi_weak(r, r_s=RS_DEFAULT):
    """Weak-field: Xi = r_s / (2r)."""
    if r <= 0:
        return 0.0
    return r_s / (2.0 * r)


def Xi_hard(r, r_s=RS_DEFAULT):
    """Hard switch at r/r_s = 100."""
    if r <= 0:
        return 0.0
    if r / r_s >= 100:
        return Xi_weak(r, r_s)
    return Xi_strong(r, r_s)


def _hermite(t):
    """C2 Hermite blend: h(t) = 3t^2 - 2t^3."""
    t = max(0.0, min(1.0, t))
    return t * t * (3.0 - 2.0 * t)


def Xi_blend(r, r_s=RS_DEFAULT, r_lo=BLEND_LO, r_hi=BLEND_HI):
    """C2 Hermite-blended Xi: smooth transition in [r_lo, r_hi] r_s."""
    if r <= 0:
        return 0.0
    ratio = r / r_s
    if ratio < r_lo:
        return Xi_strong(r, r_s)
    if ratio > r_hi:
        return Xi_weak(r, r_s)
    t = (ratio - r_lo) / (r_hi - r_lo)
    w = _hermite(t)
    return (1.0 - w) * Xi_strong(r, r_s) + w * Xi_weak(r, r_s)


# -- Derived metric functions --


def D(r, xi_fn=None, r_s=RS_DEFAULT):
    """Dilation factor D = 1/(1+Xi). Pass xi_fn or uses Xi_blend."""
    if xi_fn is None:
        xi_val = Xi_blend(r, r_s)
    elif callable(xi_fn):
        xi_val = xi_fn(r)
    else:
        xi_val = xi_fn
    return 1.0 / (1.0 + xi_val)


def s(r, xi_fn=None, r_s=RS_DEFAULT):
    """Segment factor s = 1/D = 1+Xi."""
    return 1.0 / D(r, xi_fn, r_s)


def dD_dr(xi_fn, r, r_s=RS_DEFAULT):
    """Numerical derivative dD/dr."""
    eps = max(abs(r) * 1e-7, 1e-10)
    return (D(r + eps, xi_fn, r_s) - D(r - eps, xi_fn, r_s)) / (2.0 * eps)
