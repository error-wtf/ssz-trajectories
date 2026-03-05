"""Orbit analysis: deflection, embedding statistics, bridge identity.

Authors: Carmen N. Wrede, Lino P. Casu
"""
import math

from .constants import c_light, LN_PHI
from .xi import D, s as s_fn
from .embedding import x_local, x_cumulative, count_jumps


def analyze_orbit(xi_fn, pts, r_s=1.0):
    """Analyse a non-radial null orbit.

    Parameters
    ----------
    xi_fn : callable
        Xi function.
    pts : list of (lam, r, phi)
        Orbit trajectory.
    r_s : float
        Schwarzschild radius.

    Returns
    -------
    dict with keys:
        xl_min, xl_max  - x_local range
        xl_jumps        - number of N-level crossings
        unique_N        - sorted list of visited N values
        dx_max          - max step-to-step x_local change
        xc_end          - cumulative embedding at orbit end
        xc_jumps        - jumps in cumulative embedding
        r_min, r_max    - radial range
        phi_total       - total azimuthal angle swept
        xl, xc          - full arrays
    """
    xl = [x_local(r, xi_fn, r_s) for _, r, _ in pts]
    xc = x_cumulative(pts, xi_fn)

    N_vals = [round(x) for x in xl]
    unique_N = sorted(set(N_vals))
    dx_max = max(
        (abs(xl[i + 1] - xl[i]) for i in range(len(xl) - 1)),
        default=0,
    )

    return {
        "xl_min": min(xl),
        "xl_max": max(xl),
        "xl_jumps": count_jumps(xl),
        "unique_N": unique_N,
        "dx_max": dx_max,
        "xc_end": xc[-1],
        "xc_jumps": count_jumps(xc),
        "r_min": min(r for _, r, _ in pts),
        "r_max": max(r for _, r, _ in pts),
        "phi_total": abs(pts[-1][2] - pts[0][2]),
        "xl": xl,
        "xc": xc,
    }


def deflection_angle(phi_total, b=None, r0=None, xi_fn=None):
    """Deflection angle delta = phi_total - pi.

    Parameters
    ----------
    phi_total : float
        Total azimuthal angle swept (radians).
    b, r0, xi_fn : ignored
        Kept for backward compatibility.

    Returns
    -------
    float
        Deflection angle in radians.
    """
    return phi_total - math.pi


def bridge_identity(xi_fn, r, r_s=1.0):
    """Verify bridge identity at radius r.

    Bridge: dl/dt = s(r) * dr/dt = s(r) * c * D(r)^2 = c * D(r)

    Returns
    -------
    dict with dr_dt, dl_dt, cD, relative_error
    """
    d = D(r, xi_fn, r_s)
    sv = s_fn(r, xi_fn, r_s)
    v_r = c_light * d ** 2          # dr/dt
    v_ell = sv * v_r                # dl/dt = s * dr/dt
    cD = c_light * d                # expected: c * D
    err = abs(v_ell - cD) / max(cD, 1e-30)
    return {
        "r": r,
        "D": d,
        "s": sv,
        "dr_dt": v_r,
        "dl_dt": v_ell,
        "cD": cD,
        "relative_error": err,
    }


def proper_radial_length(xi_fn, r_from, r_to, r_s=1.0, n=20000):
    """Proper radial distance: ell = integral s(r) dr."""
    from .integrator import trapz
    return trapz(lambda r: s_fn(r, xi_fn, r_s), r_from, r_to, n)


def tortoise_coordinate(xi_fn, r_from, r_to, r_s=1.0, n=20000):
    """Tortoise coordinate: r* = integral dr / D(r)^2."""
    from .integrator import trapz
    return trapz(lambda r: 1.0 / D(r, xi_fn, r_s) ** 2, r_from, r_to, n)
