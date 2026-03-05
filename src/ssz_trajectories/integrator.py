"""Numerical integrators for SSZ geodesics.

Provides:
- Generic RK4 scalar integrator
- Trapezoidal quadrature
- 2nd-order geodesic integrator for null geodesics (radial + non-radial)

Authors: Carmen N. Wrede, Lino P. Casu
"""
import math

from .constants import c_light, RS_DEFAULT
from .xi import D, dD_dr


def rk4_scalar(y, t, dt, f):
    """Single-step RK4 for scalar ODE dy/dt = f(y).

    Returns (y_new, t_new).
    """
    k1 = f(y)
    k2 = f(y + 0.5 * dt * k1)
    k3 = f(y + 0.5 * dt * k2)
    k4 = f(y + dt * k3)
    return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4), t + dt


def trapz(f, a, b, n=20000):
    """Trapezoidal integration of f from a to b."""
    if a == b:
        return 0.0
    sign = 1.0
    if b < a:
        a, b = b, a
        sign = -1.0
    h = (b - a) / n
    acc = 0.5 * (f(a) + f(b))
    for i in range(1, n):
        acc += f(a + i * h)
    return sign * acc * h


def integrate_null_radial(xi_fn, r0, tmax, dt, direction=1):
    """Integrate radial null geodesic: dr/dt = direction * c * D(r)^2.

    Parameters
    ----------
    xi_fn : callable
        Xi function (takes r, returns Xi value).
    r0 : float
        Starting radius.
    tmax : float
        Maximum coordinate time.
    dt : float
        Time step.
    direction : {+1, -1}
        +1 for outgoing, -1 for ingoing.

    Returns
    -------
    list of (t, r)
    """
    def drdt(r):
        d = D(r, xi_fn)
        return direction * c_light * d ** 2

    r, t = float(r0), 0.0
    pts = [(t, r)]
    while t < tmax and r > 0.01 * RS_DEFAULT:
        r, t = rk4_scalar(r, t, dt, drdt)
        pts.append((t, r))
    return pts


def integrate_timelike_infall(xi_fn, r0, tau_max, dtau):
    """Integrate radial timelike freefall: dr/dtau = -c*sqrt(1 - D^2).

    Parameters
    ----------
    xi_fn : callable
        Xi function.
    r0 : float
        Starting radius.
    tau_max : float
        Maximum proper time.
    dtau : float
        Proper-time step.

    Returns
    -------
    list of (tau, r)
    """
    def f(r):
        d = D(r, xi_fn)
        val = max(1.0 - d ** 2, 0.0)
        return -c_light * math.sqrt(val)

    r, tau = float(r0), 0.0
    pts = [(tau, r)]
    while tau < tau_max and r > 0.01 * RS_DEFAULT:
        r, tau = rk4_scalar(r, tau, dtau, f)
        pts.append((tau, r))
    return pts


def integrate_null_geodesic(xi_fn, b, r0, lam_max, dlam, r_s=RS_DEFAULT):
    """Integrate non-radial null geodesic using 2nd-order geodesic equation.

    Uses state (r, phi, vr) where vr = dr/dlam.  The radial acceleration
    is derived from Christoffel symbols of the SSZ metric:

        d2r/dl2 = -D'c^2/D + (D'/D)*vr^2 + c^2*b^2*D^2/r^3

    No sign-switching is needed: vr reverses naturally at the turning point.

    Parameters
    ----------
    xi_fn : callable
        Xi function (takes r, returns Xi).
    b : float
        Impact parameter (same units as r).
    r0 : float
        Starting radius.
    lam_max : float
        Maximum affine parameter.
    dlam : float
        Step size in affine parameter.
    r_s : float
        Schwarzschild radius (for boundary check).

    Returns
    -------
    list of (lam, r, phi)
    """
    D0 = D(r0, xi_fn)
    arg0 = max(1.0 / D0 ** 2 - b ** 2 / r0 ** 2, 0.0)
    vr = -D0 * c_light * math.sqrt(arg0)  # inward
    r, phi, lam = float(r0), 0.0, 0.0
    pts = [(lam, r, phi)]

    def derivs(r_, vr_, phi_):
        if r_ < 0.3 * r_s:
            return (0.0, 0.0, 0.0)
        D_ = D(r_, xi_fn)
        Dp = dD_dr(xi_fn, r_)
        dphi = c_light * b / r_ ** 2
        acc = (
            -Dp * c_light ** 2 / D_
            + (Dp / D_) * vr_ ** 2
            + c_light ** 2 * b ** 2 * D_ ** 2 / r_ ** 3
        )
        return (vr_, acc, dphi)

    while lam < lam_max and 0.3 * r_s < r < 1.1 * r0:
        k1 = derivs(r, vr, phi)
        r2 = r + 0.5 * dlam * k1[0]
        v2 = vr + 0.5 * dlam * k1[1]
        p2 = phi + 0.5 * dlam * k1[2]
        k2 = derivs(r2, v2, p2)
        r3 = r + 0.5 * dlam * k2[0]
        v3 = vr + 0.5 * dlam * k2[1]
        p3 = phi + 0.5 * dlam * k2[2]
        k3 = derivs(r3, v3, p3)
        r4 = r + dlam * k3[0]
        v4 = vr + dlam * k3[1]
        p4 = phi + dlam * k3[2]
        k4 = derivs(r4, v4, p4)
        r += dlam / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        vr += dlam / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        phi += dlam / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        lam += dlam
        pts.append((lam, r, phi))

    return pts


def find_turning_point(xi_fn, b, r_lo=0.5, r_hi=1e4):
    """Find r_tp where r/D(r) = b via bisection. Returns None if none."""
    def f(r):
        return r / D(r, xi_fn) - b

    if f(r_lo) > 0 and f(r_hi) > 0:
        return None
    if f(r_lo) < 0 and f(r_hi) < 0:
        return None
    for _ in range(80):
        rm = 0.5 * (r_lo + r_hi)
        if f(rm) * f(r_lo) <= 0:
            r_hi = rm
        else:
            r_lo = rm
    return 0.5 * (r_lo + r_hi)
