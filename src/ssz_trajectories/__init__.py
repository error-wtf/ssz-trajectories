"""SSZ Trajectories -- geodesic toolkit for Segmented Spacetime.

Authors: Carmen N. Wrede, Lino P. Casu

Provides:
- Xi functions (strong, weak, hard-switch, C2 Hermite blend)
- Metric quantities D(r), s(r), dD/dr
- Phi-Euler embedding (x_local, N-level, cumulative)
- Geodesic integrators (radial null, timelike infall, non-radial null)
- Orbit analysis (deflection, jump counting, bridge identity)
"""
from .constants import PHI, LN_PHI, c_light, RS_DEFAULT
from .xi import Xi_strong, Xi_weak, Xi_hard, Xi_blend, D, s, dD_dr
from .embedding import x_local, N_level, epsilon_residual, x_cumulative, count_jumps
from .integrator import (
    rk4_scalar,
    trapz,
    integrate_null_radial,
    integrate_timelike_infall,
    integrate_null_geodesic,
    find_turning_point,
)
from .analysis import (
    analyze_orbit,
    deflection_angle,
    bridge_identity,
    proper_radial_length,
    tortoise_coordinate,
)

__version__ = "0.1.0"
__all__ = [
    "PHI", "LN_PHI", "c_light", "RS_DEFAULT",
    "Xi_strong", "Xi_weak", "Xi_hard", "Xi_blend", "D", "s", "dD_dr",
    "x_local", "N_level", "epsilon_residual", "x_cumulative", "count_jumps",
    "rk4_scalar", "trapz",
    "integrate_null_radial", "integrate_timelike_infall",
    "integrate_null_geodesic", "find_turning_point",
    "analyze_orbit", "deflection_angle", "bridge_identity",
    "proper_radial_length", "tortoise_coordinate",
]
