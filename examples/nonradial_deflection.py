"""Example: non-radial null geodesic deflection in SSZ metric."""
import math
from ssz_trajectories import (
    Xi_blend, c_light, D,
    integrate_null_geodesic, analyze_orbit, deflection_angle,
)

print("SSZ Non-Radial Null Geodesic Deflection")
print("=" * 50)

b_values = [2, 5, 10, 50, 80, 120, 200]

print(f"\n{'b':>5} {'r_min':>8} {'r_max':>8} {'jumps':>6} {'defl(deg)':>10}")
print("-" * 45)

for b in b_values:
    r0 = max(300.0, b * 4)
    lam_max = r0 / c_light * 12
    pts = integrate_null_geodesic(Xi_blend, float(b), r0, lam_max, 2e-9)
    res = analyze_orbit(Xi_blend, pts)
    defl = deflection_angle(res["phi_total"], float(b), r0, Xi_blend)
    print(f"{b:5d} {res['r_min']:8.2f} {res['r_max']:8.1f} "
          f"{res['xl_jumps']:6d} {math.degrees(defl):10.1f}")
