"""Example: radial null outgoing geodesic in SSZ metric."""
from ssz_trajectories import (
    Xi_blend, D, c_light, integrate_null_radial, bridge_identity,
)

# Integrate outgoing null ray from r = 1.01 r_s
pts = integrate_null_radial(Xi_blend, r0=1.01, tmax=2e-6, dt=2e-9, direction=1)

print(f"Steps: {len(pts)}")
print(f"r: {pts[0][1]:.4f} -> {pts[-1][1]:.2f} r_s")
print(f"Monotone: {all(pts[i+1][1] >= pts[i][1] for i in range(len(pts)-1))}")

# Bridge identity spot checks
print("\nBridge identity dl/dt = c*D:")
print(f"{'r/r_s':>8} {'dr/dt':>14} {'dl/dt':>14} {'c*D':>14} {'error':>10}")
for r in [1.01, 2, 10, 100, 1000]:
    b = bridge_identity(Xi_blend, float(r))
    print(f"{r:8.1f} {b['dr_dt']:14.2f} {b['dl_dt']:14.2f} "
          f"{b['cD']:14.2f} {b['relative_error']:10.2e}")
