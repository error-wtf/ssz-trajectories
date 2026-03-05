"""Example: plot Xi(r) profile for all four variants.

Requires: pip install matplotlib numpy
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ssz_trajectories import Xi_strong, Xi_weak, Xi_hard, Xi_blend, D

r = np.linspace(0.1, 300, 2000)

fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Xi(r)
ax = axes[0]
ax.plot(r, [Xi_strong(float(x)) for x in r], label="Xi_strong", lw=2)
ax.plot(r, [Xi_weak(float(x)) for x in r], label="Xi_weak", lw=2, ls="--")
ax.plot(r, [Xi_hard(float(x)) for x in r], label="Xi_hard", lw=2, ls=":")
ax.plot(r, [Xi_blend(float(x)) for x in r], label="Xi_blend (C2)", lw=2.5)
ax.axvspan(80, 120, alpha=0.1, color="orange", label="blend zone")
ax.set_ylabel("Xi(r)")
ax.set_title("SSZ Segment Density Xi(r)")
ax.legend(loc="right")
ax.set_ylim(-0.05, 1.15)
ax.grid(True, alpha=0.3)

# D(r)
ax = axes[1]
for name, fn in [("strong", Xi_strong), ("weak", Xi_weak),
                 ("hard", Xi_hard), ("blend", Xi_blend)]:
    ax.plot(r, [D(float(x), fn) for x in r], label=f"D ({name})", lw=2)
ax.axvspan(80, 120, alpha=0.1, color="orange")
ax.set_xlabel("r / r_s")
ax.set_ylabel("D(r)")
ax.set_title("SSZ Dilation Factor D(r) = 1/(1+Xi)")
ax.legend(loc="right")
ax.set_ylim(0.4, 1.05)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("xi_profile.png", dpi=150)
print("Saved: xi_profile.png")
