#!/usr/bin/env python3
"""Generate all Chapter 33 figures for the SSZ book."""
import math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ssz_trajectories.constants import BLEND_LO, BLEND_HI
from ssz_trajectories import (
    PHI, LN_PHI, c_light, RS_DEFAULT,
    Xi_strong, Xi_weak, Xi_hard, Xi_blend, D, s, dD_dr,
    x_local, N_level, x_cumulative, count_jumps,
    integrate_null_radial, integrate_timelike_infall,
    integrate_null_geodesic, find_turning_point,
    analyze_orbit, deflection_angle, bridge_identity,
    proper_radial_length, tortoise_coordinate,
)

OUT = r"H:\SSZ_BOOK_PROJECT\05_OUTPUT\figures\ch33_trajectories"
OUT2 = r"H:\SSZ_BOOK_PROJECT\03_WRITING\manuscript\figures\ch33_trajectories"
os.makedirs(OUT, exist_ok=True)
os.makedirs(OUT2, exist_ok=True)

STYLE = dict(linewidth=1.8)
plt.rcParams.update({
    "figure.dpi": 180, "savefig.dpi": 180, "font.size": 10,
    "axes.titlesize": 11, "axes.labelsize": 10,
    "figure.figsize": (7, 4.5),
})

def save(fig, name):
    for d in (OUT, OUT2):
        fig.savefig(os.path.join(d, name), bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {name}")

xi = lambda r: Xi_blend(r)

# ================================================================
# Fig 33.01  Coordinate velocity vs physical velocity
# ================================================================
def fig01():
    r = np.linspace(0.5, 30, 600)
    v_coord = np.array([c_light * D(ri, xi)**2 for ri in r])
    v_phys  = np.array([c_light * D(ri, xi) for ri in r])
    fig, ax = plt.subplots()
    ax.plot(r, v_coord / c_light, label=r"$v_\mathrm{coord}/c = D^2(r)$", color="#2563eb", **STYLE)
    ax.plot(r, v_phys / c_light, label=r"$v_\mathrm{phys}/c = D(r)$  (bridge)", color="#dc2626", **STYLE)
    ax.axhline(1, ls=":", color="gray", lw=0.8)
    ax.axvline(1, ls="--", color="gray", lw=0.8, label=r"$r = r_s$")
    ax.set_xlabel(r"$r / r_s$")
    ax.set_ylabel(r"velocity / $c$")
    ax.set_title("Fig 33.1 \u2014 Coordinate vs Physical Radial Velocity (Null Geodesic)")
    ax.legend(fontsize=9)
    ax.set_xlim(0.5, 30)
    ax.set_ylim(0, 1.05)
    save(fig, "fig_33_01.png")

# ================================================================
# Fig 33.02  Physical distance and tortoise to boundary
# ================================================================
def fig02():
    r_vals = np.linspace(1.01, 20, 200)
    ell_ssz = [proper_radial_length(xi, ri, 1.0) for ri in r_vals]
    tort_ssz = [tortoise_coordinate(xi, ri, 1.0) for ri in r_vals]
    def D_gr(r):
        return math.sqrt(max(1.0 - 1.0/r, 1e-12))
    ell_gr = []
    tort_gr = []
    for ri in r_vals:
        n = 10000
        h = (ri - 1.001) / n
        e_acc, t_acc = 0.0, 0.0
        for j in range(n+1):
            rr = 1.001 + j * h
            w = 0.5 if (j == 0 or j == n) else 1.0
            dg = D_gr(rr)
            e_acc += w * (1.0 / dg) * h
            t_acc += w * (1.0 / dg**2) * h
        ell_gr.append(e_acc)
        tort_gr.append(t_acc)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
    ax1.plot(r_vals, ell_ssz, color="#2563eb", label="SSZ", **STYLE)
    ax1.plot(r_vals, ell_gr, color="#dc2626", ls="--", label="GR", **STYLE)
    ax1.set_xlabel(r"$r_0 / r_s$")
    ax1.set_ylabel(r"$\ell(r_0, r_s)$ / $r_s$")
    ax1.set_title(r"Physical distance to $r_s$")
    ax1.legend()
    ax2.plot(r_vals, tort_ssz, color="#2563eb", label="SSZ (finite)", **STYLE)
    ax2.plot(r_vals, tort_gr, color="#dc2626", ls="--", label=r"GR ($\to\infty$)", **STYLE)
    ax2.set_xlabel(r"$r_0 / r_s$")
    ax2.set_ylabel(r"$r^*(r_0, r_s)$ / $r_s$")
    ax2.set_title(r"Tortoise coordinate to $r_s$")
    ax2.set_ylim(0, 80)
    ax2.legend()
    fig.suptitle("Fig 33.2 \u2014 Finiteness at the Boundary: SSZ vs GR", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "fig_33_02.png")

# ================================================================
# Fig 33.03  Non-radial photon orbits (polar plot)
# ================================================================
def fig03():
    b_vals = [3, 5, 10, 50, 200]
    colors = ["#dc2626", "#ea580c", "#2563eb", "#16a34a", "#7c3aed"]
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(6, 6))
    for b, col in zip(b_vals, colors):
        r0 = max(b * 3, 300)
        dlam = 2e-9 if b < 20 else 5e-8
        lam_max = 2 * r0 / c_light
        pts = integrate_null_geodesic(xi, b, r0, lam_max, dlam)
        theta = [p[2] for p in pts]
        radii = [p[1] for p in pts]
        ax.plot(theta, radii, color=col, lw=1.4, label=f"b={b}")
    ax.set_rmax(350)
    ax.set_title("Fig 33.3 \u2014 Non-Radial Null Geodesics", pad=20)
    ax.legend(loc="upper right", fontsize=8, bbox_to_anchor=(1.3, 1.0))
    save(fig, "fig_33_03.png")

# ================================================================
# Fig 33.04  x_local along orbit (b=5)
# ================================================================
def fig04():
    b = 5
    r0 = 300
    dlam = 2e-9
    lam_max = 2 * r0 / c_light
    pts = integrate_null_geodesic(xi, b, r0, lam_max, dlam)
    info = analyze_orbit(xi, pts)
    lam_arr = [p[0] * 1e6 for p in pts]
    fig, ax = plt.subplots()
    ax.plot(lam_arr, info["xl"], color="#2563eb", **STYLE)
    ax.axhline(0.5, ls="--", color="#dc2626", lw=1, label="N=0/1 boundary (x=0.5)")
    ax.fill_between(lam_arr, 0, 0.5, alpha=0.08, color="blue")
    ax.fill_between(lam_arr, 0.5, 1.5, alpha=0.08, color="red")
    ax.text(lam_arr[len(lam_arr)//4], 0.2, "N = 0", fontsize=9, color="#2563eb")
    ax.text(lam_arr[len(lam_arr)//2], 0.8, "N = 1", fontsize=9, color="#dc2626")
    ax.set_xlabel(r"Affine parameter $\lambda$ ($\times 10^{-6}$)")
    ax.set_ylabel(r"$x_\mathrm{local}(r)$")
    ax.set_title(f"Fig 33.4 \u2014 Local Embedding Along Orbit (b = {b} $r_s$)")
    ax.legend(fontsize=9)
    save(fig, "fig_33_04.png")

# ================================================================
# Fig 33.05  x_local vs x_cum
# ================================================================
def fig05():
    b = 5
    r0 = 300
    dlam = 2e-9
    lam_max = 2 * r0 / c_light
    pts = integrate_null_geodesic(xi, b, r0, lam_max, dlam)
    info = analyze_orbit(xi, pts)
    lam_arr = [p[0] * 1e6 for p in pts]
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    l1, = ax1.plot(lam_arr, info["xl"], color="#2563eb", label=r"$x_\mathrm{local}$ (bounded)", **STYLE)
    l2, = ax2.plot(lam_arr, info["xc"], color="#dc2626", label=r"$x_\mathrm{cum}$ (growing)", **STYLE)
    ax1.set_xlabel(r"Affine parameter $\lambda$ ($\times 10^{-6}$)")
    ax1.set_ylabel(r"$x_\mathrm{local}$", color="#2563eb")
    ax2.set_ylabel(r"$x_\mathrm{cum}$", color="#dc2626")
    ax1.set_title(f"Fig 33.5 \u2014 Local vs Cumulative Embedding (b = {b} $r_s$)")
    ax1.legend(handles=[l1, l2], fontsize=9, loc="center left")
    save(fig, "fig_33_05.png")

# ================================================================
# Fig 33.06  Hard switch vs C2 blend
# ================================================================
def fig06():
    r = np.linspace(50, 160, 1000)
    xi_h = [Xi_hard(ri) for ri in r]
    xi_b = [Xi_blend(ri) for ri in r]
    xl_h = [x_local(ri, Xi_hard) for ri in r]
    xl_b = [x_local(ri, Xi_blend) for ri in r]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
    ax1.plot(r, xi_h, color="#dc2626", ls="--", label=r"$\Xi_\mathrm{hard}$", **STYLE)
    ax1.plot(r, xi_b, color="#2563eb", label=r"$\Xi_\mathrm{blend}$ (C\u00b2)", **STYLE)
    ax1.axvspan(BLEND_LO, BLEND_HI, alpha=0.1, color="orange", label="blend zone")
    ax1.set_xlabel(r"$r / r_s$")
    ax1.set_ylabel(r"$\Xi(r)$")
    ax1.set_title(r"$\Xi$ profiles")
    ax1.legend(fontsize=8)

    ax2.plot(r, xl_h, color="#dc2626", ls="--", label=r"$x_\mathrm{hard}$", **STYLE)
    ax2.plot(r, xl_b, color="#2563eb", label=r"$x_\mathrm{blend}$ (C\u00b2)", **STYLE)
    ax2.axvspan(BLEND_LO, BLEND_HI, alpha=0.1, color="orange")
    ax2.axhline(0.5, ls=":", color="gray", lw=0.8)
    ax2.set_xlabel(r"$r / r_s$")
    ax2.set_ylabel(r"$x_\mathrm{local}(r)$")
    ax2.set_title(r"$x_\mathrm{local}$ profiles")
    ax2.legend(fontsize=8)

    fig.suptitle("Fig 33.6 \u2014 Hard Switch vs C\u00b2 Hermite Blend", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "fig_33_06.png")

# ================================================================
# Fig 33.07 (NEW)  Xi function profiles \u2014 all 4 variants
# ================================================================
def fig07():
    r = np.linspace(0.1, 300, 2000)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    ax1.plot(r, [Xi_strong(ri) for ri in r], color="#dc2626", label=r"$\Xi_\mathrm{strong}$", **STYLE)
    ax1.plot(r, [Xi_weak(ri) for ri in r], color="#16a34a", label=r"$\Xi_\mathrm{weak}$", **STYLE)
    ax1.plot(r, [Xi_hard(ri) for ri in r], color="#ea580c", ls="--", label=r"$\Xi_\mathrm{hard}$", **STYLE)
    ax1.plot(r, [Xi_blend(ri) for ri in r], color="#2563eb", label=r"$\Xi_\mathrm{blend}$ (C\u00b2)", lw=2.2)
    ax1.set_xlabel(r"$r / r_s$")
    ax1.set_ylabel(r"$\Xi(r)$")
    ax1.set_title(r"Segment density $\Xi(r)$")
    ax1.legend(fontsize=8)
    ax1.set_xlim(0, 300)
    ax1.set_ylim(0, 1.05)

    ax2.plot(r, [D(ri, Xi_blend) for ri in r], color="#2563eb", label=r"$D(r)$ (blend)", **STYLE)
    ax2.plot(r, [s(ri, Xi_blend) for ri in r], color="#dc2626", label=r"$s(r)$ (blend)", **STYLE)
    ax2.axhline(1, ls=":", color="gray", lw=0.8)
    ax2.set_xlabel(r"$r / r_s$")
    ax2.set_ylabel("value")
    ax2.set_title(r"Dilation $D(r)$ and segment factor $s(r)$")
    ax2.legend(fontsize=8)
    ax2.set_xlim(0, 30)

    fig.suptitle("Fig 33.7 \u2014 SSZ Xi Function Profiles and Metric Quantities", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "fig_33_07.png")

# ================================================================
# Fig 33.08 (NEW)  Gravitational deflection vs impact parameter
# ================================================================
def fig08():
    b_vals = np.linspace(2.5, 118, 150)
    defl_ssz = []
    defl_gr = []
    r0 = 5000.0
    for b in b_vals:
        dlam = 2e-9 if b < 20 else 5e-8
        lam_max = 3.0 * r0 / c_light
        pts = integrate_null_geodesic(xi, b, r0, lam_max, dlam)
        phi_tot = abs(pts[-1][2] - pts[0][2])
        delta = deflection_angle(phi_tot)
        defl_ssz.append(math.degrees(delta))
        defl_gr.append(math.degrees(2.0 / b))  # GR weak-field: 2r_s/b

    fig, ax = plt.subplots()
    ax.plot(b_vals, defl_ssz, color="#2563eb",
            label=r"SSZ $\Delta\phi - \pi$", **STYLE)
    ax.plot(b_vals, defl_gr, color="#dc2626", ls="--",
            label=r"GR weak-field $2r_s/b$", **STYLE)
    ax.set_xlabel(r"Impact parameter $b\;/\;r_s$")
    ax.set_ylabel(r"Deflection $\delta$ (degrees)")
    ax.set_title(r"Fig 33.8 \u2014 Deflection $\Delta\phi - \pi$")
    ax.legend(fontsize=9)
    ax.set_xlim(2.5, 118)
    ax.set_ylim(0, None)
    save(fig, "fig_33_08.png")

# ================================================================
if __name__ == "__main__":
    print("Generating Chapter 33 figures...")
    fig01(); print("  [1/8] velocity comparison")
    fig02(); print("  [2/8] finiteness")
    fig03(); print("  [3/8] non-radial orbits")
    fig04(); print("  [4/8] x_local along orbit")
    fig05(); print("  [5/8] x_local vs x_cum")
    fig06(); print("  [6/8] hard vs blend")
    fig07(); print("  [7/8] Xi profiles")
    fig08(); print("  [8/8] deflection")
    print("Done. All 8 figures generated.")
