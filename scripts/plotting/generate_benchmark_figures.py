#!/usr/bin/env python3
"""Generate benchmark figures for GaAs/AlGaAs QW and InAs/GaSb broken-gap QW."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
FIGDIR = REPO / "docs" / "figures"
FIGDIR.mkdir(exist_ok=True)

# ─── GaAs/AlGaAs QW subband dispersion (dual-panel: VB + CB zoom) ───
eig = np.loadtxt(REPO / "docs/benchmarks/output_gaas_algaas/eigenvalues.dat",
                 comments="#")
k = eig[:, 0]
bands = eig[:, 1:]

# Material band edges (from parameters.f90)
gaas_ev = -0.8
gaas_ec = 0.719
algaas_ev = -0.959
algaas_ec = 1.018

fig, (ax_vb, ax_cb) = plt.subplots(1, 2, figsize=(11, 5))

# VB panel (energies relative to GaAs VB edge, in meV)
for i in range(bands.shape[1]):
    if bands[0, i] < 0:
        ax_vb.plot(k, (bands[:, i] - gaas_ev) * 1000, color="crimson", lw=1.3)

ax_vb.axhline(0, color="gray", ls="--", lw=0.8, label="GaAs VB edge")
ax_vb.axhline((algaas_ev - gaas_ev) * 1000, color="gray", ls=":", lw=0.8,
              label=f"Barrier VB ({(algaas_ev - gaas_ev)*1000:.0f} meV)")
ax_vb.set_xlabel(r"$k_\parallel$ (1/A)", fontsize=11)
ax_vb.set_ylabel("Energy from GaAs VB edge (meV)", fontsize=11)
ax_vb.set_title("Valence subbands", fontsize=12)
ax_vb.legend(fontsize=8)
ax_vb.invert_yaxis()

# CB panel (energies relative to GaAs CB edge, in meV)
for i in range(bands.shape[1]):
    if bands[0, i] > 0:
        ax_cb.plot(k, (bands[:, i] - gaas_ec) * 1000, color="steelblue", lw=1.3)

ax_cb.axhline(0, color="gray", ls="--", lw=0.8, label="GaAs CB edge")
ax_cb.axhline((algaas_ec - gaas_ec) * 1000, color="gray", ls=":", lw=0.8,
              label=f"Barrier CB ({(algaas_ec - gaas_ec)*1000:.0f} meV)")
ax_cb.set_xlabel(r"$k_\parallel$ (1/A)", fontsize=11)
ax_cb.set_ylabel("Energy from GaAs CB edge (meV)", fontsize=11)
ax_cb.set_title("Conduction subbands", fontsize=12)
ax_cb.legend(fontsize=8)

fig.suptitle(r"GaAs/Al$_{0.3}$Ga$_{0.7}$As QW (10 nm well) — 8-band k$\cdot$p",
             fontsize=13, y=1.02)
fig.tight_layout()
fig.savefig(FIGDIR / "benchmark_gaas_algaas_qw.png", dpi=200, bbox_inches="tight")
print("Saved benchmark_gaas_algaas_qw.png")

# ─── InAsW/GaSbW broken-gap QW subband dispersion ───
eig2 = np.loadtxt(REPO / "docs/benchmarks/output_inasw_gasbw/eigenvalues.dat",
                  comments="#")
k2 = eig2[:, 0]
bands2 = eig2[:, 1:]

fig, ax = plt.subplots(figsize=(8, 6))
for i in range(bands2.shape[1]):
    e0 = bands2[0, i]
    if e0 < -0.2:
        color = "crimson"
    elif e0 < 0.4:
        color = "darkorange"
    else:
        color = "steelblue"
    ax.plot(k2, bands2[:, i], color=color, lw=1.2)

ax.axhline(-0.03, color="green", ls="--", lw=0.8, label="GaSbW EV = -0.03")
ax.axhline(-0.172, color="purple", ls="--", lw=0.8, label="InAsW EC = -0.172")
ax.axhspan(-0.03, -0.172, alpha=0.1, color="green", label="Broken-gap overlap")
ax.set_xlabel(r"$k_\parallel$ (1/A)", fontsize=12)
ax.set_ylabel("Energy (eV)", fontsize=12)
ax.set_title("AlSbW/InAsW(15A)/GaSbW(10A) broken-gap QW — 8-band k.p", fontsize=12)
ax.legend(fontsize=8, loc="upper left")
fig.tight_layout()
fig.savefig(FIGDIR / "benchmark_inasw_gasbw_broken_gap.png", dpi=200)
print("Saved benchmark_inasw_gasbw_broken_gap.png")

# ─── Bulk GaAs band structure ───
eig3 = np.loadtxt(REPO / "tests/regression/data/bulk_gaas_kx/eigenvalues.dat",
                  comments="#")
k3 = eig3[:, 0]
bands3 = eig3[:, 1:]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5),
                                gridspec_kw={"width_ratios": [1.2, 1]})
labels = ["SO", "HH", "LH", "LH", "HH", "SO", "CB", "CB"]
colors3 = ["darkorange", "crimson", "orangered", "orangered",
           "crimson", "darkorange", "steelblue", "steelblue"]

for i in range(bands3.shape[1]):
    ax1.plot(k3, bands3[:, i], color=colors3[i], lw=1.5)

ax1.annotate(f"$E_g$ = {1.519:.3f} eV", xy=(0, 0.76), fontsize=9, ha="center")
ax1.annotate(f"$\\Delta_{{SO}}$ = {0.341:.3f} eV", xy=(0.02, -0.17), fontsize=9)
ax1.set_xlabel("k (1/A)", fontsize=12)
ax1.set_ylabel("Energy (eV)", fontsize=12)
ax1.set_title("Bulk GaAs 8-band k.p", fontsize=13)

kmask = k3 <= 0.05
for i in range(bands3.shape[1]):
    ax2.plot(k3[kmask], bands3[kmask, i], color=colors3[i], lw=1.5)
ax2.set_ylim(-0.05, 0.05)
ax2.set_xlabel("k (1/A)", fontsize=12)
ax2.set_title(r"Zoom near $\Gamma$ (VB region)", fontsize=11)

fig.tight_layout()
fig.savefig(FIGDIR / "benchmark_bulk_gaas.png", dpi=200)
print("Saved benchmark_bulk_gaas.png")

# ─── g-factor bar chart ───
fig, ax = plt.subplots(figsize=(6, 4))
methods = ["This code\n(8-band)", "Roth formula", "Experiment", "14-band\n(Winkler)"]
values = [-0.315, -0.44, -0.44, -0.44]
colors_bar = ["steelblue", "coral", "forestgreen", "mediumpurple"]
bars = ax.bar(methods, values, color=colors_bar, edgecolor="black", lw=0.8)
for bar, val in zip(bars, values):
    ax.text(bar.get_x() + bar.get_width()/2, val + 0.01,
            f"{val:.3f}", ha="center", va="bottom", fontsize=10)
ax.axhline(0, color="gray", lw=0.5)
ax.set_ylabel("g* (GaAs CB)", fontsize=12)
ax.set_title("Conduction-band g-factor: 8-band vs literature", fontsize=13)
fig.tight_layout()
fig.savefig(FIGDIR / "benchmark_gfactor_comparison.png", dpi=200)
print("Saved benchmark_gfactor_comparison.png")

print("\nAll benchmark figures generated.")
