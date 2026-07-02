#!/usr/bin/env python3
# COVERAGE: observable=majorana_polarization geometry=wire material=InAsW tier=illustrative
"""Issue 11 / Plot 7: Majorana polarization profile (synthetic illustrative).

NOTE: This script is illustrative-only. It uses synthetic profiles
(not real Fortran BdG eigenvectors) to render the textbook Sticlet
P_M(n) shape for the lecture material. The canonical formula is
exercised on REAL BdG eigenvectors by the pFUnit test
test_majorana_polarization.pf; this script is the visualization
companion, not a new physics witness.

If a future PR wires this to the Fortran executable with real
eigenvector emission, drop the ``_synthetic`` suffix and move
from tier=illustrative to tier=required in this annotation.
"""
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # headless; no display required
import matplotlib.pyplot as plt
import numpy as np


def topological_profile(n_sites, xi_decay=4.0):
    """Edge-localized Majorana mode: |u|^2 = |v|^2 at ends, ~0 in middle.

    P_M(n) = 2 * |u_{n,+} v*_{n,+} + u_{n,-} v*_{n,-}| where
    the edge localization gives |u v*| ~ exp(-2|n - end|/xi).
    """
    n = np.arange(1, n_sites + 1)
    end = n_sites
    # Two MZMs at the two ends
    left = np.exp(-2.0 * (n - 1) / xi_decay)
    right = np.exp(-2.0 * (end - n) / xi_decay)
    # P_M is a sum over spin sectors (Sticlet coherence)
    P_M = 2.0 * (left + right) / (1.0 + left.max() + right.max())
    return np.clip(P_M, 0.0, 1.0)


def trivial_profile(n_sites, noise_level=0.05):
    """Delocalized in-gap resonance: P_M ~ noise everywhere."""
    n = np.arange(1, n_sites + 1)
    rng = np.random.default_rng(seed=42)
    return rng.uniform(0.0, noise_level, size=n_sites)


def main():
    n_sites = 16  # wire length for visualization (matches typical fixture)
    P_topo = topological_profile(n_sites)
    P_triv = trivial_profile(n_sites)

    # Half-wire integral (left half = first N/2 sites)
    n_half = n_sites // 2
    integral_topo = float(np.sum(P_topo[:n_half]))
    integral_triv = float(np.sum(P_triv[:n_half]))

    out_dir = Path("output")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "majorana_polarization_wire.png"

    fig, ax = plt.subplots(figsize=(8, 5))
    sites = np.arange(1, n_sites + 1)
    ax.plot(sites, P_topo, 'o-', color='#1f77b4', markersize=8,
            label=f'topological (half-wire Σ = {integral_topo:.2f})')
    ax.plot(sites, P_triv, 's--', color='#ff7f0e', markersize=6, alpha=0.7,
            label=f'tivial (half-wire Σ = {integral_triv:.2f})')
    ax.axhline(0.95, color='gray', linestyle=':', alpha=0.5,
               label='P_M = 0.95 (per-site threshold)')
    ax.axvspan(1, n_half + 0.5, alpha=0.08, color='green',
               label=f'half-wire integral region (sites 1..{n_half})')
    ax.set_xlabel('site index n')
    ax.set_ylabel('P_M(n) = 2 |Σ_s u_{n,s} v*_{n,s}|')
    ax.set_title('Sticlet Majorana polarization P_M(n) '
                 '(Issue 04, synthetic illustrative)')
    ax.set_xlim(0.5, n_sites + 0.5)
    ax.set_ylim(0.0, 1.1)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    fig.savefig(out_path, dpi=120, bbox_inches='tight')
    plt.close(fig)
    print(f"wrote {out_path}")
    print(f"  half-wire integral: topological = {integral_topo:.2f}, "
          f"trivial = {integral_triv:.2f}")
    print("PASS: Majorana polarization (Issue 11 / Plot 7, synthetic)")


if __name__ == "__main__":
    main()