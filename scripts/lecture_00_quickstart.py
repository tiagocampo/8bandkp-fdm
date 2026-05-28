#!/usr/bin/env python3
"""Lecture 00: Quickstart -- Bulk GaAs band structure demonstration.

Runs the 8-band k.p solver for bulk GaAs at k=0 (gamma point) and along
the kx dispersion direction. Prints the 8-band spectrum with band labels
and generates an E(k) overlay plot saved to docs/lecture/figures/.
"""
import os
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parent.parent
BUILD_DIR = REPO / "build"
CONFIGS_DIR = REPO / "tests" / "regression" / "configs"
FIGURES_DIR = REPO / "docs" / "lecture" / "figures"

sys.path.insert(0, str(REPO / "tests" / "integration"))
from star_helpers import run_exe, parse_eigenvalues, compare_value

# 8-band zincblende basis labels (CLAUDE.md: bands 1-4 = valence HH/LH,
# bands 5-6 = split-off SO, bands 7-8 = conduction CB)
BAND_LABELS = [
    "HH",   # band 1
    "LH",   # band 2
    "LH",   # band 3
    "HH",   # band 4
    "SO",   # band 5
    "SO",   # band 6
    "CB",   # band 7
    "CB",   # band 8
]

# GaAs reference values (Vurgaftman 2001)
GAAS_EG = 1.519    # eV, band gap
GAAS_DELTA_SO = 0.341  # eV, spin-orbit splitting


def print_k0_spectrum(data):
    """Print the 8-band spectrum at k=0 with band labels.

    Args:
        data: list of (k, [eigenvalues]) from parse_eigenvalues.
    """
    print("\n" + "=" * 60)
    print("  Bulk GaAs 8-band spectrum at Gamma (k = 0)")
    print("=" * 60)
    if not data:
        print("  ERROR: No eigenvalue data parsed.")
        return

    _, evals = data[0]
    print(f"\n  {'Band':>4s}  {'Label':>4s}  {'Energy (eV)':>14s}")
    print(f"  {'----':>4s}  {'----':>4s}  {'------------':>14s}")
    for i, (ev, label) in enumerate(zip(evals, BAND_LABELS)):
        print(f"  {i+1:4d}  {label:>4s}  {ev:14.6f}")
    print()
    print(f"  GaAs references: Eg = {GAAS_EG} eV, Delta_SO = {GAAS_DELTA_SO} eV")
    print(f"  Expected: 2x SO at -{GAAS_DELTA_SO}, 4x VB at ~0, 2x CB at {GAAS_EG}")
    print()


def run_and_parse(config_name, work_dir):
    """Run bandStructure with the given config and return parsed eigenvalues.

    Args:
        config_name: filename in tests/regression/configs/
        work_dir: temporary working directory

    Returns:
        list of (k, [eigenvalues]) from parse_eigenvalues
    """
    config_path = CONFIGS_DIR / config_name
    rc, output_dir = run_exe(str(BUILD_DIR), "bandStructure",
                             str(config_path), work_dir)
    if rc != 0:
        print(f"  ERROR: bandStructure returned {rc} for {config_name}")
        sys.exit(1)

    eig_path = os.path.join(output_dir, "eigenvalues.dat")
    if not os.path.isfile(eig_path):
        print(f"  ERROR: eigenvalues.dat not found at {eig_path}")
        sys.exit(1)

    return parse_eigenvalues(eig_path)


def plot_dispersion(dispersion_data):
    """Generate E(k) band structure plot.

    Args:
        dispersion_data: list of (k, [eigenvalues]) from parse_eigenvalues.
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    ks = np.array([d[0] for d in dispersion_data])
    evals = np.array([d[1] for d in dispersion_data])  # shape (nk, 8)

    # Color map: valence=blue, SO=green, CB=red
    band_colors = [
        "royalblue",   # HH 1
        "steelblue",   # LH 2
        "steelblue",   # LH 3
        "royalblue",   # HH 4
        "forestgreen", # SO 5
        "forestgreen", # SO 6
        "crimson",     # CB 7
        "crimson",     # CB 8
    ]
    band_linestyles = [
        "-",   # HH
        "--",  # LH
        "--",  # LH
        "-",   # HH
        "-.",  # SO
        "-.",  # SO
        "-",   # CB
        "-",   # CB
    ]

    for i in range(evals.shape[1]):
        ax.plot(ks, evals[:, i],
                color=band_colors[i],
                ls=band_linestyles[i],
                lw=1.5,
                label=f"Band {i+1}: {BAND_LABELS[i]}")

    ax.axhline(0, color="gray", ls=":", lw=0.5)
    ax.set_xlabel(r"$k_x$ ($\AA^{-1}$)", fontsize=12)
    ax.set_ylabel("Energy (eV)", fontsize=12)
    ax.set_title("Bulk GaAs 8-band k.p Band Structure", fontsize=13)
    ax.legend(fontsize=9, ncol=2, loc="best")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "lecture_00_quickstart.png",
                dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    print("=" * 60)
    print("  Lecture 00: Quickstart -- Bulk GaAs Band Structure")
    print("=" * 60)

    # Ensure output directory exists
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # --- Part 1: k=0 spectrum ---
    print("\n[1/3] Running bandStructure at k=0 (Gamma point)...")
    with tempfile.TemporaryDirectory() as tmpdir:
        k0_data = run_and_parse("bulk_gaas_k0.toml", tmpdir)
    print_k0_spectrum(k0_data)

    # --- Part 2: kx dispersion ---
    print("[2/3] Running bandStructure kx dispersion (0 to 0.05 A^-1)...")
    with tempfile.TemporaryDirectory() as tmpdir:
        dispersion_data = run_and_parse("bulk_gaas_kx_dispersion.toml", tmpdir)
    nk = len(dispersion_data)
    print(f"  Parsed {nk} k-points with {len(dispersion_data[0][1])} bands each")

    # --- Part 3: Plot ---
    print(f"[3/3] Generating E(k) plot...")
    plot_dispersion(dispersion_data)
    print(f"  Saved: {FIGURES_DIR / 'lecture_00_quickstart.png'}")

    print("\n" + "=" * 60)
    print("  Lecture 00 complete.")
    print("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())
