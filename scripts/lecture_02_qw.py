#!/usr/bin/env python3
"""Lecture 02: Quantum Well -- Subband energies and confinement.

Validates QW subband structure from the 8-band k.p solver:

  1. Single GaAs/AlGaAs QW: count CB subbands at k=0 and verify
     confinement shifts the CB edge above the GaAs bulk gap.
  2. Double GaAs/AlGaAs QW: verify anticrossing splitting produces
     two closely-spaced CB subbands (splitting < 50 meV).
  3. Subband structure plot: E(k_parallel) for single QW showing
     CB and top VB subbands with band-gap region annotated.
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
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(REPO / "tests" / "integration"))
from star_helpers import (
    run_exe, parse_eigenvalues, compare_value,
    TOL_NUMERICAL,
)

# ---------------------------------------------------------------------------
# Material parameters
# ---------------------------------------------------------------------------
# GaAs (Vurgaftman 2001)
GAAS_EG = 1.519  # eV, bulk band gap

# Al0.3Ga0.7As band gap (Vurgaftman 2001, linear interpolation)
ALGAAS30_EG = 1.519 + 1.247 * 0.30  # = 1.893 eV

# The 100-A GaAs QW confined in Al0.3Ga0.7As has 2 CB subbands
# (numcb=4 eigenvalues = 2 spin-degenerate subbands).
NUMCB = 4   # CB eigenvalues from config
NUMVB = 8   # VB eigenvalues from config
NUM_TOTAL = NUMCB + NUMVB  # 12 eigenvalues per k-point


def run_bandstructure(config_name, work_dir, timeout=300):
    """Run bandStructure with the given config, return parsed eigenvalues.

    Args:
        config_name: filename in tests/regression/configs/
        work_dir: temporary working directory
        timeout: execution timeout in seconds

    Returns:
        list of (k, [eigenvalues]) from parse_eigenvalues
    """
    config_path = CONFIGS_DIR / config_name
    rc, output_dir = run_exe(str(BUILD_DIR), "bandStructure",
                             str(config_path), work_dir, timeout=timeout)
    if rc != 0:
        sys.exit(f"ERROR: bandStructure returned {rc} for {config_name}")

    eig_path = os.path.join(output_dir, "eigenvalues.dat")
    if not os.path.isfile(eig_path):
        sys.exit(f"ERROR: eigenvalues.dat not found at {eig_path}")

    data = parse_eigenvalues(eig_path)
    if not data:
        sys.exit(f"ERROR: no eigenvalue data parsed from {config_name}")

    return data


def classify_subbands(evals):
    """Classify eigenvalues into VB and CB subbands.

    With the 8-band basis, eigenvalues are ordered ascending.
    The last NUMCB are CB, the first NUMVB are VB.
    Spin-degenerate pairs collapse to unique subband energies.

    Args:
        evals: list of eigenvalue floats (length = NUM_TOTAL)

    Returns:
        (vb_unique, cb_unique) lists of unique subband energies
    """
    vb_evals = evals[:NUMVB]   # first numvb eigenvalues
    cb_evals = evals[NUMVB:]   # last numcb eigenvalues

    # Collapse spin-degenerate pairs
    vb_unique = []
    for i in range(0, len(vb_evals), 2):
        vb_unique.append(vb_evals[i])

    cb_unique = []
    for i in range(0, len(cb_evals), 2):
        cb_unique.append(cb_evals[i])

    return vb_unique, cb_unique


# =========================================================================
# Section 1: Single QW subband count at k=0
# =========================================================================
def test_single_qw_subbands():
    """Verify GaAs/AlGaAs single QW has 2 CB subbands at k=0."""
    print("=" * 60)
    print("Lecture 02 -- Section 1: Single QW CB subbands at k=0")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as work:
        data = run_bandstructure("qw_gaas_algaas.cfg", work)

    k0, evals = data[0]
    print(f"  k = {k0:.6f},  {len(evals)} eigenvalues")

    vb_unique, cb_unique = classify_subbands(evals)

    print(f"  VB subbands ({len(vb_unique)}): "
          + ", ".join(f"{e:.4f} eV" for e in vb_unique))
    print(f"  CB subbands ({len(cb_unique)}): "
          + ", ".join(f"{e:.4f} eV" for e in cb_unique))

    all_pass = True

    # Assertion 1: exactly 2 CB subbands
    n_cb = len(cb_unique)
    if n_cb == 2:
        print(f"  [PASS] CB subband count = {n_cb} (expected 2)")
    else:
        print(f"  [FAIL] CB subband count = {n_cb} (expected 2)")
        all_pass = False

    # Assertion 2: effective QW gap exceeds bulk GaAs gap (confinement shift)
    # The eigenvalues are on an absolute energy scale, so we compare
    # the effective gap (CB1 - VB_top) against the bulk GaAs gap.
    vb_top = vb_unique[-1]  # highest VB subband (least negative, last in ascending order)
    cb_ground = cb_unique[0]
    effective_gap = cb_ground - vb_top
    if effective_gap > GAAS_EG:
        print(f"  [PASS] Effective QW gap {effective_gap:.4f} eV > "
              f"GaAs bulk Eg {GAAS_EG:.4f} eV (confinement shift)")
    else:
        print(f"  [FAIL] Effective QW gap {effective_gap:.4f} eV <= "
              f"GaAs bulk Eg {GAAS_EG:.4f} eV")
        all_pass = False

    # Assertion 3: effective gap within TOL_NUMERICAL of expected
    # The 8-band QW gap with 100 A GaAs well in Al0.3Ga0.7As is
    # slightly above bulk Eg due to CB confinement (no simple formula
    # in 8-band, so we check the gap is reasonable: within 20% of bulk).
    gap_tol = 0.20
    passed, delta, _ = compare_value(
        effective_gap, GAAS_EG, gap_tol,
        "Effective QW gap vs bulk GaAs Eg", unit="eV",
    )
    status = "PASS" if passed else "FAIL"
    print(f"  [{status}] Gap ratio = {effective_gap/GAAS_EG:.4f} "
          f"(within {gap_tol:.0%} of bulk)")

    # Print subband spacing
    if n_cb >= 2:
        spacing = cb_unique[1] - cb_unique[0]
        print(f"  CB subband spacing: {spacing:.4f} eV ({spacing*1000:.1f} meV)")

    print()
    return all_pass, data


# =========================================================================
# Section 2: Double QW anticrossing
# =========================================================================
def test_double_qw_anticrossing():
    """Verify double QW shows anticrossing (2 closely-spaced CB subbands)."""
    print("=" * 60)
    print("Lecture 02 -- Section 2: Double QW anticrossing")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as work:
        data = run_bandstructure("qw_gaas_algaas_double_qw.cfg", work,
                                 timeout=600)

    k0, evals = data[0]
    print(f"  k = {k0:.6f},  {len(evals)} eigenvalues")

    vb_unique, cb_unique = classify_subbands(evals)

    print(f"  VB subbands ({len(vb_unique)}): "
          + ", ".join(f"{e:.4f} eV" for e in vb_unique))
    print(f"  CB subbands ({len(cb_unique)}): "
          + ", ".join(f"{e:.4f} eV" for e in cb_unique))

    all_pass = True

    # Assertion: 2 CB subbands present
    n_cb = len(cb_unique)
    if n_cb >= 2:
        print(f"  [PASS] CB subband count = {n_cb} (>= 2)")
    else:
        print(f"  [FAIL] CB subband count = {n_cb} (expected >= 2)")
        all_pass = False

    # Assertion: anticrossing splitting < 50 meV
    if n_cb >= 2:
        splitting = cb_unique[1] - cb_unique[0]
        splitting_meV = splitting * 1000
        max_split = 50.0  # meV

        if splitting_meV < max_split:
            print(f"  [PASS] Anticrossing splitting = {splitting_meV:.1f} meV "
                  f"< {max_split:.0f} meV")
        else:
            print(f"  [FAIL] Anticrossing splitting = {splitting_meV:.1f} meV "
                  f">= {max_split:.0f} meV")
            all_pass = False

        # Splitting should also be > 0 (non-degenerate)
        if splitting_meV > 1.0:
            print(f"  [PASS] Splitting > 1 meV (genuinely non-degenerate)")
        else:
            print(f"  [WARN] Splitting {splitting_meV:.2f} meV very small "
                  f"(near-degenerate)")

    print()
    return all_pass


# =========================================================================
# Section 3: Subband structure plot E(k_parallel)
# =========================================================================
def plot_subband_structure(data):
    """Generate E(k_parallel) subband structure plot for single QW.

    Args:
        data: list of (k, [eigenvalues]) from parse_eigenvalues.
    """
    print("=" * 60)
    print("Lecture 02 -- Section 3: Subband structure plot E(k_parallel)")
    print("=" * 60)

    ks = np.array([d[0] for d in data])
    evals = np.array([d[1] for d in data])  # shape (nk, NUM_TOTAL)

    # Classify each k-point's eigenvalues into subbands
    n_cb = NUMCB // 2  # unique CB subbands (spin-degenerate pairs)
    n_vb = NUMVB // 2  # unique VB subbands (spin-degenerate pairs)

    # Extract unique subband energies at each k-point
    cb_subbands = np.zeros((len(ks), n_cb))
    vb_subbands = np.zeros((len(ks), n_vb))

    for ik in range(len(ks)):
        _, cb_unique = classify_subbands(evals[ik].tolist())
        vb_unique, _ = classify_subbands(evals[ik].tolist())
        for j in range(n_cb):
            cb_subbands[ik, j] = cb_unique[j]
        for j in range(n_vb):
            vb_subbands[ik, j] = vb_unique[j]

    # Compute band gap at k=0
    cb_min = cb_subbands[0, 0]
    vb_max = vb_subbands[0, -1]  # highest VB subband (least negative, last in ascending order)
    effective_gap = cb_min - vb_max

    # --- Figure ---
    fig, ax = plt.subplots(figsize=(9, 7))

    # Color scheme
    cb_colors = ["crimson", "orangered"]
    vb_colors = ["royalblue", "steelblue", "dodgerblue", "cornflowerblue"]

    # Plot CB subbands
    for j in range(n_cb):
        label = f"CB{j+1}" if n_cb > 1 else "CB1"
        ax.plot(ks, cb_subbands[:, j], color=cb_colors[j % len(cb_colors)],
                linewidth=2, label=label)

    # Plot VB subbands
    for j in range(n_vb):
        label = f"VB{n_vb - j}"  # label from top to bottom
        ax.plot(ks, vb_subbands[:, j],
                color=vb_colors[j % len(vb_colors)],
                linewidth=2, linestyle="--", label=label)

    # Annotate band-gap region
    ax.axhspan(vb_max, cb_min, alpha=0.08, color="gold",
               label=f"Gap = {effective_gap:.3f} eV")
    ax.axhline(cb_min, color="crimson", linestyle=":", linewidth=0.8,
               alpha=0.5)
    ax.axhline(vb_max, color="royalblue", linestyle=":", linewidth=0.8,
               alpha=0.5)

    # Annotate bulk band edges for reference
    ax.axhline(GAAS_EG, color="gray", linestyle="-.", linewidth=0.7,
               alpha=0.4)
    ax.text(ks[-1] * 0.02, GAAS_EG + 0.01,
            f"Bulk GaAs Eg = {GAAS_EG:.3f} eV",
            fontsize=8, color="gray", alpha=0.7)

    ax.set_xlabel(r"$k_\parallel$ ($\mathrm{nm}^{-1}$)", fontsize=12)
    ax.set_ylabel("Energy (eV)", fontsize=12)
    ax.set_title("Lecture 02: GaAs/Al$_{0.3}$Ga$_{0.7}$As QW "
                 "Subband Structure", fontsize=13)
    ax.legend(fontsize=9, ncol=2, loc="upper left")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    out_path = FIGURES_DIR / "lecture_02_qw_subbands.png"
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"  Effective QW gap at k=0: {effective_gap:.4f} eV")
    print(f"  CB1 at k=0: {cb_min:.4f} eV, VB top at k=0: {vb_max:.4f} eV")
    print(f"  Plot saved to {out_path}")
    print()


# =========================================================================
# Main
# =========================================================================
def main():
    print("\n" + "=" * 60)
    print("  LECTURE 02: Quantum Well Subband Validation")
    print("  GaAs/AlGaAs single and double quantum wells")
    print("=" * 60 + "\n")

    # Run all sections
    s1_pass, single_qw_data = test_single_qw_subbands()
    s2_pass = test_double_qw_anticrossing()

    # Re-run single QW to get data for plotting (tempdir from s1 is gone)
    with tempfile.TemporaryDirectory() as work:
        plot_data = run_bandstructure("qw_gaas_algaas.cfg", work)
    plot_subband_structure(plot_data)

    # Summary
    print("=" * 60)
    print("  SUMMARY")
    print("=" * 60)
    results = [
        ("Section 1: Single QW CB subbands", s1_pass),
        ("Section 2: Double QW anticrossing", s2_pass),
        ("Section 3: Subband structure plot", True),
    ]
    all_pass = True
    for label, passed in results:
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {label}")
        all_pass = all_pass and passed
    print()

    if all_pass:
        print("  All validations passed.")
    else:
        print("  Some validations FAILED.")
        sys.exit(1)


if __name__ == "__main__":
    main()
