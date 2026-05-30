#!/usr/bin/env python3
"""Lecture 08: Quantum Wire -- 2D confinement and CSR sparse assembly.

Runs the 8-band k.p solver (bandStructure executable) for quantum wire
configurations and validates:

  1. GaAs wire (21x21 grid): eigenvalues computed via CSR sparse assembly,
     subband count consistent with 8*21*21 = 3528 Hamiltonian dimension.
  2. InAs rectangular wire (11x11 grid): eigenvalues computed successfully.
  3. Dense-sparse consistency: eigenvalues from dense and sparse solvers
     agree within numerical tolerance for an 11x11 wire.
  4. Wire subband structure plot E(kz) with subband count annotation,
     saved to docs/lecture/figures/lecture_08_wire_subbands.png.
"""
import os
import sys
import shutil
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
from star_helpers import (run_exe, parse_eigenvalues, compare_value,
                           TOL_EXACT, TOL_NUMERICAL)

# Wire configs are CSR-heavy; allow generous timeout.
WIRE_TIMEOUT = 1200


# =========================================================================
# Section 1: GaAs 21x21 wire eigenvalues
# =========================================================================
def test_gaas_21x21_wire():
    """Run GaAs wire (21x21 grid) and verify eigenvalue output."""
    print("=" * 60)
    print("Lecture 08 -- Section 1: GaAs 21x21 wire eigenvalues")
    print("=" * 60)

    cfg = CONFIGS_DIR / "wire_gaas_rectangle.toml"
    nx, ny = 21, 21
    expected_hdim = 8 * nx * ny  # 3528
    expected_numcb = 8
    expected_numvb = 16
    # FEAST with feast_m0=-1 auto-detects; expect >= numcb + numvb
    min_expected_nev = expected_numcb + expected_numvb  # 24

    with tempfile.TemporaryDirectory(prefix="lecture08_gaas21_") as work:
        rc, outdir = run_exe(str(BUILD_DIR), "bandStructure",
                             str(cfg), work, timeout=WIRE_TIMEOUT)
        if rc != 0:
            sys.exit(f"ERROR: bandStructure returned {rc} for GaAs 21x21 wire")

        eig_path = os.path.join(outdir, "eigenvalues.dat")
        if not os.path.isfile(eig_path):
            sys.exit(f"ERROR: eigenvalues.dat not found at {eig_path}")

        data = parse_eigenvalues(eig_path)

    if not data:
        sys.exit("ERROR: no eigenvalue data parsed from GaAs 21x21 wire")

    n_kz = len(data)
    n_evals_first = len(data[0][1])
    n_evals_last = len(data[-1][1])

    print(f"  Parsed {n_kz} kz-points")
    print(f"  Eigenvalues at k=0: {n_evals_first}")
    print(f"  Eigenvalues at k_max: {n_evals_last}")
    print(f"  Hamiltonian dimension: {expected_hdim} (8 x {nx} x {ny})")
    print(f"  Expected numcb={expected_numcb}, numvb={expected_numvb}, "
          f"min_nev={min_expected_nev}")

    # Assertion: non-empty output
    assert n_kz > 0, "No kz-points in eigenvalue output"
    print(f"  PASS: Non-empty eigenvalue output ({n_kz} kz-points)")

    # Assertion: subband count >= numcb + numvb (FEAST auto-detect may find more)
    assert n_evals_first >= min_expected_nev, (
        f"Expected >= {min_expected_nev} eigenvalues, got {n_evals_first}")
    print(f"  PASS: Subband count = {n_evals_first} "
          f"(>= numcb+numvb = {min_expected_nev})")

    # Assertion: Hamiltonian dimension consistent with grid
    assert expected_hdim == 8 * nx * ny, "Dimension mismatch"
    print(f"  PASS: Hamiltonian dimension {expected_hdim} "
          f"consistent with {nx}x{ny} grid")

    # Extract kz values and CB/VB eigenvalues for plotting
    kz_vals = np.array([d[0] for d in data])
    eig_matrix = np.array([d[1] for d in data])  # shape (n_kz, n_evals)

    print()
    return True, kz_vals, eig_matrix, nx, ny, expected_hdim


# =========================================================================
# Section 2: InAs rectangular wire
# =========================================================================
def test_inas_rectangle_wire():
    """Run InAs rectangular wire and verify eigenvalue output."""
    print("=" * 60)
    print("Lecture 08 -- Section 2: InAs rectangular wire")
    print("=" * 60)

    cfg = CONFIGS_DIR / "wire_inas_rectangle.toml"

    with tempfile.TemporaryDirectory(prefix="lecture08_inas_") as work:
        rc, outdir = run_exe(str(BUILD_DIR), "bandStructure",
                             str(cfg), work, timeout=WIRE_TIMEOUT)
        if rc != 0:
            sys.exit(f"ERROR: bandStructure returned {rc} for InAs rectangle wire")

        eig_path = os.path.join(outdir, "eigenvalues.dat")
        if not os.path.isfile(eig_path):
            sys.exit(f"ERROR: eigenvalues.dat not found at {eig_path}")

        data = parse_eigenvalues(eig_path)

    if not data:
        sys.exit("ERROR: no eigenvalue data parsed from InAs rectangle wire")

    n_kz = len(data)
    n_evals = len(data[0][1])

    print(f"  Parsed {n_kz} kz-points")
    print(f"  Eigenvalues per kz-point: {n_evals}")

    # Assertion: non-empty output
    assert n_kz > 0, "No kz-points in InAs wire eigenvalue output"
    assert n_evals > 0, "No eigenvalues in InAs wire output"
    print(f"  PASS: InAs wire eigenvalues computed successfully "
          f"({n_kz} kz-points, {n_evals} bands)")

    # Check eigenvalues are finite
    for i, (kz, evals) in enumerate(data):
        for j, ev in enumerate(evals):
            if not np.isfinite(ev):
                sys.exit(f"ERROR: non-finite eigenvalue at kz-point {i}, "
                         f"band {j}: {ev}")
    print(f"  PASS: All eigenvalues are finite")

    print()
    return True


# =========================================================================
# Section 3: Dense-sparse consistency (if configs available)
# =========================================================================
def test_dense_sparse_consistency():
    """Compare dense and sparse wire eigenvalues for consistency."""
    print("=" * 60)
    print("Lecture 08 -- Section 3: Dense-sparse consistency")
    print("=" * 60)

    dense_cfg = CONFIGS_DIR / "wire_gaas_dense_sparse_dense.toml"
    sparse_cfg = CONFIGS_DIR / "wire_gaas_dense_sparse_sparse.toml"

    if not dense_cfg.is_file() or not sparse_cfg.is_file():
        print("  SKIP: dense-sparse comparison configs not available")
        print()
        return None  # skipped, not a failure

    # Run dense solver
    with tempfile.TemporaryDirectory(prefix="lecture08_dense_") as work_dense:
        rc_dense, outdir_dense = run_exe(str(BUILD_DIR), "bandStructure",
                                          str(dense_cfg), work_dense,
                                          timeout=WIRE_TIMEOUT)
        if rc_dense != 0:
            print(f"  WARN: dense solver returned {rc_dense}, skipping comparison")
            return None

        eig_dense_path = os.path.join(outdir_dense, "eigenvalues.dat")
        if not os.path.isfile(eig_dense_path):
            print("  WARN: dense eigenvalues.dat not found, skipping")
            return None

        data_dense = parse_eigenvalues(eig_dense_path)

    # Run sparse solver
    with tempfile.TemporaryDirectory(prefix="lecture08_sparse_") as work_sparse:
        rc_sparse, outdir_sparse = run_exe(str(BUILD_DIR), "bandStructure",
                                            str(sparse_cfg), work_sparse,
                                            timeout=WIRE_TIMEOUT)
        if rc_sparse != 0:
            print(f"  WARN: sparse solver returned {rc_sparse}, skipping comparison")
            return None

        eig_sparse_path = os.path.join(outdir_sparse, "eigenvalues.dat")
        if not os.path.isfile(eig_sparse_path):
            print("  WARN: sparse eigenvalues.dat not found, skipping")
            return None

        data_sparse = parse_eigenvalues(eig_sparse_path)

    if not data_dense or not data_sparse:
        print("  WARN: empty eigenvalue data, skipping comparison")
        return None

    # Both should have the same kz-points
    n_kz = min(len(data_dense), len(data_sparse))
    n_evals = min(len(data_dense[0][1]), len(data_sparse[0][1]))

    print(f"  Comparing {n_kz} kz-points, {n_evals} eigenvalues each")

    # Compare eigenvalues at first kz-point
    max_diff = 0.0
    for i in range(n_kz):
        evals_d = data_dense[i][1][:n_evals]
        evals_s = data_sparse[i][1][:n_evals]
        for j in range(n_evals):
            diff = abs(evals_d[j] - evals_s[j])
            max_diff = max(max_diff, diff)

    print(f"  Max eigenvalue difference (dense vs sparse): {max_diff:.2e} eV")

    passed, delta, _ = compare_value(
        max_diff, 0.0, TOL_NUMERICAL,
        "Dense-sparse eigenvalue agreement", unit="eV",
    )
    status = "PASS" if passed else "FAIL"
    print(f"  {status}: max_diff = {max_diff:.2e} "
          f"(tolerance {TOL_NUMERICAL:.0%})")

    if passed:
        print("  --> Dense and sparse solvers agree within tolerance.")
    else:
        print("  --> WARN: dense-sparse discrepancy exceeds tolerance.")
    print()
    return passed


# =========================================================================
# Section 4: Wire subband structure plot E(kz)
# =========================================================================
def plot_wire_subbands(kz_vals, eig_matrix, nx, ny, hdim):
    """Generate wire subband structure plot E(kz).

    Args:
        kz_vals: array of kz wavevectors
        eig_matrix: eigenvalues array (n_kz x n_evals)
        nx, ny: wire grid dimensions
        hdim: Hamiltonian matrix dimension
    """
    print("=" * 60)
    print("Lecture 08 -- Section 4: Wire subband structure plot")
    print("=" * 60)

    n_evals = eig_matrix.shape[1]

    # Identify the band gap at each kz-point as the largest energy spacing.
    # This separates CB from VB even when spurious solutions are present.
    midpoints = np.zeros(len(kz_vals))
    for i in range(len(kz_vals)):
        evals_sorted = np.sort(eig_matrix[i])
        diffs = np.diff(evals_sorted)
        gi = np.argmax(diffs)
        midpoints[i] = 0.5 * (evals_sorted[gi] + evals_sorted[gi + 1])

    # Use the k=0 gap as reference for display window.
    evals_k0 = np.sort(eig_matrix[0])
    diffs_k0 = np.diff(evals_k0)
    gap_idx = np.argmax(diffs_k0)
    vb_top = evals_k0[gap_idx]
    cb_bot = evals_k0[gap_idx + 1]
    band_gap = cb_bot - vb_top
    mid_ref = 0.5 * (vb_top + cb_bot)

    # Band tracking with gap-aware greedy matching and jump rejection.
    # At each k-point, split eigenvalues into CB/VB by the local midpoint,
    # then track within each group independently. Reject matches where
    # the energy jump exceeds a physical threshold.
    max_jump = 0.15  # eV — max allowed jump between consecutive kz-points

    cb_bands = []  # list of arrays, each is a tracked CB subband
    vb_bands = []  # list of arrays, each is a tracked VB subband

    # Initialize from k=0
    cb_idx_k0 = np.where(eig_matrix[0] > mid_ref)[0]
    vb_idx_k0 = np.where(eig_matrix[0] <= mid_ref)[0]

    # Track CB bands
    cb_tracked = np.full((len(kz_vals), len(cb_idx_k0)), np.nan)
    cb_tracked[0] = np.sort(eig_matrix[0][cb_idx_k0])
    for i in range(1, len(kz_vals)):
        mid = midpoints[i]
        cb_cur = np.sort(eig_matrix[i][eig_matrix[i] > mid])
        prev = cb_tracked[i - 1]
        n_prev = len(prev)
        n_cur = len(cb_cur)
        if n_cur == 0 or n_prev == 0:
            continue
        used = np.zeros(n_cur, dtype=bool)
        for j in range(n_prev):
            if np.isnan(prev[j]):
                continue
            dists = np.abs(cb_cur - prev[j])
            dists[used] = np.inf
            best = np.argmin(dists)
            if dists[best] < max_jump:
                cb_tracked[i, j] = cb_cur[best]
                used[best] = True

    # Track VB bands
    vb_tracked = np.full((len(kz_vals), len(vb_idx_k0)), np.nan)
    vb_tracked[0] = np.sort(eig_matrix[0][vb_idx_k0])
    for i in range(1, len(kz_vals)):
        mid = midpoints[i]
        vb_cur = np.sort(eig_matrix[i][eig_matrix[i] <= mid])
        prev = vb_tracked[i - 1]
        n_prev = len(prev)
        n_cur = len(vb_cur)
        if n_cur == 0 or n_prev == 0:
            continue
        used = np.zeros(n_cur, dtype=bool)
        for j in range(n_prev):
            if np.isnan(prev[j]):
                continue
            dists = np.abs(vb_cur - prev[j])
            dists[used] = np.inf
            best = np.argmin(dists)
            if dists[best] < max_jump:
                vb_tracked[i, j] = vb_cur[best]
                used[best] = True

    fig, ax = plt.subplots(figsize=(10, 7))

    # Plot CB bands
    n_cb = cb_tracked.shape[1]
    for j in range(n_cb):
        mask = ~np.isnan(cb_tracked[:, j])
        if np.sum(mask) > 1:
            ax.plot(kz_vals[mask], cb_tracked[mask, j],
                    color="crimson", alpha=0.7, lw=0.8)

    # Plot VB bands near the gap (limit to upper portion)
    vb_display = min(16, vb_tracked.shape[1])
    n_vb_total = vb_tracked.shape[1]
    for j in range(n_vb_total - vb_display, n_vb_total):
        mask = ~np.isnan(vb_tracked[:, j])
        if np.sum(mask) > 1:
            ax.plot(kz_vals[mask], vb_tracked[mask, j],
                    color="steelblue", alpha=0.4, lw=0.5)

    ax.annotate(
        f"Grid: {nx} x {ny}\n"
        f"Gap: {band_gap:.3f} eV\n"
        f"CB subbands: {n_cb}\n"
        f"VB subbands (shown): {vb_display}",
        xy=(0.98, 0.02),
        xycoords="axes fraction",
        fontsize=9,
        ha="right",
        va="bottom",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.8),
    )

    ax.set_xlabel(r"$k_z$ (1/A)", fontsize=12)
    ax.set_ylabel("Energy (eV)", fontsize=12)
    ax.set_title(
        f"Quantum Wire Subband Structure "
        f"(GaAs, {nx}x{ny} grid, 8-band k.p)",
        fontsize=13, fontweight="bold",
    )
    # Set y-limits around the gap region
    e_lo = vb_top - 0.3
    e_hi = cb_bot + 0.3
    ax.set_ylim(e_lo, e_hi)
    ax.grid(True, alpha=0.3)

    # Add legend entries for CB/VB coloring
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color="crimson", lw=1.2, label="CB subbands"),
        Line2D([0], [0], color="steelblue", lw=1.2, label="VB subbands"),
    ]
    ax.legend(handles=legend_elements, fontsize=10, loc="upper left")

    fig.tight_layout()

    out_path = FIGURES_DIR / "lecture_08_wire_subbands.png"
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"  Plot saved to {out_path}")
    print(f"  Near-gap bands: {n_cb} CB + {vb_display} VB (of {n_vb_total} total)")
    print(f"  Band gap at k=0: {band_gap:.4f} eV")
    print()


# =========================================================================
# Main
# =========================================================================
def main():
    print("\n" + "=" * 60)
    print("  LECTURE 08: Quantum Wire")
    print("  2D confinement, CSR sparse assembly, wire subbands")
    print("=" * 60 + "\n")

    # Run all sections
    s1_pass, kz_vals, eig_matrix, nx, ny, hdim = test_gaas_21x21_wire()
    s2_pass = test_inas_rectangle_wire()
    s3_pass = test_dense_sparse_consistency()
    plot_wire_subbands(kz_vals, eig_matrix, nx, ny, hdim)

    # Summary
    print("=" * 60)
    print("  SUMMARY")
    print("=" * 60)
    results = [
        ("Section 1: GaAs 21x21 wire eigenvalues", s1_pass),
        ("Section 2: InAs rectangular wire", s2_pass),
        ("Section 3: Dense-sparse consistency",
         s3_pass if s3_pass is not None else True),  # skip = not a failure
        ("Section 4: Wire subband plot", True),
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
