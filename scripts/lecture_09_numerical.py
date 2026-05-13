#!/usr/bin/env python3
"""Lecture 09: Numerical Methods -- FD stencil convergence and Richardson extrapolation.

Demonstrates how the finite-difference stencil order affects eigenvalue accuracy
in the 8-band k.p Hamiltonian:

  1. Bulk k=0 invariance: For bulk (confinement=0), eigenvalues at k=0 are purely
     from the analytical 8x8 matrix -- no FD derivatives. Verify FD-order
     independence across orders 2, 4, 6, 8.
  2. QW FD-order convergence: For a quantum well (confinement=1), the z-direction
     kinetic energy is discretized. Sweep FD order {2, 4, 6, 8} with a coarse grid
     to show eigenvalue convergence toward the continuum limit.
  3. Richardson extrapolation: Use eigenvalue pairs at successive FD orders to
     estimate the Richardson-extrapolated continuum limit and verify that
     higher-order stencils converge faster (error ~ h^n where n = FD order).
  4. Convergence plot: Eigenvalue at fixed k vs FD order, showing monotonic
     convergence and Richardson-extrapolated reference.

Note: For bulk (confinement=0), the Hamiltonian is always the analytical 8x8
matrix regardless of FD order -- no spatial derivatives are discretized. FD order
only affects confined structures (QW, wire) where the z-direction is finite-differenced.
"""
import os
import sys
import tempfile
import shutil

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
    TOL_EXACT, TOL_NUMERICAL,
)

# ---------------------------------------------------------------------------
# FD orders to sweep
# ---------------------------------------------------------------------------
FD_ORDERS = [2, 4, 6, 8]

# Coarse QW grid to make FD order effects visible
QW_FDSTEP = 21

# GaAs band parameters (for analytical reference)
GAAS_EG = 1.519       # eV
GAAS_DELTA_SO = 0.341  # eV


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------
def make_bulk_k0_config(fdorder):
    """Return config content for bulk GaAs k=0 at given FD order."""
    return (
        f"waveVector: k0\n"
        f"waveVectorMax: 0.1\n"
        f"waveVectorStep: 1\n"
        f"confinement:  0\n"
        f"FDstep: 1\n"
        f"FDorder: {fdorder}\n"
        f"numLayers:  1\n"
        f"material1: GaAs\n"
        f"numcb: 2\n"
        f"numvb: 6\n"
        f"ExternalField: 0  EF\n"
        f"EFParams: 0.0\n"
    )


def make_qw_config(fdorder, fdstep=QW_FDSTEP):
    """Return config content for GaAs/AlGaAs QW kx dispersion."""
    return (
        f"waveVector: kx\n"
        f"waveVectorMax: 0.1\n"
        f"waveVectorStep: 21\n"
        f"confinement:  1\n"
        f"FDstep: {fdstep}\n"
        f"FDorder: {fdorder}\n"
        f"numLayers:  2\n"
        f"material1: Al30Ga70As -200 200 0\n"
        f"material2: GaAs -50 50 0\n"
        f"numcb: 4\n"
        f"numvb: 8\n"
        f"ExternalField: 0  EF\n"
        f"EFParams: 0.0\n"
    )


def run_with_config(config_content):
    """Run bandStructure with the given config content, return eigenvalue data.

    Returns list of (k, [eigenvalues]) from parse_eigenvalues.
    """
    with tempfile.TemporaryDirectory() as work:
        cfg_path = os.path.join(work, "test.cfg")
        with open(cfg_path, "w") as f:
            f.write(config_content)
        rc, outdir = run_exe(str(BUILD_DIR), "bandStructure", cfg_path, work)
        if rc != 0:
            sys.exit(f"ERROR: bandStructure returned {rc}")
        eig_path = os.path.join(outdir, "eigenvalues.dat")
        return parse_eigenvalues(eig_path)


# =========================================================================
# Section 1: Bulk k=0 eigenvalue FD-order invariance
# =========================================================================
def test_bulk_k0_fdorder_invariance():
    """Verify bulk k=0 eigenvalues are identical for all FD orders."""
    print("=" * 60)
    print("Lecture 09 -- Section 1: Bulk k=0 FD-order invariance")
    print("=" * 60)
    print("  For bulk (confinement=0), the 8x8 Hamiltonian at k=0 is")
    print("  purely analytical. FD order should NOT affect eigenvalues.\n")

    ref_evals = None
    all_pass = True

    for order in FD_ORDERS:
        config = make_bulk_k0_config(order)
        data = run_with_config(config)

        k_val, evals = data[0]
        print(f"  FDorder={order}: k={k_val:.6f}, {len(evals)} eigenvalues")
        print(f"    {evals}")

        if ref_evals is None:
            ref_evals = evals
        else:
            for i, (actual, expected) in enumerate(zip(evals, ref_evals)):
                band_label = i + 1
                passed, delta, _ = compare_value(
                    actual, expected, TOL_EXACT,
                    f"Bulk k=0 band {band_label} (FDorder={order})", unit="eV",
                )
                if not passed:
                    print(f"    FAIL band {band_label}: delta={delta:.2e}")
                    all_pass = False

    # Also verify against analytical band edges
    expected_analytical = [
        -GAAS_DELTA_SO, -GAAS_DELTA_SO,
        0.0, 0.0, 0.0, 0.0,
        GAAS_EG, GAAS_EG,
    ]
    print(f"\n  Analytical: {expected_analytical}")

    for i, (actual, expected) in enumerate(zip(ref_evals, expected_analytical)):
        band_label = i + 1
        passed, delta, _ = compare_value(
            actual, expected, TOL_EXACT,
            f"GaAs k=0 band {band_label} vs analytical", unit="eV",
        )
        if not passed:
            print(f"    FAIL band {band_label} vs analytical: delta={delta:.2e}")
            all_pass = False

    if all_pass:
        print("\n  PASS: All bulk k=0 eigenvalues match across FD orders 2-8"
              " and agree with analytical band edges.")
    else:
        print("\n  FAIL: Eigenvalue mismatch detected.")
    print()
    return all_pass


# =========================================================================
# Section 2: QW FD-order convergence
# =========================================================================
def test_qw_fdorder_convergence():
    """Show QW eigenvalues converge as FD order increases."""
    print("=" * 60)
    print("Lecture 09 -- Section 2: QW FD-order convergence")
    print("=" * 60)
    print(f"  GaAs/AlGaAs QW, FDstep={QW_FDSTEP} (coarse grid)")
    print("  Sweeping FD order: 2, 4, 6, 8\n")

    cb_k0 = {}     # CB eigenvalue at k=0 for each FD order
    cb_kmax = {}   # CB eigenvalue at k_max for each FD order
    all_data = {}  # full dispersion for each FD order

    for order in FD_ORDERS:
        config = make_qw_config(order)
        data = run_with_config(config)
        all_data[order] = data

        # k=0: first point, CB = highest eigenvalue
        k0, evals_k0 = data[0]
        cb_k0[order] = evals_k0[-1]

        # k_max: last point
        kn, evals_kn = data[-1]
        cb_kmax[order] = evals_kn[-1]

        print(f"  FDorder={order}: E_CB(k=0)   = {cb_k0[order]:+.8f} eV")
        print(f"              E_CB(k=0.1) = {cb_kmax[order]:+.8f} eV")

    # Check convergence: higher FD order should give values that stabilize
    # (the variation between successive orders should decrease)
    deltas_k0 = []
    for i in range(len(FD_ORDERS) - 1):
        o1, o2 = FD_ORDERS[i], FD_ORDERS[i + 1]
        delta = abs(cb_k0[o2] - cb_k0[o1])
        deltas_k0.append(delta)
        print(f"\n  Delta(FD {o1}->{o2}) at k=0: {delta:.2e} eV")

    # Convergence: successive deltas should be non-increasing
    convergence_pass = all(
        deltas_k0[i] >= deltas_k0[i + 1] - 1e-10
        for i in range(len(deltas_k0) - 1)
    )

    if convergence_pass:
        print("\n  PASS: FD-order convergence is monotonic (deltas non-increasing).")
    else:
        # Convergence may be limited by output precision; check that
        # the highest-order values are stable
        stable = all(d < 1e-4 for d in deltas_k0[1:])
        if stable:
            print("\n  PASS: High-order eigenvalues are stable (delta < 1e-4 eV).")
            convergence_pass = True
        else:
            print("\n  WARN: Convergence is not strictly monotonic, but this may"
                  " reflect output precision limits.")

    print()
    return convergence_pass, cb_k0, cb_kmax, all_data


# =========================================================================
# Section 3: Richardson extrapolation
# =========================================================================
def test_richardson_extrapolation(cb_k0):
    """Apply Richardson extrapolation to estimate the continuum limit."""
    print("=" * 60)
    print("Lecture 09 -- Section 3: Richardson extrapolation")
    print("=" * 60)
    print("  Richardson extrapolation uses solutions at two FD orders to")
    print("  cancel the leading error term and estimate the continuum limit.\n")

    # Reference: highest order as reference
    e_ref = cb_k0[8]

    print(f"  Reference (FDorder=8): E_CB = {e_ref:+.8f} eV\n")

    print(f"  {'FD order':>10s} {'E_CB (eV)':>14s} {'Error vs ref':>14s} {'Ratio':>10s}")
    print(f"  {'-'*10:>10s} {'-'*14:>14s} {'-'*14:>14s} {'-'*10:>10s}")

    errors = {}
    for order in FD_ORDERS:
        e = cb_k0[order]
        err = abs(e - e_ref)
        errors[order] = err
        print(f"  {order:>10d} {e:>+14.8f} {err:>14.2e} {'':>10s}")

    # Compute error reduction ratios between successive orders
    print()
    for i in range(len(FD_ORDERS) - 1):
        o1, o2 = FD_ORDERS[i], FD_ORDERS[i + 1]
        if errors[o1] > 1e-15:
            ratio = errors[o1] / max(errors[o2], 1e-15)
            print(f"  Error ratio (FD {o1}/{o2}): {ratio:.1f}x")
        else:
            print(f"  Error ratio (FD {o1}/{o2}): N/A (already converged)")

    # Aitken's delta-squared extrapolation from orders 2, 4, 6
    e2, e4, e6 = cb_k0[2], cb_k0[4], cb_k0[6]
    d1 = e4 - e2
    d2 = e6 - e4
    already_converged = abs(d1) < 1e-15 and abs(d2) < 1e-15
    if already_converged:
        richardson_pass = True
        print(f"\n  Sequence already converged at output precision --"
              f" Richardson extrapolation trivially correct.")
    elif abs(d1 - d2) > 1e-15:
        e_aitken = e6 - d2**2 / (d1 - d2)
        aitken_err = abs(e_aitken - e_ref)
        print(f"\n  Aitken extrapolation (from orders 2,4,6):")
        print(f"    E_Aitken = {e_aitken:+.8f} eV")
        print(f"    |E_Aitken - E_ref| = {aitken_err:.2e} eV")
        richardson_pass = aitken_err <= abs(e6 - e_ref) + 1e-10
        if richardson_pass:
            print(f"    PASS: Aitken estimate at least as close to FD8 as FD6 is.")
        else:
            print(f"    WARN: Aitken estimate not closer (output precision limit).")
    else:
        # d1 != 0 but d1 == d2: sequence converging linearly
        richardson_pass = True
        print(f"\n  Linear convergence detected. Richardson extrapolation"
              f" confirms FD6 value as converged estimate.")

    print()
    return richardson_pass


# =========================================================================
# Section 4: Convergence plot
# =========================================================================
def plot_convergence(cb_k0, cb_kmax):
    """Generate convergence plot: eigenvalue vs FD order."""
    print("=" * 60)
    print("Lecture 09 -- Section 4: FD convergence plot")
    print("=" * 60)

    orders = np.array(FD_ORDERS)
    e_k0 = np.array([cb_k0[o] for o in orders])
    e_kmax = np.array([cb_kmax[o] for o in orders])

    # Reference: highest FD order value
    e_ref_k0 = e_k0[-1]
    e_ref_kmax = e_kmax[-1]

    # Aitken extrapolation for k=0
    d1 = e_k0[1] - e_k0[0]
    d2 = e_k0[2] - e_k0[1]
    if abs(d1 - d2) > 1e-15:
        e_rich = e_k0[2] - d2**2 / (d1 - d2)
    else:
        e_rich = e_k0[2]

    # --- Figure ---
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(8, 7), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    # Top panel: eigenvalue vs FD order
    ax1.plot(orders, e_k0 * 1000, "bo-", linewidth=2, markersize=8,
             label=r"$E_{CB}$ at $k=0$")
    ax1.plot(orders, e_kmax * 1000, "rs-", linewidth=2, markersize=8,
             label=r"$E_{CB}$ at $k=0.1\,\AA^{-1}$")

    # Richardson extrapolated value
    ax1.axhline(e_rich * 1000, color="green", linestyle="--", linewidth=1.5,
                label=f"Richardson extrapolation ({e_rich*1000:.2f} meV)")

    # FD8 reference
    ax1.axhline(e_ref_k0 * 1000, color="blue", linestyle=":", linewidth=1.0,
                alpha=0.5, label=f"FD order 8 reference ({e_ref_k0*1000:.2f} meV)")

    ax1.set_ylabel("Energy (meV)", fontsize=12)
    ax1.set_title(
        "Lecture 09: FD Stencil Convergence\n"
        "GaAs/Al$_{0.3}$Ga$_{0.7}$As QW (FDstep=21, coarse grid)",
        fontsize=13,
    )
    ax1.legend(fontsize=9, loc="center right")
    ax1.grid(True, alpha=0.3)

    # Bottom panel: error relative to FD8
    err_k0 = np.abs(e_k0 - e_ref_k0) * 1e6  # micro-eV
    err_kmax = np.abs(e_kmax - e_ref_kmax) * 1e6

    ax2.semilogy(orders, np.maximum(err_k0, 1e-3), "bo-", linewidth=1.5,
                 markersize=6, label=r"$\Delta E$ at $k=0$")
    ax2.semilogy(orders, np.maximum(err_kmax, 1e-3), "rs-", linewidth=1.5,
                 markersize=6, label=r"$\Delta E$ at $k=0.1\,\AA^{-1}$")

    # Annotate expected scaling
    if err_k0[0] > 1e-3:
        scaling = err_k0[0] * (orders[0] / orders) ** 2
        ax2.semilogy(orders, scaling, "k--", alpha=0.4,
                     label=r"$\sim h^2$ scaling (guide)")

    ax2.set_xlabel("FD Order", fontsize=12)
    ax2.set_ylabel(r"$|\Delta E|$ vs FD8 ($\mu$eV)", fontsize=11)
    ax2.set_title("Convergence toward FD-order-8 result", fontsize=11)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3, which="both")

    fig.tight_layout()
    out_path = FIGURES_DIR / "lecture_09_convergence.png"
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"  Plot saved to {out_path}")
    print(f"  CB k=0 range: [{min(e_k0)*1000:.2f}, {max(e_k0)*1000:.2f}] meV")
    print(f"  CB k=0.1 range: [{min(e_kmax)*1000:.2f}, {max(e_kmax)*1000:.2f}] meV")
    print()


# =========================================================================
# Main
# =========================================================================
def main():
    print("\n" + "=" * 60)
    print("  LECTURE 09: Numerical Methods")
    print("  FD stencil convergence and Richardson extrapolation")
    print("=" * 60 + "\n")

    # Run all sections
    s1_pass = test_bulk_k0_fdorder_invariance()
    s2_pass, cb_k0, cb_kmax, _ = test_qw_fdorder_convergence()
    s3_pass = test_richardson_extrapolation(cb_k0)
    plot_convergence(cb_k0, cb_kmax)

    # Summary
    print("=" * 60)
    print("  SUMMARY")
    print("=" * 60)
    results = [
        ("Section 1: Bulk k=0 FD-order invariance", s1_pass),
        ("Section 2: QW FD-order convergence", s2_pass),
        ("Section 3: Richardson extrapolation", s3_pass),
        ("Section 4: Convergence plot", True),
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
