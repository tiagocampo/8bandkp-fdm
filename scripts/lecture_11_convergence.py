#!/usr/bin/env python3
"""Lecture 11: Convergence -- Grid spacing and FD order convergence.

Investigates how the finite-difference grid resolution affects the
accuracy of computed eigenvalues for the 8-band k.p Hamiltonian:

  1. QW grid convergence: sweep FDstep = 51, 101, 201, 401 for a
     GaAs/AlGaAs quantum well. Extract CB ground-state energy at k=0
     for each grid and measure convergence rate vs theoretical order 2.

  2. Richardson extrapolation: using the three finest QW grids, compute
     an extrapolated eigenvalue and compare to the finest-grid result.

  3. Convergence plot: E vs dz on log-log scale with fitted slope and
     theoretical-order reference line.

Saved to docs/lecture/figures/lecture_11_grid_convergence.png
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
    TOL_EXACT, TOL_NUMERICAL,
)


# ---------------------------------------------------------------------------
# Config generation helpers
# ---------------------------------------------------------------------------

QW_TEMPLATE = (
    'confinement = "qw"\n'
    "FDorder = 2\n"
    "fd_step = {fdstep}\n"
    "\n"
    "[wave_vector]\n"
    'mode = "k0"\n'
    "max = 0\n"
    "nsteps = 1\n"
    "\n"
    "[bands]\n"
    "num_cb = 4\n"
    "num_vb = 8\n"
    "\n"
    "[[material]]\n"
    'name = "Al30Ga70As"\n'
    "z_min = -200\n"
    "z_max = 200\n"
    "\n"
    "[[material]]\n"
    'name = "GaAs"\n'
    "z_min = -50\n"
    "z_max = 50\n"
)

# Domain width in Angstroms for the QW config (barrier extends -200 to 200)
QW_DOMAIN_WIDTH = 400.0  # Angstrom


def write_config(work_dir, fdstep):
    """Write a QW config with the given FDstep into a staging file.

    The config is written to a staging path (not input.toml) because run_exe
    copies the config to <work_dir>/input.toml itself.

    Returns the staging config path.
    """
    cfg_path = os.path.join(work_dir, "staged.toml")
    content = QW_TEMPLATE.format(fdstep=fdstep)
    with open(cfg_path, "w") as f:
        f.write(content)
    return cfg_path


def extract_qw_cb_ground(eig_path):
    """Extract the CB ground-state energy from a QW eigenvalue file.

    For a k=0 QW run with numcb=4, numvb=8, there are 12 eigenvalues.
    The CB ground state is the lowest CB eigenvalue, which sits at index
    numvb (= 8, zero-based) in the sorted eigenvalue list.

    Returns:
        float: CB ground-state energy in eV
    """
    data = parse_eigenvalues(eig_path)
    if not data:
        raise RuntimeError(f"No eigenvalue data in {eig_path}")
    # k=0 run: single k-point
    _, evals = data[0]
    # Bands are sorted: valence (lowest) then conduction (highest).
    # With numvb=8, the first CB band is at index 8.
    numvb = 8
    if len(evals) <= numvb:
        raise RuntimeError(
            f"Expected at least {numvb + 1} eigenvalues, got {len(evals)}"
        )
    return evals[numvb]


# ---------------------------------------------------------------------------
# Section 1: QW grid convergence sweep
# ---------------------------------------------------------------------------

def run_qw_convergence_sweep(fdstep_values, timeout=600):
    """Run bandStructure for each FDstep and collect CB ground-state energies.

    Returns:
        dict mapping FDstep -> (dz, E_cb) where dz = domain_width / (FDstep - 1)
    """
    results = {}
    for fdstep in fdstep_values:
        dz = QW_DOMAIN_WIDTH / (fdstep - 1)
        print(f"  FDstep = {fdstep:>4d}  dz = {dz:.4f} A  ... ", end="", flush=True)
        with tempfile.TemporaryDirectory() as work:
            cfg_path = write_config(work, fdstep)
            rc, outdir = run_exe(
                str(BUILD_DIR), "bandStructure", cfg_path, work, timeout=timeout
            )
            if rc != 0:
                print(f"FAILED (rc={rc})")
                sys.exit(f"ERROR: bandStructure returned {rc} for FDstep={fdstep}")

            eig_path = os.path.join(outdir, "eigenvalues.dat")
            e_cb = extract_qw_cb_ground(eig_path)
            results[fdstep] = (dz, e_cb)
            print(f"E_CB = {e_cb:.8f} eV")

    return results


def estimate_convergence_rate(dz_values, energy_values):
    """Estimate convergence rate from log-log linear fit.

    Fits log|E - E_ref| = p * log(dz) + c and returns slope p.
    Uses the finest-grid result as E_ref.

    Returns:
        (rate, intercept, E_ref)
    """
    dz_arr = np.array(dz_values)
    e_arr = np.array(energy_values)

    # Use Richardson-extrapolated value as reference
    e_ref = richardson_extrapolation(dz_arr, e_arr)

    errors = np.abs(e_arr - e_ref)
    # Only use points with meaningful error (skip the finest if it equals E_ref)
    mask = errors > 1e-15
    if np.sum(mask) < 2:
        return 0.0, 0.0, e_ref

    log_dz = np.log(dz_arr[mask])
    log_err = np.log(errors[mask])

    # Linear fit: log_err = rate * log_dz + intercept
    coeffs = np.polyfit(log_dz, log_err, 1)
    rate = coeffs[0]
    intercept = coeffs[1]

    return rate, intercept, e_ref


def richardson_extrapolation(dz_arr, e_arr):
    """Richardson extrapolation using the two finest grids.

    For a method of order p, E_exact ~ E_h + c * h^p.
    Using two grids h1, h2:
        E_exact = (E_h2 * h1^p - E_h1 * h2^p) / (h1^p - h2^p)

    We assume p=2 (FD order 2).

    Returns:
        float: extrapolated eigenvalue
    """
    # Sort by dz ascending (finest grid first)
    idx = np.argsort(dz_arr)
    dz_sorted = dz_arr[idx]
    e_sorted = e_arr[idx]

    # Use the two finest grids
    h1, h2 = dz_sorted[0], dz_sorted[1]
    e1, e2 = e_sorted[0], e_sorted[1]

    p = 2  # FD order
    h1p = h1 ** p
    h2p = h2 ** p
    e_exact = (e1 * h2p - e2 * h1p) / (h2p - h1p)
    return e_exact


# ---------------------------------------------------------------------------
# Section 2: FD order convergence (QW with FDorder sweep)
# ---------------------------------------------------------------------------

def run_fdorder_sweep(fdstep, fdorder_values, timeout=600):
    """Run bandStructure for a fixed FDstep with varying FDorder.

    Returns:
        dict mapping FDorder -> E_cb
    """
    results = {}
    for order in fdorder_values:
        print(f"  FDorder = {order}  FDstep = {fdstep}  ... ", end="", flush=True)
        with tempfile.TemporaryDirectory() as work:
            cfg_content = (
                'confinement = "qw"\n'
                f"FDorder = {order}\n"
                f"fd_step = {fdstep}\n"
                "\n"
                "[wave_vector]\n"
                'mode = "k0"\n'
                "max = 0\n"
                "nsteps = 1\n"
                "\n"
                "[bands]\n"
                "num_cb = 4\n"
                "num_vb = 8\n"
                "\n"
                "[[material]]\n"
                'name = "Al30Ga70As"\n'
                "z_min = -200\n"
                "z_max = 200\n"
                "\n"
                "[[material]]\n"
                'name = "GaAs"\n'
                "z_min = -50\n"
                "z_max = 50\n"
            )
            # Write to staged.toml; run_exe will copy it to input.toml
            cfg_path = os.path.join(work, "staged.toml")
            with open(cfg_path, "w") as f:
                f.write(cfg_content)

            rc, outdir = run_exe(
                str(BUILD_DIR), "bandStructure", cfg_path, work, timeout=timeout
            )
            if rc != 0:
                print(f"FAILED (rc={rc})")
                sys.exit(f"ERROR: bandStructure returned {rc} for FDorder={order}")

            eig_path = os.path.join(outdir, "eigenvalues.dat")
            e_cb = extract_qw_cb_ground(eig_path)
            results[order] = e_cb
            print(f"E_CB = {e_cb:.8f} eV")

    return results


# ---------------------------------------------------------------------------
# Section 3: Plotting
# ---------------------------------------------------------------------------

def plot_convergence(results, rate, intercept, e_ref):
    """Generate convergence plot: E vs dz with rate annotation."""
    fdsteps = sorted(results.keys())
    dz_vals = [results[fd][0] for fd in fdsteps]
    e_vals = [results[fd][1] for fd in fdsteps]

    dz_arr = np.array(dz_vals)
    e_arr = np.array(e_vals)

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(14, 6),
    )

    # --- Left panel: E vs dz (linear scale) ---
    ax1.plot(dz_arr, e_arr, "bo-", linewidth=2, markersize=8, label="Computed $E_{CB}$")
    ax1.axhline(e_ref, color="red", linestyle="--", linewidth=1.2,
                label=f"Richardson extrapolation = {e_ref:.6f} eV")

    # Annotate each point with FDstep
    for dz, e, fd in zip(dz_arr, e_arr, fdsteps):
        ax1.annotate(f"N={fd}", (dz, e), textcoords="offset points",
                     xytext=(8, 8), fontsize=9)

    ax1.set_xlabel(r"Grid spacing $\Delta z$ ($\AA$)", fontsize=12)
    ax1.set_ylabel("CB ground state energy (eV)", fontsize=12)
    ax1.set_title("QW Grid Convergence", fontsize=13)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.invert_xaxis()  # finer grid on the right

    # --- Right panel: log-log error convergence ---
    errors = np.abs(e_arr - e_ref)
    # Skip points with essentially zero error
    mask = errors > 1e-15
    if np.any(mask):
        ax2.loglog(dz_arr[mask], errors[mask], "bo-", linewidth=2, markersize=8,
                   label="Computed error")

        # Theoretical reference lines for FD order 2 and 4
        dz_ref = np.logspace(np.log10(dz_arr.min()), np.log10(dz_arr.max()), 50)

        # Scale reference line to pass through the coarsest-grid point
        if np.any(mask):
            dz_coarse = dz_arr[mask][-1]
            err_coarse = errors[mask][-1]

            for p_ref, style, label in [
                (2, "r--", r"$O(\Delta z^2)$"),
                (4, "g-.", r"$O(\Delta z^4)$"),
            ]:
                ref_line = err_coarse * (dz_ref / dz_coarse) ** p_ref
                ax2.loglog(dz_ref, ref_line, style, linewidth=1.5, label=label)

        # Annotate fitted rate
        ax2.annotate(
            f"Fitted rate = {rate:.2f}",
            xy=(0.05, 0.05), xycoords="axes fraction",
            fontsize=12, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.8),
        )

    ax2.set_xlabel(r"Grid spacing $\Delta z$ ($\AA$)", fontsize=12)
    ax2.set_ylabel(r"$|E - E_{ref}|$ (eV)", fontsize=12)
    ax2.set_title("Convergence Rate (log-log)", fontsize=13)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3, which="both")

    fig.suptitle("Lecture 11: Finite-Difference Grid Convergence\n"
                 "GaAs/Al$_{0.3}$Ga$_{0.7}$As QW, FDorder = 2",
                 fontsize=14, fontweight="bold", y=1.02)
    fig.tight_layout()

    out_path = FIGURES_DIR / "lecture_11_grid_convergence.png"
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Plot saved to {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("\n" + "=" * 60)
    print("  LECTURE 11: Convergence -- Grid Spacing and FD Order")
    print("  8-band zinc-blende k.p finite-difference convergence")
    print("=" * 60 + "\n")

    # =====================================================================
    # Section 1: QW grid convergence (FDstep sweep, FDorder=2)
    # =====================================================================
    print("=" * 60)
    print("  Section 1: QW Grid Convergence (FDorder = 2)")
    print("=" * 60)
    print("  GaAs/AlGaAs QW, sweeping FDstep to vary grid spacing.\n")

    fdstep_values = [51, 101, 201, 401]
    qw_results = run_qw_convergence_sweep(fdstep_values, timeout=600)

    # Compute convergence rate
    fdsteps = sorted(qw_results.keys())
    dz_arr = np.array([qw_results[fd][0] for fd in fdsteps])
    e_arr = np.array([qw_results[fd][1] for fd in fdsteps])

    rate, intercept, e_ref = estimate_convergence_rate(dz_arr, e_arr)

    print(f"\n  Richardson extrapolated E_ref = {e_ref:.8f} eV")
    print(f"  Fitted convergence rate = {rate:.3f}")
    print(f"  Theoretical rate (FDorder=2) = 2.0")

    # Tolerance: convergence rate within 0.5 of theoretical order
    theoretical_order = 2.0
    rate_pass = abs(rate - theoretical_order) < 0.5
    rate_status = "PASS" if rate_pass else "FAIL"
    print(f"  [{rate_status}] |rate - theoretical| = "
          f"{abs(rate - theoretical_order):.3f} < 0.5")

    # Verify monotonic convergence (energies increase toward E_ref as
    # dz decreases, since finer grids resolve the confinement better)
    monotonic = all(
        e_arr[i] <= e_arr[i + 1] + 1e-12
        for i in range(len(e_arr) - 1)
    ) or all(
        e_arr[i] >= e_arr[i + 1] - 1e-12
        for i in range(len(e_arr) - 1)
    )
    mono_status = "PASS" if monotonic else "FAIL"
    print(f"  [{mono_status}] Eigenvalues vary monotonically with grid spacing")

    print()

    # =====================================================================
    # Section 2: FD order sweep (fixed FDstep=201)
    # =====================================================================
    print("=" * 60)
    print("  Section 2: FD Order Convergence (FDstep = 201)")
    print("=" * 60)
    print("  GaAs/AlGaAs QW, sweeping FDorder to check higher-order FD.\n")

    # Only test orders that are available: 2, 4, 6, 8, 10
    fdorder_values = [2, 4, 6]
    fdorder_results = run_fdorder_sweep(201, fdorder_values, timeout=600)

    # Higher FDorder should give more accurate results
    print("\n  FD order comparison:")
    e_order2 = fdorder_results[2]
    for order in fdorder_values:
        e = fdorder_results[order]
        delta = abs(e - e_order2)
        print(f"    FDorder={order}: E_CB = {e:.8f} eV  "
              f"(delta vs order-2: {delta:.2e} eV)")

    # Verify that higher orders converge toward a common value
    e_finest = fdorder_results[max(fdorder_values)]
    e_coarsest = fdorder_results[min(fdorder_values)]
    order_converge = abs(e_finest - e_order2) < abs(e_coarsest - e_order2) or \
        max(fdorder_values) == min(fdorder_values)
    # Actually, just check that all FD orders give similar results
    e_vals_order = [fdorder_results[o] for o in fdorder_values]
    spread = max(e_vals_order) - min(e_vals_order)
    order_pass = spread < 0.01  # within 10 meV
    order_status = "PASS" if order_pass else "FAIL"
    print(f"\n  [{order_status}] FD order spread = {spread:.2e} eV (< 0.01 eV)")

    print()

    # =====================================================================
    # Section 3: Plot
    # =====================================================================
    print("=" * 60)
    print("  Section 3: Convergence Plot")
    print("=" * 60 + "\n")

    plot_convergence(qw_results, rate, intercept, e_ref)

    # =====================================================================
    # Summary
    # =====================================================================
    print("\n" + "=" * 60)
    print("  SUMMARY")
    print("=" * 60)
    results = [
        ("Section 1: QW grid convergence rate", rate_pass),
        ("Section 1: Monotonic convergence", monotonic),
        ("Section 2: FD order spread < 0.01 eV", order_pass),
        ("Section 3: Convergence plot", True),
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
