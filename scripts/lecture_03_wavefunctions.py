#!/usr/bin/env python3
"""Lecture 03: Wavefunctions -- Band character and probability density.

Validates eigenvector decomposition and normalization from the 8-band k.p
solver for a GaAs/AlGaAs quantum well at k=0:

  1. CB ground state: verify CB character > 90%
  2. HH ground state: verify HH character > 90%
  3. Normalization: Euclidean sum |v_i|^2 = 1 (LAPACK convention)
  4. Wavefunction plot with per-band decomposition
"""
import glob
import os
import shutil
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
# Configuration
# ---------------------------------------------------------------------------
NUMCB = 4   # CB eigenvalues from config (2 spin-degenerate subbands)
NUMVB = 8   # VB eigenvalues from config (4 spin-degenerate subbands)
NUM_TOTAL = NUMCB + NUMVB  # 12 eigenvalues per k-point
TOL_NORM = 1e-5     # normalization tolerance (g14.6 output limits precision)
TOL_CHARACTER = 0.90  # minimum band character fraction (90%)

# Band labels: 1=HH, 2=LH, 3=LH, 4=HH, 5=SO, 6=SO, 7=CB, 8=CB
BAND_NAMES = ["HH", "LH", "LH", "HH", "SO", "SO", "CB", "CB"]

# Group labels for composite character
HH_BANDS = [0, 3]   # indices 0-based: bands 1 and 4
LH_BANDS = [1, 2]   # bands 2 and 3
SO_BANDS = [4, 5]   # bands 5 and 6
CB_BANDS = [6, 7]   # bands 7 and 8

BAND_COLORS = {
    "HH": "#1f77b4",   # blue
    "LH": "#ff7f0e",   # orange
    "SO": "#2ca02c",   # green
    "CB": "#d62728",   # red
}


def run_bandstructure(config_name, work_dir, timeout=300):
    """Run bandStructure with the given config, return output_dir.

    Args:
        config_name: filename in tests/regression/configs/
        work_dir: temporary working directory
        timeout: execution timeout in seconds

    Returns:
        output_dir path
    """
    config_path = CONFIGS_DIR / config_name
    rc, output_dir = run_exe(str(BUILD_DIR), "bandStructure",
                             str(config_path), work_dir, timeout=timeout)
    if rc != 0:
        sys.exit(f"ERROR: bandStructure returned {rc} for {config_name}")
    return output_dir


def parse_eigenfunction_file(filepath):
    """Parse a QW eigenfunction file into z and per-band |psi| arrays.

    File format: fdstep rows, 9 columns:
        z  |psi_band1|  |psi_band2|  ...  |psi_band8|

    Args:
        filepath: path to eigenfunctions_k_XXXXX_ev_XXXXX.dat

    Returns:
        (z, psi_abs) where z is (fdstep,) and psi_abs is (fdstep, 8)
    """
    data = np.loadtxt(filepath, comments="#")
    z = data[:, 0]
    psi_abs = data[:, 1:]  # (fdstep, 8)
    return z, psi_abs


def parse_parts_file(filepath):
    """Parse parts.dat into a (num_evstates, 8) array of band fractions.

    File format: one row per eigenstate, 8 columns (band fractions).
    """
    data = np.loadtxt(filepath, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def compute_band_character(psi_abs):
    """Compute band character from |psi_band(z)| via Euclidean norm.

    For each band b:
        char_b = sum_z |psi_b(z)|^2

    Returns normalized fractions that sum to 1.
    The eigenvectors from LAPACK zheevx satisfy sum_i |v_i|^2 = 1.

    Args:
        psi_abs: (fdstep, 8) array of |psi_band(z)|

    Returns:
        (8,) array of band character fractions
    """
    char = np.sum(psi_abs ** 2, axis=0)  # (8,)
    total = np.sum(char)
    if total > 0:
        char = char / total
    return char


def group_character(char, band_indices):
    """Sum band character for a group of bands."""
    return sum(char[i] for i in band_indices)


def check_normalization(psi_abs):
    """Verify Euclidean norm sum_i |v_i|^2 = 1.

    LAPACK zheevx returns eigenvectors normalized so that
    sum over all 8*N components |v_i|^2 = 1.

    Args:
        psi_abs: (fdstep, 8) array of |psi_band(z)|

    Returns:
        (norm_euclidean, passed)
    """
    norm = np.sum(psi_abs ** 2)
    passed = abs(norm - 1.0) < TOL_NORM
    return norm, passed


# =========================================================================
# Section 1: CB ground state band character
# =========================================================================
def test_cb_ground_state(output_dir):
    """Verify CB ground state has dominant CB character > 90%."""
    print("=" * 60)
    print("Lecture 03 -- Section 1: CB ground state band character")
    print("=" * 60)

    all_pass = True

    # CB ground state is eigenstate index NUMVB + 1 = 9 (1-based)
    # In the eigenfunction files: ev_00009.dat
    cb_ev_idx = NUMVB + 1
    eigf_path = os.path.join(
        output_dir, f"eigenfunctions_k_00001_ev_{cb_ev_idx:05d}.dat"
    )
    if not os.path.isfile(eigf_path):
        sys.exit(f"ERROR: CB ground state eigenfunction not found: {eigf_path}")

    z, psi_abs = parse_eigenfunction_file(eigf_path)
    char = compute_band_character(psi_abs)

    cb_frac = group_character(char, CB_BANDS)
    hh_frac = group_character(char, HH_BANDS)
    lh_frac = group_character(char, LH_BANDS)
    so_frac = group_character(char, SO_BANDS)

    print(f"  CB ground state (ev {cb_ev_idx}) band decomposition:")
    print(f"    HH: {hh_frac:.6f}  ({hh_frac*100:.2f}%)")
    print(f"    LH: {lh_frac:.6f}  ({lh_frac*100:.2f}%)")
    print(f"    SO: {so_frac:.6f}  ({so_frac*100:.2f}%)")
    print(f"    CB: {cb_frac:.6f}  ({cb_frac*100:.2f}%)")

    if cb_frac > TOL_CHARACTER:
        print(f"  [PASS] CB character = {cb_frac*100:.2f}% > {TOL_CHARACTER*100:.0f}%")
    else:
        print(f"  [FAIL] CB character = {cb_frac*100:.2f}% <= {TOL_CHARACTER*100:.0f}%")
        all_pass = False

    # Also check normalization
    norm, norm_pass = check_normalization(psi_abs)
    if norm_pass:
        print(f"  [PASS] Normalization = {norm:.12f} (tol {TOL_NORM:.0e})")
    else:
        print(f"  [FAIL] Normalization = {norm:.12f} (tol {TOL_NORM:.0e})")
        all_pass = False

    print()
    return all_pass, z, psi_abs, char


# =========================================================================
# Section 2: HH ground state band character
# =========================================================================
def test_hh_ground_state(output_dir):
    """Verify HH ground state (top VB) has dominant HH character > 90%."""
    print("=" * 60)
    print("Lecture 03 -- Section 2: HH ground state band character")
    print("=" * 60)

    all_pass = True

    # The VB eigenvalues are ordered ascending; the highest VB state (least
    # negative = top of VB, HH1) is the LAST VB eigenstate = ev_00008.
    # In QW mode, il selects the top numvb valence states in ascending order,
    # so ev_00001 = deepest VB, ev_00008 = top VB (HH1 ground state).
    hh_ev_idx = NUMVB  # ev_00008
    eigf_path = os.path.join(
        output_dir, f"eigenfunctions_k_00001_ev_{hh_ev_idx:05d}.dat"
    )
    if not os.path.isfile(eigf_path):
        sys.exit(f"ERROR: HH ground state eigenfunction not found: {eigf_path}")

    z, psi_abs = parse_eigenfunction_file(eigf_path)
    char = compute_band_character(psi_abs)

    cb_frac = group_character(char, CB_BANDS)
    hh_frac = group_character(char, HH_BANDS)
    lh_frac = group_character(char, LH_BANDS)
    so_frac = group_character(char, SO_BANDS)

    print(f"  HH ground state (ev {hh_ev_idx}) band decomposition:")
    print(f"    HH: {hh_frac:.6f}  ({hh_frac*100:.2f}%)")
    print(f"    LH: {lh_frac:.6f}  ({lh_frac*100:.2f}%)")
    print(f"    SO: {so_frac:.6f}  ({so_frac*100:.2f}%)")
    print(f"    CB: {cb_frac:.6f}  ({cb_frac*100:.2f}%)")

    if hh_frac > TOL_CHARACTER:
        print(f"  [PASS] HH character = {hh_frac*100:.2f}% > {TOL_CHARACTER*100:.0f}%")
    else:
        print(f"  [FAIL] HH character = {hh_frac*100:.2f}% <= {TOL_CHARACTER*100:.0f}%")
        all_pass = False

    # Also check normalization
    norm, norm_pass = check_normalization(psi_abs)
    if norm_pass:
        print(f"  [PASS] Normalization = {norm:.12f} (tol {TOL_NORM:.0e})")
    else:
        print(f"  [FAIL] Normalization = {norm:.12f} (tol {TOL_NORM:.0e})")
        all_pass = False

    print()
    return all_pass, z, psi_abs, char


# =========================================================================
# Section 3: Wavefunction plot with per-band decomposition
# =========================================================================
def plot_wavefunctions(output_dir):
    """Generate wavefunction plot with per-band probability density.

    Shows CB ground state and HH ground state side by side,
    with stacked per-band |psi_b(z)|^2 decomposition.
    """
    print("=" * 60)
    print("Lecture 03 -- Section 3: Wavefunction visualization")
    print("=" * 60)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    states = [
        ("CB ground state", NUMVB + 1, "CB"),
        ("HH ground state", NUMVB, "HH"),
    ]

    for ax, (label, ev_idx, dominant) in zip(axes, states):
        eigf_path = os.path.join(
            output_dir, f"eigenfunctions_k_00001_ev_{ev_idx:05d}.dat"
        )
        if not os.path.isfile(eigf_path):
            print(f"  [WARN] eigenfunction file not found for {label}")
            continue

        z, psi_abs = parse_eigenfunction_file(eigf_path)
        char = compute_band_character(psi_abs)

        # Group bands for plotting: HH, LH, SO, CB
        psi_sq = psi_abs ** 2  # (fdstep, 8)
        grouped = {
            "HH": psi_sq[:, HH_BANDS].sum(axis=1),
            "LH": psi_sq[:, LH_BANDS].sum(axis=1),
            "SO": psi_sq[:, SO_BANDS].sum(axis=1),
            "CB": psi_sq[:, CB_BANDS].sum(axis=1),
        }

        # Total |psi(z)|^2
        total_psi_sq = psi_sq.sum(axis=1)

        # Plot total envelope
        ax.plot(z, total_psi_sq, "k-", linewidth=1.5, label=r"$|\psi(z)|^2$")

        # Plot per-band contributions as filled regions
        prev = np.zeros_like(z)
        for band_name in ["HH", "LH", "SO", "CB"]:
            contribution = grouped[band_name]
            ax.fill_between(z, prev, prev + contribution,
                            alpha=0.4, color=BAND_COLORS[band_name],
                            label=f"{band_name} ({char_sum(char, band_name)*100:.1f}%)")
            prev = prev + contribution

        # Annotate
        cb_c = group_character(char, CB_BANDS)
        hh_c = group_character(char, HH_BANDS)
        info = f"HH: {hh_c*100:.1f}%  CB: {cb_c*100:.1f}%"
        ax.text(0.02, 0.96, info, transform=ax.transAxes,
                fontsize=9, verticalalignment="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.7))

        ax.set_xlabel(r"$z$ ($\AA$)", fontsize=12)
        ax.set_ylabel(r"$|\psi(z)|^2$ ($\AA^{-1}$)", fontsize=12)
        ax.set_title(f"{label} (ev {ev_idx})", fontsize=12)
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(True, alpha=0.3)

        # Print band fractions for log
        print(f"  {label} (ev {ev_idx}):")
        for bn in ["HH", "LH", "SO", "CB"]:
            frac = char_sum(char, bn)
            print(f"    {bn}: {frac*100:.2f}%")

    fig.suptitle("Lecture 03: GaAs/Al$_{0.3}$Ga$_{0.7}$As QW Wavefunctions at $k=0$",
                 fontsize=14, y=1.02)
    fig.tight_layout()
    out_path = FIGURES_DIR / "lecture_03_wavefunctions.png"
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"  Plot saved to {out_path}")
    print()


def char_sum(char, band_name):
    """Sum band character fractions for a named group."""
    mapping = {
        "HH": HH_BANDS,
        "LH": LH_BANDS,
        "SO": SO_BANDS,
        "CB": CB_BANDS,
    }
    return group_character(char, mapping[band_name])


# =========================================================================
# Main
# =========================================================================
def main():
    print("\n" + "=" * 60)
    print("  LECTURE 03: Wavefunction Decomposition & Normalization")
    print("  GaAs/AlGaAs QW -- band character and probability density")
    print("=" * 60 + "\n")

    # Run bandStructure once in a temporary directory
    with tempfile.TemporaryDirectory() as work_dir:
        output_dir = run_bandstructure("qw_gaas_algaas.toml", work_dir)

        # Verify key output files exist
        eigf_dir = output_dir
        eigf_pattern = os.path.join(eigf_dir, "eigenfunctions_k_00001_ev_*.dat")
        eigf_files = sorted(glob.glob(eigf_pattern))
        if not eigf_files:
            sys.exit(f"ERROR: no eigenfunction files found in {eigf_dir}")
        print(f"  Found {len(eigf_files)} eigenfunction files at k=1")
        print()

        # Run sections
        s1_pass, _, _, _ = test_cb_ground_state(output_dir)
        s2_pass, _, _, _ = test_hh_ground_state(output_dir)
        plot_wavefunctions(output_dir)

    # Summary
    print("=" * 60)
    print("  SUMMARY")
    print("=" * 60)
    results = [
        ("Section 1: CB ground state band character", s1_pass),
        ("Section 2: HH ground state band character", s2_pass),
        ("Section 3: Wavefunction plot", True),
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
