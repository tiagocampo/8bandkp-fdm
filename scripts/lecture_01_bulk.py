#!/usr/bin/env python3
"""Lecture 01: Bulk Band Structure -- 8-band k.p Hamiltonian validation.

Validates the bulk 8-band zinc-blende Hamiltonian against analytical
results from Kane theory:

  1. k=0 eigenvalue self-check (GaAs): verifies band edge energies
     [-Delta_SO, -Delta_SO, 0, 0, 0, 0, E_g, E_g].
  2. Conduction-band effective mass (GaAs): parabolic fit near k=0
     compared to the 2-band Kane formula m*/m_e = E_g/(E_P + E_g).
  3. InAs band gap: verifies E_g = 0.417 eV at the Gamma point.
  4. Overlay plot: code E(k) vs parabolic Kane dispersion for GaAs CB,
     illustrating non-parabolicity deviation.
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
    extract_effective_mass, TOL_EXACT, TOL_ANALYTICAL, TOL_NUMERICAL,
    HBAR2_OVER_2M0,
)

# ---------------------------------------------------------------------------
# Material parameters
# ---------------------------------------------------------------------------
# GaAs (Vurgaftman 2001 / Winkler 2003)
GAAS_EG = 1.519       # eV, band gap
GAAS_DELTA_SO = 0.341  # eV, spin-orbit splitting
GAAS_EP = 28.8         # eV, Kane interband matrix element

# InAs (Vurgaftman 2001)
INAS_EG = 0.417        # eV, band gap

# 2-band Kane effective mass: m*/m_e = E_g / (E_P + E_g)
GAAS_MSTAR_KANE = GAAS_EG / (GAAS_EP + GAAS_EG)


# =========================================================================
# Section 1: k=0 eigenvalue self-check (GaAs)
# =========================================================================
def test_gaas_k0_eigenvalues():
    """Verify GaAs k=0 eigenvalues match analytical band edges."""
    print("=" * 60)
    print("Lecture 01 -- Section 1: GaAs k=0 eigenvalue self-check")
    print("=" * 60)

    cfg = CONFIGS_DIR / "bulk_gaas_k0.toml"
    with tempfile.TemporaryDirectory() as work:
        rc, outdir = run_exe(str(BUILD_DIR), "bandStructure",
                             str(cfg), work)
        if rc != 0:
            sys.exit(f"ERROR: bandStructure returned {rc} for GaAs k=0")

        eig_path = os.path.join(outdir, "eigenvalues.dat")
        data = parse_eigenvalues(eig_path)
        if not data:
            sys.exit("ERROR: no eigenvalue data parsed from GaAs k=0")

        # k=0 run: single k-point
        k_val, evals = data[0]
        print(f"  k = {k_val:.6f},  {len(evals)} eigenvalues")
        print(f"  Computed: {evals}")

        # Analytical expectations for GaAs at k=0
        # Basis: bands 1-2 = HH (E = -Delta_SO), bands 3-4 = LH (E = 0),
        #         bands 5-6 = SO (E = 0), bands 7-8 = CB (E = E_g)
        expected = [
            -GAAS_DELTA_SO, -GAAS_DELTA_SO,
            0.0, 0.0, 0.0, 0.0,
            GAAS_EG, GAAS_EG,
        ]
        print(f"  Expected: {expected}")

        all_pass = True
        for i, (actual, exp) in enumerate(zip(evals, expected)):
            band_label = i + 1
            passed, delta, _ = compare_value(
                actual, exp, TOL_EXACT,
                f"GaAs k=0 band {band_label}", unit="eV",
            )
            status = "PASS" if passed else "FAIL"
            print(f"    Band {band_label}: {actual:+.8f}  "
                  f"(expected {exp:+.3f}, delta={delta:.2e})  {status}")
            all_pass = all_pass and passed

        if all_pass:
            print("  --> All 8 eigenvalues match analytical band edges.")
        else:
            print("  --> FAIL: one or more eigenvalues outside tolerance.")
        print()
        return all_pass


# =========================================================================
# Section 2: GaAs conduction-band effective mass
# =========================================================================
def test_gaas_effective_mass():
    """Compare extracted CB effective mass with the Kane formula."""
    print("=" * 60)
    print("Lecture 01 -- Section 2: GaAs CB effective mass")
    print("=" * 60)

    cfg = CONFIGS_DIR / "bulk_gaas_kx_dispersion.toml"
    with tempfile.TemporaryDirectory() as work:
        rc, outdir = run_exe(str(BUILD_DIR), "bandStructure",
                             str(cfg), work)
        if rc != 0:
            sys.exit(f"ERROR: bandStructure returned {rc} for GaAs dispersion")

        eig_path = os.path.join(outdir, "eigenvalues.dat")
        result = extract_effective_mass(eig_path, cb_index=-1)
        if result is None:
            sys.exit("ERROR: effective mass extraction failed for GaAs")

        m_star_me, e0, r2, k_max = result
        print(f"  Extracted m* = {m_star_me:.6f} m_e")
        print(f"  Parabolic fit E0 = {e0:.6f} eV,  R^2 = {r2:.8f},  "
              f"k_max = {k_max:.4f} 1/A")
        print(f"  Kane formula m* = E_g/(E_P + E_g) = "
              f"{GAAS_MSTAR_KANE:.6f} m_e")

        # 10% tolerance: 8-band higher-order mixing causes deviation
        # from the simple 2-band Kane formula
        tol_mass = 0.10
        passed, delta, _ = compare_value(
            m_star_me, GAAS_MSTAR_KANE, tol_mass,
            "GaAs CB effective mass", unit="m_e",
        )
        status = "PASS" if passed else "FAIL"
        print(f"  {status}: delta = {delta:.4f}  (tolerance {tol_mass:.0%})")

        if passed:
            print("  --> Effective mass matches Kane formula within 10%.")
        else:
            print("  --> FAIL: effective mass outside 10% tolerance.")
        print()
        return passed, eig_path


# =========================================================================
# Section 3: InAs band gap
# =========================================================================
def test_inas_band_gap():
    """Verify InAs band gap at the Gamma point."""
    print("=" * 60)
    print("Lecture 01 -- Section 3: InAs band gap")
    print("=" * 60)

    cfg = CONFIGS_DIR / "bulk_inas_kx_dispersion.toml"
    with tempfile.TemporaryDirectory() as work:
        rc, outdir = run_exe(str(BUILD_DIR), "bandStructure",
                             str(cfg), work)
        if rc != 0:
            sys.exit(f"ERROR: bandStructure returned {rc} for InAs dispersion")

        eig_path = os.path.join(outdir, "eigenvalues.dat")
        data = parse_eigenvalues(eig_path)
        if not data:
            sys.exit("ERROR: no eigenvalue data parsed for InAs")

        # First k-point is k=0
        k0, evals = data[0]
        cb_min = evals[-1]  # highest eigenvalue = CB edge

        print(f"  InAs CB minimum at k=0: {cb_min:.6f} eV")
        print(f"  Vurgaftman 2001 E_g:     {INAS_EG:.6f} eV")

        passed, delta, _ = compare_value(
            cb_min, INAS_EG, TOL_EXACT,
            "InAs band gap", unit="eV",
        )
        status = "PASS" if passed else "FAIL"
        print(f"  {status}: delta = {delta:.2e}  (tolerance {TOL_EXACT:.0e})")

        if passed:
            print("  --> InAs band gap matches Vurgaftman 2001.")
        else:
            print("  --> FAIL: band gap outside tolerance.")
        print()
        return passed


# =========================================================================
# Section 4: Overlay plot -- code E(k) vs Kane parabolic dispersion
# =========================================================================
def plot_gaas_dispersion(gaas_eig_path):
    """Overlay code CB E(k) with parabolic Kane E(k) for GaAs."""
    print("=" * 60)
    print("Lecture 01 -- Section 4: GaAs CB dispersion overlay plot")
    print("=" * 60)

    data = parse_eigenvalues(gaas_eig_path)
    if not data:
        sys.exit("ERROR: no eigenvalue data for GaAs dispersion plot")

    ks_code = np.array([d[0] for d in data])
    es_cb = np.array([d[1][-1] for d in data])  # CB (highest eigenvalue)

    # Parabolic Kane dispersion: E(k) = E_g + hbar^2 k^2 / (2 m*)
    # where m* = Kane mass
    es_kane = GAAS_EG + HBAR2_OVER_2M0 / GAAS_MSTAR_KANE * ks_code ** 2

    # Non-parabolicity deviation
    deviation = es_cb - es_kane

    # --- Figure ---
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(8, 7), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    # Top panel: E(k) comparison
    ax1.plot(ks_code, es_cb, "b-", linewidth=2, label="8-band k.p (code)")
    ax1.plot(ks_code, es_kane, "r--", linewidth=1.5,
             label=f"Kane parabolic (m* = {GAAS_MSTAR_KANE:.4f} m$_e$)")

    # Shade the tolerance band (10% around Kane parabolic)
    tol_frac = 0.10
    band_upper = es_kane * (1 + tol_frac) if GAAS_EG > 0 else es_kane + tol_frac
    band_lower = es_kane * (1 - tol_frac) if GAAS_EG > 0 else es_kane - tol_frac
    ax1.fill_between(ks_code, band_lower, band_upper,
                     alpha=0.15, color="red", label="10% tolerance band")

    ax1.set_ylabel("Energy (eV)", fontsize=12)
    ax1.set_title("Lecture 01: GaAs Bulk CB Dispersion -- 8-band vs Kane", fontsize=13)
    ax1.legend(fontsize=10, loc="upper left")
    ax1.grid(True, alpha=0.3)

    # Bottom panel: deviation
    ax2.plot(ks_code, deviation, "k-", linewidth=1.2)
    ax2.axhline(0, color="gray", linestyle=":", linewidth=0.8)
    ax2.fill_between(ks_code, 0, deviation, alpha=0.2, color="orange")
    ax2.set_xlabel(r"$k$ (1/$\AA$)", fontsize=12)
    ax2.set_ylabel(r"$\Delta E$ (eV)", fontsize=11)
    ax2.set_title("Non-parabolicity deviation (code $-$ Kane)", fontsize=11)
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    out_path = FIGURES_DIR / "lecture_01_bulk_dispersion.png"
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"  Plot saved to {out_path}")
    print(f"  Max deviation: {np.max(np.abs(deviation)):.4f} eV "
          f"at k = {ks_code[np.argmax(np.abs(deviation))]:.4f} 1/A")
    print()


# =========================================================================
# Main
# =========================================================================
def main():
    print("\n" + "=" * 60)
    print("  LECTURE 01: Bulk Band Structure Validation")
    print("  8-band zinc-blende k.p Hamiltonian")
    print("=" * 60 + "\n")

    # Run all sections
    s1_pass = test_gaas_k0_eigenvalues()
    s2_pass, gaas_eig_path = test_gaas_effective_mass()
    s3_pass = test_inas_band_gap()

    # The GaAs eig_path from Section 2 is from a tempdir that no longer
    # exists. Re-run GaAs dispersion to get a fresh path for plotting.
    cfg = CONFIGS_DIR / "bulk_gaas_kx_dispersion.toml"
    with tempfile.TemporaryDirectory() as work:
        rc, outdir = run_exe(str(BUILD_DIR), "bandStructure",
                             str(cfg), work)
        if rc != 0:
            sys.exit(f"ERROR: bandStructure returned {rc} for GaAs plot")
        gaas_plot_path = os.path.join(outdir, "eigenvalues.dat")
        plot_gaas_dispersion(gaas_plot_path)

    # Summary
    print("=" * 60)
    print("  SUMMARY")
    print("=" * 60)
    results = [
        ("Section 1: GaAs k=0 eigenvalues", s1_pass),
        ("Section 2: GaAs effective mass", s2_pass),
        ("Section 3: InAs band gap", s3_pass),
        ("Section 4: Overlay plot", True),
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
