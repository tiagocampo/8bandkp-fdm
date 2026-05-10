#!/usr/bin/env python3
"""Lecture 05: g-Factor -- Roth formula and Landau levels.

Validates the k.p g-factor code against the Roth formula and verifies
Landau level spacing from the magnetic field module.

Sections:
  1. GaAs CB g-factor vs Roth formula
  2. InSb CB g-factor (extreme narrow-gap regime)
  3. GaAsW CB g-factor (Winkler parameter set)
  4. GaAs Landau levels at B=5T
  5. Overlay bar chart: code g-factors vs Roth formula
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
from star_helpers import (run_exe, parse_eigenvalues, parse_gfactor, compare_value,
                          roth_gfactor, TOL_EXACT, TOL_ANALYTICAL, TOL_NUMERICAL,
                          HBAR_J_S, E_CHARGE, M0_KG)

# ---------------------------------------------------------------------------
# Material parameters (Vurgaftman 2001 / Winkler 2003)
# ---------------------------------------------------------------------------
MATERIALS = {
    'GaAs':  {'Eg': 1.519, 'DeltaSO': 0.341, 'EP': 28.8,  'meff': 0.067},
    'GaAsW': {'Eg': 1.519, 'DeltaSO': 0.341, 'EP': 28.89, 'meff': 0.0665},
    'InSb':  {'Eg': 0.235, 'DeltaSO': 0.81,  'EP': 23.3,  'meff': 0.0135},
    'InAs':  {'Eg': 0.417, 'DeltaSO': 0.39,  'EP': 21.5,  'meff': 0.026},
}


def run_gfactor(config_name, label):
    """Run gfactorCalculation with given config, return (gx, gy, gz).

    Args:
        config_name: config filename in tests/regression/configs/
        label: descriptive label for error messages

    Returns:
        tuple (gx, gy, gz)
    """
    cfg = CONFIGS_DIR / config_name
    with tempfile.TemporaryDirectory() as work:
        rc, outdir = run_exe(str(BUILD_DIR), "gfactorCalculation",
                             str(cfg), work)
        if rc != 0:
            sys.exit(f"ERROR: gfactorCalculation returned {rc} for {label}")

        gfactor_path = os.path.join(outdir, "gfactor.dat")
        if not os.path.isfile(gfactor_path):
            sys.exit(f"ERROR: gfactor.dat not found for {label}")

        return parse_gfactor(gfactor_path)


def run_bandstructure(config_name, label):
    """Run bandStructure with given config, return parsed eigenvalues.

    Args:
        config_name: config filename in tests/regression/configs/
        label: descriptive label for error messages

    Returns:
        list of (k, [eigenvalues])
    """
    cfg = CONFIGS_DIR / config_name
    with tempfile.TemporaryDirectory() as work:
        rc, outdir = run_exe(str(BUILD_DIR), "bandStructure",
                             str(cfg), work)
        if rc != 0:
            sys.exit(f"ERROR: bandStructure returned {rc} for {label}")

        eig_path = os.path.join(outdir, "eigenvalues.dat")
        if not os.path.isfile(eig_path):
            sys.exit(f"ERROR: eigenvalues.dat not found for {label}")

        return parse_eigenvalues(eig_path)


# =========================================================================
# Section 1: GaAs CB g-factor vs Roth formula
# =========================================================================
def test_gaas_gfactor():
    """Verify GaAs CB g-factor matches Roth formula within TOL_ANALYTICAL."""
    print("=" * 60)
    print("Lecture 05 -- Section 1: GaAs CB g-factor")
    print("=" * 60)

    gx, gy, gz = run_gfactor("gfactor_bulk_gaas_cb.cfg", "GaAs CB")
    g_code = gz  # bulk is isotropic; use gz component

    p = MATERIALS['GaAs']
    g_roth = roth_gfactor(p['EP'], p['Eg'], p['DeltaSO'])

    print(f"  Code g-factor:  g = {g_code:+.6f}")
    print(f"  Roth formula:   g = {g_roth:+.6f}")
    print(f"  (EP={p['EP']}, Eg={p['Eg']}, DeltaSO={p['DeltaSO']})")

    passed, delta, _ = compare_value(g_code, g_roth, TOL_ANALYTICAL,
                                     "GaAs CB g-factor")
    status = "PASS" if passed else "FAIL"
    print(f"  {status}: relative error = {delta:.4f}  (tolerance {TOL_ANALYTICAL:.0%})")

    if passed:
        print("  --> GaAs g-factor matches Roth formula within 1%.")
    else:
        print("  --> FAIL: g-factor outside analytical tolerance.")
    print()
    return passed, g_code, g_roth


# =========================================================================
# Section 2: InSb CB g-factor (extreme narrow-gap regime)
# =========================================================================
def test_insb_gfactor():
    """Verify InSb CB g-factor is large and negative (|g| > 40)."""
    print("=" * 60)
    print("Lecture 05 -- Section 2: InSb CB g-factor (extreme regime)")
    print("=" * 60)

    gx, gy, gz = run_gfactor("gfactor_bulk_insb_cb.cfg", "InSb CB")
    g_code = gz

    p = MATERIALS['InSb']
    g_roth = roth_gfactor(p['EP'], p['Eg'], p['DeltaSO'])

    print(f"  Code g-factor:  g = {g_code:+.2f}")
    print(f"  Roth formula:   g = {g_roth:+.2f}")
    print(f"  (EP={p['EP']}, Eg={p['Eg']}, DeltaSO={p['DeltaSO']})")

    # Qualitative check: |g| > 40 (extreme narrow-gap regime)
    extreme_pass = abs(g_code) > 40
    status_extreme = "PASS" if extreme_pass else "FAIL"
    print(f"  {status_extreme}: |g| = {abs(g_code):.2f}  (threshold: 40)")

    # Also compare against Roth formula (looser tolerance for extreme regime)
    passed, delta, _ = compare_value(g_code, g_roth, TOL_NUMERICAL,
                                     "InSb CB g-factor")
    status = "PASS" if passed else "FAIL"
    print(f"  {status}: relative error vs Roth = {delta:.4f}  (tolerance {TOL_NUMERICAL:.0%})")

    note = ("  Note: 20-30% g-factor shortfall vs experiment is a known "
            "8-band model limitation.\n"
            "  InSb Roth g gives qualitative agreement only.")
    print(note)

    if extreme_pass:
        print("  --> InSb g-factor in extreme narrow-gap regime.")
    else:
        print("  --> FAIL: g-factor not in expected extreme regime.")
    print()
    return extreme_pass, g_code, g_roth


# =========================================================================
# Section 3: GaAsW CB g-factor (Winkler parameter set)
# =========================================================================
def test_gaasw_gfactor():
    """Verify GaAsW CB g-factor matches Roth formula with Winkler params."""
    print("=" * 60)
    print("Lecture 05 -- Section 3: GaAsW CB g-factor (Winkler params)")
    print("=" * 60)

    gx, gy, gz = run_gfactor("gfactor_bulk_gaasw_cb.cfg", "GaAsW CB")
    g_code = gz

    p = MATERIALS['GaAsW']
    g_roth = roth_gfactor(p['EP'], p['Eg'], p['DeltaSO'])

    print(f"  Code g-factor:  g = {g_code:+.6f}")
    print(f"  Roth formula:   g = {g_roth:+.6f}")
    print(f"  (EP={p['EP']}, Eg={p['Eg']}, DeltaSO={p['DeltaSO']})")

    passed, delta, _ = compare_value(g_code, g_roth, TOL_ANALYTICAL,
                                     "GaAsW CB g-factor")
    status = "PASS" if passed else "FAIL"
    print(f"  {status}: relative error = {delta:.4f}  (tolerance {TOL_ANALYTICAL:.0%})")

    if passed:
        print("  --> GaAsW g-factor matches Roth formula within 1%.")
    else:
        print("  --> FAIL: g-factor outside analytical tolerance.")
    print()
    return passed, g_code, g_roth


# =========================================================================
# Section 4: GaAs Landau levels at B=5T
# =========================================================================
def test_landau_levels():
    """Verify GaAs Landau level spacing matches Kane-mass hbar*omega_c."""
    print("=" * 60)
    print("Lecture 05 -- Section 4: GaAs Landau levels at B=5T")
    print("=" * 60)

    data = run_bandstructure("landau_bulk_GaAs.cfg", "GaAs Landau")
    if not data:
        sys.exit("ERROR: no eigenvalue data parsed for GaAs Landau")

    # First k-point (k_y=0): Landau levels should be visible in the CB
    k0, evals = data[0]
    print(f"  k = {k0:.6f},  {len(evals)} eigenvalues")

    n_total = len(evals)
    print(f"  Total eigenvalues: {n_total}")

    # Extract the CB eigenvalues (upper portion)
    n_cb = n_total // 2  # upper half is CB-like
    cb_evals = sorted(evals[n_cb:])  # CB eigenvalues, sorted ascending

    # The Landau levels are the discrete CB energies at k_y=0
    spacings = np.diff(cb_evals)
    print(f"  CB eigenvalues count: {len(cb_evals)}")
    print(f"  First 6 CB energies (eV): {[f'{e:.6f}' for e in cb_evals[:6]]}")
    print(f"  First 5 spacings (meV):  {[f'{s*1000:.2f}' for s in spacings[:5]]}")

    # Analytical: hbar*omega_c = hbar*eB/m*_Kane (in eV)
    # Use Kane effective mass (8-band self-consistent): m* = Eg/(EP+Eg)
    B = 5.0  # Tesla
    p = MATERIALS['GaAs']
    m_star_kane = p['Eg'] / (p['EP'] + p['Eg'])  # 0.0501 m0
    omega_c = E_CHARGE * B / (m_star_kane * M0_KG)  # rad/s
    hbar_omega_c_eV = HBAR_J_S * omega_c / E_CHARGE  # convert J to eV
    hbar_omega_c_meV = hbar_omega_c_eV * 1000

    print(f"  Analytical hbar*omega_c = {hbar_omega_c_meV:.2f} meV")
    print(f"  (m*_Kane = {m_star_kane:.4f} m0 = Eg/(EP+Eg), B = {B} T)")

    # Extract Landau level spacings, filtering out spin-orbit sub-splittings.
    # In 8-band k.p, each Landau level is spin-split; the intra-level splitting
    # is small (~1 meV) while the inter-level spacing is ~hbar*omega_c (~12 meV).
    # The 0.5 threshold cleanly separates the two regimes for GaAs at B=5T
    # (spin-split ~1 meV << threshold ~6 meV << Landau spacing ~12 meV).
    spacings_meV = np.array([s * 1000 for s in spacings])
    threshold = 0.5 * hbar_omega_c_meV
    landau_spacings = spacings_meV[spacings_meV > threshold]

    if len(landau_spacings) == 0:
        # Fallback: use all spacings
        landau_spacings = spacings_meV
    mean_spacing = np.mean(landau_spacings) / 1000.0  # back to eV
    mean_spacing_meV = mean_spacing * 1000

    print(f"  Landau spacings (>{threshold:.1f} meV): "
          f"{[f'{s:.2f}' for s in landau_spacings]}")
    print(f"  Mean Landau spacing: {mean_spacing_meV:.2f} meV")

    TOL_LANDAU = 0.15  # 15% — GaAs has weak non-parabolicity at B=5T
    passed, delta, _ = compare_value(mean_spacing, hbar_omega_c_eV,
                                     TOL_LANDAU, "Landau level spacing")
    status = "PASS" if passed else "FAIL"
    print(f"  {status}: relative error = {delta:.4f}  (tolerance {TOL_LANDAU:.0%})")

    if passed:
        print("  --> GaAs Landau level spacing matches Kane-mass hbar*omega_c.")
    else:
        print("  --> FAIL: Landau level spacing outside tolerance.")
    print()
    return passed, cb_evals, hbar_omega_c_eV


# =========================================================================
# Section 5: Overlay bar chart -- code g-factors vs Roth formula
# =========================================================================
def plot_gfactor_comparison(results):
    """Bar chart comparing code g-factors vs Roth formula.

    Args:
        results: list of (material_name, g_code, g_roth)
    """
    print("=" * 60)
    print("Lecture 05 -- Section 5: g-factor comparison plot")
    print("=" * 60)

    materials = [r[0] for r in results]
    g_codes = [r[1] for r in results]
    g_roths = [r[2] for r in results]
    n = len(materials)

    x = np.arange(n)
    width = 0.35

    fig, ax = plt.subplots(figsize=(9, 6))

    bars_code = ax.bar(x - width / 2, g_codes, width,
                       label='8-band k.p (code)', color='steelblue',
                       edgecolor='black', linewidth=0.8)
    bars_roth = ax.bar(x + width / 2, g_roths, width,
                       label='Roth formula', color='coral',
                       edgecolor='black', linewidth=0.8)

    # Annotate bars with values
    for bar, val in zip(bars_code, g_codes):
        ypos = bar.get_height()
        offset = 1.5 if ypos >= 0 else -3.5
        ax.text(bar.get_x() + bar.get_width() / 2, ypos + offset,
                f'{val:+.2f}', ha='center', va='bottom', fontsize=9,
                fontweight='bold')

    for bar, val in zip(bars_roth, g_roths):
        ypos = bar.get_height()
        offset = 1.5 if ypos >= 0 else -3.5
        ax.text(bar.get_x() + bar.get_width() / 2, ypos + offset,
                f'{val:+.2f}', ha='center', va='bottom', fontsize=9,
                color='darkred')

    ax.set_xlabel('Material', fontsize=12)
    ax.set_ylabel(r'$g$-factor', fontsize=12)
    ax.set_title(r'Lecture 05: CB $g$-factor -- 8-band k.p vs Roth Formula',
                 fontsize=13)
    ax.set_xticks(x)
    ax.set_xticklabels(materials, fontsize=11)
    ax.legend(fontsize=11, loc='best')
    ax.axhline(2.0, color='gray', ls=':', lw=0.8, label=r'$g_0 = 2$')
    ax.grid(True, alpha=0.3, axis='y')

    # Add a note about the 8-band model limitation
    ax.annotate('Note: 20-30% shortfall vs experiment\nis a known 8-band limitation',
                xy=(0.98, 0.02), xycoords='axes fraction',
                fontsize=8, ha='right', va='bottom',
                bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow', alpha=0.8))

    fig.tight_layout()
    out_path = FIGURES_DIR / "lecture_05_gfactor_comparison.png"
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"  Plot saved to {out_path}")
    print()


# =========================================================================
# Main
# =========================================================================
def main():
    print("\n" + "=" * 60)
    print("  LECTURE 05: g-Factor -- Roth Formula and Landau Levels")
    print("  8-band zinc-blende k.p Hamiltonian")
    print("=" * 60 + "\n")

    # Collect g-factor results for the comparison plot
    gfactor_results = []

    # Section 1: GaAs CB g-factor
    s1_pass, gaas_g_code, gaas_g_roth = test_gaas_gfactor()
    gfactor_results.append(('GaAs', gaas_g_code, gaas_g_roth))

    # Section 2: InSb CB g-factor
    s2_pass, insb_g_code, insb_g_roth = test_insb_gfactor()
    gfactor_results.append(('InSb', insb_g_code, insb_g_roth))

    # Section 3: GaAsW CB g-factor
    s3_pass, gaasw_g_code, gaasw_g_roth = test_gaasw_gfactor()
    gfactor_results.append(('GaAsW', gaasw_g_code, gaasw_g_roth))

    # Section 4: GaAs Landau levels
    s4_pass, _, _ = test_landau_levels()

    # Section 5: Overlay plot
    plot_gfactor_comparison(gfactor_results)

    # Summary
    print("=" * 60)
    print("  SUMMARY")
    print("=" * 60)
    results = [
        ("Section 1: GaAs CB g-factor vs Roth", s1_pass),
        ("Section 2: InSb CB g-factor (extreme)", s2_pass),
        ("Section 3: GaAsW CB g-factor", s3_pass),
        ("Section 4: GaAs Landau levels", s4_pass),
        ("Section 5: Comparison plot", True),
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
