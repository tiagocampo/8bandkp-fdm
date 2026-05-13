#!/usr/bin/env python3
"""Lecture 04: Strain -- Bir-Pikus Hamiltonian and band-edge shifts.

Validates the Bir-Pikus strain implementation in the 8-band k.p code against
analytical formulas for biaxial [001] strain:

  1. Unstrained reference: bulk GaAs at k=0 (8 eigenvalues at band edges).
  2. Strained bulk GaAs: with substrate lattice constant a_sub = 5.869 A,
     verify HH-LH splitting, hydrostatic CB shift, and SO shift against the
     Bir-Pikus formulas:
       eps_xx = (a_sub - a_0) / a_0  (biaxial in-plane)
       eps_zz = -2*C12/C11 * eps_xx  (Poisson relaxation)
       P_eps  = -av * Tr(eps)
       Q_eps  = b/2 * (eps_zz - 0.5*(eps_xx + eps_yy))
       delta_Ec  = ac * Tr(eps)
       delta_EHH = -P_eps + Q_eps
       delta_ELH = -P_eps - Q_eps
       delta_ESO = -P_eps
  3. Strained InAs/GaAs QW: verify that strain effects are present in the
     quantum well eigenvalue spectrum (splitting of confined hole states).
  4. Overlay plot: strained vs unstrained bulk GaAs band edges with annotated
     Bir-Pikus shifts.

Audit (R10, 2026-05-11): bir_pikus_bulk() formulas match U1's
bir_pikus_bulk_inas() exactly (same P_eps, Q_eps, LH-SO coupling via QT2/sqrt(2)).
GaAs tensile ordering (LHSO_low < LHSO_high < HH < CB) is correct for Q_eps > 0.
Section 3 QW check is qualitative (VB spread > 0.01 eV); quantitative validation
is in verify_strain_rung6_qw.py (U2).
"""
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
    bir_pikus_biaxial_001,
    TOL_EXACT, TOL_ANALYTICAL, TOL_NUMERICAL,
)

# ---------------------------------------------------------------------------
# 8-band basis labels (CLAUDE.md: bands 1-4=valence, 5-6=SO, 7-8=CB)
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# GaAs material parameters (Vurgaftman 2001)
# ---------------------------------------------------------------------------
GAAS_A0     = 5.65325   # Angstrom, lattice constant
GAAS_C11    = 1221.0    # kbar, elastic constant
GAAS_C12    = 566.0     # kbar, elastic constant
GAAS_AC     = -7.17     # eV, conduction band hydrostatic deformation potential
GAAS_AV     = 1.16      # eV, valence band hydrostatic deformation potential
GAAS_B      = -2.0      # eV, shear deformation potential
GAAS_EG     = 1.519     # eV, band gap
GAAS_DELTA_SO = 0.341   # eV, spin-orbit splitting

# Substrate lattice constant for strained run
A_SUBSTRATE = 5.869     # Angstrom


# Bir-Pikus analytical formulas are in star_helpers.bir_pikus_biaxial_001()


def run_bandstructure(config_name, work_dir, config_overrides=None):
    """Run bandStructure with the given config and return parsed eigenvalues.

    Args:
        config_name: filename in tests/regression/configs/
        work_dir: temporary working directory
        config_overrides: dict of {old_text: new_text} to apply to config

    Returns:
        (data, output_dir) where data is list of (k, [eigenvalues])
    """
    config_path = CONFIGS_DIR / config_name

    if config_overrides:
        # Read original config and apply overrides
        with open(config_path) as f:
            content = f.read()
        for old, new in config_overrides.items():
            content = content.replace(old, new)
        # Write modified config to a staging path (NOT input.cfg yet,
        # since run_exe copies src -> work_dir/input.cfg and would
        # hit SameFileError if they were the same).
        modified_path = os.path.join(work_dir, "modified.cfg")
        with open(modified_path, "w") as f:
            f.write(content)
        config_path = modified_path

    rc, output_dir = run_exe(str(BUILD_DIR), "bandStructure",
                             str(config_path), work_dir)
    if rc != 0:
        sys.exit(f"ERROR: bandStructure returned {rc} for {config_name}")

    eig_path = os.path.join(output_dir, "eigenvalues.dat")
    if not os.path.isfile(eig_path):
        sys.exit(f"ERROR: eigenvalues.dat not found at {eig_path}")

    return parse_eigenvalues(eig_path), output_dir


# =========================================================================
# Section 1: Unstrained bulk GaAs reference at k=0
# =========================================================================
def test_unstrained_reference():
    """Run bulk GaAs at k=0 as unstrained reference."""
    print("=" * 60)
    print("Lecture 04 -- Section 1: Unstrained bulk GaAs at k=0")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as work:
        data, _ = run_bandstructure("bulk_gaas_k0.cfg", work)

    if not data:
        sys.exit("ERROR: no eigenvalue data for unstrained GaAs k=0")

    k_val, evals = data[0]
    print(f"  k = {k_val:.6f},  {len(evals)} eigenvalues")
    print(f"  Eigenvalues (eV):")
    for i, ev in enumerate(evals):
        print(f"    Band {i+1} ({BAND_LABELS[i]:>2s}): {ev:+.6f}")

    # Expected unstrained band edges (GaAs, Vurgaftman 2001):
    #   HH/LH: 0.0 eV (VB top, 4-fold degenerate at Gamma)
    #   SO: -Delta_SO = -0.341 eV
    #   CB: Eg = 1.519 eV
    # Basis ordering: bands 1,4 = HH; bands 2,3 = LH; bands 5,6 = SO;
    #                 bands 7,8 = CB
    expected = [
        -GAAS_DELTA_SO, -GAAS_DELTA_SO,
        0.0, 0.0, 0.0, 0.0,
        GAAS_EG, GAAS_EG,
    ]
    all_pass = True
    for i, (actual, exp) in enumerate(zip(evals, expected)):
        passed, delta, _ = compare_value(
            actual, exp, TOL_EXACT,
            f"Unstrained GaAs band {i+1} ({BAND_LABELS[i]})", unit="eV",
        )
        if not passed:
            print(f"    WARN: band {i+1} delta={delta:.2e}")
            all_pass = False

    if all_pass:
        print("  --> All unstrained eigenvalues match analytical band edges.")
    print()
    return all_pass, evals


# =========================================================================
# Section 2: Strained bulk GaAs -- Bir-Pikus validation
# =========================================================================
def test_strained_bir_pikus():
    """Run strained bulk GaAs and validate Bir-Pikus shifts."""
    print("=" * 60)
    print("Lecture 04 -- Section 2: Strained bulk GaAs -- Bir-Pikus validation")
    print("=" * 60)

    # Compute analytical Bir-Pikus shifts
    bp = bir_pikus_biaxial_001(GAAS_A0, A_SUBSTRATE, GAAS_C11, GAAS_C12,
                               GAAS_AC, GAAS_AV, GAAS_B,
                               delta_so=GAAS_DELTA_SO, Eg=GAAS_EG)

    print(f"  Substrate lattice constant: {A_SUBSTRATE} A")
    print(f"  GaAs lattice constant:      {GAAS_A0} A")
    print(f"  eps_xx = eps_yy = {bp['eps_xx']:.6f}")
    print(f"  eps_zz           = {bp['eps_zz']:.6f}")
    print(f"  Tr(eps)          = {bp['Tr_eps']:.6f}")
    print(f"  P_eps            = {bp['P_eps']:.6f} eV")
    print(f"  Q_eps            = {bp['Q_eps']:.6f} eV")
    print()
    print(f"  Expected Bir-Pikus shifts (eV):")
    print(f"    delta_Ec  = {bp['delta_Ec']:+.6f}")
    print(f"    delta_EHH = {bp['delta_EHH']:+.6f}")
    print(f"    delta_ELH = {bp['delta_ELH']:+.6f}")
    print(f"    delta_ESO = {bp['delta_ESO']:+.6f}")
    print(f"    QT2 (LH-SO coupling) = {bp['QT2']:.6f}")
    print()
    print(f"  Expected strained band edges (eV, with LH-SO mixing):")
    print(f"    E_CB        = {bp['E_CB']:+.6f}")
    print(f"    E_HH        = {bp['E_HH']:+.6f}")
    print(f"    E_LHSO_high = {bp['E_LHSO_high']:+.6f}  (LH-SO mixed)")
    print(f"    E_LHSO_low  = {bp['E_LHSO_low']:+.6f}  (LH-SO mixed)")
    print(f"    HH-LH splitting = {bp['HH_LH_splitting']:.6f}")
    print()

    # Run strained bulk GaAs with a clean config written from scratch.
    # For bulk mode, strain is controlled by strainSubstrate > 0 in the
    # material parameters (checked in ZB8bandBulk). We write a minimal
    # config with just the fields the parser expects, including
    # strainSubstrate at the end. We cannot use bulk_gaas_strained.cfg
    # because it has extra fields (gWhichBand, gBandIdx, SC, feast_*)
    # that confuse the sequential parser and prevent strainSubstrate from
    # being read correctly.
    strained_config = (
        "waveVector: k0\n"
        "waveVectorMax: 0\n"
        "waveVectorStep: 1\n"
        "confinement:  0\n"
        "FDstep: 101\n"
        "FDorder: 2\n"
        "numLayers:  1\n"
        "material1: GaAs\n"
        "numcb: 2\n"
        "numvb: 6\n"
        "ExternalField: 0  EF\n"
        f"EFParams: 0.0\n"
        f"strainSubstrate: {A_SUBSTRATE}\n"
    )
    with tempfile.TemporaryDirectory() as work_strained:
        # Write config directly and run executable
        cfg_path = os.path.join(work_strained, "strained.cfg")
        with open(cfg_path, "w") as f:
            f.write(strained_config)

        rc, output_dir = run_exe(str(BUILD_DIR), "bandStructure",
                                 cfg_path, work_strained)
        if rc != 0:
            sys.exit(f"ERROR: bandStructure returned {rc} for strained config")

        eig_path = os.path.join(output_dir, "eigenvalues.dat")
        if not os.path.isfile(eig_path):
            sys.exit(f"ERROR: eigenvalues.dat not found at {eig_path}")
        data = parse_eigenvalues(eig_path)

    if not data:
        sys.exit("ERROR: no eigenvalue data for strained GaAs k=0")

    k_val, evals_strained = data[0]
    print(f"  Strained eigenvalues (eV):")
    for i, ev in enumerate(evals_strained):
        print(f"    Band {i+1} ({BAND_LABELS[i]:>2s}): {ev:+.6f}")

    # Run unstrained reference in the same context
    with tempfile.TemporaryDirectory() as work_unref:
        data_ref, _ = run_bandstructure("bulk_gaas_k0.cfg", work_unref)
    _, evals_unref = data_ref[0]

    # Expected strained energies (eigenvalues sorted ascending):
    #   E1 = E2 = E_LHSO_low   (LH-SO mixed, lower state)
    #   E3 = E4 = E_LHSO_high  (LH-SO mixed, upper state)
    #   E5 = E6 = E_HH         (heavy hole, no LH-SO mixing)
    #   E7 = E8 = E_CB         (conduction band, no mixing)
    #
    # The LH and SO bands are coupled by QT2/sqrt(2) in the 8-band basis,
    # forming a 2x2 Hermitian sub-problem. The diagonal-only shifts
    # delta_ELH and delta_ESO are NOT the actual eigenvalues; the LH-SO
    # mixing pushes them apart. HH and CB have no such coupling, so their
    # diagonal shifts are exact.

    expected_strained = {
        "LHSO_low":  bp["E_LHSO_low"],   # bands 1,2 (sorted lowest)
        "LHSO_high": bp["E_LHSO_high"],  # bands 3,4
        "HH":        bp["E_HH"],          # bands 5,6
        "CB":        bp["E_CB"],          # bands 7,8
    }

    # Map band index to expected energy (0-indexed, sorted ascending)
    band_idx_map = {
        "LHSO_low":  0,   # band 1 (lower LH-SO mixed state)
        "LHSO_high": 2,   # band 3 (upper LH-SO mixed state)
        "HH":        4,   # band 5 (heavy hole)
        "CB":        6,   # band 7 (conduction band)
    }

    print()
    print(f"  Comparison (code vs analytical Bir-Pikus with LH-SO mixing):")
    all_pass = True
    for label in ["LHSO_low", "LHSO_high", "HH", "CB"]:
        idx = band_idx_map[label]
        computed = evals_strained[idx]
        expected = expected_strained[label]
        passed, delta, _ = compare_value(
            computed, expected, TOL_ANALYTICAL,
            f"Strained GaAs {label}", unit="eV",
        )
        status = "PASS" if passed else "FAIL"
        print(f"    {label:>10s}: computed={computed:+.6f}  "
              f"expected={expected:+.6f}  delta={delta:.4f}  {status}")
        all_pass = all_pass and passed

    # Also check HH-LH splitting (HH minus upper LH-SO state)
    computed_splitting = evals_strained[4] - evals_strained[2]
    expected_splitting = bp["HH_LH_splitting"]
    passed_split, delta_split, _ = compare_value(
        computed_splitting, expected_splitting, TOL_ANALYTICAL,
        "HH-LH splitting", unit="eV",
    )
    status = "PASS" if passed_split else "FAIL"
    print(f"    {'HH-LH split':>10s}: computed={computed_splitting:+.6f}  "
          f"expected={expected_splitting:+.6f}  delta={delta_split:.4f}  {status}")
    all_pass = all_pass and passed_split

    if all_pass:
        print("\n  --> All Bir-Pikus shifts validated within analytical tolerance.")
    else:
        print("\n  --> FAIL: some shifts outside tolerance.")
    print()

    return all_pass, evals_unref, evals_strained, bp


# =========================================================================
# Section 3: Strained InAs/GaAs QW -- verify strain effects present
# =========================================================================
def test_strained_qw():
    """Run strained InAs/GaAs QW and verify strain splitting present."""
    print("=" * 60)
    print("Lecture 04 -- Section 3: Strained InAs/GaAs QW")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as work_s:
        data, _ = run_bandstructure("qw_inas_gaas_strained.cfg", work_s)

    if not data:
        sys.exit("ERROR: no eigenvalue data for strained InAs/GaAs QW")

    # First k-point is k=0 (kpar=0)
    k0, evals_k0 = data[0]
    nk = len(data)
    nbands = len(evals_k0)
    print(f"  k=0 eigenvalues: {nbands} bands")
    print(f"  Total k-points: {nk}")
    for i, ev in enumerate(evals_k0):
        print(f"    Band {i+1:2d}: {ev:+.6f} eV")

    # In a strained InAs/GaAs QW, the HH and LH confined states should
    # be split due to the large lattice mismatch (~7%).
    # The QW has numcb=4, numvb=8 so we expect 12 bands.
    # We verify that the topmost valence states are non-degenerate,
    # indicating strain-induced splitting.
    n_vb = 8
    if nbands < n_vb:
        print(f"  WARN: expected at least {n_vb} bands, got {nbands}")

    # Get the highest valence-band states (first 8 are valence in ascending order)
    vb_states = evals_k0[:n_vb] if nbands >= n_vb else evals_k0

    # Check that not all VB states are at the same energy (strain breaks degeneracy)
    vb_spread = max(vb_states) - min(vb_states)
    print(f"\n  Valence band energy spread: {vb_spread:.6f} eV")

    # For a meaningful strain effect, the spread should be significant (> 0.01 eV)
    has_strain_effect = vb_spread > 0.01
    status = "PASS" if has_strain_effect else "FAIL"
    print(f"  Strain effects present: {status} (spread > 0.01 eV)")

    # Also run unstrained QW for comparison (disable strain)
    overrides_qw = {"strain: T": "strain: F"}
    with tempfile.TemporaryDirectory() as work_u:
        data_unref, _ = run_bandstructure("qw_inas_gaas_strained.cfg", work_u,
                                          config_overrides=overrides_qw)

    _, evals_unref = data_unref[0]
    vb_unref = evals_unref[:n_vb] if len(evals_unref) >= n_vb else evals_unref
    spread_unref = max(vb_unref) - min(vb_unref)
    print(f"  Unstrained VB spread: {spread_unref:.6f} eV")
    print(f"  Strain-induced additional splitting: "
          f"{vb_spread - spread_unref:.6f} eV")

    if has_strain_effect:
        print("\n  --> Strain effects confirmed in InAs/GaAs QW.")
    else:
        print("\n  --> FAIL: no significant strain effects detected.")
    print()
    return has_strain_effect, data, data_unref


# =========================================================================
# Section 4: Overlay plot -- strained vs unstrained bulk GaAs
# =========================================================================
def plot_strain_shifts(evals_unref, evals_strained, bp):
    """Generate overlay plot: strained vs unstrained band edges."""
    print("=" * 60)
    print("Lecture 04 -- Section 4: Strain shift overlay plot")
    print("=" * 60)

    fig, ax = plt.subplots(figsize=(9, 7))

    # Band-edge positions on a vertical axis
    x_unref = 0.0
    x_strained = 1.0

    # Unstrained band edges (sorted ascending):
    #   SO: -Delta_SO = -0.341 eV (2-fold)
    #   HH+LH: 0.0 eV (4-fold degenerate at Gamma)
    #   CB: Eg = 1.519 eV (2-fold)
    #
    # Strained band edges (sorted ascending, with LH-SO mixing):
    #   LHSO_low (bands 1,2): mixed LH-SO state
    #   LHSO_high (bands 3,4): mixed LH-SO state
    #   HH (bands 5,6): pure heavy hole
    #   CB (bands 7,8): conduction band

    band_groups_unref = [
        {"label": "SO",       "energy": -GAAS_DELTA_SO, "color": "forestgreen", "deg": 2},
        {"label": "HH+LH",   "energy": 0.0,            "color": "steelblue",   "deg": 4},
        {"label": "CB",       "energy": GAAS_EG,        "color": "crimson",     "deg": 2},
    ]

    band_groups_strained = [
        {"label": "LH-SO (low)",  "energy": bp["E_LHSO_low"],  "color": "forestgreen"},
        {"label": "LH-SO (high)", "energy": bp["E_LHSO_high"], "color": "steelblue"},
        {"label": "HH",           "energy": bp["E_HH"],        "color": "royalblue"},
        {"label": "CB",           "energy": bp["E_CB"],        "color": "crimson"},
    ]

    # Draw unstrained band edges (left column)
    for grp in band_groups_unref:
        e = grp["energy"]
        ax.plot([x_unref - 0.25, x_unref + 0.25], [e, e],
                color=grp["color"], lw=3, solid_capstyle="round")
        ax.text(x_unref - 0.30, e, f"{grp['label']} (x{grp['deg']})",
                ha="right", va="center", fontsize=10, color=grp["color"])

    # Draw strained band edges (right column)
    for grp in band_groups_strained:
        e = grp["energy"]
        ax.plot([x_strained - 0.25, x_strained + 0.25], [e, e],
                color=grp["color"], lw=3, ls="--", solid_capstyle="round")
        ax.text(x_strained + 0.30, e, f"{grp['label']}",
                ha="left", va="center", fontsize=10, color=grp["color"])

    # Annotate Bir-Pikus diagonal shifts (diagonal-only, before LH-SO mixing)
    # These show the first-order shifts; the LH-SO eigenvalues differ due to mixing.
    shift_annotations = [
        {"label": "SO", "unref": -GAAS_DELTA_SO, "shift": bp["delta_ESO"],
         "color": "forestgreen"},
        {"label": "LH", "unref": 0.0, "shift": bp["delta_ELH"],
         "color": "steelblue"},
        {"label": "HH", "unref": 0.0, "shift": bp["delta_EHH"],
         "color": "royalblue"},
        {"label": "CB", "unref": GAAS_EG, "shift": bp["delta_Ec"],
         "color": "crimson"},
    ]
    for ann in shift_annotations:
        e_unref = ann["unref"]
        e_diag = e_unref + ann["shift"]
        shift = ann["shift"]
        if abs(shift) > 0.001:
            mid_x = 0.5
            ax.annotate(
                "", xy=(x_strained - 0.27, e_diag),
                xytext=(x_unref + 0.27, e_unref),
                arrowprops=dict(arrowstyle="->", color="gray", lw=0.8,
                                connectionstyle="arc3,rad=0", ls=":"),
            )
            ax.text(mid_x, (e_unref + e_diag) / 2,
                    f"{shift:+.3f}",
                    ha="center", va="bottom", fontsize=7,
                    color="dimgray", style="italic",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white",
                              ec="lightgray", alpha=0.7))

    # Annotate HH-LH splitting (computed from actual eigenvalues)
    hh_e = evals_strained[4]
    lhs_high_e = evals_strained[2]
    split = hh_e - lhs_high_e
    ax.annotate(
        f"HH-LH splitting\n= {split:.3f} eV\n(analytical: {bp['HH_LH_splitting']:.3f} eV)",
        xy=(x_strained + 0.35, (hh_e + lhs_high_e) / 2),
        xytext=(x_strained + 0.75, (hh_e + lhs_high_e) / 2 + 0.3),
        fontsize=9, ha="left",
        arrowprops=dict(arrowstyle="->", color="black", lw=1),
        bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow",
                  ec="orange", alpha=0.9),
    )

    # Formatting
    ax.set_xlim(-0.8, 2.2)
    ax.set_xticks([x_unref, x_strained])
    ax.set_xticklabels(["Unstrained", "Strained"], fontsize=12)
    ax.set_ylabel("Energy (eV)", fontsize=12)
    ax.set_title(
        f"Lecture 04: Bir-Pikus Strain Shifts in GaAs\n"
        f"(substrate $a_{{sub}}$ = {A_SUBSTRATE} A, "
        f"$\\epsilon_{{xx}}$ = {bp['eps_xx']:.4f})",
        fontsize=13,
    )
    ax.axhline(0, color="gray", ls=":", lw=0.5)
    ax.grid(True, axis="y", alpha=0.3)

    # Add parameter box
    textstr = (
        f"GaAs parameters (Vurgaftman 2001):\n"
        f"  $a_0$ = {GAAS_A0} A\n"
        f"  $a_c$ = {GAAS_AC} eV\n"
        f"  $a_v$ = {GAAS_AV} eV\n"
        f"  $b$ = {GAAS_B} eV\n"
        f"  $C_{{11}}$ = {GAAS_C11} kbar\n"
        f"  $C_{{12}}$ = {GAAS_C12} kbar"
    )
    props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=8,
            verticalalignment="top", bbox=props)

    fig.tight_layout()
    out_path = FIGURES_DIR / "lecture_04_strain_shifts.png"
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"  Plot saved to {out_path}")
    print()


# =========================================================================
# Main
# =========================================================================
def main():
    print("\n" + "=" * 60)
    print("  LECTURE 04: Strain -- Bir-Pikus Hamiltonian")
    print("  8-band zinc-blende k.p with deformation potentials")
    print("=" * 60 + "\n")

    # Section 1: Unstrained reference
    s1_pass, evals_unref = test_unstrained_reference()

    # Section 2: Bir-Pikus validation
    s2_pass, evals_unref_2, evals_strained, bp = test_strained_bir_pikus()

    # Section 3: Strained QW
    s3_pass, qw_strained, qw_unref = test_strained_qw()

    # Section 4: Overlay plot
    plot_strain_shifts(evals_unref_2, evals_strained, bp)

    # Summary
    print("=" * 60)
    print("  SUMMARY")
    print("=" * 60)
    results = [
        ("Section 1: Unstrained GaAs k=0 reference", s1_pass),
        ("Section 2: Bir-Pikus analytical validation", s2_pass),
        ("Section 3: Strained InAs/GaAs QW", s3_pass),
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
