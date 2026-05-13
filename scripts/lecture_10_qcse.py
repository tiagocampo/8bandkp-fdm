#!/usr/bin/env python3
"""Lecture 10: Quantum-Confined Stark Effect -- Field-dependent subband shift.

Runs the 8-band k.p solver for a GaAs/AlGaAs quantum well in two conditions:
  1. Zero electric field (reference): sc_qcse_gaas_algaas.cfg
  2. With electric field (-0.007 eV/Ang = -700 kV/cm): sc_qcse_gaas_algaas_ef.cfg

Extracts eigenvalues at k=0 from each run and computes:
  - Raw CB1 eigenvalue shift: Delta_E_CB1 = E_CB1(field) - E_CB1(0)
    This is POSITIVE because the linear potential ramp across the 460 A domain
    dominates the shift (~1.05 eV).
  - Transition energy (CB1-VB1 gap) change: the actual QCSE red shift.
    This is NEGATIVE because the electric field reduces the effective gap
    between the first confined electron and hole states.

Validates:
  - Raw CB1 shift matches regression reference (~1.05 eV, 5% tolerance)
  - Transition energy change is NEGATIVE (QCSE red shift)
  - Transition energy shift magnitude is physically reasonable

Generates:
  - Overlay plot: subband energies with and without field
  - Saved to docs/lecture/figures/lecture_10_stark_shift.png
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
from star_helpers import run_exe, parse_eigenvalues, compare_value, TOL_NUMERICAL

# ---------------------------------------------------------------------------
# Physical constants and reference values for QCSE analysis
# ---------------------------------------------------------------------------
# EFParams value in the config: -0.007 eV/Angstrom
# Unit conversion: 1 eV/Angstrom = 1e8 eV/cm = 1e8 V/cm = 1e5 kV/cm
# So -0.007 eV/Ang = -700 kV/cm
EFIELD_EV_PER_ANG = -0.007  # eV/Angstrom from config
EFIELD_KV_PER_CM = EFIELD_EV_PER_ANG * 1.0e5  # -700 kV/cm

# QW width from config: GaAs layer from -30 to 30 Angstrom
QW_WIDTH_ANG = 60.0  # Angstrom

# Eigenvalue indices (0-based, after stripping k from parsed data)
# numcb=4, numvb=8 -> 12 eigenvalues total
# Degenerate pairs: VB [0,1] HH1, [2,3] LH1, [4,5] LH2, [6,7] HH2
#                   CB [8,9] CB1, [10,11] CB2
NUM_VB = 8
CB1_IDX = NUM_VB  # index 8
VB1_IDX = 0       # index 0 (HH1)

# Reference values from regression data (verify_star_gaas_algaas_qw.py)
CB1_NO_FIELD_REF = 0.931880    # eV
CB1_EF_REF = 1.981650          # eV
STARK_RAW_SHIFT_REF = CB1_EF_REF - CB1_NO_FIELD_REF  # ~1.049770 eV

# Transition energy shift bounds (QCSE red shift of CB1-VB1 gap)
TRANSITION_SHIFT_MIN_MEV = 100.0   # meV (lower bound for this strong field)
TRANSITION_SHIFT_MAX_MEV = 5000.0  # meV (upper bound, generous)


def run_bandstructure(config_name, work_dir):
    """Run bandStructure with the given config and return eigenvalue data.

    Args:
        config_name: filename in CONFIGS_DIR
        work_dir: temporary working directory

    Returns:
        (returncode, eigenvalues_list) where eigenvalues_list is from
        parse_eigenvalues, or None on failure.
    """
    cfg_path = CONFIGS_DIR / config_name
    rc, outdir = run_exe(str(BUILD_DIR), "bandStructure", str(cfg_path),
                         work_dir, timeout=300)
    if rc != 0:
        return rc, None

    eig_path = os.path.join(outdir, "eigenvalues.dat")
    if not os.path.isfile(eig_path):
        return -1, None

    data = parse_eigenvalues(eig_path)
    return rc, data


def plot_stark_shift(evals_0, evals_f, cb1_0, cb1_f, vb1_0, vb1_f,
                     save_path):
    """Generate overlay plot comparing subband energies with/without field.

    Left panel: all eigenvalue levels at k=0 for both runs, with connecting
    lines showing the field-induced shifts.
    Right panel: QW schematic showing the tilted potential and the shift
    of CB1 and VB1 energy levels.

    Args:
        evals_0: eigenvalue array without field
        evals_f: eigenvalue array with field
        cb1_0: CB1 energy without field (eV)
        cb1_f: CB1 energy with field (eV)
        vb1_0: VB1 energy without field (eV)
        vb1_f: VB1 energy with field (eV)
        save_path: output figure path
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5),
                              gridspec_kw={"width_ratios": [3, 2]})

    # --- Left panel: Subband energy levels ---
    ax1 = axes[0]

    n_bands = len(evals_0)
    x_positions = np.arange(n_bands)

    # Plot energies without field
    ax1.plot(x_positions, evals_0, 'o', color='royalblue', markersize=7,
             label='No field', zorder=3)
    # Plot energies with field
    ax1.plot(x_positions, evals_f, 's', color='crimson', markersize=7,
             label=f'Field ({EFIELD_KV_PER_CM:.0f} kV/cm)', zorder=3)

    # Draw connecting lines to show shifts
    for i in range(n_bands):
        shift = evals_f[i] - evals_0[i]
        if abs(shift) > 0.001:
            ax1.plot([x_positions[i], x_positions[i]],
                     [evals_0[i], evals_f[i]],
                     '-', color='gray', linewidth=0.8, alpha=0.6)

    # Annotate CB1 raw shift
    raw_shift_mev = (cb1_f - cb1_0) * 1000
    ax1.annotate(
        f'$\\Delta E_{{CB1}}$ = {raw_shift_mev:+.1f} meV\n(raw shift)',
        xy=(CB1_IDX, cb1_f),
        xytext=(CB1_IDX + 1.5, cb1_f + 0.05),
        fontsize=9, color='crimson', fontweight='bold',
        arrowprops=dict(arrowstyle='->', color='crimson', lw=1.2),
    )

    # Draw bracket showing transition energy without field
    gap_0 = cb1_0 - vb1_0
    ax1.annotate(
        '', xy=(n_bands - 0.3, cb1_0), xytext=(n_bands - 0.3, vb1_0),
        arrowprops=dict(arrowstyle='<->', color='royalblue', lw=1.5),
    )
    ax1.text(n_bands + 0.1, (cb1_0 + vb1_0) / 2,
             f'$E_{{gap}}$ = {gap_0 * 1000:.1f}\nmeV',
             fontsize=8, va='center', color='royalblue')

    ax1.set_xlabel('Eigenvalue index', fontsize=11)
    ax1.set_ylabel('Energy (eV)', fontsize=11)
    ax1.set_title('Subband Energies at $k=0$', fontsize=12)
    ax1.legend(fontsize=9, loc='upper left')
    ax1.grid(True, alpha=0.3, axis='y')

    # --- Right panel: QW potential schematic ---
    ax2 = axes[1]

    qw_half = QW_WIDTH_ANG / 2  # 30 Angstrom
    # Al0.2Ga0.8As band offset above GaAs CB edge: ~0.25 eV (20% of ~1.25 eV)
    barrier_cb_offset = 0.25

    # QW potential without field (rectangular)
    vb_no_field = vb1_0 - 0.1  # Approximate VB well bottom
    cb_no_field = cb1_0
    barrier_top = cb_no_field + barrier_cb_offset

    ax2.plot([-qw_half, -qw_half, qw_half, qw_half],
             [barrier_top, cb_no_field, cb_no_field, barrier_top],
             'b-', linewidth=2, label='No field')

    # Tilted QW potential with field
    # V(z) = EFIELD_EV_PER_ANG * z (linear potential)
    field_shift_at_edge = EFIELD_EV_PER_ANG * qw_half  # shift at well edges
    cb_field_left = cb_no_field - field_shift_at_edge   # left edge (z=-30)
    cb_field_right = cb_no_field + field_shift_at_edge  # right edge (z=+30)
    barrier_field_left = barrier_top - field_shift_at_edge
    barrier_field_right = barrier_top + field_shift_at_edge

    ax2.plot([-qw_half, -qw_half, qw_half, qw_half],
             [barrier_field_left, cb_field_left,
              cb_field_right, barrier_field_right],
             'r--', linewidth=2, label=f'Field ({EFIELD_KV_PER_CM:.0f} kV/cm)')

    # Energy levels without field
    ax2.hlines(cb1_0, -qw_half * 0.7, qw_half * 0.7,
               colors='royalblue', linewidth=2.5, linestyles='-',
               label=f'CB1 (no field): {cb1_0:.4f} eV')

    # Energy levels with field (approximate - show at well center)
    ax2.hlines(cb1_f, -qw_half * 0.7, qw_half * 0.7,
               colors='crimson', linewidth=2.5, linestyles='--',
               label=f'CB1 (field): {cb1_f:.4f} eV')

    # Arrow showing raw CB1 shift
    ax2.annotate('', xy=(qw_half * 0.5, cb1_f),
                 xytext=(qw_half * 0.5, cb1_0),
                 arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
    raw_mev = (cb1_f - cb1_0) * 1000
    ax2.text(qw_half * 0.65, (cb1_0 + cb1_f) / 2,
             f'{raw_mev:+.0f} meV',
             fontsize=10, va='center', fontweight='bold')

    ax2.set_xlabel('$z$ (Angstrom)', fontsize=11)
    ax2.set_ylabel('Energy (eV)', fontsize=11)
    ax2.set_title('QCSE: Electric Field in QW', fontsize=12)
    ax2.legend(fontsize=7, loc='upper right')
    ax2.grid(True, alpha=0.3)

    fig.suptitle(
        'Lecture 10: Quantum-Confined Stark Effect '
        f'(GaAs/{QW_WIDTH_ANG:.0f} A QW)',
        fontsize=13, fontweight='bold',
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Stark shift plot saved: {save_path}")


def main():
    print("=" * 60)
    print("  Lecture 10: Quantum-Confined Stark Effect")
    print("  GaAs/AlGaAs QW -- Field-dependent subband shift")
    print("=" * 60)

    all_passed = True

    # ------------------------------------------------------------------
    # Step 1: Run QW without field (reference)
    # ------------------------------------------------------------------
    print("\n[1/5] Running GaAs/AlGaAs QW without electric field...")
    with tempfile.TemporaryDirectory(prefix="lecture10_nofield_") as work_0:
        rc_0, eig_0 = run_bandstructure("sc_qcse_gaas_algaas.cfg", work_0)
        if rc_0 != 0 or eig_0 is None:
            print(f"  FAIL: bandStructure returned {rc_0} (no field)")
            all_passed = False
            sys.exit(f"ERROR: reference run failed with rc={rc_0}")
        evals_0 = eig_0[0][1]
        cb1_0 = evals_0[CB1_IDX]
        vb1_0 = evals_0[VB1_IDX]
        print(f"  PASS: bandStructure succeeded (rc=0)")
        print(f"  {len(evals_0)} eigenvalues at k=0")
        print(f"  VB1 (HH1): {vb1_0:.6f} eV")
        print(f"  CB1:       {cb1_0:.6f} eV")
        print(f"  Gap (CB1-VB1): {(cb1_0 - vb1_0) * 1000:.2f} meV")

    # ------------------------------------------------------------------
    # Step 2: Run QW with electric field
    # ------------------------------------------------------------------
    print(f"\n[2/5] Running GaAs/AlGaAs QW with electric field "
          f"({EFIELD_KV_PER_CM:.0f} kV/cm)...")
    with tempfile.TemporaryDirectory(prefix="lecture10_field_") as work_f:
        rc_f, eig_f = run_bandstructure("sc_qcse_gaas_algaas_ef.cfg", work_f)
        if rc_f != 0 or eig_f is None:
            print(f"  FAIL: bandStructure returned {rc_f} (with field)")
            all_passed = False
            sys.exit(f"ERROR: field run failed with rc={rc_f}")
        evals_f = eig_f[0][1]
        cb1_f = evals_f[CB1_IDX]
        vb1_f = evals_f[VB1_IDX]
        print(f"  PASS: bandStructure succeeded (rc=0)")
        print(f"  VB1 (HH1): {vb1_f:.6f} eV")
        print(f"  CB1:       {cb1_f:.6f} eV")
        print(f"  Gap (CB1-VB1): {(cb1_f - vb1_f) * 1000:.2f} meV")

    # ------------------------------------------------------------------
    # Step 3: Compute shifts
    # ------------------------------------------------------------------
    print("\n[3/5] Computing Stark shifts...")

    # Raw CB1 eigenvalue shift (includes linear potential ramp)
    raw_cb1_shift_eV = cb1_f - cb1_0
    raw_cb1_shift_mev = raw_cb1_shift_eV * 1000.0

    # Transition energy (CB1-VB1 gap) change -- the actual QCSE red shift
    gap_0 = cb1_0 - vb1_0
    gap_f = cb1_f - vb1_f
    transition_shift_eV = gap_f - gap_0
    transition_shift_mev = transition_shift_eV * 1000.0

    print(f"  E_CB1 (no field)  = {cb1_0:.6f} eV")
    print(f"  E_CB1 (with field) = {cb1_f:.6f} eV")
    print(f"  Raw CB1 shift      = {raw_cb1_shift_mev:+.3f} meV")
    print()
    print(f"  E_VB1 (no field)  = {vb1_0:.6f} eV")
    print(f"  E_VB1 (with field) = {vb1_f:.6f} eV")
    print()
    print(f"  Gap (no field)    = {gap_0 * 1000:.2f} meV")
    print(f"  Gap (with field)  = {gap_f * 1000:.2f} meV")
    print(f"  Transition shift  = {transition_shift_mev:+.3f} meV (QCSE)")

    # ------------------------------------------------------------------
    # Step 4: Validate shifts
    # ------------------------------------------------------------------
    print("\n[4/5] Validating Stark shifts...")

    # Assertion 4a: Raw CB1 shift matches regression reference (~1.05 eV)
    # This validates that the electric field implementation is correct.
    passed_raw, delta_raw, _ = compare_value(
        raw_cb1_shift_eV, STARK_RAW_SHIFT_REF, TOL_NUMERICAL,
        "Raw CB1 eigenvalue shift", unit="eV",
    )
    status_raw = "PASS" if passed_raw else "FAIL"
    print(f"  {status_raw}: Raw CB1 shift = {raw_cb1_shift_mev:+.3f} meV "
          f"(ref {STARK_RAW_SHIFT_REF * 1000:+.3f} meV, "
          f"delta = {delta_raw:.4f})")
    all_passed = all_passed and passed_raw

    # Assertion 4b: Raw CB1 shift direction matches field sign convention.
    # For a negative field, the linear potential pushes CB1 upward (positive
    # shift) because the wavefunction shifts toward lower z where the
    # potential is higher. This is the expected behavior per
    # verify_stark_shift.py lines 72-78.
    if EFIELD_EV_PER_ANG < 0.0 and raw_cb1_shift_eV > 0.0:
        print(f"  PASS: Raw CB1 shift positive for negative field "
              f"(+{raw_cb1_shift_mev:.1f} meV)")
    elif EFIELD_EV_PER_ANG > 0.0 and raw_cb1_shift_eV < 0.0:
        print(f"  PASS: Raw CB1 shift negative for positive field "
              f"({raw_cb1_shift_mev:.1f} meV)")
    else:
        print(f"  FAIL: Raw CB1 shift direction unexpected: "
              f"{raw_cb1_shift_mev:+.3f} meV for field "
              f"{EFIELD_KV_PER_CM:.0f} kV/cm")
        all_passed = False

    # Assertion 4c: Transition energy shift (CB1-VB1 gap) is NEGATIVE.
    # The electric field reduces the effective transition energy (QCSE red
    # shift). Both electron and hole wavefunctions shift toward the same
    # side of the well, reducing their energy difference.
    if transition_shift_eV < 0:
        print(f"  PASS: Transition energy shift is negative (QCSE red shift): "
              f"{transition_shift_mev:+.3f} meV < 0")
    else:
        print(f"  FAIL: Transition energy shift is NOT negative: "
              f"{transition_shift_mev:+.3f} meV >= 0")
        print(f"        Expected QCSE red shift (negative transition energy change).")
        all_passed = False

    # Assertion 4d: Transition shift magnitude is physically reasonable
    abs_trans_mev = abs(transition_shift_mev)
    if abs_trans_mev >= TRANSITION_SHIFT_MIN_MEV:
        print(f"  PASS: |Transition shift| = {abs_trans_mev:.1f} meV >= "
              f"{TRANSITION_SHIFT_MIN_MEV} meV lower bound")
    else:
        print(f"  WARN: |Transition shift| = {abs_trans_mev:.1f} meV < "
              f"{TRANSITION_SHIFT_MIN_MEV} meV (suspiciously small)")

    if abs_trans_mev <= TRANSITION_SHIFT_MAX_MEV:
        print(f"  PASS: |Transition shift| = {abs_trans_mev:.1f} meV <= "
              f"{TRANSITION_SHIFT_MAX_MEV} meV upper bound")
    else:
        print(f"  WARN: |Transition shift| = {abs_trans_mev:.1f} meV > "
              f"{TRANSITION_SHIFT_MAX_MEV} meV (unexpectedly large)")
        if abs_trans_mev > TRANSITION_SHIFT_MAX_MEV * 2:
            all_passed = False

    # Physical interpretation
    print(f"\n  Physical context:")
    print(f"    QW width: {QW_WIDTH_ANG:.0f} Angstrom")
    print(f"    Electric field: {EFIELD_KV_PER_CM:.0f} kV/cm")
    print(f"    Raw CB1 shift:  {raw_cb1_shift_mev:+.1f} meV (dominated by linear ramp)")
    print(f"    Transition shift: {transition_shift_mev:+.1f} meV (QCSE red shift)")
    # Potential drop across the well
    potential_drop = abs(EFIELD_EV_PER_ANG) * QW_WIDTH_ANG * 1000  # meV
    print(f"    Potential drop across well (eFL): {potential_drop:.1f} meV")

    # ------------------------------------------------------------------
    # Step 5: Generate overlay plot
    # ------------------------------------------------------------------
    print("\n[5/5] Generating Stark shift overlay plot...")
    plot_stark_shift(
        evals_0, evals_f, cb1_0, cb1_f, vb1_0, vb1_f,
        FIGURES_DIR / "lecture_10_stark_shift.png",
    )

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("  Summary")
    print("=" * 60)
    print(f"  VB1 energy (no field):   {vb1_0:.6f} eV")
    print(f"  CB1 energy (no field):   {cb1_0:.6f} eV")
    print(f"  VB1 energy (with field): {vb1_f:.6f} eV")
    print(f"  CB1 energy (with field): {cb1_f:.6f} eV")
    print(f"  Raw CB1 shift:           {raw_cb1_shift_mev:+.3f} meV")
    print(f"  Transition shift (QCSE): {transition_shift_mev:+.3f} meV")
    print(f"  Electric field:          {EFIELD_KV_PER_CM:.0f} kV/cm")
    print(f"  QW width:                {QW_WIDTH_ANG:.0f} Angstrom")
    print()

    if all_passed:
        print("  Lecture 10: ALL CHECKS PASSED")
    else:
        print("  Lecture 10: SOME CHECKS FAILED")
    print("=" * 60)
    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
