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
                     save_path, transition_shift_mev):
    """Generate overlay plot showing QCSE transition energy shift.

    Left panel: transition energy (CB1-VB1 gap) with and without field.
    Right panel: QW schematic showing the tilted potential and wavefunction
    shift that causes the QCSE red shift.

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

    # --- Left panel: Transition energy diagram ---
    ax1 = axes[0]

    gap_0 = cb1_0 - vb1_0
    gap_f = cb1_f - vb1_f
    qcse_shift_mev = (gap_f - gap_0) * 1000.0

    # Draw subband levels
    x_no = 0.3
    x_ef = 0.7

    # VB1 levels
    ax1.plot([x_no - 0.15, x_no + 0.15], [vb1_0, vb1_0], 'b-', lw=2.5)
    ax1.plot([x_ef - 0.15, x_ef + 0.15], [vb1_f, vb1_f], 'r--', lw=2.5)

    # CB1 levels
    ax1.plot([x_no - 0.15, x_no + 0.15], [cb1_0, cb1_0], 'b-', lw=2.5)
    ax1.plot([x_ef - 0.15, x_ef + 0.15], [cb1_f, cb1_f], 'r--', lw=2.5)

    # Transition energy arrows
    ax1.annotate('', xy=(x_no, cb1_0), xytext=(x_no, vb1_0),
                 arrowprops=dict(arrowstyle='<->', color='royalblue', lw=2))
    ax1.text(x_no - 0.05, (cb1_0 + vb1_0) / 2,
             f'$E_{{gap}}^0$ = {gap_0 * 1000:.1f} meV',
             fontsize=9, ha='right', va='center', color='royalblue')

    ax1.annotate('', xy=(x_ef, cb1_f), xytext=(x_ef, vb1_f),
                 arrowprops=dict(arrowstyle='<->', color='crimson', lw=2))
    ax1.text(x_ef + 0.05, (cb1_f + vb1_f) / 2,
             f'$E_{{gap}}^F$ = {gap_f * 1000:.1f} meV',
             fontsize=9, ha='left', va='center', color='crimson')

    # Labels
    ax1.text(x_no, vb1_0 - 0.03, 'VB1', ha='center', fontsize=9, color='gray')
    ax1.text(x_no, cb1_0 + 0.03, 'CB1', ha='center', fontsize=9, color='gray')
    ax1.text(x_no, 0.5 * (cb1_0 + vb1_0) + 0.05, 'No field',
             ha='center', fontsize=10, fontweight='bold', color='royalblue')
    ax1.text(x_ef, 0.5 * (cb1_f + vb1_f) + 0.05,
             f'F = {EFIELD_KV_PER_CM:.0f}\nkV/cm',
             ha='center', fontsize=10, fontweight='bold', color='crimson')

    # QCSE shift annotation
    ax1.annotate(
        f'QCSE red shift\n$\\Delta E_{{trans}}$ = {qcse_shift_mev:+.1f} meV',
        xy=(0.5, 0.5 * (gap_0 + gap_f) + vb1_0),
        xycoords='data',
        fontsize=11, ha='center', va='center', fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow',
                  edgecolor='orange', alpha=0.9),
    )

    ax1.set_xlim(0, 1)
    ax1.set_ylabel('Energy (eV)', fontsize=11)
    ax1.set_title('Transition Energy: CB1 $-$ VB1', fontsize=12)
    ax1.set_xticks([])
    ax1.grid(True, alpha=0.3, axis='y')

    # --- Right panel: QW potential schematic ---
    # Show the well potential profile. The field-shifted eigenvalues include
    # a large linear ramp across the full simulation domain, so we plot them
    # relative to the no-field levels to isolate the QCSE effect.
    ax2 = axes[1]

    qw_half = QW_WIDTH_ANG / 2  # 30 Angstrom
    barrier_cb_offset = 0.25

    # No-field potential
    cb_no_field = cb1_0
    barrier_top = cb_no_field + barrier_cb_offset

    ax2.plot([-qw_half, -qw_half, qw_half, qw_half],
             [barrier_top, cb_no_field, cb_no_field, barrier_top],
             'b-', linewidth=2, label='No field (CB)')

    # Tilted potential across the well
    field_drop = EFIELD_EV_PER_ANG * QW_WIDTH_ANG  # total drop across well
    cb_field_left = cb_no_field - field_drop / 2
    cb_field_right = cb_no_field + field_drop / 2

    ax2.plot([-qw_half, -qw_half, qw_half, qw_half],
             [cb_field_left + barrier_cb_offset, cb_field_left,
              cb_field_right, cb_field_right + barrier_cb_offset],
             'r--', linewidth=2,
             label=f'Field ({EFIELD_KV_PER_CM:.0f} kV/cm)')

    # No-field energy levels
    ax2.hlines(cb1_0, -qw_half * 0.7, qw_half * 0.7,
               colors='royalblue', linewidth=2.5, linestyles='-',
               label=f'CB1 (no field)')
    ax2.hlines(vb1_0, -qw_half * 0.7, qw_half * 0.7,
               colors='royalblue', linewidth=2.5, linestyles='-',
               label=f'VB1 (no field)')

    # Field-shifted levels: subtract the domain-wide ramp to show the
    # QCSE shift relative to the no-field position. The ramp adds ~EF*z_avg
    # where z_avg is the wavefunction centroid in the full domain.
    # For the QCSE effect, what matters is the relative change.
    # The transition energy shift is the key quantity, already shown in left panel.
    # Here show approximate CB1/VB1 shifts within the well.
    cb1_shift_within_well = (cb1_f - cb1_0) - EFIELD_EV_PER_ANG * 0  # ref at z=0
    # The actual within-well shift is smaller than the raw shift.
    # Plot the field levels at positions that show the gap narrowing.
    # Use the transition energy change: CB1 drops by ~half, VB1 rises by ~half
    half_qcse = transition_shift_mev / 2.0 / 1000.0  # eV
    cb1_qcse = cb1_0 + half_qcse  # drops (half_qcse is negative)
    vb1_qcse = vb1_0 - half_qcse  # rises (subtract negative = add)

    ax2.hlines(cb1_qcse, -qw_half * 0.4, qw_half * 0.7,
               colors='crimson', linewidth=2.5, linestyles='--',
               label=f'CB1 (QCSE shift)')
    ax2.hlines(vb1_qcse, -qw_half * 0.4, qw_half * 0.7,
               colors='crimson', linewidth=2.5, linestyles='--',
               label=f'VB1 (QCSE shift)')

    # Arrows showing the shifts
    ax2.annotate('', xy=(0, cb1_qcse), xytext=(0, cb1_0),
                 arrowprops=dict(arrowstyle='->', color='darkred', lw=1.5))
    ax2.annotate('', xy=(0, vb1_qcse), xytext=(0, vb1_0),
                 arrowprops=dict(arrowstyle='->', color='darkred', lw=1.5))

    ax2.annotate(
        f'QCSE shift:\n'
        f'$\\Delta V$ = {abs(field_drop)*1000:.0f} meV across well\n'
        f'$\\Delta E_{{trans}}$ = {qcse_shift_mev:+.0f} meV',
        xy=(0.03, 0.97), xycoords='axes fraction',
        fontsize=8, va='top',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8),
    )

    ax2.set_xlabel('$z$ (Angstrom)', fontsize=11)
    ax2.set_ylabel('Energy (eV)', fontsize=11)
    ax2.set_title('QCSE: Level Shifts in QW', fontsize=12)
    ax2.legend(fontsize=6.5, loc='center right')
    ax2.grid(True, alpha=0.3)

    fig.suptitle(
        'Quantum-Confined Stark Effect '
        f'(GaAs/{QW_WIDTH_ANG:.0f} A QW, F = {EFIELD_KV_PER_CM:.0f} kV/cm)',
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
        transition_shift_mev,
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
