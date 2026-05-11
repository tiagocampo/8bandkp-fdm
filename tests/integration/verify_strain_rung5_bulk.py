#!/usr/bin/env python3
"""Rung 5 — Bulk Strain Verification (U1).

Validates the full bulk strain chain: strainSubstrate -> apply_bp_strain_inline
produces correct strained eigenvalues, HH-LH splitting, and additive modification.

Requirements: R1, R2, R3, R11, R13, R14.

Usage:
    verify_strain_rung5_bulk.py <build_dir> <source_dir>
"""

# COVERAGE: observable=strain_shift geometry=bulk material=InAs/GaAs ref=Vurgaftman2001
# COVERAGE: observable=HH_LH_splitting geometry=bulk material=InAs/GaAs ref=Vurgaftman2001

import os
import shutil
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from star_helpers import (
    run_exe, parse_eigenvalues, compare_value,
    TOL_EXACT, TOL_ANALYTICAL, TOL_NUMERICAL,
)

# ---------------------------------------------------------------------------
# InAs material parameters (Vurgaftman 2001, parameters.f90)
# ---------------------------------------------------------------------------
INAS_A0    = 6.0583    # Angstrom
INAS_C11   = 832.9     # kbar
INAS_C12   = 452.6     # kbar
INAS_AC    = -5.08     # eV
INAS_AV    = 1.00      # eV
INAS_B     = -1.8      # eV  (b_dp)
INAS_D     = -3.6      # eV  (d_dp)
INAS_EG    = 0.417     # eV
INAS_DELTA = 0.390     # eV  (Delta_SO)

# GaAs substrate lattice constant
GAAS_A0_SUB = 5.65325  # Angstrom

# Tolerances
TOL_STRAIN = 0.01      # 1% for analytical Bir-Pikus comparison
TOL_ADDITIVE = 1e-5    # eigenvalue solver precision for additive check (CB, HH)


def bir_pikus_bulk_inas(a0, a_sub, C11, C12, ac, av, b_dp, delta_so, Eg):
    """Compute Bir-Pikus band-edge shifts for biaxial [001] strain.

    Returns dict with all band-edge shifts in eV, including the LH-SO
    coupled eigenvalues from the off-diagonal QT2 coupling.

    Sign convention matches compute_bp_scalar in strain_solver.f90.
    """
    eps_xx = (a_sub - a0) / a0
    eps_yy = eps_xx
    eps_zz = -2.0 * C12 / C11 * eps_xx

    Tr_eps = eps_xx + eps_yy + eps_zz
    P_eps = -av * Tr_eps
    Q_eps = b_dp / 2.0 * (eps_zz - 0.5 * (eps_xx + eps_yy))

    delta_Ec  = ac * Tr_eps
    delta_EHH = -P_eps + Q_eps
    delta_ELH = -P_eps - Q_eps
    delta_ESO = -P_eps

    # LH-SO off-diagonal coupling
    QT2 = 2.0 * Q_eps
    coupling = QT2 / np.sqrt(2.0)

    # Strained band edges (diagonal only, before LH-SO mixing)
    E_SO_unstrained = -delta_so
    E_HH_unstrained = 0.0
    E_LH_unstrained = 0.0
    E_CB_unstrained = Eg

    E_HH_strained = E_HH_unstrained + delta_EHH
    E_CB_strained = E_CB_unstrained + delta_Ec

    # LH-SO coupled 2x2 system
    a = E_LH_unstrained + delta_ELH
    b = E_SO_unstrained + delta_ESO
    c = abs(coupling)
    mid = (a + b) / 2.0
    half_diff = (a - b) / 2.0
    root = np.sqrt(half_diff**2 + c**2)
    E_LHSO_low  = mid - root
    E_LHSO_high = mid + root

    HH_LH_splitting = E_HH_strained - E_LHSO_high

    return {
        "eps_xx": eps_xx, "eps_yy": eps_yy, "eps_zz": eps_zz,
        "Tr_eps": Tr_eps, "P_eps": P_eps, "Q_eps": Q_eps,
        "delta_Ec": delta_Ec, "delta_EHH": delta_EHH,
        "delta_ELH": delta_ELH, "delta_ESO": delta_ESO,
        "E_CB": E_CB_strained, "E_HH": E_HH_strained,
        "E_LHSO_low": E_LHSO_low, "E_LHSO_high": E_LHSO_high,
        "HH_LH_splitting": HH_LH_splitting,
        "QT2": QT2, "coupling": coupling,
    }


# Inline configs for InAs bulk
CONFIG_UNSTRAINED = (
    "waveVector: k0\n"
    "waveVectorMax: 0\n"
    "waveVectorStep: 1\n"
    "confinement:  0\n"
    "FDstep: 101\n"
    "FDorder: 2\n"
    "numLayers:  1\n"
    "material1: InAs\n"
    "numcb: 2\n"
    "numvb: 6\n"
    "ExternalField: 0  EF\n"
    "EFParams: 0.0\n"
)

CONFIG_STRAINED = (
    "waveVector: k0\n"
    "waveVectorMax: 0\n"
    "waveVectorStep: 1\n"
    "confinement:  0\n"
    "FDstep: 101\n"
    "FDorder: 2\n"
    "numLayers:  1\n"
    "material1: InAs\n"
    "numcb: 2\n"
    "numvb: 6\n"
    "ExternalField: 0  EF\n"
    "EFParams: 0.0\n"
    f"strainSubstrate: {GAAS_A0_SUB}\n"
)


def run_with_config(build_dir, config_str, label):
    """Run bandStructure with an inline config string. Returns eigenvalues."""
    work = tempfile.mkdtemp(prefix=f"rung5_{label}_")
    try:
        # Write to a staging name; run_exe copies src -> work_dir/input.cfg
        cfg_path = os.path.join(work, "strain.cfg")
        with open(cfg_path, "w") as f:
            f.write(config_str)
        rc, output_dir = run_exe(build_dir, "bandStructure", cfg_path, work)
        if rc != 0:
            print(f"  FATAL: bandStructure returned {rc} for {label}")
            sys.exit(1)
        eig_path = os.path.join(output_dir, "eigenvalues.dat")
        if not os.path.isfile(eig_path):
            print(f"  FATAL: eigenvalues.dat not found for {label}")
            sys.exit(1)
        data = parse_eigenvalues(eig_path)
        if not data:
            print(f"  FATAL: no eigenvalue data for {label}")
            sys.exit(1)
        return data[0][1]  # eigenvalues at first k-point
    finally:
        shutil.rmtree(work, ignore_errors=True)


def test_r1_strained_eigenvalues(build_dir, bp):
    """R1: Strained bulk eigenvalues match analytical Bir-Pikus with LH-SO mixing."""
    print("  [R1] Strained eigenvalue match (LH-SO mixed)")
    evals = run_with_config(build_dir, CONFIG_STRAINED, "strained")

    # Ascending order with LH-SO mixing:
    #   bands 0,1 = LHSO_low; bands 2,3 = HH; bands 4,5 = LHSO_high; bands 6,7 = CB
    # Under compressive strain, HH shifts below the upper LH-SO mixed state.
    expected_map = {
        "LHSO_low":  (0, bp["E_LHSO_low"]),
        "HH":        (2, bp["E_HH"]),
        "LHSO_high": (4, bp["E_LHSO_high"]),
        "CB":        (6, bp["E_CB"]),
    }

    all_pass = True
    for label, (idx, expected) in expected_map.items():
        computed = evals[idx]
        passed, delta, _ = compare_value(
            computed, expected, TOL_STRAIN,
            f"  {label}", unit="eV",
        )
        status = "PASS" if passed else "FAIL"
        print(f"    {label:>10s}: computed={computed:+.6f}  "
              f"expected={expected:+.6f}  delta={delta:.4e}  {status}")
        all_pass = all_pass and passed

    return all_pass, evals


def test_r2_hh_lh_splitting(evals, bp):
    """R2: HH-LH splitting matches analytical value with LH-SO mixing."""
    print("  [R2] HH-LH splitting (with LH-SO mixing)")
    # Ascending: LHSO_low(0,1), HH(2,3), LHSO_high(4,5), CB(6,7)
    computed_splitting = evals[2] - evals[4]  # HH - LHSO_high
    expected_splitting = bp["HH_LH_splitting"]
    passed, delta, _ = compare_value(
        computed_splitting, expected_splitting, TOL_STRAIN,
        "  HH-LH splitting", unit="eV",
    )
    status = "PASS" if passed else "FAIL"
    print(f"    computed={computed_splitting:+.6f}  "
          f"expected={expected_splitting:+.6f}  delta={delta:.4e}  {status}")
    return passed


def test_r3_additive_modification(build_dir, bp):
    """R3: Strained minus unstrained equals Bir-Pikus shift exactly."""
    print("  [R3] Strain as additive modification")
    evals_unref = run_with_config(build_dir, CONFIG_UNSTRAINED, "unstrained")
    evals_strained = run_with_config(build_dir, CONFIG_STRAINED, "strained")

    all_pass = True

    # Unstrained (ascending): SO(0,1)=-Delta, HH/LH(2-5)=0.0, CB(6,7)=Eg
    # Strained (ascending):   LHSO_low(0,1), HH(2,3), LHSO_high(4,5), CB(6,7)

    # CB: no LH-SO mixing, shift should be exact (both at index 6)
    cb_shift = evals_strained[6] - evals_unref[6]
    expected_cb = bp["delta_Ec"]
    passed_cb, delta_cb, _ = compare_value(
        cb_shift, expected_cb, TOL_ADDITIVE,
        "  CB shift", unit="eV",
    )
    status_cb = "PASS" if passed_cb else "FAIL"
    print(f"    CB shift: {cb_shift:+.10f}  expected={expected_cb:+.6f}  "
          f"delta={delta_cb:.2e}  {status_cb}")
    all_pass = all_pass and passed_cb

    # HH: no LH-SO mixing. Strained HH at index 2, unstrained HH at index 4 (both 0.0)
    hh_shift = evals_strained[2] - evals_unref[4]
    expected_hh = bp["delta_EHH"]
    passed_hh, delta_hh, _ = compare_value(
        hh_shift, expected_hh, TOL_ADDITIVE,
        "  HH shift", unit="eV",
    )
    status_hh = "PASS" if passed_hh else "FAIL"
    print(f"    HH shift: {hh_shift:+.10f}  expected={expected_hh:+.6f}  "
          f"delta={delta_hh:.2e}  {status_hh}")
    all_pass = all_pass and passed_hh

    # LH+SO block: off-diagonal QT2 coupling preserves the trace.
    # Strained LHSO = indices 0,1 (LHSO_low) + 4,5 (LHSO_high)
    # Unstrained LH+SO = indices 0,1 (SO) + 2,3 (LH)
    lhsstrained_sum = (evals_strained[0] + evals_strained[1] +
                       evals_strained[4] + evals_strained[5])
    expected_block_sum = 2 * bp["E_LHSO_low"] + 2 * bp["E_LHSO_high"]
    passed_block, delta_block, _ = compare_value(
        lhsstrained_sum, expected_block_sum, TOL_ADDITIVE,
        "  LH+SO block sum", unit="eV",
    )
    status_block = "PASS" if passed_block else "FAIL"
    print(f"    LH+SO block sum: {lhsstrained_sum:+.10f}  "
          f"expected={expected_block_sum:+.6f}  delta={delta_block:.2e}  {status_block}")
    all_pass = all_pass and passed_block

    return all_pass


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        sys.exit(1)

    build_dir = os.path.abspath(sys.argv[1])

    exe_path = os.path.join(build_dir, "src", "bandStructure")
    if not os.path.isfile(exe_path):
        print(f"ERROR: executable not found: {exe_path}")
        sys.exit(1)

    print("=" * 60)
    print("  RUNG 5: Bulk Strain Verification (InAs/GaAs)")
    print("=" * 60)
    print()

    # Compute analytical reference values
    bp = bir_pikus_bulk_inas(INAS_A0, GAAS_A0_SUB, INAS_C11, INAS_C12,
                              INAS_AC, INAS_AV, INAS_B, INAS_DELTA, INAS_EG)

    print(f"  InAs a0 = {INAS_A0} A, GaAs substrate a_sub = {GAAS_A0_SUB} A")
    print(f"  eps_xx = eps_yy = {bp['eps_xx']:.6f}")
    print(f"  eps_zz           = {bp['eps_zz']:.6f}")
    print(f"  Tr(eps)          = {bp['Tr_eps']:.6f}")
    print(f"  P_eps            = {bp['P_eps']:+.6f} eV")
    print(f"  Q_eps            = {bp['Q_eps']:+.6f} eV")
    print(f"  delta_Ec         = {bp['delta_Ec']:+.6f} eV")
    print(f"  delta_EHH        = {bp['delta_EHH']:+.6f} eV")
    print(f"  delta_ELH        = {bp['delta_ELH']:+.6f} eV")
    print(f"  delta_ESO        = {bp['delta_ESO']:+.6f} eV")
    print(f"  E_CB             = {bp['E_CB']:+.6f} eV")
    print(f"  E_HH             = {bp['E_HH']:+.6f} eV")
    print(f"  E_LHSO_low       = {bp['E_LHSO_low']:+.6f} eV")
    print(f"  E_LHSO_high      = {bp['E_LHSO_high']:+.6f} eV")
    print(f"  HH-LH splitting  = {bp['HH_LH_splitting']:+.6f} eV")
    print()

    all_pass = True

    # R1
    r1_pass, evals_strained = test_r1_strained_eigenvalues(build_dir, bp)
    all_pass = all_pass and r1_pass
    print()

    # R2
    r2_pass = test_r2_hh_lh_splitting(evals_strained, bp)
    all_pass = all_pass and r2_pass
    print()

    # R3
    r3_pass = test_r3_additive_modification(build_dir, bp)
    all_pass = all_pass and r3_pass
    print()

    # Summary
    print("=" * 60)
    if all_pass:
        print("  PASS: all bulk strain verification checks passed")
        sys.exit(0)
    else:
        print("  FAIL: one or more bulk strain checks failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
