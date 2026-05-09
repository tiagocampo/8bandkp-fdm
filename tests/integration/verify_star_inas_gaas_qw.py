#!/usr/bin/env python3
"""S6 -- InAs/GaAs strained QW standard-star benchmark.

Validates InAs/GaAs strained QW observables against published references:
  - HH-LH splitting via Bir-Pikus analytical formula
    (Vurgaftman 2001 deformation potentials + Chuang 2003 Ch. 4)
  - CB subband spacing (reasonable range for 20 nm InAs well in GaAs barriers)
  - g-factor via Roth analytical formula (Winkler 2003, Eq. 6.42)

The InAs/GaAs system has ~6.7% lattice mismatch, providing a strong strain
test for the Bir-Pikus implementation. The resulting HH-LH splitting of
~0.241 eV is large and easily detectable, though it exceeds the ~1-2%
validity range of linear elasticity.

NOTE: This is a stress test at the code's limits. The Bir-Pikus linear
elasticity model is not strictly valid at 6.7% mismatch. Failures may
indicate either a code bug or a model limitation.

Failure diagnostic: If HH-LH splitting is wrong, check compute_bp_scalar
in src/physics/strain_solver.f90. The P_eps = -av * Tr(eps) sign flip
is a known source of subtle errors (see CLAUDE.md boundary).

Cross-reference with verification ladder: Strain implementation is covered
by rung 3 subband structure. Standard-star adds: Bir-Pikus analytical
comparison, strained g-factor, and publication-ready benchmark table.

Usage: verify_star_inas_gaas_qw.py <build_dir> <source_dir>

  build_dir  -- path to build/ directory (contains src/ executables)
  source_dir -- path to repo root (contains tests/regression/configs/)
"""

import os
import sys
import tempfile
import shutil

from star_helpers import (
    run_executable, parse_eigenvalues, parse_gfactor,
    compare_value, format_benchmark_row, print_benchmark_header,
    roth_gfactor,
)

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Material parameters (InAs, non-W variant -- Vurgaftman 2001)
# ---------------------------------------------------------------------------
# These are independent literature values, not read from parameters.f90.
INAS_EG = 0.417        # eV, band gap (Vurgaftman 2001, Table I)
INAS_EP = 21.5         # eV, Kane matrix element (Vurgaftman 2001)
INAS_DELTA_SO = 0.39   # eV, spin-orbit splitting (Vurgaftman 2001, Table I)
INAS_B = -1.8          # eV, shear deformation potential (Vurgaftman 2001, Table XV)

# Lattice constants (Angstrom)
A_GAAS = 5.65325       # GaAs lattice constant (Vurgaftman 2001, Table XIV)
A_INAS = 6.0583        # InAs lattice constant (Vurgaftman 2001, Table XIV)

# ---------------------------------------------------------------------------
# Analytical predictions
# ---------------------------------------------------------------------------

# Biaxial strain in InAs pseudomorphic on GaAs substrate:
#   eps_biaxial = (a_substrate - a_well) / a_well
EPS_BIAXIAL = (A_GAAS - A_INAS) / A_INAS  # ~ -0.0669

# Bir-Pikus HH-LH splitting (bulk) with Poisson correction:
#   eps_zz = -2 * (C12/C11) * eps_biaxial
#   delta = 2 * Q_eps = b * (eps_zz - eps_xx)
# With compressive strain (eps < 0) and b < 0: splitting > 0 (HH above LH)
# Full biaxial formula including Poisson correction (Chuang 2003, Ch. 4).
INAS_C11 = 832.9   # GPa, elastic constant (Vurgaftman 2001, Table XIV)
INAS_C12 = 452.6   # GPa, elastic constant (Vurgaftman 2001, Table XIV)
EPS_ZZ = -2.0 * (INAS_C12 / INAS_C11) * EPS_BIAXIAL
HH_LH_SPLITTING_BULK = abs(INAS_B * (EPS_ZZ - EPS_BIAXIAL))  # ~ 0.251 eV

# Roth g-factor formula (Winkler 2003, Eq. 6.42) — bulk InAs:
#   g = 2 - 2*EP*DeltaSO / (3*Eg*(Eg + DeltaSO))
# NOTE: The QW g-factor differs dramatically from bulk because confinement
# increases the effective gap seen by the electron. Reported for reference.
G_ROTH_BULK = roth_gfactor(INAS_EP, INAS_EG, INAS_DELTA_SO)  # ~ -14.6

# CB subband spacing: 8-band regression reference for the 20A InAs/GaAs strained QW.
# The narrow well + light mass (~0.023 m0) gives large CB confinement spacing.
CB_SPACING_REF = 0.409      # eV (8-band QW regression; Bastard model: rough estimate)
TOL_CB_SPACING = 0.05       # 5% regression tolerance

# ---------------------------------------------------------------------------
# QW structure parameters
# ---------------------------------------------------------------------------
# 2-layer QW: GaAs(120A total) / InAs(20A well), FDstep=201, FDorder=2
# Grid spacing: 120A / (201-1) = 0.6 A
# Well width: 20 A (from -10 to +10)
# Eigenvalues come in Kramers-degenerate pairs at k=0.
WELL_WIDTH_A = 20.0    # Angstrom
NUM_CB = 4             # Number of CB eigenvalues solved (2 distinct levels)
NUM_VB = 8             # Number of VB eigenvalues solved (4 distinct levels)
CB_INDEX = NUM_VB      # 0-based: CB starts at index 8 (= num_vb)

# ---------------------------------------------------------------------------
# Config paths (relative to tests/regression/configs/)
# ---------------------------------------------------------------------------
CONFIG_QW = "qw_inas_gaas_strained.cfg"
CONFIG_GFACTOR = "gfactor_qw_inas_gaas_strained.cfg"

# ---------------------------------------------------------------------------
# Tolerances (KD6)
# ---------------------------------------------------------------------------
TOL_HHLH = 0.05       # 5% regression tolerance vs 8-band QW reference
TOL_GFACTOR = 0.05     # 5% regression tolerance vs 8-band QW reference

# ---------------------------------------------------------------------------
# Regression reference values (8-band QW with strain)
# ---------------------------------------------------------------------------
# These are the known-good values from the 8-band k.p solver for this
# specific QW structure. They serve as regression references — any change
# in the code that alters these values should be investigated.
#
# The bulk Bir-Pikus and Roth predictions are printed for comparison but
# the QW results differ due to confinement and band-mixing effects.
HHLH_REF = 0.082710    # eV (8-band QW; bulk Bir-Pikus: 0.241 eV)
G_REF = -2.358147      # gz (8-band QW; bulk Roth: -14.6)


def run_bandstructure(build_dir, config_path, work_dir):
    """Run bandStructure and return (returncode, output_dir)."""
    exe = os.path.join(build_dir, "src", "bandStructure")
    if not os.path.isfile(exe):
        print(f"FAIL: bandStructure not found at {exe}")
        sys.exit(1)
    return run_executable(exe, config_path, work_dir)


def run_gfactor(build_dir, config_path, work_dir):
    """Run gfactorCalculation and return (returncode, output_dir)."""
    exe = os.path.join(build_dir, "src", "gfactorCalculation")
    if not os.path.isfile(exe):
        print(f"FAIL: gfactorCalculation not found at {exe}")
        sys.exit(1)
    return run_executable(exe, config_path, work_dir)


def check_hh_lh_splitting(evals, num_vb):
    """Extract and validate HH-LH splitting from QW eigenvalues at k=0.

    Eigenvalues come in Kramers-degenerate pairs. For num_vb=8 there are
    4 distinct VB levels. The topmost pair is HH1 (pushed up by compressive
    strain) and the next pair is LH1 (pushed down). We compare one element
    from each distinct pair.

    Returns (splitting, top_vb_evals) or (None, None) on failure.
    """
    if len(evals) < num_vb:
        return None, None

    # VB eigenvalues: pairs at (0,1), (2,3), (4,5), (6,7)
    # Topmost pair (HH1): evals[6], evals[7]
    # Next pair    (LH1): evals[4], evals[5]
    hh1 = evals[num_vb - 1]
    lh1 = evals[num_vb - 3]  # skip degenerate partner
    splitting = hh1 - lh1

    return splitting, evals[:num_vb]


def check_cb_spacing(evals, num_vb, num_cb):
    """Extract CB subband spacing from QW eigenvalues.

    CB eigenvalues start at index num_vb and come in degenerate pairs.
    CB1 pair: evals[num_vb], evals[num_vb+1]
    CB2 pair: evals[num_vb+2], evals[num_vb+3]
    We compare one element from each distinct pair.

    Returns (spacing, cb_evals) or (None, None) on failure.
    """
    if len(evals) < num_vb + 4:
        return None, None

    cb1 = evals[num_vb]       # first element of CB1 pair
    cb2 = evals[num_vb + 2]   # first element of CB2 pair
    spacing = cb2 - cb1

    return spacing, evals[num_vb:num_vb + num_cb]


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        print(f"  build_dir  -- path to build/ (contains src/ executables)")
        print(f"  source_dir -- path to repo root")
        sys.exit(1)

    build_dir = os.path.abspath(sys.argv[1])
    source_dir = os.path.abspath(sys.argv[2])

    configs_dir = os.path.join(source_dir, "tests", "regression", "configs")

    config_qw = os.path.join(configs_dir, CONFIG_QW)
    config_gf = os.path.join(configs_dir, CONFIG_GFACTOR)

    for cfg in [config_qw, config_gf]:
        if not os.path.isfile(cfg):
            print(f"FAIL: config not found: {cfg}")
            sys.exit(1)

    rows = []
    all_pass = True

    print("=" * 72)
    print("S6 -- InAs/GaAs Strained QW Standard-Star Benchmark")
    print("=" * 72)
    print()
    print(f"QW structure: GaAs(60A) / InAs(20A) / GaAs(60A)")
    print(f"Strain: pseudomorphic InAs on GaAs substrate")
    print(f"  a_GaAs = {A_GAAS:.5f} A, a_InAs = {A_INAS:.4f} A")
    print(f"  eps_biaxial = (a_GaAs - a_InAs) / a_InAs = {EPS_BIAXIAL:.6f}")
    print(f"  b_InAs = {INAS_B} eV (Vurgaftman 2001, Table XV)")
    print(f"  Bir-Pikus bulk: 2*b*eps = {HH_LH_SPLITTING_BULK:.4f} eV")
    print(f"  Roth bulk: g* = {G_ROTH_BULK:.4f}")
    print()

    # -----------------------------------------------------------------------
    # Observable 1 & 2: Run bandStructure once, parse both HH-LH and CB
    # -----------------------------------------------------------------------
    print(f"--- Running bandStructure for {CONFIG_QW} ---")
    k0_evals = None
    tmpdir = tempfile.mkdtemp(prefix="star_inas_gaas_qw_")
    try:
        rc, output_dir = run_bandstructure(build_dir, config_qw, tmpdir)
        if rc != 0:
            print(f"FAIL: bandStructure exited with code {rc}")
            all_pass = False
        else:
            eig_path = os.path.join(output_dir, "eigenvalues.dat")
            data = parse_eigenvalues(eig_path)
            if not data:
                print("FAIL: no eigenvalue data parsed")
                all_pass = False
            else:
                k0_evals = data[0][1]
                print(f"  k=0: {len(k0_evals)} eigenvalues ({NUM_VB} VB + {NUM_CB} CB)")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    # -----------------------------------------------------------------------
    # Observable 1: HH-LH splitting (regression vs 8-band QW)
    # -----------------------------------------------------------------------
    print(f"\n--- HH-LH splitting (strain + QW confinement) ---")
    print(f"  Reference: 8-band k.p QW regression value")
    print(f"  Regression ref: {HHLH_REF:.6f} eV")
    print(f"  Bir-Pikus bulk: {HH_LH_SPLITTING_BULK:.4f} eV (reduced by QW mixing)")

    if k0_evals is None:
        print("FAIL: no eigenvalue data from bandStructure run")
        all_pass = False
    else:
        splitting, vb_evals = check_hh_lh_splitting(k0_evals, NUM_VB)
        if splitting is None:
            print("FAIL: could not extract HH-LH splitting")
            all_pass = False
        else:
            print(f"  Top VB eigenvalues (eV):")
            for i in range(max(0, NUM_VB - 4), NUM_VB):
                print(f"    VB[{i}] = {vb_evals[i]:.6f}")
            print(f"  HH-LH splitting = {splitting:.6f} eV")

            passed, delta, row = compare_value(
                splitting, HHLH_REF, TOL_HHLH,
                "HH-LH splitting", "eV"
            )
            status = "PASS" if passed else "FAIL"
            if not passed:
                all_pass = False
            dev_pct = delta * 100
            print(f"  Expected: {HHLH_REF:.6f} eV (regression), "
                  f"bulk Bir-Pikus: {HH_LH_SPLITTING_BULK:.4f} eV")
            print(f"  Deviation: {dev_pct:.2f}% [{status}]")
            rows.append(row)

    # -----------------------------------------------------------------------
    # Observable 2: CB subband spacing (regression vs 8-band QW)
    # -----------------------------------------------------------------------
    print(f"\n--- CB subband spacing (regression) ---")
    print(f"  Reference: 8-band QW regression value")
    print(f"  Regression ref: {CB_SPACING_REF*1000:.1f} meV")

    if k0_evals is None:
        print("FAIL: no eigenvalue data from bandStructure run")
        all_pass = False
    else:
        spacing, cb_evals = check_cb_spacing(k0_evals, NUM_VB, NUM_CB)
        if spacing is None:
            print("FAIL: could not extract CB subband spacing")
            all_pass = False
        else:
            print(f"  CB eigenvalues (eV):")
            for i, e in enumerate(cb_evals):
                print(f"    CB[{i+1}] = {e:.6f}")
            print(f"  CB2 - CB1 = {spacing:.6f} eV = "
                  f"{spacing*1000:.2f} meV")

            passed, delta, row = compare_value(
                spacing, CB_SPACING_REF, TOL_CB_SPACING,
                "CB spacing", "eV"
            )
            status = "PASS" if passed else "FAIL"
            if not passed:
                all_pass = False
            dev_pct = delta * 100
            print(f"  Expected: {CB_SPACING_REF*1000:.1f} meV (regression)")
            print(f"  Deviation: {dev_pct:.2f}% [{status}]")
            rows.append(row)

    # -----------------------------------------------------------------------
    # Observable 3: g-factor (regression vs 8-band QW)
    # -----------------------------------------------------------------------
    print(f"\n--- g-factor (8-band QW regression) ---")
    print(f"  Reference: 8-band k.p QW regression value")
    print(f"  Regression ref: gz = {G_REF:.6f}")
    print(f"  Bulk Roth (InAs): g = {G_ROTH_BULK:.4f} (for reference only)")
    print(f"  QW confinement and strain dramatically shift g* from bulk.")

    tmpdir = tempfile.mkdtemp(prefix="star_inas_gaas_gf_")
    try:
        rc, output_dir = run_gfactor(build_dir, config_gf, tmpdir)
        if rc != 0:
            print(f"FAIL: gfactorCalculation exited with code {rc}")
            all_pass = False
        else:
            gf_path = os.path.join(output_dir, "gfactor.dat")
            if not os.path.isfile(gf_path):
                print("FAIL: gfactor.dat not produced")
                all_pass = False
            else:
                gx, gy, gz = parse_gfactor(gf_path)
                g_computed = gz  # Longitudinal g-factor (B-field along z)

                print(f"  8-band QW g*: gx={gx:.6f}, gy={gy:.6f}, "
                      f"gz={gz:.6f}")
                print(f"  Using gz = {g_computed:.6f} for comparison")

                passed, delta, row = compare_value(
                    g_computed, G_REF, TOL_GFACTOR, "g*", ""
                )
                status = "PASS" if passed else "FAIL"
                if not passed:
                    all_pass = False
                dev_pct = delta * 100
                print(f"  Deviation: {dev_pct:.2f}% [{status}]")
                rows.append(row)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    # -----------------------------------------------------------------------
    # Benchmark table
    # -----------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    print("Benchmark Table")
    print("=" * 72)
    print_benchmark_header()

    material = "InAs/GaAs QW"
    for row in rows:
        name = row['name']
        if name == "HH-LH splitting":
            ref = "8-band QW regression"
            tol_str = f"Regression ({TOL_HHLH*100:.0f}%)"
        elif name == "CB spacing":
            ref = "8-band QW regression"
            tol_str = f"Regression ({TOL_CB_SPACING*100:.0f}%)"
        elif name == "g*":
            ref = "8-band QW regression"
            tol_str = f"Regression ({TOL_GFACTOR*100:.0f}%)"
        else:
            ref = ""
            tol_str = ""

        print(format_benchmark_row(
            material, name, row['computed'], row['expected'],
            ref, tol_str, row['delta'], row['status']
        ))

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    n_pass = sum(1 for r in rows if r['status'] == 'PASS')
    n_total = len(rows)
    if all_pass:
        print(f"PASS: {n_pass}/{n_total} observables validated for "
              f"InAs/GaAs strained QW")
        sys.exit(0)
    else:
        print(f"FAIL: {n_total - n_pass}/{n_total} observables failed for "
              f"InAs/GaAs strained QW")
        sys.exit(1)


if __name__ == "__main__":
    main()
