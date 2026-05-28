#!/usr/bin/env python3
"""S3 -- InSb bulk standard-star benchmark.

Validates InSb bulk observables against published literature:
  - Eg = 0.235 eV          (Vurgaftman 2001, Table I)
  - DeltaSO = 0.810 eV     (Vurgaftman 2001, Table I)
  - m*_e via Kane formula   (2-band: Eg/(EP+Eg))
  - g* via Roth formula     (Winkler 2003, Eq. 6.42)

InSb has extreme spin-orbit coupling (DeltaSO = 0.81 eV) and a narrow gap
(Eg = 0.235 eV), making it the most stringent g-factor test in the suite.
The Roth prediction g ~ -49 is far from the free-electron value g = 2.

Cross-reference with verification ladder:
  - Eg and DeltaSO covered by rung 1 (R1).
  - Effective mass covered by rung 2 (R6-R7).
  Standard-star adds: g-factor validation and publication-ready benchmark table.

Usage: verify_star_insb_bulk.py <build_dir> <source_dir>

  build_dir  -- path to build/ directory (contains src/ executables)
  source_dir -- path to repo root (contains tests/regression/configs/)
"""

# COVERAGE: observable=Eg geometry=bulk material=InSb ref=Vurgaftman2001
# COVERAGE: observable=Delta_SO geometry=bulk material=InSb ref=Vurgaftman2001
# COVERAGE: observable=m*_e geometry=bulk material=InSb ref=Kane
# COVERAGE: observable=g*_cb geometry=bulk material=InSb ref=Winkler2003
import os
import sys
import tempfile
import shutil

try:
    import numpy as np
except ImportError:
    print("FAIL: numpy is required (pip install numpy)")
    sys.exit(1)

from star_helpers import (
    run_exe, parse_eigenvalues, parse_gfactor,
    compare_value, format_benchmark_row, print_benchmark_header,
    HBAR2_OVER_2M0, roth_gfactor, TOL_EXACT,
    extract_effective_mass,
)

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# InSb material parameters
# Sources: Vurgaftman 2001 (Eg, DeltaSO), parameters.f90 (EP)
# ---------------------------------------------------------------------------
EG_INSB = 0.235       # eV, Vurgaftman 2001 Table I
DELTA_SO = 0.810      # eV, Vurgaftman 2001 Table I
EP_INSB = 23.3        # eV, parameters.f90

# 2-band Kane formula: m*_kane = Eg / (EP + Eg) (informational)
M_STAR_KANE = EG_INSB / (EP_INSB + EG_INSB)  # ~0.00997 m0

# 8-band model reference (validated against kdotpy)
M_STAR_8BAND = 0.0069  # m0

# Roth g-factor formula (Winkler 2003, Eq. 6.42):
#   g_roth = 2 - 2*EP*DeltaSO / (3*Eg*(Eg + DeltaSO))
G_ROTH = roth_gfactor(EP_INSB, EG_INSB, DELTA_SO)

# ---------------------------------------------------------------------------
# Tolerances (KD6)
# ---------------------------------------------------------------------------
TOL_MASS = 0.05        # 5%: 8-band reference comparison
TOL_GFACTOR = 0.01     # 1%: Roth g-factor vs 8-band

# CB eigenvalue index (0-based): bands 1-4 valence, 5-6 SO, 7-8 CB
CB_INDEX = 6


def check_eigenvalues_at_k0(evals):
    """Check Eg and DeltaSO at k=0 against parameters.f90.

    At k=0, eigenvalues in ascending order:
      [-DeltaSO, -DeltaSO, 0, 0, 0, 0, Eg, Eg]

    Returns list of (passed, row_dict) tuples.
    """
    results = []

    # Eg: eigenvalue at CB_INDEX
    eg_computed = evals[CB_INDEX]
    passed_eg, delta_eg, row_eg = compare_value(
        eg_computed, EG_INSB, TOL_EXACT, "Eg", "eV"
    )
    results.append((passed_eg, row_eg))

    # DeltaSO: eigenvalue at index 0 (should be -DeltaSO)
    delta_so_computed = -evals[0]
    passed_so, delta_so, row_so = compare_value(
        delta_so_computed, DELTA_SO, TOL_EXACT, "DeltaSO", "eV"
    )
    results.append((passed_so, row_so))

    return results


def _extract_mass(eig_path):
    """Extract CB effective mass via shared adaptive parabolic fit."""
    result = extract_effective_mass(eig_path, cb_index=CB_INDEX)
    if result is None:
        return None
    return result[0]  # m_star


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        print(f"  build_dir  -- path to build/ (contains src/ executables)")
        print(f"  source_dir -- path to repo root")
        sys.exit(1)

    build_dir = os.path.abspath(sys.argv[1])
    source_dir = os.path.abspath(sys.argv[2])

    configs_dir = os.path.join(source_dir, "tests", "regression", "configs")

    config_k0 = os.path.join(configs_dir, "bulk_insb_k0.toml")
    config_disp = os.path.join(configs_dir, "bulk_insb_kx_dispersion.toml")
    config_gf = os.path.join(configs_dir, "gfactor_bulk_insb_cb.toml")

    for cfg in [config_k0, config_disp, config_gf]:
        if not os.path.isfile(cfg):
            print(f"FAIL: config not found: {cfg}")
            sys.exit(1)

    rows = []
    all_pass = True

    print("=" * 72)
    print("S3 -- InSb Bulk Standard-Star Benchmark")
    print("=" * 72)

    # -----------------------------------------------------------------------
    # Observable 1 & 2: Eg and DeltaSO at k=0
    # -----------------------------------------------------------------------
    print(f"\n--- Eg and DeltaSO (k=0) ---")
    tmpdir = tempfile.mkdtemp(prefix="star_insb_k0_")
    try:
        rc, output_dir = run_exe(build_dir, 'bandStructure', config_k0, tmpdir)
        if rc != 0:
            print(f"FAIL: bandStructure exited with code {rc}")
            all_pass = False
        else:
            eig_path = os.path.join(output_dir, "eigenvalues.dat")
            data = parse_eigenvalues(eig_path)
            if not data:
                print("FAIL: no eigenvalue data at k=0")
                all_pass = False
            else:
                evals = data[0][1]
                print(f"  Eigenvalues at k=0: {['%.6f' % e for e in evals]}")

                eig_results = check_eigenvalues_at_k0(evals)
                for passed, row in eig_results:
                    status = "PASS" if passed else "FAIL"
                    if not passed:
                        all_pass = False
                    name = row['name']
                    print(f"  {name}: computed={row['computed']:.6f}, "
                          f"expected={row['expected']:.6f}, "
                          f"delta={row['delta']:.2e} [{status}]")
                    rows.append(row)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    # -----------------------------------------------------------------------
    # Observable 3: Effective mass via parabolic fitting
    # -----------------------------------------------------------------------
    print(f"\n--- Effective mass (Kane) ---")
    tmpdir = tempfile.mkdtemp(prefix="star_insb_disp_")
    try:
        rc, output_dir = run_exe(build_dir, 'bandStructure', config_disp, tmpdir)
        if rc != 0:
            print(f"FAIL: bandStructure exited with code {rc}")
            all_pass = False
        else:
            eig_path = os.path.join(output_dir, "eigenvalues.dat")
            m_star = _extract_mass(eig_path)
            if m_star is None:
                print("FAIL: could not extract effective mass")
                all_pass = False
            else:
                print(f"  m* (numerical): {m_star:.6f} m0")
                print(f"  m* (8-band ref): {M_STAR_8BAND:.6f} m0 "
                      f"(Kane: {M_STAR_KANE:.6f})")
                passed, delta, row = compare_value(
                    m_star, M_STAR_8BAND, TOL_MASS, "m*_e", "m0"
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
    # Observable 4: g-factor via Roth formula
    # -----------------------------------------------------------------------
    print(f"\n--- g-factor (Roth vs 8-band) ---")
    print(f"  Roth prediction: g = 2 - 2*EP*DeltaSO/(3*Eg*(Eg+DeltaSO))")
    print(f"                 = 2 - 2*{EP_INSB}*{DELTA_SO}/"
          f"(3*{EG_INSB}*{EG_INSB + DELTA_SO:.3f})")
    print(f"                 = {G_ROTH:.4f}")

    tmpdir = tempfile.mkdtemp(prefix="star_insb_gf_")
    try:
        rc, output_dir = run_exe(build_dir, 'gfactorCalculation', config_gf, tmpdir)
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
                g_computed = gz  # B-field along z for bulk

                print(f"  8-band g*: gx={gx:.6f}, gy={gy:.6f}, gz={gz:.6f}")
                print(f"  Using gz = {g_computed:.6f} for comparison")

                passed, delta, row = compare_value(
                    g_computed, G_ROTH, TOL_GFACTOR, "g*", ""
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

    material = "InSb"
    for row in rows:
        tol_str = f"{row['tolerance']:.0e}" if row['tolerance'] < 0.001 else \
                  f"{row['tolerance']*100:.0f}%"
        print(format_benchmark_row(
            material, row['name'], row['computed'], row['expected'],
            "Vurgaftman 2001" if row['name'] in ("Eg", "DeltaSO") else
            "2-band Kane" if row['name'] == "m*_e" else
            "Winkler 2003 Eq.6.42",
            tol_str, row['delta'], row['status']
        ))

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    n_pass = sum(1 for r in rows if r['status'] == 'PASS')
    n_total = len(rows)
    if all_pass:
        print(f"PASS: {n_pass}/{n_total} observables validated for InSb bulk")
        sys.exit(0)
    else:
        print(f"FAIL: {n_total - n_pass}/{n_total} observables failed for InSb bulk")
        sys.exit(1)


if __name__ == "__main__":
    main()
