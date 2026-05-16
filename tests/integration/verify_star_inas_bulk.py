#!/usr/bin/env python3
"""S2 -- InAs bulk standard-star benchmark.

Validates InAs bulk observables against published literature references:
  - Eg at k=0 (Vurgaftman 2001)
  - Effective mass via parabolic fitting vs 2-band Kane formula
  - g-factor vs Roth analytical formula (Winkler 2003, InAsW variant)

Cross-reference with verification ladder: Eg covered by rung 1 (R1),
effective mass covered by rung 2 (R6-R7). Standard-star adds: g-factor
validation with Roth analytical comparison and publication-ready benchmark
table output.

Usage: verify_star_inas_bulk.py <build_dir> <source_dir>

  build_dir  -- path to build/ directory (contains src/bandStructure, src/gfactorCalculation)
  source_dir -- path to repo root (contains tests/regression/configs/)
"""

# COVERAGE: observable=Eg geometry=bulk material=InAs ref=Vurgaftman2001
# COVERAGE: observable=m*_e geometry=bulk material=InAs ref=Kane
# COVERAGE: observable=g*_cb geometry=bulk material=InAsW ref=Winkler2003
import os
import shutil
import sys
import tempfile

try:
    import numpy as np
except ImportError:
    print("FAIL: numpy is required (pip install numpy)")
    sys.exit(1)

from star_helpers import (
    run_executable,
    parse_eigenvalues,
    parse_gfactor,
    compare_value,
    format_benchmark_row,
    print_benchmark_header,
    HBAR2_OVER_2M0, roth_gfactor,
)

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Material parameters
# ---------------------------------------------------------------------------
# InAs (non-W) -- Vurgaftman 2001, Table I
# Used for Eg and effective mass benchmarks.
INAS_EG = 0.417       # eV, band gap
INAS_EP = 21.5        # eV, Kane matrix element
INAS_DELTA_SO = 0.390 # eV, spin-orbit splitting

# InAsW -- Winkler 2003 parameter set
# Used for g-factor benchmark (existing g-factor config specifies InAsW).
INASW_EG = 0.418       # eV
INASW_EP = 22.2        # eV
INASW_DELTA_SO = 0.38  # eV

# Kane effective mass prediction: m* = Eg / (EP + Eg)
M_STAR_KANE = INAS_EG / (INAS_EP + INAS_EG)  # 0.0190 m0

# Roth g-factor (Winkler 2003, Eq. 6.42):
#   g = 2 - 2*EP*DeltaSO / (3*Eg*(Eg + DeltaSO))
G_ROTH = roth_gfactor(INASW_EP, INASW_EG, INASW_DELTA_SO)

# ---------------------------------------------------------------------------
# Config paths (relative to tests/regression/configs/)
# ---------------------------------------------------------------------------
CONFIG_EG = "bulk_inas_k0.cfg"
CONFIG_DISPERSION = "bulk_inas_kx_dispersion.cfg"
CONFIG_GFACTOR = "gfactor_bulk_inasw_cb.cfg"

# ---------------------------------------------------------------------------
# Tolerances (KD6)
# ---------------------------------------------------------------------------
TOL_EXACT = 1e-12      # machine precision for Eg at k=0
TOL_KANE = 0.10        # 10% for Kane effective mass comparison
TOL_ROTH = 0.01        # ~1% for Roth g-factor vs 8-band


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def extract_cb_dispersion(data):
    """Extract CB band E(k) from parsed eigenvalue data.

    Returns (k_array, E_array) as numpy arrays, using the highest
    eigenvalue (CB top) at each k-point.
    """
    ks = []
    es = []
    for k, evals in data:
        ks.append(k)
        es.append(evals[-1])  # highest eigenvalue = CB
    ks = np.array(ks)
    es = np.array(es)
    idx = np.argsort(ks)
    return ks[idx], es[idx]


def extract_mass_numerical_derivative(k_vals, e_vals):
    """Extract effective mass from the first nonzero k-point.

    Uses the two-point formula: c2 = (E(k1) - E0) / k1^2
    where E(k) = E0 + c2*k^2 near Gamma.

    Returns (m_star, c2, E0) or None if insufficient data.
    """
    mask = k_vals >= 0
    k_pos = k_vals[mask]
    e_pos = e_vals[mask]

    if len(k_pos) < 2:
        return None

    e0 = e_pos[0]
    k1 = k_pos[1]
    e1 = e_pos[1]

    if k1 <= 0:
        return None

    c2 = (e1 - e0) / k1**2
    if c2 <= 0:
        return None

    m_star = HBAR2_OVER_2M0 / c2
    return m_star, c2, e0


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        sys.exit(1)

    build_dir = os.path.abspath(sys.argv[1])
    source_dir = os.path.abspath(sys.argv[2])

    exe_band = os.path.join(build_dir, "src", "bandStructure")
    exe_gfactor = os.path.join(build_dir, "src", "gfactorCalculation")
    config_dir = os.path.join(source_dir, "tests", "regression", "configs")

    # Check executables
    for exe, label in [(exe_band, "bandStructure"), (exe_gfactor, "gfactorCalculation")]:
        if not os.path.isfile(exe):
            print(f"FAIL: {label} executable not found at {exe}")
            sys.exit(1)

    # Check configs
    for cfg in [CONFIG_EG, CONFIG_DISPERSION, CONFIG_GFACTOR]:
        path = os.path.join(config_dir, cfg)
        if not os.path.isfile(path):
            print(f"FAIL: config not found: {path}")
            sys.exit(1)

    all_pass = True
    rows = []

    print("=" * 72)
    print("S2 -- InAs Bulk Standard-Star Benchmark")
    print("=" * 72)

    # ------------------------------------------------------------------
    # Observable 1: Eg at k=0 (exact, self-consistency check)
    # Reference: Vurgaftman 2001, Table I
    # ------------------------------------------------------------------
    print(f"\n{'─' * 72}")
    print("Observable 1: Eg (band gap at k=0)")
    print(f"{'─' * 72}")

    tmpdir = tempfile.mkdtemp(prefix="star_inas_eg_")
    try:
        config_path = os.path.join(config_dir, CONFIG_EG)
        rc, output_dir = run_executable(exe_band, config_path, tmpdir)
        if rc != 0:
            print(f"  FAIL: bandStructure returned exit code {rc}")
            all_pass = False
        else:
            eig_path = os.path.join(output_dir, "eigenvalues.dat")
            data = parse_eigenvalues(eig_path)

            if not data:
                print("  FAIL: no eigenvalue data parsed")
                all_pass = False
            else:
                k0_evals = data[0][1]
                # CB eigenvalue at index 6 (0-based) for bulk 8-band
                # Eigenvalues in ascending order: SO, SO, HH/LH x4, CB, CB
                eg_computed = k0_evals[6]

                print(f"  Computed Eg: {eg_computed:.12f} eV")
                print(f"  Expected Eg: {INAS_EG:.12f} eV (Vurgaftman 2001, Table I)")

                passed, delta, row = compare_value(
                    eg_computed, INAS_EG, TOL_EXACT,
                    "Eg", "eV"
                )

                bench_row = format_benchmark_row(
                    "InAs", "Eg", eg_computed, INAS_EG,
                    "Vurgaftman 2001, Table I",
                    f"Exact ({TOL_EXACT:.0e})", delta,
                    "PASS" if passed else "FAIL"
                )
                rows.append(bench_row)

                if passed:
                    print(f"  PASS: |delta| = {delta:.2e} <= {TOL_EXACT:.0e}")
                else:
                    print(f"  FAIL: |delta| = {delta:.2e} > {TOL_EXACT:.0e}")
                    all_pass = False
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Observable 2: Effective mass via parabolic fitting (Kane model)
    # Reference: 2-band Kane formula m* = Eg/(EP+Eg)
    # Uses InAs (non-W) with EP=21.5, Eg=0.417
    # ------------------------------------------------------------------
    print(f"\n{'─' * 72}")
    print("Observable 2: m*_e (Kane effective mass)")
    print(f"{'─' * 72}")

    tmpdir = tempfile.mkdtemp(prefix="star_inas_mstar_")
    try:
        config_path = os.path.join(config_dir, CONFIG_DISPERSION)
        rc, output_dir = run_executable(exe_band, config_path, tmpdir)
        if rc != 0:
            print(f"  FAIL: bandStructure returned exit code {rc}")
            all_pass = False
        else:
            eig_path = os.path.join(output_dir, "eigenvalues.dat")
            data = parse_eigenvalues(eig_path)

            if not data:
                print("  FAIL: no eigenvalue data parsed")
                all_pass = False
            else:
                k_vals, e_vals = extract_cb_dispersion(data)

                result = extract_mass_numerical_derivative(k_vals, e_vals)
                if result is None:
                    print("  FAIL: could not extract effective mass")
                    all_pass = False
                else:
                    m_star_computed, c2, e0 = result

                    print(f"  Computed m*: {m_star_computed:.4f} m0 "
                          f"(d2E/dk2 = {c2:.4f} eV*A^2)")
                    print(f"  Kane m*:    {M_STAR_KANE:.4f} m0 "
                          f"(Eg/(EP+Eg) = {INAS_EG}/({INAS_EP}+{INAS_EG}))")

                    passed, delta, _ = compare_value(
                        m_star_computed, M_STAR_KANE, TOL_KANE,
                        "m*_e", "m0"
                    )

                    bench_row = format_benchmark_row(
                        "InAs", "m*_e (Kane)", m_star_computed, M_STAR_KANE,
                        "2-band Kane formula",
                        f"Analytical ({TOL_KANE*100:.0f}%)", delta,
                        "PASS" if passed else "FAIL"
                    )
                    rows.append(bench_row)

                    if passed:
                        print(f"  PASS: |delta| = {delta:.4f} ({delta*100:.1f}%) "
                              f"<= {TOL_KANE*100:.0f}%")
                    else:
                        print(f"  FAIL: |delta| = {delta:.4f} ({delta*100:.1f}%) "
                              f"> {TOL_KANE*100:.0f}%")
                        all_pass = False
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Observable 3: g-factor (Roth vs 8-band)
    # Uses InAsW variant (EP=22.2, Eg=0.418, DeltaSO=0.38)
    # Roth formula: g = 2 - 2*EP*DeltaSO / (3*Eg*(Eg+DeltaSO))
    # Reference: Winkler 2003, Eq. (6.42)
    # ------------------------------------------------------------------
    print(f"\n{'─' * 72}")
    print("Observable 3: g* (Landau g-factor)")
    print(f"{'─' * 72}")
    print(f"  InAsW parameters: EP={INASW_EP}, Eg={INASW_EG}, "
          f"DeltaSO={INASW_DELTA_SO}")

    tmpdir = tempfile.mkdtemp(prefix="star_inas_gfactor_")
    try:
        config_path = os.path.join(config_dir, CONFIG_GFACTOR)
        rc, output_dir = run_executable(exe_gfactor, config_path, tmpdir)
        if rc != 0:
            print(f"  FAIL: gfactorCalculation returned exit code {rc}")
            all_pass = False
        else:
            gf_path = os.path.join(output_dir, "gfactor.dat")
            if not os.path.isfile(gf_path):
                print(f"  FAIL: gfactor.dat not found in {output_dir}")
                all_pass = False
            else:
                gx, gy, gz = parse_gfactor(gf_path)

                print(f"  Computed g*: gx={gx:.4f}, gy={gy:.4f}, gz={gz:.4f}")
                print(f"  Roth g*:     {G_ROTH:.4f}")
                print(f"    = 2 - 2*{INASW_EP}*{INASW_DELTA_SO} "
                      f"/ (3*{INASW_EG}*({INASW_EG}+{INASW_DELTA_SO}))")
                print(f"    = 2 - {2*INASW_EP*INASW_DELTA_SO:.4f} "
                      f"/ {3*INASW_EG*(INASW_EG+INASW_DELTA_SO):.6f}")

                # Compare gz (longitudinal g-factor) against Roth
                passed, delta, _ = compare_value(
                    gz, G_ROTH, TOL_ROTH,
                    "g*", ""
                )

                bench_row = format_benchmark_row(
                    "InAsW", "g* (Roth)", gz, G_ROTH,
                    "Winkler 2003, Eq. (6.42)",
                    f"Analytical ({TOL_ROTH*100:.0f}%)", delta,
                    "PASS" if passed else "FAIL"
                )
                rows.append(bench_row)

                if passed:
                    print(f"  PASS: |delta| = {delta:.4f} ({delta*100:.1f}%) "
                          f"<= {TOL_ROTH*100:.0f}%")
                else:
                    print(f"  FAIL: |delta| = {delta:.4f} ({delta*100:.1f}%) "
                          f"> {TOL_ROTH*100:.0f}%")
                    all_pass = False
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Summary and benchmark table
    # ------------------------------------------------------------------
    print(f"\n{'=' * 72}")
    print("Benchmark Table")
    print("=" * 72)
    print_benchmark_header()
    for row in rows:
        print(row)

    print(f"\n{'=' * 72}")
    print("Summary")
    print("=" * 72)

    if all_pass:
        print("PASS: all InAs bulk standard-star benchmarks passed")
        sys.exit(0)
    else:
        print("FAIL: one or more InAs bulk benchmarks failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
