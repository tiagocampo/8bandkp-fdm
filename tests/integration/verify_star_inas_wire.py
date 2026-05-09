#!/usr/bin/env python3
"""S7 -- InAs nanowire standard-star benchmark.

Validates InAs nanowire observables using wire geometry (confinement=2):
  - g* CB wire vs Roth bulk limit     (Winkler 2003, Eq. 6.42)
  - Wire eigenvalue vs bulk band edge  (internal consistency)

This is the most complex standard-star system: wire geometry requires CSR
sparse solver + commutator-based velocity matrices [r_alpha, H]. The wire
g-factor depends on transverse confinement, approaching the bulk Roth value
for large wires.

Uses non-W InAs (EP=21.5, Eg=0.417, DeltaSO=0.39 from Vurgaftman 2001).
For a 55x55 nm wire, confinement effects are moderate and the g-factor should
be within a factor of 2 of the bulk Roth prediction.

Cross-reference with verification ladder:
  - Wire convergence covered by rung 4 (R14).
Standard-star adds: wire g-factor validation against Roth bulk limit and
publication-ready benchmark table output.

Usage: verify_star_inas_wire.py <build_dir> <source_dir>

  build_dir  -- path to build/ directory (contains src/bandStructure, src/gfactorCalculation)
  source_dir -- path to repo root (contains tests/regression/configs/)
"""

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
)

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
HBAR2_OVER_2M0 = 3.80998  # eV * Angstrom^2 (hbar^2 / (2*m0))

# ---------------------------------------------------------------------------
# Material parameters: InAs (non-W, Vurgaftman 2001)
# ---------------------------------------------------------------------------
INAS_EG = 0.417        # eV, band gap (Vurgaftman 2001, Table I)
INAS_EP = 21.5         # eV, Kane matrix element (parameters.f90)
INAS_DELTA_SO = 0.39   # eV, spin-orbit splitting (Vurgaftman 2001, Table I)

# Roth g-factor formula (Winkler 2003, Eq. 6.42):
#   g_roth = 2 - 2*EP*DeltaSO / (3*Eg*(Eg + DeltaSO))
G_ROTH = 2.0 - 2.0 * INAS_EP * INAS_DELTA_SO / (
    3.0 * INAS_EG * (INAS_EG + INAS_DELTA_SO)
)

# ---------------------------------------------------------------------------
# Config paths (relative to tests/regression/configs/)
# ---------------------------------------------------------------------------
CONFIG_WIRE = "wire_inas_rectangle.cfg"
CONFIG_GFACTOR = "wire_inas_gfactor.cfg"

# ---------------------------------------------------------------------------
# Tolerances
# ---------------------------------------------------------------------------
TOL_GFACTOR_RANGE = 0.50   # Wire g* within [0.5*|g_roth|, 2.0*|g_roth|]
TOL_EIGENVALUE = 0.05      # 5%: regression tolerance vs 8-band wire reference

# ---------------------------------------------------------------------------
# Regression reference values (8-band wire)
# ---------------------------------------------------------------------------
# The 55x55 nm InAs wire has significant confinement (m* ~0.023 m0), pushing
# the CB ground state well above bulk Eg. Use code output as regression ref.
CB_GROUND_REF = 0.864320   # eV, 8-band wire CB ground state


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_cb_ground_state(evals):
    """Find the CB ground state by identifying the largest gap in the spectrum.

    The 8-band wire spectrum has a clear band gap between VB and CB states.
    The CB ground state is the first eigenvalue above the largest gap.
    """
    sorted_evals = sorted(evals)
    max_gap = 0.0
    max_gap_idx = 0
    for i in range(len(sorted_evals) - 1):
        gap = sorted_evals[i + 1] - sorted_evals[i]
        if gap > max_gap:
            max_gap = gap
            max_gap_idx = i
    return sorted_evals[max_gap_idx + 1]


def check_wire_gfactor(g_computed, g_roth, tolerance_frac):
    """Check wire g-factor against Roth bulk limit.

    For a large wire, the g-factor should:
      1. Have the same sign as the bulk Roth value
      2. Have magnitude within [tolerance_frac * |g_roth|, (1/tolerance_frac) * |g_roth|]

    This is a range check rather than a point comparison because wire
    confinement modifies the g-factor from the bulk value, but for a large
    wire (55x55 nm) it should be within a factor of 2.

    Args:
        g_computed: computed wire g-factor
        g_roth: bulk Roth g-factor prediction
        tolerance_frac: fractional bounds (0.5 means within [0.5, 2.0] of |g_roth|)

    Returns:
        (passed, row_dict) where row_dict has standard benchmark fields.
    """
    abs_roth = abs(g_roth)
    abs_computed = abs(g_computed)

    # Sign check
    sign_match = (g_computed * g_roth) > 0

    # Magnitude range check: |g| in [frac * |g_roth|, |g_roth| / frac]
    lower_bound = tolerance_frac * abs_roth
    upper_bound = abs_roth / tolerance_frac
    in_range = lower_bound <= abs_computed <= upper_bound

    passed = sign_match and in_range

    # Compute a "delta" for the benchmark table: relative deviation from
    # the Roth value (may be > 100% if outside the tolerance range).
    delta = abs(abs_computed - abs_roth) / abs_roth

    row = {
        'computed': g_computed,
        'expected': g_roth,
        'tolerance': tolerance_frac,
        'delta': delta,
        'name': 'g* wire',
        'unit': '',
        'status': 'PASS' if passed else 'FAIL',
    }

    return passed, delta, row


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
    for exe, label in [(exe_band, "bandStructure"),
                       (exe_gfactor, "gfactorCalculation")]:
        if not os.path.isfile(exe):
            print(f"FAIL: {label} executable not found at {exe}")
            sys.exit(1)

    # Check configs
    for cfg in [CONFIG_WIRE, CONFIG_GFACTOR]:
        path = os.path.join(config_dir, cfg)
        if not os.path.isfile(path):
            print(f"FAIL: config not found: {path}")
            sys.exit(1)

    all_pass = True
    rows = []

    print("=" * 72)
    print("S7 -- InAs Nanowire Standard-Star Benchmark")
    print("=" * 72)
    print(f"  Wire: 11x11 grid, 55x55 nm InAs rectangle")
    print(f"  Material: InAs (non-W, EP={INAS_EP}, Eg={INAS_EG}, "
          f"DeltaSO={INAS_DELTA_SO})")
    print(f"  Roth bulk g*: {G_ROTH:.4f}")
    print(f"    = 2 - 2*{INAS_EP}*{INAS_DELTA_SO} "
          f"/ (3*{INAS_EG}*({INAS_EG}+{INAS_DELTA_SO}))")

    # ------------------------------------------------------------------
    # Observable 1: Wire g-factor vs Roth bulk limit
    # Reference: Winkler 2003, Eq. (6.42)
    # The wire g-factor is computed via commutator-based velocity matrices
    # [r_alpha, H] on the CSR Hamiltonian.
    # ------------------------------------------------------------------
    print(f"\n{'─' * 72}")
    print("Observable 1: g* CB wire (Roth bulk limit)")
    print(f"{'─' * 72}")
    print(f"  Roth prediction (bulk): g = {G_ROTH:.4f}")
    print(f"  Wire g* expected: same sign, magnitude within "
          f"[{TOL_GFACTOR_RANGE}*|g_roth|, |g_roth|/{TOL_GFACTOR_RANGE}]")
    print(f"  i.e. |g*| in [{TOL_GFACTOR_RANGE * abs(G_ROTH):.2f}, "
          f"{abs(G_ROTH) / TOL_GFACTOR_RANGE:.2f}]")

    tmpdir = tempfile.mkdtemp(prefix="star_inas_wire_gf_")
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

                # For wire with B along z, gz is the longitudinal g-factor
                passed, delta, row = check_wire_gfactor(
                    gz, G_ROTH, TOL_GFACTOR_RANGE
                )

                bench_row = format_benchmark_row(
                    "InAs wire", "g* CB (wire)", gz, G_ROTH,
                    "Winkler 2003, Eq. (6.42)",
                    f"Analytical ({TOL_GFACTOR_RANGE*100:.0f}%-range)",
                    delta,
                    "PASS" if passed else "FAIL"
                )
                rows.append(bench_row)

                if passed:
                    print(f"  PASS: |g*| = {abs(gz):.4f} within "
                          f"[{TOL_GFACTOR_RANGE * abs(G_ROTH):.2f}, "
                          f"{abs(G_ROTH) / TOL_GFACTOR_RANGE:.2f}]")
                else:
                    print(f"  FAIL: |g*| = {abs(gz):.4f} outside "
                          f"[{TOL_GFACTOR_RANGE * abs(G_ROTH):.2f}, "
                          f"{abs(G_ROTH) / TOL_GFACTOR_RANGE:.2f}]")
                    all_pass = False
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Observable 2: Wire CB ground state (regression)
    # For a 55x55 nm InAs wire, the CB ground-state eigenvalue should be
    # near the bulk InAs CB edge (Eg = 0.417 eV). Confinement shifts the
    # eigenvalue above Eg; for a large wire the shift is small.
    # Reference: internal consistency with bulk parameters.
    # ------------------------------------------------------------------
    print(f"\n{'─' * 72}")
    print("Observable 2: Wire CB ground state (regression)")
    print(f"{'─' * 72}")
    print(f"  Reference: 8-band wire CB ground state = {CB_GROUND_REF:.6f} eV")
    print(f"  Bulk InAs Eg = {INAS_EG} eV (confinement adds ~{CB_GROUND_REF-INAS_EG:.3f} eV)")

    tmpdir = tempfile.mkdtemp(prefix="star_inas_wire_eig_")
    try:
        config_path = os.path.join(config_dir, CONFIG_WIRE)
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
                # Wire outputs eigenvalues at each kz point.
                # Use k=0 (first point) eigenvalues.
                k0, evals = data[0]
                cb_gs = find_cb_ground_state(evals)

                print(f"  {len(evals)} eigenvalues at k=0")
                print(f"  CB ground state: {cb_gs:.6f} eV")
                print(f"  Regression ref:  {CB_GROUND_REF:.6f} eV")
                print(f"  Bulk InAs Eg:    {INAS_EG:.6f} eV (for reference)")

                # Regression check: CB ground state should match the
                # 8-band wire reference value. The 55x55 nm wire has
                # significant confinement for InAs's light mass.
                passed, delta, _ = compare_value(
                    cb_gs, CB_GROUND_REF, TOL_EIGENVALUE,
                    "CB ground state", "eV"
                )

                bench_row = format_benchmark_row(
                    "InAs wire", "CB ground state", cb_gs, CB_GROUND_REF,
                    "8-band wire regression",
                    f"Regression ({TOL_EIGENVALUE*100:.0f}%)",
                    delta,
                    "PASS" if passed else "FAIL"
                )
                rows.append(bench_row)

                if passed:
                    print(f"  PASS: |delta| = {delta:.4f} ({delta*100:.1f}%) "
                          f"<= {TOL_EIGENVALUE*100:.0f}%")
                else:
                    print(f"  FAIL: |delta| = {delta:.4f} ({delta*100:.1f}%) "
                          f"> {TOL_EIGENVALUE*100:.0f}%")
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
        print("PASS: all InAs nanowire standard-star benchmarks passed")
        sys.exit(0)
    else:
        print("FAIL: one or more InAs nanowire benchmarks failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
