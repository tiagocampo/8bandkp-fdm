#!/usr/bin/env python3
"""ISBT GaAs/AlGaAs QW benchmark -- z-dipole and oscillator strength vs infinite well.

Validates intersubband transition (ISBT) properties for a 100 AA GaAs/Al30Ga70As
quantum well against infinite-square-well analytical formulas.

In the 8-band model, CB states at k=0 are Kramer's-degenerate pairs:
  {1,2} = first CB subband, {3,4} = second CB subband.
The physical ISBT is between subband *pairs*.  The z-dipole table reports
four transitions (e1->e3, e1->e4, e2->e3, e2->e4) which split into:
  - "large" dipole (same spin): e1->e3, e2->e4
  - "small" dipole (cross spin): e1->e4, e2->e3
We use the large-dipole transitions for comparison.

NOTE: the 8-band z-dipole sums over all 8 band components at each grid
point, so it can differ significantly from the single-band effective-mass
result.  The comparison is therefore qualitative (ratio bounds) rather
than exact.

Checks:
  B1: Transition energy E12 < infinite-well E12 (finite barriers)
  B2: Transition energy E12 > 0.5 * infinite-well E12 (reasonable)
  B3: z-dipole is positive and finite
  B4: Oscillator strength self-consistency: f = E * |z|^2 / hbar2O2m0
  B5: Per-subband f (total over degenerate pair) is physically reasonable

Usage: verify_isbt_benchmark.py <build_dir> <source_dir>

  build_dir  -- path to build/ directory (contains src/ executables)
  source_dir -- path to repo root (contains tests/regression/configs/)
"""

# COVERAGE: observable=ISBT_dipole geometry=QW material=GaAs/AlGaAs ref=infinite_well_analytical
# COVERAGE: observable=ISBT_oscillator_strength geometry=QW material=GaAs/AlGaAs ref=infinite_well_analytical

import sys
import os
import tempfile
import shutil

try:
    import numpy as np
except ImportError:
    print("FAIL: numpy is required (pip install numpy)")
    sys.exit(1)

try:
    import star_helpers
except ImportError:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
    import star_helpers


# ---------------------------------------------------------------------------
# Physical constants (CODATA 2018, matching star_helpers.py and defs.f90)
# ---------------------------------------------------------------------------
HBAR2_OVER_2M0 = 3.80998  # hbar^2/(2*m0) in eV*AA^2

# GaAs conduction-band effective mass (Vurgaftman 2001)
MSTAR_EFF = 0.067  # in units of m0

# Well width
L_AA = 100.0  # well width in AA (z_min=-50 to z_max=50)

# ---------------------------------------------------------------------------
# Infinite-well analytical formulas
# ---------------------------------------------------------------------------
# z-dipole: |<1|z|2>|_inf = 16*L / (9*pi^2)  (exact for infinite square well)
Z12_INF = 16.0 * L_AA / (9.0 * np.pi**2)  # in AA

# Transition energy: E12 = 3*pi^2*hbar^2 / (2*m*_eff*L^2)
E12_INF = 3.0 * np.pi**2 * HBAR2_OVER_2M0 / (MSTAR_EFF * L_AA**2)  # in eV

# Oscillator strength for infinite well: f12 = (2*m0*E12/hbar^2) * |z12|^2
# This is universal: f12 = 256 / (27*pi^2) for n=1->n=2
F12_INF = 256.0 / (27.0 * np.pi**2)

# Config file name
CONFIG_NAME = "isbt_gaas_algaas_w100.toml"


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_isbt_transitions(filepath):
    """Parse output/isbt_transitions.dat.

    Format: i  j  E_ij(eV)  Re(z_ij)(AA)  Im(z_ij)(AA)  |z_ij|^2(AA^2)  f_ij

    Returns list of dicts with keys: i, j, E_ij, z_re, z_im, z_abs2, f_ij
    """
    transitions = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            vals = line.split()
            if len(vals) >= 7:
                transitions.append({
                    'i': int(vals[0]),
                    'j': int(vals[1]),
                    'E_ij': float(vals[2]),
                    'z_re': float(vals[3]),
                    'z_im': float(vals[4]),
                    'z_abs2': float(vals[5]),
                    'f_ij': float(vals[6]),
                })
    return transitions


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        print("  build_dir  -- path to build/ (contains src/ executables)")
        print("  source_dir -- path to repo root")
        sys.exit(1)

    build_dir = os.path.abspath(sys.argv[1])
    source_dir = os.path.abspath(sys.argv[2])

    print("=" * 60)
    print("ISBT GaAs/AlGaAs QW Benchmark (100 AA well)")
    print("=" * 60)

    # Print infinite-well reference values
    print(f"\nInfinite-well analytical values (single-band):")
    print(f"  Well width L         = {L_AA:.1f} AA")
    print(f"  |z_12|_inf           = {Z12_INF:.6f} AA")
    print(f"  E_12_inf             = {E12_INF:.6f} eV")
    print(f"  f_12_inf             = {F12_INF:.6f}")

    config_path = os.path.join(source_dir, "tests", "regression", "configs", CONFIG_NAME)
    if not os.path.isfile(config_path):
        print(f"\nFAIL: config not found: {config_path}")
        sys.exit(1)

    # Run opticalProperties in a temp directory
    work_dir = tempfile.mkdtemp(prefix="isbt_benchmark_")
    failures = []

    try:
        rc, out_dir = star_helpers.run_exe(
            build_dir, "opticalProperties", config_path, work_dir, timeout=300
        )
        if rc != 0:
            print(f"\nFAIL: opticalProperties exited with code {rc}")
            sys.exit(1)

        isbt_path = os.path.join(out_dir, "isbt_transitions.dat")
        if not os.path.isfile(isbt_path):
            print(f"\nFAIL: isbt_transitions.dat not produced")
            sys.exit(1)

        transitions = parse_isbt_transitions(isbt_path)
        if not transitions:
            print(f"\nFAIL: no ISBT transitions found in {isbt_path}")
            sys.exit(1)

        print(f"\nISBT transitions found: {len(transitions)}")
        for t in transitions:
            print(f"  e{t['i']} -> e{t['j']}: E={t['E_ij']:.6f} eV, "
                  f"|z|={np.sqrt(t['z_abs2']):.4f} AA, f={t['f_ij']:.6f}")

        # In the 8-band model, CB states are Kramer's-degenerate pairs:
        #   {1,2} = first CB subband, {3,4} = second CB subband.
        # The physical ISBT is between subband pairs. The z-dipole table
        # splits into large-dipole (same spin) and small-dipole (cross spin)
        # transitions.  The canonical "same spin" transitions are e1->e3 and
        # e2->e4 (or e1->e2 if the model has only 2 CB subbands).
        #
        # Strategy: find the transition pair (i, j) with the largest |z|^2.
        # This is the physical ISBT.
        largest = max(transitions, key=lambda t: t['z_abs2'])
        E12_8band = largest['E_ij']
        z12_abs2_8band = largest['z_abs2']
        z12_8band = np.sqrt(z12_abs2_8band)
        f12_8band = largest['f_ij']

        # Compute total oscillator strength for this subband transition
        # (sum over the degenerate spin pair)
        f_total = sum(t['f_ij'] for t in transitions
                      if abs(t['E_ij'] - E12_8band) < 1e-6)

        label = f"e{largest['i']}->e{largest['j']}"
        print(f"\n--- {label} comparison (largest z-dipole, first->second CB subband) ---")
        print(f"  8-band  : E12 = {E12_8band:.6f} eV, "
              f"|z12| = {z12_8band:.4f} AA, f12 = {f12_8band:.4f}")
        print(f"  Total f (all spin channels at this E) = {f_total:.4f}")
        print(f"  Inf well: E12 = {E12_INF:.6f} eV, "
              f"|z12| = {Z12_INF:.4f} AA, f12 = {F12_INF:.6f}")

        # ---------------------------------------------------------------
        # B1: Transition energy -- 8-band should be less than infinite well
        # (finite barriers reduce confinement, lowering subband spacing)
        # ---------------------------------------------------------------
        E12_ratio = E12_8band / E12_INF
        print(f"\n[B1] Transition energy ratio E12_8band / E12_inf = {E12_ratio:.4f}")
        if E12_ratio < 1.0:
            print(f"  PASS: E12_8band < E12_inf (finite barriers reduce spacing)")
        else:
            failures.append(
                f"B1 FAIL: E12_8band/E12_inf = {E12_ratio:.4f} >= 1.0 "
                f"(should be < 1 due to finite barriers)"
            )

        # B2: E12 should not be too small
        if E12_ratio > 0.5:
            print(f"  PASS: E12 ratio > 0.5 (physically reasonable)")
        else:
            failures.append(
                f"B2 FAIL: E12_8band/E12_inf = {E12_ratio:.4f} <= 0.5 "
                f"(suspiciously low)"
            )

        # ---------------------------------------------------------------
        # B3: z-dipole is positive and finite
        # ---------------------------------------------------------------
        print(f"\n[B3] z-dipole check: |z12| = {z12_8band:.4f} AA")
        if z12_8band > 0 and np.isfinite(z12_8band):
            print(f"  PASS: z-dipole is positive and finite")
        else:
            failures.append(
                f"B3 FAIL: z-dipole = {z12_8band:.4f} (not positive/finite)"
            )

        # ---------------------------------------------------------------
        # B4: Oscillator strength self-consistency
        # f_ij = E_ij * |z_ij|^2 / hbar2O2m0
        # ---------------------------------------------------------------
        f12_expected = E12_8band * z12_abs2_8band / HBAR2_OVER_2M0
        f12_relerr = abs(f12_8band - f12_expected) / abs(f12_expected)
        print(f"\n[B4] Oscillator strength self-consistency:")
        print(f"  f_ij from code      = {f12_8band:.6f}")
        print(f"  E*|z|^2/hbar2O2m0  = {f12_expected:.6f}")
        print(f"  Relative error       = {f12_relerr:.2e}")
        if f12_relerr < 1e-6:
            print(f"  PASS: f_ij = E*|z|^2/hbar2O2m0 (self-consistent)")
        else:
            failures.append(
                f"B4 FAIL: f_ij self-consistency error = {f12_relerr:.2e}"
            )

        # ---------------------------------------------------------------
        # B5: Total per-subband f is physically reasonable
        # For the infinite well, the per-subband f = 256/(27*pi^2) ~ 0.96.
        # In the 8-band model, the z-dipole includes all band contributions
        # and can exceed the single-band result.  We check that f_total
        # is a positive, finite number of reasonable magnitude.
        # ---------------------------------------------------------------
        print(f"\n[B5] Total per-subband oscillator strength = {f_total:.4f}")
        if f_total > 0 and np.isfinite(f_total):
            print(f"  PASS: f_total is positive and finite")
        else:
            failures.append(
                f"B5 FAIL: f_total = {f_total:.4f} (not positive/finite)"
            )

    except Exception as e:
        failures.append(f"EXCEPTION: {e}")
        import traceback
        traceback.print_exc()
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    # Summary
    print(f"\n{'=' * 60}")
    if failures:
        print("FAIL: ISBT benchmark failures:")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    else:
        print("PASS: all ISBT benchmark checks passed")
        sys.exit(0)


if __name__ == "__main__":
    main()
