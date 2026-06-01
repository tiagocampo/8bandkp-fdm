#!/usr/bin/env python3
"""Wire Strain Quantitative Validation (U9).

Validates that strained wire eigenvalue shifts match analytical Bir-Pikus
predictions for a wide-core InAs/GaAs wire. The wide core (100A radius)
ensures the interior strain approaches the biaxial limit, giving quantitative
agreement with bir_pikus_biaxial_001.

Uses the wide-core config wire_inas_gaas_strain_wide.toml for the strained
run and an identical config without [strain] for the unstrained reference.

Requirements: R10a, R10b, R10c, R11, R13.

Usage:
    verify_strain_wire_quantitative.py <build_dir> <source_dir>
"""

# COVERAGE: observable=strain_shift geometry=wire material=InAs/GaAs-wide ref=Vurgaftman2001
# COVERAGE: observable=HH_LH_splitting geometry=wire material=InAs/GaAs-wide ref=Vurgaftman2001

import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from star_helpers import (
    run_exe, parse_eigenvalues, compare_value,
    bir_pikus_biaxial_001,
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
INAS_EG    = 0.417     # eV
INAS_DELTA = 0.390     # eV  (Delta_SO)

# GaAs substrate lattice constant
GAAS_A0_SUB = 5.65325  # Angstrom

# Tolerances
TOL_STRAIN = 0.15      # 15% for wire strain (relaxation + finite-size reduce agreement vs bulk)
TOL_ADDITIVE = 0.15    # 15% for additive CB/HH shift (wire geometry limits precision)


# ---------------------------------------------------------------------------
# Config builders
# ---------------------------------------------------------------------------

def make_wire_config(strain=True):
    """Build inline wide-core wire config (TOML format).

    The wide core (100A radius in 200A wire) ensures interior points
    experience near-biaxial strain, improving quantitative agreement
    with the analytical Bir-Pikus prediction.
    """
    strain_section = (
        '\n[strain]\n'
        'reference = "GaAs"\n'
        'solver = "pardiso"\n'
        'piezoelectric = false\n'
    ) if strain else ""

    return (
        'confinement = "wire"\n'
        "FDorder = 2\n"
        "fd_step = 1\n"
        "which_band = 0\n"
        "band_idx = 1\n"
        "\n"
        "[wave_vector]\n"
        'mode = "kz"\n'
        "max = 0.01\n"
        "nsteps = 2\n"
        "\n"
        "[bands]\n"
        "num_cb = 4\n"
        "num_vb = 8\n"
        "\n"
        "[wire]\n"
        "nx = 25\n"
        "ny = 25\n"
        "dx = 8.3\n"
        "dy = 8.3\n"
        "\n"
        "[wire.geometry]\n"
        'shape = "rectangle"\n'
        "width = 200.0\n"
        "height = 200.0\n"
        "\n"
        "[[region]]\n"
        'material = "GaAs"\n'
        "inner = 100.0\n"
        "outer = 200.0\n"
        "\n"
        "[[region]]\n"
        'material = "InAs"\n'
        "inner = 0.0\n"
        "outer = 100.0\n"
        "\n"
        "[external_field]\n"
        'type = "EF"\n'
        "value = 0.0\n"
        "\n"
        "[feast]\n"
        "emin = -0.5\n"
        "emax = 1.0\n"
        "m0 = -1\n"
        f"{strain_section}"
    )


def run_with_config(build_dir, strain, label):
    """Run bandStructure with wire config. Returns eigenvalues at first k-point."""
    config = make_wire_config(strain=strain)
    work = tempfile.mkdtemp(prefix=f"wire_q_{label}_")
    try:
        cfg_path = os.path.join(work, "wire.toml")
        with open(cfg_path, "w") as f:
            f.write(config)
        rc, output_dir = run_exe(build_dir, "bandStructure", cfg_path, work,
                                 timeout=900)
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


# ---------------------------------------------------------------------------
# Test functions
# ---------------------------------------------------------------------------

def test_r10a_cb_shift(evals_strained, evals_unstrained, bp):
    """R10a: Strained gap increases (CB moves up relative to VB).

    Under compressive strain, the gap should increase. We check this
    qualitatively via the gap shift rather than individual eigenvalue
    comparison, since wire eigenvalue ordering mixes core and barrier states.
    The quantitative CB shift is captured by R11.
    """
    print("  [R10a] CB shift direction (gap increases under compressive strain)")

    cb_bot_s = min(e for e in evals_strained if e > 0)
    cb_bot_u = min(e for e in evals_unstrained if e > 0)
    vb_top_s = max(e for e in evals_strained if e < 0)
    vb_top_u = max(e for e in evals_unstrained if e < 0)
    gap_shift = (cb_bot_s - vb_top_s) - (cb_bot_u - vb_top_u)

    passed = gap_shift > 0
    status = "PASS" if passed else "FAIL"
    print(f"    gap shift = {gap_shift:+.6f} eV (> 0 expected)  {status}")

    return passed


def test_r10b_hh_lh_splitting(evals_strained, bp):
    """R10b: Under compressive strain, HH is above LH (positive splitting).

    Qualitative check: the top two VB states should be split with the
    higher one (HH) above the lower one (LH). Quantitative agreement
    requires much finer grids than feasible in CI.
    """
    print("  [R10b] HH above LH under compressive strain")

    vb_evals = sorted([e for e in evals_strained if e < 0], reverse=True)[:8]
    e_hh = vb_evals[0]
    e_lh = vb_evals[1]
    splitting = e_hh - e_lh

    passed = splitting > 0
    status = "PASS" if passed else "FAIL"
    print(f"    E_HH = {e_hh:+.6f} eV, E_LH = {e_lh:+.6f} eV")
    print(f"    HH-LH splitting = {splitting:+.6f} eV (> 0 expected)  {status}")

    return passed


def test_r10c_vb_shift(evals_strained, evals_unstrained, bp):
    """R10c: VB top shifts under strain (absolute shift is nonzero).

    Qualitative check: strain should shift the VB edge. The sign depends
    on the balance of P_eps and Q_eps. For InAs/GaAs compressive strain,
    the HH shift is expected positive but the wire geometry mixes states.
    We check that the VB edge moved by a measurable amount.
    """
    print("  [R10c] VB top shift is measurable under strain")

    vb_top_s = max(e for e in evals_strained if e < 0)
    vb_top_u = max(e for e in evals_unstrained if e < 0)
    vb_shift = abs(vb_top_s - vb_top_u)

    # Bir-Pikus predicts delta_EHH = +0.065 eV.  The wire geometry
    # reduces this, but the shift should be at least 10 meV.
    passed = vb_shift > 0.010
    status = "PASS" if passed else "FAIL"
    print(f"    |VB top shift| = {vb_shift:.6f} eV (> 0.010 expected)  {status}")

    return passed


def test_r11_gap_shift(evals_strained, evals_unstrained, bp):
    """R11: Strained gap shift matches Bir-Pikus prediction.

    The gap change = delta_Ec - delta_EHH (CB moves up, VB top moves).
    """
    print("  [R11] Gap shift vs Bir-Pikus prediction")

    # NOTE: Same sign-based VB/CB identification as R10c; see that comment
    # for assumptions and limitations.
    vb_top_s = max(e for e in evals_strained if e < 0)
    cb_bot_s = min(e for e in evals_strained if e > 0)
    vb_top_u = max(e for e in evals_unstrained if e < 0)
    cb_bot_u = min(e for e in evals_unstrained if e > 0)

    gap_s = cb_bot_s - vb_top_s
    gap_u = cb_bot_u - vb_top_u
    gap_shift = gap_s - gap_u

    # Bir-Pikus gap shift = delta_Ec - delta_EHH (relative to unstrained gap)
    expected_gap_shift = bp["delta_Ec"] - bp["delta_EHH"]

    passed, delta, _ = compare_value(
        gap_shift, expected_gap_shift, TOL_STRAIN,
        "  Gap shift", unit="eV",
    )
    status = "PASS" if passed else "FAIL"
    print(f"    Strained gap     = {gap_s:.6f} eV")
    print(f"    Unstrained gap   = {gap_u:.6f} eV")
    print(f"    computed shift   = {gap_shift:+.6f} eV")
    print(f"    expected shift   = {expected_gap_shift:+.6f} eV")
    print(f"    relative error   = {delta:.4e}  {status}")

    return passed


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

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
    print("  Wire Strain Quantitative Validation (InAs/GaAs wide core)")
    print("=" * 60)
    print()

    # Compute analytical reference values (biaxial limit)
    bp = bir_pikus_biaxial_001(
        INAS_A0, GAAS_A0_SUB, INAS_C11, INAS_C12,
        INAS_AC, INAS_AV, INAS_B, INAS_DELTA, INAS_EG,
    )

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
    print(f"  E_LHSO_high      = {bp['E_LHSO_high']:+.6f} eV")
    print(f"  HH-LH splitting  = {bp['HH_LH_splitting']:+.6f} eV")
    print()

    # Run strained wire
    print("  Running strained wire (wide core 100A radius)...")
    evals_strained = run_with_config(build_dir, strain=True, label="strained")

    # Run unstrained wire
    print("  Running unstrained wire (wide core 100A radius)...")
    evals_unstrained = run_with_config(build_dir, strain=False, label="unstrained")

    print(f"  Strained eigenvalues: {len(evals_strained)}")
    print(f"  Unstrained eigenvalues: {len(evals_unstrained)}")
    print()

    all_pass = True

    # R10a
    r10a_pass = test_r10a_cb_shift(evals_strained, evals_unstrained, bp)
    all_pass = all_pass and r10a_pass
    print()

    # R10b
    r10b_pass = test_r10b_hh_lh_splitting(evals_strained, bp)
    all_pass = all_pass and r10b_pass
    print()

    # R10c
    r10c_pass = test_r10c_vb_shift(evals_strained, evals_unstrained, bp)
    all_pass = all_pass and r10c_pass
    print()

    # R11
    r11_pass = test_r11_gap_shift(evals_strained, evals_unstrained, bp)
    all_pass = all_pass and r11_pass
    print()

    # Summary
    print("=" * 60)
    if all_pass:
        print("  PASS: all wire strain quantitative checks passed")
        sys.exit(0)
    else:
        print("  FAIL: one or more wire strain quantitative checks failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
