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
TOL_STRAIN = 0.10      # 10% for wire strain (relaxation reduces agreement vs bulk)
TOL_ADDITIVE = 1e-3    # 0.1% for additive CB/HH shift (no LH-SO mixing)


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
        "nx = 40\n"
        "ny = 40\n"
        "dx = 5.0\n"
        "dy = 5.0\n"
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
        "emin = -1.5\n"
        "emax = 2.0\n"
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
                                 timeout=600)
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
    """R10a: CB eigenvalue shift matches Bir-Pikus delta_Ec.

    The CB has no LH-SO mixing, so the shift should be a clean additive
    modification: delta_Ec = ac * Tr(eps). For the wide core wire, the
    interior strain is close to biaxial, so we expect ~10% agreement.
    """
    print("  [R10a] CB eigenvalue shift vs Bir-Pikus delta_Ec")

    cb_strained = evals_strained[-1]   # highest eigenvalue (CB top)
    cb_unstrained = evals_unstrained[-1]
    cb_shift = cb_strained - cb_unstrained

    expected_shift = bp["delta_Ec"]

    passed, delta, _ = compare_value(
        cb_shift, expected_shift, TOL_STRAIN,
        "  CB shift", unit="eV",
    )
    status = "PASS" if passed else "FAIL"
    print(f"    computed shift  = {cb_shift:+.6f} eV")
    print(f"    expected delta_Ec = {expected_shift:+.6f} eV")
    print(f"    relative error  = {delta:.4e}  {status}")

    return passed


def test_r10b_hh_lh_splitting(evals_strained, bp):
    """R10b: HH-LH splitting in strained wire matches Bir-Pikus prediction.

    Under compressive strain (InAs on GaAs), HH shifts above LH. The
    wide-core wire interior approximates bulk biaxial strain, so the
    splitting should be close to the analytical value.
    """
    print("  [R10b] HH-LH splitting vs Bir-Pikus prediction")

    # FEAST returns eigenvalues sorted ascending (lowest first).
    # Taking evals_strained[:8] would pick the 8 DEEPEST VB states (near -1.5 eV),
    # not the top VB states. Instead, select the top 8 VB eigenvalues: the 8
    # largest negative eigenvalues (closest to zero from below), sorted descending.
    # Under compressive strain, HH shifts above LH, so:
    #   vb_evals[0] = highest VB state (HH top)
    #   vb_evals[1] = next VB state (LH top)
    #   vb_evals[2..3] = second LH pair, vb_evals[4..7] = SO + deeper states
    vb_evals = sorted([e for e in evals_strained if e < 0], reverse=True)[:8]

    # HH-LH splitting = energy difference between top HH and top LH
    e_hh = vb_evals[0]
    e_lh = vb_evals[1]

    splitting = e_hh - e_lh
    expected_splitting = bp["HH_LH_splitting"]

    passed, delta, _ = compare_value(
        splitting, expected_splitting, TOL_STRAIN,
        "  HH-LH splitting", unit="eV",
    )
    status = "PASS" if passed else "FAIL"
    print(f"    E_HH (VB top)   = {e_hh:+.6f} eV")
    print(f"    E_LH            = {e_lh:+.6f} eV")
    print(f"    computed split  = {splitting:+.6f} eV")
    print(f"    expected split  = {expected_splitting:+.6f} eV")
    print(f"    relative error  = {delta:.4e}  {status}")

    return passed


def test_r10c_vb_shift(evals_strained, evals_unstrained, bp):
    """R10c: VB top shift matches Bir-Pikus delta_EHH.

    The HH has no LH-SO mixing, so its shift should be exactly the
    Bir-Pikus delta_EHH = -P_eps + Q_eps.
    """
    print("  [R10c] VB top shift vs Bir-Pikus delta_EHH")

    # NOTE: Identifying VB top as max(e < 0) assumes all VB eigenvalues are
    # negative and all CB eigenvalues are positive. For InAs/GaAs with the
    # solver's energy reference (EV=0 for unstrained), the unstrained CB
    # bottom starts at Eg=0.417 eV, and the strain-induced CB shift is < 0.5 eV,
    # so the sign boundary remains valid. If extending to narrow-gap materials
    # with large strain, verify this assumption or use num_cb/num_vb from config.
    vb_top_s = max(e for e in evals_strained if e < 0)
    vb_top_u = max(e for e in evals_unstrained if e < 0)
    vb_shift = vb_top_s - vb_top_u

    expected_shift = bp["delta_EHH"]

    passed, delta, _ = compare_value(
        vb_shift, expected_shift, TOL_STRAIN,
        "  VB top shift", unit="eV",
    )
    status = "PASS" if passed else "FAIL"
    print(f"    computed shift    = {vb_shift:+.6f} eV")
    print(f"    expected delta_EHH = {expected_shift:+.6f} eV")
    print(f"    relative error    = {delta:.4e}  {status}")

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
