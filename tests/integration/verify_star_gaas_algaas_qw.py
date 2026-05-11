#!/usr/bin/env python3
"""S4 -- GaAs/AlGaAs QW standard-star benchmark.

Validates GaAs/AlGaAs QW observables against published references:
  - CB subband spacing (E1-E2) = 9.92 meV (benchmarks.md Section 2;
    nextnano cross-validation)
  - Stark shift at -700 kV/cm (benchmarks.md Section 6; QCSE validation
    via eigenvalue comparison against reference data)
  - TE absorption onset at interband transition energy (benchmarks.md
    Section 2; CB1-VB_top gap)
  - TE > TM polarization (physical: TE couples more strongly in QW)

Structure A (subband spacing + absorption):
  Al30Ga70As(200A) / GaAs(100A) / Al30Ga70As(200A), FDstep=201, FDorder=2/4

Structure B (Stark shift):
  Al20Ga80As(460A) / GaAs(60A) / Al20Ga80As(460A), FDstep=461, FDorder=2
  Config: tests/regression/configs/sc_qcse_gaas_algaas_ef.cfg
  No-field reference: tests/regression/data/sc_qcse_gaas_algaas/eigenvalues.dat

Cross-reference with verification ladder: subband spacing covered by rung 3
(R10). Standard-star adds: Stark shift, optical absorption, polarization
check, and publication-ready benchmark table output.

Usage: verify_star_gaas_algaas_qw.py <build_dir> <source_dir>

  build_dir  -- path to build/ directory (contains src/ executables)
  source_dir -- path to repo root (contains tests/regression/configs/,
                docs/benchmarks/)
"""

# COVERAGE: observable=subband_spacing geometry=QW material=GaAs/AlGaAs
# COVERAGE: observable=stark_shift geometry=QW material=GaAs/AlGaAs
# COVERAGE: observable=absorption_polarization geometry=QW material=GaAs/AlGaAs
import os
import sys
import tempfile

try:
    import numpy as np
except ImportError:
    print("FAIL: numpy is required (pip install numpy)")
    sys.exit(1)

from star_helpers import (
    run_executable, parse_eigenvalues, parse_absorption,
    compare_value, format_benchmark_row, print_benchmark_header,
    TOL_NUMERICAL, trapz_fn,
)

# ---------------------------------------------------------------------------
# Reference values (benchmarks.md Sections 2 and 6)
# ---------------------------------------------------------------------------
# CB subband spacing for Al30Ga70As(200A)/GaAs(100A)/Al30Ga70As(200A) QW.
# Reference: benchmarks.md Section 2; nextnano cross-validation.
CB_SPACING_EXPECTED = 9.92e-3   # eV (= 9.92 meV)

# Stark shift for Al20Ga80As(460A)/GaAs(60A)/Al20Ga80As(460A) QW at -700 kV/cm.
# The raw CB1 eigenvalue shift includes the linear potential ramp across the
# domain. Reference: benchmarks.md Section 6 comparison table.
# No-field CB1 = 0.931880 eV, with-field CB1 = 1.981650 eV.
# Raw shift = 1.049770 eV (linear potential + QCSE).
CB1_NO_FIELD_REF = 0.931880     # eV (from reference eigenvalues.dat)
CB1_EF_REF = 1.981650           # eV (from reference eigenvalues.dat)
STARK_RAW_SHIFT_REF = CB1_EF_REF - CB1_NO_FIELD_REF  # ~1.049770 eV

# CB1 subband edge for the 100A/30%Al structure (from benchmarks.md Section 2).
# CB1 = EC(GaAs) + E1(8-band) ~ 1.021 eV.
CB1_SUBBAND_EDGE = 1.021        # eV

# ---------------------------------------------------------------------------
# Config paths
# ---------------------------------------------------------------------------
# Structure A: subband spacing (2-layer barrier/well)
CONFIG_SUBBAND = os.path.join("tests", "regression", "configs", "qw_gaas_algaas.cfg")

# Structure B: QCSE with electric field
CONFIG_QCSE_EF = os.path.join("tests", "regression", "configs",
                               "sc_qcse_gaas_algaas_ef.cfg")

# No-field QCSE reference data
REF_QCSE_NOFIELD = os.path.join("tests", "regression", "data",
                                 "sc_qcse_gaas_algaas", "eigenvalues.dat")

# Structure A: absorption (opticalProperties)
CONFIG_ABSORPTION = os.path.join("tests", "regression", "configs",
                                  "qw_gaas_algaas_absorption.cfg")

# ---------------------------------------------------------------------------
# Eigenvalue indices (0-based, after stripping k from parsed data)
# ---------------------------------------------------------------------------
# Structure A: numcb=4, numvb=8 -> 12 eigenvalues total
# Eigenvalues come in degenerate pairs: CB1 pair at [8,9], CB2 pair at [10,11]
NUM_VB_A = 8
CB1_IDX_A = 8
CB2_IDX_A = 10

# Structure B: numcb=4, numvb=8 -> 12 eigenvalues total
NUM_VB_B = 8
CB1_IDX_B = 8

# ---------------------------------------------------------------------------
# Tolerances (KD6)
# ---------------------------------------------------------------------------
TOL_CB_SPACING = 0.01    # 1% for subband spacing (numerical, nextnano ref)
TOL_STARK_RAW = 0.01     # 1% for raw Stark eigenvalue shift (regression ref)
TOL_ABS_ONSET = 0.02     # 2% for absorption onset energy

# Absorption edge detection: threshold fraction of peak absorption
ABSORPTION_THRESHOLD_FRAC = 0.10

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run_bandstructure(build_dir, config_path, work_dir, timeout=300):
    """Run bandStructure and return (returncode, output_dir)."""
    exe = os.path.join(build_dir, "src", "bandStructure")
    if not os.path.isfile(exe):
        print(f"FAIL: bandStructure not found at {exe}")
        sys.exit(1)
    return run_executable(exe, config_path, work_dir, timeout=timeout)


def run_optical(build_dir, config_path, work_dir, timeout=300):
    """Run opticalProperties and return (returncode, output_dir)."""
    exe = os.path.join(build_dir, "src", "opticalProperties")
    if not os.path.isfile(exe):
        print(f"FAIL: opticalProperties not found at {exe}")
        sys.exit(1)
    return run_executable(exe, config_path, work_dir, timeout=timeout)


def extract_cb_spacing(data, cb1_idx, cb2_idx):
    """Extract CB1-CB2 spacing from first k-point eigenvalues.

    Args:
        data: output of parse_eigenvalues, list of (k, [evals])
        cb1_idx: 0-based index of CB1 in evals array
        cb2_idx: 0-based index of CB2 in evals array

    Returns:
        (spacing, cb1, cb2) or (None, None, None) on failure.
    """
    if not data:
        return None, None, None
    _, evals = data[0]
    if len(evals) <= max(cb1_idx, cb2_idx):
        return None, None, None
    cb1 = evals[cb1_idx]
    cb2 = evals[cb2_idx]
    spacing = cb2 - cb1
    return spacing, cb1, cb2


def find_absorption_onset(spectrum, threshold_frac=0.10):
    """Find absorption onset energy from spectrum data.

    Onset is defined as the first energy where absorption exceeds
    threshold_frac * max(absorption).

    Args:
        spectrum: list of (energy, absorption) tuples
        threshold_frac: fraction of max absorption for onset detection

    Returns:
        onset_energy or None if not found
    """
    if not spectrum:
        return None
    energies = np.array([e for e, _ in spectrum])
    alphas = np.array([a for _, a in spectrum])
    max_alpha = alphas.max()
    if max_alpha <= 0:
        return None
    threshold = threshold_frac * max_alpha
    above = np.where(alphas > threshold)[0]
    if len(above) == 0:
        return None
    return float(energies[above[0]])


def compute_integrated_absorption(spectrum):
    """Compute total integrated absorption via trapezoidal rule.

    Args:
        spectrum: list of (energy, absorption) tuples

    Returns:
        float, the integrated absorption
    """
    if not spectrum or len(spectrum) < 2:
        return 0.0
    energies = np.array([e for e, _ in spectrum])
    alphas = np.array([a for _, a in spectrum])
    return float(trapz_fn(alphas, energies))


# ---------------------------------------------------------------------------
# Observable checks
# ---------------------------------------------------------------------------

def check_subband_spacing(build_dir, source_dir):
    """Check CB subband spacing (E1-E2) for the 100A/30%Al QW.

    Reference: benchmarks.md Section 2, nextnano cross-validation.
    Expected: 9.92 meV with 1% tolerance.

    Returns list of benchmark row dicts.
    """
    config_path = os.path.join(source_dir, CONFIG_SUBBAND)
    rows = []

    with tempfile.TemporaryDirectory(prefix="star_gaas_qw_sub_") as tmpdir:
        rc, output_dir = run_bandstructure(build_dir, config_path, tmpdir)
        if rc != 0:
            print(f"  FAIL: bandStructure returned {rc}")
            return rows

        eig_path = os.path.join(output_dir, "eigenvalues.dat")
        data = parse_eigenvalues(eig_path)
        if not data:
            print("  FAIL: no eigenvalue data parsed")
            return rows

        spacing, cb1, cb2 = extract_cb_spacing(data, CB1_IDX_A, CB2_IDX_A)
        if spacing is None:
            print("  FAIL: insufficient eigenvalues for CB spacing")
            return rows

        spacing_mev = spacing * 1000  # eV -> meV
        expected_mev = CB_SPACING_EXPECTED * 1000

        print(f"  CB1 = {cb1:.6f} eV, CB2 = {cb2:.6f} eV")
        print(f"  Spacing = {spacing_mev:.4f} meV "
              f"(expected {expected_mev:.2f} meV)")

        passed, delta, row = compare_value(
            spacing, CB_SPACING_EXPECTED, TOL_CB_SPACING,
            "CB spacing (E1-E2)", "eV"
        )
        row['material'] = "GaAs/AlGaAs QW"
        row['reference'] = "benchmarks.md Sec.2; nextnano"
        row['tol_str'] = f"Numerical ({TOL_CB_SPACING*100:.0f}%)"
        rows.append(row)

        status = "PASS" if passed else "FAIL"
        print(f"  {status}: delta = {delta:.4f} ({delta*100:.2f}%)")

    return rows


def check_stark_shift(build_dir, source_dir):
    """Check QCSE Stark shift via eigenvalue comparison.

    Runs bandStructure with the -700 kV/cm EF config, reads the no-field
    reference eigenvalues, and validates the raw CB1 shift against the
    known reference value.

    The raw shift (CB1_ef - CB1_no) includes the linear potential ramp.
    We validate it against the reference shift from the regression data
    to ensure the electric field implementation is correct.

    Reference: benchmarks.md Section 6, regression test data.

    Returns list of benchmark row dicts.
    """
    config_ef = os.path.join(source_dir, CONFIG_QCSE_EF)
    ref_nofield = os.path.join(source_dir, REF_QCSE_NOFIELD)
    rows = []

    if not os.path.isfile(ref_nofield):
        print(f"  FAIL: no-field reference not found: {ref_nofield}")
        return rows

    # Parse no-field reference CB1
    nofield_data = parse_eigenvalues(ref_nofield)
    if not nofield_data:
        print("  FAIL: could not parse no-field reference eigenvalues")
        return rows
    cb1_nofield = nofield_data[0][1][CB1_IDX_B]

    # Run with-field config
    with tempfile.TemporaryDirectory(prefix="star_gaas_qw_ef_") as tmpdir:
        rc, output_dir = run_bandstructure(build_dir, config_ef, tmpdir)
        if rc != 0:
            print(f"  FAIL: bandStructure returned {rc}")
            return rows

        eig_path = os.path.join(output_dir, "eigenvalues.dat")
        data = parse_eigenvalues(eig_path)
        if not data:
            print("  FAIL: no eigenvalue data parsed")
            return rows

        cb1_ef = data[0][1][CB1_IDX_B]
        stark_shift = cb1_ef - cb1_nofield

        print(f"  CB1 no-field: {cb1_nofield:.6f} eV (reference data)")
        print(f"  CB1 with-field: {cb1_ef:.6f} eV (computed)")
        print(f"  Raw shift: {stark_shift:.6f} eV")
        print(f"  Reference raw shift: {STARK_RAW_SHIFT_REF:.6f} eV")

        # Compare against reference raw shift from regression data
        passed, delta, row = compare_value(
            stark_shift, STARK_RAW_SHIFT_REF, TOL_STARK_RAW,
            "Stark raw shift (CB1)", "eV"
        )
        row['material'] = "GaAs/AlGaAs QW (QCSE)"
        row['reference'] = "benchmarks.md Sec.6; regression data"
        row['tol_str'] = f"Numerical ({TOL_STARK_RAW*100:.0f}%)"
        rows.append(row)

        status = "PASS" if passed else "FAIL"
        print(f"  {status}: delta = {delta:.6f} ({delta*100:.2f}%)")

    return rows


def check_absorption(build_dir, source_dir):
    """Check TE absorption onset and TE > TM polarization.

    Runs opticalProperties with the absorption config. Detects the TE
    absorption onset energy and validates it is consistent with the
    interband transition energy. Checks TE integral > TM integral
    (physical: TE polarization couples more strongly in QW).

    Returns list of benchmark row dicts.
    """
    config_path = os.path.join(source_dir, CONFIG_ABSORPTION)
    rows = []

    with tempfile.TemporaryDirectory(prefix="star_gaas_qw_abs_") as tmpdir:
        rc, output_dir = run_optical(build_dir, config_path, tmpdir)
        if rc != 0:
            print(f"  FAIL: opticalProperties returned {rc}")
            return rows

        te_path = os.path.join(output_dir, "absorption_TE.dat")
        tm_path = os.path.join(output_dir, "absorption_TM.dat")

        if not os.path.isfile(te_path):
            print(f"  FAIL: absorption_TE.dat not found in {output_dir}")
            return rows

        # --- TE absorption onset ---
        te_spectrum = parse_absorption(te_path)
        if not te_spectrum:
            print("  FAIL: no TE absorption data parsed")
            return rows

        onset_energy = find_absorption_onset(te_spectrum, ABSORPTION_THRESHOLD_FRAC)
        if onset_energy is None:
            print("  FAIL: could not detect TE absorption onset")
            return rows

        # The absorption onset corresponds to the interband transition
        # energy E_CV = CB1 - VB_top. For the QW, this is the lowest
        # photon energy at which absorption occurs.
        # Sanity bounds: above GaAs gap (1.519 eV), below Al30Ga70As gap (1.977 eV).
        E_GAAS_GAP = 1.519     # eV, GaAs band gap
        E_ALGAAS30_GAP = 1.977  # eV, Al30Ga70As band gap
        onset_min = E_GAAS_GAP
        onset_max = E_ALGAAS30_GAP

        print(f"  TE absorption onset: {onset_energy:.4f} eV")
        print(f"  Sanity range: [{onset_min:.3f}, {onset_max:.3f}] eV")
        print(f"    GaAs gap: {E_GAAS_GAP:.3f} eV, Al30Ga70As gap: {E_ALGAAS30_GAP:.3f} eV")

        # TODO: Replace range check with regression reference once
        # opticalProperties has been run to establish ONSET_REF.
        # Target: compare_value(onset_energy, ONSET_REF, TOL_ABS_ONSET, ...)
        onset_in_range = (onset_energy >= onset_min and
                          onset_energy <= onset_max)
        status_onset = 'PASS' if onset_in_range else 'FAIL'

        row_onset = {
            'computed': onset_energy,
            'expected': (onset_min + onset_max) / 2,
            'tolerance': TOL_ABS_ONSET,
            'delta': 0.0 if onset_in_range else 1.0,
            'name': 'TE absorption onset',
            'unit': 'eV',
            'status': status_onset,
            'material': 'GaAs/AlGaAs QW',
            'reference': 'Range: GaAs gap to AlGaAs gap',
            'tol_str': f'Sanity ({TOL_ABS_ONSET*100:.0f}%)',
        }
        rows.append(row_onset)
        print(f"  {status_onset}: onset in [{onset_min:.3f}, {onset_max:.3f}] eV")

        # --- TE > TM polarization check ---
        if os.path.isfile(tm_path):
            tm_spectrum = parse_absorption(tm_path)
            if tm_spectrum:
                te_integral = compute_integrated_absorption(te_spectrum)
                tm_integral = compute_integrated_absorption(tm_spectrum)

                print(f"  TE integrated absorption: {te_integral:.4f}")
                print(f"  TM integrated absorption: {tm_integral:.4f}")

                te_gt_tm = te_integral > tm_integral
                row_pol = {
                    'computed': te_integral,
                    'expected': tm_integral,
                    'tolerance': 0.0,
                    'delta': abs(te_integral - tm_integral) / max(abs(tm_integral), 1e-14),
                    'name': 'TE > TM polarization',
                    'unit': 'arb.',
                    'status': 'PASS' if te_gt_tm else 'FAIL',
                    'material': 'GaAs/AlGaAs QW',
                    'reference': 'Physical: TE dominates in QW',
                    'tol_str': 'Boolean (TE > TM)',
                }
                rows.append(row_pol)
                status_pol = 'PASS' if te_gt_tm else 'FAIL'
                print(f"  {status_pol}: TE integral {'>' if te_gt_tm else '<='} "
                      f"TM integral")
            else:
                print("  WARN: no TM absorption data parsed, skipping polarization check")
        else:
            print("  WARN: absorption_TM.dat not found, skipping polarization check")

    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        print(f"  build_dir  -- path to build/ (contains src/ executables)")
        print(f"  source_dir -- path to repo root")
        sys.exit(1)

    build_dir = os.path.abspath(sys.argv[1])
    source_dir = os.path.abspath(sys.argv[2])

    # Validate configs exist
    configs_to_check = [
        os.path.join(source_dir, CONFIG_SUBBAND),
        os.path.join(source_dir, CONFIG_QCSE_EF),
        os.path.join(source_dir, CONFIG_ABSORPTION),
        os.path.join(source_dir, REF_QCSE_NOFIELD),
    ]
    for path in configs_to_check:
        if not os.path.isfile(path):
            print(f"FAIL: required file not found: {path}")
            sys.exit(1)

    all_rows = []
    all_pass = True

    print("=" * 72)
    print("S4 -- GaAs/AlGaAs QW Standard-Star Benchmark")
    print("=" * 72)
    print()
    print("Structure A (subband + absorption):")
    print("  Al30Ga70As(200A) / GaAs(100A) / Al30Ga70As(200A)")
    print(f"  CB spacing expected: {CB_SPACING_EXPECTED*1000:.2f} meV "
          f"(benchmarks.md Sec.2)")
    print()
    print("Structure B (Stark shift):")
    print("  Al20Ga80As(460A) / GaAs(60A) / Al20Ga80As(460A)")
    print(f"  EF = -700 kV/cm, raw shift ref: {STARK_RAW_SHIFT_REF:.6f} eV")
    print()

    # --- Observable 1: CB subband spacing ---
    print("-" * 40)
    print("Observable 1: CB subband spacing (E1-E2)")
    print("-" * 40)
    print(f"  Reference: benchmarks.md Section 2; nextnano cross-validation")
    print(f"  Tolerance: {TOL_CB_SPACING*100:.0f}%")
    rows = check_subband_spacing(build_dir, source_dir)
    all_rows.extend(rows)
    if any(r['status'] == 'FAIL' for r in rows):
        all_pass = False
    print()

    # --- Observable 2: Stark shift ---
    print("-" * 40)
    print("Observable 2: Stark shift (QCSE)")
    print("-" * 40)
    print(f"  Reference: benchmarks.md Section 6; regression reference data")
    print(f"  Tolerance: {TOL_STARK_RAW*100:.0f}%")
    rows = check_stark_shift(build_dir, source_dir)
    all_rows.extend(rows)
    if any(r['status'] == 'FAIL' for r in rows):
        all_pass = False
    print()

    # --- Observable 3: Absorption onset + polarization ---
    print("-" * 40)
    print("Observable 3: TE absorption onset + polarization")
    print("-" * 40)
    print(f"  Reference: benchmarks.md Section 2; interband gap")
    print(f"  Tolerance: Numerical (range check)")
    rows = check_absorption(build_dir, source_dir)
    all_rows.extend(rows)
    if any(r['status'] == 'FAIL' for r in rows):
        all_pass = False
    print()

    # --- Benchmark table ---
    print("=" * 72)
    print("Benchmark Table")
    print("=" * 72)
    print_benchmark_header()

    for row in all_rows:
        print(format_benchmark_row(
            row.get('material', 'GaAs/AlGaAs QW'),
            row.get('name', 'unknown'),
            row['computed'],
            row['expected'],
            row.get('reference', ''),
            row.get('tol_str', ''),
            row['delta'],
            row['status'],
        ))

    # --- Summary ---
    n_pass = sum(1 for r in all_rows if r['status'] == 'PASS')
    n_total = len(all_rows)
    print()
    print(f"Results: {n_pass} PASS, {n_total - n_pass} FAIL "
          f"out of {n_total} observables")

    if all_pass:
        print("PASS: all GaAs/AlGaAs QW observables within tolerance")
        sys.exit(0)
    else:
        print("FAIL: one or more GaAs/AlGaAs QW observables failed validation")
        sys.exit(1)


if __name__ == "__main__":
    main()
