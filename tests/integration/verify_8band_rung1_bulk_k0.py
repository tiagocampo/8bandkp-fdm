#!/usr/bin/env python3
"""Rung 1 — Bulk k=0 Structural Validation (U1).

Validates bulk eigenvalues, degeneracies, basis ordering, T_d symmetry,
and eigenfunction normalization at k=0 for 5 materials:
GaAs, InAs, InSb, GaSb, GaSbW.

Requirements: R1, R2, R3, R4, R5.

Usage:
    verify_8band_rung1_bulk_k0.py <build_dir> <test_dir>

    build_dir  — path to the build directory (contains src/bandStructure)
    test_dir   — path to the repo root (contains tests/regression/configs/)

The script discovers configs named bulk_<material>_k0.toml in
<test_dir>/tests/regression/configs/ and runs bandStructure for each one
in a temporary directory.
"""

# COVERAGE: observable=Eg geometry=bulk material=GaAs ref=Vurgaftman2001
# COVERAGE: observable=Delta_SO geometry=bulk material=GaAs ref=Vurgaftman2001
# COVERAGE: observable=Eg geometry=bulk material=InAs ref=Vurgaftman2001
# COVERAGE: observable=Delta_SO geometry=bulk material=InAs ref=Vurgaftman2001
# COVERAGE: observable=Eg geometry=bulk material=InSb ref=Vurgaftman2001
# COVERAGE: observable=Delta_SO geometry=bulk material=InSb ref=Vurgaftman2001
# COVERAGE: observable=Eg geometry=bulk material=GaSb ref=Vurgaftman2001
# COVERAGE: observable=Delta_SO geometry=bulk material=GaSb ref=Vurgaftman2001
# COVERAGE: observable=Eg geometry=bulk material=GaSbW ref=Winkler2003
# COVERAGE: observable=Delta_SO geometry=bulk material=GaSbW ref=Winkler2003
import os
import sys
import tempfile
import shutil
import subprocess

# ---------------------------------------------------------------------------
# Material database — values from parameters.f90
# ---------------------------------------------------------------------------
# Bulk mode (ZB8bandBulk) places Eg and DeltaSO directly on the Hamiltonian
# diagonal.  At k=0 the eigenvalues are exactly:
#   [-DeltaSO, -DeltaSO, 0, 0, 0, 0, Eg, Eg]  (ascending order)
# Bulk mode does NOT use EV/EC, so GaSb non-W (which lacks EV/EC) works.
# ---------------------------------------------------------------------------
MATERIALS = {
    "GaAs": {
        "Eg": 1.519,
        "DeltaSO": 0.341,
    },
    "InAs": {
        "Eg": 0.417,
        "DeltaSO": 0.390,
    },
    "InSb": {
        "Eg": 0.235,
        "DeltaSO": 0.810,
    },
    "GaSb": {
        "Eg": 0.812,
        "DeltaSO": 0.760,
    },
    "GaSbW": {
        "Eg": 0.812,
        "DeltaSO": 0.760,
    },
}

# Tolerances
TOL_MACHINE = 1e-12   # R1: machine-precision eigenvalue comparison
TOL_DEGEN = 1e-14     # R2: degeneracy tolerance
TOL_NORM = 1e-10      # R5: eigenfunction normalization
TOL_CHARACTER = 0.90   # R3: minimum fractional band character


def parse_eigenvalues(filepath):
    """Parse eigenvalues file, returning list of eigenvalues for the first k-point.

    Returns a list of floats (eigenvalues in ascending order).
    """
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            vals = [float(x) for x in line.split()]
            # First value is |k|, rest are eigenvalues
            return vals[1:]
    return None


def parse_parts(filepath):
    """Parse parts.dat file, returning list of 8-element lists.

    Each row is one eigenstate, columns 1-8 are per-band weights.
    For bulk mode, raw weights are |psi_band|^2 and should sum to 1.
    """
    rows = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            vals = [float(x) for x in line.split()]
            if len(vals) >= 8:
                rows.append(vals[:8])
    return rows


def run_bandstructure(exe_path, config_path, work_dir):
    """Run bandStructure executable in work_dir with the given config.

    Copies config to <work_dir>/input.toml, runs the executable,
    returns (returncode, eigenvalues_path, parts_path).
    """
    # Copy config to input.toml in the work directory
    dst_cfg = os.path.join(work_dir, "input.toml")
    shutil.copy2(config_path, dst_cfg)

    # Ensure output directory exists
    output_dir = os.path.join(work_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    # Run the executable
    result = subprocess.run(
        [exe_path],
        cwd=work_dir,
        capture_output=True,
        text=True,
        timeout=120,
    )

    eigenvalues_path = os.path.join(output_dir, "eigenvalues.dat")
    parts_path = os.path.join(output_dir, "parts.dat")

    return result.returncode, eigenvalues_path, parts_path


def check_eigenvalues(material, evals):
    """R1: Check eigenvalues match expected band edges at machine precision.

    Expected ascending order: [-DeltaSO, -DeltaSO, 0, 0, 0, 0, Eg, Eg]
    """
    ref = MATERIALS[material]
    expected = [
        -ref["DeltaSO"], -ref["DeltaSO"],
        0.0, 0.0, 0.0, 0.0,
        ref["Eg"], ref["Eg"],
    ]
    failures = []

    if len(evals) != 8:
        failures.append(f"Expected 8 eigenvalues, got {len(evals)}")
        return failures

    for i, (got, exp) in enumerate(zip(evals, expected)):
        if abs(exp) > 1e-6:
            rel_err = abs(got - exp) / abs(exp)
            if rel_err > TOL_MACHINE:
                failures.append(
                    f"Eigenvalue[{i}]: got {got:.10e}, expected {exp:.6f}, "
                    f"rel_err = {rel_err:.2e}"
                )
        else:
            # Near zero: use absolute tolerance
            if abs(got) > TOL_MACHINE:
                failures.append(
                    f"Eigenvalue[{i}]: got {got:.10e}, expected 0.0, "
                    f"abs_err = {abs(got):.2e}"
                )

    return failures


def check_degeneracies(evals):
    """R2: Verify eigenvalue degeneracies.

    Expected: SO pair (indices 0,1), HH/LH quartet (indices 2,3,4,5),
    CB pair (indices 6,7).
    """
    failures = []

    if len(evals) != 8:
        failures.append("Cannot check degeneracies: need 8 eigenvalues")
        return failures

    # SO pair (indices 0, 1)
    if abs(evals[0] - evals[1]) > TOL_DEGEN:
        failures.append(
            f"SO degeneracy broken: E[0]={evals[0]:.10e}, E[1]={evals[1]:.10e}, "
            f"diff={abs(evals[0]-evals[1]):.2e}"
        )

    # HH/LH quartet (indices 2, 3, 4, 5)
    for i in range(2, 6):
        if abs(evals[i] - evals[2]) > TOL_DEGEN:
            failures.append(
                f"HH/LH quartet degeneracy broken: E[2]={evals[2]:.10e}, "
                f"E[{i}]={evals[i]:.10e}, diff={abs(evals[i]-evals[2]):.2e}"
            )

    # CB pair (indices 6, 7)
    if abs(evals[6] - evals[7]) > TOL_DEGEN:
        failures.append(
            f"CB degeneracy broken: E[6]={evals[6]:.10e}, E[7]={evals[7]:.10e}, "
            f"diff={abs(evals[6]-evals[7]):.2e}"
        )

    return failures


def check_td_symmetry(evals, material):
    """R4: Verify T_d symmetry at k=0 (indirect check).

    At k=0 in a zinc-blende lattice, the Hamiltonian has T_d symmetry.
    Any spurious k-dependent coupling would cause eigenvalue deviations from
    the exact band edges.  Since R1 already checks eigenvalues match the
    diagonal to machine precision, this is implicitly satisfied.
    We verify there are exactly 3 distinct eigenvalue levels.
    """
    failures = []

    if len(evals) != 8:
        failures.append("Cannot check T_d symmetry: need 8 eigenvalues")
        return failures

    # Collect distinct energy levels
    levels = [evals[0], evals[2], evals[6]]  # SO, HH/LH, CB representatives
    # Check they are all distinct (within machine precision)
    for i in range(len(levels)):
        for j in range(i + 1, len(levels)):
            if abs(levels[i] - levels[j]) < TOL_MACHINE:
                failures.append(
                    f"T_d symmetry: energy levels {i} and {j} are unexpectedly "
                    f"degenerate: {levels[i]:.10e} vs {levels[j]:.10e}"
                )

    # Check we have exactly 3 distinct levels (SO, HH/LH, CB)
    unique_levels = []
    for e in evals:
        found = False
        for u in unique_levels:
            if abs(e - u) < TOL_DEGEN:
                found = True
                break
        if not found:
            unique_levels.append(e)

    if len(unique_levels) != 3:
        failures.append(
            f"T_d symmetry: expected 3 distinct eigenvalue levels, got "
            f"{len(unique_levels)}: {[f'{u:.6f}' for u in unique_levels]}"
        )

    return failures


def check_basis_ordering(parts):
    """R3: Verify basis ordering via eigenfunction weights.

    Band ordering: 1=HH, 2=LH, 3=LH, 4=HH, 5=SO, 6=SO, 7=CB, 8=CB.
    Grouped: HH = cols 0+3, LH = cols 1+2, SO = cols 4+5, CB = cols 6+7.

    Expected at k=0:
      States 0,1 (SO): SO weight > TOL_CHARACTER
      States 2,3,4,5 (HH/LH): HH+LH weight > TOL_CHARACTER
      States 6,7 (CB): CB weight > TOL_CHARACTER
    """
    failures = []

    if len(parts) < 8:
        failures.append(
            f"Cannot check basis ordering: need 8 eigenstates in parts.dat, "
            f"got {len(parts)}"
        )
        return failures

    for state_idx in range(8):
        row = parts[state_idx]
        row_sum = sum(row)
        if row_sum == 0:
            failures.append(
                f"State {state_idx}: parts row sums to zero"
            )
            continue
        # Normalize to fractional character
        frac = [r / row_sum for r in row]

        hh_weight = frac[0] + frac[3]
        lh_weight = frac[1] + frac[2]
        so_weight = frac[4] + frac[5]
        cb_weight = frac[6] + frac[7]

        if state_idx < 2:
            # SO states
            if so_weight < TOL_CHARACTER:
                failures.append(
                    f"State {state_idx} (expected SO): SO weight = {so_weight:.6f} "
                    f"< {TOL_CHARACTER}"
                )
        elif state_idx < 6:
            # HH/LH states
            vvb_weight = hh_weight + lh_weight
            if vvb_weight < TOL_CHARACTER:
                failures.append(
                    f"State {state_idx} (expected HH/LH): HH+LH weight = "
                    f"{vvb_weight:.6f} < {TOL_CHARACTER}"
                )
        else:
            # CB states
            if cb_weight < TOL_CHARACTER:
                failures.append(
                    f"State {state_idx} (expected CB): CB weight = {cb_weight:.6f} "
                    f"< {TOL_CHARACTER}"
                )

    return failures


def check_normalization(parts):
    """R5: Verify eigenfunction normalization (unit sum of squares).

    For bulk mode, parts(j, band) = |A(band, j)|^2 and each row should
    sum to 1.0.
    """
    failures = []

    if len(parts) < 8:
        failures.append(
            f"Cannot check normalization: need 8 eigenstates, got {len(parts)}"
        )
        return failures

    for state_idx in range(8):
        row_sum = sum(parts[state_idx])
        if abs(row_sum - 1.0) > TOL_NORM:
            failures.append(
                f"State {state_idx}: normalization = {row_sum:.10e}, "
                f"deviation from 1.0 = {abs(row_sum - 1.0):.2e} "
                f"(tolerance {TOL_NORM:.0e})"
            )

    return failures


def validate_material(material, exe_path, config_path):
    """Run all R1-R5 checks for one material. Returns (passed, failures)."""
    failures = []

    with tempfile.TemporaryDirectory(prefix=f"rung1_{material}_") as tmpdir:
        returncode, eigenvalues_path, parts_path = run_bandstructure(
            exe_path, config_path, tmpdir
        )

        if returncode != 0:
            failures.append(
                f"bandStructure exited with code {returncode}"
            )
            return False, failures

        # Parse eigenvalues
        if not os.path.exists(eigenvalues_path):
            failures.append("eigenvalues.dat not found in output/")
            return False, failures

        evals = parse_eigenvalues(eigenvalues_path)
        if evals is None:
            failures.append("eigenvalues.dat is empty or could not be parsed")
            return False, failures

        print(f"  Eigenvalues: {['%.6f' % e for e in evals]}")

        # R1: Eigenvalue accuracy
        r1_failures = check_eigenvalues(material, evals)
        if r1_failures:
            failures.extend([f"[R1] {f}" for f in r1_failures])
        else:
            print(f"  [R1] Eigenvalue accuracy: PASS (machine precision)")

        # R2: Degeneracies
        r2_failures = check_degeneracies(evals)
        if r2_failures:
            failures.extend([f"[R2] {f}" for f in r2_failures])
        else:
            print(f"  [R2] Degeneracies: PASS (SO 2x, HH/LH 4x, CB 2x)")

        # R4: T_d symmetry (indirect)
        r4_failures = check_td_symmetry(evals, material)
        if r4_failures:
            failures.extend([f"[R4] {f}" for f in r4_failures])
        else:
            print(f"  [R4] T_d symmetry: PASS (3 distinct levels)")

        # R3 + R5: Eigenfunction checks (parts.dat)
        if os.path.exists(parts_path):
            parts = parse_parts(parts_path)
            if parts:
                # R3: Basis ordering
                r3_failures = check_basis_ordering(parts)
                if r3_failures:
                    failures.extend([f"[R3] {f}" for f in r3_failures])
                else:
                    print(f"  [R3] Basis ordering: PASS")

                # R5: Normalization
                r5_failures = check_normalization(parts)
                if r5_failures:
                    failures.extend([f"[R5] {f}" for f in r5_failures])
                else:
                    print(f"  [R5] Normalization: PASS (unit norm)")
            else:
                failures.append("[R3/R5] parts.dat is empty")
        else:
            print(f"  [R3/R5] parts.dat not found — skipping basis ordering "
                  f"and normalization checks")

    return len(failures) == 0, failures


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <repo_dir>")
        print(f"  build_dir  — path to build/ (contains src/bandStructure)")
        print(f"  repo_dir   — path to repo root")
        sys.exit(1)

    build_dir = os.path.abspath(sys.argv[1])
    test_dir = os.path.abspath(sys.argv[2])

    exe_path = os.path.join(build_dir, "src", "bandStructure")
    if not os.path.isfile(exe_path):
        print(f"ERROR: executable not found: {exe_path}")
        sys.exit(1)

    config_dir = os.path.join(test_dir, "tests", "regression", "configs")

    all_pass = True
    total = 0

    for material in MATERIALS:
        config_name = f"bulk_{material.lower()}_k0.toml"
        config_path = os.path.join(config_dir, config_name)

        if not os.path.isfile(config_path):
            print(f"\n{'='*60}")
            print(f"MATERIAL: {material}")
            print(f"  SKIP: config not found: {config_path}")
            continue

        total += 1
        ref = MATERIALS[material]
        print(f"\n{'='*60}")
        print(f"MATERIAL: {material}")
        print(f"  Expected: Eg={ref['Eg']}, DeltaSO={ref['DeltaSO']}")
        print(f"  Config: {config_path}")

        passed, failures = validate_material(material, exe_path, config_path)

        if passed:
            print(f"  >>> ALL CHECKS PASSED for {material}")
        else:
            all_pass = False
            print(f"  >>> FAILURES for {material}:")
            for f in failures:
                print(f"      {f}")

    # Summary
    print(f"\n{'='*60}")
    print(f"SUMMARY: {total} materials tested")
    if all_pass:
        print("PASS: all bulk k=0 structural validation checks passed")
        sys.exit(0)
    else:
        print("FAIL: one or more materials failed validation")
        sys.exit(1)


if __name__ == "__main__":
    main()
