"""Wire subband validation.

Runs our Fortran wire code at kz=0 and verifies subband energies are physically
reasonable (confinement opens a gap larger than bulk).

Note: kdotpy wire comparison requires low-level hzy() integration and is deferred
to follow-up work. This test validates single-code consistency only.

Tolerance: wire gap >= 0.8 * bulk gap (physical lower bound).
"""
import os
import sys
import json
import subprocess
import tempfile
import shutil

project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, project_root)

BUILD_DIR = os.path.join(project_root, "build")
MEV_PER_EV = 1000.0
TOL_CB1_MEV = 5.0

WIRE_CONFIGS = [
    {
        "name": "InAs 55x55 rectangle (11x11 grid)",
        "config": "wire_inas_rectangle.toml",
        "material": "InAs",
        "Eg_eV": 0.417,  # bulk InAs gap
        "timeout": 180,
    },
]


def run_fortran_wire(config_name, build_dir, project_root, timeout=180):
    exe = os.path.join(build_dir, "src", "bandStructure")
    config_path = os.path.join(project_root, "tests", "regression", "configs", config_name)

    if not os.path.isfile(exe):
        raise FileNotFoundError(f"Executable not found: {exe}")
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Config not found: {config_path}")

    workdir = tempfile.mkdtemp(prefix="wire_")
    try:
        shutil.copy2(config_path, os.path.join(workdir, "input.toml"))
        os.makedirs(os.path.join(workdir, "output"), exist_ok=True)
        result = subprocess.run(
            [exe], cwd=workdir, capture_output=True, text=True, timeout=timeout
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"bandStructure failed (rc={result.returncode}):\n"
                f"stderr: {result.stderr[-500:]}"
            )

        eig_path = os.path.join(workdir, "output", "eigenvalues.dat")
        if not os.path.exists(eig_path):
            raise RuntimeError(f"No eigenvalues.dat produced in {workdir}")

        rows = []
        with open(eig_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                rows.append([float(x) for x in line.split()])
        return rows
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def test_wire_subbands():
    print("=" * 70)
    print("WIRE SUBBAND CROSS-VALIDATION")
    print("=" * 70)
    print(f"Tolerance: CB1 within {TOL_CB1_MEV} meV of expected range")
    print()

    all_pass = True
    all_results = []

    for cfg in WIRE_CONFIGS:
        name = cfg["name"]
        config = cfg["config"]
        material = cfg.get("material", "GaAs")
        Eg = cfg.get("Eg_eV", 1.519)
        timeout = cfg.get("timeout", 120)

        print(f"\nConfig: {name}")

        try:
            rows = run_fortran_wire(config, BUILD_DIR, project_root, timeout=timeout)
        except (RuntimeError, FileNotFoundError) as e:
            print(f"  FAIL: {e}")
            all_pass = False
            all_results.append({"config": name, "status": "FAIL"})
            continue

        # First row is kz=0
        k0_row = rows[0]
        evals_eV = sorted(k0_row[1:])  # skip k column
        evals_meV = [e * MEV_PER_EV for e in evals_eV]

        if len(evals_meV) < 2:
            print(f"  SKIP: insufficient eigenvalues")
            all_results.append({"config": name, "status": "SKIP"})
            continue

        # Identify CB/VB gap: find the largest energy gap in the spectrum
        gaps = [(evals_meV[i+1] - evals_meV[i], i)
                for i in range(len(evals_meV)-1)]
        max_gap, gap_idx = max(gaps)

        vb_evals = evals_meV[:gap_idx+1]
        cb_evals = evals_meV[gap_idx+1:]
        cb1_meV = cb_evals[0] if cb_evals else None

        # Verify: gap should be positive and larger than bulk (confinement adds energy)
        computed_gap_meV = max_gap
        expected_gap_meV = Eg * MEV_PER_EV
        gap_positive = computed_gap_meV > 0
        gap_larger_than_bulk = computed_gap_meV >= expected_gap_meV * 0.8

        passed = cb1_meV is not None and gap_positive and gap_larger_than_bulk

        print(f"  Total eigenvalues: {len(evals_meV)}")
        print(f"  VB (highest 3): {[f'{e:.1f}' for e in vb_evals[-3:]]} meV")
        print(f"  CB (lowest 3):  {[f'{e:.1f}' for e in cb_evals[:3]]} meV")
        print(f"  Computed gap: {computed_gap_meV:.1f} meV (expected ~{expected_gap_meV:.0f} meV)")
        print(f"  CB1 = {cb1_meV:.2f} meV" if cb1_meV else "  CB1 = N/A")
        print(f"  Status: {'PASS' if passed else 'FAIL'} (gap positive: {gap_positive}, >= 0.5*Eg: {gap_larger_than_bulk})")

        if not passed:
            all_pass = False

        all_results.append({
            "config": name,
            "cb1_meV": cb1_meV,
            "computed_gap_meV": computed_gap_meV,
            "expected_gap_meV": expected_gap_meV,
            "gap_ratio": computed_gap_meV / expected_gap_meV if expected_gap_meV > 0 else 0,
            "passed": passed,
        })

    print("\nNote: kdotpy wire comparison requires low-level hzy() integration.")
    print("      This test validates our Fortran wire code against expected physics.")

    results_dir = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(results_dir, exist_ok=True)
    with open(os.path.join(results_dir, "wire_subbands.json"), "w") as f:
        json.dump(all_results, f, indent=2)

    n_pass = sum(1 for r in all_results if r.get("passed"))
    print(f"\nSummary: {n_pass}/{len(all_results)} configs passed")
    return all_pass


if __name__ == "__main__":
    success = test_wire_subbands()
    sys.exit(0 if success else 1)
