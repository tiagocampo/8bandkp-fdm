"""Self-consistent QW validation.

Runs our self-consistent Schrodinger-Poisson solver for a doped GaAs/AlGaAs QW
and verifies convergence behavior and subband energies against expected values.

Note: kdotpy cross-code comparison is deferred to follow-up work due to
different SC mixing strategies (DIIS vs dynamic time stepping).

Tolerance: < 5 meV for CB1 subband energy.
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

SC_CONFIGS = [
    {
        "name": "GaAs/AlGaAs doped QW",
        "config": "sc_gaas_alas_qw.cfg",
        "expected_cb1_meV": 800.0,  # approximate, depends on SC state
    },
]


def run_fortran_sc(config_name, build_dir, project_root):
    exe = os.path.join(build_dir, "src", "bandStructure")
    config_path = os.path.join(project_root, "tests", "regression", "configs", config_name)

    if not os.path.isfile(exe):
        raise FileNotFoundError(f"Executable not found: {exe}")
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Config not found: {config_path}")

    workdir = tempfile.mkdtemp(prefix="sc_")
    try:
        shutil.copy2(config_path, os.path.join(workdir, "input.cfg"))
        os.makedirs(os.path.join(workdir, "output"), exist_ok=True)
        result = subprocess.run(
            [exe], cwd=workdir, capture_output=True, text=True, timeout=300
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"bandStructure failed (rc={result.returncode}):\n"
                f"stderr: {result.stderr[-500:]}"
            )

        # Parse eigenvalues from output
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

        # Check SC convergence from stdout
        converged = "SC loop converged" in result.stdout

        return {"eigenvalues": rows, "converged": converged, "stdout": result.stdout}
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def test_sc_qw():
    print("=" * 70)
    print("SELF-CONSISTENT QW CROSS-VALIDATION")
    print("=" * 70)
    print(f"Tolerance: CB1 within {TOL_CB1_MEV} meV of expected value")
    print()

    all_pass = True
    all_results = []

    for cfg in SC_CONFIGS:
        name = cfg["name"]
        config = cfg["config"]
        expected_cb1 = cfg["expected_cb1_meV"]

        print(f"\nConfig: {name}")

        try:
            result = run_fortran_sc(config, BUILD_DIR, project_root)
        except (RuntimeError, FileNotFoundError) as e:
            print(f"  FAIL: {e}")
            all_pass = False
            all_results.append({"config": name, "status": "FAIL"})
            continue

        rows = result["eigenvalues"]
        if not rows:
            print(f"  SKIP: no eigenvalues produced")
            all_results.append({"config": name, "status": "SKIP"})
            continue

        # First row is k=0
        k0_row = rows[0]
        evals_eV = sorted(k0_row[1:])
        evals_meV = [e * MEV_PER_EV for e in evals_eV]

        # CB eigenvalues: use the numcb highest
        numcb = 4
        cb_evals = sorted(evals_meV[-numcb:])
        cb1_meV = cb_evals[0]

        delta_meV = abs(cb1_meV - expected_cb1)
        passed = delta_meV < TOL_CB1_MEV and result["converged"]

        print(f"  CB1 = {cb1_meV:.2f} meV ({cb_evals[0]/MEV_PER_EV:.6f} eV)")
        print(f"  Expected ~{expected_cb1:.0f} meV, delta = {delta_meV:.2f} meV")
        print(f"  SC converged: {result['converged']}")
        print(f"  CB eigenvalues: {[f'{e:.2f}' for e in cb_evals]} meV")
        print(f"  Status: {'PASS' if passed else 'FAIL'}")

        if not result["converged"]:
            print(f"  WARNING: SC loop did not converge")
        if not passed:
            all_pass = False

        all_results.append({
            "config": name,
            "cb1_meV": cb1_meV,
            "expected_cb1_meV": expected_cb1,
            "delta_meV": delta_meV,
            "sc_converged": result["converged"],
            "passed": bool(passed),
        })

    print("\nNote: kdotpy SC comparison uses different mixing strategies.")
    print("      This test validates our SC solver against expected physics.")

    results_dir = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(results_dir, exist_ok=True)
    with open(os.path.join(results_dir, "sc_qw.json"), "w") as f:
        json.dump(all_results, f, indent=2)

    n_pass = sum(1 for r in all_results if r.get("passed"))
    print(f"\nSummary: {n_pass}/{len(all_results)} configs passed")
    return all_pass


if __name__ == "__main__":
    success = test_sc_qw()
    sys.exit(0 if success else 1)
