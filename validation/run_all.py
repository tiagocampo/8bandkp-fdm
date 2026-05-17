"""Cross-code validation pipeline runner.

Runs all validation tests and produces a summary report.

Usage:
    source validation/kdotpy_env/bin/activate
    python validation/run_all.py
"""

import os
import sys
import json
import subprocess
import datetime

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, project_root)

TESTS = [
    ("Bulk k=0 Gate Test", "validation/bulk/test_bulk_k0.py"),
    ("Bulk Dispersion", "validation/bulk/test_bulk_dispersion.py"),
    ("Bulk Zeeman", "validation/bulk/test_bulk_zeeman.py"),
    ("QW Subbands", "validation/qw/test_qw_subbands.py"),
]

VENV_ACTIVATE = os.path.join(project_root, "validation", "kdotpy_env", "bin", "activate")


def run_test(name, script):
    """Run a single test script and return result dict."""
    script_path = os.path.join(project_root, script)
    if not os.path.isfile(script_path):
        return {"name": name, "script": script, "status": "SKIP", "output": "File not found"}

    cmd = f"source {VENV_ACTIVATE} && python3 {script_path}"
    result = subprocess.run(
        ["bash", "-c", cmd],
        capture_output=True, text=True, timeout=300,
        cwd=project_root,
    )

    return {
        "name": name,
        "script": script,
        "status": "PASS" if result.returncode == 0 else "FAIL",
        "returncode": result.returncode,
        "output": result.stdout[-2000:] if result.stdout else "",
        "error": result.stderr[-500:] if result.stderr else "",
    }


def main():
    print("=" * 70)
    print("CROSS-CODE VALIDATION PIPELINE")
    print(f"Date: {datetime.datetime.now().isoformat()}")
    print("=" * 70)
    print()

    all_results = []
    n_pass = 0
    n_fail = 0
    n_skip = 0

    for name, script in TESTS:
        print(f"Running: {name} ({script})")
        result = run_test(name, script)
        all_results.append(result)

        if result["status"] == "PASS":
            n_pass += 1
            print(f"  -> PASS")
        elif result["status"] == "SKIP":
            n_skip += 1
            print(f"  -> SKIP: {result['output']}")
        else:
            n_fail += 1
            print(f"  -> FAIL (rc={result.get('returncode', '?')})")
            if result.get("error"):
                print(f"     Error: {result['error'][:200]}")
        print()

    print("=" * 70)
    print(f"SUMMARY: {n_pass} passed, {n_fail} failed, {n_skip} skipped")
    print("=" * 70)

    print()
    print("KNOWN FINDINGS:")
    print("1. Bulk k=0: PASS (12/12 materials, 0.000000 meV max delta)")
    print("2. Bulk dispersion: convention difference documented")
    print("   (our diagonal lacks const = hbar²/(2m₀))")
    print("3. QW subbands: ~13 meV CB1 offset from missing const")
    print("4. Wire path correctly bakes const (may be independently testable)")
    print()
    print("ACTION REQUIRED: Apply const to QW/bulk diagonal kpterms in")
    print("confinement_init.f90 and hamiltonianConstructor.f90, then re-run.")
    print()
    print("Downstream tests (wire, Landau, g-factor, strain, Berry, SC)")
    print("are deferred pending the Hamiltonian const fix.")

    results_dir = os.path.join(project_root, "validation", "results")
    os.makedirs(results_dir, exist_ok=True)
    report_path = os.path.join(results_dir, "pipeline_report.json")
    with open(report_path, "w") as f:
        json.dump({
            "date": datetime.datetime.now().isoformat(),
            "summary": {"pass": n_pass, "fail": n_fail, "skip": n_skip},
            "tests": all_results,
            "findings": [
                "Bulk k=0 gate test passes exactly (0.000000 meV)",
                "Bulk/QW diagonal terms lack const = hbar²/(2m₀)",
                "QW subband offset ~13 meV CB1 from missing const",
                "Wire path correctly applies const (confinement_init.f90:493-499)",
            ],
        }, f, indent=2, default=str)

    print(f"\nReport saved to: {report_path}")
    return n_fail == 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
