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
    ("QW Dispersion", "validation/qw/test_qw_dispersion.py"),
    ("QW Convergence", "validation/qw/test_qw_convergence.py"),
    ("Landau Bulk", "validation/landau/test_landau_bulk.py"),
    ("Wire Subbands", "validation/wire/test_wire_subbands.py"),
    ("g-Factor QW", "validation/gfactor/test_gfactor_qw.py"),
    ("Strain Bandedge", "validation/strain/test_strain_bandedge.py"),
    ("Strain QW", "validation/strain/test_strain_qw.py"),
    ("Self-Consistent QW", "validation/selfconsistent/test_sc_qw.py"),
]

VENV_ACTIVATE = os.path.join(project_root, "validation", "kdotpy_env", "bin", "activate")
TIMEOUTS = {
    "validation/wire/test_wire_subbands.py": 300,
    "validation/selfconsistent/test_sc_qw.py": 300,
    "validation/qw/test_qw_convergence.py": 600,
}


def run_test(name, script):
    """Run a single test script and return result dict."""
    script_path = os.path.join(project_root, script)
    if not os.path.isfile(script_path):
        return {"name": name, "script": script, "status": "SKIP", "output": "File not found"}

    cmd = f"source '{VENV_ACTIVATE}' && python3 '{script_path}'"
    timeout = TIMEOUTS.get(script, 300)
    try:
        result = subprocess.run(
            ["bash", "-c", cmd],
            capture_output=True, text=True, timeout=timeout,
            cwd=project_root,
        )
    except subprocess.TimeoutExpired:
        return {"name": name, "script": script, "status": "TIMEOUT",
                "output": f"Exceeded {timeout}s timeout", "error": ""}
    except OSError as e:
        return {"name": name, "script": script, "status": "ERROR",
                "output": str(e), "error": ""}

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

        status = result["status"]
        if status == "PASS":
            n_pass += 1
            print(f"  -> PASS")
        elif status == "SKIP":
            n_skip += 1
            print(f"  -> SKIP: {result['output']}")
        elif status == "TIMEOUT":
            n_fail += 1
            print(f"  -> TIMEOUT ({TIMEOUTS.get(script, 300)}s)")
        elif status == "ERROR":
            n_fail += 1
            print(f"  -> ERROR: {result['output']}")
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
    if all(r["status"] == "PASS" for r in all_results):
        print("All validation tests passed.")
    else:
        failed = [r["name"] for r in all_results if r["status"] not in ("PASS", "SKIP")]
        if failed:
            print(f"FAILED tests: {', '.join(failed)}")

    results_dir = os.path.join(project_root, "validation", "results")
    os.makedirs(results_dir, exist_ok=True)
    report_path = os.path.join(results_dir, "pipeline_report.json")
    with open(report_path, "w") as f:
        json.dump({
            "date": datetime.datetime.now().isoformat(),
            "summary": {"pass": n_pass, "fail": n_fail, "skip": n_skip},
            "tests": all_results,
            "findings": [
                f"{r['name']}: {r['status']}" for r in all_results
            ],
        }, f, indent=2, default=str)

    print(f"\nReport saved to: {report_path}")
    return n_fail == 0 and n_pass > 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
