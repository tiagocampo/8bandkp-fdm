"""Bulk k=0 cross-code validation — the gate test.

Compares all 8 eigenvalues at Gamma for all available III-V materials
between our Fortran code and kdotpy. This validates the entire parameter
mapping pipeline.

Tolerance: 0.01 meV (parameter mapping should be exact at k=0).
"""

import os
import sys
import json
import tempfile

# Add project root to path
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, project_root)

from validation.shared.param_mapper import map_material, list_materials
from validation.shared.kdotpy_runner import run_bulk_single
from validation.shared.fortran_runner import run_bulk
from validation.shared.comparison import compare_eigenvalues, TOL_EXACT, write_json_report, print_summary

BUILD_DIR = os.path.join(project_root, "build")

# Materials to test (skip Vacuum, Al, and W-variants for initial validation)
BULK_MATERIALS = [
    "GaAs", "InAs", "InSb", "AlAs", "GaSb", "AlSb", "InP",
    "GaAsW", "InAsW", "InSbW",
    "Al20Ga80As", "Al30Ga70As",
]


def test_all_bulk_k0():
    """Run bulk k=0 comparison for all materials."""
    results = []
    failures = []

    print("=" * 70)
    print("BULK k=0 CROSS-CODE VALIDATION (Gate Test)")
    print("=" * 70)
    print(f"Tolerance: {TOL_EXACT} meV")
    print(f"Materials: {len(BULK_MATERIALS)}")
    print()

    print("| Material | CB (meV) | VB top (meV) | SO (meV) | Max delta | Status |")
    print("|----------|----------|-------------|----------|-----------|--------|")

    for mat_name in BULK_MATERIALS:
        try:
            # Run kdotpy
            kdotpy_mat = map_material(mat_name)
            kdotpy_evals = run_bulk_single(kdotpy_mat)

            # Run Fortran
            fortran_results = run_bulk(BUILD_DIR, mat_name, [(0.0, 0.0, 0.0)])
            fortran_eV = fortran_results[0][1]

            # Compare
            comp = compare_eigenvalues(fortran_eV, kdotpy_evals, tolerance_meV=TOL_EXACT)

            # Extract key energies for summary
            cb_meV = comp["per_band"][-1]["ours_meV"]
            vb_meV = comp["per_band"][5]["ours_meV"]  # VB top
            so_meV = comp["per_band"][0]["ours_meV"]   # SO bottom

            status = "PASS" if comp["passed"] else "FAIL"
            print(f"| {mat_name} | {cb_meV:.2f} | {vb_meV:.2f} | "
                  f"{so_meV:.2f} | {comp['max_delta_meV']:.6f} | {status} |")

            if not comp["passed"]:
                failures.append(mat_name)
                # Print per-band details for failures
                for b in comp["per_band"]:
                    if not b["passed"]:
                        print(f"  -> Band {b['band']}: ours={b['ours_meV']:.6f} "
                              f"kdotpy={b['kdotpy_meV']:.6f} delta={b['delta_meV']:.6f}")

            results.append({
                "material": mat_name,
                "test": "bulk_k0",
                "comparison": comp,
            })

        except Exception as e:
            print(f"| {mat_name} | ERROR | ERROR | ERROR | ERROR | ERROR |")
            print(f"  -> {e}")
            failures.append(mat_name)

    # Save results
    results_dir = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(results_dir, exist_ok=True)
    write_json_report(results, os.path.join(results_dir, "bulk_k0.json"))

    print()
    print(f"Summary: {len(results) - len(failures)}/{len(results)} materials passed")
    if failures:
        print(f"FAILURES: {failures}")
        return False
    return True


if __name__ == "__main__":
    success = test_all_bulk_k0()
    sys.exit(0 if success else 1)
