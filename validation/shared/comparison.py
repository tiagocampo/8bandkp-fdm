"""Comparison and reporting utilities for cross-code validation.

Compares eigenvalues between our Fortran code (eV) and kdotpy (meV),
handles unit conversion and band-index alignment.
"""

import json
import os
import sys

import numpy as np


# Tolerance tiers in meV
TOL_EXACT = 0.01       # Parameter mapping validation
TOL_BULK_DISPERSION = 0.1
TOL_QW_SUBBAND = 1.0
TOL_QW_CONVERGED = 0.5
TOL_WIRE = 2.0
TOL_LANDAU = 1.0
TOL_GFACTOR_PCT = 5.0   # percent
TOL_BHZ_PCT = 10.0      # percent
TOL_STRAIN = 1.0
TOL_BERRY_PCT = 10.0    # percent
TOL_SC_SUBBAND = 5.0
TOL_SC_CHARGE_PCT = 10.0

EV_TO_MEV = 1000.0


def compare_eigenvalues(fortran_eV, kdotpy_meV, tolerance_meV=TOL_EXACT):
    """Compare eigenvalues from both codes.

    Both arrays must have the same length. Fortran eigenvalues are sorted
    ascending (our output format). kdotpy eigenvalues are sorted ascending.

    Args:
        fortran_eV: array-like, eigenvalues from our code in eV
        kdotpy_meV: array-like, eigenvalues from kdotpy in meV
        tolerance_meV: maximum allowed absolute difference in meV

    Returns:
        dict with keys:
          passed: bool
          max_delta_meV: float
          per_band: list of {band, ours_meV, kdotpy_meV, delta_meV, passed}
    """
    ours_meV = np.array(fortran_eV) * EV_TO_MEV
    theirs_meV = np.array(kdotpy_meV)

    n = min(len(ours_meV), len(theirs_meV))
    per_band = []
    for i in range(n):
        delta = abs(float(ours_meV[i]) - float(theirs_meV[i]))
        per_band.append({
            "band": i,
            "ours_meV": float(ours_meV[i]),
            "kdotpy_meV": float(theirs_meV[i]),
            "delta_meV": delta,
            "passed": delta <= tolerance_meV,
        })

    max_delta = max(b["delta_meV"] for b in per_band) if per_band else 0.0
    passed = all(b["passed"] for b in per_band)

    return {
        "passed": passed,
        "max_delta_meV": max_delta,
        "per_band": per_band,
    }


def format_report(comparison, material, test_name, tolerance_meV):
    """Format a comparison result as a markdown table row.

    Args:
        comparison: result from compare_eigenvalues
        material: material name
        test_name: test identifier
        tolerance_meV: tolerance used

    Returns:
        list of markdown table row strings
    """
    rows = []
    status = "PASS" if comparison["passed"] else "FAIL"
    for b in comparison["per_band"]:
        rows.append(
            f"| {material} | {test_name} band {b['band']} | "
            f"{b['ours_meV']:.4f} | {b['kdotpy_meV']:.4f} | "
            f"{b['delta_meV']:.6f} | {tolerance_meV} | "
            f"{'PASS' if b['passed'] else 'FAIL'} |"
        )
    return rows


def write_json_report(results, filepath):
    """Write comparison results to JSON file.

    Args:
        results: list of comparison dicts
        filepath: output path
    """
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        json.dump(results, f, indent=2)


def print_summary(results):
    """Print a summary table of all comparison results.

    Args:
        results: list of dicts with 'material', 'test', 'comparison' keys
    """
    print("\n| Material | Test | Max Delta (meV) | Status |")
    print("|----------|------|-----------------|--------|")
    for r in results:
        comp = r["comparison"]
        status = "PASS" if comp["passed"] else "FAIL"
        print(f"| {r['material']} | {r['test']} | "
              f"{comp['max_delta_meV']:.6f} | {status} |")

    total = len(results)
    passed = sum(1 for r in results if r["comparison"]["passed"])
    print(f"\nSummary: {passed}/{total} tests passed")
    return passed == total
