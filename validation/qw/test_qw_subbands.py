"""QW subband cross-code validation.

Compares QW subband energies at k_par=0 between our Fortran code and kdotpy.
After the const fix (confinement_init.f90 applies const to gamma/A profiles
and hamiltonianConstructor.f90 applies const in the bulk Hamiltonian), the
remaining discrepancy is from grid resolution differences.

Tolerance: 2 meV for CB1 subband energy (accounts for FD vs plane-wave
discretization differences at finite grid resolution).
"""

import os
import sys
import json

import numpy as np

project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, project_root)

from validation.shared.param_mapper import map_material, FORTRAN_MATERIALS
from validation.shared.kdotpy_runner import run_qw
from validation.shared.fortran_runner import run_qw as fortran_run_qw

BUILD_DIR = os.path.join(project_root, "build")
ANG_TO_NM = 10.0
TOL_QW_CB = 2.0  # meV
TOL_QW_VB = 3.0  # meV (VB more sensitive to grid)

QW_TESTS = [
    {"barrier": "Al20Ga80As", "well": "GaAs", "l_well_nm": 10.0,
     "l_barrier_nm": 15.0, "fdstep": 201, "zres": 0.25},
    {"barrier": "Al20Ga80As", "well": "GaAs", "l_well_nm": 7.0,
     "l_barrier_nm": 15.0, "fdstep": 201, "zres": 0.25},
    {"barrier": "Al30Ga70As", "well": "GaAs", "l_well_nm": 10.0,
     "l_barrier_nm": 15.0, "fdstep": 201, "zres": 0.25},
]


def get_unique_evals(evals, threshold=0.5):
    """Remove Kramers degeneracy from eigenvalue list."""
    unique = []
    for e in sorted(evals):
        if not unique or abs(e - unique[-1]) > threshold:
            unique.append(e)
    return unique


def classify_bands(unique_evals, well_material):
    """Split eigenvalues into CB and VB using material-aware mid-gap."""
    well_params = FORTRAN_MATERIALS[well_material]
    ec_meV = well_params["EC"] * 1000.0
    eg_meV = well_params["Eg"] * 1000.0
    ev_meV = ec_meV - eg_meV
    midgap = (ec_meV + ev_meV) / 2.0

    cb = [e for e in unique_evals if e > midgap]
    vb = sorted([e for e in unique_evals if e <= midgap], reverse=True)
    return cb, vb


def test_qw_subbands():
    """Compare QW subband energies between codes."""
    results = []
    all_pass = True

    print("=" * 70)
    print("QW SUBBAND CROSS-CODE VALIDATION")
    print("=" * 70)
    print(f"Tolerance: {TOL_QW_CB} meV for CB1")
    print()

    # Check kdotpy availability
    try:
        from kdotpy.config import initialize_config
        initialize_config()
    except ImportError:
        print("SKIP: kdotpy not available (activate kdotpy_env first)")
        return False

    print("| Config | CB1 ours | CB1 kd | Delta | VB1 ours | VB1 kd | Delta | Status |")
    print("|--------|----------|--------|-------|----------|--------|-------|--------|")

    for cfg in QW_TESTS:
        barrier = cfg["barrier"]
        well = cfg["well"]
        l_well_nm = cfg["l_well_nm"]
        l_barrier_nm = cfg["l_barrier_nm"]
        fdstep = cfg["fdstep"]
        zres = cfg["zres"]

        l_well_ang = l_well_nm * ANG_TO_NM
        l_barrier_ang = l_barrier_nm * ANG_TO_NM
        total_ang = 2 * l_barrier_ang + l_well_ang

        label = f"{well}/{barrier} {l_well_nm}nm"

        try:
            f_results = fortran_run_qw(
                BUILD_DIR, barrier, well,
                l_well_ang, l_barrier_ang, total_ang,
                [(0.0, 0.0)], fdstep=fdstep, numcb=6, numvb=12,
            )
            f_evals = f_results[0][1]
            f_unique = get_unique_evals([e * 1000 for e in f_evals])
        except (FileNotFoundError, RuntimeError, OSError) as e:
            print(f"| {label} | ERROR ({type(e).__name__}: {e}) | ERROR | - | ERROR | ERROR | - | ERROR |")
            all_pass = False
            continue

        try:
            barrier_mat = map_material(barrier, qw_mode=True)
            well_mat = map_material(well, qw_mode=True)
            kd_results = run_qw(
                barrier_mat, well_mat, l_well_nm, l_barrier_nm,
                [(0.0, 0.0)], zres=zres, neig=50, energy=700.0,
            )
            kd_evals = sorted(kd_results[0])
            kd_unique = get_unique_evals(kd_evals)
        except (ImportError, RuntimeError, ValueError, OSError) as e:
            print(f"| {label} | {f_unique[-1]:.2f} | ERROR ({type(e).__name__}: {e}) | - | "
                  f"{f_unique[0] if f_unique else 'N/A'} | ERROR | - | ERROR |")
            all_pass = False
            continue

        f_cb, f_vb = classify_bands(f_unique, well)
        kd_cb, kd_vb = classify_bands(kd_unique, well)

        if not f_cb or not kd_cb:
            print(f"| {label} | MISSING | MISSING | - | MISSING | MISSING | - | ERROR |")
            all_pass = False
            continue

        cb_delta = abs(f_cb[0] - kd_cb[0])
        vb_delta = abs(f_vb[0] - kd_vb[0]) if (f_vb and kd_vb) else 0.0
        cb_passed = cb_delta < TOL_QW_CB
        vb_passed = vb_delta < TOL_QW_VB if (f_vb and kd_vb) else True

        if not cb_passed:
            all_pass = False
        if not vb_passed:
            all_pass = False

        status = "PASS" if (cb_passed and vb_passed) else "FAIL"
        print(f"| {label} | {f_cb[0]:.2f} | {kd_cb[0]:.2f} | "
              f"{cb_delta:.2f} | {f_vb[0]:.2f} | {kd_vb[0]:.2f} | "
              f"{vb_delta:.2f} | {status} |")

        results.append({
            "config": label,
            "barrier": barrier,
            "well": well,
            "l_well_nm": l_well_nm,
            "cb1_fortran_meV": float(f_cb[0]),
            "cb1_kdotpy_meV": float(kd_cb[0]),
            "cb1_delta_meV": float(cb_delta),
            "vb1_fortran_meV": float(f_vb[0]) if f_vb else None,
            "vb1_kdotpy_meV": float(kd_vb[0]) if kd_vb else None,
            "vb1_delta_meV": float(vb_delta),
            "passed": bool(cb_passed and vb_passed),
            "cb_passed": bool(cb_passed),
            "vb_passed": bool(vb_passed),
        })

    results_dir = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(results_dir, exist_ok=True)
    with open(os.path.join(results_dir, "qw_subbands.json"), "w") as f:
        json.dump(results, f, indent=2)

    passed_count = sum(1 for r in results if r.get("passed", False))
    print(f"\nSummary: {passed_count}/{len(results)} passed (CB1 < {TOL_QW_CB} meV)")
    return all_pass


if __name__ == "__main__":
    success = test_qw_subbands()
    sys.exit(0 if success else 1)
