"""QW convergence study with Richardson extrapolation.

Runs both our Fortran code and kdotpy at multiple grid resolutions,
extracts CB1 at k_par=0, performs Richardson extrapolation to the
continuum limit, and compares the extrapolated values.

Tolerance: < 0.5 meV for Richardson-extrapolated CB1 energy.
"""

import os
import sys
import json

import numpy as np

project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, project_root)

from validation.shared.param_mapper import map_material
from validation.shared.kdotpy_runner import run_qw
from validation.shared.fortran_runner import run_qw as fortran_run_qw
from validation.qw.test_qw_subbands import get_unique_evals, classify_bands
from validation.shared.param_mapper import FORTRAN_MATERIALS

BUILD_DIR = os.path.join(project_root, "build")
ANG_TO_NM = 10.0
MEV_PER_EV = 1000.0
TOL_RICHARDSON_MEV = 1.0  # 1 meV for extrapolated CB1 (FD vs plane-wave)

QW_CONFIG = {
    "barrier": "Al20Ga80As", "well": "GaAs",
    "l_well_nm": 10.0, "l_barrier_nm": 15.0,
}

FORTRAN_RESOLUTIONS = [
    {"fdstep": 81, "label": "N=81"},
    {"fdstep": 121, "label": "N=121"},
    {"fdstep": 161, "label": "N=161"},
    {"fdstep": 201, "label": "N=201"},
    {"fdstep": 301, "label": "N=301"},
]

KDOTPY_RESOLUTIONS = [
    {"zres": 0.5, "label": "dz=0.50nm"},
    {"zres": 0.25, "label": "dz=0.25nm"},
    {"zres": 0.125, "label": "dz=0.125nm"},
    {"zres": 0.1, "label": "dz=0.10nm"},
]


def extract_cb1_k0(evals_eV, well_material):
    """Extract CB1 energy at k=0 from eigenvalue list (in eV)."""
    evals_meV = [e * MEV_PER_EV for e in evals_eV]
    unique = get_unique_evals(sorted(evals_meV))
    cb, vb = classify_bands(unique, well_material)
    return cb[0] if cb else None


def richardson_extrapolate(h_vals, E_vals, order=2):
    """Richardson extrapolation to h=0.

    Args:
        h_vals: grid spacing values (ascending or descending)
        E_vals: observable values at each grid spacing
        order: convergence order (2 for FD order 2, typically p=2)

    Returns:
        Extrapolated value at h=0
    """
    h = np.array(h_vals, dtype=float)
    E = np.array(E_vals, dtype=float)

    if len(h) < 2:
        return E[0]

    # Richardson extrapolation: E_ext = E(h1) + (E(h1) - E(h2)) / ((h2/h1)^p - 1)
    # Use the two finest grids
    h1, h2 = h[-1], h[-2]
    E1, E2 = E[-1], E[-2]
    p = order
    E_ext = E1 + (E1 - E2) / ((h2 / h1)**p - 1)
    return E_ext


def test_qw_convergence():
    """Run convergence study and Richardson extrapolation for both codes."""
    barrier = QW_CONFIG["barrier"]
    well = QW_CONFIG["well"]
    l_well_nm = QW_CONFIG["l_well_nm"]
    l_barrier_nm = QW_CONFIG["l_barrier_nm"]
    l_well_ang = l_well_nm * ANG_TO_NM
    l_barrier_ang = l_barrier_nm * ANG_TO_NM
    total_ang = 2 * l_barrier_ang + l_well_ang

    print("=" * 70)
    print("QW CONVERGENCE STUDY WITH RICHARDSON EXTRAPOLATION")
    print("=" * 70)
    print(f"QW: {well}/{barrier} {l_well_nm}nm")
    print(f"Tolerance: extrapolated CB1 delta < {TOL_RICHARDSON_MEV} meV")
    print()

    # --- Fortran convergence ---
    print("Fortran FD convergence:")
    print(f"  {'Resolution':>12s}  {'CB1 (meV)':>12s}  {'h (A)':>10s}")
    f_cb1_list = []
    f_h_list = []

    for res in FORTRAN_RESOLUTIONS:
        fdstep = res["fdstep"]
        h_ang = total_ang / (fdstep - 1)
        try:
            f_results = fortran_run_qw(
                BUILD_DIR, barrier, well,
                l_well_ang, l_barrier_ang, total_ang,
                [(0.0, 0.0)], fdstep=fdstep, numcb=6, numvb=12,
            )
            cb1 = extract_cb1_k0(f_results[0][1], well)
            if cb1 is not None:
                f_cb1_list.append(cb1)
                f_h_list.append(h_ang)
                print(f"  {res['label']:>12s}  {cb1:>12.3f}  {h_ang:>10.3f}")
        except Exception as e:
            print(f"  {res['label']:>12s}  ERROR: {e}")

    # --- kdotpy convergence ---
    print("\nkdotpy plane-wave convergence:")
    print(f"  {'Resolution':>12s}  {'CB1 (meV)':>12s}  {'dz (nm)':>10s}")
    kd_cb1_list = []
    kd_h_list = []

    barrier_mat = map_material(barrier, qw_mode=True)
    well_mat = map_material(well, qw_mode=True)

    for res in KDOTPY_RESOLUTIONS:
        zres = res["zres"]
        try:
            kd_results = run_qw(
                barrier_mat, well_mat, l_well_nm, l_barrier_nm,
                [(0.0, 0.0)], zres=zres, neig=50, energy=700.0,
            )
            evals_meV = sorted(kd_results[0])
            unique = get_unique_evals(evals_meV)
            cb, vb = classify_bands(unique, well)
            if cb:
                kd_cb1_list.append(cb[0])
                kd_h_list.append(zres)
                print(f"  {res['label']:>12s}  {cb[0]:>12.3f}  {zres:>10.4f}")
        except Exception as e:
            print(f"  {res['label']:>12s}  ERROR: {e}")

    # --- Richardson extrapolation ---
    f_extrap = None
    kd_extrap = None

    if len(f_cb1_list) >= 2:
        f_extrap = richardson_extrapolate(f_h_list, f_cb1_list, order=2)
        print(f"\nFortran Richardson extrapolated CB1: {f_extrap:.3f} meV")

    if len(kd_cb1_list) >= 2:
        kd_extrap = richardson_extrapolate(kd_h_list, kd_cb1_list, order=2)
        print(f"kdotpy Richardson extrapolated CB1:  {kd_extrap:.3f} meV")

    # --- Compare extrapolated values ---
    passed = False
    if f_extrap is not None and kd_extrap is not None:
        delta = abs(f_extrap - kd_extrap)
        passed = delta < TOL_RICHARDSON_MEV
        print(f"\nExtrapolated CB1 delta: {delta:.3f} meV "
              f"(tolerance: {TOL_RICHARDSON_MEV} meV)")
        print(f"Overall: {'PASS' if passed else 'FAIL'}")
    else:
        print("\nERROR: insufficient data for Richardson extrapolation")

    results = {
        "fortran_resolution_sweep": [
            {"label": r["label"], "cb1_meV": float(v), "h_ang": float(h)}
            for r, v, h in zip(FORTRAN_RESOLUTIONS, f_cb1_list, f_h_list)
        ],
        "kdotpy_resolution_sweep": [
            {"label": r["label"], "cb1_meV": float(v), "zres_nm": float(h)}
            for r, v, h in zip(KDOTPY_RESOLUTIONS, kd_cb1_list, kd_h_list)
        ],
        "fortran_richardson_cb1": float(f_extrap) if f_extrap else None,
        "kdotpy_richardson_cb1": float(kd_extrap) if kd_extrap else None,
        "delta_meV": float(abs(f_extrap - kd_extrap)) if (f_extrap and kd_extrap) else None,
        "passed": bool(passed),
    }

    results_dir = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(results_dir, exist_ok=True)
    with open(os.path.join(results_dir, "qw_convergence.json"), "w") as f:
        json.dump(results, f, indent=2)

    return passed


if __name__ == "__main__":
    success = test_qw_convergence()
    sys.exit(0 if success else 1)
