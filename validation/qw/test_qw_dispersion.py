"""QW in-plane dispersion cross-code validation.

Compares QW subband dispersion E(k_par) between our Fortran code and kdotpy
for a GaAs/AlGaAs QW. Extracts CB1 effective mass from parabolic curvature
and verifies agreement.

Tolerance: < 2% for effective mass, < 2 meV per k-point CB1 eigenvalue.
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
from validation.qw.test_qw_subbands import classify_bands, get_unique_evals

BUILD_DIR = os.path.join(project_root, "build")
ANG_TO_NM = 10.0
MEV_PER_EV = 1000.0

TOL_MASS_PCT = 0.02   # 2% for effective mass
TOL_CB1_MEV = 2.0      # 2 meV per k-point CB1 eigenvalue (FD vs plane-wave)

QW_CONFIGS = [
    {"barrier": "Al20Ga80As", "well": "GaAs", "l_well_nm": 10.0,
     "l_barrier_nm": 15.0, "fdstep": 201, "zres": 0.25,
     "k_max_inv_nm": 0.2, "n_k": 21},
]


def extract_cb1_dispersion(evals_all, well_material, threshold=0.5):
    """Extract CB1 energy at each k-point from eigenvalue lists.

    Returns list of CB1 energies in meV (one per k-point).
    """
    cb1_list = []
    for evals in evals_all:
        unique = get_unique_evals(sorted(evals), threshold)
        cb, vb = classify_bands(unique, well_material)
        if cb:
            cb1_list.append(cb[0])
        else:
            cb1_list.append(None)
    return cb1_list


def parabolic_mass(k_inv_nm, e_meV):
    """Extract effective mass from parabolic fit to E(k) near k=0.

    Fits E(k) = E0 + (hbar^2 / 2m0) * k^2 / m* using first N points.
    Returns m* in units of m0.
    """
    # hbar^2/(2m0) = 3.80998 eV*A^2 = 38.0998 meV*nm^2
    HBAR2_OVER_2M0_MEV_NM2 = 38.0998

    n_fit = min(5, len(k_inv_nm))
    k_fit = np.array(k_inv_nm[:n_fit])
    e_fit = np.array(e_meV[:n_fit])

    k2 = k_fit**2
    coeffs = np.polyfit(k2, e_fit, 1)
    c2 = coeffs[0]  # meV * nm^2

    if abs(c2) < 1e-15:
        return None

    return HBAR2_OVER_2M0_MEV_NM2 / c2


def test_qw_dispersion():
    """Compare QW dispersion between codes."""
    all_pass = True
    all_results = []

    print("=" * 70)
    print("QW DISPERSION CROSS-CODE VALIDATION")
    print("=" * 70)
    print(f"Tolerance: CB1 per k-point < {TOL_CB1_MEV} meV, "
          f"effective mass < {TOL_MASS_PCT*100:.0f}%")
    print()

    # Check kdotpy availability
    try:
        from kdotpy.config import initialize_config
        initialize_config()
    except ImportError:
        print("SKIP: kdotpy not available (activate kdotpy_env first)")
        return False

    for cfg in QW_CONFIGS:
        barrier = cfg["barrier"]
        well = cfg["well"]
        l_well_nm = cfg["l_well_nm"]
        l_barrier_nm = cfg["l_barrier_nm"]
        fdstep = cfg["fdstep"]
        zres = cfg["zres"]
        k_max = cfg["k_max_inv_nm"]
        n_k = cfg["n_k"]

        l_well_ang = l_well_nm * ANG_TO_NM
        l_barrier_ang = l_barrier_nm * ANG_TO_NM
        total_ang = 2 * l_barrier_ang + l_well_ang

        label = f"{well}/{barrier} {l_well_nm}nm"
        print(f"\nConfig: {label}")

        k_points_nm = [(k_max * i / (n_k - 1), 0.0) for i in range(n_k)]
        k_mags_nm = [kx for kx, ky in k_points_nm]

        # Fortran: k in 1/Angstrom
        k_points_ang = [(kx / ANG_TO_NM, ky / ANG_TO_NM) for kx, ky in k_points_nm]

        # Run Fortran
        try:
            f_results = fortran_run_qw(
                BUILD_DIR, barrier, well,
                l_well_ang, l_barrier_ang, total_ang,
                k_points_ang, fdstep=fdstep, numcb=6, numvb=12,
            )
            f_cb1 = extract_cb1_dispersion(
                [[e * MEV_PER_EV for e in row[1]] for row in f_results],
                well,
            )
        except (FileNotFoundError, RuntimeError, OSError) as e:
            print(f"  Fortran ERROR ({type(e).__name__}): {e}")
            all_pass = False
            continue

        # Run kdotpy
        try:
            barrier_mat = map_material(barrier, qw_mode=True)
            well_mat = map_material(well, qw_mode=True)
            kd_results = run_qw(
                barrier_mat, well_mat, l_well_nm, l_barrier_nm,
                k_points_nm, zres=zres, neig=50, energy=700.0,
            )
            kd_cb1 = extract_cb1_dispersion(kd_results, well)
        except (ImportError, RuntimeError, ValueError, OSError) as e:
            print(f"  kdotpy ERROR ({type(e).__name__}): {e}")
            all_pass = False
            continue

        # Compare CB1 at each k-point
        if None in f_cb1 or None in kd_cb1:
            print(f"  ERROR: could not extract CB1 from eigenvalues")
            all_pass = False
            continue

        max_delta = 0.0
        n_pass_k = 0
        for i, (f_val, kd_val) in enumerate(zip(f_cb1, kd_cb1)):
            delta = abs(f_val - kd_val)
            max_delta = max(max_delta, delta)
            if delta < TOL_CB1_MEV:
                n_pass_k += 1

        # Extract effective masses
        f_mass = parabolic_mass(k_mags_nm, f_cb1)
        kd_mass = parabolic_mass(k_mags_nm, kd_cb1)

        mass_pass = False
        if f_mass and kd_mass and f_mass > 0 and kd_mass > 0:
            mass_dev = abs(f_mass - kd_mass) / kd_mass
            mass_pass = mass_dev < TOL_MASS_PCT
            print(f"  m* Fortran: {f_mass:.4f} m0")
            print(f"  m* kdotpy:  {kd_mass:.4f} m0")
            print(f"  Mass deviation: {mass_dev*100:.2f}% "
                  f"({'PASS' if mass_pass else 'FAIL'})")
        else:
            print(f"  WARNING: could not extract effective mass")
            all_pass = False

        print(f"  CB1 max delta: {max_delta:.2f} meV")
        print(f"  k-points: {n_pass_k}/{n_k} within {TOL_CB1_MEV} meV")

        cb1_pass = n_pass_k == n_k
        config_pass = cb1_pass and mass_pass

        if not config_pass:
            all_pass = False

        print(f"  Overall: {'PASS' if config_pass else 'FAIL'}")

        all_results.append({
            "config": label,
            "f_mass": float(f_mass) if f_mass else None,
            "kd_mass": float(kd_mass) if kd_mass else None,
            "mass_deviation_pct": float(abs(f_mass - kd_mass) / kd_mass * 100) if (f_mass and kd_mass and kd_mass > 0) else None,
            "mass_pass": bool(mass_pass),
            "cb1_max_delta_meV": float(max_delta),
            "cb1_pass": bool(cb1_pass),
            "passed": bool(config_pass),
        })

    results_dir = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(results_dir, exist_ok=True)
    with open(os.path.join(results_dir, "qw_dispersion.json"), "w") as f:
        json.dump(all_results, f, indent=2)

    n_pass = sum(1 for r in all_results if r.get("passed"))
    print(f"\nSummary: {n_pass}/{len(all_results)} configs passed")
    return all_pass


if __name__ == "__main__":
    success = test_qw_dispersion()
    sys.exit(0 if success else 1)
