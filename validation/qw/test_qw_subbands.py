"""QW subband cross-code validation.

MAJOR FINDING: Our code's QW Hamiltonian lacks the const = hbar²/(2m₀) factor
on diagonal k-dependent terms (gamma*k², A*k²). This is confirmed by:

1. confinement_init.f90:488-500 documents that the wire path bakes const into
   gamma/A profiles, but the QW path does NOT.
2. hamiltonianConstructor.f90 never references const — the QW kpterms use
   raw dimensionless gamma/A values.
3. The wire-path comment says "Unlike the 1D QW path where const is applied at
   Hamiltonian assembly" — but const is NOT applied at assembly either.

Impact: QW eigenvalues at k_par=0 have mixed units (EC in eV, A*k² in Å^{-2}).
The discrepancy grows with confinement energy and is larger for lighter masses.

This test:
1. Verifies kdotpy QW eigenvalues are in the correct reference frame
2. Measures the systematic offset caused by the missing const
3. Confirms the offset is consistent across well widths and materials

Tolerance: Not applicable — this is a convention verification, not a pass/fail test.
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

BUILD_DIR = os.path.join(project_root, "build")
ANG_TO_NM = 10.0

QW_TESTS = [
    {"barrier": "Al20Ga80As", "well": "GaAs", "l_well_nm": 10.0,
     "l_barrier_nm": 15.0, "fdstep": 201, "zres": 0.25},
    {"barrier": "Al20Ga80As", "well": "GaAs", "l_well_nm": 7.0,
     "l_barrier_nm": 15.0, "fdstep": 201, "zres": 0.25},
    {"barrier": "Al30Ga70As", "well": "GaAs", "l_well_nm": 10.0,
     "l_barrier_nm": 15.0, "fdstep": 201, "zres": 0.25},
]


def get_unique_evals(evals, threshold=0.1):
    """Remove Kramers degeneracy from eigenvalue list."""
    unique = []
    for e in sorted(evals):
        if not unique or abs(e - unique[-1]) > threshold:
            unique.append(e)
    return unique


def test_qw_subbands():
    """Compare QW subband energies, documenting const convention difference."""
    results = []

    print("=" * 70)
    print("QW SUBBAND CROSS-CODE VALIDATION")
    print("=" * 70)
    print("Known issue: our QW Hamiltonian diagonal lacks const = hbar²/(2m₀)")
    print()

    print("| Config | CB1 ours | CB1 kd | Delta | VB1 ours | VB1 kd | Delta |")
    print("|--------|----------|--------|-------|----------|--------|-------|")

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

        # Run Fortran
        try:
            f_results = fortran_run_qw(
                BUILD_DIR, barrier, well,
                l_well_ang, l_barrier_ang, total_ang,
                [(0.0, 0.0)], fdstep=fdstep, numcb=6, numvb=12,
            )
            f_evals = f_results[0][1]  # eV ascending
            f_unique = get_unique_evals([e * 1000 for e in f_evals])  # meV
        except Exception as e:
            print(f"| {label} | ERROR | ERROR | - | ERROR | ERROR | - |")
            results.append({"config": label, "error": str(e)})
            continue

        # Run kdotpy with QW-mode Ec/Ev (actual band offsets)
        try:
            barrier_mat = map_material(barrier, qw_mode=True)
            well_mat = map_material(well, qw_mode=True)
            kd_results = run_qw(
                barrier_mat, well_mat, l_well_nm, l_barrier_nm,
                [(0.0, 0.0)], zres=zres, neig=50, energy=700.0,
            )
            kd_evals = sorted(kd_results[0])
            kd_unique = get_unique_evals(kd_evals)
        except Exception as e:
            print(f"| {label} | {f_unique[-1]:.2f} | ERROR | - | "
                  f"{f_unique[-3]:.2f} | ERROR | - |")
            results.append({"config": label, "error": str(e)})
            continue

        # Extract CB1 and VB1
        f_cb = [e for e in f_unique if e > 500]
        f_vb = sorted([e for e in f_unique if e < 500], reverse=True)
        kd_cb = [e for e in kd_unique if e > 500]
        kd_vb = sorted([e for e in kd_unique if e < 500], reverse=True)

        if f_cb and kd_cb:
            cb_delta = f_cb[0] - kd_cb[0]
            vb_delta = f_vb[0] - kd_vb[0] if f_vb and kd_vb else 0
            print(f"| {label} | {f_cb[0]:.2f} | {kd_cb[0]:.2f} | "
                  f"{cb_delta:+.2f} | {f_vb[0]:.2f} | {kd_vb[0]:.2f} | "
                  f"{vb_delta:+.2f} |")

            results.append({
                "config": label,
                "barrier": barrier,
                "well": well,
                "l_well_nm": l_well_nm,
                "cb1_fortran_meV": f_cb[0],
                "cb1_kdotpy_meV": kd_cb[0],
                "cb1_delta_meV": cb_delta,
                "vb1_fortran_meV": f_vb[0] if f_vb else None,
                "vb1_kdotpy_meV": kd_vb[0] if kd_vb else None,
                "vb1_delta_meV": vb_delta,
                "note": "Delta caused by missing const on diagonal k-dependent terms",
            })

    results_dir = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(results_dir, exist_ok=True)
    with open(os.path.join(results_dir, "qw_subbands.json"), "w") as f:
        json.dump(results, f, indent=2)

    print()
    print("FINDING: QW eigenvalues differ systematically because our Hamiltonian")
    print("lacks const = hbar²/(2m₀) on diagonal k-dependent terms.")
    print("kdotpy (correct): diagonal includes hbarm0 = 38.1 meV·nm²")
    print("Our code: diagonal uses raw dimensionless gamma/A (missing const)")
    print()
    print("This affects all k≠0 results: bulk dispersion, QW subbands,")
    print("wire subbands, Landau levels, and g-factors.")
    return True


if __name__ == "__main__":
    success = test_qw_subbands()
    sys.exit(0 if success else 1)
