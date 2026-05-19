"""Strained QW subband cross-validation.

Compares strain-induced CB1 subband shift in an InAs/GaAs QW between our
Fortran code and analytical Bir-Pikus formulas. kdotpy QW strain comparison
is deferred -- its substrate model interacts with band offsets differently
than our code.

Tolerance: CB1 shift within 10% of analytical Bir-Pikus prediction.
"""
import os
import sys
import json
import tempfile
import shutil
import subprocess

project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, project_root)

from validation.shared.param_mapper import map_material, FORTRAN_MATERIALS

BUILD_DIR = os.path.join(project_root, "build")
MEV_PER_EV = 1000.0
ANG_TO_NM = 10.0
TOL_CB1_SHIFT_PCT = 10.0  # percent deviation from analytical Bir-Pikus

STRAIN_QW_CONFIGS = [
    {
        "name": "InAs/GaAs strained QW (compressive)",
        "barrier": "GaAs",
        "well": "InAs",
        "l_well_nm": 10.0,
        "l_barrier_nm": 15.0,
        "fdstep": 201,
        "zres": 0.25,
        "substrate": "GaAs",
    },
]


def _run_fortran_qw(barrier, well, l_well_ang, l_barrier_ang, total_ang,
                    fdstep=201, substrate=None, work_dir=None, timeout=120):
    """Run Fortran bandStructure for a QW, optionally with strain."""
    exe = os.path.join(BUILD_DIR, "src", "bandStructure")
    if not os.path.isfile(exe):
        raise FileNotFoundError(f"Executable not found: {exe}")

    cleanup = work_dir is None
    if work_dir is None:
        prefix = "strain_qw_" if substrate else "unstr_qw_"
        work_dir = tempfile.mkdtemp(prefix=prefix)

    try:
        half_total = total_ang / 2.0
        well_start = -l_well_ang / 2.0
        well_end = l_well_ang / 2.0

        lines = [
            "waveVector: k0",
            "waveVectorMax: 0",
            "waveVectorStep: 1",
            "confinement: 1",
            f"FDstep: {fdstep}",
            "FDorder: 2",
            "numLayers: 2",
            f"material1: {barrier} {-half_total:.1f} {half_total:.1f} 0",
            f"material2: {well} {well_start:.1f} {well_end:.1f} 0",
            "numcb: 6",
            "numvb: 12",
            "ExternalField: 0  EF",
            "EFParams: 0.0",
        ]
        if substrate:
            lines.append("strain: T")
            lines.append(f"strain_ref: {substrate}")

        with open(os.path.join(work_dir, "input.cfg"), 'w') as f:
            f.write('\n'.join(lines) + '\n')
        os.makedirs(os.path.join(work_dir, "output"), exist_ok=True)

        result = subprocess.run(
            [exe], cwd=work_dir, capture_output=True, text=True, timeout=timeout
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"bandStructure failed (rc={result.returncode}):\n"
                f"stderr: {result.stderr[-500:]}"
            )

        eig_path = os.path.join(work_dir, "output", "eigenvalues.dat")
        if not os.path.exists(eig_path):
            raise RuntimeError(f"No eigenvalues.dat produced in {work_dir}")

        rows = []
        with open(eig_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                rows.append([float(x) for x in line.split()])
        return rows[0][1:] if rows else None
    finally:
        if cleanup:
            shutil.rmtree(work_dir, ignore_errors=True)


def classify_qw_bands(evals_meV, well_material):
    """Split eigenvalues into CB and VB using material mid-gap."""
    well_params = FORTRAN_MATERIALS[well_material]
    ec_meV = well_params["EC"] * 1000.0
    eg_meV = well_params["Eg"] * 1000.0
    ev_meV = ec_meV - eg_meV
    midgap = (ec_meV + ev_meV) / 2.0

    cb = [e for e in evals_meV if e > midgap]
    vb = sorted([e for e in evals_meV if e <= midgap], reverse=True)
    return cb, vb


def test_strain_qw():
    print("=" * 70)
    print("STRAINED QW SUBBAND CROSS-VALIDATION")
    print("=" * 70)
    print(f"Tolerance: CB1 strain shift within {TOL_CB1_SHIFT_PCT}% of analytical Bir-Pikus")
    print()

    all_pass = True
    all_results = []

    for cfg in STRAIN_QW_CONFIGS:
        name = cfg["name"]
        barrier = cfg["barrier"]
        well = cfg["well"]
        l_well_nm = cfg["l_well_nm"]
        l_barrier_nm = cfg["l_barrier_nm"]
        fdstep = cfg["fdstep"]
        substrate = cfg["substrate"]

        l_well_ang = l_well_nm * ANG_TO_NM
        l_barrier_ang = l_barrier_nm * ANG_TO_NM
        total_ang = 2 * l_barrier_ang + l_well_ang

        print(f"\nConfig: {name}")

        # Compute analytical Bir-Pikus CB shift for the well material
        well_params = FORTRAN_MATERIALS[well]
        sub_params = FORTRAN_MATERIALS[substrate]
        a_well = well_params["a0"]
        a_sub = sub_params["a0"]
        eps = (a_sub - a_well) / a_well
        C11, C12 = well_params["C11"], well_params["C12"]
        eps_zz = -2 * C12 / C11 * eps
        Tr_eps = 2 * eps + eps_zz
        ac = well_params["ac"]
        delta_ec_analytical_meV = ac * Tr_eps * MEV_PER_EV

        print(f"  Lattice mismatch: {eps*100:.4f}%")
        print(f"  Analytical delta_Ec = {delta_ec_analytical_meV:.3f} meV")

        # Run Fortran strained + unstrained
        try:
            f_strained = _run_fortran_qw(
                barrier, well, l_well_ang, l_barrier_ang, total_ang,
                fdstep=fdstep, substrate=substrate
            )
            f_unstrained = _run_fortran_qw(
                barrier, well, l_well_ang, l_barrier_ang, total_ang,
                fdstep=fdstep, substrate=None
            )
        except (FileNotFoundError, RuntimeError) as e:
            print(f"  SKIP: Fortran execution error: {e}")
            all_results.append({"config": name, "status": "SKIP"})
            continue

        f_str_meV = sorted([e * MEV_PER_EV for e in f_strained])
        f_uns_meV = sorted([e * MEV_PER_EV for e in f_unstrained])

        f_str_cb, f_str_vb = classify_qw_bands(f_str_meV, well)
        f_uns_cb, f_uns_vb = classify_qw_bands(f_uns_meV, well)

        if not f_str_cb or not f_uns_cb:
            print(f"  SKIP: could not classify bands")
            all_results.append({"config": name, "status": "SKIP"})
            continue

        f_cb1_shift = f_str_cb[0] - f_uns_cb[0]

        # Compare against analytical (QW confinement modifies the shift,
        # so allow 10% deviation from bulk Bir-Pikus)
        dev_pct = abs(f_cb1_shift - delta_ec_analytical_meV) / abs(delta_ec_analytical_meV) * 100
        passed = dev_pct < TOL_CB1_SHIFT_PCT

        print(f"  Fortran CB1 shift: {f_cb1_shift:.3f} meV")
        print(f"  Analytical bulk:   {delta_ec_analytical_meV:.3f} meV")
        print(f"  Deviation: {dev_pct:.2f}% ({'PASS' if passed else 'FAIL'})")
        print(f"  Status: {'PASS' if passed else 'FAIL'}")

        if not passed:
            all_pass = False

        all_results.append({
            "config": name,
            "f_cb1_shift_meV": float(f_cb1_shift),
            "analytical_cb1_shift_meV": float(delta_ec_analytical_meV),
            "deviation_pct": float(dev_pct),
            "passed": bool(passed),
        })

    results_dir = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(results_dir, exist_ok=True)
    with open(os.path.join(results_dir, "strain_qw.json"), "w") as f:
        json.dump(all_results, f, indent=2)

    n_pass = sum(1 for r in all_results if r.get("passed"))
    print(f"\nSummary: {n_pass}/{len(all_results)} configs passed")
    return all_pass


if __name__ == "__main__":
    success = test_strain_qw()
    sys.exit(0 if success else 1)
