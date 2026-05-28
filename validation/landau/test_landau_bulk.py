"""Bulk Landau level validation.

Compares our Fortran Landau level eigenvalues against analytical formulas.

Our code uses confinement=3 (Landau gauge discretization) with landau_nx
and landau_width parameters. Analytical reference: E_n = E_C + hbar*omega_c*(n+1/2).

Note: kdotpy's bulk-ll mode requires SymbolicHamiltonian setup that differs
fundamentally from our FD approach. This test validates against analytical
formulas; kdotpy LL comparison is deferred to follow-up work.

Tolerance: < 5% for CB ground state LL energy.
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
TOL_REL = 0.05  # 5% for LL energy

# Physical constants (SI)
HBAR_JS = 1.054571817e-34
EC_C = 1.602176634e-19
M0_KG = 9.1093837015e-31

# Material database: m_star from 8-band model, EC = EV + Eg
MATERIALS = {
    "GaAs": {"m_star": 0.0317, "EV": -0.80, "Eg": 1.519, "config": "landau_bulk_GaAs.toml"},
    "InAs": {"m_star": 0.0123, "EV": -0.59, "Eg": 0.417, "config": "landau_bulk_InAs.toml"},
}


def hbar_omega_c_eV(B_T, m_star):
    omega_c = EC_C * B_T / (m_star * M0_KG)
    return HBAR_JS * omega_c / EC_C


def run_fortran_landau(config_name, build_dir, project_root):
    exe = os.path.join(build_dir, "src", "bandStructure")
    config_path = os.path.join(project_root, "tests", "regression", "configs", config_name)

    if not os.path.isfile(exe):
        raise FileNotFoundError(f"Executable not found: {exe}")
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Config not found: {config_path}")

    workdir = tempfile.mkdtemp(prefix="landau_")
    try:
        shutil.copy2(config_path, os.path.join(workdir, "input.toml"))
        os.makedirs(os.path.join(workdir, "output"), exist_ok=True)
        result = subprocess.run(
            [exe], cwd=workdir, capture_output=True, text=True, timeout=120
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"bandStructure failed (rc={result.returncode}):\n"
                f"stderr: {result.stderr[-500:]}"
            )

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
        return rows
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def extract_cb_ll(rows, numcb=4):
    if not rows:
        return None
    k0_row = rows[0]
    cb_eigs = sorted(k0_row[-numcb:])
    return cb_eigs[0]  # ground state


def test_landau_bulk():
    print("=" * 70)
    print("BULK LANDAU LEVEL CROSS-VALIDATION")
    print("=" * 70)
    print(f"Tolerance: CB LL0 relative error < {TOL_REL*100:.0f}%")
    print()

    # This test validates against analytical formulas, not kdotpy.
    # No kdotpy availability check needed.

    all_pass = True
    all_results = []

    for mat_name, mat_data in MATERIALS.items():
        m_star = mat_data["m_star"]
        EC_material = mat_data["EV"] + mat_data["Eg"]
        config = mat_data["config"]

        print(f"\nMaterial: {mat_name} (m* = {m_star} m0)")

        try:
            rows = run_fortran_landau(config, BUILD_DIR, project_root)
        except (RuntimeError, FileNotFoundError) as e:
            print(f"  FAIL: {e}")
            all_pass = False
            all_results.append({"material": mat_name, "status": "FAIL"})
            continue

        # Parse B-field from config
        config_path = os.path.join(project_root, "tests", "regression", "configs", config)
        B_mag = None
        with open(config_path) as f:
            for line in f:
                if line.strip().startswith("b_field:"):
                    b_vals = [float(x) for x in line.split(":")[1].split()]
                    B_mag = (b_vals[0]**2 + b_vals[1]**2 + b_vals[2]**2)**0.5
                    break

        if B_mag is None:
            print(f"  SKIP: B-field not found in config")
            all_results.append({"material": mat_name, "status": "SKIP"})
            continue

        # Analytical LL energy
        hw = hbar_omega_c_eV(B_mag, m_star)
        E0_analytical = EC_material + hw / 2.0

        # Extract numerical CB ground state
        cb0 = extract_cb_ll(rows)
        if cb0 is None:
            print(f"  SKIP: could not extract CB eigenvalues")
            all_results.append({"material": mat_name, "status": "SKIP"})
            continue

        rel_err = abs(cb0 - E0_analytical) / abs(E0_analytical)
        passed = rel_err < TOL_REL

        print(f"  B = {B_mag:.1f} T")
        print(f"  hbar*omega_c = {hw*MEV_PER_EV:.2f} meV")
        print(f"  LL0 analytical: {E0_analytical:.6f} eV")
        print(f"  LL0 numerical:  {cb0:.6f} eV")
        print(f"  Relative error: {rel_err:.4f} ({rel_err*100:.1f}%)")
        print(f"  Status: {'PASS' if passed else 'FAIL'}")

        if not passed:
            all_pass = False

        all_results.append({
            "material": mat_name,
            "B_T": B_mag,
            "hbar_omega_c_meV": hw * MEV_PER_EV,
            "LL0_analytical_eV": E0_analytical,
            "LL0_numerical_eV": cb0,
            "rel_error": rel_err,
            "passed": passed,
        })

        # ky degeneracy check
        if len(rows) >= 3:
            numcb = 4
            k0_cb = sorted(rows[0][-numcb:])
            kmid_cb = sorted(rows[len(rows)//2][-numcb:])
            max_deg_err = max(abs(a - b) for a, b in zip(k0_cb, kmid_cb))
            deg_tol = 1e-4  # eV
            deg_pass = max_deg_err < deg_tol
            print(f"  ky degeneracy: max dev = {max_deg_err:.2e} eV "
                  f"({'PASS' if deg_pass else 'FAIL'})")
            if not deg_pass:
                all_pass = False

    results_dir = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(results_dir, exist_ok=True)
    with open(os.path.join(results_dir, "landau_bulk.json"), "w") as f:
        json.dump(all_results, f, indent=2)

    n_pass = sum(1 for r in all_results if r.get("passed"))
    n_total = len(all_results)
    print(f"\nSummary: {n_pass}/{n_total} materials passed")
    return all_pass


if __name__ == "__main__":
    success = test_landau_bulk()
    sys.exit(0 if success else 1)
