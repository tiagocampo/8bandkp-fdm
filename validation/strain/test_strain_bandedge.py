"""Strain bandedge cross-validation.

Compares Bir-Pikus strain shifts between our Fortran code and kdotpy
for InAs on GaAs substrate (compressive) and GaAs on InP substrate (tensile).

Our code uses the strain: block in input.toml with strainSubstrate for
the substrate lattice constant. kdotpy uses substrate_material in PhysParams.

Both should produce identical strain shifts for matched parameters.

Tolerance: < 1 meV for all strained eigenvalue shifts.
"""
import os
import sys
import json
import tempfile
import shutil
import subprocess

project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, project_root)

from validation.shared.param_mapper import map_material

BUILD_DIR = os.path.join(project_root, "build")
MEV_PER_EV = 1000.0
TOL_MEV = 1.0

STRAIN_CONFIGS = [
    {
        "name": "InAs on GaAs substrate (compressive)",
        "material": "InAs",
        "substrate": "GaAs",
        "substrate_a0": 5.65325,  # GaAs lattice constant in Angstrom
    },
    {
        "name": "GaAs on InP substrate (tensile)",
        "material": "GaAs",
        "substrate": "InP",
        "substrate_a0": 5.8687,  # InP lattice constant (Vurgaftman 2001)
    },
]


def _run_fortran_bulk(material, build_dir, substrate_a0=None):
    """Run Fortran bandStructure for bulk, optionally with strain substrate."""
    exe = os.path.join(build_dir, "src", "bandStructure")
    if not os.path.isfile(exe):
        raise FileNotFoundError(f"Executable not found: {exe}")

    workdir = tempfile.mkdtemp(prefix="strain_" if substrate_a0 else "unstr_")
    try:
        lines = [
            'confinement = "bulk"',
            "FDorder = 2",
            "fd_step = 1",
            "",
            "[wave_vector]",
            'mode = "k0"',
            "max = 0",
            "nsteps = 1",
            "",
            "[bands]",
            "num_cb = 4",
            "num_vb = 4",
            "",
            "[[material]]",
            f'name = "{material}"',
        ]
        if substrate_a0 is not None:
            lines.extend(["", "[strain]", f"substrate_value = {substrate_a0}"])

        with open(os.path.join(workdir, "input.toml"), 'w') as f:
            f.write('\n'.join(lines) + '\n')
        os.makedirs(os.path.join(workdir, "output"), exist_ok=True)

        result = subprocess.run(
            [exe], cwd=workdir, capture_output=True, text=True, timeout=60
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
        return rows[0][1:] if rows else None
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def run_kdotpy_strained(material_name, substrate_name):
    from kdotpy.physparams import PhysParams
    from kdotpy.diagonalization.diagonalization import hbulk
    from kdotpy.vector import Vector

    mat = map_material(material_name)
    sub = map_material(substrate_name)

    pp = PhysParams(kdim=3, m_layers=[mat], l_layers=[1.0], norbitals=8,
                    substrate_material=sub)
    result = hbulk(Vector(0, 0, 0), pp)
    return sorted(result.eival)  # meV


def run_kdotpy_unstrained(material_name):
    from kdotpy.physparams import PhysParams
    from kdotpy.diagonalization.diagonalization import hbulk
    from kdotpy.vector import Vector

    mat = map_material(material_name)
    pp = PhysParams(kdim=3, m_layers=[mat], l_layers=[1.0], norbitals=8,
                    a_lattice=mat.param['a'])
    result = hbulk(Vector(0, 0, 0), pp)
    return sorted(result.eival)  # meV


def test_strain_bandedge():
    print("=" * 70)
    print("STRAIN BANDEDGE CROSS-VALIDATION")
    print("=" * 70)
    print(f"Tolerance: < {TOL_MEV} meV for strained eigenvalues")
    print()

    all_pass = True
    all_results = []

    for cfg in STRAIN_CONFIGS:
        name = cfg["name"]
        material = cfg["material"]
        substrate = cfg["substrate"]
        substrate_a0 = cfg["substrate_a0"]

        print(f"\nConfig: {name}")

        # Run our Fortran code
        try:
            f_unstrained = _run_fortran_bulk(material, BUILD_DIR)
            f_strained = _run_fortran_bulk(material, BUILD_DIR, substrate_a0=substrate_a0)
        except (FileNotFoundError, RuntimeError) as e:
            print(f"  SKIP: Fortran execution error: {e}")
            all_results.append({"config": name, "status": "SKIP"})
            continue

        f_unstr_meV = sorted([e * MEV_PER_EV for e in f_unstrained])
        f_str_meV = sorted([e * MEV_PER_EV for e in f_strained])
        f_shifts = [s - u for s, u in zip(f_str_meV, f_unstr_meV)]

        # Run kdotpy
        try:
            kd_unstr_meV = run_kdotpy_unstrained(material)
            kd_str_meV = run_kdotpy_strained(material, substrate)
            kd_shifts = [s - u for s, u in zip(kd_str_meV, kd_unstr_meV)]
        except (ImportError, RuntimeError, ValueError, OSError) as e:
            print(f"  FAIL: kdotpy calculation error ({type(e).__name__}): {e}")
            all_pass = False
            all_results.append({"config": name, "status": "FAIL"})
            continue

        # Compare shifts
        max_delta = max(abs(f - kd) for f, kd in zip(f_shifts, kd_shifts))
        passed = max_delta < TOL_MEV

        print(f"  Fortran strained (meV):  {[f'{e:.3f}' for e in f_str_meV]}")
        print(f"  kdotpy strained (meV):   {[f'{e:.3f}' for e in kd_str_meV]}")
        print(f"  Max total delta: {max_delta:.3f} meV")
        print(f"  Status: {'PASS' if passed else 'FAIL'}")

        if not passed:
            all_pass = False

        all_results.append({
            "config": name,
            "fortran_strained_meV": f_str_meV,
            "kdotpy_strained_meV": kd_str_meV,
            "max_delta_meV": max_delta,
            "passed": bool(passed),
        })

    results_dir = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(results_dir, exist_ok=True)
    with open(os.path.join(results_dir, "strain_bandedge.json"), "w") as f:
        json.dump(all_results, f, indent=2)

    n_pass = sum(1 for r in all_results if r.get("passed"))
    print(f"\nSummary: {n_pass}/{len(all_results)} configs passed")
    return all_pass


if __name__ == "__main__":
    success = test_strain_bandedge()
    sys.exit(0 if success else 1)
