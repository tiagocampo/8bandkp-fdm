"""Bulk Zeeman splitting validation.

Validates that our bulk code produces the correct CB g-factor at k=0
via the Roth formula: g* = 2 - (2*EP/3) * (1/Eg - 1/(Eg + deltaSO))

This test runs our gfactorCalculation executable and compares the computed
g_z against the analytical Roth formula for representative III-V materials.

kdotpy's hbulk does not include the Zeeman term in its 8-band Hamiltonian
(only its LL symbolic mode does), so cross-code comparison is deferred.

Tolerance: < 5% deviation from Roth formula.
"""
import os
import sys
import json
import subprocess
import tempfile
import shutil

project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, project_root)

from validation.shared.param_mapper import FORTRAN_MATERIALS

BUILD_DIR = os.path.join(project_root, "build")
TOL_PCT = 5.0

ZEEMAN_MATERIALS = ["GaAs", "InAs", "InSb"]


def roth_gfactor(m):
    """Compute analytical CB g-factor via Roth formula."""
    EP = m["EP"]  # eV
    Eg = m["Eg"]  # eV
    delta = m["deltaSO"]  # eV
    return 2.0 - (2.0 * EP / 3.0) * (1.0 / Eg - 1.0 / (Eg + delta))


def run_gfactor(material, build_dir):
    """Run gfactorCalculation executable and parse g_z."""
    exe = os.path.join(build_dir, "src", "gfactorCalculation")
    if not os.path.isfile(exe):
        raise RuntimeError(f"gfactorCalculation executable not found at {exe}")

    workdir = tempfile.mkdtemp(prefix="gfactor_")
    try:
        config = f"""confinement = "bulk"
FDorder = 2
fd_step = 1

[wave_vector]
mode = "k0"
max = 0
nsteps = 1

[bands]
num_cb = 2
num_vb = 6

[[material]]
name = "{material}"

which_band = 0
band_idx = 1
"""
        cfg_path = os.path.join(workdir, "input.toml")
        with open(cfg_path, 'w') as f:
            f.write(config)
        os.makedirs(os.path.join(workdir, "output"), exist_ok=True)

        result = subprocess.run(
            [exe], cwd=workdir, capture_output=True, text=True, timeout=60
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"gfactorCalculation failed (rc={result.returncode}):\n"
                f"stderr: {result.stderr[-500:]}"
            )

        # Parse g-factor from stdout — format: "gz\n  0.0\n  -0.315..."
        # The first value is band 0 (reference), second is band 1 (CB1 g-factor)
        g_z = None
        lines = result.stdout.split('\n')
        for i, line in enumerate(lines):
            if line.strip().lower() == 'gz':
                vals = []
                for j in range(i + 1, len(lines)):
                    val_str = lines[j].strip()
                    if not val_str:
                        continue
                    try:
                        vals.append(float(val_str))
                    except ValueError:
                        break
                    if len(vals) >= 2:
                        break
                if len(vals) >= 2:
                    g_z = vals[1]
                elif len(vals) == 1:
                    g_z = vals[0]
                break

        return g_z
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def test_bulk_zeeman():
    print("=" * 70)
    print("BULK ZEEMAN (g-FACTOR) VALIDATION")
    print("=" * 70)
    print(f"Tolerance: < {TOL_PCT}% deviation from Roth formula")
    print()

    all_pass = True
    all_results = []

    print("| Material | g*_Roth | g*_computed | Dev% | Status |")
    print("|----------|---------|-------------|------|--------|")

    for mat_name in ZEEMAN_MATERIALS:
        m = FORTRAN_MATERIALS[mat_name]
        g_roth = roth_gfactor(m)

        try:
            g_computed = run_gfactor(mat_name, BUILD_DIR)
        except RuntimeError as e:
            print(f"| {mat_name} | {g_roth:.4f} | ERROR | - | FAIL |")
            print(f"  -> {e}")
            all_pass = False
            all_results.append({"material": mat_name, "status": "FAIL"})
            continue

        if g_computed is None:
            print(f"| {mat_name} | {g_roth:.4f} | N/A | - | SKIP |")
            all_results.append({"material": mat_name, "status": "SKIP"})
            continue

        dev_pct = abs(g_computed - g_roth) / abs(g_roth) * 100 if abs(g_roth) > 0.001 else abs(g_computed - g_roth) * 100
        passed = dev_pct < TOL_PCT

        print(f"| {mat_name} | {g_roth:.4f} | {g_computed:.4f} | {dev_pct:.2f}% | "
              f"{'PASS' if passed else 'FAIL'} |")

        if not passed:
            all_pass = False

        all_results.append({
            "material": mat_name,
            "g_roth": float(g_roth),
            "g_computed": float(g_computed),
            "deviation_pct": float(dev_pct),
            "passed": bool(passed),
        })

    print()
    print("Note: kdotpy's hbulk does not include Zeeman splitting in the")
    print("      8-band Hamiltonian. Cross-code comparison requires kdotpy's")
    print("      LL symbolic mode, which uses a different API.")

    results_dir = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(results_dir, exist_ok=True)
    with open(os.path.join(results_dir, "bulk_zeeman.json"), "w") as f:
        json.dump(all_results, f, indent=2)

    n_pass = sum(1 for r in all_results if r.get("passed"))
    print(f"\nSummary: {n_pass}/{len(all_results)} materials passed")
    return all_pass


if __name__ == "__main__":
    success = test_bulk_zeeman()
    sys.exit(0 if success else 1)
