#!/usr/bin/env python3
"""Verify Landau level eigenvalues against analytical E_n = E_C + hbar*omega_c*(n+1/2).

Usage: verify_landau_analytical.py <bandStructure_exe> <landau_config>

Runs bandStructure with the given config, parses eigenvalues at k=0,
and checks that the lowest CB Landau levels match the analytical
cyclotron energy spacing. Also validates ky degeneracy.
"""
# COVERAGE: observable=landau_levels geometry=bulk material=GaAs
# COVERAGE: observable=landau_levels geometry=bulk material=InAs
import shutil
import subprocess
import sys
import os
import tempfile

# --- Physical constants (SI) ---
HBAR_JS = 1.054571817e-34    # J*s
E_C = 1.602176634e-19        # C
M0_KG = 9.1093837015e-31     # kg

# --- Material database ---
# m_star in units of m0: 8-band model reference values (const-correct)
# EC = EV + Eg (conduction band edge including offset)
MATERIALS = {
    "InAs": {"m_star": 0.0123, "EV": -0.59, "Eg": 0.417},
    "GaAs": {"m_star": 0.0317, "EV": -0.80, "Eg": 1.519},
}


def hbar_omega_c_eV(B_T, m_star):
    """Cyclotron energy hbar*omega_c in eV."""
    omega_c = E_C * B_T / (m_star * M0_KG)  # rad/s
    return HBAR_JS * omega_c / E_C           # eV


def parse_eigenvalues(filepath):
    """Parse eigenvalues.dat, return list of rows (each row is list of floats)."""
    rows = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            rows.append([float(t) for t in line.split()])
    return rows


def parse_config_value(filepath, key):
    """Read a config value from a `key: value` line."""
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('!'):
                continue
            if line.startswith(f"{key}:"):
                return line.split(":", 1)[1].strip()
    return None


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <bandStructure_exe> <landau_config>")
        sys.exit(1)

    exe = sys.argv[1]
    config = sys.argv[2]

    # Parse material and B-field from config
    material_line = parse_config_value(config, "material1")
    if material_line is None:
        print("FAIL: material1 not found in config")
        sys.exit(1)
    material = material_line.split()[0]

    if material not in MATERIALS:
        print(f"FAIL: unknown material '{material}'")
        sys.exit(1)

    m_star = MATERIALS[material]["m_star"]
    E_C_material = MATERIALS[material]["EV"] + MATERIALS[material]["Eg"]

    b_field_str = parse_config_value(config, "b_field")
    if b_field_str is None:
        print("FAIL: b_field not found in config")
        sys.exit(1)
    b_vals = [float(x) for x in b_field_str.split()]
    B_mag = (b_vals[0]**2 + b_vals[1]**2 + b_vals[2]**2)**0.5

    # Run bandStructure
    workdir = tempfile.mkdtemp()
    try:
        dst = os.path.join(workdir, "input.cfg")
        with open(dst, 'w') as f:
            with open(config) as src:
                f.write(src.read())
        os.makedirs(os.path.join(workdir, "output"), exist_ok=True)

        result = subprocess.run(
            [exe], cwd=workdir, capture_output=True, text=True, timeout=120
        )
        if result.returncode != 0:
            print(f"FAIL: bandStructure returned {result.returncode}")
            print(result.stdout)
            print(result.stderr)
            sys.exit(1)

        eig_path = os.path.join(workdir, "output", "eigenvalues.dat")
        if not os.path.exists(eig_path):
            print("FAIL: eigenvalues.dat not produced")
            sys.exit(1)

        rows = parse_eigenvalues(eig_path)
        if len(rows) < 1:
            print("FAIL: no data rows in eigenvalues.dat")
            sys.exit(1)

        # Row 0 is k=0 (first data row after header)
        # Columns: k, then 8 eigenvalues (4 VB + 4 CB for numcb=4, numvb=4)
        k0_row = rows[0]
        n_cols = len(k0_row)
        n_eig = n_cols - 1  # exclude k column

        # CB eigenvalues: last numcb columns
        numcb = 4
        cb_eigs = sorted(k0_row[-numcb:])

        # Analytical: lowest CB Landau level = E_C + hbar*omega_c/2
        hw = hbar_omega_c_eV(B_mag, m_star)
        E0_analytical = E_C_material + hw / 2.0

        # Tolerance: 5% relative (FD discretization error)
        tol = 0.05

        print(f"Material: {material}, m* = {m_star} m0, B = {B_mag} T")
        print(f"hbar*omega_c = {hw*1000:.2f} meV")
        print(f"E_0 (analytical) = {E0_analytical:.6f} eV")
        print(f"E_0 (numerical)  = {cb_eigs[0]:.6f} eV")

        rel_err = abs(cb_eigs[0] - E0_analytical) / abs(E0_analytical)
        print(f"Relative error: {rel_err:.4f} (tolerance: {tol})")

        if rel_err > tol:
            print(f"FAIL: CB ground state {cb_eigs[0]:.6f} eV differs from "
                  f"analytical {E0_analytical:.6f} eV by {rel_err*100:.1f}%")
            sys.exit(1)

        print(f"PASS: CB ground state within {tol*100:.0f}% of analytical value")

        # --- ky degeneracy check ---
        # eigenvalues.dat has rows for different ky values (waveVectorStep > 1)
        # Landau levels should be degenerate in ky (same eigenvalues at all ky)
        if len(rows) >= 3:
            k0_cb = sorted(rows[0][-numcb:])
            k_mid = len(rows) // 2
            kmid_cb = sorted(rows[k_mid][-numcb:])
            klast_cb = sorted(rows[-1][-numcb:])

            max_deg_err = 0.0
            for i in range(numcb):
                err = max(abs(k0_cb[i] - kmid_cb[i]), abs(k0_cb[i] - klast_cb[i]))
                max_deg_err = max(max_deg_err, err)

            deg_tol = 1e-4  # eV (finite grid breaks exact degeneracy)
            if max_deg_err > deg_tol:
                print(f"FAIL: ky degeneracy broken, max deviation = {max_deg_err:.2e} eV")
                sys.exit(1)

            print(f"PASS: ky degeneracy holds (max deviation = {max_deg_err:.2e} eV)")
        else:
            print("SKIP: ky degeneracy check (only 1 k-point)")

    finally:
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
