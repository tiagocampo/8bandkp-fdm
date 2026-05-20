#!/usr/bin/env python3
"""Exciton smoke test (R11a).

Single-resolution test verifying the exciton module runs correctly and
produces a physically plausible binding energy (1-20 meV for GaAs/AlGaAs QW).

COVERAGE: observable=exciton_Eb geometry=QW material=GaAs/AlGaAs tier=smoke

Usage:
    python3 test_exciton_smoke.py <build_dir> <source_dir>
"""

import os
import re
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from star_helpers import run_exe

# Binding energy bounds for GaAs/AlGaAs QW (meV)
EB_MIN_MEV = 1.0
EB_MAX_MEV = 20.0


def parse_exciton_stdout(stdout_text):
    """Parse exciton binding energy from bandStructure stdout."""
    result = {}

    m_eb = re.search(r'E_binding\s*=\s*([\d.eE+-]+)', stdout_text)
    if m_eb:
        result['E_binding_meV'] = float(m_eb.group(1))

    m_lam = re.search(r'lambda_opt\s*=\s*([\d.eE+-]+)', stdout_text)
    if m_lam:
        result['lambda_opt_AA'] = float(m_lam.group(1))

    if 'E_binding_meV' not in result:
        m_alt = re.search(r'Exciton binding energy:\s*([\d.eE+-]+)', stdout_text)
        if m_alt:
            result['E_binding_meV'] = float(m_alt.group(1))

    if 'lambda_opt_AA' not in result:
        m_alt2 = re.search(r'Variational parameter:\s*([\d.eE+-]+)', stdout_text)
        if m_alt2:
            result['lambda_opt_AA'] = float(m_alt2.group(1))

    return result if result else None


def parse_exciton_file(filepath):
    """Parse output/exciton.dat file."""
    import numpy as np
    try:
        data = np.loadtxt(filepath, comments='#')
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.shape[1] >= 2:
            return {
                'lambda_opt_AA': float(data[0, 0]),
                'E_binding_meV': float(data[0, 1]),
                'mu_over_m0': float(data[0, 2]) if data.shape[1] > 2 else None,
                'eps_r': float(data[0, 3]) if data.shape[1] > 3 else None,
            }
    except Exception:
        pass
    return None


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        sys.exit(1)

    build_dir = sys.argv[1]
    source_dir = sys.argv[2]
    config_path = os.path.join(
        source_dir, "tests", "regression", "configs", "exciton_gaas_algaas.cfg"
    )

    if not os.path.isfile(config_path):
        print(f"FAIL: Config not found: {config_path}")
        sys.exit(1)

    print("=" * 60)
    print("  EXCITON SMOKE TEST (R11a)")
    print("  GaAs/AlGaAs QW, FDstep=101, FDorder=2")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as work_dir:
        rc, output_dir = run_exe(build_dir, "bandStructure", config_path, work_dir, timeout=120)
        if rc != 0:
            print(f"FAIL: bandStructure returned {rc}")
            sys.exit(1)

        # Check exciton output file
        exciton_path = os.path.join(output_dir, "exciton.dat")
        file_result = None
        if os.path.isfile(exciton_path):
            file_result = parse_exciton_file(exciton_path)
            print(f"  exciton.dat parsed: {file_result}")

        # Also parse stdout from the run
        stdout_path = os.path.join(work_dir, "stdout.txt") if os.path.exists(
            os.path.join(work_dir, "stdout.txt")) else None
        stdout_result = None

        # Re-run to capture stdout (or read from log)
        # Actually, run_exe uses capture_output=True, so we need a different approach
        # Let's just use the file output which is more reliable

        if file_result is None:
            print("FAIL: Could not parse exciton output")
            sys.exit(1)

        eb = file_result['E_binding_meV']
        print(f"\n  Binding energy: {eb:.4f} meV")
        print(f"  Expected range: [{EB_MIN_MEV}, {EB_MAX_MEV}] meV")

        # Assert binding energy in physical range
        if eb < EB_MIN_MEV or eb > EB_MAX_MEV:
            print(f"  FAIL: Binding energy {eb:.4f} meV outside [{EB_MIN_MEV}, {EB_MAX_MEV}]")
            sys.exit(1)
        print(f"  PASS: Binding energy in physical range")

        # Check variational parameters
        lam = file_result.get('lambda_opt_AA')
        if lam is not None:
            print(f"  lambda_opt: {lam:.4f} A")
            if lam <= 0:
                print("  FAIL: lambda_opt must be positive")
                sys.exit(1)
            print("  PASS: lambda_opt is positive")

        mu = file_result.get('mu_over_m0')
        if mu is not None:
            print(f"  mu/m0: {mu:.6f}")
            if mu <= 0 or mu > 1.0:
                print(f"  WARN: mu/m0 = {mu:.6f} outside typical (0, 1)")

        eps_r = file_result.get('eps_r')
        if eps_r is not None:
            print(f"  eps_r: {eps_r:.4f}")
            if eps_r < 1.0:
                print(f"  FAIL: eps_r = {eps_r:.4f} < 1 (unphysical)")
                sys.exit(1)
            print("  PASS: eps_r >= 1")

        print("\n  All checks passed.")
        sys.exit(0)


if __name__ == "__main__":
    main()
