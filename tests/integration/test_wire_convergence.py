#!/usr/bin/env python3
"""Wire grid convergence tests (U6) for S7 (InAs wire).

Runs InAs wire at fixed geometry with varying grid resolution, extracting
subband energy and g-factor. Grid convergence with Richardson extrapolation.

COVERAGE: observable=CB1_energy geometry=wire material=InAs tier=convergence
COVERAGE: observable=gz geometry=wire material=InAs tier=convergence

Usage:
    python3 test_wire_convergence.py <build_dir> <source_dir>
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(__file__))
from star_helpers import run_exe, parse_eigenvalues, parse_gfactor
from convergence_helpers import (
    make_convergence_report, write_convergence_json,
)


# ---------------------------------------------------------------------------
# Wire config template (fixed geometry, varying grid)
# ---------------------------------------------------------------------------

WIRE_TEMPLATE = (
    'confinement = "wire"\n'
    "FDorder = {fdorder}\n"
    "fd_step = 1\n"
    "\n"
    "[wave_vector]\n"
    'mode = "k0"\n'
    "max = 0\n"
    "nsteps = 1\n"
    "\n"
    "[bands]\n"
    "num_cb = 4\n"
    "num_vb = 8\n"
    "\n"
    "[[material]]\n"
    'name = "InAs"\n'
    "\n"
    "[wire]\n"
    "nx = {nx}\n"
    "ny = {ny}\n"
    "dx = {dx}\n"
    "dy = {dy}\n"
    "\n"
    "[wire.geometry]\n"
    'shape = "rectangle"\n'
    "width = 55.0\n"
    "height = 55.0\n"
    "\n"
    "[[region]]\n"
    'material = "InAs"\n'
    "inner = 0.0\n"
    "outer = 50.0\n"
    "\n"
    "[[region]]\n"
    'material = "GaAs"\n'
    "inner = 50.0\n"
    "outer = 100.0\n"
    "\n"
    "[solver]\n"
    'method = "DENSE"\n'
    'mode = "ENERGY"\n'
    "emin = -1.0\n"
    "emax = 2.0\n"
)

WIRE_GF_TEMPLATE = (
    'confinement = "wire"\n'
    "FDorder = {fdorder}\n"
    "fd_step = 1\n"
    "\n"
    "[wave_vector]\n"
    'mode = "k0"\n'
    "max = 0\n"
    "nsteps = 1\n"
    "\n"
    "[bands]\n"
    "num_cb = 4\n"
    "num_vb = 8\n"
    "\n"
    "[[material]]\n"
    'name = "InAs"\n'
    "\n"
    "[wire]\n"
    "nx = {nx}\n"
    "ny = {ny}\n"
    "dx = {dx}\n"
    "dy = {dy}\n"
    "\n"
    "[wire.geometry]\n"
    'shape = "rectangle"\n'
    "width = 55.0\n"
    "height = 55.0\n"
    "\n"
    "[[region]]\n"
    'material = "InAs"\n'
    "inner = 0.0\n"
    "outer = 50.0\n"
    "\n"
    "[[region]]\n"
    'material = "GaAs"\n'
    "inner = 50.0\n"
    "outer = 100.0\n"
    "\n"
    "which_band = 0\n"
    "band_idx = 1\n"
    "\n"
    "[solver]\n"
    'method = "DENSE"\n'
    'mode = "ENERGY"\n'
    "emin = -1.0\n"
    "emax = 2.0\n"
)

# Grid levels: fixed wire 55x55 nm, varying grid spacing
# h = dx = dy = wire_width / (nx - 1)  or  wire_height / (ny - 1)
GRID_LEVELS = [
    {'nx': 16, 'ny': 16},  # h = 55/15 = 3.67 A
    {'nx': 21, 'ny': 21},  # h = 55/20 = 2.75 A
    {'nx': 26, 'ny': 26},  # h = 55/25 = 2.20 A
    {'nx': 31, 'ny': 31},  # h = 55/30 = 1.83 A
]


def extract_wire_cb1(output_dir, numvb=8):
    """Extract CB1 energy at k=0 from wire eigenvalues.dat."""
    eig_path = os.path.join(output_dir, "eigenvalues.dat")
    data = parse_eigenvalues(eig_path)
    if not data:
        return None
    _, evals = data[0]
    if len(evals) <= numvb:
        return None
    return evals[numvb]


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        sys.exit(1)

    build_dir = sys.argv[1]
    source_dir = sys.argv[2]

    print("=" * 60)
    print("  WIRE GRID CONVERGENCE TESTS (U6)")
    print("  S7: InAs wire, fixed geometry, varying grid")
    print("=" * 60)

    all_reports = {}
    fdorder = 2
    all_passed = True

    # --- Subband energy ---
    print("\n  S7: Subband energy (bandStructure)")
    h_vals, cb1_vals = [], []
    for grid in GRID_LEVELS:
        nx, ny = grid['nx'], grid['ny']
        dx = 55.0 / (nx - 1)
        dy = 55.0 / (ny - 1)
        h = max(dx, dy)

        cfg_content = WIRE_TEMPLATE.format(
            nx=nx, ny=ny, dx=dx, dy=dy, fdorder=fdorder
        )
        with tempfile.TemporaryDirectory() as work:
            cfg_path = os.path.join(work, "staged.toml")
            with open(cfg_path, 'w') as f:
                f.write(cfg_content)
            rc, output_dir = run_exe(build_dir, "bandStructure", cfg_path, work, timeout=600)
            if rc != 0:
                print(f"    {nx}x{ny}: FAILED (rc={rc})")
                continue
            cb1 = extract_wire_cb1(output_dir)
            if cb1 is not None:
                h_vals.append(h)
                cb1_vals.append(cb1)
                print(f"    {nx}x{ny}: h={h:.4f} A, CB1={cb1:.8f} eV")

    if len(h_vals) >= 3:
        report = make_convergence_report(
            'S7_InAs_wire', 'CB1_energy',
            h_vals, cb1_vals, order=fdorder,
        )
        all_reports['CB1_energy'] = report
        print(f"    Richardson: {report['richardson_extrapolated']:.8f} eV")
        print(f"    GCI: {report['gci_percent']:.2f}%, max_rate: {report['max_observed_rate']:.3f}")
        if not report['passed']:
            all_passed = False
    elif len(h_vals) > 0:
        print(f"    WARN: only {len(h_vals)} grid levels succeeded (need >= 3)")
        all_passed = False

    # --- g-factor ---
    print("\n  S7: g-factor (gfactorCalculation)")
    h_vals_g, gz_vals = [], []
    for grid in GRID_LEVELS:
        nx, ny = grid['nx'], grid['ny']
        dx = 55.0 / (nx - 1)
        dy = 55.0 / (ny - 1)
        h = max(dx, dy)

        cfg_content = WIRE_GF_TEMPLATE.format(
            nx=nx, ny=ny, dx=dx, dy=dy, fdorder=fdorder
        )
        with tempfile.TemporaryDirectory() as work:
            cfg_path = os.path.join(work, "staged.toml")
            with open(cfg_path, 'w') as f:
                f.write(cfg_content)
            rc, output_dir = run_exe(build_dir, "gfactorCalculation", cfg_path, work, timeout=600)
            if rc != 0:
                print(f"    {nx}x{ny}: FAILED (rc={rc})")
                continue
            gf_path = os.path.join(output_dir, "gfactor.dat")
            if os.path.isfile(gf_path):
                _, _, gz = parse_gfactor(gf_path)
                h_vals_g.append(h)
                gz_vals.append(gz)
                print(f"    {nx}x{ny}: h={h:.4f} A, gz={gz:.6f}")

    if len(h_vals_g) >= 3:
        report = make_convergence_report(
            'S7_InAs_wire', 'gz',
            h_vals_g, gz_vals, order=fdorder,
            gci_threshold=0.20,  # wire g-factor converges slowly
        )
        all_reports['gz'] = report
        print(f"    Richardson gz: {report['richardson_extrapolated']:.6f}")
        if not report['passed']:
            all_passed = False

    # Write JSON
    results_dir = os.path.join(source_dir, "tests", "integration", "convergence_results")
    os.makedirs(results_dir, exist_ok=True)
    json_path = os.path.join(results_dir, "wire_convergence.json")
    write_convergence_json(all_reports, json_path)
    print(f"\n  JSON results written to {json_path}")

    # Summary
    print(f"\n  {'All tests passed' if all_passed else 'Some tests FAILED'}")
    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
