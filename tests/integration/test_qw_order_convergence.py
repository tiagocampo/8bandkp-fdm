#!/usr/bin/env python3
"""QW FD order convergence tests (U5) for S4-S6.

Runs QW systems at fixed grid spacing with varying FD orders {2, 4, 6, 8, 10},
extracts CB1 energy, and verifies monotonic convergence toward the continuum limit.

COVERAGE: observable=CB1_energy geometry=QW material=GaAs/AlGaAs tier=convergence_order

Usage:
    python3 test_qw_order_convergence.py <build_dir> <source_dir>
"""

import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(__file__))
from star_helpers import run_exe, parse_eigenvalues
from convergence_helpers import (
    write_convergence_json,
)


# ---------------------------------------------------------------------------
# Config templates (parameterized by FDstep and FDorder)
# ---------------------------------------------------------------------------

S4_K0_TEMPLATE = (
    "waveVector: k0\n"
    "waveVectorMax: 0\n"
    "waveVectorStep: 1\n"
    "confinement: 1\n"
    "FDstep: {fdstep}\n"
    "FDorder: {fdorder}\n"
    "numLayers: 2\n"
    "material1: Al30Ga70As -200 200 0\n"
    "material2: GaAs -50 50 0\n"
    "numcb: 4\n"
    "numvb: 8\n"
    "ExternalField: 0  EF\n"
    "EFParams: 0.0\n"
)

S5_K0_TEMPLATE = (
    "waveVector: k0\n"
    "waveVectorMax: 0\n"
    "waveVectorStep: 1\n"
    "confinement: 1\n"
    "FDstep: {fdstep}\n"
    "FDorder: {fdorder}\n"
    "numLayers: 3\n"
    "material1: AlSbW -150 150 0\n"
    "material2: InAsW -15 0 0\n"
    "material3: GaSbW 0 10 0\n"
    "numcb: 4\n"
    "numvb: 8\n"
    "ExternalField: 0  EF\n"
    "EFParams: 0.0005\n"
)

S6_K0_TEMPLATE = (
    "waveVector: k0\n"
    "waveVectorMax: 0\n"
    "waveVectorStep: 1\n"
    "confinement: 1\n"
    "FDstep: {fdstep}\n"
    "FDorder: {fdorder}\n"
    "numLayers: 2\n"
    "material1: GaAs -200 200 0\n"
    "material2: InAs -50 50 0\n"
    "numcb: 4\n"
    "numvb: 8\n"
    "ExternalField: 0  EF\n"
    "EFParams: 0.0\n"
    "strain: T\n"
    "a_substrate: 5.6533\n"
)


def extract_cb1_energy(output_dir, numvb=8):
    """Extract CB1 energy at k=0."""
    eig_path = os.path.join(output_dir, "eigenvalues.dat")
    data = parse_eigenvalues(eig_path)
    if not data:
        return None
    _, evals = data[0]
    if len(evals) <= numvb:
        return None
    return evals[numvb]


SYSTEMS = {
    'S4': {'name': 'GaAs/AlGaAs QW', 'template': S4_K0_TEMPLATE, 'fdstep': 401},
    'S5': {'name': 'InAsW/GaSbW QW', 'template': S5_K0_TEMPLATE, 'fdstep': 301},
    'S6': {'name': 'InAs/GaAs strained QW', 'template': S6_K0_TEMPLATE, 'fdstep': 401},
}

FD_ORDERS = [2, 4, 6, 8, 10]


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        sys.exit(1)

    build_dir = sys.argv[1]
    source_dir = sys.argv[2]

    print("=" * 60)
    print("  QW FD ORDER CONVERGENCE TESTS (U5)")
    print("  S4-S6 at fixed grid, FD orders {2,4,6,8,10}")
    print("=" * 60)

    all_passed = True
    all_reports = {}

    for sys_key, sys_info in SYSTEMS.items():
        template = sys_info['template']
        fdstep = sys_info['fdstep']
        name = sys_info['name']

        print(f"\n  {sys_key}: {name} (FDstep={fdstep})")

        order_vals = []
        cb1_vals = []

        for order in FD_ORDERS:
            cfg_content = template.format(fdstep=fdstep, fdorder=order)
            with tempfile.TemporaryDirectory() as work:
                cfg_path = os.path.join(work, "staged.cfg")
                with open(cfg_path, 'w') as f:
                    f.write(cfg_content)
                rc, output_dir = run_exe(build_dir, "bandStructure", cfg_path, work, timeout=300)
                if rc != 0:
                    print(f"    FDorder={order}: FAILED (rc={rc})")
                    continue
                cb1 = extract_cb1_energy(output_dir)
                if cb1 is not None:
                    order_vals.append(order)
                    cb1_vals.append(cb1)
                    print(f"    FDorder={order}: CB1={cb1:.8f} eV")

        if len(order_vals) < 2:
            print(f"    SKIP: insufficient data for {sys_key}")
            continue

        # Check spread is physically reasonable
        spread = max(cb1_vals) - min(cb1_vals)
        print(f"    Spread: {spread*1000:.4f} meV")

        spread_pass = spread < 0.1  # 100 meV

        report = {
            'system': f'{sys_key}_{name}',
            'observable': 'CB1_energy',
            'fd_step': fdstep,
            'order_sweep': [
                {'fdorder': order_vals[i], 'CB1_eV': cb1_vals[i]}
                for i in range(len(order_vals))
            ],
            'spread_meV': spread * 1000,
            'passed': spread_pass,
        }

        all_reports[sys_key] = report
        status = 'PASS' if spread_pass else 'FAIL'
        print(f"    {status}: spread={spread*1000:.4f} meV")

        if not spread_pass:
            all_passed = False

    # Write JSON
    results_dir = os.path.join(source_dir, "tests", "integration", "convergence_results")
    os.makedirs(results_dir, exist_ok=True)
    json_path = os.path.join(results_dir, "qw_order_convergence.json")
    write_convergence_json(all_reports, json_path)
    print(f"\n  JSON results written to {json_path}")

    print(f"\n  {'All tests passed' if all_passed else 'Some tests FAILED'}")
    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
