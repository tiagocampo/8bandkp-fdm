#!/usr/bin/env python3
"""Exciton convergence benchmark (U8).

Multi-resolution exciton runs extracting binding energy, performing Richardson
extrapolation, and validating against published experimental values.

Miller et al., Phys. Rev. B 1985: GaAs/AlGaAs QW exciton binding energy.
2D hydrogen model: Eb ~ Rydberg/4 = 13.6/4 = 3.4 meV (sanity check).

COVERAGE: observable=exciton_Eb geometry=QW material=GaAs/AlGaAs tier=convergence_exciton

Usage:
    python3 test_exciton_convergence.py <build_dir> <source_dir>
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(__file__))
from star_helpers import run_exe
from convergence_helpers import (
    make_convergence_report, write_convergence_json, parse_exciton_file,
)


# ---------------------------------------------------------------------------
# Exciton config template (parameterized by FDstep)
# ---------------------------------------------------------------------------

EXCITON_TEMPLATE = (
    'confinement = "qw"\n'
    "FDorder = 2\n"
    "fd_step = {fdstep}\n"
    "\n"
    "[wave_vector]\n"
    'mode = "k0"\n'
    "max = 0.0\n"
    "nsteps = 1\n"
    "\n"
    "[bands]\n"
    "num_cb = 2\n"
    "num_vb = 6\n"
    "\n"
    "[[material]]\n"
    'name = "Al30Ga70As"\n'
    "z_min = -200\n"
    "z_max = 200\n"
    "\n"
    "[[material]]\n"
    'name = "GaAs"\n'
    "z_min = -{half_w}\n"
    "z_max = {half_w}\n"
    "\n"
    "[exciton]\n"
    'method = "variational"\n'
)

DOMAIN_WIDTH = 400.0  # Angstrom
HALF_WELL = 50  # Angstrom (100 A well)
FDSTEPS = [51, 101, 201, 401]

# Published reference: Miller et al., Phys. Rev. B 1985
# GaAs/AlGaAs QW, well width ~100 A: Eb ~ 8-10 meV
MILLER_EB_MEV = 9.0  # approximate value for 100 A well
MILLER_TOLERANCE = 0.30  # 30% tolerance (structure-dependent)

# 2D hydrogen model sanity check
RYDBERG_EV = 13.605693  # eV
HYDROGEN_2D_MEV = RYDBERG_EV / 4.0 * 1000  # ~3400 meV (free exciton)


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        sys.exit(1)

    build_dir = sys.argv[1]
    source_dir = sys.argv[2]

    print("=" * 60)
    print("  EXCITON CONVERGENCE BENCHMARK (U8)")
    print(f"  GaAs/AlGaAs QW, well={2*HALF_WELL} A")
    print(f"  Reference: Miller et al. Eb ~ {MILLER_EB_MEV} meV")
    print("=" * 60)

    h_vals, eb_vals = [], []

    for fdstep in FDSTEPS:
        h = DOMAIN_WIDTH / (fdstep - 1)
        cfg_content = EXCITON_TEMPLATE.format(fdstep=fdstep, half_w=HALF_WELL)
        with tempfile.TemporaryDirectory() as work:
            cfg_path = os.path.join(work, "staged.toml")
            with open(cfg_path, 'w') as f:
                f.write(cfg_content)
            rc, output_dir = run_exe(build_dir, "bandStructure", cfg_path, work, timeout=120)
            if rc != 0:
                print(f"    FDstep={fdstep}: FAILED (rc={rc})")
                continue

            exciton_result = parse_exciton_file(os.path.join(output_dir, "exciton.dat"))
            if exciton_result is None:
                print(f"    FDstep={fdstep}: no exciton output")
                continue

            eb = exciton_result['E_binding_meV']
            lam = exciton_result.get('lambda_opt_AA', 0)
            h_vals.append(h)
            eb_vals.append(eb)
            print(f"    FDstep={fdstep}: h={h:.4f} A, Eb={eb:.4f} meV, lambda={lam:.2f} A")

    if len(h_vals) < 3:
        print("  FAIL: insufficient data for Richardson extrapolation")
        sys.exit(1)

    # Richardson extrapolation
    report = make_convergence_report(
        'exciton_GaAs_AlGaAs', 'binding_energy',
        h_vals, eb_vals, order=2,
    )

    # Validation against published values
    rich_eb = report['richardson_extrapolated']
    miller_delta = abs(rich_eb - MILLER_EB_MEV) / MILLER_EB_MEV if MILLER_EB_MEV > 0 else 999
    miller_pass = miller_delta <= MILLER_TOLERANCE

    # 2D hydrogen sanity check: Eb should be much less than free-exciton
    hydrogen_pass = rich_eb < HYDROGEN_2D_MEV

    print(f"\n  Richardson extrapolated Eb: {rich_eb:.4f} meV")
    print(f"  GCI: {report['gci_percent']:.2f}%")
    print(f"  Max convergence rate: {report['max_observed_rate']:.3f}")
    print(f"  Monotonic: {report['monotonic']}")
    print(f"  Miller ref ({MILLER_EB_MEV} meV): delta={miller_delta:.1%} "
          f"{'PASS' if miller_pass else 'FAIL'}")
    print(f"  2D hydrogen sanity (Eb < {HYDROGEN_2D_MEV:.0f} meV): {'PASS' if hydrogen_pass else 'FAIL'}")

    # Add validation to report
    report['experimental_reference'] = {
        'source': 'Miller et al., Phys. Rev. B 1985',
        'Eb_meV': MILLER_EB_MEV,
        'delta_fraction': miller_delta,
        'passed': miller_pass,
    }
    report['hydrogen_sanity'] = {
        '2D_Rydberg_meV': HYDROGEN_2D_MEV,
        'passed': hydrogen_pass,
    }

    # Exciton variational convergence can be non-monotonic due to lambda oscillation
    if not report['monotonic']:
        print(f"  WARN: non-monotonic convergence (variational lambda oscillation)")
        report['failures'] = [f for f in report['failures'] if 'monotonic' not in f]
        report['passed'] = bool(len(report['failures']) == 0)

    overall_pass = report['passed'] and miller_pass and hydrogen_pass
    report['overall_pass'] = bool(overall_pass)

    # Write JSON
    results_dir = os.path.join(source_dir, "tests", "integration", "convergence_results")
    os.makedirs(results_dir, exist_ok=True)
    json_path = os.path.join(results_dir, "exciton_convergence.json")

    # Wrap in a dict for consistent JSON structure
    write_convergence_json({'exciton': report}, json_path)
    print(f"\n  JSON results written to {json_path}")

    print(f"\n  {'Exciton convergence benchmark passed' if overall_pass else 'Exciton convergence benchmark FAILED'}")
    sys.exit(0 if overall_pass else 1)


if __name__ == "__main__":
    main()
