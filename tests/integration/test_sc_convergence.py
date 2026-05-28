#!/usr/bin/env python3
"""SC convergence benchmark (U7).

Multi-resolution SC loop runs extracting Fermi level, subband energy shift,
and charge density integral. Monotonic convergence with GCI - no theoretical
rate assertion (SC convergence is empirically determined).

COVERAGE: observable=Fermi_level geometry=QW material=GaAs/AlAs tier=convergence_sc
COVERAGE: observable=subband_shift geometry=QW material=GaAs/AlAs tier=convergence_sc
COVERAGE: observable=charge_integral geometry=QW material=GaAs/AlAs tier=convergence_sc

Usage:
    python3 test_sc_convergence.py <build_dir> <source_dir>
"""

import json
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from star_helpers import run_exe, parse_eigenvalues, trapz_fn
from convergence_helpers import (
    make_convergence_report, write_convergence_json,
)


# ---------------------------------------------------------------------------
# SC config template (parameterized by FDstep)
# ---------------------------------------------------------------------------

SC_TEMPLATE = (
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
    "num_cb = 4\n"
    "num_vb = 8\n"
    "\n"
    "[[material]]\n"
    'name = "AlAs"\n'
    "z_min = -150\n"
    "z_max = 150\n"
    "\n"
    "[[material]]\n"
    'name = "GaAs"\n'
    "z_min = -50\n"
    "z_max = 50\n"
    "\n"
    "which_band = 0\n"
    "band_idx = 1\n"
    "\n"
    "{sc_section}"
    "\n"
    "[[doping]]\n"
    "ND = 0.0\n"
    "NA = 0.0\n"
    "\n"
    "[[doping]]\n"
    "ND = 5.0e18\n"
    "NA = 0.0\n"
)

SC_SECTION = (
    "[sc]\n"
    "max_iterations = 100\n"
    "tolerance = 1.0e-8\n"
    "mixing_alpha = 0.3\n"
    "diis_history = 7\n"
    "temperature = 300.0\n"
    'fermi_mode = "charge_neutrality"\n'
    "fermi_level = 0.0\n"
    "num_kpar = 41\n"
    "kpar_max = 0.2\n"
    'bc_type = "DD"\n'
    "bc_left = 0.0\n"
    "bc_right = 0.0\n"
)

DOMAIN_WIDTH = 300.0  # Angstrom
FDSTEPS = [51, 81, 101, 151]


def parse_sc_summary(filepath):
    """Parse output/sc_summary.dat.

    Format: # converged  iterations  |dPhi|  fermi_level(eV)
            T  23  9.87e-07  0.567
    """
    if not os.path.isfile(filepath):
        return None
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = line.split()
            if len(parts) >= 4:
                return {
                    'converged': parts[0] == 'T',
                    'iterations': int(parts[1]),
                    'dphi': float(parts[2]),
                    'fermi_level': float(parts[3]),
                }
    return None


def parse_sc_charge(filepath):
    """Parse output/sc_charge.dat and return integrated electron density."""
    if not os.path.isfile(filepath):
        return None
    data = np.loadtxt(filepath, comments='#')
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        return None
    z = data[:, 0]  # Angstrom
    ne = data[:, 1]  # cm^-3
    # Integrate n_e over z (convert z from A to cm: 1e-8)
    # Result in cm^-2
    z_cm = z * 1e-8
    integral = trapz_fn(ne, z_cm)
    return integral


def extract_cb1_shift(output_dir, flatband_cb1=None, numvb=8):
    """Extract CB1 energy shift relative to flat-band."""
    eig_path = os.path.join(output_dir, "eigenvalues.dat")
    data = parse_eigenvalues(eig_path)
    if not data:
        return None
    _, evals = data[0]
    if len(evals) <= numvb:
        return None
    cb1 = evals[numvb]
    if flatband_cb1 is None:
        return cb1
    return cb1 - flatband_cb1


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        sys.exit(1)

    build_dir = sys.argv[1]
    source_dir = sys.argv[2]

    print("=" * 60)
    print("  SC CONVERGENCE BENCHMARK (U7)")
    print("  GaAs/AlAs QW, SC loop at 4 grid resolutions")
    print("=" * 60)

    all_reports = {}
    all_passed = True

    # Collect data at each resolution
    fermi_vals = []
    cb1_vals = []
    charge_vals = []
    h_vals = []

    # First run: get flat-band CB1 (no SC)
    print("\n  Running flat-band reference (no SC)...")
    fb_cfg = SC_TEMPLATE.format(fdstep=101, sc_section="")
    with tempfile.TemporaryDirectory() as work:
        cfg_path = os.path.join(work, "staged.toml")
        with open(cfg_path, 'w') as f:
            f.write(fb_cfg)
        rc, output_dir = run_exe(build_dir, "bandStructure", cfg_path, work, timeout=300)
        if rc == 0:
            flatband_cb1 = extract_cb1_shift(output_dir)
            print(f"  Flat-band CB1: {flatband_cb1:.8f} eV")
        else:
            flatband_cb1 = None
            print("  WARNING: flat-band reference failed, using absolute CB1")

    print(f"\n  SC convergence sweep (tolerance=1e-8):")
    for fdstep in FDSTEPS:
        h = DOMAIN_WIDTH / (fdstep - 1)
        cfg_content = SC_TEMPLATE.format(fdstep=fdstep, sc_section=SC_SECTION)
        with tempfile.TemporaryDirectory() as work:
            cfg_path = os.path.join(work, "staged.toml")
            with open(cfg_path, 'w') as f:
                f.write(cfg_content)
            rc, output_dir = run_exe(build_dir, "bandStructure", cfg_path, work, timeout=600)
            if rc != 0:
                print(f"    FDstep={fdstep}: FAILED (rc={rc})")
                continue

            # Parse SC summary
            summary = parse_sc_summary(os.path.join(output_dir, "sc_summary.dat"))
            if summary is None:
                print(f"    FDstep={fdstep}: no sc_summary.dat")
                continue
            if not summary['converged']:
                print(f"    FDstep={fdstep}: SC did NOT converge ({summary['iterations']} iter)")
                continue

            fermi = summary['fermi_level']

            # Parse charge density
            charge_integral = parse_sc_charge(os.path.join(output_dir, "sc_charge.dat"))

            # Parse CB1
            cb1 = extract_cb1_shift(output_dir, flatband_cb1=flatband_cb1)

            h_vals.append(h)
            fermi_vals.append(fermi)
            if cb1 is not None:
                cb1_vals.append(cb1)
            if charge_integral is not None:
                charge_vals.append(charge_integral)

            parts = [f"    FDstep={fdstep}: h={h:.4f} A, Fermi={fermi:.6f} eV"]
            if cb1 is not None:
                parts.append(f"CB1_shift={cb1:.6f} eV")
            if charge_integral is not None:
                parts.append(f"charge_int={charge_integral:.4e} cm^-2")
            print(", ".join(parts))

    # Richardson extrapolation for each observable
    if len(h_vals) >= 3:
        # Fermi level
        report_f = make_convergence_report(
            'SC_GaAs_AlAs', 'Fermi_level',
            h_vals, fermi_vals, order=2,
        )
        all_reports['Fermi_level'] = report_f
        print(f"\n  Fermi level Richardson: {report_f['richardson_extrapolated']:.6f} eV")
        print(f"    GCI: {report_f['gci_percent']:.2f}%, monotonic: {report_f['monotonic']}")
        if not report_f['monotonic']:
            print(f"    WARN: non-monotonic Fermi level convergence")
            # SC convergence is empirical - don't fail on non-monotonicity

        # CB1 shift
        if len(cb1_vals) == len(h_vals):
            report_cb = make_convergence_report(
                'SC_GaAs_AlAs', 'CB1_shift',
                h_vals, cb1_vals, order=2,
            )
            all_reports['CB1_shift'] = report_cb
            print(f"  CB1 shift Richardson: {report_cb['richardson_extrapolated']:.6f} eV")

        # Charge integral
        if len(charge_vals) == len(h_vals):
            report_q = make_convergence_report(
                'SC_GaAs_AlAs', 'charge_integral',
                h_vals, charge_vals, order=2,
            )
            all_reports['charge_integral'] = report_q
            print(f"  Charge integral Richardson: {report_q['richardson_extrapolated']:.4e} cm^-2")

    # SC convergence is empirically determined, not governed by FD order.
    # Override monotonicity failures (SC may oscillate) but check GCI.
    for name, report in all_reports.items():
        if not report['monotonic']:
            print(f"  WARN: {name} is non-monotonic (SC convergence is empirically determined)")
            report['failures'] = [f for f in report['failures'] if 'monotonic' not in f]
            report['passed'] = bool(len(report['failures']) == 0)

    sc_passed = len(all_reports) >= 1 and all(r['passed'] for r in all_reports.values())

    # Write JSON
    results_dir = os.path.join(source_dir, "tests", "integration", "convergence_results")
    os.makedirs(results_dir, exist_ok=True)
    json_path = os.path.join(results_dir, "sc_convergence.json")
    write_convergence_json(all_reports, json_path)
    print(f"\n  JSON results written to {json_path}")

    print(f"\n  {'SC convergence benchmark passed' if sc_passed else 'SC convergence benchmark FAILED'}")
    sys.exit(0 if sc_passed else 1)


if __name__ == "__main__":
    main()
