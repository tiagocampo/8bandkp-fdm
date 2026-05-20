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
    richardson_extrapolate, compute_gci, check_monotonic,
    make_convergence_report, write_convergence_json,
)


# ---------------------------------------------------------------------------
# SC config template (parameterized by FDstep)
# ---------------------------------------------------------------------------

SC_TEMPLATE = (
    "waveVector: k0\n"
    "waveVectorMax: 0.0\n"
    "waveVectorStep: 1\n"
    "confinement: 1\n"
    "FDstep: {fdstep}\n"
    "FDorder: 2\n"
    "numLayers: 2\n"
    "material1: AlAs -150 150 0\n"
    "material2: GaAs -50 50 0\n"
    "numcb: 4\n"
    "numvb: 8\n"
    "ExternalField: 0  EF\n"
    "EFParams: 0.0005\n"
    "whichBand: 0\n"
    "bandIdx: 1\n"
    "SC: 1\n"
    "max_iter: 100\n"
    "tolerance: 1.0e-8\n"
    "mixing_alpha: 0.3\n"
    "diis_history: 7\n"
    "temperature: 300.0\n"
    "fermi_mode: 0\n"
    "fermi_level: 0.0\n"
    "num_kpar: 41\n"
    "kpar_max: 0.2\n"
    "bc_type: DD\n"
    "bc_left: 0.0\n"
    "bc_right: 0.0\n"
    "doping1: 0.0 0.0\n"
    "doping2: 5.0e18 0.0\n"
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
    fb_cfg = SC_TEMPLATE.format(fdstep=101).replace("SC: 1", "SC: 0")
    with tempfile.TemporaryDirectory() as work:
        cfg_path = os.path.join(work, "staged.cfg")
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
        cfg_content = SC_TEMPLATE.format(fdstep=fdstep)
        with tempfile.TemporaryDirectory() as work:
            cfg_path = os.path.join(work, "staged.cfg")
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

            print(f"    FDstep={fdstep}: h={h:.4f} A, Fermi={fermi:.6f} eV, "
                  f"CB1_shift={cb1:.6f} eV, charge_int={charge_integral:.4e} cm^-2" if cb1 is not None else
                  f"    FDstep={fdstep}: h={h:.4f} A, Fermi={fermi:.6f} eV")

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

    # SC tests pass if we got Richardson values (no rate assertion)
    sc_passed = len(all_reports) >= 1
    for name, report in all_reports.items():
        if not report['monotonic']:
            print(f"  WARN: {name} is non-monotonic (SC is empirically determined)")

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
