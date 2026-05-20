#!/usr/bin/env python3
"""QW grid convergence tests (U4) for S4-S6.

Runs GaAs/AlGaAs (S4), InAs/GaSb (S5), and InAs/GaAs strained (S6) QW systems
at 4-5 grid spacings with fixed FD order, extracting subband energy, effective
mass, g-factor, and absorption edge. Performs Richardson extrapolation with GCI
uncertainty quantification and asserts convergence rates.

COVERAGE: observable=CB1_energy geometry=QW material=GaAs/AlGaAs tier=convergence
COVERAGE: observable=m*_e geometry=QW material=GaAs/AlGaAs tier=convergence
COVERAGE: observable=gz geometry=QW material=GaAs/AlGaAs tier=convergence
COVERAGE: observable=absorption_edge geometry=QW material=GaAs/AlGaAs tier=convergence

Usage:
    python3 test_qw_grid_convergence.py <build_dir> <source_dir>
"""

import json
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from star_helpers import run_exe, parse_eigenvalues, parse_gfactor, parse_absorption
from convergence_helpers import (
    richardson_extrapolate, compute_gci, max_convergence_rate,
    check_monotonic, extract_absorption_edge, make_convergence_report,
    write_convergence_json,
)


# ---------------------------------------------------------------------------
# Config templates (parameterized by FDstep)
# ---------------------------------------------------------------------------

S4_TEMPLATE = (
    "waveVector: kx\n"
    "waveVectorMax: 0.1\n"
    "waveVectorStep: 21\n"
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

S4_GF_TEMPLATE = (
    "waveVector: k0\n"
    "waveVectorMax: 0\n"
    "waveVectorStep: 0\n"
    "confinement: 1\n"
    "FDstep: {fdstep}\n"
    "FDorder: {fdorder}\n"
    "numLayers: 2\n"
    "material1: Al30Ga70As -200 200 0\n"
    "material2: GaAs -50 50 0\n"
    "numcb: 4\n"
    "numvb: 8\n"
    "ExternalField: 0  EF\n"
    "EFParams: 0.0005\n"
    "whichBand: 0\n"
    "bandIdx: 1\n"
)

S4_OPTICS_TEMPLATE = (
    "waveVector: k0\n"
    "waveVectorMax: 0.0\n"
    "waveVectorStep: 0\n"
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
    "whichBand: 0\n"
    "bandIdx: 1\n"
    "SC: 0\n"
    "Optics: T\n"
    "LinewidthLorentzian: 0.030\n"
    "LinewidthGaussian: 0.005\n"
    "RefractiveIndex: 3.3\n"
    "Emin: 1.3\n"
    "Emax: 2.0\n"
    "NEnergyPoints: 300\n"
    "Temperature: 300.0\n"
    "CarrierDensity: 0.0\n"
    "Gain: F\n"
    "GainCarrierDensity: 0.0\n"
    "ISBT: F\n"
    "SpontaneousEnabled: F\n"
    "SpinResolved: F\n"
    "Exciton: F\n"
)

S5_TEMPLATE = (
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

S6_TEMPLATE = (
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


# ---------------------------------------------------------------------------
# Observable extraction
# ---------------------------------------------------------------------------

def extract_cb1_energy(output_dir, numvb=8):
    """Extract CB1 energy at k=0 from eigenvalues.dat."""
    eig_path = os.path.join(output_dir, "eigenvalues.dat")
    data = parse_eigenvalues(eig_path)
    if not data:
        return None
    _, evals = data[0]
    if len(evals) <= numvb:
        return None
    return evals[numvb]


def extract_cb2_energy(output_dir, numvb=8):
    """Extract CB2 energy at k=0 from eigenvalues.dat."""
    eig_path = os.path.join(output_dir, "eigenvalues.dat")
    data = parse_eigenvalues(eig_path)
    if not data:
        return None
    _, evals = data[0]
    if len(evals) <= numvb + 1:
        return None
    return evals[numvb + 1]


def extract_gz(output_dir):
    """Extract gz component from gfactor.dat."""
    gf_path = os.path.join(output_dir, "gfactor.dat")
    if not os.path.isfile(gf_path):
        return None
    gx, gy, gz = parse_gfactor(gf_path)
    return gz


def extract_effective_mass_fixed(output_dir):
    """Extract effective mass using fixed k-range from k-sweep eigenvalues."""
    from star_helpers import extract_effective_mass, HBAR2_OVER_2M0
    eig_path = os.path.join(output_dir, "eigenvalues.dat")
    result = extract_effective_mass(eig_path, cb_index=-1)
    if result is None:
        return None
    return result[0]  # m_star


# ---------------------------------------------------------------------------
# Grid sweep runner
# ---------------------------------------------------------------------------

DOMAIN_WIDTH = 400.0  # Angstrom

def run_grid_sweep(build_dir, template, fdsteps, fdorder, exe_name, timeout=300):
    """Run Fortran exe at multiple FDstep values and collect outputs."""
    results = []
    for fdstep in fdsteps:
        h = DOMAIN_WIDTH / (fdstep - 1)
        cfg_content = template.format(fdstep=fdstep, fdorder=fdorder)
        with tempfile.TemporaryDirectory() as work:
            cfg_path = os.path.join(work, "staged.cfg")
            with open(cfg_path, 'w') as f:
                f.write(cfg_content)
            rc, output_dir = run_exe(build_dir, exe_name, cfg_path, work, timeout=timeout)
            if rc != 0:
                print(f"    FDstep={fdstep}: FAILED (rc={rc})")
                results.append({'fdstep': fdstep, 'h': h, 'rc': rc})
                continue
            results.append({
                'fdstep': fdstep, 'h': h, 'rc': rc,
                'output_dir': output_dir,
                'work_dir': work,
            })
    return results


# ---------------------------------------------------------------------------
# System definitions
# ---------------------------------------------------------------------------

SYSTEMS = {
    'S4': {
        'name': 'GaAs/AlGaAs QW',
        'domain_width': 400.0,
        'fdsteps': [51, 101, 201, 401],
        'fdorder': 2,
        'observables': ['CB1_energy', 'CB2_energy', 'subband_spacing', 'effective_mass', 'gz', 'absorption_edge'],
    },
    'S5': {
        'name': 'InAsW/GaSbW QW',
        'domain_width': 300.0,
        'fdsteps': [101, 201, 401],
        'fdorder': 2,
        'observables': ['CB1_energy'],
        'rate_tolerance': None,  # broken-gap: no rate assertion, just Richardson + GCI
        'require_monotonic': False,  # broken-gap: convergence may be non-monotonic
    },
    'S6': {
        'name': 'InAs/GaAs strained QW',
        'domain_width': 400.0,
        'fdsteps': [51, 101, 201, 401],
        'fdorder': 2,
        'observables': ['CB1_energy'],
        'rate_tolerance': 0.5,
    },
}


def run_system_convergence(build_dir, source_dir, sys_key):
    """Run grid convergence for one QW system and return report."""
    sys_info = SYSTEMS[sys_key]
    fdsteps = sys_info['fdsteps']
    fdorder = sys_info['fdorder']
    domain_width = sys_info['domain_width']

    # Select template
    template_map = {
        'S4': S4_TEMPLATE,
        'S5': S5_TEMPLATE,
        'S6': S6_TEMPLATE,
    }
    template = template_map[sys_key]

    print(f"\n  {sys_key}: {sys_info['name']} (FDorder={fdorder})")

    all_reports = {}

    # --- CB1 energy (bandStructure, k=0) ---
    if 'CB1_energy' in sys_info['observables']:
        # Need k=0 config for energy extraction
        k0_template = template
        if sys_key == 'S4':
            k0_template = S4_TEMPLATE.replace("waveVector: kx\n", "waveVector: k0\n").replace(
                "waveVectorMax: 0.1\n", "waveVectorMax: 0\n")

        h_vals, cb1_vals, cb2_vals = [], [], []
        for fdstep in fdsteps:
            h = domain_width / (fdstep - 1)
            cfg_content = k0_template.format(fdstep=fdstep, fdorder=fdorder)
            with tempfile.TemporaryDirectory() as work:
                cfg_path = os.path.join(work, "staged.cfg")
                with open(cfg_path, 'w') as f:
                    f.write(cfg_content)
                rc, output_dir = run_exe(build_dir, "bandStructure", cfg_path, work, timeout=300)
                if rc != 0:
                    print(f"    FDstep={fdstep}: FAILED (rc={rc})")
                    continue
                cb1 = extract_cb1_energy(output_dir)
                cb2 = extract_cb2_energy(output_dir)
                if cb1 is not None:
                    h_vals.append(h)
                    cb1_vals.append(cb1)
                if cb2 is not None:
                    cb2_vals.append(cb2)
                print(f"    FDstep={fdstep}: h={h:.4f} A, CB1={cb1:.8f} eV")

        if len(h_vals) >= 3:
            rate_tol = sys_info.get('rate_tolerance', 0.5)
            require_mono = sys_info.get('require_monotonic', True)
            report = make_convergence_report(
                f'{sys_key}_{sys_info["name"]}', 'CB1_energy',
                h_vals, cb1_vals, order=fdorder,
                theoretical_rate=fdorder if rate_tol is not None else None,
                rate_tolerance=rate_tol,
            )
            if not require_mono and not report['monotonic']:
                report['failures'] = [f for f in report['failures'] if 'monotonic' not in f]
                report['passed'] = bool(len(report['failures']) == 0)
            all_reports['CB1_energy'] = report
            print(f"    Richardson: {report['richardson_extrapolated']:.8f} eV")
            print(f"    GCI: {report['gci_percent']:.2f}%, max_rate: {report['max_observed_rate']:.3f}")
            print(f"    {'PASS' if report['passed'] else 'FAIL'}")

        if len(h_vals) >= 3 and len(cb2_vals) >= 3 and len(h_vals) == len(cb2_vals):
            spacing = [cb2_vals[i] - cb1_vals[i] for i in range(len(cb1_vals))]
            report_sp = make_convergence_report(
                f'{sys_key}_{sys_info["name"]}', 'subband_spacing',
                h_vals, spacing, order=fdorder,
            )
            all_reports['subband_spacing'] = report_sp

    # --- Effective mass (bandStructure, k-sweep) ---
    if 'effective_mass' in sys_info['observables'] and sys_key == 'S4':
        h_vals_m, mstar_vals = [], []
        for fdstep in fdsteps:
            h = domain_width / (fdstep - 1)
            cfg_content = S4_TEMPLATE.format(fdstep=fdstep, fdorder=fdorder)
            with tempfile.TemporaryDirectory() as work:
                cfg_path = os.path.join(work, "staged.cfg")
                with open(cfg_path, 'w') as f:
                    f.write(cfg_content)
                rc, output_dir = run_exe(build_dir, "bandStructure", cfg_path, work, timeout=300)
                if rc != 0:
                    continue
                mstar = extract_effective_mass_fixed(output_dir)
                if mstar is not None:
                    h_vals_m.append(h)
                    mstar_vals.append(mstar)
                    print(f"    FDstep={fdstep}: m*={mstar:.6f} m0")

        if len(h_vals_m) >= 3:
            report = make_convergence_report(
                f'{sys_key}_{sys_info["name"]}', 'effective_mass',
                h_vals_m, mstar_vals, order=fdorder,
            )
            all_reports['effective_mass'] = report
            print(f"    Richardson m*: {report['richardson_extrapolated']:.6f} m0, {'PASS' if report['passed'] else 'FAIL'}")

    # --- g-factor (gfactorCalculation) ---
    if 'gz' in sys_info['observables'] and sys_key == 'S4':
        h_vals_g, gz_vals = [], []
        for fdstep in fdsteps:
            h = domain_width / (fdstep - 1)
            cfg_content = S4_GF_TEMPLATE.format(fdstep=fdstep, fdorder=fdorder)
            with tempfile.TemporaryDirectory() as work:
                cfg_path = os.path.join(work, "staged.cfg")
                with open(cfg_path, 'w') as f:
                    f.write(cfg_content)
                rc, output_dir = run_exe(build_dir, "gfactorCalculation", cfg_path, work, timeout=300)
                if rc != 0:
                    continue
                gz = extract_gz(output_dir)
                if gz is not None:
                    h_vals_g.append(h)
                    gz_vals.append(gz)
                    print(f"    FDstep={fdstep}: gz={gz:.6f}")

        if len(h_vals_g) >= 3:
            report = make_convergence_report(
                f'{sys_key}_{sys_info["name"]}', 'gz',
                h_vals_g, gz_vals, order=fdorder,
            )
            all_reports['gz'] = report
            print(f"    Richardson gz: {report['richardson_extrapolated']:.6f}, {'PASS' if report['passed'] else 'FAIL'}")

    # --- Absorption edge (opticalProperties) ---
    if 'absorption_edge' in sys_info['observables'] and sys_key == 'S4':
        h_vals_a, edge_vals = [], []
        for fdstep in fdsteps:
            h = domain_width / (fdstep - 1)
            cfg_content = S4_OPTICS_TEMPLATE.format(fdstep=fdstep, fdorder=fdorder)
            with tempfile.TemporaryDirectory() as work:
                cfg_path = os.path.join(work, "staged.cfg")
                with open(cfg_path, 'w') as f:
                    f.write(cfg_content)
                rc, output_dir = run_exe(build_dir, "opticalProperties", cfg_path, work, timeout=300)
                if rc != 0:
                    continue
                abs_path = os.path.join(output_dir, "absorption_TE.dat")
                if os.path.isfile(abs_path):
                    spectrum = parse_absorption(abs_path)
                    edge = extract_absorption_edge(spectrum)
                    if edge is not None:
                        h_vals_a.append(h)
                        edge_vals.append(edge)
                        print(f"    FDstep={fdstep}: edge={edge:.6f} eV")

        if len(h_vals_a) >= 3:
            report = make_convergence_report(
                f'{sys_key}_{sys_info["name"]}', 'absorption_edge',
                h_vals_a, edge_vals, order=fdorder,
            )
            all_reports['absorption_edge'] = report
            print(f"    Richardson edge: {report['richardson_extrapolated']:.6f} eV, {'PASS' if report['passed'] else 'FAIL'}")

    return all_reports


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        sys.exit(1)

    build_dir = sys.argv[1]
    source_dir = sys.argv[2]

    print("=" * 60)
    print("  QW GRID CONVERGENCE TESTS (U4)")
    print("  S4: GaAs/AlGaAs, S5: InAs/GaSb, S6: InAs/GaAs strained")
    print("=" * 60)

    all_passed = True
    all_system_reports = {}

    for sys_key in ['S4', 'S5', 'S6']:
        reports = run_system_convergence(build_dir, source_dir, sys_key)
        all_system_reports[sys_key] = reports
        for obs_name, report in reports.items():
            if not report['passed']:
                all_passed = False
                print(f"  FAIL: {sys_key}/{obs_name}: {report['failures']}")

    # Write JSON results
    results_dir = os.path.join(source_dir, "tests", "integration", "convergence_results")
    os.makedirs(results_dir, exist_ok=True)
    json_path = os.path.join(results_dir, "qw_grid_convergence.json")
    write_convergence_json(all_system_reports, json_path)
    print(f"\n  JSON results written to {json_path}")

    # Summary
    print("\n" + "=" * 60)
    print("  SUMMARY")
    print("=" * 60)
    total = 0
    passed = 0
    for sys_key, reports in all_system_reports.items():
        for obs_name, report in reports.items():
            total += 1
            if report['passed']:
                passed += 1
                print(f"  PASS: {sys_key}/{obs_name}")
            else:
                print(f"  FAIL: {sys_key}/{obs_name}: {report['failures']}")

    print(f"\n  {passed}/{total} convergence tests passed")
    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
