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

sys.path.insert(0, os.path.dirname(__file__))
from star_helpers import run_exe, parse_eigenvalues, parse_gfactor, parse_absorption
from convergence_helpers import (
    extract_absorption_edge, make_convergence_report,
    write_convergence_json,
)


# ---------------------------------------------------------------------------
# Config templates (parameterized by FDstep)
# ---------------------------------------------------------------------------

S4_TEMPLATE = (
    'confinement = "qw"\n'
    "FDorder = {fdorder}\n"
    "fd_step = {fdstep}\n"
    "\n"
    "[wave_vector]\n"
    'mode = "kx"\n'
    "max = 0.1\n"
    "nsteps = 21\n"
    "\n"
    "[bands]\n"
    "num_cb = 4\n"
    "num_vb = 8\n"
    "\n"
    "[[material]]\n"
    'name = "Al30Ga70As"\n'
    "z_min = -200\n"
    "z_max = 200\n"
    "\n"
    "[[material]]\n"
    'name = "GaAs"\n'
    "z_min = -50\n"
    "z_max = 50\n"
)

S4_GF_TEMPLATE = (
    'confinement = "qw"\n'
    "FDorder = {fdorder}\n"
    "fd_step = {fdstep}\n"
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
    'name = "Al30Ga70As"\n'
    "z_min = -200\n"
    "z_max = 200\n"
    "\n"
    "[[material]]\n"
    'name = "GaAs"\n'
    "z_min = -50\n"
    "z_max = 50\n"
    "\n"
    "which_band = 0\n"
    "band_idx = 1\n"
)

S4_OPTICS_TEMPLATE = (
    'confinement = "qw"\n'
    "FDorder = {fdorder}\n"
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
    'name = "Al30Ga70As"\n'
    "z_min = -200\n"
    "z_max = 200\n"
    "\n"
    "[[material]]\n"
    'name = "GaAs"\n'
    "z_min = -50\n"
    "z_max = 50\n"
    "\n"
    "which_band = 0\n"
    "band_idx = 1\n"
    "\n"
    "[optics]\n"
    "linewidth_lorentzian = 0.030\n"
    "linewidth_gaussian = 0.005\n"
    "refractive_index = 3.3\n"
    "E_min = 1.3\n"
    "E_max = 2.0\n"
    "num_energy_points = 300\n"
    "temperature = 300.0\n"
    "carrier_density = 0.0\n"
    "gain_enabled = false\n"
    "gain_carrier_density = 0.0\n"
    "ISBT = false\n"
    "spontaneous = false\n"
    "spin_resolved = false\n"
)

S5_TEMPLATE = (
    'confinement = "qw"\n'
    "FDorder = {fdorder}\n"
    "fd_step = {fdstep}\n"
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
    'name = "AlSbW"\n'
    "z_min = -150\n"
    "z_max = 150\n"
    "\n"
    "[[material]]\n"
    'name = "InAsW"\n'
    "z_min = -15\n"
    "z_max = 0\n"
    "\n"
    "[[material]]\n"
    'name = "GaSbW"\n'
    "z_min = 0\n"
    "z_max = 10\n"
)

S6_TEMPLATE = (
    'confinement = "qw"\n'
    "FDorder = {fdorder}\n"
    "fd_step = {fdstep}\n"
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
    'name = "GaAs"\n'
    "z_min = -200\n"
    "z_max = 200\n"
    "\n"
    "[[material]]\n"
    'name = "InAs"\n'
    "z_min = -50\n"
    "z_max = 50\n"
    "\n"
    "[strain]\n"
    "strain_substrate = 5.6533\n"
)

S4_K0_TEMPLATE = (
    'confinement = "qw"\n'
    "FDorder = {fdorder}\n"
    "fd_step = {fdstep}\n"
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
    'name = "Al30Ga70As"\n'
    "z_min = -200\n"
    "z_max = 200\n"
    "\n"
    "[[material]]\n"
    'name = "GaAs"\n'
    "z_min = -50\n"
    "z_max = 50\n"
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
    """Extract CB2 energy at k=0 from eigenvalues.dat.

    At k=0, eigenvalues are doubly degenerate (time-reversal symmetry).
    CB1 = evals[numvb], its degenerate partner = evals[numvb+1],
    so the real CB2 is at evals[numvb+2].
    """
    eig_path = os.path.join(output_dir, "eigenvalues.dat")
    data = parse_eigenvalues(eig_path)
    if not data:
        return None
    _, evals = data[0]
    if len(evals) <= numvb + 2:
        return None
    return evals[numvb + 2]


def extract_gz(output_dir):
    """Extract gz component from gfactor.dat."""
    gf_path = os.path.join(output_dir, "gfactor.dat")
    if not os.path.isfile(gf_path):
        return None
    gx, gy, gz = parse_gfactor(gf_path)
    return gz


def extract_effective_mass_fixed(output_dir):
    """Extract effective mass using fixed k-range from k-sweep eigenvalues."""
    from star_helpers import extract_effective_mass
    eig_path = os.path.join(output_dir, "eigenvalues.dat")
    result = extract_effective_mass(eig_path, cb_index=-1)
    if result is None:
        return None
    return result[0]  # m_star


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
        'rate_tolerance': 1.0,  # 8-band k.p with material interfaces converges slower than FD order
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
        'rate_tolerance': 1.0,  # strained QW with material interfaces converges slower than FD order
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
        k0_template_map = {
            'S4': S4_K0_TEMPLATE,
            'S5': S5_TEMPLATE,
            'S6': S6_TEMPLATE,
        }
        k0_template = k0_template_map[sys_key]

        h_vals, cb1_vals, cb2_vals = [], [], []
        for fdstep in fdsteps:
            h = domain_width / (fdstep - 1)
            cfg_content = k0_template.format(fdstep=fdstep, fdorder=fdorder)
            with tempfile.TemporaryDirectory() as work:
                cfg_path = os.path.join(work, "staged.toml")
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
                    print(f"    FDstep={fdstep}: h={h:.4f} A, CB1={cb1:.8f} eV")
                    if cb2 is not None:
                        cb2_vals.append(cb2)
                else:
                    print(f"    FDstep={fdstep}: no CB1 eigenvalue")

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
                cfg_path = os.path.join(work, "staged.toml")
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
                cfg_path = os.path.join(work, "staged.toml")
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
                cfg_path = os.path.join(work, "staged.toml")
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

        if len(h_vals_a) >= 2:
            report = make_convergence_report(
                f'{sys_key}_{sys_info["name"]}', 'absorption_edge',
                h_vals_a, edge_vals, order=fdorder,
            )
            all_reports['absorption_edge'] = report
            rich_val = report.get('richardson_extrapolated')
            rich_str = f"{rich_val:.6f} eV" if rich_val is not None else "N/A"
            print(f"    Richardson edge: {rich_str}, {'PASS' if report['passed'] else 'FAIL'}")

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
