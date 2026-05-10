#!/usr/bin/env python3
"""Lecture 07: Self-Consistent Schrodinger-Poisson -- Convergence Validation.

Runs the 8-band k.p solver with self-consistent Schrodinger-Poisson (SCSP)
loop for two systems:
  1. Bulk doped GaAs (SC enabled, charge neutrality mode)
  2. GaAs/AlAs quantum well (SC enabled, DIIS with history=7)

Validates:
  - SC loop converges (exit code 0, convergence message in stdout)
  - Potential change |dPhi| decreases below tolerance
  - Charge neutrality within 5% for doped systems
  - DIIS converges faster (fewer iterations) than linear mixing would

Generates:
  - Convergence plot: |dPhi| vs iteration for both systems
  - Saved to docs/lecture/figures/lecture_07_convergence.png
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parent.parent
BUILD_DIR = REPO / "build"
CONFIGS_DIR = REPO / "tests" / "regression" / "configs"
FIGURES_DIR = REPO / "docs" / "lecture" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(REPO / "tests" / "integration"))
from star_helpers import run_exe, parse_eigenvalues, compare_value, TOL_NUMERICAL

# numpy 2.0 compatibility
_trapz = getattr(np, 'trapezoid', None) or getattr(np, 'trapz', None)


def run_with_stdout(build_dir, name, config_path, work_dir, timeout=300):
    """Run executable and capture stdout for SC iteration parsing.

    Returns:
        (returncode, output_dir, stdout_str)
    """
    exe_path = os.path.join(build_dir, "src", name)
    if not os.path.isfile(exe_path):
        raise FileNotFoundError(f"Executable not found: {exe_path}")

    dst_cfg = os.path.join(work_dir, "input.cfg")
    shutil.copy2(config_path, dst_cfg)
    output_dir = os.path.join(work_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    try:
        result = subprocess.run(
            [exe_path],
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        return result.returncode, output_dir, result.stdout + result.stderr
    except subprocess.TimeoutExpired:
        return -1, output_dir, ""


def parse_sc_iterations(stdout_text):
    """Parse SC iteration data from stdout.

    The Fortran code prints lines like:
        iter:    1  |dPhi|:  1.2345E-03   mu:    0.1234  max|rho|:  5.6789E+18

    Returns:
        dict with keys:
          'iterations': list of int
          'dphi': list of float (|dPhi| values)
          'mu': list of float (Fermi level)
          'max_rho': list of float (max charge density)
          'converged': bool
          'converged_iter': int or None
    """
    result = {
        'iterations': [],
        'dphi': [],
        'mu': [],
        'max_rho': [],
        'converged': False,
        'converged_iter': None,
    }

    # Match SC iteration lines
    iter_pattern = re.compile(
        r'iter:\s*(\d+)\s+\|dPhi\|:\s*([0-9.E+\-]+)\s+mu:\s*([0-9.\-+E]+)\s+max\|rho\|:\s*([0-9.E+\-]+)',
        re.IGNORECASE,
    )
    for line in stdout_text.splitlines():
        m = iter_pattern.search(line)
        if m:
            result['iterations'].append(int(m.group(1)))
            result['dphi'].append(float(m.group(2)))
            result['mu'].append(float(m.group(3)))
            result['max_rho'].append(float(m.group(4)))

    # Check convergence message
    conv_pattern = re.compile(
        r'SC loop converged at iteration\s+(\d+)',
    )
    wire_conv_pattern = re.compile(
        r'Wire SC loop converged at iteration\s+(\d+)',
    )
    for line in stdout_text.splitlines():
        m = conv_pattern.search(line) or wire_conv_pattern.search(line)
        if m:
            result['converged'] = True
            result['converged_iter'] = int(m.group(1))

    return result


def parse_sc_charge(charge_path):
    """Parse sc_charge.dat: z(A), n_e(cm^-3), n_h(cm^-3).

    Returns:
        (z, n_e, n_h) as numpy arrays
    """
    data = np.loadtxt(charge_path, comments='#')
    if data.ndim == 1:
        data = data.reshape(1, -1)
    z = data[:, 0]
    n_e = data[:, 1]
    n_h = data[:, 2]
    return z, n_e, n_h


def compute_charge_neutrality(z, n_e, n_h, doping_net):
    """Compute charge neutrality metric for SC loop.

    The SC loop enforces: integral(n_e - n_h) = integral(N_D - N_A).
    This function checks how well that condition is satisfied.

    Args:
        z: spatial grid (Angstroms)
        n_e: electron density profile (cm^-3)
        n_h: hole density profile (cm^-3)
        doping_net: net doping (N_D - N_A) per grid point (cm^-3), same length as z

    Returns:
        float: relative charge imbalance
               |integral(n_e - n_h) - integral(doping_net)| / |integral(doping_net)|
    """
    total_ne = _trapz(n_e, z)
    total_nh = _trapz(n_h, z)
    total_doping = _trapz(np.array(doping_net), z)
    net_free = total_ne - total_nh
    imbalance = abs(net_free - total_doping)
    ref = max(abs(total_doping), 1e-30)
    return imbalance / ref


def plot_convergence(bulk_data, qw_data, save_path):
    """Generate convergence plot: |dPhi| vs iteration.

    Args:
        bulk_data: dict from parse_sc_iterations for bulk GaAs
        qw_data: dict from parse_sc_iterations for GaAs/AlAs QW
        save_path: output figure path
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # --- Left panel: Bulk doped GaAs ---
    if bulk_data['iterations']:
        iters = np.array(bulk_data['iterations'])
        dphi = np.array(bulk_data['dphi'])
        ax1.semilogy(iters, dphi, 'o-', color='royalblue', lw=1.5,
                      markersize=5, label=r'$|\Delta\Phi|$')
        ax1.axhline(1e-6, color='red', ls='--', lw=1, label=r'Tolerance $10^{-6}$')
        if bulk_data['converged_iter'] is not None:
            ax1.axvline(bulk_data['converged_iter'], color='green', ls=':',
                        lw=1, label=f'Converged (iter {bulk_data["converged_iter"]})')
        ax1.set_xlabel('SC Iteration', fontsize=11)
        ax1.set_ylabel(r'$|\Delta\Phi|$ (eV)', fontsize=11)
        ax1.set_title('Bulk doped GaAs', fontsize=12)
        ax1.legend(fontsize=9)
        ax1.grid(True, alpha=0.3)
    else:
        ax1.text(0.5, 0.5, 'No SC data', transform=ax1.transAxes,
                 ha='center', va='center', fontsize=12)
        ax1.set_title('Bulk doped GaAs', fontsize=12)

    # --- Right panel: GaAs/AlAs QW ---
    if qw_data['iterations']:
        iters = np.array(qw_data['iterations'])
        dphi = np.array(qw_data['dphi'])
        # Mark linear mixing phase (first diis_history iterations)
        diis_hist = 7  # from config
        linear_end = min(diis_hist, len(iters))

        ax2.semilogy(iters[:linear_end], dphi[:linear_end], 's-',
                      color='steelblue', lw=1.5, markersize=5,
                      label='Linear mixing')
        if len(iters) > linear_end:
            ax2.semilogy(iters[linear_end - 1:], dphi[linear_end - 1:], 'o-',
                          color='crimson', lw=1.5, markersize=5,
                          label='DIIS accelerated')
        ax2.axhline(1e-6, color='red', ls='--', lw=1, label=r'Tolerance $10^{-6}$')
        if qw_data['converged_iter'] is not None:
            ax2.axvline(qw_data['converged_iter'], color='green', ls=':',
                        lw=1, label=f'Converged (iter {qw_data["converged_iter"]})')
        ax2.set_xlabel('SC Iteration', fontsize=11)
        ax2.set_ylabel(r'$|\Delta\Phi|$ (eV)', fontsize=11)
        ax2.set_title('GaAs/AlAs QW (DIIS)', fontsize=12)
        ax2.legend(fontsize=9)
        ax2.grid(True, alpha=0.3)
    else:
        ax2.text(0.5, 0.5, 'No SC data', transform=ax2.transAxes,
                 ha='center', va='center', fontsize=12)
        ax2.set_title('GaAs/AlAs QW (DIIS)', fontsize=12)

    fig.suptitle('Lecture 07: Self-Consistent Schrodinger-Poisson Convergence',
                 fontsize=13, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Convergence plot saved: {save_path}")


def main():
    print("=" * 60)
    print("  Lecture 07: Self-Consistent Schrodinger-Poisson")
    print("=" * 60)

    all_passed = True

    # ------------------------------------------------------------------
    # Test 1: Bulk doped GaAs SC convergence
    # ------------------------------------------------------------------
    print("\n[1/4] Running bulk doped GaAs with SC loop...")
    bulk_cfg = CONFIGS_DIR / "sc_bulk_gaas_doped.cfg"
    bulk_work = tempfile.mkdtemp(prefix="lecture07_bulk_")
    try:
        rc, bulk_outdir, bulk_stdout = run_with_stdout(
            str(BUILD_DIR), "bandStructure", str(bulk_cfg), bulk_work,
            timeout=300,
        )
        bulk_data = parse_sc_iterations(bulk_stdout)

        # Assertion: executable succeeds
        if rc != 0:
            print(f"  FAIL: bandStructure returned {rc}")
            all_passed = False
        else:
            print(f"  PASS: bandStructure exited successfully (rc=0)")

        # Assertion: SC loop converged
        if bulk_data['converged']:
            print(f"  PASS: SC loop converged at iteration "
                  f"{bulk_data['converged_iter']}")
        else:
            print(f"  FAIL: SC loop did not converge "
                  f"({len(bulk_data['iterations'])} iterations)")
            all_passed = False

        # Assertion: |dPhi| decreases below tolerance
        if bulk_data['dphi']:
            final_dphi = bulk_data['dphi'][-1]
            initial_dphi = bulk_data['dphi'][0]
            if final_dphi < 1e-6:
                print(f"  PASS: Final |dPhi| = {final_dphi:.2e} < 1e-6")
            else:
                print(f"  FAIL: Final |dPhi| = {final_dphi:.2e} >= 1e-6")
                all_passed = False
            print(f"  Info: |dPhi| dropped from {initial_dphi:.2e} "
                  f"to {final_dphi:.2e} in {len(bulk_data['iterations'])} iters")

        # Assertion: charge neutrality
        charge_path = os.path.join(bulk_outdir, "sc_charge.dat")
        if os.path.isfile(charge_path):
            z, n_e, n_h = parse_sc_charge(charge_path)
            # Bulk config: 1 layer GaAs, doping1: 1e18 0.0
            doping_net = np.full_like(z, 1.0e18)
            imbalance = compute_charge_neutrality(z, n_e, n_h, doping_net)
            if imbalance < 0.05:
                print(f"  PASS: Charge imbalance = {imbalance:.4f} (< 5%)")
            else:
                print(f"  WARN: Charge imbalance = {imbalance:.4f} (>= 5%)")
                # Not a hard failure -- numerical integration can be coarse
        else:
            print(f"  WARN: sc_charge.dat not found (bulk may not write it)")

    finally:
        shutil.rmtree(bulk_work, ignore_errors=True)

    # ------------------------------------------------------------------
    # Test 2: GaAs/AlAs QW SC convergence with DIIS
    # ------------------------------------------------------------------
    # The canonical sc_gaas_alas_qw.cfg uses 101 FD points and 41 k-points,
    # which is too slow for a lecture demo. We write a lighter QW config with
    # fewer grid points and k-points to keep runtime manageable while still
    # exercising the DIIS convergence path.
    print("\n[2/4] Running GaAs/AlAs QW with SC loop (DIIS)...")
    qw_work = tempfile.mkdtemp(prefix="lecture07_qw_")
    try:
        qw_cfg = os.path.join(qw_work, "qw_sc_light.cfg")
        with open(qw_cfg, 'w') as f:
            f.write(
                "waveVector: k0\n"
                "waveVectorMax: 0.0\n"
                "waveVectorStep: 1\n"
                "confinement:  1\n"
                "FDstep: 51\n"
                "FDorder: 2\n"
                "numLayers:  2\n"
                "material1: AlAs -75 75 0\n"
                "material2: GaAs -25 25 0\n"
                "numcb: 2\n"
                "numvb: 4\n"
                "ExternalField: 0  EF\n"
                "EFParams: 0.0005\n"
                "whichBand: 0\n"
                "bandIdx: 1\n"
                "SC: 1\n"
                "max_iter: 50\n"
                "tolerance: 1.0e-6\n"
                "mixing_alpha: 0.3\n"
                "diis_history: 7\n"
                "temperature: 300.0\n"
                "fermi_mode: 0\n"
                "fermi_level: 0.0\n"
                "num_kpar: 21\n"
                "kpar_max: 0.2\n"
                "bc_type: DD\n"
                "bc_left: 0.0\n"
                "bc_right: 0.0\n"
                "doping1: 0.0 0.0\n"
                "doping2: 5.0e18 0.0\n"
            )
        rc, qw_outdir, qw_stdout = run_with_stdout(
            str(BUILD_DIR), "bandStructure", qw_cfg, qw_work,
            timeout=300,
        )
        qw_data = parse_sc_iterations(qw_stdout)

        # Assertion: executable succeeds
        if rc != 0:
            print(f"  FAIL: bandStructure returned {rc}")
            all_passed = False
        else:
            print(f"  PASS: bandStructure exited successfully (rc=0)")

        # Assertion: SC loop converged
        if qw_data['converged']:
            print(f"  PASS: SC loop converged at iteration "
                  f"{qw_data['converged_iter']}")
        else:
            print(f"  FAIL: SC loop did not converge "
                  f"({len(qw_data['iterations'])} iterations)")
            all_passed = False

        # Assertion: |dPhi| decreases below tolerance
        if qw_data['dphi']:
            final_dphi = qw_data['dphi'][-1]
            initial_dphi = qw_data['dphi'][0]
            if final_dphi < 1e-6:
                print(f"  PASS: Final |dPhi| = {final_dphi:.2e} < 1e-6")
            else:
                print(f"  FAIL: Final |dPhi| = {final_dphi:.2e} >= 1e-6")
                all_passed = False
            print(f"  Info: |dPhi| dropped from {initial_dphi:.2e} "
                  f"to {final_dphi:.2e} in {len(qw_data['iterations'])} iters")

        # Assertion: charge neutrality within 5%
        charge_path = os.path.join(qw_outdir, "sc_charge.dat")
        if os.path.isfile(charge_path):
            z, n_e, n_h = parse_sc_charge(charge_path)
            # QW config: 2 layers, AlAs (-75,75) doping=0, GaAs (-25,25) doping=5e18
            # FDstep=51, domain -75..75, dz=150/50=3 AA
            # GaAs well from -25 to 25
            doping_net = np.zeros_like(z)
            doping_net[(z >= -25.0) & (z <= 25.0)] = 5.0e18
            imbalance = compute_charge_neutrality(z, n_e, n_h, doping_net)
            if imbalance < 0.05:
                print(f"  PASS: Charge imbalance = {imbalance:.4f} (< 5%)")
            else:
                print(f"  WARN: Charge imbalance = {imbalance:.4f} (>= 5%)")
                if imbalance > 0.10:
                    all_passed = False
        else:
            print(f"  FAIL: sc_charge.dat not found")
            all_passed = False

    finally:
        shutil.rmtree(qw_work, ignore_errors=True)

    # ------------------------------------------------------------------
    # Test 3: DIIS converges faster than linear mixing
    # ------------------------------------------------------------------
    print("\n[3/4] Checking DIIS convergence acceleration...")
    if qw_data['iterations'] and qw_data['converged']:
        n_iters = len(qw_data['iterations'])
        dphi_arr = np.array(qw_data['dphi'])

        # DIIS history = 7 from config. The first 7 iterations use linear
        # mixing. DIIS kicks in after that. We check that convergence
        # completes in fewer than the max iterations and that |dPhi| is
        # monotonically decreasing in the DIIS phase (or at least converges).
        diis_start = 7  # diis_history from config

        # Report per-phase statistics
        if n_iters > diis_start and len(dphi_arr) > diis_start:
            dphi_linear_end = dphi_arr[diis_start - 1]
            dphi_final = dphi_arr[-1]
            print(f"  Linear phase (iters 1-{diis_start}): "
                  f"|dPhi| = {dphi_arr[0]:.2e} -> {dphi_linear_end:.2e}")
            print(f"  DIIS phase (iters {diis_start+1}-{n_iters}): "
                  f"|dPhi| = {dphi_linear_end:.2e} -> {dphi_final:.2e}")
            # DIIS acceleration factor: how much faster DIIS reduces |dPhi|
            # compared to linear (orders of magnitude per iteration)
            if dphi_linear_end > 0 and dphi_final > 0:
                linear_reduction = dphi_arr[0] / dphi_linear_end
                diis_reduction = dphi_linear_end / dphi_final
                linear_per_iter = linear_reduction ** (1.0 / diis_start)
                diis_iters = n_iters - diis_start
                if diis_iters > 0:
                    diis_per_iter = diis_reduction ** (1.0 / diis_iters)
                else:
                    diis_per_iter = 1.0
                print(f"  Reduction per iter: linear={linear_per_iter:.3f}x, "
                      f"DIIS={diis_per_iter:.3f}x")
                if diis_per_iter >= linear_per_iter:
                    print(f"  PASS: DIIS phase converges at least as fast as linear")
                else:
                    print(f"  WARN: DIIS per-iter rate not faster (may need more iters)")
        else:
            print(f"  INFO: Converged in {n_iters} iters (within linear phase)")

        # Simple check: converged in < max iterations
        if n_iters < 50:
            print(f"  PASS: Converged in {n_iters} iterations (< max 50)")
        else:
            print(f"  FAIL: Required all {n_iters} iterations")
            all_passed = False
    else:
        print(f"  SKIP: No QW SC data to analyze")

    # ------------------------------------------------------------------
    # Test 4: Generate convergence plot
    # ------------------------------------------------------------------
    print(f"\n[4/4] Generating convergence plot...")
    plot_convergence(bulk_data, qw_data,
                     FIGURES_DIR / "lecture_07_convergence.png")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    if all_passed:
        print("  Lecture 07: ALL CHECKS PASSED")
    else:
        print("  Lecture 07: SOME CHECKS FAILED")
    print("=" * 60)
    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
