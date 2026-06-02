#!/usr/bin/env python3
"""SC charge neutrality hard check (Phase 21, Issue #02).

Runs bandStructure with SC enabled on a GaAs/AlAs QW, parses the charge
density output, and asserts charge neutrality as a hard conservation law.

The self-consistent solver with fermi_mode = "charge_neutrality" adjusts
the Fermi level so that int(n_e) - int(n_h) = int(ND - NA). This is a
conservation law enforced by bisection, not an approximation. We check
that the residual is below 1%.

COVERAGE: observable=charge_neutrality geometry=QW material=GaAs/AlAs ref=conservation_law

Usage:
    python3 verify_sc_benchmark.py <build_dir> <source_dir>
"""

import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from star_helpers import run_exe, trapz_fn


# ---------------------------------------------------------------------------
# SC config — GaAs/AlAs QW, FDstep=101, ND=5e18 cm^-3 in well, tolerance=1e-8
# ---------------------------------------------------------------------------

SC_CONFIG = (
    'confinement = "qw"\n'
    "FDorder = 2\n"
    "fd_step = 101\n"
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
    "\n"
    "[[doping]]\n"
    "ND = 0.0\n"
    "NA = 0.0\n"
    "\n"
    "[[doping]]\n"
    "ND = 5.0e18\n"
    "NA = 0.0\n"
)

# Physical parameters for reference calculation
ND_WELL = 5.0e18       # cm^-3 (donor concentration in GaAs well)
L_WELL_A = 100.0       # Angstrom (well width: -50 to +50)
L_WELL_CM = L_WELL_A * 1e-8  # cm
ND_INTEGRAL = ND_WELL * L_WELL_CM  # areal density in cm^-2

# Tolerance for charge neutrality check
CHARGE_TOL = 0.01  # 1%


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_sc_summary(filepath):
    """Parse output/sc_summary.dat.

    Format: # converged  iterations  |dPhi|  fermi_level(eV)
            T  23  9.87e-07  0.567

    Returns dict or None.
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
    """Parse output/sc_charge.dat and return (z_A, n_e, n_h) arrays.

    File format (QW): z(A) n_e(cm^-3) n_h(cm^-3)
    """
    if not os.path.isfile(filepath):
        return None
    data = np.loadtxt(filepath, comments='#')
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        return None
    z_A = data[:, 0]    # Angstrom
    n_e = data[:, 1]    # cm^-3
    n_h = data[:, 2] if data.shape[1] >= 3 else None
    return z_A, n_e, n_h


def integrate_charge(z_A, n_e):
    """Integrate n_e(z) dz to get areal density in cm^-2.

    z in Angstrom, n_e in cm^-3. Convert z to cm for integration.
    """
    z_cm = z_A * 1e-8
    return trapz_fn(n_e, z_cm)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        sys.exit(1)

    build_dir = sys.argv[1]
    source_dir = sys.argv[2]

    print("=" * 60)
    print("  SC CHARGE NEUTRALITY HARD CHECK (Phase 21, Issue #02)")
    print("  GaAs/AlAs QW, SC with charge_neutrality fermi mode")
    print("=" * 60)

    checks_passed = 0
    checks_failed = 0

    # ------------------------------------------------------------------
    # Run bandStructure with SC enabled, parse output inside tempdir scope
    # ------------------------------------------------------------------
    print("\n  Running bandStructure with SC loop...")

    with tempfile.TemporaryDirectory() as work:
        cfg_path = os.path.join(work, "staged.toml")
        with open(cfg_path, 'w') as f:
            f.write(SC_CONFIG)

        rc, output_dir = run_exe(build_dir, "bandStructure", cfg_path, work,
                                 timeout=600)

        if rc != 0:
            print(f"  FAIL: bandStructure exited with rc={rc}")
            checks_failed += 1
            print_summary(checks_passed, checks_failed)
            sys.exit(1)

        print(f"  bandStructure completed (rc=0)")

        # --------------------------------------------------------------
        # Check 1: SC convergence
        # --------------------------------------------------------------
        print("\n  --- Check 1: SC convergence ---")
        summary_path = os.path.join(output_dir, "sc_summary.dat")
        summary = parse_sc_summary(summary_path)

        if summary is None:
            print(f"  FAIL: could not parse {summary_path}")
            checks_failed += 1
        elif not summary['converged']:
            print(f"  FAIL: SC did NOT converge ({summary['iterations']} iterations)")
            checks_failed += 1
        else:
            print(f"  PASS: SC converged in {summary['iterations']} iterations "
                  f"(|dPhi|={summary['dphi']:.2e}, E_F={summary['fermi_level']:.6f} eV)")
            checks_passed += 1

        # --------------------------------------------------------------
        # Check 2: Charge neutrality conservation law
        #   int(n_e) - int(n_h) = int(ND - NA) = int(ND)  [NA=0]
        # The Fermi bisection enforces this exactly; residual is numerical.
        # --------------------------------------------------------------
        print("\n  --- Check 2: Charge neutrality "
              "|int(n_e) - int(n_h) - int(N_D)| / int(N_D) < 1% ---")
        charge_path = os.path.join(output_dir, "sc_charge.dat")
        parsed = parse_sc_charge(charge_path)

        if parsed is None:
            print(f"  FAIL: could not parse {charge_path}")
            checks_failed += 1
        else:
            z_A, n_e, n_h = parsed

            # Integrate electron and hole densities
            ne_integral = integrate_charge(z_A, n_e)

            if n_h is not None:
                nh_integral = integrate_charge(z_A, n_h)
            else:
                nh_integral = 0.0

            # Net electron areal density
            net_electron = ne_integral - nh_integral

            # Compute reference doping integral on the actual grid.
            # The grid maps z_min=-50, z_max=50 for GaAs onto discrete points
            # with spacing dz = totalSize / (fd_step - 1) = 300 / 100 = 3.0 A.
            # The doping layer covers the well material's grid indices.
            # Build a doping profile matching the Fortran code's map_layer_to_grid:
            # last-material-wins, so GaAs (layer 2) indices overwrite AlAs (layer 1).
            dz_A = z_A[1] - z_A[0] if len(z_A) > 1 else 1.0
            # Well grid points: where |z| <= 50 A (matches Fortran nint mapping)
            well_mask = (z_A >= -50.0 - dz_A) & (z_A <= 50.0 + dz_A)
            # More precise: Fortran nint gives indices covering z=-51..+51 for
            # the well, so match by checking which grid points the Fortran would
            # assign to layer 2.
            nd_profile = np.where(well_mask, ND_WELL, 0.0)
            nd_integral = integrate_charge(z_A, nd_profile)

            print(f"  int(N_D) dz    = {nd_integral:.4e} cm^-2  "
                  f"(grid: {int(np.sum(well_mask))} points, dz={dz_A:.1f} A)")
            print(f"  int(n_e) dz    = {ne_integral:.4e} cm^-2")
            print(f"  int(n_h) dz    = {nh_integral:.4e} cm^-2")
            print(f"  int(n_e)-int(n_h) = {net_electron:.4e} cm^-2")

            if abs(nd_integral) < 1e-20:
                print("  FAIL: reference ND integral is zero")
                checks_failed += 1
            else:
                rel_err = abs(net_electron - nd_integral) / abs(nd_integral)
                print(f"  Relative error: {rel_err:.4e} ({rel_err*100:.4f}%)")

                if rel_err < CHARGE_TOL:
                    print(f"  PASS: charge neutrality satisfied "
                          f"(rel_err < {CHARGE_TOL*100:.0f}%)")
                    checks_passed += 1
                else:
                    print(f"  FAIL: charge neutrality VIOLATED "
                          f"(rel_err={rel_err*100:.2f}% >= {CHARGE_TOL*100:.0f}%)")
                    checks_failed += 1

        # --------------------------------------------------------------
        # Check 3: SC tolerance reached (|dPhi| < requested tolerance)
        # --------------------------------------------------------------
        print("\n  --- Check 3: |dPhi| < requested tolerance (1e-8) ---")
        if summary is not None:
            requested_tol = 1.0e-8
            if summary['dphi'] < requested_tol:
                print(f"  PASS: |dPhi| = {summary['dphi']:.2e} < {requested_tol:.0e}")
                checks_passed += 1
            else:
                print(f"  FAIL: |dPhi| = {summary['dphi']:.2e} >= {requested_tol:.0e}")
                checks_failed += 1
        else:
            print("  SKIP: no summary data available")
            checks_failed += 1

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print_summary(checks_passed, checks_failed)
    sys.exit(0 if checks_failed == 0 else 1)


def print_summary(passed, failed):
    total = passed + failed
    print("\n" + "=" * 60)
    print(f"  Results: {passed}/{total} checks passed, {failed} failed")
    if failed == 0:
        print("  PASS: SC charge neutrality hard check")
    else:
        print("  FAIL: SC charge neutrality hard check")
    print("=" * 60)


if __name__ == "__main__":
    main()
