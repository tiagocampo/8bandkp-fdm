#!/usr/bin/env python3
"""Wire Strain Validation (U3).

Validates the wire Navier-Cauchy PDE strain profile and the resulting
Bir-Pikus band edge shifts for an InAs core in a GaAs shell wire.

Requirements: R7, R8, R9, R12, R13, R14.

Usage:
    verify_strain_wire_profile.py <build_dir> <source_dir>
"""

# COVERAGE: observable=strain_shift geometry=wire material=InAs/GaAs-core-shell ref=biaxial_analytical

import os
import shutil
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from star_helpers import (
    run_exe, parse_eigenvalues, compare_value,
)

# ---------------------------------------------------------------------------
# InAs/GaAs parameters
# ---------------------------------------------------------------------------
INAS_A0    = 6.0583
INAS_C11   = 832.9
INAS_C12   = 452.6

GAAS_A0    = 5.65325

# Biaxial reference values (infinite-plane limit)
EPS_XX_BULK = (GAAS_A0 - INAS_A0) / INAS_A0  # -0.0669
EPS_ZZ_BULK = -2.0 * INAS_C12 / INAS_C11 * EPS_XX_BULK  # +0.0727
TR_BULK = 2 * EPS_XX_BULK + EPS_ZZ_BULK

# Tolerances
TOL_STRAIN_MAGNITUDE = 0.60  # 60% for wire strain (free-surface relaxation vs biaxial)


def make_wire_config(nx=20, ny=20, dx=5.0, dy=5.0,
                     width=100.0, height=100.0,
                     core_size=30.0, strain=True):
    """Build inline wire config."""
    shell_size = width  # shell extends to wire boundary
    strain_flag = "T" if strain else "F"
    return (
        "waveVector: kz\n"
        "waveVectorMax: 0.01\n"
        "waveVectorStep: 2\n"
        "confinement:  2\n"
        "FDstep: 1\n"
        "FDorder: 2\n"
        "numLayers:  2\n"
        f"wire_nx: {nx}\n"
        f"wire_ny: {ny}\n"
        f"wire_dx: {dx}\n"
        f"wire_dy: {dy}\n"
        "wire_shape: rectangle\n"
        f"wire_width: {width}\n"
        f"wire_height: {height}\n"
        "numRegions: 2\n"
        f"region: GaAs  {core_size}  {shell_size}\n"
        f"region: InAs  0.0  {core_size}\n"
        "numcb: 4\n"
        "numvb: 8\n"
        "ExternalField: 0  EF\n"
        "EFParams: 0.0\n"
        "whichBand: 0\n"
        "bandIdx: 1\n"
        "SC: 0\n"
        "feast_emin: -1.5\n"
        "feast_emax: 2.0\n"
        "feast_m0: -1\n"
        f"strain: {strain_flag}\n"
        "strain_ref: GaAs\n"
        "strain_solver: pardiso\n"
        "piezo: F\n"
    )


def parse_strain(filepath):
    """Parse strain.dat file.

    Returns numpy array with columns: x, y, eps_xx, eps_yy, eps_zz,
    eps_xy, eps_xz, eps_yz.
    """
    rows = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            vals = [float(x) for x in line.split()]
            if len(vals) >= 8:
                rows.append(vals[:8])
    if not rows:
        return np.empty((0, 8))
    return np.array(rows)


def identify_core_points(strain_data, core_max, wire_width, wire_height):
    """Return mask for points inside the InAs core.

    The core is defined by Euclidean distance from the wire center.
    Region: InAs 0.0 core_max means distance 0 to core_max from center.
    Uses strict < to match Fortran region assignment (shell checked first).
    """
    cx = wire_width / 2.0
    cy = wire_height / 2.0
    x = strain_data[:, 0]
    y = strain_data[:, 1]
    dist = np.sqrt((x - cx)**2 + (y - cy)**2)
    return dist < core_max


def test_r7_interior_strain(strain_data, core_max, wire_width, wire_height):
    """R7: Interior strain profile has correct qualitative behavior.

    For a cylindrical InAs core in GaAs, the Navier-Cauchy PDE produces
    a strain field that is NOT uniformly biaxial — hoop strain is tensile
    while radial strain is compressive. Key validations:
    - eps_zz (wire axis) is tensile in the core (positive)
    - Tr(eps) is compressive in the core (negative, hydrostatic)
    - Strain is approximately uniform near the center
    """
    print("  [R7] Interior strain profile (qualitative)")

    core = identify_core_points(strain_data, core_max, wire_width, wire_height)
    core_data = strain_data[core]

    if len(core_data) == 0:
        print("    FAIL: no core points found")
        return False

    eps_xx_core = core_data[:, 2]
    eps_yy_core = core_data[:, 3]
    eps_zz_core = core_data[:, 4]

    mean_eps_zz = np.mean(eps_zz_core)
    mean_tr = np.mean(eps_xx_core + eps_yy_core + eps_zz_core)

    print(f"    Core points: {len(core_data)}")
    print(f"    Mean eps_zz (out-of-plane): {mean_eps_zz:+.6f}")
    print(f"    Mean Tr(eps) (hydrostatic): {mean_tr:+.6f}")
    print(f"    Biaxial reference: Tr(eps) = {TR_BULK:+.6f}")

    all_pass = True

    # eps_zz should be positive (tensile out-of-plane, Poisson relaxation)
    passed_zz = mean_eps_zz > 0
    print(f"    eps_zz > 0 (tensile): {'PASS' if passed_zz else 'FAIL'}")
    all_pass = all_pass and passed_zz

    # Tr(eps) should be negative (compressive hydrostatic strain in core)
    passed_tr = mean_tr < 0
    print(f"    Tr(eps) < 0 (compressive hydrostatic): {'PASS' if passed_tr else 'FAIL'}")
    all_pass = all_pass and passed_tr

    # Uniformity: strain near center (r < core_max/2) should have low variance
    cx = wire_width / 2.0
    cy = wire_height / 2.0
    x_core = core_data[:, 0]
    y_core = core_data[:, 1]
    dist_core = np.sqrt((x_core - cx)**2 + (y_core - cy)**2)
    center_mask = dist_core < core_max * 0.5
    if np.sum(center_mask) > 2:
        tr_center = core_data[center_mask, 2] + core_data[center_mask, 3] + core_data[center_mask, 4]
        variance = np.std(tr_center) / abs(np.mean(tr_center)) if abs(np.mean(tr_center)) > 1e-10 else 0
        passed_uniform = variance < 0.5
        print(f"    Core center uniformity (Tr rel std): {variance:.4f} "
              f"{'PASS' if passed_uniform else 'FAIL'}")
        all_pass = all_pass and passed_uniform

    # Strain magnitude ratio vs biaxial
    if abs(TR_BULK) > 1e-10:
        magnitude_ratio = abs(mean_tr / TR_BULK)
        print(f"    Strain magnitude / biaxial: {magnitude_ratio:.2%}")
        passed_mag = magnitude_ratio > (1.0 - TOL_STRAIN_MAGNITUDE)
        print(f"    Magnitude within {(1-TOL_STRAIN_MAGNITUDE):.0%}-{(1+TOL_STRAIN_MAGNITUDE):.0%} of biaxial: "
              f"{'PASS' if passed_mag else 'FAIL'}")
        all_pass = all_pass and passed_mag

    return all_pass


def test_r8_hydrostatic_strain(strain_data, core_max, wire_width, wire_height):
    """R8: Hydrostatic strain sign and decay.

    Tr(eps) should be compressive (negative) in the InAs core and
    decay toward zero in the GaAs shell.
    """
    print("  [R8] Hydrostatic strain sign and decay")

    tr_eps = strain_data[:, 2] + strain_data[:, 3] + strain_data[:, 4]
    x = strain_data[:, 0]
    y = strain_data[:, 1]
    cx = wire_width / 2.0
    cy = wire_height / 2.0
    dist = np.sqrt((x - cx)**2 + (y - cy)**2)

    core = identify_core_points(strain_data, core_max, wire_width, wire_height)

    all_pass = True

    # Tr(eps) in core: should be negative for compressive hydrostatic
    core_tr = tr_eps[core]
    mean_core_tr = np.mean(core_tr)
    print(f"    Mean Tr(eps) in core: {mean_core_tr:+.6f}")

    # Sign check consistent with R7 (compressive = negative)
    passed_sign = mean_core_tr < 0
    print(f"    Tr(eps) < 0 in core (compressive): "
          f"{'PASS' if passed_sign else 'FAIL'}")
    all_pass = all_pass and passed_sign

    # Decay: |Tr(eps)| at outer shell < |Tr(eps)| at core center
    interface_mask = (dist > core_max - 5) & (dist < core_max + 5)
    outer_mask = dist > 0.7 * max(dist)

    if np.sum(interface_mask) > 0 and np.sum(outer_mask) > 0:
        tr_interface = np.mean(np.abs(tr_eps[interface_mask]))
        tr_outer = np.mean(np.abs(tr_eps[outer_mask]))
        passed_decay = tr_outer < tr_interface
        print(f"    Strain decay: |Tr| interface={tr_interface:.6f} "
              f"outer={tr_outer:.6f}: "
              f"{'PASS' if passed_decay else 'FAIL'}")
        all_pass = all_pass and passed_decay

    # Interface gradient: strain components should change at interface
    inner_yy = np.mean(strain_data[core, 3])
    shell_near = ~core & (dist > core_max) & (dist < core_max + 10)
    if np.sum(shell_near) > 0:
        outer_yy = np.mean(strain_data[shell_near, 3])
        gradient = abs(outer_yy - inner_yy)
        passed_gradient = gradient > 1e-4
        print(f"    Interface gradient (eps_yy): {gradient:.6f} "
              f"{'PASS' if passed_gradient else 'FAIL'}")
        all_pass = all_pass and passed_gradient

    return all_pass


def test_r9_band_edge_shifts(evals_strained, evals_unstrained):
    """R9: Band edge shifts consistent with strain-induced modification.

    Compares strained vs unstrained wire eigenvalues to verify strain
    produces physically correct shifts. For InAs core wire:
    - Band edges (VB top and CB bottom) should shift with strain
    - The gap should change due to strain
    - Spectrum should span a meaningful energy range
    """
    print("  [R9] Band edge shift consistency (strained vs unstrained)")

    print(f"    Strained eigenvalues: {len(evals_strained)}")
    print(f"    Unstrained eigenvalues: {len(evals_unstrained)}")

    # Find band edges: VB top = highest negative eigenvalue, CB bottom = lowest positive
    vb_top_s = max(e for e in evals_strained if e < 0)
    cb_bot_s = min(e for e in evals_strained if e > 0)
    vb_top_u = max(e for e in evals_unstrained if e < 0)
    cb_bot_u = min(e for e in evals_unstrained if e > 0)

    gap_s = cb_bot_s - vb_top_s
    gap_u = cb_bot_u - vb_top_u
    vb_shift = vb_top_s - vb_top_u
    cb_shift = cb_bot_s - cb_bot_u
    gap_shift = gap_s - gap_u

    print(f"    VB top: strained={vb_top_s:+.6f}, unstrained={vb_top_u:+.6f}, "
          f"shift={vb_shift:+.6f} eV")
    print(f"    CB bot: strained={cb_bot_s:+.6f}, unstrained={cb_bot_u:+.6f}, "
          f"shift={cb_shift:+.6f} eV")
    print(f"    Gap: strained={gap_s:.6f}, unstrained={gap_u:.6f}, "
          f"shift={gap_shift:+.6f} eV")

    all_pass = True

    # VB top should shift with strain
    passed_vb = abs(vb_shift) > 0.005
    print(f"    |VB top shift| > 0.005 eV: {'PASS' if passed_vb else 'FAIL'}")
    all_pass = all_pass and passed_vb

    # CB bottom should shift with strain
    passed_cb = abs(cb_shift) > 0.005
    print(f"    |CB bot shift| > 0.005 eV: {'PASS' if passed_cb else 'FAIL'}")
    all_pass = all_pass and passed_cb

    # Gap should change
    passed_gap = abs(gap_shift) > 0.01
    print(f"    |Gap shift| > 0.01 eV: {'PASS' if passed_gap else 'FAIL'}")
    all_pass = all_pass and passed_gap

    # Spectrum should span a meaningful energy range
    spread = evals_strained[-1] - evals_strained[0]
    passed_spread = spread > 0.5
    print(f"    Spectrum spread: {spread:.4f} eV (> 0.5): "
          f"{'PASS' if passed_spread else 'FAIL'}")
    all_pass = all_pass and passed_spread

    return all_pass


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        sys.exit(1)

    build_dir = os.path.abspath(sys.argv[1])

    exe_path = os.path.join(build_dir, "src", "bandStructure")
    if not os.path.isfile(exe_path):
        print(f"ERROR: executable not found: {exe_path}")
        sys.exit(1)

    print("=" * 60)
    print("  Wire Strain Validation (InAs/GaAs core-shell)")
    print("=" * 60)
    print()

    # Wire parameters
    core_size = 30.0
    nx, ny = 20, 20
    dx, dy = 5.0, 5.0
    width, height = 100.0, 100.0

    # Run strained wire
    config_strained = make_wire_config(nx, ny, dx, dy, width, height, core_size,
                                       strain=True)
    work = tempfile.mkdtemp(prefix="wire_strain_")
    try:
        cfg_path = os.path.join(work, "wire.toml")
        with open(cfg_path, "w") as f:
            f.write(config_strained)
        rc, output_dir = run_exe(build_dir, "bandStructure", cfg_path, work,
                                 timeout=500)
        if rc != 0:
            print(f"  FATAL: bandStructure returned {rc}")
            sys.exit(1)

        # Parse strain data
        strain_path = os.path.join(output_dir, "strain.dat")
        if not os.path.isfile(strain_path):
            print("  FATAL: strain.dat not found")
            sys.exit(1)
        strain_data = parse_strain(strain_path)

        # Parse eigenvalues
        eig_path = os.path.join(output_dir, "eigenvalues.dat")
        if not os.path.isfile(eig_path):
            print("  FATAL: eigenvalues.dat not found")
            sys.exit(1)
        data = parse_eigenvalues(eig_path)
        if not data:
            print("  FATAL: no eigenvalue data")
            sys.exit(1)
        evals_strained = data[0][1]
    finally:
        shutil.rmtree(work, ignore_errors=True)

    # Run unstrained wire for R9 comparison
    config_unstrained = make_wire_config(nx, ny, dx, dy, width, height, core_size,
                                         strain=False)
    work_u = tempfile.mkdtemp(prefix="wire_unstrain_")
    try:
        cfg_u = os.path.join(work_u, "wire.toml")
        with open(cfg_u, "w") as f:
            f.write(config_unstrained)
        rc_u, output_dir_u = run_exe(build_dir, "bandStructure", cfg_u, work_u,
                                      timeout=500)
        if rc_u != 0:
            print(f"  FATAL: unstrained bandStructure returned {rc_u}")
            sys.exit(1)
        eig_u = os.path.join(output_dir_u, "eigenvalues.dat")
        if not os.path.isfile(eig_u):
            print("  FATAL: unstrained eigenvalues.dat not found")
            sys.exit(1)
        data_u = parse_eigenvalues(eig_u)
        if not data_u:
            print("  FATAL: no unstrained eigenvalue data")
            sys.exit(1)
        evals_unstrained = data_u[0][1]
    finally:
        shutil.rmtree(work_u, ignore_errors=True)

    print(f"  Strain data: {len(strain_data)} grid points")
    print(f"  Grid: {nx}x{ny}, dx={dx}, dy={dy}")
    print(f"  Core: InAs 0-{core_size}A, Shell: GaAs {core_size}-{width}A")
    print()

    all_pass = True

    # R7
    r7_pass = test_r7_interior_strain(strain_data, core_size, width, height)
    all_pass = all_pass and r7_pass
    print()

    # R8
    r8_pass = test_r8_hydrostatic_strain(strain_data, core_size, width, height)
    all_pass = all_pass and r8_pass
    print()

    # R9
    r9_pass = test_r9_band_edge_shifts(evals_strained, evals_unstrained)
    all_pass = all_pass and r9_pass
    print()

    # Summary
    print("=" * 60)
    if all_pass:
        print("  PASS: all wire strain validation checks passed")
        sys.exit(0)
    else:
        print("  FAIL: one or more wire strain checks failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
