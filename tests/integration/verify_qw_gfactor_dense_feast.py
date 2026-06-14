#!/usr/bin/env python3
"""Verify QW g-factor: FEAST+ENERGY (full-spectrum window) vs DENSE+FULL.

This is the FEAST-enable regression gate for issue #08 (ADR 0005).

The QW g-factor is computed via Lowdin (quasi-degenerate) partitioning,
which projects out the non-degenerate subspace and therefore needs the
COMPLETE spectrum (all 8N eigenstates). FEAST is a partial-spectrum
contour solver: if its energy window only covers a narrow band, the
Lowdin projection silently uses a truncated subspace and the g-factor
is wrong.

Before issue #08 this combination was rejected up front by a stopgap
guard in validate_semantic (defs.f90). With the dispersion-aware window
authority (issue #03) and ENERGY mode (the Gamma-point solver no longer
forced to a mode that truncates under FEAST), QW g-factor runs under
FEAST+ENERGY with the window authority's Gershgorin-envelope bound set
to cover the ENTIRE spectral range, so FEAST returns all 8N eigenvalues
and Lowdin sees the full spectrum.

This script runs the gfactorCalculation executable on the SAME QW
config twice - once method=DENSE, once method=FEAST - and asserts:

  1. The g-factor tensor (gx, gy, gz) matches within 1e-6 between the
     two backends. A truncated FEAST spectrum would move the Lowdin
     result, so this guards against silent truncation.
  2. BOTH backends return the FULL 8N spectrum. We count the
     eigenfunction files written to output/ (eigenfunctions_k_00001_ev_*.dat)
     and assert the count == 8 * fd_step under both backends. If FEAST
     truncated the spectrum it would write fewer files.

The config is a GaSbW/InAsW/AlSbW QW (the same Winkler-parameter system
used by the regression_gfactor_cb golden test) with a small grid
(fd_step = 21 -> N = 168) so the run is fast while still exercising a
non-trivial 8-band QW g-factor. The FEAST run uses the AUTO window
(emin = emax = 0) so the window authority's Gershgorin-envelope bound
is exercised end-to-end through the Gamma-point path.

Usage: verify_qw_gfactor_dense_feast.py <build_dir> <source_dir>

# COVERAGE: observable=g*_cb geometry=qw material=GaSbW/InAsW backend=FEAST-vs-DENSE
"""

import glob
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))

from star_helpers import run_exe, parse_gfactor

TOL = 1.0e-6  # g-factor relative tolerance for FEAST-vs-DENSE agreement

FD_STEP = 21
N_TOTAL = 8 * FD_STEP  # full spectrum count the Gamma solver must return


# Shared QW config body (GaSbW/InAsW/AlSbW, Winkler parameters). The only
# field that varies between the two runs is [solver].method (+ mode/window).
# Both runs are at k=0 (mode = "k0", nsteps = 0), the Gamma-point g-factor.
QW_BODY = """\
confinement = "qw"
FDorder = 2
fd_step = {fd_step}
which_band = 0
band_idx = 1

[wave_vector]
mode = "k0"
max = 0.1
nsteps = 0

[bands]
num_cb = 8
num_vb = 24

[[material]]
name = "AlSbW"
z_min = -250
z_max = 250

[[material]]
name = "GaSbW"
z_min = -135
z_max = 135

[[material]]
name = "InAsW"
z_min = -35
z_max = 35
"""

# DENSE backend: FULL mode (LAPACK zheev returns all eigenvalues).
QW_DENSE = QW_BODY.format(fd_step=FD_STEP) + """\
[solver]
method = "DENSE"
mode = "FULL"
"""

# FEAST backend: ENERGY mode with AUTO window (emin = emax = 0). The
# Gamma-point solver routes this through apply_solver_window's single-H
# Gershgorin-envelope bound, which covers the entire spectrum so FEAST
# returns all 8N eigenvalues. m0 = N_total lets FEAST size its subspace
# to hold the full spectrum (the info=3 retry loop also grows M0 to N).
QW_FEAST = QW_BODY.format(fd_step=FD_STEP) + """\
[solver]
method = "FEAST"
mode = "ENERGY"
emin = 0.0
emax = 0.0
m0 = {m0}
""".format(m0=N_TOTAL)


def _run_gfactor(build_dir, work, config_text, label):
    """Write config, run gfactorCalculation, return (g, n_eig) or None.

    g = (gx, gy, gz) parsed from output/gfactor.dat.
    n_eig = count of eigenfunctions_k_00001_ev_*.dat files (full spectrum).
    """
    cfg_path = os.path.join(work, f"{label}.toml")
    with open(cfg_path, "w") as f:
        f.write(config_text)
    rc, outdir = run_exe(build_dir, "gfactorCalculation", cfg_path, work,
                         timeout=300)
    if rc != 0:
        print(f"  FAIL: gfactorCalculation exited {rc} ({label})")
        return None
    gf_path = os.path.join(outdir, "gfactor.dat")
    if not os.path.isfile(gf_path):
        print(f"  FAIL: gfactor.dat not produced ({label})")
        return None
    g = parse_gfactor(gf_path)
    eig_files = glob.glob(os.path.join(
        outdir, "eigenfunctions_k_00001_ev_*.dat"))
    return g, len(eig_files)


def _compare_gfactors(g_dense, g_feast):
    """Compare (gx, gy, gz) between backends; return True if all match."""
    labels = ("gx", "gy", "gz")
    worst = 0.0
    for i, name in enumerate(labels):
        a, b = g_dense[i], g_feast[i]
        # Relative tolerance with absolute floor for near-zero values.
        atol = 1.0e-12
        denom = max(abs(a), atol)
        err = abs(a - b) / denom
        if err > worst:
            worst = err
        if err > TOL:
            print(f"  FAIL: {name} DENSE={a:.10e} FEAST={b:.10e} "
                  f"rel_err={err:.2e} (tol {TOL:.0e})")
            return False
        print(f"  PASS: {name} DENSE={a:.10e} FEAST={b:.10e} "
              f"rel_err={err:.2e}")
    print(f"  g-factor worst relative diff across (gx,gy,gz): {worst:.2e}")
    return True


def verify_qw_gfactor_dense_vs_feast(build_dir, source_dir):
    print("\n" + "=" * 64)
    print("QW g-factor: FEAST+ENERGY (full-spectrum window) vs DENSE+FULL")
    print(f"  config: GaSbW/InAsW/AlSbW QW, fd_step={FD_STEP}, N={N_TOTAL}")
    print("=" * 64)

    with tempfile.TemporaryDirectory(prefix="qw_gf_ds_") as work:
        print("\n  [1/2] DENSE+FULL baseline run...")
        dense = _run_gfactor(build_dir, work, QW_DENSE, "dense")
        if dense is None:
            return False
        g_dense, n_dense = dense
        print(f"        g = {g_dense}, eigenfunctions written = {n_dense}")

        print("\n  [2/2] FEAST+ENERGY (AUTO window) run...")
        feast = _run_gfactor(build_dir, work, QW_FEAST, "feast")
        if feast is None:
            return False
        g_feast, n_feast = feast
        print(f"        g = {g_feast}, eigenfunctions written = {n_feast}")

        # --- Full-spectrum assertion (8N) under BOTH backends ------------
        # The Lowdin projection needs the complete spectrum. If FEAST's
        # window did not cover the full spectral range it would return
        # fewer than 8N eigenvalues and write fewer eigenfunction files.
        print(f"\n  Full-spectrum count check (expected {N_TOTAL} = 8*fd_step):")
        ok = True
        if n_dense != N_TOTAL:
            print(f"  FAIL: DENSE wrote {n_dense} eigenfunctions, "
                  f"expected {N_TOTAL}")
            ok = False
        else:
            print(f"  PASS: DENSE wrote {n_dense} eigenfunctions "
                  f"(== 8*fd_step)")
        if n_feast != N_TOTAL:
            print(f"  FAIL: FEAST wrote {n_feast} eigenfunctions, "
                  f"expected {N_TOTAL} -- spectrum TRUNCATED")
            ok = False
        else:
            print(f"  PASS: FEAST wrote {n_feast} eigenfunctions "
                  f"(== 8*fd_step, full spectrum)")

        # --- g-factor agreement -----------------------------------------
        print(f"\n  g-factor agreement (tol {TOL:.0e} relative):")
        if not _compare_gfactors(g_dense, g_feast):
            ok = False

        return ok


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        sys.exit(2)
    build_dir = os.path.abspath(sys.argv[1])
    source_dir = os.path.abspath(sys.argv[2])

    exe = os.path.join(build_dir, "src", "gfactorCalculation")
    if not os.path.isfile(exe):
        print(f"SKIP: gfactorCalculation not found at {exe}")
        sys.exit(0)

    if not verify_qw_gfactor_dense_vs_feast(build_dir, source_dir):
        print("\nFAIL: QW g-factor FEAST-vs-DENSE mismatch or spectrum "
              "truncated (issue #08)")
        sys.exit(1)

    print("\nPASS: QW g-factor under FEAST+ENERGY matches DENSE+FULL and "
          "returns the full 8N spectrum (issue #08)")
    sys.exit(0)


if __name__ == "__main__":
    main()
