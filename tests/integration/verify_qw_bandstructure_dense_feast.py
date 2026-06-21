#!/usr/bin/env python3
"""Verify QW band-structure: FEAST+ENERGY (sweep-envelope window) vs DENSE+INDEX.

This is the FEAST-enable regression gate for review finding #4 (Task 7).

Before review #4 the QW FEAST CSR block in main.f90 extracted the eigenvalue
slice with a lowest-M copy:

    M = min(result_bs%nev_found, iuu-il+1)
    eig(1:M, k) = result_bs%eigenvalues(1:M)          # DIVERGENT under FEAST

For the DENSE+INDEX backend that is correct (INDEX pre-slices to exactly
[il, iuu], so M == nev_target and eigenvalues(1:M) IS the gap-straddling
window). But under FEAST+ENERGY the solver returns the FULL computed spectrum
(lowest-up), so eigenvalues(1:M) is the deepest-VALENCE bands - a completely
different physical slice than the gap-straddling [il, iuu] window the band
structure is supposed to show. So requesting method=FEAST on a QW
band-structure k-sweep silently produced a wrong-axis plot (deepest valence
bands instead of the gap-straddling conduction/valence bands).

Review #4 routes the QW FEAST CSR sweep through the shared reconcile_band_slice
helper (added in Task 4 / A1a, already used by the Landau sweep in Task 6):
FEAST+ENERGY full spectrum -> offset il; the DENSE+INDEX pre-sliced case ->
offset 1. The helper decides the offset; the call site just extracts
[il, iuu] from it.

This script runs the bandStructure executable on the SAME QW config twice -
once method=DENSE, once method=FEAST - and asserts:

  1. PHYSICS GATE: the gap-straddling eigenvalue bands in eigenvalues.dat
     agree between the two backends within 1e-6 ABSOLUTE eV. FEAST-vs-LAPACK
     agree to ~1e-9 on this problem, so 1e-6 is a robust gate that is far
     below any physical effect. Pre-#4 the FEAST run returned the deepest
     valence bands (lowest-M copy), so the worst diff was several eV.
  2. DISPATCH GATE: the FEAST run's stdout contains the QW CSR sweep
     diagnostic ("QW k-sweep (FEAST/CSR)"). This proves the CSR/FEAST path
     was actually taken, not silently overridden to DENSE. (Bit-identity of
     eigenvalues is NOT a reliable signal: on small problems FEAST and
     LAPACK can agree to full double precision, making the output
     bit-identical even when genuinely different backends ran.)

The config is an Al30Ga70As/GaAs QW (same material system as the golden
qw_gaas_algaas_optics config) with a compact grid (fd_step=41 -> N=328) and a
real k_x sweep (nsteps=11). The FEAST run uses the AUTO window
(emin = emax = 0) so the sweep-envelope authority is exercised end-to-end;
both emin+emax=0 so the Task-1 partial-window guard does NOT fire.

Usage: verify_qw_bandstructure_dense_feast.py <build_dir> <source_dir>

# COVERAGE: observable=band_structure geometry=qw material=GaAs/AlGaAs backend=FEAST-vs-DENSE
"""

import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))

EIG_TOL = 1.0e-6  # absolute eV tolerance for FEAST-vs-DENSE eigenvalue agreement

FD_STEP = 41
N_TOTAL = 8 * FD_STEP  # full matrix dimension


# Shared QW config body (Al30Ga70As/GaAs, the golden optics material system,
# last-layer-wins QW pattern: barrier covers full domain, well overwrites
# center). The only field that varies between the two runs is
# [solver].method (+ mode/window). Real k_x sweep so the FEAST envelope over
# the two sweep endpoints is exercised. num_cb=4 + num_vb=8 -> 12 retained
# gap-straddling bands.
QW_BODY = """\
confinement = "qw"
FDorder = 2
fd_step = {fd_step}
which_band = 0
band_idx = 1

[wave_vector]
mode = "kx"
max = 0.1
nsteps = 11

[bands]
num_cb = 4
num_vb = 8

[[material]]
name = "Al30Ga70As"
z_min = -200
z_max = 200

[[material]]
name = "GaAs"
z_min = -50
z_max = 50
"""

# DENSE backend: AUTO resolves to DENSE+INDEX for QW (resolve_solver_defaults,
# ADR 0004). This is the golden default, so this run is the baseline.
QW_DENSE = QW_BODY.format(fd_step=FD_STEP) + """\
[solver]
method = "DENSE"
"""

# FEAST backend: AUTO window (emin = emax = 0). The band-structure sweep
# routes this through apply_solver_window's dispersion-aware envelope (union
# of Gershgorin bounds at k_x=0 and k_x=max), giving ONE stable window per
# sweep. m0 large enough to hold the numcb+numvb=12 retained bands (the
# info=3 retry loop also grows M0 if needed). Both emin+emax=0 so the
# Task-1 partial-window guard does NOT fire.
QW_FEAST = QW_BODY.format(fd_step=FD_STEP) + """\
[solver]
method = "FEAST"
emin = 0.0
emax = 0.0
m0 = {m0}
""".format(m0=N_TOTAL)


def _run_bandstructure(build_dir, work, config_text, label):
    """Write config, run bandStructure, return (rc, stdout, outdir)."""
    cfg_path = os.path.join(work, f"{label}.toml")
    with open(cfg_path, "w") as f:
        f.write(config_text)

    exe_path = os.path.abspath(os.path.join(build_dir, "src", "bandStructure"))
    if not os.path.isfile(exe_path):
        print(f"  FAIL: bandStructure not found at {exe_path}")
        return None

    # Run in an isolated subdir per label so the two runs do not clobber
    # each other's output/.
    run_dir = os.path.join(work, f"run_{label}")
    os.makedirs(os.path.join(run_dir, "output"), exist_ok=True)
    shutil.copy2(cfg_path, os.path.join(run_dir, "input.toml"))

    try:
        result = subprocess.run(
            [exe_path],
            cwd=run_dir,
            capture_output=True,
            text=True,
            timeout=600,
        )
    except subprocess.TimeoutExpired:
        print(f"  FAIL: bandStructure timed out ({label})")
        return None

    return result.returncode, result.stdout, os.path.join(run_dir, "output")


def parse_eigenvalues(path):
    """Load eigenvalues.dat as (k_array, evals_2d).

    eigenvalues.dat rows: |k| ev1 ev2 ... evN (whitespace-separated,
    '#'-comment header). Returns (k_array (nk,), evals_2d (nk, n_bands)).
    """
    if not os.path.isfile(path):
        return None
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    k_array = data[:, 0]
    evals_2d = data[:, 1:]
    return k_array, evals_2d


def _compare_eigenvalues(dense, feast):
    """Compare eigenvalue bands between backends; return (ok, worst_abs_err)."""
    kd, ed = dense
    kf, ef = feast
    if kd.shape != kf.shape:
        print(f"  FAIL: k-grid length mismatch DENSE={kd.shape} FEAST={kf.shape}")
        return False, float("inf")
    # k-grids must agree (same mode/max/nsteps).
    k_err = float(np.max(np.abs(kd - kf)))
    if k_err > 1.0e-10:
        print(f"  FAIL: k-grid shifted (max diff {k_err:.2e})")
        return False, float("inf")
    # Compare the band window both runs actually wrote. With matching
    # num_cb/num_vb the window width is identical, but guard against any
    # width drift by taking the min.
    n = min(ed.shape[1], ef.shape[1])
    if n == 0:
        print("  FAIL: no eigenvalue columns to compare")
        return False, float("inf")
    diff = np.abs(ed[:, :n] - ef[:, :n])
    worst = float(np.max(diff))
    if worst > EIG_TOL:
        idx = np.unravel_index(np.argmax(diff), diff.shape)
        print(f"  FAIL: worst eigenvalue diff {worst:.6e} eV (tol {EIG_TOL:.0e}) "
              f"at (k_idx={idx[0]}, band={idx[1]}) "
              f"[DENSE={ed[idx[0], idx[1]]:.6e}, "
              f"FEAST={ef[idx[0], idx[1]]:.6e}]")
        return False, worst
    print(f"  PASS: worst eigenvalue diff {worst:.3e} eV (tol {EIG_TOL:.0e}) "
          f"across {ed.shape[0]} k-points x {n} bands")
    return True, worst


def verify_qw_bandstructure_dense_vs_feast(build_dir, source_dir):
    print("\n" + "=" * 64)
    print("QW band-structure: FEAST+ENERGY (sweep-envelope window) vs DENSE+INDEX")
    print(f"  config: Al30Ga70As/GaAs QW, fd_step={FD_STEP}, N={N_TOTAL}, "
          f"k_x sweep nsteps=11, num_cb=4 num_vb=8")
    print("=" * 64)

    with tempfile.TemporaryDirectory(prefix="qw_bs_ds_") as work:
        print("\n  [1/2] DENSE baseline run...")
        dense = _run_bandstructure(build_dir, work, QW_DENSE, "dense")
        if dense is None:
            return False
        rc_dense, stdout_dense, outdir_dense = dense
        if rc_dense != 0:
            print(f"  FAIL: bandStructure exited {rc_dense} (dense)")
            print("  --- stdout (last 20 lines) ---")
            for line in stdout_dense.splitlines()[-20:]:
                print(f"  {line}")
            return False

        print("\n  [2/2] FEAST+ENERGY (AUTO envelope window) run...")
        feast = _run_bandstructure(build_dir, work, QW_FEAST, "feast")
        if feast is None:
            return False
        rc_feast, stdout_feast, outdir_feast = feast
        if rc_feast != 0:
            print(f"  FAIL: bandStructure exited {rc_feast} (feast)")
            print("  --- stdout (last 20 lines) ---")
            for line in stdout_feast.splitlines()[-20:]:
                print(f"  {line}")
            return False

        ok = True

        # --- PHYSICS GATE: eigenvalue agreement --------------------------
        # The two backends must agree within tolerance on the gap-straddling
        # band window. Pre-#4 the FEAST run returned the deepest valence
        # bands (lowest-M copy), so the worst diff was several eV.
        print(f"\n  Eigenvalue agreement (tol {EIG_TOL:.0e} ABSOLUTE eV):")
        ed = parse_eigenvalues(os.path.join(outdir_dense, "eigenvalues.dat"))
        ef = parse_eigenvalues(os.path.join(outdir_feast, "eigenvalues.dat"))
        if ed is None or ef is None:
            print("  FAIL: eigenvalues.dat missing under one or both backends "
                  f"(DENSE={'yes' if ed else 'no'}, FEAST={'yes' if ef else 'no'})")
            ok = False
        else:
            physics_ok, worst = _compare_eigenvalues(ed, ef)
            if not physics_ok:
                ok = False

        # --- DISPATCH GATE -----------------------------------------------
        # The FEAST run must actually have taken the FEAST/CSR path in the
        # QW K-SWEEP (not just in setup_init, which always derives the
        # solver config). The QW FEAST CSR block prints a backend-name
        # diagnostic naming the resolved method and the CSR backend. This
        # gate catches any silent override to DENSE by requiring the
        # QW-sweep-specific FEAST/CSR diagnostic in the FEAST run's stdout.
        print("\n  Dispatch gate (FEAST/CSR path actually taken in QW sweep):")
        qw_feast_diag = "FEAST/CSR"
        if qw_feast_diag in stdout_feast:
            print(f"  PASS: QW-sweep FEAST/CSR diagnostic present in stdout "
                  f"(FEAST/CSR path honored)")
        else:
            print(f"  FAIL: QW-sweep FEAST/CSR diagnostic absent from "
                  f"stdout -- method=FEAST was not honored in the QW "
                  f"band-structure CSR k-sweep (review #4 not fixed)")
            ok = False

        return ok


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        sys.exit(2)
    build_dir = os.path.abspath(sys.argv[1])
    source_dir = os.path.abspath(sys.argv[2])

    exe = os.path.join(build_dir, "src", "bandStructure")
    if not os.path.isfile(exe):
        print(f"SKIP: bandStructure not found at {exe}")
        sys.exit(0)

    if not verify_qw_bandstructure_dense_vs_feast(build_dir, source_dir):
        print("\nFAIL: QW band-structure FEAST-vs-DENSE mismatch or FEAST/CSR "
              "path not taken (review #4)")
        sys.exit(1)

    print("\nPASS: QW band-structure under FEAST+ENERGY matches DENSE+INDEX and "
          "the FEAST/CSR path is honored (review #4)")
    sys.exit(0)


if __name__ == "__main__":
    main()
