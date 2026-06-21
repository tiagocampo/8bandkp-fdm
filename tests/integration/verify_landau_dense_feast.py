#!/usr/bin/env python3
"""Verify Landau k-solve: FEAST+ENERGY (sweep-envelope window) vs DENSE+INDEX.

This is the FEAST-enable regression gate for issue #10 (ADR 0005).

Before issue #10 the Landau k-sweep block in main.f90 HARDCODED the dense
backend path and ASSUMED INDEX-sized results:

    M = landau_result%nev_found       # INDEX: == iuu-il+1 (safe)
    eig(1:M, k) = landau_result%eigenvalues(1:M)
    eigv(:, 1:M, k) = landau_result%eigenvectors(:, 1:M)

so requesting method=FEAST on Landau either silently used a default
[-1,+1] eV FEAST window (far too narrow - the Landau Hamiltonian at finite
B has dense LL structure spanning a wide band) AND overran the `eig` /
`eigv` arrays (sized to iuu-il+1) when FEAST returned more in-window
states. The #01 implementer flagged both modes segfaulting.

Issue #10 routes the Landau k-sweep through the window authority's
dispersion-aware envelope variant (asw_envelope, issue #03): the
Gershgorin bounds at the sweep's two endpoints (k=0 and k=max) are
unioned into ONE stable window per sweep - the same shape as the
band-structure k-sweep in main.f90 and the QW optics k_par sweep in
main_optics.f90 (#09). It also fixes the result-handling to be
mode-robust: under ENERGY/FEAST it extracts the global-index [il, iuu]
slice from FEAST's full in-window sorted spectrum (the same
reconciliation #09 did for QW optics); under INDEX/DENSE it stays
byte-identical to the prior code.

This script runs the bandStructure executable on the SAME Landau config
twice - once method=DENSE (AUTO->DENSE+INDEX, the golden default), once
method=FEAST+ENERGY (AUTO sweep-envelope window) - and asserts:

  1. PHYSICS GATE: the eigenvalues.dat output matches between the two
     backends within 1e-6 relative. Divergent eigenvalues would mean the
     FEAST window truncated the spectrum, the result-handling extracted
     the wrong slice, or the band indexing shifted between INDEX and
     ENERGY mode.
  2. DISPATCH GATE: the FEAST run's stdout contains the Landau-loop
     backend diagnostic. At HEAD (pre-#10) the Landau loop printed no
     backend-name diagnostic, so the FEAST-tagged run could silently run
     a broken path; this gate catches that. (Bit-identity of eigenvalues
     is NOT a reliable signal: on small problems FEAST and LAPACK can
     agree to full double precision, making the eigenvalues bit-identical
     even when genuinely different backends ran - so we require BOTH the
     physics agreement AND the dispatch diagnostic.)

The config is the golden landau_InAs.toml shape (InAs, nx=100, Bz=5T,
ky-sweep nsteps=10) - a genuine in-plane k-dispersion sweep at fixed B
(NOT the B-sweep; the B-sweep / fan diagram stays application-layer per
ADR 0003). The FEAST run uses the AUTO window (emin = emax = 0) so the
envelope authority is exercised end-to-end.

Usage: verify_landau_dense_feast.py <build_dir> <source_dir>

# COVERAGE: observable=landau_levels geometry=bulk material=InAs backend=FEAST-vs-DENSE
"""

import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))

TOL = 1.0e-6  # eigenvalue relative tolerance for FEAST-vs-DENSE agreement

NX = 100
N_TOTAL = 8 * NX  # full matrix dimension

# Shared Landau config body (InAs, the golden landau_InAs.toml material
# system). The only field that varies between the two runs is
# [solver].method (+ mode/window). Real ky-sweep so the FEAST envelope
# over the two sweep endpoints is exercised. This is an in-plane
# k-dispersion sweep at FIXED B (the B-sweep fan diagram is a separate
# application-layer block, untouched by issue #10 per ADR 0003).
LANDAU_BODY = """\
confinement = "landau"
FDorder = 2
fd_step = {fd_step}

[wave_vector]
mode = "ky"
max = 0.05
nsteps = 10

[bands]
num_cb = 4
num_vb = 4

[[material]]
name = "InAs"

[landau]
nx = {nx}
width = 2000.0

[b_field]
components = [0.0, 0.0, 5.0]
"""

# DENSE backend: AUTO resolves to DENSE+INDEX for Landau
# (resolve_solver_defaults, ADR 0004). This is the golden default, so
# this run is the baseline.
LANDAU_DENSE = LANDAU_BODY.format(fd_step=NX, nx=NX) + """\
[solver]
method = "DENSE"
"""

# FEAST backend: AUTO window (emin = emax = 0). The Landau k-sweep routes
# this through apply_solver_window's dispersion-aware envelope (union of
# Gershgorin bounds at k=0 and k=max), giving ONE stable window per
# sweep. m0 large enough to hold the numcb+numvb=8 retained bands (the
# info=3 retry loop also grows M0 if needed).
LANDAU_FEAST = LANDAU_BODY.format(fd_step=NX, nx=NX) + """\
[solver]
method = "FEAST"
emin = 0.0
emax = 0.0
m0 = {m0}
""".format(m0=N_TOTAL)


def _run_landau(build_dir, work, config_text, label):
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


def _load_eigenvalues(path):
    """Load eigenvalues.dat as (k_values, eig_array).

    eigenvalues.dat rows: |k| eig1 eig2 ... eigN (space-separated, '#'-comment header).
    Returns (k_vec (nsteps,), eig (nsteps, nev)).
    """
    if not os.path.isfile(path):
        return None
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    k_vec = data[:, 0]
    eigs = data[:, 1:]
    return k_vec, eigs


def _compare_eigenvalues(dense, feast):
    """Compare eigenvalues between backends; return (ok, worst_rel_err)."""
    kd, ed = dense
    kf, ef = feast
    if ed.shape != ef.shape:
        print(f"  FAIL: eigenvalue shape mismatch DENSE={ed.shape} "
              f"FEAST={ef.shape}")
        return False, float("inf")
    # k-grids must agree (same mode/max/nsteps).
    k_err = float(np.max(np.abs(kd - kf)))
    if k_err > 1.0e-10:
        print(f"  FAIL: k-grid shifted (max diff {k_err:.2e})")
        return False, float("inf")
    # Relative tolerance with absolute floor for near-zero eigenvalues.
    denom = np.maximum(np.abs(ed), 1.0e-10)
    rel = np.abs(ed - ef) / denom
    worst = float(np.max(rel))
    if worst > TOL:
        # Locate the worst eigenvalue for diagnostics.
        idx = np.unravel_index(np.argmax(rel), rel.shape)
        print(f"  FAIL: worst rel diff {worst:.2e} (tol {TOL:.0e}) at "
              f"(k_idx={idx[0]}, band={idx[1]}) "
              f"[DENSE={ed[idx]:.6e}, FEAST={ef[idx]:.6e}]")
        return False, worst
    print(f"  PASS: worst rel diff {worst:.2e} across {ed.size} eigenvalues "
          f"({ed.shape[0]} k-points x {ed.shape[1]} bands)")
    return True, worst


def verify_landau_dense_vs_feast(build_dir, source_dir):
    print("\n" + "=" * 64)
    print("Landau k-solve: FEAST+ENERGY (sweep-envelope window) vs DENSE+INDEX")
    print(f"  config: InAs, nx={NX}, N={N_TOTAL}, Bz=5T, ky-sweep nsteps=10")
    print("=" * 64)

    with tempfile.TemporaryDirectory(prefix="landau_ds_") as work:
        print("\n  [1/2] DENSE baseline run...")
        dense = _run_landau(build_dir, work, LANDAU_DENSE, "dense")
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
        feast = _run_landau(build_dir, work, LANDAU_FEAST, "feast")
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
        # The two backends must agree within tolerance: divergent
        # eigenvalues would mean the FEAST window truncated the spectrum,
        # the result-handling extracted the wrong slice, or the band
        # indexing shifted between INDEX and ENERGY mode.
        print(f"\n  Eigenvalue agreement (tol {TOL:.0e} relative):")
        ed = _load_eigenvalues(os.path.join(outdir_dense, "eigenvalues.dat"))
        ef = _load_eigenvalues(os.path.join(outdir_feast, "eigenvalues.dat"))
        if ed is None or ef is None:
            print("  FAIL: eigenvalues.dat missing under one or both backends "
                  f"(DENSE={'yes' if ed else 'no'}, FEAST={'yes' if ef else 'no'})")
            ok = False
        else:
            physics_ok, worst = _compare_eigenvalues(ed, ef)
            if not physics_ok:
                ok = False

        # --- DISPATCH GATE -----------------------------------------------
        # The FEAST run must actually have taken the FEAST path in the
        # LANDAU K-SWEEP (not just in setup_init, which always derives
        # the solver config). The post-#10 Landau loop prints a
        # backend-name diagnostic naming the resolved method. At HEAD
        # (pre-#10) the Landau loop printed no such line, so the
        # FEAST-tagged run could silently run a broken path. This gate
        # catches that override by requiring the Landau-loop-specific
        # FEAST diagnostic in the FEAST run's stdout.
        print("\n  Dispatch gate (FEAST path actually taken in Landau loop):")
        landau_feast_diag = "Landau k-sweep (FEAST"
        if landau_feast_diag in stdout_feast:
            print(f"  PASS: Landau-loop FEAST diagnostic present in stdout "
                  f"(FEAST path honored)")
        else:
            print(f"  FAIL: Landau-loop FEAST diagnostic absent from "
                  f"stdout -- method=FEAST was not honored in the Landau "
                  f"k-sweep (issue #10 not fixed)")
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

    if not verify_landau_dense_vs_feast(build_dir, source_dir):
        print("\nFAIL: Landau FEAST-vs-DENSE mismatch or FEAST path not "
              "taken (issue #10)")
        sys.exit(1)

    print("\nPASS: Landau k-solve under FEAST+ENERGY matches DENSE+INDEX and "
          "the FEAST path is honored (issue #10)")
    sys.exit(0)


if __name__ == "__main__":
    main()
