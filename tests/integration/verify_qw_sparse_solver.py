#!/usr/bin/env python3
"""Verify QW solver mode dispatch: INDEX (default) vs INDEX with explicit [solver].

Builds two QW configs — one without [solver] (smart default: DENSE INDEX),
one with explicit [solver] method="DENSE" mode="INDEX" — runs bandStructure
for each, and confirms eigenvalues match exactly.

This validates that the [solver] TOML section correctly configures the
polymorphic eigensolver and produces identical results to the smart defaults.

Usage: verify_qw_sparse_solver.py <build_dir> <source_dir>

# COVERAGE: observable=qw_solver_dispatch_index geometry=qw material=GaAs/AlGaAs
"""

import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from star_helpers import run_exe, parse_eigenvalues


# ── QW config WITHOUT [solver] (smart default: DENSE INDEX) ─────────

QW_DEFAULT_CONFIG = """\
confinement = "qw"
FDorder = 2
fd_step = 51

[wave_vector]
mode = "k0"
max = 0
nsteps = 1

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

# ── QW config WITH explicit [solver] DENSE INDEX ────────────────────

QW_EXPLICIT_CONFIG = """\
confinement = "qw"
FDorder = 2
fd_step = 51

[wave_vector]
mode = "k0"
max = 0
nsteps = 1

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

[solver]
method = "DENSE"
mode = "INDEX"
"""

# ── QW multi-kpoint FEAST sweep config (exercises fast path) ─────────
# nsteps=3 means k-points 2,3 take the value-only fast path
# (qw_ws%initialized == .true. after k=1). emin=emax=0 selects the AUTO
# sweep-envelope window (Gershgorin bounds at the two sweep endpoints,
# unioned) so the QW band-structure path extracts the full-spectrum
# gap-straddling [il,iu] bands via reconcile_band_slice. A narrow user
# window can't honor global-index extraction (review #4); the envelope is
# also safe from the contour-edge-on-band FEAST divide-by-zero.

QW_FASTPATH_CONFIG = """\
confinement = "qw"
FDorder = 2
fd_step = 41

[wave_vector]
mode = "kx"
max = 0.1
nsteps = 3

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

[solver]
method = "FEAST"
mode = "ENERGY"
emin = 0.0
emax = 0.0
m0 = 328
"""


def verify_qw_fastpath(build_dir, source_dir):
    """Run a 3-kpoint QW FEAST sweep; k-points 2,3 take the value-only
    fast path (qw_ws%initialized). Assert the eigenvalue count is stable
    across all k-points and no NaN/Inf appears."""
    import math
    print("\n" + "=" * 60)
    print("QW FEAST fast-path: multi-kpoint sweep stability")
    print("=" * 60)
    with tempfile.TemporaryDirectory(prefix="qw_fastpath_") as work:
        cfg_path = os.path.join(work, "cfg.toml")
        with open(cfg_path, "w") as f:
            f.write(QW_FASTPATH_CONFIG)
        rc, outdir = run_exe(build_dir, "bandStructure", cfg_path, work)
        if rc != 0:
            print("  FAIL: bandStructure exited nonzero (fast-path)")
            return False
        eig_path = os.path.join(outdir, "eigenvalues.dat")
        if not os.path.exists(eig_path):
            print("  FAIL: eigenvalues.dat not produced (fast-path)")
            return False
        rows = parse_eigenvalues(eig_path)
        if not rows:
            print("  FAIL: could not parse eigenvalues (fast-path)")
            return False
        # Every k-point must produce the same number of eigenvalues: the
        # fast path must not drop/gain entries vs the slow path.
        counts = {len(r[1]) for r in rows}
        if len(counts) != 1:
            print(f"  FAIL: inconsistent eigenvalue counts across k-points: {counts}")
            return False
        for ki, (_k, evals) in enumerate(rows):
            for v in evals:
                if math.isnan(v) or math.isinf(v):
                    print(f"  FAIL: NaN/Inf eigenvalue at k-point {ki}")
                    return False
        print(f"  PASS: QW FEAST fast-path sweep stable ({len(rows)} k-points, "
              f"{counts.pop()} eigenvalues each)")
        return True


def verify_qw_fastpath_vs_dense(build_dir, source_dir):
    """Run the same QW sweep with method=FEAST (fast path active on
    k-points 2,3) and method=DENSE (no CSR fast path), then assert each
    k-point's retained eigenvalues match to 1e-8. This is the rigorous
    two-paths-same-result equivalence gate."""
    import math
    print("\n" + "=" * 60)
    print("QW fast-path equivalence: FEAST (fast) vs DENSE reference")
    print("=" * 60)
    with tempfile.TemporaryDirectory(prefix="qw_equiv_") as work:
        # FEAST run (uses the fast path for k >= 2)
        cfgf = os.path.join(work, "feast.toml")
        with open(cfgf, "w") as f:
            f.write(QW_FASTPATH_CONFIG)
        rc, ofeast = run_exe(build_dir, "bandStructure", cfgf, work)
        if rc != 0:
            print("  FAIL: FEAST run exited nonzero")
            return False
        # DENSE reference (same energy window; no CSR fast path)
        dcfg = QW_FASTPATH_CONFIG.replace('method = "FEAST"', 'method = "DENSE"')
        cfgd = os.path.join(work, "dense.toml")
        with open(cfgd, "w") as f:
            f.write(dcfg)
        rc, odense = run_exe(build_dir, "bandStructure", cfgd, work)
        if rc != 0:
            print("  FAIL: DENSE reference run exited nonzero")
            return False
        feast_rows = parse_eigenvalues(os.path.join(ofeast, "eigenvalues.dat"))
        dense_rows = parse_eigenvalues(os.path.join(odense, "eigenvalues.dat"))
        if not feast_rows or not dense_rows:
            print("  FAIL: could not parse eigenvalues")
            return False
        if len(feast_rows) != len(dense_rows):
            print(f"  FAIL: k-point count mismatch {len(feast_rows)} vs {len(dense_rows)}")
            return False
        worst = 0.0
        for ki, (fr, dr) in enumerate(zip(feast_rows, dense_rows)):
            fe, de = fr[1], dr[1]
            m = min(len(fe), len(de))
            if m == 0:
                print(f"  FAIL: k={ki} produced zero eigenvalues")
                return False
            for a, b in zip(fe[:m], de[:m]):
                d = abs(a - b)
                if d > worst:
                    worst = d
                if d > 1.0e-8:
                    print(f"  FAIL: k={ki} eigenvalue mismatch {a:.10e} vs {b:.10e} "
                          f"(diff={d:.2e})")
                    return False
        print(f"  PASS: QW fast-path FEAST matches DENSE reference to 1e-8 "
              f"(worst diff {worst:.2e})")
        return True


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        sys.exit(1)

    build_dir = os.path.abspath(sys.argv[1])
    source_dir = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else os.getcwd()
    exe = os.path.join(build_dir, "src", "bandStructure")
    if not os.path.isfile(exe):
        print(f"SKIP: bandStructure not found at {exe}")
        sys.exit(0)

    print("=" * 60)
    print("QW solver dispatch: default vs explicit [solver] DENSE INDEX")
    print("=" * 60)

    failures = []

    # ── Run with smart default (no [solver]) ──────────────────────
    print("\n  Running QW without [solver] (smart default: DENSE INDEX)...")
    with tempfile.TemporaryDirectory(prefix="qw_default_") as work:
        cfg_path = os.path.join(work, "config.toml")
        with open(cfg_path, "w") as f:
            f.write(QW_DEFAULT_CONFIG)

        rc, outdir = run_exe(build_dir, "bandStructure", cfg_path, work)
        if rc != 0:
            print(f"  FAIL: default run returned {rc}")
            sys.exit(1)

        eig_path = os.path.join(outdir, "eigenvalues.dat")
        if not os.path.isfile(eig_path):
            print("  FAIL: eigenvalues.dat not produced (default)")
            sys.exit(1)

        default_rows = parse_eigenvalues(eig_path)

    if not default_rows:
        print("  FAIL: could not parse default eigenvalues")
        sys.exit(1)

    default_k0, default_evals = default_rows[0]
    print(f"    Eigenvalues at k=0: {len(default_evals)}")

    # ── Run with explicit [solver] ────────────────────────────────
    print("  Running QW with explicit [solver] DENSE INDEX...")
    with tempfile.TemporaryDirectory(prefix="qw_explicit_") as work:
        cfg_path = os.path.join(work, "config.toml")
        with open(cfg_path, "w") as f:
            f.write(QW_EXPLICIT_CONFIG)

        rc, outdir = run_exe(build_dir, "bandStructure", cfg_path, work)
        if rc != 0:
            print(f"  FAIL: explicit run returned {rc}")
            sys.exit(1)

        eig_path = os.path.join(outdir, "eigenvalues.dat")
        if not os.path.isfile(eig_path):
            print("  FAIL: eigenvalues.dat not produced (explicit)")
            sys.exit(1)

        explicit_rows = parse_eigenvalues(eig_path)

    if not explicit_rows:
        print("  FAIL: could not parse explicit eigenvalues")
        sys.exit(1)

    explicit_k0, explicit_evals = explicit_rows[0]
    print(f"    Eigenvalues at k=0: {len(explicit_evals)}")

    # ── Compare eigenvalues ───────────────────────────────────────
    if len(default_evals) != len(explicit_evals):
        failures.append(f"eigenvalue count mismatch: default={len(default_evals)}, "
                       f"explicit={len(explicit_evals)}")
    else:
        max_diff = max(abs(d - e) for d, e in zip(default_evals, explicit_evals))
        print(f"    Max absolute difference: {max_diff:.2e}")

        if max_diff > 1e-12:
            failures.append(f"eigenvalue mismatch: max_diff={max_diff:.2e}")
        else:
            print("    PASS: eigenvalues identical")

    # ── Fast-path checks (multi-kpoint FEAST sweep) ───────────────
    if not verify_qw_fastpath(build_dir, source_dir):
        sys.exit(1)
    if not verify_qw_fastpath_vs_dense(build_dir, source_dir):
        sys.exit(1)

    # ── Summary ────────────────────────────────────────────────────
    print()
    if failures:
        print(f"FAIL: {len(failures)} issue(s):")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    else:
        print("PASS: QW solver dispatch produces identical eigenvalues")
        sys.exit(0)


if __name__ == "__main__":
    main()
