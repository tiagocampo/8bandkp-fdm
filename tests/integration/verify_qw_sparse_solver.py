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


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        sys.exit(1)

    build_dir = os.path.abspath(sys.argv[1])
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
