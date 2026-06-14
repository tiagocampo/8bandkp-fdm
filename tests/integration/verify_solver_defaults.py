#!/usr/bin/env python3
"""Verify that omitting [solver] produces correct smart defaults for all confinement modes.

Smart defaults (from simulation_setup.f90):
  - bulk:  DENSE + FULL
  - qw:    DENSE + INDEX
  - wire:  FEAST + ENERGY
  - landau: DENSE + INDEX

The test creates minimal TOML configs WITHOUT a [solver] section for each mode,
runs the bandStructure executable, and verifies eigenvalues are produced and
physically reasonable.

Usage: verify_solver_defaults.py <build_dir> <source_dir>

# COVERAGE: observable=eigensolver_defaults geometry=bulk material=GaAs
# COVERAGE: observable=eigensolver_defaults geometry=qw material=GaAs
# COVERAGE: observable=eigensolver_defaults geometry=wire material=GaAs
# COVERAGE: observable=eigensolver_defaults geometry=landau material=InAs
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from star_helpers import run_exe, parse_eigenvalues


# ── Inline TOML configs without [solver] section ──────────────────────

BULK_GAAS_K0 = """\
confinement = "bulk"
FDorder = 2
fd_step = 101

[wave_vector]
mode = "k0"
max = 0
nsteps = 1

[bands]
num_cb = 2
num_vb = 6

[[material]]
name = "GaAs"
"""

QW_GAAS_K0 = """\
confinement = "qw"
FDorder = 2
fd_step = 101

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

WIRE_GAAS_K0 = """\
confinement = "wire"
FDorder = 2
fd_step = 1
which_band = 0
band_idx = 1

[wave_vector]
mode = "k0"
max = 0
nsteps = 1

[bands]
num_cb = 4
num_vb = 8

[wire]
nx = 11
ny = 11
dx = 3.0
dy = 3.0

[wire.geometry]
shape = "rectangle"
width = 33.0
height = 33.0

[[region]]
material = "GaAs"
inner = 0.0
outer = 100.0
"""

LANDAU_INAS_K0 = """\
confinement = "landau"
FDorder = 2
fd_step = 100

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
nx = 100
width = 2000.0

[b_field]
components = [0.0, 0.0, 5.0]
"""


def run_config(build_dir, config_text, label):
    """Write config to temp file, run bandStructure, return parsed eigenvalues."""
    with tempfile.TemporaryDirectory(prefix=f"solver_defaults_{label}_") as work:
        cfg_path = os.path.join(work, "config.toml")
        with open(cfg_path, "w") as f:
            f.write(config_text)

        rc, outdir = run_exe(build_dir, "bandStructure", cfg_path, work)
        if rc != 0:
            print(f"  FAIL: {label} — bandStructure returned {rc}")
            return None

        eig_path = os.path.join(outdir, "eigenvalues.dat")
        if not os.path.isfile(eig_path):
            print(f"  FAIL: {label} — eigenvalues.dat not produced")
            return None

        rows = parse_eigenvalues(eig_path)
        if not rows:
            print(f"  FAIL: {label} — could not parse eigenvalues")
            return None

        return rows


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        sys.exit(1)

    build_dir = os.path.abspath(sys.argv[1])
    # source_dir not directly used but kept for consistent interface
    exe = os.path.join(build_dir, "src", "bandStructure")
    if not os.path.isfile(exe):
        print(f"SKIP: bandStructure not found at {exe}")
        sys.exit(0)

    print("=" * 60)
    print("Solver defaults: verify smart dispatch without [solver] section")
    print("=" * 60)

    failures = []

    # ── Bulk: DENSE + FULL ─────────────────────────────────────────
    print("\n  Bulk GaAs (expected: DENSE + FULL)...")
    rows = run_config(build_dir, BULK_GAAS_K0, "bulk")
    if rows is None:
        failures.append("bulk: execution or parsing failed")
    else:
        k0, evals = rows[0]
        n_evals = len(evals)
        print(f"    k=0 eigenvalues: {n_evals}")
        if n_evals != 8:
            failures.append(f"bulk: expected 8 eigenvalues (FULL mode), got {n_evals}")
        else:
            # Check GaAs gap is reasonable (~1.42 eV at Gamma)
            evals_sorted = sorted(evals)
            # Bands 1-6 = valence (top ~0 eV), 7-8 = conduction
            gap = evals_sorted[6] - evals_sorted[5]
            print(f"    Band gap: {gap:.4f} eV (expected ~1.42)")
            if abs(gap - 1.42) > 0.5:
                failures.append(f"bulk: gap {gap:.4f} eV too far from expected ~1.42 eV")
            else:
                print("    PASS")

    # ── QW: DENSE + INDEX ──────────────────────────────────────────
    print("\n  QW GaAs/AlGaAs (expected: DENSE + INDEX)...")
    rows = run_config(build_dir, QW_GAAS_K0, "qw")
    if rows is None:
        failures.append("qw: execution or parsing failed")
    else:
        k0, evals = rows[0]
        n_evals = len(evals)
        print(f"    k=0 eigenvalues: {n_evals}")
        # INDEX mode returns num_cb + num_vb = 4 + 8 = 12
        expected = 4 + 8
        if n_evals < expected:
            failures.append(f"qw: expected >= {expected} eigenvalues (INDEX mode), got {n_evals}")
        else:
            print("    PASS")

    # ── Wire: FEAST + ENERGY ───────────────────────────────────────
    print("\n  Wire GaAs 11x11 (expected: FEAST + ENERGY)...")
    rows = run_config(build_dir, WIRE_GAAS_K0, "wire")
    if rows is None:
        failures.append("wire: execution or parsing failed")
    else:
        k0, evals = rows[0]
        n_evals = len(evals)
        print(f"    k=0 eigenvalues: {n_evals}")
        # FEAST ENERGY mode: should find at least num_cb + num_vb eigenvalues
        expected = 4 + 8
        if n_evals < expected:
            failures.append(f"wire: expected >= {expected} eigenvalues (FEAST ENERGY), got {n_evals}")
        else:
            print("    PASS")

    # ── Landau: DENSE + INDEX ──────────────────────────────────────
    print("\n  Landau InAs B=5T (expected: DENSE + INDEX)...")
    rows = run_config(build_dir, LANDAU_INAS_K0, "landau")
    if rows is None:
        failures.append("landau: execution or parsing failed")
    else:
        k0, evals = rows[0]
        n_evals = len(evals)
        print(f"    k=0 eigenvalues: {n_evals}")
        # INDEX mode returns num_cb + num_vb = 4 + 4 = 8
        expected = 4 + 4
        if n_evals < expected:
            failures.append(f"landau: expected >= {expected} eigenvalues (INDEX mode), got {n_evals}")
        else:
            print("    PASS")

    # ── Summary ────────────────────────────────────────────────────
    print()
    if failures:
        print(f"FAIL: {len(failures)} issue(s):")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    else:
        print("PASS: all 4 confinement modes produce correct results without [solver]")
        sys.exit(0)


if __name__ == "__main__":
    main()
