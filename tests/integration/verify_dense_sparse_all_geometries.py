#!/usr/bin/env python3
"""Verify eigensolver dispatch produces bit-identical eigenvalues across
all four geometries (bulk, QW, wire, Landau).

This extends the dense-vs-sparse verification-ladder rung. The original
rung covered:
  - wire: test_wire_dense_sparse_consistency.sh (DENSE+ENERGY vs FEAST)
  - QW:   verify_qw_sparse_solver.py (DENSE+INDEX vs FEAST+ENERGY)

This script extends coverage to the remaining geometries:

  bulk:   DENSE+FULL  vs  DENSE+INDEX   (8x8 — FEAST is invalid for N=8;
          asserts mode-dispatch consistency: FULL and INDEX return the
          same eigenvalues for the same matrix)
  landau: AUTO (resolves DENSE+INDEX) vs explicit DENSE+INDEX
          (asserts the AUTO resolver dispatches to the correct backend
          and mode for Landau, producing identical eigenvalues)

For QW and wire, the genuine dense-vs-sparse (CSR) cross-check is
already covered by the existing tests cited above; this script focuses
on the geometries not yet covered. Once issue #10 (FEAST-enable Landau)
lands, the Landau case will gain a DENSE-vs-FEAST assertion.

The QW/wire cases are included here as well, re-running the DENSE-vs-FEAST
comparison so this script is a single-entry-point regression for the
full dispatch surface across all four geometries.

Bit-identical assertion (within published tolerance, 1e-8).

Usage: verify_dense_sparse_all_geometries.py <build_dir> <source_dir>

# COVERAGE: observable=eigensolver_dispatch_consistency geometry=bulk material=GaAs
# COVERAGE: observable=dense_sparse_agreement geometry=qw material=GaAs/AlGaAs
# COVERAGE: observable=dense_sparse_agreement geometry=wire material=GaAs
# COVERAGE: observable=eigensolver_dispatch_consistency geometry=landau material=GaAs
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))

from star_helpers import run_exe, parse_eigenvalues

TOL = 1.0e-8  # golden-regression tolerance (bit-identical within published tol)


# -- Bulk (8x8): DENSE+FULL vs DENSE+INDEX -----------------------------
# FEAST is invalid for the 8x8 bulk matrix (contour integration on N=8
# is degenerate). The regression here is mode-dispatch consistency:
# FULL (zheev, all eigenvalues) and INDEX (zheevx, first 8) must agree.
BULK_FULL = """\
confinement = "bulk"
FDorder = 2
fd_step = 1

[wave_vector]
mode = "k0"
max = 0.0
nsteps = 1

[bands]
num_cb = 4
num_vb = 4

[[material]]
name = "GaAs"

[solver]
method = "DENSE"
mode = "FULL"
"""

BULK_INDEX = """\
confinement = "bulk"
FDorder = 2
fd_step = 1

[wave_vector]
mode = "k0"
max = 0.0
nsteps = 1

[bands]
num_cb = 4
num_vb = 4

[[material]]
name = "GaAs"

[solver]
method = "DENSE"
mode = "INDEX"
"""


# -- QW (8N x 8N): DENSE+INDEX (default) vs FEAST+AUTO-envelope (k=0) ----
# Both backends extract the gap-straddling [il, iu] bands via the QW
# band-structure path (Task 7 / reconcile_band_slice). DENSE defaults to
# DENSE+INDEX, which pre-slices via zheevx range 'I' so its lowest-M copy
# IS exactly [il, iu]. FEAST uses the AUTO envelope (emin = emax = 0) so
# apply_solver_window's dispersion-aware authority picks the window; m0 =
# 8*fd_step = 328 is large enough for FEAST to converge on the wide
# envelope. num_cb=4 + num_vb=8 -> 12 gap-straddling bands retained.
# This mirrors verify_qw_bandstructure_dense_feast.py. Note: a narrow user
# ENERGY window (e.g. [-1, 2]) is incompatible with the QW path's
# global-index extraction (reconcile_band_slice needs nev_found >= iu to
# offset into the full spectrum; a narrow window returns fewer), so the
# two backends cannot be compared in ENERGY mode on a narrow window.
QW_DENSE = """\
confinement = "qw"
FDorder = 2
fd_step = 41

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
"""

QW_FEAST = """\
confinement = "qw"
FDorder = 2
fd_step = 41

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
method = "FEAST"
emin = 0.0
emax = 0.0
m0 = 328
"""


# -- Wire (CSR): DENSE+ENERGY vs FEAST+ENERGY --------------------------
# Reuses the proven wire_gaas_dense_sparse parameter set (24 eigenvalues,
# window [-1.5, 2.0]). Both backends operate on the same CSR Hamiltonian.
# nsteps=2 (the nsteps=1 wire path has a pre-existing LAPACK info=1 bug
# unrelated to this issue; nsteps>=2 is the working path).
WIRE_DENSE = """\
confinement = "wire"
FDorder = 2
fd_step = 1
which_band = 0
band_idx = 1

[wave_vector]
mode = "kz"
max = 0.1
nsteps = 2

[bands]
num_cb = 8
num_vb = 16

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


[solver]
method = "DENSE"
mode = "ENERGY"
emin = -1.5
emax = 2.0
"""

WIRE_FEAST = """\
confinement = "wire"
FDorder = 2
fd_step = 1
which_band = 0
band_idx = 1

[wave_vector]
mode = "kz"
max = 0.1
nsteps = 2

[bands]
num_cb = 8
num_vb = 16

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


[solver]
method = "FEAST"
mode = "ENERGY"
emin = -1.5
emax = 2.0
"""


# -- Landau (8nx x 8nx): AUTO vs explicit DENSE+INDEX ------------------
# The Landau setup path hardcodes solve_dense (issue #10 will FEAST-enable
# it). The regression here is AUTO-mode resolution: the resolver must pick
# DENSE+INDEX for Landau, matching the explicit DENSE+INDEX result.
LANDAU_AUTO = """\
confinement = "landau"
FDorder = 2
fd_step = 1

[wave_vector]
mode = "k0"
max = 0.0
nsteps = 1

[bands]
num_cb = 4
num_vb = 4

[[material]]
name = "GaAs"

[landau]
nx = 60
width = 2000.0

[b_field]
components = [0.0, 0.0, 5.0]
"""

LANDAU_DENSE_INDEX = """\
confinement = "landau"
FDorder = 2
fd_step = 1

[wave_vector]
mode = "k0"
max = 0.0
nsteps = 1

[bands]
num_cb = 4
num_vb = 4

[[material]]
name = "GaAs"

[landau]
nx = 60
width = 2000.0

[b_field]
components = [0.0, 0.0, 5.0]

[solver]
method = "DENSE"
mode = "INDEX"
"""


def _run_and_parse(build_dir, work, config_text, label):
    """Write config, run bandStructure, return parsed eigenvalues."""
    cfg_path = os.path.join(work, f"{label}.toml")
    with open(cfg_path, "w") as f:
        f.write(config_text)
    rc, outdir = run_exe(build_dir, "bandStructure", cfg_path, work)
    if rc != 0:
        print(f"  FAIL: bandStructure exited {rc} ({label})")
        return None
    eig_path = os.path.join(outdir, "eigenvalues.dat")
    if not os.path.isfile(eig_path):
        print(f"  FAIL: eigenvalues.dat not produced ({label})")
        return None
    rows = parse_eigenvalues(eig_path)
    if not rows:
        print(f"  FAIL: could not parse eigenvalues ({label})")
        return None
    return rows


def _compare(label, rows_a, rows_b, sublabel_a, sublabel_b):
    """Compare two eigenvalue-row sets; return True if all match."""
    if len(rows_a) != len(rows_b):
        print(f"  FAIL: k-point count mismatch {len(rows_a)} ({sublabel_a}) vs "
              f"{len(rows_b)} ({sublabel_b}) ({label})")
        return False
    worst = 0.0
    for ki, (ra, rb) in enumerate(zip(rows_a, rows_b)):
        ea, eb = ra[1], rb[1]
        m = min(len(ea), len(eb))
        if m == 0:
            print(f"  FAIL: k={ki} produced zero eigenvalues ({label})")
            return False
        for a, b in zip(ea[:m], eb[:m]):
            d = abs(a - b)
            if d > worst:
                worst = d
            if d > TOL:
                print(f"  FAIL: {label} k={ki} eigenvalue mismatch "
                      f"{a:.10e} ({sublabel_a}) vs {b:.10e} ({sublabel_b}) "
                      f"(diff={d:.2e})")
                return False
    print(f"  PASS: {label} {sublabel_a} vs {sublabel_b} match "
          f"(worst diff {worst:.2e}, {len(rows_a)} k-point(s))")
    return True


def verify_geometry(build_dir, label, cfg_a, cfg_b, sublabel_a, sublabel_b):
    """Run one geometry's two-way comparison."""
    print("\n" + "=" * 60)
    print(f"Eigensolver dispatch consistency: {label} ({sublabel_a} vs {sublabel_b})")
    print("=" * 60)
    with tempfile.TemporaryDirectory(prefix=f"ds_{label}_") as work:
        rows_a = _run_and_parse(build_dir, work, cfg_a, f"{label}_A")
        if rows_a is None:
            return False
        rows_b = _run_and_parse(build_dir, work, cfg_b, f"{label}_B")
        if rows_b is None:
            return False
        return _compare(label, rows_a, rows_b, sublabel_a, sublabel_b)


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        sys.exit(1)

    build_dir = os.path.abspath(sys.argv[1])
    exe = os.path.join(build_dir, "src", "bandStructure")
    if not os.path.isfile(exe):
        print(f"SKIP: bandStructure not found at {exe}")
        sys.exit(0)

    geometries = [
        ("bulk",   BULK_FULL,        BULK_INDEX,        "DENSE+FULL",  "DENSE+INDEX"),
        ("qw",     QW_DENSE,         QW_FEAST,          "DENSE+INDEX", "FEAST+AUTO"),
        ("wire",   WIRE_DENSE,       WIRE_FEAST,        "DENSE+ENERGY","FEAST+ENERGY"),
        ("landau", LANDAU_AUTO,      LANDAU_DENSE_INDEX,"AUTO",        "DENSE+INDEX"),
    ]

    all_passed = True
    for label, cfg_a, cfg_b, sl_a, sl_b in geometries:
        if not verify_geometry(build_dir, label, cfg_a, cfg_b, sl_a, sl_b):
            all_passed = False

    print()
    if all_passed:
        print("PASS: eigensolver dispatch produces bit-identical eigenvalues "
              "across all 4 geometries (bulk/QW/wire/Landau)")
        sys.exit(0)
    print("FAIL: eigenvalue mismatch detected across geometries")
    sys.exit(1)


if __name__ == "__main__":
    main()
