#!/usr/bin/env python3
"""Verify the dispersion-aware energy-window authority (issue #03).

The window authority (apply_solver_window in eigensolver.f90) is the single
source of [emin, emax] for every FEAST consumer reached via the standard
init / k-sweep paths. When the user leaves the window at the parser's auto
sentinel (emin = emax = 0), the authority builds the QW Hamiltonian CSR at
BOTH sweep endpoints (k = 0 and k = k_max), computes the Gershgorin bound
at each, and UNIONS them into ONE stable window for the whole sweep.

This test exercises that envelope path end-to-end through the bandStructure
executable: a GaAs/AlGaAs QW with FEAST + ENERGY + auto window over a wide
in-plane wave-vector sweep (k_max = 0.15 1/A). We assert FEAST returns the
REQUESTED number of eigenvalues (num_cb + num_vb) at the LAST k-point
(k = k_max), i.e. the dispersion extremum where a single-k = 0 window
would have under-delivered before issue #03.

Contrast with the pre-#03 behavior (documented, not asserted here): the old
main.f90 k-sweep block built only the k = k_max CSR and took that single
endpoint's Gershgorin bound — which under-delivered eigenvalues at the
opposite endpoint (k = 0) whenever the dispersion was asymmetric, and
vice-versa for the init path which used k = 0 only. The envelope unions
both, so the window is guaranteed to bracket the spectrum at every k-point
of the sweep. This regression test would FAIL on the pre-#03 code path
whenever the k = 0 bound alone did not cover the k_max eigenvalues.

Usage: verify_qw_window_envelope.py <build_dir> <source_dir>

# COVERAGE: observable=energy_window_envelope geometry=qw material=GaAs/AlGaAs
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from star_helpers import run_exe, parse_eigenvalues


# Wide-sweep QW config: FEAST + ENERGY + AUTO window (emin = emax = 0).
#
# The wave-vector sweep spans k = 0 .. 0.1 1/A (3 k-points). At HEAD (before
# issue #03) the k-sweep envelope block was unreachable (its guard
# `ecfg_bs%emin >= ecfg_bs%emax` was already stale after issue #02 wired
# derive_eigensolver to propagate the window), so the sweep fell through to
# the eigensolver_config type default window [-1, +1] eV. That default does
# NOT bracket this QW's spectrum: the lowest subbands sit near -1.65 eV at
# k = 0 and disperse down to -2.2 eV at k = 0.1 — entirely outside [-1, +1].
# Under the type-default window FEAST returns 0 eigenvalues at every k-point
# (the under-delivery the envelope is meant to fix).
#
# With the #03 envelope live (guard corrected to the [0,0] AUTO sentinel),
# the sweep derives the union of the k = 0 and k = k_max Gershgorin bounds —
# ~[-3.8, +4.9] eV for this config — which brackets the whole sweep, so
# FEAST delivers the requested 8 bands at every k-point.
#
# The narrow 30 A GaAs well + fd_step = 21 keeps the total eigenvalue count
# in the envelope window within FEAST's auto-grown m0 (the valence manifold
# of a QW is dense; a wider well or larger grid pushes the in-window count
# past m0 and trips info=3, which is an m0-tuning matter, not a window
# matter). num_cb = 2, num_vb = 6 -> 8 retained bands; m0 = 64 lets FEAST
# grow the subspace as needed for the ~168 in-window eigenvalues.

QW_WIDE_SWEEP_AUTO_WINDOW = """\
confinement = "qw"
FDorder = 2
fd_step = 21

[wave_vector]
mode = "kx"
max = 0.1
nsteps = 3

[bands]
num_cb = 2
num_vb = 6

[[material]]
name = "Al30Ga70As"
z_min = -200
z_max = 200

[[material]]
name = "GaAs"
z_min = -30
z_max = 30

[solver]
method = "FEAST"
mode = "ENERGY"
emin = 0.0
emax = 0.0
m0 = 64
"""

# Expected eigenvalue count at every k-point: num_cb + num_vb = 2 + 6 = 8.
# The envelope window is GUARANTEED to bracket the whole sweep, so FEAST
# must deliver exactly this count at every k-point, including the last.
EXPECTED_NEV = 8


def verify_envelope_covers_full_sweep(build_dir, source_dir):
    """Run a wide-k QW FEAST sweep with AUTO window and assert FEAST
    delivers the requested eigenvalue count at the LAST k-point."""
    print("\n" + "=" * 64)
    print("QW window envelope: wide-k sweep, AUTO window, full delivery")
    print("=" * 64)
    with tempfile.TemporaryDirectory(prefix="qw_window_env_") as work:
        cfg_path = os.path.join(work, "cfg.toml")
        with open(cfg_path, "w") as f:
            f.write(QW_WIDE_SWEEP_AUTO_WINDOW)
        rc, outdir = run_exe(build_dir, "bandStructure", cfg_path, work,
                             timeout=300)
        if rc != 0:
            print(f"  FAIL: bandStructure exited nonzero (rc={rc})")
            return False

        eig_path = os.path.join(outdir, "eigenvalues.dat")
        if not os.path.exists(eig_path):
            print("  FAIL: eigenvalues.dat not produced")
            return False
        rows = parse_eigenvalues(eig_path)
        if not rows:
            print("  FAIL: could not parse eigenvalues")
            return False

        n_k = len(rows)
        print(f"  Parsed {n_k} k-points from eigenvalues.dat")
        if n_k < 2:
            print("  FAIL: sweep must have >= 2 k-points to test the envelope")
            return False

        # The acceptance criterion: FEAST delivers EXPECTED_NEV NON-ZERO
        # eigenvalues at the LAST k-point (k = k_max, the dispersion
        # extremum). The envelope's whole point is that the window derived
        # from the two endpoints brackets the spectrum at the extremum.
        #
        # Crucially, we count NON-ZERO eigenvalues: FEAST under-delivery
        # zero-fills the missing bands (the "missing bands zero-filled"
        # warning path in main.f90), so a raw column count of 8 would hide
        # an under-delivery of 6 real + 2 zero-filled. Real semiconductor
        # subband eigenvalues are never exactly 0.0 in this QW (the
        # spectrum sits near -1 eV), so 0.0 is an unambiguous fill marker.
        last_k, last_evals = rows[-1]
        n_last_real = sum(1 for v in last_evals if v != 0.0)
        n_last = len(last_evals)
        print(f"  Last k-point (|k|={last_k:.4f} 1/A): {n_last_real} real / "
              f"{n_last} total eigenvalues (expected {EXPECTED_NEV} real)")
        if n_last_real != EXPECTED_NEV:
            print(f"  FAIL: under-delivery at k_max — envelope did not cover "
                  f"the sweep extremum ({n_last_real} real != {EXPECTED_NEV})")
            return False

        # Stronger: the real-eigenvalue count must be STABLE across ALL
        # k-points. A per-k moving window or an under-covering window would
        # show a real-count drop somewhere mid-sweep. This is the ADR-0005
        # "one stable window per sweep" invariant made observable.
        real_counts = {sum(1 for v in r[1] if v != 0.0) for r in rows}
        if len(real_counts) != 1:
            print(f"  FAIL: real-eigenvalue count varies across sweep: "
                  f"{sorted(real_counts)}")
            return False
        stable_count = real_counts.pop()
        if stable_count != EXPECTED_NEV:
            print(f"  FAIL: stable real count {stable_count} != "
                  f"expected {EXPECTED_NEV}")
            return False

        # Sanity: no NaN / Inf (FEAST divide-by-zero on a band-edge window
        # would manifest here).
        import math
        for ki, (_k, evals) in enumerate(rows):
            for v in evals:
                if math.isnan(v) or math.isinf(v):
                    print(f"  FAIL: NaN/Inf eigenvalue at k-point {ki}")
                    return False

        print(f"  PASS: envelope delivered {stable_count} eigenvalues stably "
              f"across all {n_k} k-points of the wide sweep")
        return True


def main():
    if len(sys.argv) != 3:
        print("Usage: verify_qw_window_envelope.py <build_dir> <source_dir>")
        sys.exit(2)
    build_dir = os.path.abspath(sys.argv[1])
    source_dir = os.path.abspath(sys.argv[2])

    if not verify_envelope_covers_full_sweep(build_dir, source_dir):
        sys.exit(1)

    print("\nPASS: QW window envelope delivers the full eigenvalue set "
          "across a wide k-sweep (issue #03)")
    sys.exit(0)


if __name__ == "__main__":
    main()
