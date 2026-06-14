#!/usr/bin/env python3
"""Verify QW optics: FEAST+ENERGY (sweep-envelope window) vs DENSE+INDEX.

This is the FEAST-enable regression gate for issue #09 (ADR 0005).

Before issue #09 the QW optics block in main_optics.f90 HARDCODED the
dense backend:

    qw_cfg%method = 'DENSE'        # ignores cfg%solver%method
    qw_cfg%mode   = EIGEN_MODE_INDEX
    qw_cfg%il/iu/nev = ...         # hand-assembled

so requesting method=FEAST on QW optics was silently overridden to
DENSE - a silent-correction violation of the project's validation rule
(one executable, two behaviors: wire optics honored [solver].method,
QW optics did not).

Issue #09 routes QW optics through the single derivation seam
(derive_eigensolver, issue #02) so it honors [solver].method, and makes
FEAST work for the optics k_par sweep via the window authority's
dispersion-aware envelope variant (asw_envelope, issue #03): the
Gershgorin bounds at the sweep's two endpoints (k_par=0 and k_par=max)
are unioned into ONE stable window per sweep - the same shape as the
band-structure k-sweep in main.f90.

This script runs the opticalProperties executable on the SAME QW config
twice - once method=DENSE, once method=FEAST - and asserts:

  1. PHYSICS GATE: the output spectra (absorption_TE, absorption_TM, and
     whichever of gain/spontaneous/ISBT the config enables) match
     between the two backends within 1e-6 relative. Divergent optics
     would mean the FEAST window truncated the spectrum or the band
     indexing shifted between INDEX and ENERGY mode.
  2. DISPATCH GATE: the FEAST run's stdout contains the optics-loop
     backend diagnostic ("QW optics eigensolver: FEAST backend"). At
     HEAD (pre-#09) the optics loop hardcoded DENSE and printed no such
     line, so the FEAST-tagged run silently ran DENSE - this gate
     catches that override. (Bit-identity of spectra is NOT a reliable
     signal: on small problems FEAST and LAPACK can agree to full double
     precision, making the spectra bit-identical even when genuinely
     different backends ran.)

The config is an Al30Ga70As/GaAs QW (same material system as the golden
qw_gaas_algaas_optics config) with a compact grid (fd_step=41 -> N=328)
and a real k_par sweep (nsteps=21) so the envelope over the two sweep
endpoints is exercised. The FEAST run uses the AUTO window
(emin = emax = 0) so the envelope authority is exercised end-to-end.

Usage: verify_qw_optics_dense_feast.py <build_dir> <source_dir>

# COVERAGE: observable=absorption spectra geometry=qw material=GaAs/AlGaAs backend=FEAST-vs-DENSE
"""

import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))

TOL = 1.0e-6  # spectra relative tolerance for FEAST-vs-DENSE agreement

FD_STEP = 41
N_TOTAL = 8 * FD_STEP  # full matrix dimension

# Spectra files to compare (only those the config enables will be checked).
SPECTRA_FILES = [
    "absorption_TE.dat",
    "absorption_TM.dat",
    "absorption_ISBT.dat",
    "gain_TE.dat",
    "gain_TM.dat",
    "spontaneous_TE.dat",
    "spontaneous_TM.dat",
]


# Shared QW config body (Al30Ga70As/GaAs, the golden optics material
# system). The only field that varies between the two runs is
# [solver].method (+ mode/window). Real k_par sweep so the FEAST
# envelope over the two sweep endpoints is exercised. gain + ISBT
# enabled so multiple spectra channels are compared.
QW_BODY = """\
confinement = "qw"
FDorder = 2
fd_step = {fd_step}
which_band = 0
band_idx = 1

[wave_vector]
mode = "kx"
max = 0.1
nsteps = 21

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

[optics]
linewidth_lorentzian = 0.030
linewidth_gaussian = 0.005
refractive_index = 3.3
E_min = 1.3
E_max = 2.0
num_energy_points = 200
temperature = 300.0
carrier_density = 0.0
gain_enabled = true
gain_carrier_density = 3.0e12
ISBT = true
spontaneous = true
spin_resolved = false
"""

# DENSE backend: AUTO resolves to DENSE+INDEX for QW (resolve_solver_defaults,
# ADR 0004). This is the golden default, so this run is the baseline.
QW_DENSE = QW_BODY.format(fd_step=FD_STEP) + """\
[solver]
method = "DENSE"
"""

# FEAST backend: AUTO window (emin = emax = 0). The optics sweep routes
# this through apply_solver_window's dispersion-aware envelope (union of
# Gershgorin bounds at k_par=0 and k_par=max), giving ONE stable window
# per sweep. m0 large enough to hold the numcb+numvb=12 retained bands
# (the info=3 retry loop also grows M0 if needed).
QW_FEAST = QW_BODY.format(fd_step=FD_STEP) + """\
[solver]
method = "FEAST"
emin = 0.0
emax = 0.0
m0 = {m0}
""".format(m0=N_TOTAL)


def _run_optics(build_dir, work, config_text, label):
    """Write config, run opticalProperties, return (rc, stdout, outdir)."""
    cfg_path = os.path.join(work, f"{label}.toml")
    with open(cfg_path, "w") as f:
        f.write(config_text)
    dst_cfg = os.path.join(work, f"input_{label}.toml")
    shutil.copy2(cfg_path, dst_cfg)

    exe_path = os.path.abspath(os.path.join(build_dir, "src", "opticalProperties"))
    if not os.path.isfile(exe_path):
        print(f"  FAIL: opticalProperties not found at {exe_path}")
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
        print(f"  FAIL: opticalProperties timed out ({label})")
        return None

    return result.returncode, result.stdout, os.path.join(run_dir, "output")


def _load_spectrum(path):
    """Load a spectrum file as (energy, value) arrays, or None if absent."""
    if not os.path.isfile(path):
        return None
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data[:, 0], data[:, 1]


def _compare_spectrum(name, dense, feast):
    """Compare one spectrum between backends; return (ok, worst_rel_err)."""
    ed, vd = dense
    ef, vf = feast
    if len(ed) != len(ef):
        print(f"  FAIL: {name} energy-grid length mismatch "
              f"DENSE={len(ed)} FEAST={len(ef)}")
        return False, float("inf")
    # Energy grids must agree (same E_min/E_max/num_energy_points).
    e_err = float(np.max(np.abs(ed - ef)))
    if e_err > 1.0e-10:
        print(f"  FAIL: {name} energy grid shifted (max diff {e_err:.2e})")
        return False, float("inf")
    # Relative tolerance with absolute floor for near-zero values.
    denom = np.maximum(np.abs(vd), 1.0e-12)
    rel = np.abs(vd - vf) / denom
    worst = float(np.max(rel))
    # Only flag points where the dense value is non-negligible; near-zero
    # regions are dominated by broadening tail noise.
    mask = np.abs(vd) > 1.0e-8 * max(float(np.max(np.abs(vd))), 1.0)
    if np.any(mask):
        worst_sig = float(np.max(rel[mask]))
    else:
        worst_sig = worst
    if worst_sig > TOL:
        print(f"  FAIL: {name} worst rel diff {worst_sig:.2e} (tol {TOL:.0e}) "
              f"[DENSE max={float(np.max(np.abs(vd))):.4e}, "
              f"FEAST max={float(np.max(np.abs(vf))):.4e}]")
        return False, worst_sig
    print(f"  PASS: {name} worst rel diff {worst_sig:.2e} "
          f"(max |alpha|={float(np.max(np.abs(vd))):.4e})")
    return True, worst_sig


def verify_qw_optics_dense_vs_feast(build_dir, source_dir):
    print("\n" + "=" * 64)
    print("QW optics: FEAST+ENERGY (sweep-envelope window) vs DENSE+INDEX")
    print(f"  config: Al30Ga70As/GaAs QW, fd_step={FD_STEP}, N={N_TOTAL}, "
          f"k_par sweep nsteps=21")
    print("=" * 64)

    with tempfile.TemporaryDirectory(prefix="qw_opt_ds_") as work:
        print("\n  [1/2] DENSE baseline run...")
        dense = _run_optics(build_dir, work, QW_DENSE, "dense")
        if dense is None:
            return False
        rc_dense, stdout_dense, outdir_dense = dense
        if rc_dense != 0:
            print(f"  FAIL: opticalProperties exited {rc_dense} (dense)")
            return False

        print("\n  [2/2] FEAST+ENERGY (AUTO envelope window) run...")
        feast = _run_optics(build_dir, work, QW_FEAST, "feast")
        if feast is None:
            return False
        rc_feast, stdout_feast, outdir_feast = feast
        if rc_feast != 0:
            print(f"  FAIL: opticalProperties exited {rc_feast} (feast)")
            return False

        ok = True

        # --- PHYSICS GATE: spectra agreement ------------------------------
        # Compare every spectrum channel the config enables. The two
        # backends must agree within tolerance: divergent optics would mean
        # the FEAST window truncated the spectrum or the band indexing
        # shifted between INDEX and ENERGY mode.
        print(f"\n  Spectra agreement (tol {TOL:.0e} relative):")
        compared = 0
        worst_overall = 0.0
        for fname in SPECTRA_FILES:
            d = _load_spectrum(os.path.join(outdir_dense, fname))
            f = _load_spectrum(os.path.join(outdir_feast, fname))
            if d is None and f is None:
                # Config did not enable this channel on either run: skip.
                continue
            if d is None or f is None:
                print(f"  FAIL: {fname} produced under one backend but not "
                      f"the other (DENSE={'yes' if d else 'no'}, "
                      f"FEAST={'yes' if f else 'no'})")
                ok = False
                continue
            spectrum_ok, worst = _compare_spectrum(fname, d, f)
            if not spectrum_ok:
                ok = False
            compared += 1
            if worst > worst_overall:
                worst_overall = worst

        if compared == 0:
            print("  FAIL: no spectra files found to compare")
            ok = False
        else:
            print(f"\n  Compared {compared} spectrum file(s); "
                  f"worst rel diff across all = {worst_overall:.2e}")

        # --- DISPATCH GATE ------------------------------------------------
        # The FEAST run must actually have taken the FEAST path in the
        # OPTICS LOOP (not just in setup_init, which always derives the
        # solver config). The post-#09 optics loop prints a backend
        # diagnostic line naming the resolved method. At HEAD (pre-#09) the
        # optics loop hardcoded DENSE and printed no such line, so the
        # FEAST-tagged run silently ran DENSE. This gate catches that
        # override by requiring the optics-loop-specific FEAST diagnostic
        # in the FEAST run's stdout. (We do NOT use bit-identity of the
        # spectra as the signal: on small problems FEAST and LAPACK can
        # agree to full double precision, making the spectra bit-identical
        # even when genuinely different backends ran.)
        print("\n  Dispatch gate (FEAST path actually taken in optics loop):")
        optics_feast_diag = "QW optics eigensolver: FEAST backend"
        if optics_feast_diag in stdout_feast:
            print(f"  PASS: optics-loop FEAST diagnostic present in stdout "
                  f"(FEAST path honored)")
        else:
            print(f"  FAIL: optics-loop FEAST diagnostic absent from "
                  f"stdout -- method=FEAST was silently overridden to "
                  f"DENSE in the optics loop (issue #09 not fixed)")
            ok = False

        return ok


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        sys.exit(2)
    build_dir = os.path.abspath(sys.argv[1])
    source_dir = os.path.abspath(sys.argv[2])

    exe = os.path.join(build_dir, "src", "opticalProperties")
    if not os.path.isfile(exe):
        print(f"SKIP: opticalProperties not found at {exe}")
        sys.exit(0)

    if not verify_qw_optics_dense_vs_feast(build_dir, source_dir):
        print("\nFAIL: QW optics FEAST-vs-DENSE mismatch or FEAST path not "
              "taken (issue #09)")
        sys.exit(1)

    print("\nPASS: QW optics under FEAST+ENERGY matches DENSE+INDEX and the "
          "FEAST path is honored (issue #09)")
    sys.exit(0)


if __name__ == "__main__":
    main()
