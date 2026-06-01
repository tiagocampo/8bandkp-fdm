#!/usr/bin/env python3
"""Rung 7 -- g-factor Validation (R7.1, R7.2).

Validates g-factor calculation against the Roth analytical formula
(Winkler 2003, Eq. 6.42) for bulk materials and a QW system.

R7.1: Bulk CB g-factor vs Roth formula (4 materials):
  - GaAs, InAsW, InSb, GaSb
  - Tolerance: 2% (Roth is 2-band approximation vs full 8-band)

R7.2: QW g-factor physical consistency checks (1 material):
  - AlSbW/GaSbW/InAsW QW
  - Checks: sign consistency with bulk Roth, |g_z| < |g_x| = |g_y|,
    in-plane g close to Roth bulk (5%), confinement reduces |g_z|

Usage:
    verify_8band_rung7_gfactor.py <build_dir> <source_dir>

    build_dir  -- path to build/ directory (contains src/gfactorCalculation)
    source_dir -- path to repo root (contains tests/regression/configs/)
"""
# COVERAGE: observable=g*_cb geometry=bulk material=GaAs ref=Roth1959
# COVERAGE: observable=g*_cb geometry=bulk material=InAsW ref=Roth1959
# COVERAGE: observable=g*_cb geometry=bulk material=InSb ref=Roth1959
# COVERAGE: observable=g*_cb geometry=bulk material=GaSb ref=Roth1959
# COVERAGE: observable=g*_cb geometry=QW material=InAsW/GaSbW/AlSbW ref=Roth1959
import os
import sys
import tempfile
import shutil

sys.path.insert(
    0,
    os.path.join(os.path.dirname(__file__), "..", "..", "tests", "integration"),
)
from star_helpers import (
    run_exe,
    roth_gfactor,
    parse_gfactor,
    compare_value,
    print_benchmark_header,
    format_benchmark_row,
)

# ---------------------------------------------------------------------------
# Material database -- EP/Eg/DeltaSO from parameters.f90
# Roth g-factor: g = 2 - 2*EP*DeltaSO / (3*Eg*(Eg + DeltaSO))
# ---------------------------------------------------------------------------
BULK_MATERIALS = {
    "GaAs": {
        "config": "gfactor_bulk_gaas_cb.toml",
        "EP": 28.8,
        "Eg": 1.519,
        "DeltaSO": 0.341,
        "tolerance": 0.02,  # 2%
    },
    "InAsW": {
        "config": "gfactor_bulk_inasw_cb.toml",
        "EP": 22.2,
        "Eg": 0.418,
        "DeltaSO": 0.38,
        "tolerance": 0.02,  # 2%
    },
    "InSb": {
        "config": "gfactor_bulk_insb_cb.toml",
        "EP": 23.3,
        "Eg": 0.235,
        "DeltaSO": 0.81,
        "tolerance": 0.02,  # 2%
    },
    "GaSb": {
        "config": "gfactor_bulk_gasb_cb.toml",
        "EP": 27.0,
        "Eg": 0.812,
        "DeltaSO": 0.76,
        "tolerance": 0.02,  # 2%
    },
}

# R7.2: QW g-factor physical consistency.
# The QW has AlSbW/GaSbW/InAsW layers; the CB g-factor is dominated by the
# InAsW well material.  Confinement modifies g-factors:
#   - In-plane (gx = gy): probes 2D dispersion, close to bulk Roth for
#     narrow-gap QWs where non-parabolicity dominates.
#   - Out-of-plane (gz): probes quantized confinement direction, typically
#     reduced in magnitude relative to in-plane due to subband quantization.
QW_MATERIAL = {
    "config": "gfactor_qw_cb.toml",
    "well_EP": 22.2,        # InAsW
    "well_Eg": 0.418,       # InAsW bulk gap
    "well_DeltaSO": 0.38,   # InAsW
    "inplane_tolerance": 0.05,  # 5% for in-plane g vs Roth bulk
}


def run_gfactor_test(build_dir, config_path, work_dir):
    """Run gfactorCalculation and return parsed g-factor tuple.

    Returns:
        (gx, gy, gz) on success, None on failure.
    """
    rc, output_dir = run_exe(build_dir, "gfactorCalculation", config_path, work_dir)
    if rc != 0:
        print(f"  FAIL: gfactorCalculation returned exit code {rc}")
        return None

    gf_path = os.path.join(output_dir, "gfactor.dat")
    if not os.path.isfile(gf_path):
        print(f"  FAIL: gfactor.dat not found in {output_dir}")
        return None

    try:
        gx, gy, gz = parse_gfactor(gf_path)
    except RuntimeError as e:
        print(f"  FAIL: could not parse gfactor.dat: {e}")
        return None

    return gx, gy, gz


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        print("  build_dir  -- path to build/ (contains src/gfactorCalculation)")
        print("  source_dir -- path to repo root")
        sys.exit(1)

    build_dir = os.path.abspath(sys.argv[1])
    source_dir = os.path.abspath(sys.argv[2])
    config_dir = os.path.join(source_dir, "tests", "regression", "configs")

    exe_path = os.path.join(build_dir, "src", "gfactorCalculation")
    if not os.path.isfile(exe_path):
        print(f"FAIL: gfactorCalculation not found at {exe_path}")
        sys.exit(1)

    failures = []
    all_pass = True

    # =================================================================
    # R7.1: Bulk CB g-factor vs Roth formula
    # =================================================================
    print("=" * 72)
    print("Rung 7 -- g-factor Validation")
    print("=" * 72)
    print()
    print("-" * 72)
    print("R7.1: Bulk CB g-factor vs Roth formula")
    print("-" * 72)

    print_benchmark_header()

    for mat_name, mat_data in BULK_MATERIALS.items():
        ep = mat_data["EP"]
        eg = mat_data["Eg"]
        dso = mat_data["DeltaSO"]
        tol = mat_data["tolerance"]
        g_roth = roth_gfactor(ep, eg, dso)

        config_path = os.path.join(config_dir, mat_data["config"])
        if not os.path.isfile(config_path):
            msg = f"{mat_name}: config not found: {config_path}"
            print(f"  FAIL: {msg}")
            failures.append(msg)
            all_pass = False
            continue

        tmpdir = tempfile.mkdtemp(prefix=f"rung7_{mat_name}_")
        try:
            result = run_gfactor_test(build_dir, config_path, tmpdir)
            if result is None:
                failures.append(f"{mat_name}: gfactorCalculation failed")
                all_pass = False
                continue

            gx, gy, gz = result

            # For bulk, all components should be equal (isotropic).
            # Use gz for comparison.
            passed, delta, row = compare_value(
                gz, g_roth, tol,
                f"g* ({mat_name})", ""
            )

            bench_row = format_benchmark_row(
                mat_name, "g*_cb (Roth)", gz, g_roth,
                "Winkler 2003 Eq. 6.42",
                f"Analytical ({tol*100:.0f}%)", delta,
                "PASS" if passed else "FAIL",
            )
            print(bench_row)

            if passed:
                print(f"  PASS: {mat_name} g*={gz:.4f}, Roth={g_roth:.4f}, "
                      f"delta={delta:.2e}")
            else:
                msg = (f"{mat_name}: g*={gz:.4f}, Roth={g_roth:.4f}, "
                       f"delta={delta:.4f} ({delta*100:.1f}%) >= {tol*100:.0f}%")
                print(f"  FAIL: {msg}")
                failures.append(msg)
                all_pass = False

        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    # =================================================================
    # R7.2: QW g-factor physical consistency checks
    # =================================================================
    print()
    print("-" * 72)
    print("R7.2: QW CB g-factor (AlSbW/GaSbW/InAsW)")
    print("-" * 72)

    qw_data = QW_MATERIAL
    g_roth_bulk = roth_gfactor(qw_data["well_EP"], qw_data["well_Eg"],
                               qw_data["well_DeltaSO"])
    tol_ip = qw_data["inplane_tolerance"]

    print(f"  QW config: {qw_data['config']}")
    print(f"  Well material: InAsW (EP={qw_data['well_EP']}, "
          f"Eg={qw_data['well_Eg']}, DeltaSO={qw_data['well_DeltaSO']})")
    print(f"  Roth bulk g* (InAsW): {g_roth_bulk:.4f}")

    config_path = os.path.join(config_dir, qw_data["config"])
    if not os.path.isfile(config_path):
        msg = f"QW: config not found: {config_path}"
        print(f"  FAIL: {msg}")
        failures.append(msg)
        all_pass = False
    else:
        tmpdir = tempfile.mkdtemp(prefix="rung7_qw_")
        try:
            result = run_gfactor_test(build_dir, config_path, tmpdir)
            if result is None:
                failures.append("QW: gfactorCalculation failed")
                all_pass = False
            else:
                gx, gy, gz = result
                print(f"  Computed g*: gx={gx:.4f}, gy={gy:.4f}, gz={gz:.4f}")

                # Check 1: In-plane isotropy (gx == gy for [001]-grown QW)
                iso_diff = abs(gx - gy)
                iso_tol = 1e-6
                if iso_diff < iso_tol:
                    print(f"  [R7.2a] In-plane isotropy: PASS "
                          f"(|gx-gy| = {iso_diff:.2e})")
                else:
                    msg = (f"QW in-plane isotropy broken: "
                           f"gx={gx:.4f}, gy={gy:.4f}, diff={iso_diff:.2e}")
                    print(f"  [R7.2a] FAIL: {msg}")
                    failures.append(msg)
                    all_pass = False

                # Check 2: Sign consistency with Roth bulk (should be negative
                # for narrow-gap CB)
                if gz < 0 and gx < 0:
                    print(f"  [R7.2b] Sign consistency: PASS "
                          f"(all components negative, matching Roth bulk)")
                else:
                    msg = (f"QW sign mismatch: gx={gx:.4f}, gy={gy:.4f}, "
                           f"gz={gz:.4f} (expected all negative)")
                    print(f"  [R7.2b] FAIL: {msg}")
                    failures.append(msg)
                    all_pass = False

                # Check 3: In-plane g vs Roth bulk (within tolerance).
                # For narrow-gap QWs, gx probes the in-plane k.p dispersion
                # which is close to the bulk value.
                passed_ip, delta_ip, _ = compare_value(
                    gx, g_roth_bulk, tol_ip,
                    "g*_ip (QW)", ""
                )
                if passed_ip:
                    print(f"  [R7.2c] In-plane g vs Roth bulk: PASS "
                          f"(gx={gx:.4f}, Roth={g_roth_bulk:.4f}, "
                          f"delta={delta_ip:.2e})")
                else:
                    msg = (f"QW in-plane g: gx={gx:.4f}, "
                           f"Roth={g_roth_bulk:.4f}, "
                           f"delta={delta_ip:.4f} ({delta_ip*100:.1f}%) "
                           f">= {tol_ip*100:.0f}%")
                    print(f"  [R7.2c] FAIL: {msg}")
                    failures.append(msg)
                    all_pass = False

                # Check 4: Confinement reduction of |gz| relative to |gx|.
                # For a [001] QW, quantization along z reduces the
                # out-of-plane g-factor magnitude relative to in-plane.
                if abs(gz) < abs(gx):
                    reduction = (abs(gx) - abs(gz)) / abs(gx) * 100
                    print(f"  [R7.2d] Confinement reduction |gz|<|gx|: PASS "
                          f"(|gz|/|gx| = {abs(gz)/abs(gx):.4f}, "
                          f"reduction = {reduction:.1f}%)")
                else:
                    msg = (f"QW confinement: |gz|={abs(gz):.4f} >= "
                           f"|gx|={abs(gx):.4f} (expected |gz| < |gx|)")
                    print(f"  [R7.2d] FAIL: {msg}")
                    failures.append(msg)
                    all_pass = False

        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    # =================================================================
    # Summary
    # =================================================================
    print()
    print("=" * 72)
    print("Summary")
    print("=" * 72)

    if failures:
        print(f"FAIL: {len(failures)} check(s) failed:")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    else:
        n_bulk = len(BULK_MATERIALS)
        print(f"PASS: all g-factor validation checks passed")
        print(f"  R7.1: {n_bulk} bulk materials tested vs Roth formula")
        print(f"  R7.2: 1 QW system tested (isotropy, sign, in-plane vs Roth, "
              f"confinement)")
        sys.exit(0)


if __name__ == "__main__":
    main()
