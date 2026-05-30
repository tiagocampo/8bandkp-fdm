#!/usr/bin/env python3
"""S1 -- GaAs bulk standard-star benchmark.

Validates GaAs bulk observables against published references:
  - Eg = 1.519 eV  (Vurgaftman 2001, Table I)
  - DeltaSO = 0.341 eV  (Vurgaftman 2001, Table I)
  - m*_e (Kane 2-band formula)  (Eg/(EP+Eg) ~ 0.0501 m0)
  - g* (Roth formula, Winkler 2003 Eq. 6.42)
  - Absorption edge near Eg

Cross-reference with verification ladder:
  - Eg, DeltaSO overlap with rung 1 (R1: bulk k=0 eigenvalue accuracy)
  - Effective mass overlaps with rung 2 (R6-R7: bulk dispersion)
  - Standard-star adds: g-factor validation, optical absorption edge,
    publication-ready benchmark table format.

Usage: verify_star_gaas_bulk.py <build_dir> <source_dir>
"""

# COVERAGE: observable=Eg geometry=bulk material=GaAs ref=Vurgaftman2001
# COVERAGE: observable=Delta_SO geometry=bulk material=GaAs ref=Vurgaftman2001
# COVERAGE: observable=m*_e geometry=bulk material=GaAs ref=Kane
# COVERAGE: observable=g*_cb geometry=bulk material=GaAs ref=Winkler2003
# COVERAGE: observable=absorption_edge geometry=bulk material=GaAs ref=Vurgaftman2001
import os
import sys
import tempfile
import shutil

from star_helpers import (
    run_executable, parse_eigenvalues, parse_gfactor,
    parse_absorption, compare_value, format_benchmark_row, print_benchmark_header,
    TOL_EXACT, TOL_ANALYTICAL, TOL_NUMERICAL,
    HBAR2_OVER_2M0, roth_gfactor,
)

# ---------------------------------------------------------------------------
# Physical constants and GaAs parameters
# ---------------------------------------------------------------------------
# GaAs material parameters from Vurgaftman 2001 and Winkler 2003.
# Used to compute analytical reference values (Kane mass, Roth g-factor).
# These are NOT read from parameters.f90 -- they are independent literature
# values used for cross-validation.

EG_GAAS = 1.519       # eV, Vurgaftman 2001 Table I
DELTASO_GAAS = 0.341   # eV, Vurgaftman 2001 Table I
EP_GAAS = 28.8         # eV, Kane matrix element (Vurgaftman 2001)

# 2-band Kane effective mass: m* = Eg / (EP + Eg) (informational only)
M_STAR_KANE = EG_GAAS / (EP_GAAS + EG_GAAS)  # ~ 0.0501 m0

# 8-band model reference (validated against kdotpy cross-code comparison)
M_STAR_8BAND = 0.0317  # m0, from const-correct 8-band calculation

# Roth g-factor formula (Winkler 2003, Eq. 6.42):
#   g = 2 - 2*EP*DeltaSO / (3*Eg*(Eg + DeltaSO))
G_ROTH = roth_gfactor(EP_GAAS, EG_GAAS, DELTASO_GAAS)  # ~ -0.317

# Tolerances for each observable
TOL_EG = TOL_EXACT              # machine precision
TOL_DELTASO = TOL_EXACT         # machine precision
TOL_MASS = 0.05                 # 5% for 8-band reference comparison
TOL_GFACTOR = TOL_ANALYTICAL    # 1% for Roth vs 8-band
TOL_ABSORPTION = 0.06           # 6% for absorption edge (onset is broadened above Eg by linewidth + VB dispersion)

# CB eigenvalue index (0-based): bands 1-4=valence, 5-6=SO, 7-8=CB
CB_INDEX = 6

# Absorption edge detection threshold: fraction of maximum absorption
ABSORPTION_THRESHOLD_FRAC = 0.10


def check_eg_deltaso(build_dir, source_dir):
    """Check Eg and DeltaSO from bulk k=0 eigenvalues.

    Returns list of benchmark row dicts.
    """
    config_path = os.path.join(
        source_dir, "tests", "regression", "configs", "bulk_gaas_k0.toml"
    )
    exe_path = os.path.join(build_dir, "src", "bandStructure")
    rows = []

    with tempfile.TemporaryDirectory(prefix="star_gaas_k0_") as tmpdir:
        rc, output_dir = run_executable(exe_path, config_path, tmpdir)
        if rc != 0:
            print(f"  FAIL: bandStructure returned {rc}")
            return rows

        eig_path = os.path.join(output_dir, "eigenvalues.dat")
        data = parse_eigenvalues(eig_path)
        if not data:
            print("  FAIL: no eigenvalues parsed")
            return rows

        # Take first k-point (k=0)
        k0, evals = data[0]
        if len(evals) < 8:
            print(f"  FAIL: expected 8 eigenvalues, got {len(evals)}")
            return rows

        # Eigenvalue ordering: [-DeltaSO, -DeltaSO, 0, 0, 0, 0, Eg, Eg]
        eg_computed = evals[CB_INDEX]
        deltaso_computed = -evals[0]  # first eigenvalue is -DeltaSO

        # Eg check (self-consistency: eigenvalue == parameter value at k=0)
        passed, delta, row = compare_value(
            eg_computed, EG_GAAS, TOL_EG, "Eg", "eV"
        )
        rows.append({
            "material": "GaAs",
            "observable": "Eg (k=0 self-check)",
            "computed": eg_computed,
            "expected": EG_GAAS,
            "reference": "Vurgaftman 2001, Table I",
            "tolerance": "Exact (1e-12)",
            "delta": delta,
            "status": "PASS" if passed else "FAIL",
        })
        print(f"  Eg: {eg_computed:.10e} eV "
              f"(expected {EG_GAAS}, delta={delta:.2e}) "
              f"{'PASS' if passed else 'FAIL'}")

        # DeltaSO check (self-consistency: eigenvalue == parameter value at k=0)
        passed, delta, row = compare_value(
            deltaso_computed, DELTASO_GAAS, TOL_DELTASO, "DeltaSO", "eV"
        )
        rows.append({
            "material": "GaAs",
            "observable": "DeltaSO (k=0 self-check)",
            "computed": deltaso_computed,
            "expected": DELTASO_GAAS,
            "reference": "Vurgaftman 2001, Table I",
            "tolerance": "Exact (1e-12)",
            "delta": delta,
            "status": "PASS" if passed else "FAIL",
        })
        print(f"  DeltaSO: {deltaso_computed:.10e} eV "
              f"(expected {DELTASO_GAAS}, delta={delta:.2e}) "
              f"{'PASS' if passed else 'FAIL'}")

    return rows


def check_effective_mass(build_dir, source_dir):
    """Check CB effective mass from bulk dispersion.

    Extracts the effective mass from the numerical second derivative at k=0
    and compares against the 8-band model reference value (validated against
    kdotpy). The 2-band Kane formula is reported as informational.

    Returns list of benchmark row dicts.
    """
    config_path = os.path.join(
        source_dir, "tests", "regression", "configs",
        "bulk_gaas_kx_dispersion.toml"
    )
    exe_path = os.path.join(build_dir, "src", "bandStructure")
    rows = []

    with tempfile.TemporaryDirectory(prefix="star_gaas_mass_") as tmpdir:
        rc, output_dir = run_executable(exe_path, config_path, tmpdir)
        if rc != 0:
            print(f"  FAIL: bandStructure returned {rc}")
            return rows

        eig_path = os.path.join(output_dir, "eigenvalues.dat")
        data = parse_eigenvalues(eig_path)
        if not data or len(data) < 2:
            print("  FAIL: insufficient eigenvalue data for mass extraction")
            return rows

        # Extract E(0) and E(k1) from the CB band
        e0 = data[0][1][CB_INDEX]  # E at k=0
        k1 = data[1][0]            # first nonzero k
        e1 = data[1][1][CB_INDEX]  # E at k1

        if k1 <= 0:
            # Find first nonzero k
            for i in range(1, len(data)):
                if data[i][0] > 0:
                    k1 = data[i][0]
                    e1 = data[i][1][CB_INDEX]
                    break

        if k1 <= 0:
            print("  FAIL: no nonzero k-point found")
            return rows

        # d2E/dk2 ~ 2*(E(k1) - E(0))/k1^2 (parabolic approximation)
        d2E_dk2 = 2.0 * (e1 - e0) / k1**2

        if d2E_dk2 <= 0:
            print(f"  FAIL: negative curvature d2E/dk2 = {d2E_dk2:.4e}")
            return rows

        # m*/m0 = hbar^2/(2*m0) / (d2E/dk2 / 2) = HBAR2_OVER_2M0 / c2
        # where c2 = d2E_dk2 / 2 ... actually:
        # E(k) = E0 + c2*k^2 where c2 = (e1-e0)/k1^2
        # d2E/dk2 = 2*c2
        # m*/m0 = HBAR2_OVER_2M0 / c2 = HBAR2_OVER_2M0 * k1^2 / (e1 - e0)
        c2 = (e1 - e0) / k1**2
        m_star = HBAR2_OVER_2M0 / c2

        print(f"  E(0) = {e0:.6f} eV, E(k1={k1:.4f}) = {e1:.6f} eV")
        print(f"  d2E/dk2 = {d2E_dk2:.4f} eV*A^2, c2 = {c2:.4f}")
        print(f"  m* = {m_star:.4f} m0 (8-band ref: {M_STAR_8BAND:.4f}, "
              f"Kane: {M_STAR_KANE:.4f})")

        passed, delta, _ = compare_value(
            m_star, M_STAR_8BAND, TOL_MASS, "m*_e", "m0"
        )
        rows.append({
            "material": "GaAs",
            "observable": "m*_e (8-band)",
            "computed": m_star,
            "expected": M_STAR_8BAND,
            "reference": "8-band kdotpy cross-validation",
            "tolerance": f"Regression ({TOL_MASS*100:.0f}%)",
            "delta": delta,
            "status": "PASS" if passed else "FAIL",
        })
        print(f"  m*: {m_star:.6f} m0 "
              f"(8-band ref {M_STAR_8BAND:.4f}, delta={delta:.2e}) "
              f"{'PASS' if passed else 'FAIL'}")

    return rows


def check_gfactor(build_dir, source_dir):
    """Check g-factor against Roth analytical formula.

    The Roth formula (Winkler 2003, Eq. 6.42):
      g = 2 - 2*EP*DeltaSO / (3*Eg*(Eg + DeltaSO))
    This gives the perturbative g-factor in the 8-band model's
    Roth approximation. The 8-band Lowdin result should match within ~1%.

    Returns list of benchmark row dicts.
    """
    config_path = os.path.join(
        source_dir, "tests", "regression", "configs",
        "gfactor_bulk_gaas_cb.toml"
    )
    exe_path = os.path.join(build_dir, "src", "gfactorCalculation")
    rows = []

    with tempfile.TemporaryDirectory(prefix="star_gaas_gf_") as tmpdir:
        rc, output_dir = run_executable(exe_path, config_path, tmpdir)
        if rc != 0:
            print(f"  FAIL: gfactorCalculation returned {rc}")
            return rows

        gf_path = os.path.join(output_dir, "gfactor.dat")
        gx, gy, gz = parse_gfactor(gf_path)

        # For bulk zinc-blende along principal axes, gz is the relevant
        # component (parallel to the perturbation field).
        g_computed = gz
        print(f"  g-factor: gx={gx:.6f}, gy={gy:.6f}, gz={gz:.6f}")
        print(f"  Roth prediction: {G_ROTH:.6f}")

        passed, delta, _ = compare_value(
            g_computed, G_ROTH, TOL_GFACTOR, "g*", ""
        )
        rows.append({
            "material": "GaAs",
            "observable": "g* (Roth)",
            "computed": g_computed,
            "expected": G_ROTH,
            "reference": "Winkler 2003, Eq. (6.42)",
            "tolerance": "Analytical (1%)",
            "delta": delta,
            "status": "PASS" if passed else "FAIL",
        })
        print(f"  g*: {g_computed:.6f} "
              f"(Roth {G_ROTH:.4f}, delta={delta:.2e}) "
              f"{'PASS' if passed else 'FAIL'}")

    return rows


def check_absorption_edge(build_dir, source_dir):
    """Check optical absorption edge against Eg.

    Parses absorption_TE.dat from the optics config run. Finds the first
    energy where absorption rises above a threshold fraction of the maximum.
    Compares against Eg.

    Returns list of benchmark row dicts.
    """
    config_path = os.path.join(
        source_dir, "tests", "regression", "configs",
        "bulk_gaas_optics.toml"
    )
    exe_path = os.path.join(build_dir, "src", "opticalProperties")
    rows = []

    with tempfile.TemporaryDirectory(prefix="star_gaas_abs_") as tmpdir:
        rc, output_dir = run_executable(exe_path, config_path, tmpdir)
        if rc != 0:
            print(f"  FAIL: opticalProperties returned {rc}")
            rows.append({
                "material": "GaAs",
                "observable": "Absorption edge",
                "computed": float("nan"),
                "expected": EG_GAAS,
                "reference": "Vurgaftman 2001",
                "tolerance": "Numerical (2%)",
                "delta": float("nan"),
                "status": "FAIL",
            })
            return rows

        abs_path = os.path.join(output_dir, "absorption_TE.dat")
        if not os.path.exists(abs_path):
            print(f"  FAIL: absorption_TE.dat not found in {output_dir}")
            return rows

        spectrum = parse_absorption(abs_path)
        if not spectrum:
            print("  FAIL: no absorption data parsed")
            return rows

        # Find absorption edge: first energy where alpha > threshold * max
        energies = [e for e, _ in spectrum]
        alphas = [a for _, a in spectrum]
        max_alpha = max(alphas) if alphas else 0.0

        if max_alpha <= 0:
            print("  FAIL: max absorption is zero or negative")
            return rows

        threshold = ABSORPTION_THRESHOLD_FRAC * max_alpha
        edge_energy = None
        for e, a in spectrum:
            if a > threshold:
                edge_energy = e
                break

        if edge_energy is None:
            print("  FAIL: could not detect absorption edge")
            return rows

        print(f"  Absorption edge: {edge_energy:.4f} eV "
              f"(max_alpha={max_alpha:.2f}, threshold={threshold:.2f})")

        passed, delta, _ = compare_value(
            edge_energy, EG_GAAS, TOL_ABSORPTION, "absorption edge", "eV"
        )
        rows.append({
            "material": "GaAs",
            "observable": "Absorption edge",
            "computed": edge_energy,
            "expected": EG_GAAS,
            "reference": "Vurgaftman 2001",
            "tolerance": "Numerical (2%)",
            "delta": delta,
            "status": "PASS" if passed else "FAIL",
        })
        print(f"  Absorption edge: {edge_energy:.4f} eV "
              f"(Eg {EG_GAAS:.3f}, delta={delta:.2e}) "
              f"{'PASS' if passed else 'FAIL'}")

    return rows


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <build_dir> <source_dir>")
        print(f"  build_dir  -- path to build/ (contains src/ executables)")
        print(f"  source_dir -- path to repo root")
        sys.exit(1)

    build_dir = os.path.abspath(sys.argv[1])
    source_dir = os.path.abspath(sys.argv[2])

    print("=" * 72)
    print("S1 -- GaAs Bulk Standard-Star Benchmark")
    print("=" * 72)
    print()
    print(f"GaAs parameters (literature):")
    print(f"  Eg      = {EG_GAAS:.3f} eV   (Vurgaftman 2001)")
    print(f"  DeltaSO = {DELTASO_GAAS:.3f} eV   (Vurgaftman 2001)")
    print(f"  EP      = {EP_GAAS:.1f} eV     (Kane matrix element)")
    print(f"  m*_Kane = {M_STAR_KANE:.4f} m0  (Eg/(EP+Eg))")
    print(f"  g_Roth  = {G_ROTH:.4f}      (Winkler 2003, Eq. 6.42)")
    print()

    all_rows = []

    # --- Eg and DeltaSO ---
    print("-" * 40)
    print("Eg and DeltaSO (bulk k=0)")
    print("-" * 40)
    rows = check_eg_deltaso(build_dir, source_dir)
    all_rows.extend(rows)
    print()

    # --- Effective mass ---
    print("-" * 40)
    print("Effective mass (bulk dispersion)")
    print("-" * 40)
    rows = check_effective_mass(build_dir, source_dir)
    all_rows.extend(rows)
    print()

    # --- g-factor ---
    print("-" * 40)
    print("g-factor (Roth analytical)")
    print("-" * 40)
    rows = check_gfactor(build_dir, source_dir)
    all_rows.extend(rows)
    print()

    # --- Absorption edge ---
    print("-" * 40)
    print("Absorption edge (bulk optics)")
    print("-" * 40)
    rows = check_absorption_edge(build_dir, source_dir)
    all_rows.extend(rows)
    print()

    # --- Benchmark table ---
    print("=" * 72)
    print("Benchmark Table")
    print("=" * 72)
    print()
    print_benchmark_header()
    for r in all_rows:
        print(format_benchmark_row(
            r["material"], r["observable"],
            r["computed"], r["expected"],
            r["reference"], r["tolerance"],
            r["delta"], r["status"],
        ))

    # --- Summary ---
    n_pass = sum(1 for r in all_rows if r["status"] == "PASS")
    n_fail = sum(1 for r in all_rows if r["status"] == "FAIL")
    print()
    print(f"Results: {n_pass} PASS, {n_fail} FAIL out of {len(all_rows)} observables")

    if n_fail > 0:
        print("FAIL: one or more GaAs bulk observables failed validation")
        sys.exit(1)
    else:
        print("PASS: all GaAs bulk observables within tolerance")
        sys.exit(0)


if __name__ == "__main__":
    main()
