#!/usr/bin/env python3
"""Lecture 06: Optical Properties -- Absorption, polarization, ISBT validation.

Runs the 8-band k.p optical solver (opticalProperties executable) for three
configurations and validates:

  1. Bulk GaAs absorption onset near Eg = 1.519 eV (within 5% numerical
     tolerance). Onset is detected via the inflection point (maximum d(alpha)/dE)
     of the absorption edge, which is robust against Lorentzian/Gaussian
     broadening tails below the band gap.
  2. QW GaAs/AlGaAs interband absorption: TE polarization dominant over TM
     for interband transitions (TE/TM ratio > 1 at peak absorption).
     Uses qw_optics_commutator.cfg which has a proper kx sweep with pure
     interband absorption (qw_gaas_algaas_optics.cfg is k0-only).
  3. QW intersubband (ISBT) absorption peak position validation.
  4. Overlay plot: absorption spectra with band-gap onset annotated and
     TE vs TM comparison, saved to docs/lecture/figures/.
"""
import os
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parent.parent
BUILD_DIR = REPO / "build"
CONFIGS_DIR = REPO / "tests" / "regression" / "configs"
FIGURES_DIR = REPO / "docs" / "lecture" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(REPO / "tests" / "integration"))
from star_helpers import (run_exe, parse_eigenvalues, parse_absorption,
                          compare_value, TOL_EXACT, TOL_ANALYTICAL,
                          TOL_NUMERICAL)

# ---------------------------------------------------------------------------
# Material parameters
# ---------------------------------------------------------------------------
# GaAs (Vurgaftman 2001 / Winkler 2003)
GAAS_EG = 1.519       # eV, band gap
GAAS_DELTA_SO = 0.341  # eV, spin-orbit splitting

# =========================================================================
# Section 1: Bulk GaAs absorption onset
# =========================================================================
def _find_inflection_onset(energies, absorption, e_center=None, e_window=0.3):
    """Find absorption onset using a threshold fraction within a search window.

    Locates the energy where absorption first exceeds 5% of the peak value
    within a restricted energy window centered on e_center. This avoids
    confusion from high-energy oscillations (van Hove singularities) and
    low-energy broadening tails.

    Args:
        energies: array of energy values (eV).
        absorption: array of absorption values (cm^-1).
        e_center: center of search window (eV). If None, uses midpoint.
        e_window: half-width of search window (eV).

    Returns:
        (onset_energy, onset_idx, peak_absorption_in_window) or None on failure.
    """
    if len(energies) < 5:
        return None

    # Determine search window
    if e_center is None:
        e_center = (energies[0] + energies[-1]) / 2.0
    e_lo = e_center - e_window
    e_hi = e_center + e_window

    # Restrict to the search window
    mask = (energies >= e_lo) & (energies <= e_hi)
    if np.sum(mask) < 5:
        return None

    e_win = energies[mask]
    a_win = absorption[mask]

    # Onset = energy where absorption first exceeds 5% of the peak
    # within the search window. This is robust against broadening tails
    # and high-energy oscillations.
    a_peak = a_win.max()
    if a_peak <= 0:
        return None

    threshold = 0.05 * a_peak
    onset_idx_win = np.where(a_win > threshold)[0]
    if len(onset_idx_win) == 0:
        return None

    onset_energy = e_win[onset_idx_win[0]]
    # Map back to the original index
    onset_idx = np.where(mask)[0][onset_idx_win[0]]

    return onset_energy, onset_idx, a_peak


def test_bulk_absorption_onset():
    """Verify bulk GaAs absorption onset near Eg = 1.519 eV."""
    print("=" * 60)
    print("Lecture 06 -- Section 1: Bulk GaAs absorption onset")
    print("=" * 60)

    cfg = CONFIGS_DIR / "bulk_gaas_optics.cfg"
    with tempfile.TemporaryDirectory() as work:
        rc, outdir = run_exe(str(BUILD_DIR), "opticalProperties",
                             str(cfg), work)
        if rc != 0:
            sys.exit(f"ERROR: opticalProperties returned {rc} for bulk GaAs")

        te_path = os.path.join(outdir, "absorption_TE.dat")
        if not os.path.isfile(te_path):
            sys.exit(f"ERROR: absorption_TE.dat not found at {te_path}")

        data = parse_absorption(te_path)
        if not data:
            sys.exit("ERROR: no absorption data parsed from bulk GaAs TE")

    energies = np.array([d[0] for d in data])
    absorption = np.array([d[1] for d in data])

    print(f"  Parsed {len(data)} energy points: "
          f"E = [{energies[0]:.3f}, {energies[-1]:.3f}] eV")
    print(f"  Absorption range: [{absorption.min():.4e}, "
          f"{absorption.max():.4e}] cm^-1")

    # Find onset via inflection point, searching in a window around the
    # expected band gap. The 8-band joint-DOS produces oscillations at
    # higher energies, so we restrict the search to avoid those.
    result = _find_inflection_onset(energies, absorption,
                                    e_center=GAAS_EG, e_window=0.3)
    if result is None:
        sys.exit("ERROR: inflection-point onset detection failed")

    onset_energy, infl_idx, a_peak_win = result
    print(f"  Peak absorption in Eg window: {a_peak_win:.2f} cm^-1")
    print(f"  5% threshold onset: E = {onset_energy:.4f} eV")
    print(f"  Expected Eg:        {GAAS_EG:.4f} eV")

    passed, delta, _ = compare_value(
        onset_energy, GAAS_EG, TOL_NUMERICAL,
        "Bulk GaAs absorption onset (inflection)", unit="eV",
    )
    status = "PASS" if passed else "FAIL"
    print(f"  {status}: delta = {delta:.4f}  (tolerance {TOL_NUMERICAL:.0%})")

    if passed:
        print("  --> Absorption onset matches GaAs band gap within 5%.")
    else:
        print("  --> FAIL: onset energy outside 5% tolerance.")
    print()
    return passed, energies, absorption


# =========================================================================
# Section 2: QW interband absorption -- TE vs TM polarization
# =========================================================================
def test_qw_te_tm_polarization():
    """Verify TE polarization dominant for QW interband transitions."""
    print("=" * 60)
    print("Lecture 06 -- Section 2: QW TE vs TM polarization")
    print("=" * 60)

    # qw_gaas_algaas_optics.cfg is k0-only (waveVectorStep=0) and produces
    # zero absorption for QW interband transitions. The _full variant has a
    # k-sweep but enables gain+ISBT which distorts the TE/TM ratio.
    # Use qw_optics_commutator.cfg which has a proper kx sweep with pure
    # interband absorption (no gain, no ISBT) for a clean TE vs TM comparison.
    cfg = CONFIGS_DIR / "qw_optics_commutator.cfg"
    with tempfile.TemporaryDirectory() as work:
        rc, outdir = run_exe(str(BUILD_DIR), "opticalProperties",
                             str(cfg), work)
        if rc != 0:
            sys.exit(f"ERROR: opticalProperties returned {rc} for QW GaAs/AlGaAs")

        te_path = os.path.join(outdir, "absorption_TE.dat")
        tm_path = os.path.join(outdir, "absorption_TM.dat")
        if not os.path.isfile(te_path):
            sys.exit(f"ERROR: absorption_TE.dat not found at {te_path}")
        if not os.path.isfile(tm_path):
            sys.exit(f"ERROR: absorption_TM.dat not found at {tm_path}")

        te_data = parse_absorption(te_path)
        tm_data = parse_absorption(tm_path)
        if not te_data:
            sys.exit("ERROR: no TE absorption data parsed")
        if not tm_data:
            sys.exit("ERROR: no TM absorption data parsed")

    te_energies = np.array([d[0] for d in te_data])
    te_alpha = np.array([d[1] for d in te_data])
    tm_energies = np.array([d[0] for d in tm_data])
    tm_alpha = np.array([d[1] for d in tm_data])

    te_peak = te_alpha.max()
    tm_peak = tm_alpha.max()
    print(f"  TE peak absorption: {te_peak:.4e} cm^-1")
    print(f"  TM peak absorption: {tm_peak:.4e} cm^-1")

    if tm_peak > 0:
        ratio = te_peak / tm_peak
    else:
        ratio = float("inf")

    print(f"  TE/TM peak ratio:  {ratio:.2f}")

    # For interband transitions in a QW, TE (in-plane) should dominate
    # over TM (out-of-plane). Check ratio > 1.
    passed = ratio > 1.0
    status = "PASS" if passed else "FAIL"
    print(f"  {status}: TE/TM ratio {'>' if passed else '<='} 1  "
          f"(interband TE dominant)")

    if passed:
        print("  --> TE polarization correctly dominant for interband transitions.")
    else:
        print("  --> FAIL: TM unexpectedly dominant over TE.")
    print()
    return passed, te_energies, te_alpha, tm_energies, tm_alpha


# =========================================================================
# Section 3: ISBT absorption peak
# =========================================================================
def test_isbt_peak():
    """Verify QW intersubband absorption peak position."""
    print("=" * 60)
    print("Lecture 06 -- Section 3: ISBT absorption peak")
    print("=" * 60)

    cfg = CONFIGS_DIR / "qw_gaas_algaas_isbt.cfg"
    with tempfile.TemporaryDirectory() as work:
        rc, outdir = run_exe(str(BUILD_DIR), "opticalProperties",
                             str(cfg), work)
        if rc != 0:
            sys.exit(f"ERROR: opticalProperties returned {rc} for ISBT config")

        isbt_path = os.path.join(outdir, "absorption_ISBT.dat")
        if not os.path.isfile(isbt_path):
            sys.exit(f"ERROR: absorption_ISBT.dat not found at {isbt_path}")

        isbt_data = parse_absorption(isbt_path)
        if not isbt_data:
            sys.exit("ERROR: no ISBT absorption data parsed")

    energies = np.array([d[0] for d in isbt_data])
    alpha = np.array([d[1] for d in isbt_data])

    print(f"  Parsed {len(isbt_data)} energy points: "
          f"E = [{energies[0]:.4f}, {energies[-1]:.4f}] eV")

    peak_idx = np.argmax(alpha)
    peak_energy = energies[peak_idx]
    peak_alpha = alpha[peak_idx]
    print(f"  ISBT peak: {peak_energy:.4f} eV  (alpha = {peak_alpha:.4e} cm^-1)")

    # The ISBT config uses a GaAs/Al30Ga70As QW with FDstep=201, FDorder=4.
    # Subband spacing depends on well width and confinement. For a 100 A GaAs
    # well, the e1-e2 ISBT transition is approximately 0.08-0.15 eV.
    # Validate that (1) the peak is in the expected ISBT energy range and
    # (2) the peak absorption is well above numerical noise (> 1 cm^-1).
    isbt_e_min, isbt_e_max = 0.03, 0.30  # eV, generous ISBT range
    min_alpha = 1.0  # cm^-1, above numerical noise floor
    in_range = isbt_e_min <= peak_energy <= isbt_e_max
    above_noise = peak_alpha > min_alpha

    if in_range and above_noise:
        print(f"  ISBT peak validated: E={peak_energy:.4f} eV in [{isbt_e_min}, {isbt_e_max}], "
              f"alpha={peak_alpha:.4e} cm^-1 > {min_alpha} cm^-1")
        passed = True
    elif not in_range:
        print(f"  FAIL: ISBT peak at {peak_energy:.4f} eV outside expected "
              f"range [{isbt_e_min}, {isbt_e_max}] eV")
        passed = False
    else:
        print(f"  FAIL: ISBT peak alpha={peak_alpha:.4e} cm^-1 below "
              f"noise threshold {min_alpha} cm^-1")
        passed = False

    status = "PASS" if passed else "FAIL"
    print(f"  {status}: ISBT peak at {peak_energy:.4f} eV")
    print()
    return passed, energies, alpha


# =========================================================================
# Section 4: Overlay plot
# =========================================================================
def plot_absorption_spectra(bulk_e, bulk_alpha,
                            qw_te_e, qw_te_alpha, qw_tm_e, qw_tm_alpha,
                            isbt_e, isbt_alpha):
    """Generate overlay absorption plot with band-gap onset annotation."""
    print("=" * 60)
    print("Lecture 06 -- Section 4: Absorption overlay plot")
    print("=" * 60)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # --- Panel 1: Bulk GaAs absorption with onset ---
    ax1 = axes[0]
    ax1.plot(bulk_e, bulk_alpha, "b-", linewidth=1.5, label=r"$\alpha_{\rm TE}$")
    ax1.axvline(GAAS_EG, color="red", ls="--", lw=1.2,
                label=f"$E_g$ = {GAAS_EG:.3f} eV")

    # Annotate onset via inflection point
    alpha_max_bulk = bulk_alpha.max()
    result = _find_inflection_onset(bulk_e, bulk_alpha,
                                    e_center=GAAS_EG, e_window=0.3)
    if result is not None:
        onset_e, infl_idx, _ = result
        ax1.annotate(
            f"Onset = {onset_e:.3f} eV",
            xy=(onset_e, bulk_alpha[infl_idx]),
            xytext=(onset_e + 0.2, alpha_max_bulk * 0.5),
            fontsize=9,
            arrowprops=dict(arrowstyle="->", color="gray", lw=0.8),
            color="gray",
        )

    ax1.set_xlabel("Energy (eV)", fontsize=11)
    ax1.set_ylabel(r"Absorption (cm$^{-1}$)", fontsize=11)
    ax1.set_title("Bulk GaAs Absorption", fontsize=12)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    # --- Panel 2: QW TE vs TM ---
    ax2 = axes[1]
    ax2.plot(qw_te_e, qw_te_alpha, "b-", linewidth=1.5,
             label=r"TE ($\alpha_{\rm TE}$)")
    ax2.plot(qw_tm_e, qw_tm_alpha, "r--", linewidth=1.5,
             label=r"TM ($\alpha_{\rm TM}$)")

    te_peak = qw_te_alpha.max()
    tm_peak = qw_tm_alpha.max()
    if tm_peak > 0:
        ratio = te_peak / tm_peak
    else:
        ratio = float("inf")
    ax2.set_xlabel("Energy (eV)", fontsize=11)
    ax2.set_ylabel(r"Absorption (cm$^{-1}$)", fontsize=11)
    ax2.set_title(f"QW GaAs/AlGaAs  TE/TM = {ratio:.1f}", fontsize=12)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    # --- Panel 3: ISBT ---
    ax3 = axes[2]
    ax3.plot(isbt_e, isbt_alpha, "g-", linewidth=1.5,
             label=r"ISBT ($\alpha_{\rm TM}$)")

    if isbt_alpha.max() > 0:
        peak_idx = np.argmax(isbt_alpha)
        peak_e = isbt_e[peak_idx]
        ax3.axvline(peak_e, color="gray", ls=":", lw=1.0,
                    label=f"Peak = {peak_e:.4f} eV")

    ax3.set_xlabel("Energy (eV)", fontsize=11)
    ax3.set_ylabel(r"Absorption (cm$^{-1}$)", fontsize=11)
    ax3.set_title("QW Intersubband (ISBT)", fontsize=12)
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3)

    fig.suptitle("Lecture 06: Optical Properties -- Absorption Spectra",
                 fontsize=14, fontweight="bold", y=1.02)
    fig.tight_layout()

    out_path = FIGURES_DIR / "lecture_06_absorption.png"
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"  Plot saved to {out_path}")
    print()


# =========================================================================
# Main
# =========================================================================
def main():
    print("\n" + "=" * 60)
    print("  LECTURE 06: Optical Properties Validation")
    print("  Absorption, polarization, intersubband transitions")
    print("=" * 60 + "\n")

    # Run all sections
    s1_pass, bulk_e, bulk_alpha = test_bulk_absorption_onset()
    s2_pass, qw_te_e, qw_te_alpha, qw_tm_e, qw_tm_alpha = \
        test_qw_te_tm_polarization()
    s3_pass, isbt_e, isbt_alpha = test_isbt_peak()

    # Generate overlay plot
    plot_absorption_spectra(bulk_e, bulk_alpha,
                            qw_te_e, qw_te_alpha, qw_tm_e, qw_tm_alpha,
                            isbt_e, isbt_alpha)

    # Summary
    print("=" * 60)
    print("  SUMMARY")
    print("=" * 60)
    results = [
        ("Section 1: Bulk GaAs absorption onset", s1_pass),
        ("Section 2: QW TE/TM polarization", s2_pass),
        ("Section 3: ISBT peak detection", s3_pass),
        ("Section 4: Overlay plot", True),
    ]
    all_pass = True
    for label, passed in results:
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {label}")
        all_pass = all_pass and passed
    print()

    if all_pass:
        print("  All validations passed.")
    else:
        print("  Some validations FAILED.")
        sys.exit(1)


if __name__ == "__main__":
    main()
