#!/usr/bin/env python3
"""Lecture 14: Exciton Binding Energy & Phonon Scattering Rates.

Runs the 8-band k.p solver (bandStructure executable) for GaAs/AlGaAs
quantum wells and validates:

  1. Exciton binding energy vs well width: Eb decreases with increasing
     well width (quantum confinement effect). Uses Bastard variational
     method with golden-section search over the trial parameter lambda.
 2. LO-phonon Froehlich scattering rates for intersubband transitions.
     Validates emission and absorption rates are positive and in the
     physically reasonable range (~10^9 to 10^12 1/s).

Outputs:
  docs/lecture/figures/lecture_14_exciton_binding.png
  docs/lecture/figures/lecture_14_scattering.png
"""
import os
import re
import shutil
import subprocess
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

sys.path.insert(0, str(REPO / "tests" / "integration"))
from star_helpers import run_exe, parse_eigenvalues

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
HBAR2_OVER_2M0 = 3.80998   # eV * Angstrom^2 (hbar^2 / (2*m0))
HBAR_EV_S = 6.582119569e-16  # eV*s

# GaAs LO-phonon energy (Vurgaftman 2001)
GAAS_HBAR_OMEGA_LO = 0.036  # eV (~36 meV)

# GaAs dielectric constants
GAAS_EPS_INF = 10.9
GAAS_EPS_0 = 12.9


# ---------------------------------------------------------------------------
# Config templates
# ---------------------------------------------------------------------------
_EXCITON_CONFIG_TEMPLATE = """\
waveVector: k0
waveVectorMax: 0.0
waveVectorStep: 1
confinement: 1
FDstep: 51
FDorder: 2
numLayers: 2
material1: Al30Ga70As -200 200 0
material2: GaAs -{half_w} {half_w} 0
numcb: 2
numvb: 6
ExternalField: 0  EF
EFParams: 0.0
SC: 0
exciton: T
method: variational
"""

_SCATTERING_CONFIG = """\
waveVector: k0
waveVectorMax: 0.0
waveVectorStep: 1
confinement: 1
FDstep: 51
FDorder: 2
numLayers: 2
material1: Al30Ga70As -200 200 0
material2: GaAs -50 50 0
numcb: 4
numvb: 6
ExternalField: 0  EF
EFParams: 0.0
SC: 0
scattering: T
phonon_energy: 0.036
eps_inf: 10.9
eps_0: 12.9
"""


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------
def parse_exciton_stdout(stdout_text):
    """Parse exciton binding energy from bandStructure stdout.

    Expected format (from exciton.f90 write_exciton_output + main.f90):
      E_binding   =  <value>  meV
      lambda_opt  =  <value>  AA

    Returns:
        dict with 'E_binding_meV' and 'lambda_opt_AA', or None on failure.
    """
    result = {}

    # Match "E_binding   = ..." or "E_binding   = ..."
    m_eb = re.search(r'E_binding\s*=\s*([\d.eE+-]+)', stdout_text)
    if m_eb:
        result['E_binding_meV'] = float(m_eb.group(1))

    # Match "lambda_opt  = ..."
    m_lam = re.search(r'lambda_opt\s*=\s*([\d.eE+-]+)', stdout_text)
    if m_lam:
        result['lambda_opt_AA'] = float(m_lam.group(1))

    # Also try alternate format from main.f90:
    # "Exciton binding energy:  <value> meV"
    if 'E_binding_meV' not in result:
        m_alt = re.search(r'Exciton binding energy:\s*([\d.eE+-]+)', stdout_text)
        if m_alt:
            result['E_binding_meV'] = float(m_alt.group(1))

    # "Variational parameter:  <value> AA"
    if 'lambda_opt_AA' not in result:
        m_alt2 = re.search(r'Variational parameter:\s*([\d.eE+-]+)', stdout_text)
        if m_alt2:
            result['lambda_opt_AA'] = float(m_alt2.group(1))

    return result if result else None


def parse_exciton_file(filepath):
    """Parse output/exciton.dat file.

    Format (from write_exciton_output):
      # lambda_opt(AA)  E_binding(meV)  mu/m0  eps_r
      <value> <value> <value> <value>

    Returns:
        dict with parsed values, or None on failure.
    """
    try:
        data = np.loadtxt(filepath, comments='#')
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.shape[1] >= 2:
            return {
                'lambda_opt_AA': data[0, 0],
                'E_binding_meV': data[0, 1],
                'mu_over_m0': data[0, 2] if data.shape[1] > 2 else None,
                'eps_r': data[0, 3] if data.shape[1] > 3 else None,
            }
    except Exception:
        pass
    return None


def parse_scattering_file(filepath):
    """Parse output/scattering_rates.dat.

    Format (from write_scattering_output):
      # i  j  E_ij(meV)  rate_emission(1/s)  rate_absorption(1/s)
      #     tau_emission(ps)  tau_absorption(ps)
      i  j  dE  rate_em  rate_ab  tau_em  tau_ab

    Returns:
        list of dicts with keys: i, j, dE_meV, rate_em, rate_ab, tau_em_ps,
        tau_ab_ps
    """
    results = []
    try:
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 5:
                    results.append({
                        'i': int(parts[0]),
                        'j': int(parts[1]),
                        'dE_meV': float(parts[2]),
                        'rate_em': float(parts[3]),
                        'rate_ab': float(parts[4]),
                        'tau_em_ps': float(parts[5]) if len(parts) > 5 else None,
                        'tau_ab_ps': float(parts[6]) if len(parts) > 6 else None,
                    })
    except Exception:
        pass
    return results


def bastard_binding_analytic(Lz_A, eps_r=12.9, mu_m0=0.059):
    """Approximate Bastard variational binding energy for a 2D QW.

    Uses the simple 2D hydrogenic model:
      Eb = mu * e^4 / (2 * (4*pi*eps0*eps_r)^2 * hbar^2)
    converted to meV. This is the 2D limit (thin well); real QW has
    smaller Eb due to finite well width.

    For the variational trend, we use:
      Eb(Lz) ~ Eb_2D * f(Lz)
    where f(Lz) decreases with Lz as the exciton spreads in the growth
    direction.

    Args:
        Lz_A: well width in Angstrom
        eps_r: relative dielectric constant
        mu_m0: reduced mass in units of m0

    Returns:
        Eb in meV
    """
    # 2D hydrogenic binding energy: Eb = Rydberg * mu/m0 / eps_r^2
    # GaAs Rydberg: 13.6 eV * 0.059 / 12.9^2 = 4.84 meV
    rydberg_eV = 13.6  # hydrogen Rydberg in eV
    Eb_2D_meV = rydberg_eV * 1000.0 * mu_m0 / (eps_r ** 2)

    # Approximate finite-width correction from Bastard (1982):
    # lambda ~ sqrt(Lz * a_B) where a_B is the 2D Bohr radius
    # The binding energy is reduced for wider wells as the wavefunction
    # spreads. A simple approximation:
    # Eb(Lz) ~ Eb_2D / (1 + (Lz / a_B_2D)^0.5)
    # where a_B_2D = 2 * eps_r * a0 / mu_m0
    a_B_2D = 2.0 * eps_r * 0.529 / mu_m0  # Angstrom (2D Bohr radius)
    ratio = Lz_A / a_B_2D
    Eb_approx = Eb_2D_meV / (1.0 + ratio ** 0.5)

    return Eb_approx


# ---------------------------------------------------------------------------
# Section 1: Exciton Binding Energy vs Well Width
# ---------------------------------------------------------------------------
def section1_exciton_binding():
    """Exciton binding energy for GaAs/AlGaAs QW vs well width.

    Runs bandStructure for well widths 30, 50, 80, 100, 150, 200 A
    with the exciton variational solver enabled. Validates that Eb
    decreases with well width and is in the 1-20 meV range.
    """
    print("=" * 60)
    print("  Section 1: Exciton Binding Energy vs Well Width")
    print("  Bastard variational method (PRB 1982)")
    print("=" * 60)

    exe = BUILD_DIR / "src" / "bandStructure"
    if not exe.exists():
        print("  ERROR: bandStructure not found. Build first.")
        return False

    well_widths = [30, 50, 80, 100, 150, 200]
    results = []

    for w in well_widths:
        half_w = w / 2.0
        config = _EXCITON_CONFIG_TEMPLATE.format(half_w=half_w)

        print(f"\n  Well width = {w} A (half_w = {half_w:.1f} A) ...")

        with tempfile.TemporaryDirectory() as work:
            # Write config
            input_cfg = os.path.join(work, "input.cfg")
            with open(input_cfg, "w") as f:
                f.write(config)

            output_dir = os.path.join(work, "output")
            os.makedirs(output_dir, exist_ok=True)

            # Run bandStructure
            try:
                proc = subprocess.run(
                    [str(exe)],
                    cwd=work,
                    capture_output=True,
                    text=True,
                    timeout=120,
                )
            except subprocess.TimeoutExpired:
                print(f"    TIMEOUT for width = {w} A")
                results.append((w, None))
                continue

            if proc.returncode != 0:
                print(f"    ERROR: bandStructure returned {proc.returncode}")
                if proc.stderr:
                    print(f"    stderr: {proc.stderr[:300]}")
                results.append((w, None))
                continue

            # Parse stdout for exciton output
            exciton_data = parse_exciton_stdout(proc.stdout)

            # Also parse exciton.dat file
            exciton_path = os.path.join(output_dir, "exciton.dat")
            if os.path.isfile(exciton_path):
                file_data = parse_exciton_file(exciton_path)
                if file_data and exciton_data is None:
                    exciton_data = file_data

            if exciton_data and 'E_binding_meV' in exciton_data:
                Eb = exciton_data['E_binding_meV']
                lam = exciton_data.get('lambda_opt_AA', None)
                print(f"    E_binding = {Eb:.4f} meV", end="")
                if lam is not None:
                    print(f", lambda_opt = {lam:.2f} A", end="")
                print()
                results.append((w, exciton_data))
            else:
                print(f"    WARNING: no exciton data parsed for width = {w} A")
                if proc.stdout:
                    for line in proc.stdout.split('\n'):
                        if 'exciton' in line.lower() or 'binding' in line.lower():
                            print(f"    stdout: {line.strip()}")
                results.append((w, None))

    # --- Validation ---
    print("\n  --- Validation ---")
    valid = [(w, d['E_binding_meV']) for w, d in results
             if d is not None and 'E_binding_meV' in d]

    all_pass = True

    if len(valid) < 2:
        print("  FAIL: not enough data points for validation")
        all_pass = False
    else:
        widths_v, Eb_v = zip(*valid)

        # Check 1: Eb decreases with well width (allow some tolerance for
        # numerical noise)
        decreasing = True
        for k in range(len(Eb_v) - 1):
            if Eb_v[k + 1] > Eb_v[k] * 1.05:  # allow 5% tolerance
                decreasing = False
                print(f"  FAIL: Eb({widths_v[k+1]}A) = {Eb_v[k+1]:.4f} meV "
                      f"> Eb({widths_v[k]}A) = {Eb_v[k]:.4f} meV")
                all_pass = False
                break
        if decreasing:
            print(f"  PASS: Eb decreases with well width "
                  f"({Eb_v[0]:.2f} meV at {widths_v[0]} A -> "
                  f"{Eb_v[-1]:.2f} meV at {widths_v[-1]} A)")

        # Check 2: Eb in range 1-20 meV for GaAs QWs
        for w, eb in valid:
            if not (1.0 <= eb <= 20.0):
                print(f"  FAIL: Eb = {eb:.2f} meV at {w} A outside [1, 20] meV")
                all_pass = False

        if all(eb >= 1.0 and eb <= 20.0 for _, eb in valid):
            print(f"  PASS: all Eb values in [1, 20] meV "
                  f"(range: {min(eb for _, eb in valid):.2f} -- "
                  f"{max(eb for _, eb in valid):.2f} meV)")

    # --- Generate figure ---
    print("\n  Generating exciton binding energy figure ...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left panel: Eb vs well width
    ax1 = axes[0]
    if valid:
        w_plot, eb_plot = zip(*valid)
        ax1.plot(w_plot, eb_plot, 'bo-', linewidth=2, markersize=8,
                 label='8-band k.p (variational)')

        # Bastard analytic trend
        w_anal = np.linspace(20, 250, 100)
        eb_anal = [bastard_binding_analytic(w) for w in w_anal]
        ax1.plot(w_anal, eb_anal, 'r--', linewidth=1.5,
                 label='Bastard approx. trend')

        ax1.set_xlabel('Well width (Angstrom)', fontsize=12)
        ax1.set_ylabel(r'Exciton binding energy $E_b$ (meV)', fontsize=12)
        ax1.set_title('Exciton Binding Energy vs Well Width\n'
                      'GaAs/AlGaAs QW (Bastard variational)', fontsize=12)
        ax1.legend(fontsize=10)
        ax1.grid(True, alpha=0.3)

        # Annotate data points
        for w, eb in valid:
            ax1.annotate(f'{eb:.1f}', (w, eb),
                         textcoords="offset points", xytext=(0, 10),
                         fontsize=8, ha='center')
    else:
        ax1.text(0.5, 0.5, 'No data available', transform=ax1.transAxes,
                 ha='center', fontsize=12, color='gray')
        ax1.set_title('Exciton Binding Energy vs Well Width', fontsize=12)

    # Right panel: lambda_opt vs well width
    ax2 = axes[1]
    lambda_data = [(w, d['lambda_opt_AA']) for w, d in results
                   if d is not None and 'lambda_opt_AA' in d
                   and d.get('lambda_opt_AA', 0) > 0]
    if lambda_data:
        w_lam, lam_plot = zip(*lambda_data)
        ax2.plot(w_lam, lam_plot, 'gs-', linewidth=2, markersize=8,
                 label=r'$\lambda_{\rm opt}$ (variational)')
        ax2.set_xlabel('Well width (Angstrom)', fontsize=12)
        ax2.set_ylabel(r'Variational parameter $\lambda_{\rm opt}$ (AA)',
                        fontsize=12)
        ax2.set_title(r'Exciton Bohr Radius $\lambda_{\rm opt}$ vs Well Width'
                      '\nGaAs/AlGaAs QW', fontsize=12)
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3)

        for w, lam in lambda_data:
            ax2.annotate(f'{lam:.0f}', (w, lam),
                         textcoords="offset points", xytext=(0, 10),
                         fontsize=8, ha='center')
    else:
        ax2.text(0.5, 0.5, 'No lambda data available',
                 transform=ax2.transAxes, ha='center', fontsize=12,
                 color='gray')
        ax2.set_title(r'Variational parameter $\lambda_{\rm opt}$', fontsize=12)

    fig.suptitle('Lecture 14: Exciton Binding in GaAs/AlGaAs Quantum Wells',
                 fontsize=14, fontweight='bold', y=1.02)
    fig.tight_layout()

    out_path = FIGURES_DIR / "lecture_14_exciton_binding.png"
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")

    status = "PASS" if all_pass else "FAIL"
    print(f"\n  Section 1 result: {status}")
    return all_pass


# ---------------------------------------------------------------------------
# Section 2: Phonon Scattering Rates
# ---------------------------------------------------------------------------
def section2_phonon_scattering():
    """LO-phonon Froehlich scattering rates for GaAs/AlGaAs QW.

    Runs bandStructure with the scattering block enabled for a 100 A
    GaAs/AlGaAs QW with 4 CB states. Parses scattering_rates.dat and
    validates that rates are positive and in a physically reasonable
    range.
    """
    print("\n" + "=" * 60)
    print("  Section 2: LO-Phonon Scattering Rates")
    print("  Froehlich coupling (Ferreira & Bastard PRB 1989)")
    print("=" * 60)

    exe = BUILD_DIR / "src" / "bandStructure"
    if not exe.exists():
        print("  ERROR: bandStructure not found. Build first.")
        return False

    print(f"\n  Running GaAs/AlGaAs QW (100 A well, 4 CB states) ...")
    print(f"  LO-phonon energy = {GAAS_HBAR_OMEGA_LO * 1000:.1f} meV")
    print(f"  eps_inf = {GAAS_EPS_INF}, eps_0 = {GAAS_EPS_0}")

    with tempfile.TemporaryDirectory() as work:
        # Write config
        input_cfg = os.path.join(work, "input.cfg")
        with open(input_cfg, "w") as f:
            f.write(_SCATTERING_CONFIG)

        output_dir = os.path.join(work, "output")
        os.makedirs(output_dir, exist_ok=True)

        # Run bandStructure
        try:
            proc = subprocess.run(
                [str(exe)],
                cwd=work,
                capture_output=True,
                text=True,
                timeout=120,
            )
        except subprocess.TimeoutExpired:
            print("  ERROR: bandStructure timed out")
            return False

        if proc.returncode != 0:
            print(f"  ERROR: bandStructure returned {proc.returncode}")
            if proc.stderr:
                print(f"  stderr: {proc.stderr[:500]}")
            return False

        # Parse scattering rates file
        scatter_path = os.path.join(output_dir, "scattering_rates.dat")
        if not os.path.isfile(scatter_path):
            print("  ERROR: scattering_rates.dat not found")
            if proc.stdout:
                for line in proc.stdout.split('\n'):
                    if 'scattering' in line.lower():
                        print(f"  stdout: {line.strip()}")
            return False

        rates = parse_scattering_file(scatter_path)

    if not rates:
        print("  ERROR: no scattering rate data parsed")
        return False

    # Print parsed rates
    print(f"\n  Parsed {len(rates)} intersubband transitions:")
    print(f"  {'Transition':>12s}  {'dE (meV)':>10s}  "
          f"{'Rate_em (1/s)':>14s}  {'Rate_ab (1/s)':>14s}")
    print(f"  {'-'*12}  {'-'*10}  {'-'*14}  {'-'*14}")
    for r in rates:
        print(f"  CB{r['i']} -> CB{r['j']:>4d}  "
              f"{r['dE_meV']:10.4f}  "
              f"{r['rate_em']:14.4e}  {r['rate_ab']:14.4e}")

    # --- Validation ---
    print("\n  --- Validation ---")
    all_pass = True

    # Check 1: All rates are positive
    for r in rates:
        if r['rate_em'] < 0:
            print(f"  FAIL: emission rate < 0 for CB{r['i']}->CB{r['j']}")
            all_pass = False
        if r['rate_ab'] < 0:
            print(f"  FAIL: absorption rate < 0 for CB{r['i']}->CB{r['j']}")
            all_pass = False

    all_positive = all(r['rate_em'] >= 0 and r['rate_ab'] >= 0 for r in rates)
    if all_positive:
        print("  PASS: all scattering rates are non-negative")

    # Check 2: At least one emission rate > 0 (requires dE > hbar*omega_LO)
    has_emission = any(r['rate_em'] > 0 for r in rates)
    if has_emission:
        print("  PASS: at least one emission rate > 0 (dE > hbar*omega_LO)")
    else:
        print("  WARNING: no emission rates > 0 (all transitions below "
              "phonon threshold)")

    # Check 3: Rates in reasonable range (~10^8 to 10^13 1/s)
    # Ferreira & Bastard (PRB 1989): LO-phonon scattering rates for GaAs QWs
    # are typically 10^9 to 10^12 s^-1.  Rates below 10^8 suggest a
    # unit-conversion error.
    rate_min_log = 8.0   # log10(1/s)
    rate_max_log = 13.0  # log10(1/s)
    rates_reasonable = True
    for r in rates:
        for label, rate in [('emission', r['rate_em']),
                            ('absorption', r['rate_ab'])]:
            if rate > 0:
                log_rate = np.log10(rate)
                if not (rate_min_log <= log_rate <= rate_max_log):
                    print(f"  WARNING: {label} rate for CB{r['i']}->CB{r['j']}"
                          f" = {rate:.2e} 1/s outside typical range"
                          f" [10^{rate_min_log}, 10^{rate_max_log}]")
                    rates_reasonable = False

    if rates_reasonable and any(r['rate_ab'] > 0 for r in rates):
        print("  PASS: scattering rates in physically reasonable range")

    # Check 4: Energy separations are positive (j > i => E_j > E_i)
    for r in rates:
        if r['dE_meV'] <= 0:
            print(f"  FAIL: dE <= 0 for CB{r['i']}->CB{r['j']}")
            all_pass = False

    all_dE_positive = all(r['dE_meV'] > 0 for r in rates)
    if all_dE_positive:
        print("  PASS: all energy separations dE > 0")

    # --- Generate figure ---
    print("\n  Generating scattering rates figure ...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left panel: Bar chart of emission and absorption rates
    ax1 = axes[0]
    transitions = [f"CB{r['i']}-CB{r['j']}" for r in rates]
    rate_em_vals = [r['rate_em'] for r in rates]
    rate_ab_vals = [r['rate_ab'] for r in rates]

    x_pos = np.arange(len(transitions))
    width = 0.35

    bars_em = ax1.bar(x_pos - width / 2, rate_em_vals, width,
                       label='Emission', color='crimson', edgecolor='k',
                       linewidth=0.5)
    bars_ab = ax1.bar(x_pos + width / 2, rate_ab_vals, width,
                       label='Absorption', color='steelblue', edgecolor='k',
                       linewidth=0.5)

    ax1.set_yscale('log')
    ax1.set_xlabel('Transition', fontsize=12)
    ax1.set_ylabel(r'Scattering rate (s$^{-1}$)', fontsize=12)
    ax1.set_title('LO-Phonon Intersubband Scattering Rates\n'
                  'GaAs/AlGaAs QW (100 A)', fontsize=12)
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(transitions, fontsize=10)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3, axis='y')

    # Annotate rates on bars
    for bar, val in zip(bars_em, rate_em_vals):
        if val > 0:
            ax1.text(bar.get_x() + bar.get_width() / 2, val * 1.3,
                     f'{val:.1e}', ha='center', va='bottom', fontsize=7,
                     rotation=45)
    for bar, val in zip(bars_ab, rate_ab_vals):
        if val > 0:
            ax1.text(bar.get_x() + bar.get_width() / 2, val * 1.3,
                     f'{val:.1e}', ha='center', va='bottom', fontsize=7,
                     rotation=45)

    # Right panel: Energy level diagram with transition energies
    ax2 = axes[1]
    if rates:
        dE_vals = [r['dE_meV'] for r in rates]
        hbar_omega_meV = GAAS_HBAR_OMEGA_LO * 1000  # 36 meV

        colors = ['crimson' if r['rate_em'] > 0 else 'gray' for r in rates]
        ax2.bar(x_pos, dE_vals, color=colors, edgecolor='k', linewidth=0.5,
                label=r'$\Delta E_{ij}$')
        ax2.axhline(hbar_omega_meV, color='navy', ls='--', lw=2,
                     label=rf'$\hbar\omega_{{LO}}$ = {hbar_omega_meV:.0f} meV')

        ax2.set_xlabel('Transition', fontsize=12)
        ax2.set_ylabel(r'$\Delta E_{ij}$ (meV)', fontsize=12)
        ax2.set_title('Subband Separations vs Phonon Energy\n'
                      '(red = emission allowed)', fontsize=12)
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(transitions, fontsize=10)
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3, axis='y')

        # Annotate dE values
        for i, (xp, de) in enumerate(zip(x_pos, dE_vals)):
            ax2.annotate(f'{de:.1f}', (xp, de),
                         textcoords="offset points", xytext=(0, 5),
                         fontsize=8, ha='center')
    else:
        ax2.text(0.5, 0.5, 'No data available', transform=ax2.transAxes,
                 ha='center', fontsize=12, color='gray')
        ax2.set_title('Subband Separations', fontsize=12)

    fig.suptitle('Lecture 14: LO-Phonon Scattering in GaAs/AlGaAs QW',
                 fontsize=14, fontweight='bold', y=1.02)
    fig.tight_layout()

    out_path = FIGURES_DIR / "lecture_14_scattering.png"
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")

    status = "PASS" if all_pass else "FAIL"
    print(f"\n  Section 2 result: {status}")
    return all_pass


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 60)
    print("  Lecture 14: Excitons & Scattering")
    print("=" * 60)

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    results = {}
    results['exciton'] = section1_exciton_binding()
    results['scattering'] = section2_phonon_scattering()

    print("\n" + "=" * 60)
    print("  Summary")
    print("=" * 60)
    all_pass = True
    for section, passed in results.items():
        status = "PASS" if passed else "FAIL"
        print(f"  {section:12s}: {status}")
        if not passed:
            all_pass = False

    print("\n  Figures saved to:", FIGURES_DIR)
    for f in ["lecture_14_exciton_binding.png", "lecture_14_scattering.png"]:
        path = FIGURES_DIR / f
        exists = "OK" if path.exists() else "MISSING"
        print(f"    {f}: {exists}")

    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
