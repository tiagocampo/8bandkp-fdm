#!/usr/bin/env python3
"""Lecture 13: Topological Superconductivity -- Chern numbers, Z2 invariants,
Majorana modes, Landau levels, spectral function, phase diagrams, and
Hall conductance.

Consolidates seven sections into one lecture-pair script:
  - Section 1: QWZ Chern numbers (verify_qwz_chern.py)
  - Section 2: BHZ Z2 invariant (verify_bhz_z2.py)
  - Section 3: BdG Majorana phase transition (sweep_rashba_bdg.py)
  - Section 4: Landau level quantization (verify_landau_levels.py)
  - Section 5: Spectral function A(k,E) for bulk GaAs
  - Section 6: Z2 phase diagram (BHZ analytic sweep)
  - Section 7: Hall conductance via Kubo formula (QWZ model)

Outputs:
  docs/lecture/figures/lecture_13_chern.png
  docs/lecture/figures/lecture_13_z2.png
  docs/lecture/figures/lecture_13_majorana.png
  docs/lecture/figures/lecture_13_landau.png
  docs/lecture/figures/lecture_13_spectral.png
  docs/lecture/figures/lecture_13_phase_diagram.png
  docs/lecture/figures/lecture_13_conductance.png
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
from star_helpers import (run_exe, parse_eigenvalues, parse_topology_result,
                          HBAR_J_S, E_CHARGE, M0_KG)

MATERIAL_DB = {
    'InAs': {'meff': 0.026, 'gamma1': 19.0, 'Eg': 0.417},
    'GaAs': {'meff': 0.067, 'gamma1': 6.95, 'Eg': 1.519},
}


# ---------------------------------------------------------------------------
# Section 1 -- QWZ Chern numbers
# ---------------------------------------------------------------------------
QWZ_CONFIGS = [
    ("topology_qwz_chern_u-0.8.cfg", -0.8, +1),
    ("topology_qwz_chern_u0.5.cfg",   0.5, -1),
    ("topology_qwz_chern_u2.5.cfg",   2.5,  0),
]


def run_qwz_chern(config_name, expected_C, work_dir):
    """Run topologicalAnalysis with a QWZ config, return parsed Chern number.

    Returns:
        (chern_int, topo_dict) or (None, {}) on failure.
    """
    config_path = CONFIGS_DIR / config_name
    rc, output_dir = run_exe(str(BUILD_DIR), "topologicalAnalysis",
                             str(config_path), work_dir)
    if rc != 0:
        print(f"    ERROR: topologicalAnalysis returned {rc}")
        return None, {}

    topo_file = os.path.join(output_dir, "topology_result.dat")
    if not os.path.isfile(topo_file):
        print(f"    ERROR: topology_result.dat not found")
        return None, {}

    topo = parse_topology_result(topo_file)

    # Parse Chern number from "# Chern number: N"
    content = Path(topo_file).read_text()
    match = re.search(r"# Chern number:\s*(-?\d+)", content)
    if match:
        return int(match.group(1)), topo
    return None, topo


def compute_berry_curvature_qwz(u, nk=100):
    """Compute Berry curvature F(kx, ky) for the QWZ model on an nk x nk grid.

    Uses the lower band of H = [[mz, sin_kx], [sin_ky, -mz]]
    with mz = u + cos(kx) + cos(ky).
    """
    dk = 2.0 * np.pi / nk
    kx = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    ky = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    KX, KY = np.meshgrid(kx, ky, indexing='ij')

    mz = u + np.cos(KX) + np.cos(KY)
    sin_kx = np.sin(KX)
    sin_ky = np.sin(KY)

    # Berry curvature for 2-level system:
    # F = (1/2) * mz / (mz^2 + sin_kx^2 + sin_ky^2)^(3/2)
    #   * (cos_kx * sin_ky - sin_kx * cos_ky) ... but simpler form:
    # F = (1/4) * sin(mz_hat . (dm/dkx x dm/dky)) where m = (sin_ky, sin_kx, mz)
    # For 2-band: F = (1/2) * d_hat . (dd/dkx x dd/dky) with d = (sin_ky, sin_kx, mz)
    # ddkx = (0, cos_kx, -sin_kx), ddky = (cos_ky, 0, -sin_ky)
    # cross = (-cos_kx * (-sin_ky) - 0, -sin_kx * (-sin_ky) - (-sin_kx)*cos_ky, ...)
    # Actually: cross(ddkx, ddky)_z = cos_kx*cos_ky - 0 = cos_kx*cos_ky ... no
    # Let d = (sin_ky, sin_kx, mz)
    # dd/dkx = (0, cos_kx, -sin_kx)
    # dd/dky = (cos_ky, 0, -sin_ky)
    # dd/dkx x dd/dky = |i  j  k|
    #                    |0  cos_kx  -sin_kx|
    #                    |cos_ky  0  -sin_ky|
    # i: cos_kx*(-sin_ky) - (-sin_kx)*0 = -cos_kx*sin_ky
    # j: -sin_kx*(-sin_ky) - 0*(-sin_kx) = sin_kx*sin_ky  ... wait no
    # j: (0*(-sin_ky) - (-sin_kx)*cos_ky) = sin_kx*cos_ky
    # k: 0*0 - cos_kx*cos_ky = -cos_kx*cos_ky
    # So cross = (-cos_kx*sin_ky, sin_kx*cos_ky, -cos_kx*cos_ky) ... hmm
    # d . cross = sin_ky*(-cos_kx*sin_ky) + sin_kx*(sin_kx*cos_ky)
    #           + mz*(-cos_kx*cos_ky)
    # Actually let me use the standard formula:
    # F = (1/2) * d_hat . (dd/dkx x dd/dky) / |d|^3 ... no
    # For a 2-band Hamiltonian H = d . sigma, the Berry curvature is:
    # F = (1/2) * d_hat . (dd/dkx x dd/dky) where d_hat = d/|d|

    d_norm = np.sqrt(mz**2 + sin_kx**2 + sin_ky**2)
    d_norm = np.where(d_norm < 1e-30, 1e-30, d_norm)

    dd_dkx = np.stack([np.zeros_like(KX), np.cos(KX), -np.sin(KX)], axis=-1)
    dd_dky = np.stack([np.cos(KY), np.zeros_like(KY), -np.sin(KY)], axis=-1)

    cross = np.cross(dd_dkx, dd_dky)
    d_vec = np.stack([sin_ky, sin_kx, mz], axis=-1)

    F = 0.5 * np.sum(d_vec * cross, axis=-1) / (d_norm**3)
    return KX, KY, F


def section1_chern():
    """Section 1: QWZ Chern number verification and Berry curvature plot."""
    print("\n" + "=" * 60)
    print("  Section 1: QWZ Chern Numbers")
    print("  Qi, Wu & Zhang, Phys. Rev. B 74, 085308 (2006)")
    print("=" * 60)

    exe_path = BUILD_DIR / "src" / "topologicalAnalysis"
    if not exe_path.exists():
        print("  ERROR: topologicalAnalysis not found. Build first.")
        return False

    results = {}
    all_pass = True
    for config_name, u, expected_C in QWZ_CONFIGS:
        print(f"\n  Running {config_name} (u={u:+.1f}, expected C={expected_C})...")
        with tempfile.TemporaryDirectory() as tmpdir:
            C, topo = run_qwz_chern(config_name, expected_C, tmpdir)

        if C is None:
            print(f"    FAIL: Could not parse Chern number")
            all_pass = False
            results[u] = None
        elif C == expected_C:
            print(f"    PASS: C = {C}")
            results[u] = C
        else:
            print(f"    FAIL: C = {C}, expected {expected_C}")
            all_pass = False
            results[u] = C

    # --- Generate Berry curvature heatmap for u=-0.8 (C=+1 topological) ---
    print("\n  Computing Berry curvature heatmap (u=-0.8)...")
    KX, KY, F = compute_berry_curvature_qwz(u=-0.8, nk=200)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Berry curvature heatmap
    ax = axes[0]
    vmax = np.percentile(np.abs(F), 95)
    im = ax.pcolormesh(KX, KY, F, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                       shading='auto')
    ax.set_xlabel(r'$k_x$', fontsize=12)
    ax.set_ylabel(r'$k_y$', fontsize=12)
    ax.set_title(r'Berry curvature $F(k)$, QWZ $u=-0.8$ ($C=+1$)', fontsize=12)
    ax.set_aspect('equal')
    fig.colorbar(im, ax=ax, label=r'$F(k)$')

    # Phase diagram bar chart
    ax2 = axes[1]
    u_vals = [-0.8, 0.5, 2.5]
    u_strs = [f"{u:+.1f}" for u in u_vals]
    C_vals = [results.get(u, 0) for u in u_vals]
    expected = [+1, -1, 0]
    bar_colors = ['crimson' if c == 1 else 'steelblue' if c == -1 else 'forestgreen'
                  for c in expected]
    bars = ax2.bar(u_strs, C_vals, color=bar_colors, width=0.5, edgecolor='k', lw=0.5)
    ax2.axhline(0, color='gray', ls='--', lw=0.6)
    for bar, c in zip(bars, C_vals):
        if c is not None:
            yoff = 0.05 if c >= 0 else -0.15
            ax2.text(bar.get_x() + bar.get_width() / 2, c + yoff, str(c),
                     ha='center', va='bottom' if c >= 0 else 'top', fontsize=13)
    ax2.set_xlabel(r'Mass parameter $u$', fontsize=12)
    ax2.set_ylabel(r'Chern number $C$', fontsize=12)
    ax2.set_title('QWZ Phase Diagram (nk=50, Fortran FHS)', fontsize=12)
    ax2.set_yticks([-2, -1, 0, 1, 2])
    ax2.grid(True, alpha=0.3, axis='y')

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "lecture_13_chern.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {FIGURES_DIR / 'lecture_13_chern.png'}")

    status = "PASS" if all_pass else "FAIL"
    print(f"\n  Section 1 result: {status}")
    return all_pass


# ---------------------------------------------------------------------------
# Section 2 -- BHZ Z2 invariant
# ---------------------------------------------------------------------------
BHZ_WIDTH_MIN = 40
BHZ_WIDTH_MAX = 100
BHZ_WIDTH_STEP = 2


def run_bhz_z2(exe, width_angstrom, work_dir):
    """Run topologicalAnalysis for a given wire width, return Z2 invariant.

    Modifies the trivial config template with new wire_width.
    """
    config_path = CONFIGS_DIR / "topology_bhz_z2_trivial.cfg"
    config = config_path.read_text()
    config = config.replace("wire_width: 58.0", f"wire_width: {width_angstrom:.1f}")
    config = config.replace("wire_height: 58.0", f"wire_height: {width_angstrom:.1f}")

    input_cfg = os.path.join(work_dir, "input.cfg")
    with open(input_cfg, "w") as f:
        f.write(config)

    output_dir = os.path.join(work_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    try:
        result = subprocess.run(
            [str(exe)],
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=120,
        )
    except subprocess.TimeoutExpired:
        print(f"    TIMEOUT after 120s for width={width_angstrom}A")
        return None
    if result.returncode != 0:
        print(f"    ERROR: topologicalAnalysis returned {result.returncode} "
              f"for width={width_angstrom}A")
        if result.stderr:
            print(f"    stderr: {result.stderr[:500]}")
        if result.stdout:
            print(f"    stdout: {result.stdout[:500]}")
        return None

    topo_file = os.path.join(output_dir, "topology_result.dat")
    if not os.path.isfile(topo_file):
        print(f"    ERROR: topology_result.dat not found for width={width_angstrom}A")
        if result.stdout:
            print(f"    stdout: {result.stdout[:500]}")
        return None

    content = Path(topo_file).read_text()
    match = re.search(r'Z2 invariant:\s*(\d+)', content)
    if match:
        return int(match.group(1))
    print(f"    WARNING: Z2 invariant not found in topology_result.dat for width={width_angstrom}A")
    return None


def section2_z2():
    """Section 2: BHZ Z2 invariant sweep across wire widths."""
    print("\n" + "=" * 60)
    print("  Section 2: BHZ Z2 Invariant vs Wire Width")
    print("=" * 60)

    exe = BUILD_DIR / "src" / "topologicalAnalysis"
    if not exe.exists():
        print("  ERROR: topologicalAnalysis not found.")
        return False

    widths = list(range(BHZ_WIDTH_MIN, BHZ_WIDTH_MAX + 1, BHZ_WIDTH_STEP))
    results = []

    for w in widths:
        with tempfile.TemporaryDirectory() as tmpdir:
            z2 = run_bhz_z2(exe, w, tmpdir)
        results.append((w, z2))
        status = f"Z2={z2}" if z2 is not None else "FAILED"
        print(f"  width={w:3d}A: {status}")

    # Validate: trivial below 60A, topological at/above 80A.
    # Avoid testing exactly at the critical width (~70A) where numerical
    # noise can flip the Z2 invariant.
    all_pass = True
    for w, z2 in results:
        if z2 is None:
            all_pass = False
            continue
        if w <= 60:
            expected_z2 = 0
        elif w >= 80:
            expected_z2 = 1
        else:
            continue  # skip transition region
        if z2 != expected_z2:
            all_pass = False

    # --- Generate phase diagram ---
    print("\n  Generating Z2 phase diagram...")
    valid = [(w, z2) for w, z2 in results if z2 is not None]
    if valid:
        w_v, z2_v = zip(*valid)

        fig, ax = plt.subplots(figsize=(10, 5))

        ax.step(w_v, z2_v, where='post', linewidth=2, marker='o',
                markersize=6, color='navy')

        ax.fill_between(w_v, 0, z2_v, step='post', alpha=0.15, color='steelblue')
        ax.axvline(x=70, color='red', ls='--', lw=1.5, label='Critical width (70 A)')

        ax.set_xlabel('Wire Width (Angstrom)', fontsize=12)
        ax.set_ylabel(r'$\mathbb{Z}_2$ Invariant', fontsize=12)
        ax.set_title(r'BHZ Model: $\mathbb{Z}_2$ Invariant vs Wire Width', fontsize=13)
        ax.set_ylim(-0.1, 1.5)
        ax.set_xlim(BHZ_WIDTH_MIN - 2, BHZ_WIDTH_MAX + 2)
        ax.set_yticks([0, 1])
        ax.set_yticklabels([r'0 (Trivial)', r'1 (Topological)'])
        ax.legend(loc='center right', fontsize=10)
        ax.grid(True, alpha=0.3)

        ax.annotate('Trivial phase\n(M = +10 meV)', xy=(55, 0.15),
                     fontsize=10, ha='center', color='navy')
        ax.annotate('Topological phase\n(M = -10 meV)', xy=(85, 0.85),
                     fontsize=10, ha='center', color='darkorange')

        fig.tight_layout()
        fig.savefig(FIGURES_DIR / "lecture_13_z2.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {FIGURES_DIR / 'lecture_13_z2.png'}")

    status = "PASS" if all_pass else "FAIL"
    print(f"\n  Section 2 result: {status}")
    return all_pass


# ---------------------------------------------------------------------------
# Section 3 -- BdG Majorana phase transition (QW mode)
# ---------------------------------------------------------------------------

# InAs/GaAs QW: 200 A InAs well in 500 A GaAs barriers, N=101 FD points.
# Chemical potential mu = -0.1413 eV tuned to the CB minimum of the QW.
# g_factor = 15 (InAs), delta_0 = 0.2 meV (Al proximity s-wave).
# BdG matrix is 16N x 16N = 1616 x 1616 (dense LAPACK, ~0.3 s per run).
# Intra-band s-wave pairing: Delta = delta_0 * (i sigma_y x I_4).
# Expected B_crit ~ 0.25 T from Gamma_c = delta_0 / (g * mu_B).
_BDG_QW_CONFIG = """\
waveVector: k0
waveVectorMax: 0.0
waveVectorStep: 1
confinement: 1
FDstep: 101
FDorder: 2
numLayers: 2
material1: GaAs -250 250 0
material2: InAs -100 100 0
numcb: 2
numvb: 6
ExternalField: 0  EF
EFParams: 0.0
SC: 0
topology: T
mode: bdg
compute_chern: F
compute_hall: F
qwz_u: 0.0
compute_z2: F
extract_edge_states: F
edge_E_window: 0.01
compute_ldos: F
ldos_eta: 0.001
ldos_E_range: -0.1
ldos_E_range: 0.1
ldos_num_E: 20
bdg: T
mu: {mu}
delta_0: {delta_0}
g_factor: {g_factor}
B_vec: 0.0 0.0 {B}
gauge: landau
kz: {kz}
"""


def _run_bdg_qw(kz_val, B_val, work_dir, mu=-0.1413, delta_0=0.0002,
                g_factor=15.0):
    """Run topologicalAnalysis in BdG QW mode.

    Parameters:
        kz_val: in-plane momentum (1/A)
        B_val: magnetic field along z (T)
        work_dir: temporary working directory
        mu: chemical potential (eV)
        delta_0: s-wave pairing amplitude (eV)
        g_factor: effective g-factor

    Returns:
        dict with keys: eigenvalues (eV), min_gap (eV), profile (z, rho) or None
    """
    exe = BUILD_DIR / "src" / "topologicalAnalysis"
    config = _BDG_QW_CONFIG.format(mu=mu, delta_0=delta_0, g_factor=g_factor,
                                   B=B_val, kz=kz_val)

    input_cfg = os.path.join(work_dir, "input.cfg")
    with open(input_cfg, "w") as f:
        f.write(config)

    output_dir = os.path.join(work_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    try:
        result = subprocess.run(
            [str(exe)],
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=120,
        )
    except subprocess.TimeoutExpired:
        return None

    if result.returncode != 0:
        return None

    out = {"eigenvalues": None, "min_gap": None, "profile": None}

    # Parse eigenvalues
    eig_path = os.path.join(output_dir, "bdg_eigenvalues.dat")
    if os.path.isfile(eig_path):
        evals = []
        with open(eig_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    evals.append(float(parts[1]))
        out["eigenvalues"] = np.array(evals)

    # Parse min gap from topology_result.dat
    topo_path = os.path.join(output_dir, "topology_result.dat")
    if os.path.isfile(topo_path):
        content = Path(topo_path).read_text()
        match = re.search(r"# Min gap \(eV\):\s*([0-9.e+-]+)", content)
        if match:
            out["min_gap"] = float(match.group(1))

    # Parse lowest-state spatial profile
    prof_path = os.path.join(output_dir, "bdg_lowest_state_profile.dat")
    if os.path.isfile(prof_path):
        z_vals, rho_vals = [], []
        with open(prof_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    z_vals.append(float(parts[0]))
                    rho_vals.append(float(parts[1]))
        if z_vals:
            out["profile"] = (np.array(z_vals), np.array(rho_vals))

    return out


def section3_majorana():
    """Section 3: BdG Majorana phase transition (QW, 8-band k.p).

    Oreg-Lutchyn model via 8-band k.p: InAs/GaAs QW with proximity-induced
    s-wave pairing (Delta = 0.2 meV, g = 15). Chemical potential tuned to
    CB edge. The topological phase transition occurs at B_crit ~ 0.25 T
    where the Zeeman energy overcomes the SC gap.

    Generates 2x2 figure:
      (a) Trivial dispersion at B=0: SC gap Delta_0 at Fermi level
      (b) Topological dispersion at B > B_crit: reopened topological gap
      (c) Gap vs B: V-shaped closing/reopening at B_crit ~ 0.25 T
      (d) Near-zero BdG levels: gap closing at topological transition
    """
    print("\n" + "=" * 60)
    print("  Section 3: BdG Majorana Phase Transition (QW mode)")
    print("  InAs/GaAs QW, 8-band k.p + s-wave pairing")
    print("  Oreg-Lutchyn model: Delta=0.2meV, g=15, B_crit~0.25T")
    print("=" * 60)

    exe = BUILD_DIR / "src" / "topologicalAnalysis"
    if not exe.exists():
        print("  ERROR: topologicalAnalysis not found.")
        return False

    mu_eV = -0.1413      # eV — CB minimum of 200A InAs QW
    delta_0_eV = 0.0002  # eV = 0.2 meV (Al proximity s-wave pairing)
    g_factor = 15.0      # InAs effective g-factor
    B_crit_approx = 0.25 # T — approximate critical field
    B_topo = 0.6         # T — well into topological phase
    Delta_meV = delta_0_eV * 1000

    # Analytical B_crit for single-band model: B_c = Delta_0 / (g * mu_B)
    mu_B_eV_per_T = 5.7884e-5  # eV/T
    B_crit_theory = delta_0_eV / (g_factor * mu_B_eV_per_T)
    print(f"  Theory: B_crit = Delta/(g*mu_B) = {B_crit_theory:.3f} T")

    # --- (a,b) Band dispersion: E vs k_par at B=0 and B>B_crit ---
    k_par_vals = np.linspace(-0.04, 0.04, 15)
    print(f"\n  Dispersion sweep: {len(k_par_vals)} k_par points at B=0 and B={B_topo} T ...")

    dispersion_trivial = []
    dispersion_topo = []

    for kz in k_par_vals:
        with tempfile.TemporaryDirectory() as tmpdir:
            res = _run_bdg_qw(kz, 0.0, tmpdir, mu=mu_eV, delta_0=delta_0_eV,
                              g_factor=g_factor)
        dispersion_trivial.append((kz, res["eigenvalues"] if res else None))

        with tempfile.TemporaryDirectory() as tmpdir:
            res = _run_bdg_qw(kz, B_topo, tmpdir, mu=mu_eV, delta_0=delta_0_eV,
                              g_factor=g_factor)
        dispersion_topo.append((kz, res["eigenvalues"] if res else None))

    n_ok_triv = sum(1 for _, e in dispersion_trivial if e is not None)
    n_ok_topo = sum(1 for _, e in dispersion_topo if e is not None)
    print(f"  Trivial dispersion: {n_ok_triv}/{len(k_par_vals)} points OK")
    print(f"  Topological dispersion: {n_ok_topo}/{len(k_par_vals)} points OK")

    # --- (c,d) Gap + near-zero levels vs B sweep ---
    B_sweep = np.concatenate([
        np.linspace(0.01, 0.15, 6),
        np.linspace(0.15, 0.40, 15),
        np.linspace(0.40, 1.0, 7),
    ])
    print(f"\n  Gap sweep: {len(B_sweep)} B values from {B_sweep[0]:.2f} to {B_sweep[-1]:.2f} T ...")

    gap_data = []
    near_zero_data = []

    for B in B_sweep:
        with tempfile.TemporaryDirectory() as tmpdir:
            res = _run_bdg_qw(0.0, B, tmpdir, mu=mu_eV, delta_0=delta_0_eV,
                              g_factor=g_factor)
        if res and res["min_gap"] is not None:
            gap_data.append((B, res["min_gap"]))
        else:
            gap_data.append((B, None))

        if res and res["eigenvalues"] is not None:
            evals = res["eigenvalues"]
            idx_sorted = np.argsort(np.abs(evals))
            n_extract = min(8, len(evals))
            near_zero_data.append((B, evals[idx_sorted[:n_extract]]))
        else:
            near_zero_data.append((B, None))

        g_meV = res["min_gap"] * 1000 if res and res["min_gap"] is not None else None
        status = f"{g_meV:.3f} meV" if g_meV is not None else "FAILED"
        print(f"  B={B:5.3f} T: min_gap = {status}")

    # Find B_crit from gap minimum
    valid_gaps = [(B, g) for B, g in gap_data if g is not None and g > 0]
    B_crit_est = None
    if valid_gaps:
        B_arr_g, g_arr = zip(*valid_gaps)
        idx_min = np.argmin(g_arr)
        B_crit_est = B_arr_g[idx_min]
        print(f"\n  Measured B_crit ~ {B_crit_est:.3f} T "
              f"(gap minimum = {g_arr[idx_min]*1000:.3f} meV)")

    # --- Validation ---
    all_pass = True

    # Trivial gap at B ~ 0 should be approximately Delta_0
    gap_b0 = next((g for B, g in gap_data if B < 0.05 and g is not None), None)
    if gap_b0 is not None and 0.1e-3 < gap_b0 < 1.0e-3:
        print(f"\n  Trivial gap at B~0: {gap_b0*1000:.3f} meV ≈ Delta_0 — PASS")
    else:
        v = 'None' if gap_b0 is None else f'{gap_b0*1000:.3f} meV'
        print(f"\n  Trivial gap at B~0: {v} — expected ~0.2 meV — FAIL")
        all_pass = False

    # Topological gap at B >> B_crit should be > 0 (reopened)
    gap_topo = next((g for B, g in gap_data if B > B_topo - 0.05 and g is not None), None)
    if gap_topo is not None and gap_topo > 0:
        print(f"  Topological gap at B={B_topo} T: {gap_topo*1000:.3f} meV — PASS")
    else:
        v = 'None' if gap_topo is None else f'{gap_topo*1000:.3f} meV'
        print(f"  Topological gap at B={B_topo} T: {v} — FAIL")
        all_pass = False

    # B_crit should be close to theoretical prediction (within 50%)
    if B_crit_est is not None:
        rel_err = abs(B_crit_est - B_crit_theory) / B_crit_theory
        if rel_err < 0.5:
            print(f"  B_crit = {B_crit_est:.3f} T vs theory {B_crit_theory:.3f} T "
                  f"(err={rel_err:.0%}) — PASS")
        else:
            print(f"  B_crit = {B_crit_est:.3f} T vs theory {B_crit_theory:.3f} T "
                  f"(err={rel_err:.0%}) — FAIL")
            all_pass = False

    # Particle-hole symmetry check
    for label, disp in [("trivial", dispersion_trivial), ("topological", dispersion_topo)]:
        for kz, evals in disp:
            if evals is not None and len(evals) > 0:
                evals_s = np.sort(evals)
                n2 = len(evals_s) // 2
                asym = np.max(np.abs(evals_s[:n2] + evals_s[-1:-n2-1:-1]))
                if asym > 1e-6:
                    print(f"  WARNING: {label} at kz={kz:.3f} PH asymmetry = {asym:.2e} eV")
                break

    # --- Generate 2x2 figure ---
    print("\n  Generating Majorana phase transition figure ...")
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    E_ZOOM = 2.0  # meV — energy window for dispersion (near gap)

    def _plot_dispersion(ax, disp_data, title, color):
        """Plot BdG eigenvalues vs k_par, zoomed to near-gap region."""
        for kz, evals in disp_data:
            if evals is None:
                continue
            evals_meV = evals * 1000
            mask = np.abs(evals_meV) < E_ZOOM
            ax.plot([kz] * mask.sum(), evals_meV[mask], '.', color=color,
                    ms=2, alpha=0.7)
        ax.axhline(0, color='k', ls='-', lw=0.5, alpha=0.3)
        ax.set_xlabel(r'$k_\parallel$ (1/A)', fontsize=11)
        ax.set_ylabel('Energy (meV)', fontsize=11)
        ax.set_title(title, fontsize=12)
        ax.set_ylim(-E_ZOOM, E_ZOOM)
        ax.grid(True, alpha=0.2)

    # (a) Trivial dispersion — SC gap Delta_0 at k=0
    _plot_dispersion(axes[0, 0], dispersion_trivial,
                     r'Trivial phase ($B = 0$, gap $\approx \Delta_0$)',
                     'steelblue')

    # (b) Topological dispersion — reopened gap
    _plot_dispersion(axes[0, 1], dispersion_topo,
                     rf'Topological ($B = {B_topo}$ T, gap reopened)',
                     'crimson')

    # (c) Gap vs B — V-shaped topological transition
    ax_gap = axes[1, 0]
    B_valid = np.array([B for B, g in gap_data if g is not None])
    g_valid = np.array([g * 1000 for B, g in gap_data if g is not None])
    if len(B_valid) > 0:
        ax_gap.plot(B_valid, g_valid, 'o-', ms=4, color='navy', lw=1.5)
    if B_crit_est is not None:
        ax_gap.axvline(B_crit_est, color='red', ls='--', lw=1.5,
                       label=rf'$B_{{crit}} \approx {B_crit_est:.2f}$ T')
    ax_gap.axhline(Delta_meV, color='gray', ls=':', lw=1,
                   label=rf'$\Delta_0 = {Delta_meV:.1f}$ meV')
    if len(B_valid) > 0 and B_crit_est is not None:
        ax_gap.axvspan(B_crit_est, max(B_valid), alpha=0.08, color='gold',
                       label='Topological')
        ax_gap.axvspan(min(B_valid), B_crit_est, alpha=0.08, color='steelblue',
                       label='Trivial')
    ax_gap.set_xlabel('Magnetic field $B$ (T)', fontsize=11)
    ax_gap.set_ylabel('Min spectral gap (meV)', fontsize=11)
    ax_gap.set_title('Topological Phase Transition (gap closes at $B_{crit}$)',
                     fontsize=12)
    ax_gap.annotate(r'InAs/GaAs QW: $g=15$, $\Delta_0=0.2$ meV'
                    '\n' r'$B_{crit} \approx \Delta_0 / (g\mu_B) \approx 0.23$ T',
                    xy=(0.97, 0.97), xycoords='axes fraction',
                    fontsize=8, ha='right', va='top',
                    bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow',
                              ec='gray', alpha=0.9))
    ax_gap.legend(fontsize=9, loc='center right')
    ax_gap.grid(True, alpha=0.3)

    # (d) Near-zero BdG levels vs B
    ax_fan = axes[1, 1]
    for B, ev in near_zero_data:
        if ev is not None:
            ev_meV = ev * 1000
            mask = np.abs(ev_meV) < E_ZOOM
            ax_fan.scatter([B] * mask.sum(), ev_meV[mask], s=20, c='steelblue',
                           alpha=0.7, edgecolors='none', zorder=3)
    ax_fan.axhline(0, color='k', ls='-', lw=0.5, alpha=0.3)
    if B_crit_est is not None:
        ax_fan.axvline(B_crit_est, color='red', ls='--', lw=1, alpha=0.7,
                       label=rf'$B_{{crit}}$')
        ax_fan.annotate('Gap\ncloses', xy=(B_crit_est, 0),
                        xytext=(B_crit_est + 0.15, 0.5), fontsize=9,
                        arrowprops=dict(arrowstyle='->', color='red'),
                        color='red')
    ax_fan.set_xlabel('Magnetic field $B$ (T)', fontsize=11)
    ax_fan.set_ylabel('BdG eigenvalues (meV)', fontsize=11)
    ax_fan.set_title('Near-Zero BdG Levels', fontsize=12)
    ax_fan.set_ylim(-E_ZOOM, E_ZOOM)
    ax_fan.legend(fontsize=9)
    ax_fan.grid(True, alpha=0.3)

    fig.suptitle('Majorana Phase Transition in InAs/GaAs QW\n'
                 r'(8-band k.p, $\Delta = 0.2$ meV, $g = 15$, intra-band s-wave)',
                 fontsize=13, fontweight='bold', y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "lecture_13_majorana.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {FIGURES_DIR / 'lecture_13_majorana.png'}")

    status = "PASS" if all_pass else "FAIL"
    print(f"\n  Section 3 result: {status}")
    return all_pass


# ---------------------------------------------------------------------------
# Section 4 -- Landau levels
# ---------------------------------------------------------------------------


def compute_cyclotron_energy(m_eff, B):
    """Compute hbar*omega_c in meV. m_eff in units of m0, B in Tesla."""
    omega_c = E_CHARGE * B / (m_eff * M0_KG)
    hbar_omega_J = HBAR_J_S * omega_c
    return hbar_omega_J * 6.241509e21  # J -> meV


def run_landau(config_name, work_dir):
    """Run bandStructure with a Landau config, return eigenvalues in meV.

    BUGFIX: reads eigenvalues.dat (not band_results.dat).
    """
    config_path = CONFIGS_DIR / config_name
    rc, output_dir = run_exe(str(BUILD_DIR), "bandStructure",
                             str(config_path), work_dir)
    if rc != 0:
        print(f"    ERROR: bandStructure returned {rc}")
        return None

    # BUGFIX: the correct output file is eigenvalues.dat
    eig_path = os.path.join(output_dir, "eigenvalues.dat")
    if not os.path.isfile(eig_path):
        print(f"    WARNING: eigenvalues.dat not found")
        return None

    data = parse_eigenvalues(eig_path)
    if not data:
        return None

    # Return eigenvalues at first k-point, converted to meV
    _, evals = data[0]
    return [e * 1000.0 for e in evals]  # eV -> meV


def run_landau_bsweep(work_dir):
    """Run bandStructure with landau_bulk_InAs_Bsweep.cfg for Zeeman fan.

    Returns:
        list of (B, eigenvalues_meV) tuples, or None on failure.
    """
    config_path = CONFIGS_DIR / "landau_bulk_InAs_Bsweep.cfg"
    rc, output_dir = run_exe(str(BUILD_DIR), "bandStructure",
                             str(config_path), work_dir, timeout=300)
    if rc != 0:
        print(f"    ERROR: bandStructure returned {rc} for B-sweep")
        return None

    eig_path = os.path.join(output_dir, "eigenvalues.dat")
    if not os.path.isfile(eig_path):
        print(f"    WARNING: eigenvalues.dat not found for B-sweep")
        return None

    data = parse_eigenvalues(eig_path)
    # Each line: B_value eval1 eval2 ...
    return [(d[0], [e * 1000.0 for e in d[1]]) for d in data]


def section4_landau():
    """Section 4: Landau level quantization and Zeeman fan diagram."""
    print("\n" + "=" * 60)
    print("  Section 4: Landau Level Quantization")
    print("=" * 60)

    material = 'InAs'
    B = 5.0
    meff, cb_edge = MATERIAL_DB[material]['meff'], MATERIAL_DB[material]['Eg']
    hbar_omega = compute_cyclotron_energy(meff, B)
    print(f"  {material}: m* = {meff:.3f} m0, CB edge = {cb_edge*1000:.0f} meV")
    print(f"  hbar*omega_c = {hbar_omega:.2f} meV")

    # Run single-B Landau config
    eigenvalues = None
    with tempfile.TemporaryDirectory() as tmpdir:
        eigenvalues = run_landau("landau_bulk_InAs.cfg", tmpdir)

    if eigenvalues is not None:
        print(f"  Parsed {len(eigenvalues)} eigenvalues (meV):")
        for i, ev in enumerate(eigenvalues[:8]):
            print(f"    {i:2d}: {ev:10.3f} meV")
    else:
        print("  WARNING: No eigenvalues parsed from landau_bulk_InAs.cfg")

    # Analytical Landau levels
    N_MAX = 6
    n_analytical = list(range(N_MAX + 1))
    e_analytical = [cb_edge * 1000.0 + hbar_omega * (n + 0.5) for n in n_analytical]

    # Verify spacing
    all_pass = True
    if eigenvalues is not None and len(eigenvalues) >= 2:
        # Filter conduction band eigenvalues (above CB edge)
        cb_meV = cb_edge * 1000.0
        cond_evals = sorted([e for e in eigenvalues if e > cb_meV * 0.5])
        if len(cond_evals) >= 2:
            spacing = cond_evals[1] - cond_evals[0]
            expected_spacing = hbar_omega
            rel_err = abs(spacing - expected_spacing) / expected_spacing
            print(f"\n  CB Landau spacing: {spacing:.2f} meV "
                  f"(expected hbar*omega_c = {expected_spacing:.2f} meV, "
                  f"rel_err = {rel_err:.2%})")
            # Allow 20% tolerance for 8-band non-parabolicity
            if rel_err > 0.20:
                all_pass = False

    # Run B-sweep for Zeeman fan
    bsweep_data = None
    print("\n  Running Landau B-sweep...")
    with tempfile.TemporaryDirectory() as tmpdir:
        bsweep_data = run_landau_bsweep(tmpdir)

    if bsweep_data is not None:
        print(f"  Parsed {len(bsweep_data)} B-points")
    else:
        print("  WARNING: B-sweep did not produce data")

    # --- Generate combined figure ---
    print("\n  Generating Landau level + Zeeman fan figure...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left panel: Landau level comparison
    ax = axes[0]
    ax.hlines(e_analytical, -0.5, 2.5, colors='forestgreen', linestyles='--',
              lw=2.5,
              label=rf'Analytical $E_n = CB + \hbar\omega_c\,(n+\frac{{1}}{{2}})$',
              zorder=3)

    if eigenvalues is not None:
        cb_meV = cb_edge * 1000.0
        cond_evals = sorted([e for e in eigenvalues if e > cb_meV * 0.5])
        n_comp = list(range(len(cond_evals)))
        ax.scatter(n_comp, cond_evals[:N_MAX + 1], color='royalblue', s=120,
                   zorder=5, label='Computed (8-band k.p)',
                   edgecolors='k', linewidths=0.7, marker='o')
    else:
        ax.text(0.5, 0.5, 'No computed data', transform=ax.transAxes,
                ha='center', fontsize=12, color='gray')

    for n, E in zip(n_analytical, e_analytical):
        ax.annotate(f'n={n}\n{E:.1f} meV', xy=(2.15, E), fontsize=7,
                     color='forestgreen', va='center')

    info_text = (
        f'{material} at B={B} T:\n'
        f'  m* = {meff:.3f} m0\n'
        f'  hbar*omega_c = {hbar_omega:.2f} meV'
    )
    ax.annotate(info_text, xy=(0.55, 0.55), xycoords='axes fraction',
                 fontsize=9, color='darkblue',
                 bbox=dict(boxstyle='round,pad=0.5', facecolor='lavender',
                           edgecolor='steelblue', lw=1.5),
                 va='top', ha='left')

    ax.set_xlabel(r'Landau level index $n$', fontsize=12)
    ax.set_ylabel('Energy (meV)', fontsize=12)
    ax.set_title(f'Landau Levels: 8-band k.p vs Effective Mass\n'
                 f'{material}, B = {B} T', fontsize=12)
    ax.set_xlim(-0.5, 3.2)
    ax.set_xticks(range(N_MAX + 1))
    ax.grid(True, alpha=0.3, axis='y')
    ax.legend(fontsize=9, loc='upper left')

    # Right panel: Zeeman fan diagram
    ax2 = axes[1]
    if bsweep_data and len(bsweep_data) > 1:
        B_arr = np.array([d[0] for d in bsweep_data])
        eval_arr = np.array([d[1] for d in bsweep_data])

        # Plot each band as a function of B
        n_bands = eval_arr.shape[1]
        band_colors = plt.cm.viridis(np.linspace(0.1, 0.9, n_bands))
        for i in range(n_bands):
            ax2.plot(B_arr, eval_arr[:, i], '-', color=band_colors[i], lw=1.2)

        # Overlay analytical E_n = CB + (n+1/2)*hbar*omega_c(B)
        B_anal = np.linspace(B_arr.min(), B_arr.max(), 100)
        for n in range(min(N_MAX + 1, n_bands)):
            hbar_omega_B = np.array([compute_cyclotron_energy(meff, b) for b in B_anal])
            E_n = cb_edge * 1000.0 + hbar_omega_B * (n + 0.5)
            ax2.plot(B_anal, E_n, '--', color='red', lw=0.8, alpha=0.5)

        ax2.set_xlabel('Magnetic field $B$ (T)', fontsize=12)
        ax2.set_ylabel('Energy (meV)', fontsize=12)
        ax2.set_title(f'Zeeman Fan Diagram: {material}\n'
                       f'(dashed red = analytical Landau levels)', fontsize=12)
        ax2.grid(True, alpha=0.3)
    else:
        ax2.text(0.5, 0.5, 'B-sweep data not available',
                 transform=ax2.transAxes, ha='center', fontsize=12, color='gray')
        ax2.set_title('Zeeman Fan Diagram (no data)', fontsize=12)

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "lecture_13_landau.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {FIGURES_DIR / 'lecture_13_landau.png'}")

    status = "PASS" if all_pass else "FAIL"
    print(f"\n  Section 4 result: {status}")
    return all_pass


# ---------------------------------------------------------------------------
# Section 5 -- Spectral function A(k, E) for bulk GaAs
# ---------------------------------------------------------------------------

_SPECTRAL_BULK_CONFIG = """\
waveVector: k0
waveVectorMax: 0.0
waveVectorStep: 1
confinement: 0
FDstep: 1
FDorder: 2
numLayers: 1
material1: GaAs
numcb: 2
numvb: 6
ExternalField: 0  EF
EFParams: 0.0
SC: 0
topology: T
mode: spectral
compute_chern: F
compute_hall: F
qwz_u: 0.0
compute_z2: F
extract_edge_states: F
edge_E_window: 0.01
compute_ldos: F
ldos_eta: 0.001
ldos_E_range: -0.1
ldos_E_range: 0.1
ldos_num_E: 20
compute_spectral: T
spectral_k_grid: -0.05 0.05 51
spectral_E_grid: -2.0 2.0 101 0.002
"""


def _parse_spectral_function(filepath):
    """Parse spectral_function.dat into (k_array, E_array, A_2d).

    Returns:
        (k_arr, E_arr, A_2d) or (None, None, None) on failure.
    """
    nk = nE = None
    k_list, E_list, A_list = [], [], []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                # Parse header for nk and nE
                m_nk = re.search(r'nk=(\d+)', line)
                if m_nk:
                    nk = int(m_nk.group(1))
                m_nE = re.search(r'nE=(\d+)', line)
                if m_nE:
                    nE = int(m_nE.group(1))
                continue
            parts = line.split()
            if len(parts) >= 3:
                k_list.append(float(parts[0]))
                E_list.append(float(parts[1]))
                A_list.append(float(parts[2]))

    if nk is None or nE is None or len(k_list) == 0:
        return None, None, None

    k_arr = np.array(k_list)
    E_arr = np.array(E_list)
    A_flat = np.array(A_list)

    # Reshape: Fortran writes k outer loop, E inner loop
    A_2d = A_flat.reshape(nk, nE)

    # Extract unique k and E arrays
    k_unique = k_arr.reshape(nk, nE)[:, 0]
    E_unique = E_arr.reshape(nk, nE)[0, :]

    return k_unique, E_unique, A_2d


def section5_spectral():
    """Section 5: Spectral function A(k,E) for bulk GaAs.

    Computes A(k,E) = -(1/pi) Im Tr G(k,E) on a 51x101 grid using the
    8-band bulk Hamiltonian. The spectral function should show peaks along
    the band dispersion, with valence bands near E~0 and the conduction
    band near E~1.519 eV (GaAs band gap).
    """
    print("\n" + "=" * 60)
    print("  Section 5: Spectral Function A(k,E) for Bulk GaAs")
    print("=" * 60)

    exe = BUILD_DIR / "src" / "topologicalAnalysis"
    if not exe.exists():
        print("  ERROR: topologicalAnalysis not found.")
        return False

    nk, nE = 51, 101
    print(f"  Grid: nk={nk}, nE={nE}, k=[-0.05, 0.05] 1/A, E=[-2.0, 2.0] eV")
    print(f"  eta = 0.002 eV")

    # Run topologicalAnalysis
    with tempfile.TemporaryDirectory() as tmpdir:
        input_cfg = os.path.join(tmpdir, "input.cfg")
        with open(input_cfg, "w") as f:
            f.write(_SPECTRAL_BULK_CONFIG)

        output_dir = os.path.join(tmpdir, "output")
        os.makedirs(output_dir, exist_ok=True)

        try:
            result = subprocess.run(
                [str(exe)],
                cwd=tmpdir,
                capture_output=True,
                text=True,
                timeout=120,
            )
        except subprocess.TimeoutExpired:
            print("  ERROR: topologicalAnalysis timed out after 120s")
            return False

        if result.returncode != 0:
            print(f"  ERROR: topologicalAnalysis returned {result.returncode}")
            if result.stderr:
                print(f"  stderr: {result.stderr[:500]}")
            return False

        spectral_path = os.path.join(output_dir, "spectral_function.dat")
        if not os.path.isfile(spectral_path):
            print("  ERROR: spectral_function.dat not found")
            return False

        k_arr, E_arr, A_2d = _parse_spectral_function(spectral_path)

    if k_arr is None:
        print("  ERROR: Failed to parse spectral_function.dat")
        return False

    print(f"  Parsed spectral function: shape = {A_2d.shape}")
    print(f"  k range: [{k_arr[0]:.4f}, {k_arr[-1]:.4f}] 1/A")
    print(f"  E range: [{E_arr[0]:.4f}, {E_arr[-1]:.4f}] eV")
    print(f"  A(k,E) range: [{A_2d.min():.4f}, {A_2d.max():.4f}] 1/eV")

    # --- Validation: check peaks near band edges ---
    all_pass = True
    Eg_GaAs = 1.519  # eV

    # Find the energy of maximum spectral weight in the CB region
    E_cb_min_idx = np.argmin(np.abs(E_arr - Eg_GaAs))
    cb_region = A_2d[:, max(0, E_cb_min_idx - 5):E_cb_min_idx + 5]
    if cb_region.size > 0:
        cb_peak_E = E_arr[max(0, E_cb_min_idx - 5) + np.unravel_index(
            np.argmax(cb_region), cb_region.shape)[1]]
        cb_err = abs(cb_peak_E - Eg_GaAs)
        if cb_err < 0.3:
            print(f"  CB spectral peak at E={cb_peak_E:.3f} eV "
                  f"(expected ~{Eg_GaAs:.3f}, err={cb_err:.3f} eV) -- PASS")
        else:
            print(f"  CB spectral peak at E={cb_peak_E:.3f} eV "
                  f"(expected ~{Eg_GaAs:.3f}, err={cb_err:.3f} eV) -- FAIL")
            all_pass = False
    else:
        print("  WARNING: CB region empty, skipping CB peak validation")

    # Check VB region (E < 0.5 eV) has spectral weight
    vb_mask = E_arr < 0.5
    if np.any(vb_mask):
        vb_weight = A_2d[:, vb_mask].sum()
        total_weight = A_2d.sum()
        vb_frac = vb_weight / total_weight if total_weight > 0 else 0
        if vb_frac > 0.1:
            print(f"  VB spectral fraction: {vb_frac:.1%} -- PASS")
        else:
            print(f"  VB spectral fraction: {vb_frac:.1%} -- FAIL (expected > 10%)")
            all_pass = False
    else:
        print("  WARNING: no VB energy points found")

    # Check total spectral weight is positive
    if A_2d.min() >= 0:
        print("  Spectral function non-negative -- PASS")
    else:
        print(f"  Spectral function has negative values (min={A_2d.min():.4e}) -- WARNING")

    # --- Generate 2D heatmap ---
    print("\n  Generating spectral function heatmap...")
    fig, ax = plt.subplots(figsize=(10, 7))

    im = ax.pcolormesh(k_arr, E_arr, A_2d.T, cmap='hot_r', shading='auto')
    ax.axhline(0, color='cyan', ls='--', lw=0.8, alpha=0.7, label='VB top (E=0)')
    ax.axhline(Eg_GaAs, color='lime', ls='--', lw=0.8, alpha=0.7,
               label=f'CB edge (E={Eg_GaAs:.3f} eV)')

    ax.set_xlabel(r'$k$ (1/A)', fontsize=12)
    ax.set_ylabel(r'$E$ (eV)', fontsize=12)
    ax.set_title(r'Spectral Function $A(k,E)$ for Bulk GaAs (8-band k.p)'
                 '\n' r'$\eta = 2$ meV, nk=51, nE=101', fontsize=12)
    ax.legend(fontsize=9, loc='upper left')
    fig.colorbar(im, ax=ax, label=r'$A(k,E)$ (1/eV)')

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "lecture_13_spectral.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {FIGURES_DIR / 'lecture_13_spectral.png'}")

    status = "PASS" if all_pass else "FAIL"
    print(f"\n  Section 5 result: {status}")
    return all_pass


# ---------------------------------------------------------------------------
# Section 6 -- Z2 Phase Diagram (BHZ analytic sweep)
# ---------------------------------------------------------------------------

_SWEEP_BHZ_CONFIG = """\
waveVector: k0
waveVectorMax: 0.0
waveVectorStep: 1
confinement: 0
FDstep: 1
FDorder: 2
numLayers: 1
material1: GaAs
numcb: 2
numvb: 6
ExternalField: 0  EF
EFParams: 0.0
SC: 0
topology: T
mode: sweep
compute_chern: F
compute_hall: F
qwz_u: 0.0
compute_z2: F
extract_edge_states: F
edge_E_window: 0.01
compute_ldos: F
ldos_eta: 0.001
ldos_E_range: -0.1
ldos_E_range: 0.1
ldos_num_E: 10
compute_gap_sweep: T
gap_sweep_B: -0.02 0.02 41
gap_sweep_mu: -0.01 0.03 41
sweep_model: bhz_analytic
"""


def _parse_z2_phase_diagram(filepath):
    """Parse z2_phase_diagram.dat into (B_arr, mu_arr, z2_2d, gap_2d).

    Returns:
        (B_unique, mu_unique, z2_2d, gap_2d) or Nones on failure.
    """
    nB = nMu = None
    B_list, mu_list, z2_list, gap_list = [], [], [], []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                m_nB = re.search(r'nB=(\d+)', line)
                if m_nB:
                    nB = int(m_nB.group(1))
                m_nMu = re.search(r'nMu=(\d+)', line)
                if m_nMu:
                    nMu = int(m_nMu.group(1))
                continue
            parts = line.split()
            if len(parts) >= 4:
                B_list.append(float(parts[0]))
                mu_list.append(float(parts[1]))
                z2_list.append(float(parts[2]))
                gap_list.append(float(parts[3]))

    if nB is None or nMu is None or len(B_list) == 0:
        return None, None, None, None

    B_arr = np.array(B_list)
    mu_arr = np.array(mu_list)
    z2_flat = np.array(z2_list)
    gap_flat = np.array(gap_list)

    # Layout: B outer loop, mu inner loop -> reshape (nB, nMu)
    z2_2d = z2_flat.reshape(nB, nMu).T   # shape (nMu, nB)
    gap_2d = gap_flat.reshape(nB, nMu).T

    B_unique = B_arr.reshape(nB, nMu)[:, 0]
    mu_unique = mu_arr.reshape(nB, nMu)[0, :]

    return B_unique, mu_unique, z2_2d, gap_2d


def _parse_z2_transitions(filepath):
    """Parse z2_transitions.dat into list of (B, mu) tuples."""
    transitions = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                transitions.append((float(parts[0]), float(parts[1])))
    return transitions


def section6_phase_diagram():
    """Section 6: Z2 Phase Diagram via BHZ analytic sweep.

    Sweeps (B, mu) parameter space using the BHZ analytic model.
    The Z2 invariant changes from 0 (trivial) to 1 (topological) across
    a critical line in (B, mu) space, where the bulk gap closes.
    """
    print("\n" + "=" * 60)
    print("  Section 6: Z2 Phase Diagram (BHZ Analytic Sweep)")
    print("=" * 60)

    exe = BUILD_DIR / "src" / "topologicalAnalysis"
    if not exe.exists():
        print("  ERROR: topologicalAnalysis not found.")
        return False

    nB, nMu = 41, 41
    print(f"  Grid: nB={nB}, nMu={nMu}")
    print(f"  B range: [-0.02, 0.02], mu range: [-0.01, 0.03] eV")
    print(f"  BHZ: M=10 meV, transition at B = mu - M")

    # Run topologicalAnalysis
    with tempfile.TemporaryDirectory() as tmpdir:
        input_cfg = os.path.join(tmpdir, "input.cfg")
        with open(input_cfg, "w") as f:
            f.write(_SWEEP_BHZ_CONFIG)

        output_dir = os.path.join(tmpdir, "output")
        os.makedirs(output_dir, exist_ok=True)

        try:
            result = subprocess.run(
                [str(exe)],
                cwd=tmpdir,
                capture_output=True,
                text=True,
                timeout=120,
            )
        except subprocess.TimeoutExpired:
            print("  ERROR: topologicalAnalysis timed out after 120s")
            return False

        if result.returncode != 0:
            print(f"  ERROR: topologicalAnalysis returned {result.returncode}")
            if result.stderr:
                print(f"  stderr: {result.stderr[:500]}")
            return False

        phase_path = os.path.join(output_dir, "z2_phase_diagram.dat")
        trans_path = os.path.join(output_dir, "z2_transitions.dat")

        if not os.path.isfile(phase_path):
            print("  ERROR: z2_phase_diagram.dat not found")
            return False

        B_arr, mu_arr, z2_2d, gap_2d = _parse_z2_phase_diagram(phase_path)
        transitions = _parse_z2_transitions(trans_path) if os.path.isfile(trans_path) else []

    if B_arr is None:
        print("  ERROR: Failed to parse z2_phase_diagram.dat")
        return False

    print(f"  Parsed phase diagram: z2 shape = {z2_2d.shape}")
    n_topo = int(z2_2d.sum())
    n_trivial = z2_2d.size - n_topo
    print(f"  Topological (Z2=1): {n_topo} points, Trivial (Z2=0): {n_trivial} points")
    if transitions:
        print(f"  Phase transitions detected: {len(transitions)}")
        for B_t, mu_t in transitions:
            print(f"    B={B_t:.4f} T, mu={mu_t:.6f} eV")

    # --- Validation ---
    all_pass = True

    # Must have both phases present
    if n_topo > 0 and n_trivial > 0:
        print(f"  Both trivial and topological phases present -- PASS")
    else:
        print(f"  Missing phase: topo={n_topo}, trivial={n_trivial} -- FAIL")
        all_pass = False

    # Phase transitions must be detected
    if len(transitions) > 0:
        print(f"  Phase boundary detected ({len(transitions)} points) -- PASS")
    else:
        print("  No phase transitions detected -- FAIL")
        all_pass = False

    # Gap must close somewhere (min gap near zero)
    min_gap = gap_2d.min()
    if min_gap < 0.01:  # eV
        print(f"  Minimum gap in phase diagram: {min_gap:.6f} eV (closes at transition) -- PASS")
    else:
        print(f"  Minimum gap in phase diagram: {min_gap:.6f} eV (expected ~0) -- FAIL")
        all_pass = False

    # --- Generate phase diagram figure ---
    print("\n  Generating Z2 phase diagram...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left panel: gap magnitude colormap with Z2 overlay
    ax1 = axes[0]
    B_mesh, mu_mesh = np.meshgrid(B_arr, mu_arr)
    im = ax1.pcolormesh(B_mesh, mu_mesh * 1000, gap_2d * 1000,
                        cmap='viridis', shading='auto')
    # Contour for Z2 boundary
    ax1.contour(B_mesh, mu_mesh * 1000, z2_2d, levels=[0.5],
                colors='red', linewidths=2)
    fig.colorbar(im, ax=ax1, label='Gap (meV)')
    ax1.set_xlabel('Magnetic field $B$ (T)', fontsize=12)
    ax1.set_ylabel(r'Chemical potential $\mu$ (meV)', fontsize=12)
    ax1.set_title(r'Z$_2$ Phase Diagram: Gap Magnitude' '\n(red = phase boundary)',
                  fontsize=12)

    # Right panel: Z2 invariant map
    ax2 = axes[1]
    cmap_z2 = plt.cm.colors.ListedColormap(['steelblue', 'gold'])
    bounds = [-0.5, 0.5, 1.5]
    norm = plt.cm.colors.BoundaryNorm(bounds, cmap_z2.N)
    im2 = ax2.pcolormesh(B_mesh, mu_mesh * 1000, z2_2d, cmap=cmap_z2,
                         norm=norm, shading='auto')
    ax2.contour(B_mesh, mu_mesh * 1000, z2_2d, levels=[0.5],
                colors='red', linewidths=2)
    cbar2 = fig.colorbar(im2, ax=ax2, ticks=[0, 1])
    cbar2.set_ticklabels([r'Trivial ($Z_2=0$)', r'Topological ($Z_2=1$)'])
    if transitions:
        B_trans = [t[0] for t in transitions]
        mu_trans = [t[1] * 1000 for t in transitions]
        ax2.plot(B_trans, mu_trans, 'r.', ms=8, label='Transition points')
        ax2.legend(fontsize=9)
    ax2.set_xlabel('Magnetic field $B$ (T)', fontsize=12)
    ax2.set_ylabel(r'Chemical potential $\mu$ (meV)', fontsize=12)
    ax2.set_title(r'Z$_2$ Invariant Map (BHZ Analytic)', fontsize=12)

    fig.suptitle('BHZ Topological Phase Diagram', fontsize=14, fontweight='bold',
                 y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "lecture_13_phase_diagram.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {FIGURES_DIR / 'lecture_13_phase_diagram.png'}")

    status = "PASS" if all_pass else "FAIL"
    print(f"\n  Section 6 result: {status}")
    return all_pass


# ---------------------------------------------------------------------------
# Section 7 -- Hall Conductance via Kubo Formula (QWZ model)
# ---------------------------------------------------------------------------

_CONDUCTANCE_CONFIG_TEMPLATE = """\
waveVector: k0
waveVectorMax: 0.0
waveVectorStep: 1
confinement: 0
FDstep: 1
FDorder: 2
numLayers: 1
material1: GaAs
numcb: 2
numvb: 6
ExternalField: 0  EF
EFParams: 0.0
topology: T
mode: conductance
compute_chern: F
compute_hall: F
qwz_u: {u}
compute_z2: F
extract_edge_states: F
edge_E_window: 0.01
compute_ldos: F
ldos_eta: 0.001
ldos_E_range: -0.1
ldos_E_range: 0.1
ldos_num_E: 200
compute_conductance: T
conductance_method: kubo_chern
berry_nk: 50
landauer_energy: 0.0
"""


def _run_conductance(u_val, work_dir):
    """Run topologicalAnalysis in conductance mode for a given QWZ u.

    Returns:
        dict with chern_number, conductance_xy, or None on failure.
    """
    exe = BUILD_DIR / "src" / "topologicalAnalysis"
    config = _CONDUCTANCE_CONFIG_TEMPLATE.format(u=u_val)

    input_cfg = os.path.join(work_dir, "input.cfg")
    with open(input_cfg, "w") as f:
        f.write(config)

    output_dir = os.path.join(work_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    try:
        result = subprocess.run(
            [str(exe)],
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=120,
        )
    except subprocess.TimeoutExpired:
        return None

    if result.returncode != 0:
        return None

    topo_path = os.path.join(output_dir, "topology_result.dat")
    if not os.path.isfile(topo_path):
        return None

    content = Path(topo_path).read_text()

    out = {"chern_number": None, "conductance_xy": None}

    m_chern = re.search(r"# Chern number:\s*(-?\d+)", content)
    if m_chern:
        out["chern_number"] = int(m_chern.group(1))

    m_cond = re.search(r"# Conductance xy \(e\^2/h\):\s*([0-9.e+-]+)", content)
    if m_cond:
        out["conductance_xy"] = float(m_cond.group(1))

    # Fallback: parse Hall conductance
    if out["conductance_xy"] is None:
        m_hall = re.search(r"# Hall conductance \(e\^2/h\):\s*([0-9.e+-]+)", content)
        if m_hall:
            out["conductance_xy"] = float(m_hall.group(1))

    return out


def section7_conductance():
    """Section 7: Hall Conductance via Kubo Formula (QWZ model).

    Computes the quantized Hall conductance sigma_xy = C * e^2/h for the
    QWZ model at u=-0.8 (C=+1, topological) and u=2.5 (C=0, trivial).
    Also computes Berry curvature on a fine Python grid and generates
    a combined figure.
    """
    print("\n" + "=" * 60)
    print("  Section 7: Hall Conductance via Kubo Formula (QWZ)")
    print("=" * 60)

    exe = BUILD_DIR / "src" / "topologicalAnalysis"
    if not exe.exists():
        print("  ERROR: topologicalAnalysis not found.")
        return False

    # Test points: (u, expected_Chern)
    test_points = [
        (-0.8, +1),
        (2.5, 0),
    ]

    results = {}
    all_pass = True

    for u_val, expected_C in test_points:
        print(f"\n  Running u={u_val:+.1f} (expected C={expected_C})...")
        with tempfile.TemporaryDirectory() as tmpdir:
            res = _run_conductance(u_val, tmpdir)

        if res is None or res["chern_number"] is None:
            print(f"    FAIL: Could not parse results for u={u_val}")
            all_pass = False
            results[u_val] = {"C": None, "sigma_xy": None, "expected_C": expected_C}
            continue

        C = res["chern_number"]
        sigma_xy = res["conductance_xy"]
        print(f"    Chern number C = {C}")
        print(f"    sigma_xy = {sigma_xy} e^2/h" if sigma_xy is not None
              else "    sigma_xy = None")

        # Validate Chern number
        if C == expected_C:
            print(f"    Chern number -- PASS")
        else:
            print(f"    Chern number -- FAIL (expected {expected_C})")
            all_pass = False

        # Validate conductance quantization: sigma_xy = C * e^2/h
        # (In units of e^2/h, sigma_xy should equal C exactly)
        if sigma_xy is not None:
            expected_sigma = float(C)
            err = abs(sigma_xy - expected_sigma)
            if err < 0.01:
                print(f"    Conductance quantization: sigma_xy={sigma_xy:.4f} "
                      f"= {C}*e^2/h (err={err:.4f}) -- PASS")
            else:
                print(f"    Conductance quantization: sigma_xy={sigma_xy:.4f} "
                      f"expected {expected_sigma:.4f} (err={err:.4f}) -- FAIL")
                all_pass = False
        else:
            print(f"    Conductance not parsed -- FAIL")
            all_pass = False

        results[u_val] = {"C": C, "sigma_xy": sigma_xy, "expected_C": expected_C}

    # --- Compute Berry curvature on fine grid for u=-0.8 ---
    print("\n  Computing Berry curvature on fine grid (u=-0.8, nk=200)...")
    KX, KY, F = compute_berry_curvature_qwz(u=-0.8, nk=200)
    C_python = np.sum(F) * (2 * np.pi / 200)**2 / (2 * np.pi)
    print(f"  Python Berry curvature integral: C = {C_python:.4f}")

    # --- Generate combined figure ---
    print("\n  Generating conductance figure...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left panel: Berry curvature heatmap
    ax1 = axes[0]
    vmax = np.percentile(np.abs(F), 95)
    im = ax1.pcolormesh(KX, KY, F, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                        shading='auto')
    ax1.set_xlabel(r'$k_x$', fontsize=12)
    ax1.set_ylabel(r'$k_y$', fontsize=12)
    ax1.set_title(r'Berry curvature $F(k)$, QWZ $u=-0.8$'
                  '\n' rf'$C = \int F\,dk_x\,dk_y / 2\pi = {C_python:.2f}$',
                  fontsize=12)
    ax1.set_aspect('equal')
    fig.colorbar(im, ax=ax1, label=r'$F(k)$')

    # Right panel: Bar chart of Hall conductance
    ax2 = axes[1]
    u_vals_plot = list(results.keys())
    u_strs = [f"u={u:+.1f}" for u in u_vals_plot]
    sigma_vals = [results[u]["sigma_xy"] if results[u]["sigma_xy"] is not None else 0
                  for u in u_vals_plot]
    C_vals = [results[u]["C"] if results[u]["C"] is not None else 0
              for u in u_vals_plot]
    expected_C_vals = [results[u]["expected_C"] for u in u_vals_plot]

    bar_colors = ['crimson' if c == 1 else 'forestgreen' for c in expected_C_vals]
    x_pos = np.arange(len(u_vals_plot))
    bars = ax2.bar(x_pos - 0.15, sigma_vals, width=0.3, color=bar_colors,
                   edgecolor='k', lw=0.5, label=r'$\sigma_{xy}$ (computed)')
    bars2 = ax2.bar(x_pos + 0.15, expected_C_vals, width=0.3, color=bar_colors,
                    edgecolor='k', lw=0.5, alpha=0.4, label=r'$C$ (expected)')

    for i, (s, c) in enumerate(zip(sigma_vals, C_vals)):
        ax2.text(i - 0.15, s + 0.05, f'{s:.2f}', ha='center', fontsize=10)
        ax2.text(i + 0.15, c + 0.05, f'{c}', ha='center', fontsize=10)

    ax2.axhline(0, color='gray', ls='--', lw=0.6)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(u_strs, fontsize=11)
    ax2.set_ylabel(r'$\sigma_{xy}$ ($e^2/h$)', fontsize=12)
    ax2.set_title(r'Hall Conductance: $\sigma_{xy} = C \cdot e^2/h$'
                  '\n(Kubo-Chern method, nk=50)', fontsize=12)
    ax2.set_yticks([-1, 0, 1])
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3, axis='y')

    fig.suptitle('Quantized Hall Conductance in QWZ Model', fontsize=14,
                 fontweight='bold', y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "lecture_13_conductance.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {FIGURES_DIR / 'lecture_13_conductance.png'}")

    status = "PASS" if all_pass else "FAIL"
    print(f"\n  Section 7 result: {status}")
    return all_pass


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 60)
    print("  Lecture 13: Topological Superconductivity")
    print("=" * 60)

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    results = {}
    results['chern'] = section1_chern()
    results['z2'] = section2_z2()
    results['majorana'] = section3_majorana()
    results['landau'] = section4_landau()
    results['spectral'] = section5_spectral()
    results['phase_diagram'] = section6_phase_diagram()
    results['conductance'] = section7_conductance()

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
    for f in ["lecture_13_chern.png", "lecture_13_z2.png",
              "lecture_13_majorana.png", "lecture_13_landau.png",
              "lecture_13_spectral.png", "lecture_13_phase_diagram.png",
              "lecture_13_conductance.png"]:
        path = FIGURES_DIR / f
        exists = "OK" if path.exists() else "MISSING"
        print(f"    {f}: {exists}")

    print("\n  Lecture 13 complete." if all_pass else "\n  Lecture 13 had failures.")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
