#!/usr/bin/env python3
"""
generate_all_figures.py
=======================
Run Fortran executables for each test configuration and produce
publication-quality matplotlib figures (PNG, 300 dpi).

Usage
-----
    python scripts/plotting/generate_all_figures.py [--skip-build] [--only FIG ...]

Options
-------
    --skip-build    Skip the cmake build step (assume executables exist).
    --only FIG ...  Only generate the named figures (without extension).

Output
------
    docs/figures/*.png
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
except ImportError:
    print("ERROR: matplotlib and numpy are required.  pip install matplotlib numpy")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Paths (relative to repo root, resolved at import time)
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
BUILD_DIR = REPO_ROOT / "build"
EXE_BAND = BUILD_DIR / "src" / "bandStructure"
EXE_GFACTOR = BUILD_DIR / "src" / "gfactorCalculation"
CONFIG_DIR = REPO_ROOT / "tests" / "regression" / "configs"
FIGURE_DIR = REPO_ROOT / "docs" / "figures"
INPUT_CFG = REPO_ROOT / "input.cfg"

# ---------------------------------------------------------------------------
# Band names and colours (basis ordering: 1-4 valence, 5-6 SO, 7-8 CB)
# ---------------------------------------------------------------------------
BAND_NAMES = [
    "HH (+3/2)",
    "LH (+1/2)",
    "LH (-1/2)",
    "HH (-3/2)",
    "SO (+1/2)",
    "SO (-1/2)",
    "CB (+1/2)",
    "CB (-1/2)",
]

BAND_COLORS = [
    "#d62728",  # HH red
    "#1f77b4",  # LH blue
    "#2ca02c",  # LH green
    "#9467bd",  # HH purple
    "#ff7f0e",  # SO orange
    "#e377c2",  # SO pink
    "#17becf",  # CB cyan
    "#bcbd22",  # CB olive
]

# ---------------------------------------------------------------------------
# Matplotlib style
# ---------------------------------------------------------------------------

def setup_style() -> None:
    """Configure matplotlib for clean academic plots."""
    rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["DejaVu Serif", "Computer Modern", "serif"],
            "font.size": 11,
            "axes.linewidth": 0.8,
            "axes.labelsize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 9,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.05,
            "mathtext.fontset": "cm",
        }
    )


# ===========================================================================
# Utility helpers
# ===========================================================================


def ensure_dirs() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)


def build_fortran() -> None:
    """Configure (if needed) and build the Fortran executables."""
    if EXE_BAND.exists() and EXE_GFACTOR.exists():
        print("[build] Executables already exist, skipping build.")
        return
    print("[build] Configuring and building Fortran code ...")
    if not BUILD_DIR.exists():
        mkl_dir = os.environ.get("MKLROOT", "")
        cmake_cmd = [
            "cmake",
            "-G",
            "Ninja",
            "-B",
            str(BUILD_DIR),
        ]
        if mkl_dir:
            cmake_cmd.append(f"-DMKL_DIR={mkl_dir}/lib/cmake/mkl")
        subprocess.check_call(cmake_cmd, cwd=str(REPO_ROOT))
    subprocess.check_call(["cmake", "--build", str(BUILD_DIR)], cwd=str(REPO_ROOT))
    print("[build] Done.")


def run_executable(
    exe: Path,
    config_path: Path,
    cwd: Path,
    label: str = "",
    timeout: int = 300,
) -> subprocess.CompletedProcess:
    """Copy config to input.cfg, run executable, return CompletedProcess."""
    # Copy config file
    shutil.copy2(str(config_path), str(cwd / "input.cfg"))
    print(f"  [{label}] Running {exe.name} ...")
    t0 = time.time()
    result = subprocess.run(
        [str(exe)],
        cwd=str(cwd),
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    elapsed = time.time() - t0
    if result.returncode != 0:
        print(f"  [{label}] FAILED (rc={result.returncode}) in {elapsed:.1f}s")
        print(f"    stderr: {result.stderr[:500]}")
    else:
        print(f"  [{label}] OK in {elapsed:.1f}s")
    return result


def parse_eigenvalues(output_dir: Path) -> Tuple[np.ndarray, np.ndarray]:
    """
    Parse output/eigenvalues.dat.

    Returns
    -------
    k_vals : 1-D array of |k|
    eig : 2-D array (n_bands, n_kpoints)
    """
    path = output_dir / "eigenvalues.dat"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")
    data = np.loadtxt(str(path), comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    k_vals = data[:, 0]
    eig = data[:, 1:].T  # (n_bands, n_kpoints)
    return k_vals, eig


def parse_parts(output_dir: Path) -> np.ndarray:
    """
    Parse output/parts.dat.

    Returns
    -------
    parts : 2-D array (n_eigenvalues, 8)
    """
    path = output_dir / "parts.dat"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")
    return np.loadtxt(str(path))


def parse_potential_profile(output_dir: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Parse output/potential_profile.dat (QW mode).

    Returns
    -------
    z, EV, EV_DeltaSO, EC : 1-D arrays
    """
    path = output_dir / "potential_profile.dat"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")
    data = np.loadtxt(str(path))
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]


def parse_sc_potential(output_dir: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Parse output/sc_potential_profile.dat."""
    path = output_dir / "sc_potential_profile.dat"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")
    data = np.loadtxt(str(path))
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]


def parse_sc_charge(output_dir: Path, is_2d: bool = False):
    """
    Parse output/sc_charge.dat.

    For QW (1-D): returns (z, n_e, n_h)
    For wire (2-D): returns (x, y, n_e, n_h) grids
    """
    path = output_dir / "sc_charge.dat"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")
    if is_2d:
        # 2D format with blank-line separators for gnuplot splot
        xs, ys, ne, nh = [], [], [], []
        block: List[List[float]] = []
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    if block:
                        for row in block:
                            xs.append(row[0])
                            ys.append(row[1])
                            ne.append(row[2])
                            nh.append(row[3])
                        block = []
                    continue
                block.append([float(v) for v in line.split()])
            if block:
                for row in block:
                    xs.append(row[0])
                    ys.append(row[1])
                    ne.append(row[2])
                    nh.append(row[3])
        return np.array(xs), np.array(ys), np.array(ne), np.array(nh)
    else:
        data = np.loadtxt(str(path), comments="#")
        return data[:, 0], data[:, 1], data[:, 2]


def parse_gfactor(output_dir: Path) -> Tuple[float, float, float]:
    """Parse output/gfactor.dat -> (gx, gy, gz)."""
    path = output_dir / "gfactor.dat"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")
    data = np.loadtxt(str(path))
    return float(data[0]), float(data[1]), float(data[2])


def parse_eigenfunctions_bulk(output_dir: Path, k_idx: int, n_ev: int) -> np.ndarray:
    """
    Parse bulk eigenfunctions for a given k-point.

    Returns
    -------
    components : (n_ev, 8) array of |psi_i|^2 per band
    """
    result = []
    for ev in range(1, n_ev + 1):
        fname = output_dir / f"eigenfunctions_k_{k_idx:05d}_ev_{ev:05d}.dat"
        if not fname.exists():
            continue
        data = np.loadtxt(str(fname))
        result.append(data)
    if not result:
        return np.zeros((0, 8))
    return np.array(result)


def parse_eigenfunctions_qw(
    output_dir: Path, k_idx: int, n_ev: int, n_z: int
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Parse QW eigenfunctions for a given k-point.

    Returns
    -------
    z : (n_z,) positions
    wf : (n_ev, n_z, 8) |psi| per band
    """
    z_arr = None
    wfs = []
    for ev in range(1, n_ev + 1):
        fname = output_dir / f"eigenfunctions_k_{k_idx:05d}_ev_{ev:05d}.dat"
        if not fname.exists():
            continue
        data = np.loadtxt(str(fname))
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if z_arr is None:
            z_arr = data[:, 0]
        wfs.append(data[:, 1:])  # (n_z, 8)
    if not wfs or z_arr is None:
        return np.array([]), np.zeros((0, 0, 8))
    return z_arr, np.array(wfs)  # (n_ev, n_z, 8)


def parse_wire_eigenfunction_2d(
    output_dir: Path, k_idx: int, ev_idx: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Parse wire 2-D eigenfunction.

    Returns
    -------
    x, y : 1-D coordinate arrays
    psi2 : 2-D |psi(x,y)|^2
    """
    fname = output_dir / f"eigenfunctions_k_{k_idx:05d}_ev_{ev_idx:05d}.dat"
    if not fname.exists():
        raise FileNotFoundError(f"{fname} not found")
    xs, ys, psi2 = [], [], []
    row_block: List[List[float]] = []
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            if not line:
                if row_block:
                    for row in row_block:
                        xs.append(row[0])
                        ys.append(row[1])
                        psi2.append(row[2])
                    row_block = []
                continue
            row_block.append([float(v) for v in line.split()])
        if row_block:
            for row in row_block:
                xs.append(row[0])
                ys.append(row[1])
                psi2.append(row[2])
    xs, ys, psi2 = np.array(xs), np.array(ys), np.array(psi2)
    if len(xs) == 0:
        return xs, ys, psi2
    x_unique = np.unique(xs)
    y_unique = np.unique(ys)
    psi2_grid = psi2.reshape(len(y_unique), len(x_unique))
    return x_unique, y_unique, psi2_grid


def parse_convergence_from_stdout(stdout: str) -> Tuple[List[int], List[float]]:
    """
    Parse SC convergence history from stdout lines like:
      iter:   1  |dPhi|:  1.2345E-03   ...
    """
    iters: List[int] = []
    deltas: List[float] = []
    for line in stdout.splitlines():
        line = line.strip()
        if "iter:" in line and "|dPhi|:" in line:
            try:
                i_pos = line.index("iter:") + len("iter:")
                d_pos = line.index("|dPhi|:") + len("|dPhi|:")
                iter_val = int(line[i_pos:].split()[0])
                rest = line[d_pos:].split()
                delta_val = float(rest[0])
                iters.append(iter_val)
                deltas.append(delta_val)
            except (ValueError, IndexError):
                continue
    return iters, deltas


# ===========================================================================
# Figure generators
# ===========================================================================


def fig_bulk_gaas_bands(output_dir: Path) -> None:
    """bulk_gaas_bands.png: 8-band E(k) dispersion from bulk GaAs."""
    print("[figure] bulk_gaas_bands")
    cfg = CONFIG_DIR / "bulk_gaas_kx.cfg"
    run_executable(EXE_BAND, cfg, REPO_ROOT, label="bulk_gaas_kx")
    k_vals, eig = parse_eigenvalues(output_dir)

    fig, ax = plt.subplots(figsize=(5, 4.5))
    n_bands = eig.shape[0]
    for i in range(n_bands):
        ax.plot(k_vals, eig[i], color=BAND_COLORS[i], linewidth=1.2)
    ax.set_xlabel(r"$k$ (1/A)")
    ax.set_ylabel(r"$E$ (eV)")
    ax.set_title("Bulk GaAs 8-band dispersion")
    ax.axhline(0, color="grey", linewidth=0.4, linestyle="--")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "bulk_gaas_bands.png")
    plt.close(fig)
    print("  -> docs/figures/bulk_gaas_bands.png")


def fig_bulk_gaas_parts(output_dir: Path) -> None:
    """bulk_gaas_parts.png: band character at Gamma from bulk."""
    print("[figure] bulk_gaas_parts")
    parts = parse_parts(output_dir)
    if parts.size == 0:
        # Re-run if needed (uses output from bulk_gaas_bands run)
        cfg = CONFIG_DIR / "bulk_gaas_kx.cfg"
        run_executable(EXE_BAND, cfg, REPO_ROOT, label="bulk_gaas_kx")
        parts = parse_parts(output_dir)

    n_ev = parts.shape[0]
    fig, ax = plt.subplots(figsize=(6, 4))
    bar_width = 0.8
    x = np.arange(n_ev)
    bottom = np.zeros(n_ev)
    for b in range(8):
        ax.bar(
            x,
            parts[:, b],
            bar_width,
            bottom=bottom,
            label=BAND_NAMES[b],
            color=BAND_COLORS[b],
            linewidth=0,
        )
        bottom += parts[:, b]
    ax.set_xlabel("Eigenstate index")
    ax.set_ylabel("Band character")
    ax.set_title("Bulk GaAs band decomposition at $\\Gamma$")
    ax.legend(ncol=2, fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "bulk_gaas_parts.png")
    plt.close(fig)
    print("  -> docs/figures/bulk_gaas_parts.png")


def fig_qw_alsbw_gasbw_inasw_bands(output_dir: Path) -> None:
    """qw_alsbw_gasbw_inasw_bands.png: E(k_parallel) subbands from QW."""
    print("[figure] qw_alsbw_gasbw_inasw_bands")
    cfg = CONFIG_DIR / "qw_alsbw_gasbw_inasw.cfg"
    run_executable(EXE_BAND, cfg, REPO_ROOT, label="qw_alsbw_gasbw_inasw")
    k_vals, eig = parse_eigenvalues(output_dir)

    fig, ax = plt.subplots(figsize=(5.5, 5))
    n_bands = eig.shape[0]
    # Separate CB and VB by sign
    for i in range(n_bands):
        energy_mid = np.mean(eig[i])
        if energy_mid > 0:
            color = "#17becf"
            alpha = 0.85
        else:
            color = "#d62728"
            alpha = 0.85
        ax.plot(k_vals, eig[i], color=color, linewidth=0.9, alpha=alpha)
    ax.set_xlabel(r"$k_{\parallel}$ (1/A)")
    ax.set_ylabel(r"$E$ (eV)")
    ax.set_title("AlSb/GaSb/InAs QW subbands")
    ax.axhline(0, color="grey", linewidth=0.4, linestyle="--")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "qw_alsbw_gasbw_inasw_bands.png")
    plt.close(fig)
    print("  -> docs/figures/qw_alsbw_gasbw_inasw_bands.png")


def fig_qw_potential_profile(output_dir: Path) -> None:
    """qw_potential_profile.png: potential profile from QW."""
    print("[figure] qw_potential_profile")
    # Uses output from the QW run if potential_profile.dat exists
    try:
        z, EV, EV_SO, EC = parse_potential_profile(output_dir)
    except FileNotFoundError:
        cfg = CONFIG_DIR / "qw_alsbw_gasbw_inasw.cfg"
        run_executable(EXE_BAND, cfg, REPO_ROOT, label="qw_alsbw_gasbw_inasw")
        z, EV, EV_SO, EC = parse_potential_profile(output_dir)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(z, EV, color="#d62728", linewidth=1.5, label="$E_V$")
    ax.plot(z, EV_SO, color="#ff7f0e", linewidth=1.5, label="$E_{\\Delta SO}$")
    ax.plot(z, EC, color="#17becf", linewidth=1.5, label="$E_C$")
    ax.fill_between(z, EV.min() - 0.1, EV, alpha=0.06, color="#d62728")
    ax.fill_between(z, EC, EC.max() + 0.1, alpha=0.06, color="#17becf")
    ax.set_xlabel(r"$z$ (A)")
    ax.set_ylabel("Energy (eV)")
    ax.set_title("QW band-edge profile (AlSb/GaSb/InAs)")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "qw_potential_profile.png")
    plt.close(fig)
    print("  -> docs/figures/qw_potential_profile.png")


def fig_qw_wavefunctions(output_dir: Path) -> None:
    """qw_wavefunctions.png: |psi(z)|^2 per band for lowest CB states."""
    print("[figure] qw_wavefunctions")
    # Use eigenfunctions from the QW run
    n_z = 101  # FDstep from config
    try:
        z, wf = parse_eigenfunctions_qw(output_dir, k_idx=1, n_ev=4, n_z=n_z)
    except FileNotFoundError:
        cfg = CONFIG_DIR / "qw_alsbw_gasbw_inasw.cfg"
        run_executable(EXE_BAND, cfg, REPO_ROOT, label="qw_alsbw_gasbw_inasw")
        z, wf = parse_eigenfunctions_qw(output_dir, k_idx=1, n_ev=4, n_z=n_z)

    if z.size == 0:
        print("  WARNING: no wavefunction data found, skipping.")
        return

    n_show = min(wf.shape[0], 4)
    fig, axes = plt.subplots(1, n_show, figsize=(3.5 * n_show, 4), sharey=True)
    if n_show == 1:
        axes = [axes]
    for idx in range(n_show):
        ax = axes[idx]
        psi2_total = np.sum(wf[idx] ** 2, axis=1)
        ax.plot(z, psi2_total, color="#17becf", linewidth=1.2)
        ax.fill_between(z, 0, psi2_total, alpha=0.15, color="#17becf")
        ax.set_xlabel(r"$z$ (A)")
        ax.set_title(f"State {idx + 1}")
    axes[0].set_ylabel(r"$|\psi(z)|^2$")
    fig.suptitle("QW probability density (AlSb/GaSb/InAs, $k_{\\parallel}=0$)", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "qw_wavefunctions.png")
    plt.close(fig)
    print("  -> docs/figures/qw_wavefunctions.png")


def fig_qw_parts(output_dir: Path) -> None:
    """qw_parts.png: integrated band character bar chart from QW."""
    print("[figure] qw_parts")
    try:
        parts = parse_parts(output_dir)
    except FileNotFoundError:
        cfg = CONFIG_DIR / "qw_alsbw_gasbw_inasw.cfg"
        run_executable(EXE_BAND, cfg, REPO_ROOT, label="qw_alsbw_gasbw_inasw")
        parts = parse_parts(output_dir)

    if parts.size == 0:
        print("  WARNING: no parts data, skipping.")
        return

    # Show up to 8 states (first few are usually CB)
    n_show = min(parts.shape[0], 8)
    parts = parts[:n_show]

    fig, ax = plt.subplots(figsize=(7, 4))
    x = np.arange(n_show)
    bottom = np.zeros(n_show)
    # Group into: HH (0,3), LH (1,2), SO (4,5), CB (6,7)
    groups = [
        ("HH", [0, 3], "#d62728"),
        ("LH", [1, 2], "#1f77b4"),
        ("SO", [4, 5], "#ff7f0e"),
        ("CB", [6, 7], "#17becf"),
    ]
    for name, indices, color in groups:
        vals = sum(parts[:, i] for i in indices)
        ax.bar(x, vals, 0.7, bottom=bottom, label=name, color=color, linewidth=0)
        bottom += vals
    ax.set_xlabel("Eigenstate index")
    ax.set_ylabel("Band character")
    ax.set_title("QW band decomposition (AlSb/GaSb/InAs, $k_{\\parallel}=0$)")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "qw_parts.png")
    plt.close(fig)
    print("  -> docs/figures/qw_parts.png")


def fig_gfactor_components(output_dir: Path) -> None:
    """gfactor_components.png: g-factor tensor components for CB and VB."""
    print("[figure] gfactor_components")
    gx_vals: Dict[str, float] = {}
    gy_vals: Dict[str, float] = {}
    gz_vals: Dict[str, float] = {}

    gfactor_configs = [
        ("CB bulk GaAs", "gfactor_bulk_gaas_cb.cfg"),
        ("VB bulk GaAs", "gfactor_bulk_gaas_vb.cfg"),
        ("CB QW", "gfactor_qw_cb.cfg"),
        ("VB QW", "gfactor_qw_vb.cfg"),
    ]

    for label, cfg_name in gfactor_configs:
        cfg_path = CONFIG_DIR / cfg_name
        if not cfg_path.exists():
            print(f"  WARNING: {cfg_name} not found, skipping.")
            continue
        result = run_executable(EXE_GFACTOR, cfg_path, REPO_ROOT, label=label)
        if result.returncode != 0:
            print(f"  WARNING: {label} failed, skipping.")
            continue
        try:
            gx, gy, gz = parse_gfactor(output_dir)
            gx_vals[label] = gx
            gy_vals[label] = gy
            gz_vals[label] = gz
        except FileNotFoundError:
            print(f"  WARNING: gfactor.dat not found for {label}, skipping.")

    if not gx_vals:
        print("  WARNING: no g-factor data produced, skipping figure.")
        return

    labels = list(gx_vals.keys())
    x_pos = np.arange(len(labels))
    width = 0.25

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(x_pos - width, [gx_vals[l] for l in labels], width, label=r"$g_x$", color="#1f77b4")
    ax.bar(x_pos, [gy_vals[l] for l in labels], width, label=r"$g_y$", color="#ff7f0e")
    ax.bar(x_pos + width, [gz_vals[l] for l in labels], width, label=r"$g_z$", color="#2ca02c")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, rotation=15, ha="right")
    ax.set_ylabel("g-factor")
    ax.set_title("Landau g-factor tensor components")
    ax.legend(loc="best")
    ax.axhline(0, color="grey", linewidth=0.5)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "gfactor_components.png")
    plt.close(fig)
    print("  -> docs/figures/gfactor_components.png")


def fig_gfactor_zeeman(output_dir: Path) -> None:
    """gfactor_zeeman.png: Zeeman splitting diagram (schematic)."""
    print("[figure] gfactor_zeeman")
    # Re-run CB bulk GaAs to get gfactor
    cfg_path = CONFIG_DIR / "gfactor_bulk_gaas_cb.cfg"
    if cfg_path.exists():
        result = run_executable(EXE_GFACTOR, cfg_path, REPO_ROOT, label="gfactor_cb_bulk")
        try:
            gx, gy, gz = parse_gfactor(output_dir)
        except FileNotFoundError:
            gx, gy, gz = -0.5, -0.5, -0.5  # fallback schematic
    else:
        gx, gy, gz = -0.5, -0.5, -0.5

    g = gz  # use z-component for the schematic

    fig, ax = plt.subplots(figsize=(5, 4))
    B_values = np.linspace(0, 10, 100)  # Tesla
    E_up = 0.5 * g * B_values   # m_j = +1/2
    E_dn = -0.5 * g * B_values  # m_j = -1/2

    ax.plot(B_values, E_up, color="#d62728", linewidth=1.5, label=r"$m_J = +\frac{1}{2}$")
    ax.plot(B_values, E_dn, color="#1f77b4", linewidth=1.5, label=r"$m_J = -\frac{1}{2}$")
    ax.axhline(0, color="grey", linewidth=0.4, linestyle="--")
    ax.set_xlabel("Magnetic field $B$ (T)")
    ax.set_ylabel("Zeeman energy (meV)")
    ax.set_title(f"Zeeman splitting ($g_z = {g:.3f}$)")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "gfactor_zeeman.png")
    plt.close(fig)
    print("  -> docs/figures/gfactor_zeeman.png")


def fig_sc_potential(output_dir: Path) -> None:
    """sc_potential.png: band-edge profile from QW heterostructure.

    Uses the W-variant QW (AlSbW/GaSbW/InAsW) which has well-defined
    band offsets. If the SC run produced ``sc_potential_profile.dat`` with
    meaningful band bending it will be used; otherwise the base profile
    from ``potential_profile.dat`` is shown.
    """
    print("[figure] sc_potential")

    # Remove any stale sc_potential_profile.dat that may have been produced
    # by a prior wire SC run (5-column 2D format). The QW run below will
    # produce a fresh potential_profile.dat with 4-column 1D format.
    sc_path = output_dir / "sc_potential_profile.dat"
    if sc_path.exists():
        sc_path.unlink()
        print("  Removed stale sc_potential_profile.dat")

    # Run QW first to get the heterostructure profile
    cfg = CONFIG_DIR / "qw_alsbw_gasbw_inasw.cfg"
    run_executable(EXE_BAND, cfg, REPO_ROOT, label="qw_potential")

    # potential_profile.dat should now be 4-column 1D (QW format)
    try:
        z, EV, EV_SO, EC = parse_potential_profile(output_dir)
        title = "Band-edge profile (AlSbW/GaSbW/InAsW QW)"
    except (FileNotFoundError, ValueError) as exc:
        print(f"  WARNING: could not parse potential profile ({exc}), skipping.")
        return

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(z, EV, color="#d62728", linewidth=1.5, label="$E_V$")
    ax.plot(z, EV_SO, color="#ff7f0e", linewidth=1.5, label="$E_{\\Delta SO}$")
    ax.plot(z, EC, color="#17becf", linewidth=1.5, label="$E_C$")
    ax.fill_between(z, EV.min() - 0.1, EV, alpha=0.06, color="#d62728")
    ax.fill_between(z, EC, EC.max() + 0.1, alpha=0.06, color="#17becf")
    ax.set_xlabel(r"$z$ (A)")
    ax.set_ylabel("Energy (eV)")
    ax.set_title(title)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "sc_potential.png")
    plt.close(fig)
    print("  -> docs/figures/sc_potential.png")


def fig_sc_charge_density(output_dir: Path) -> None:
    """sc_charge_density.png: n(z) charge density from SC run."""
    print("[figure] sc_charge_density")
    # Assumes SC run already happened (from fig_sc_potential)
    try:
        z, n_e, n_h = parse_sc_charge(output_dir, is_2d=False)
    except FileNotFoundError:
        print("  WARNING: sc_charge.dat not found. Run SC first.")
        return

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(z, n_e, color="#17becf", linewidth=1.5, label="$n_e(z)$")
    ax.plot(z, n_h, color="#d62728", linewidth=1.5, label="$n_h(z)$")
    ax.set_xlabel(r"$z$ (A)")
    ax.set_ylabel(r"Carrier density (cm$^{-3}$)")
    ax.set_title("Self-consistent charge density (GaAs/AlAs QW)")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "sc_charge_density.png")
    plt.close(fig)
    print("  -> docs/figures/sc_charge_density.png")


def fig_sc_convergence(output_dir: Path) -> None:
    """sc_convergence.png: SC loop convergence history."""
    print("[figure] sc_convergence")
    # Re-run SC to capture stdout
    cfg = CONFIG_DIR / "sc_gaas_alas_qw.cfg"
    result = run_executable(EXE_BAND, cfg, REPO_ROOT, label="sc_convergence", timeout=600)
    if result.returncode != 0:
        print("  WARNING: SC run failed, skipping convergence plot.")
        return

    iters, deltas = parse_convergence_from_stdout(result.stdout)
    if not iters:
        print("  WARNING: no convergence data parsed from stdout, skipping.")
        return

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.semilogy(iters, deltas, "o-", color="#1f77b4", linewidth=1.2, markersize=4)
    ax.set_xlabel("Iteration")
    ax.set_ylabel(r"$|\Delta\phi|$ (eV)")
    ax.set_title("SC loop convergence (GaAs/AlAs QW)")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "sc_convergence.png")
    plt.close(fig)
    print("  -> docs/figures/sc_convergence.png")


def fig_wire_subbands(output_dir: Path) -> None:
    """wire_subbands.png: valence subbands E(k_z) for GaAs rectangular wire."""
    print("[figure] wire_subbands")
    cfg = CONFIG_DIR / "wire_gaas_rectangle.cfg"
    run_executable(EXE_BAND, cfg, REPO_ROOT, label="wire_gaas", timeout=300)
    try:
        k_vals, eig = parse_eigenvalues(output_dir)
    except FileNotFoundError:
        print("  WARNING: eigenvalues.dat not found for wire, skipping.")
        return

    # GaAs valence band edge (all eigenvalues are deep VB states)
    EV_GAAS = -0.80

    fig, ax = plt.subplots(figsize=(5, 4.5))
    n_bands = eig.shape[0]
    vb_colors = plt.cm.Blues_r(np.linspace(0.3, 0.9, n_bands))
    for i in range(n_bands):
        ax.plot(k_vals, eig[i], color=vb_colors[i], linewidth=0.9)
    ax.axhline(EV_GAAS, color="grey", linewidth=0.8, linestyle="--",
               label=f"$E_V$ = {EV_GAAS:.2f} eV")
    ax.set_xlabel(r"$k_z$ (1/A)")
    ax.set_ylabel(r"$E$ (eV)")
    ax.set_title("GaAs rectangular wire valence subbands")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "wire_subbands.png")
    plt.close(fig)
    print("  -> docs/figures/wire_subbands.png")


def _find_ev_idx_nearest_to_energy(output_dir: Path, target_eV: float) -> int:
    """Return the eigenfunction index (1-based) whose eigenvalue at k=1 is
    closest to *target_eV* by reading the comment header of each file."""
    for ev_idx in range(1, 100):
        fname = output_dir / f"eigenfunctions_k_00001_ev_{ev_idx:05d}.dat"
        if not fname.exists():
            break
    n_ev = ev_idx - 1  # last existing index
    if n_ev == 0:
        return 1
    # Read eigenvalues from eigenvalues.dat row 1 (k=1)
    try:
        _, eig = parse_eigenvalues(output_dir)
        e_at_k1 = eig[:, 0]  # first k-point
    except (FileNotFoundError, IndexError):
        return n_ev  # fall back to highest VB state
    best_idx = int(np.argmin(np.abs(e_at_k1 - target_eV))) + 1  # 1-based
    return best_idx


def fig_wire_density_2d(output_dir: Path) -> None:
    """wire_density_2d.png: 2D |psi(x,y)|^2 for the VB-edge state."""
    print("[figure] wire_density_2d")

    # GaAs valence band edge
    EV_GAAS = -0.80

    # Find eigenfunction closest to VB edge (not the deepest state)
    ev_idx = _find_ev_idx_nearest_to_energy(output_dir, EV_GAAS)
    print(f"  Using eigenfunction ev_idx={ev_idx} (closest to EV={EV_GAAS} eV)")

    try:
        x, y, psi2 = parse_wire_eigenfunction_2d(output_dir, k_idx=1, ev_idx=ev_idx)
    except FileNotFoundError:
        cfg = CONFIG_DIR / "wire_gaas_rectangle.cfg"
        run_executable(EXE_BAND, cfg, REPO_ROOT, label="wire_gaas", timeout=300)
        ev_idx = _find_ev_idx_nearest_to_energy(output_dir, EV_GAAS)
        x, y, psi2 = parse_wire_eigenfunction_2d(output_dir, k_idx=1, ev_idx=ev_idx)

    if x.size == 0:
        print("  WARNING: no wire eigenfunction data, skipping.")
        return

    # Read the energy of the chosen state for the title
    state_energy = ""
    ev_file = output_dir / f"eigenfunctions_k_00001_ev_{ev_idx:05d}.dat"
    if ev_file.exists():
        with open(ev_file) as f:
            for line in f:
                if line.startswith("# E ="):
                    state_energy = line.replace("# E =", "").strip()
                    break

    fig, ax = plt.subplots(figsize=(5, 4.5))
    im = ax.pcolormesh(x, y, psi2, shading="auto", cmap="viridis")
    ax.set_xlabel(r"$x$ (A)")
    ax.set_ylabel(r"$y$ (A)")
    title = r"$|\psi(x,y)|^2$ VB-edge state"
    if state_energy:
        title += f" (E = {state_energy} eV)"
    ax.set_title(title)
    ax.set_aspect("equal")
    fig.colorbar(im, ax=ax, label=r"$|\psi|^2$")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "wire_density_2d.png")
    plt.close(fig)
    print("  -> docs/figures/wire_density_2d.png")


def fig_convergence_fd_order(output_dir: Path) -> None:
    """convergence_fd_order.png: convergence vs FD order (2,4,6,8).

    Uses a QW config because FD order only affects the finite-difference
    discretisation of the z-derivative, which is absent in bulk mode.
    """
    print("[figure] convergence_fd_order")
    # Use a QW config — FD order matters for confined structures
    base_cfg = CONFIG_DIR / "qw_alsbw_gasbw_inasw.cfg"
    if not base_cfg.exists():
        print("  WARNING: qw_alsbw_gasbw_inasw.cfg not found, skipping.")
        return

    base_text = base_cfg.read_text()

    fd_orders = [2, 4, 6, 8]
    energies_at_gamma: Dict[int, np.ndarray] = {}

    for order in fd_orders:
        # Modify FDorder line in config text; also reduce numcb/numvb for speed
        lines = base_text.splitlines()
        modified_lines = []
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("FDorder"):
                modified_lines.append(f"FDorder: {order}")
            elif stripped.startswith("waveVectorMax"):
                modified_lines.append("waveVectorMax: 0.0")
            elif stripped.startswith("waveVectorStep"):
                modified_lines.append("waveVectorStep: 1")
            elif stripped.startswith("waveVector:"):
                modified_lines.append("waveVector: k0")
            elif stripped.startswith("numcb"):
                modified_lines.append("numcb: 8")
            elif stripped.startswith("numvb"):
                modified_lines.append("numvb: 8")
            else:
                modified_lines.append(line)
        modified_text = "\n".join(modified_lines) + "\n"

        # Write temp config
        tmp_cfg = REPO_ROOT / "input.cfg"
        tmp_cfg.write_text(modified_text)

        print(f"  [FDorder={order}] Running QW ...")
        result = subprocess.run(
            [str(EXE_BAND)],
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
            timeout=120,
        )
        if result.returncode != 0:
            print(f"  [FDorder={order}] FAILED, skipping.")
            continue

        try:
            k_vals, eig = parse_eigenvalues(output_dir)
            # Take eigenvalues at k=0 (first and only k-point)
            energies_at_gamma[order] = eig[:, 0]
        except FileNotFoundError:
            print(f"  [FDorder={order}] eigenvalues.dat not found, skipping.")
            continue

    if len(energies_at_gamma) < 2:
        print("  WARNING: not enough FD order results for convergence plot, skipping.")
        return

    fig, ax = plt.subplots(figsize=(6, 4.5))
    orders = sorted(energies_at_gamma.keys())
    # Reference = highest order
    ref_order = max(orders)
    ref_e = energies_at_gamma[ref_order]
    for band_idx in range(min(ref_e.shape[0], 8)):
        diffs = []
        valid_orders = []
        for order in orders:
            if order == ref_order:
                continue  # skip reference — diff is zero, breaks log scale
            e = energies_at_gamma[order]
            if band_idx < e.shape[0]:
                diffs.append(abs(e[band_idx] - ref_e[band_idx]))
                valid_orders.append(order)
        if len(valid_orders) >= 2:
            ax.semilogy(valid_orders, diffs, "o-", linewidth=1.2, markersize=5,
                        label=f"Band {band_idx + 1}", color=BAND_COLORS[band_idx])

    ax.set_xlabel("FD order")
    ax.set_ylabel(r"$|E_n - E_n^{\rm ref}|$ (eV)")
    ax.set_title(f"FD order convergence (ref: order {ref_order})")
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "convergence_fd_order.png")
    plt.close(fig)
    print("  -> docs/figures/convergence_fd_order.png")


def fig_timing_dense_vs_sparse(output_dir: Path) -> None:
    """timing_dense_vs_sparse.png: timing comparison dense vs sparse."""
    print("[figure] timing_dense_vs_sparse")
    # Run QW (dense) and wire (sparse) and compare timing from stdout.
    # This is approximate -- we measure wall-clock time.

    results: Dict[str, float] = {}

    # Dense: QW with AlSb/GaSb/InAs
    cfg_qw = CONFIG_DIR / "qw_alsbw_gasbw_inasw.cfg"
    if cfg_qw.exists():
        t0 = time.time()
        r = run_executable(EXE_BAND, cfg_qw, REPO_ROOT, label="QW dense")
        t_dense = time.time() - t0
        if r.returncode == 0:
            results["QW (dense)"] = t_dense

    # Sparse: wire
    cfg_wire = CONFIG_DIR / "wire_gaas_rectangle.cfg"
    if cfg_wire.exists():
        t0 = time.time()
        r = run_executable(EXE_BAND, cfg_wire, REPO_ROOT, label="Wire (sparse)")
        t_sparse = time.time() - t0
        if r.returncode == 0:
            results["Wire (sparse)"] = t_sparse

    if not results:
        print("  WARNING: no timing data collected, skipping figure.")
        return

    fig, ax = plt.subplots(figsize=(5, 4))
    labels = list(results.keys())
    times = [results[l] for l in labels]
    colors = ["#1f77b4", "#ff7f0e"][: len(labels)]
    bars = ax.bar(labels, times, color=colors, width=0.5, linewidth=0)
    ax.set_ylabel("Wall-clock time (s)")
    ax.set_title("Dense vs sparse solver timing")
    # Add time labels on bars
    for bar, t in zip(bars, times):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.05 * max(times),
            f"{t:.1f}s",
            ha="center",
            va="bottom",
            fontsize=10,
        )
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "timing_dense_vs_sparse.png")
    plt.close(fig)
    print("  -> docs/figures/timing_dense_vs_sparse.png")


# ===========================================================================
# Main
# ===========================================================================

ALL_FIGURES = {
    "bulk_gaas_bands": fig_bulk_gaas_bands,
    "bulk_gaas_parts": fig_bulk_gaas_parts,
    "qw_alsbw_gasbw_inasw_bands": fig_qw_alsbw_gasbw_inasw_bands,
    "qw_potential_profile": fig_qw_potential_profile,
    "qw_wavefunctions": fig_qw_wavefunctions,
    "qw_parts": fig_qw_parts,
    "gfactor_components": fig_gfactor_components,
    "gfactor_zeeman": fig_gfactor_zeeman,
    "sc_potential": fig_sc_potential,
    "sc_charge_density": fig_sc_charge_density,
    "sc_convergence": fig_sc_convergence,
    "wire_subbands": fig_wire_subbands,
    "wire_density_2d": fig_wire_density_2d,
    "convergence_fd_order": fig_convergence_fd_order,
    "timing_dense_vs_sparse": fig_timing_dense_vs_sparse,
}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate all publication figures from Fortran output.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            f"""\
            Available figures:
              {chr(10).join(f'  - {name}' for name in ALL_FIGURES)}

            Output directory: {FIGURE_DIR}
            """
        ),
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Skip cmake build step (assume executables exist).",
    )
    parser.add_argument(
        "--only",
        nargs="+",
        metavar="FIG",
        help="Only generate the named figures (without .png extension).",
    )
    args = parser.parse_args()

    setup_style()
    ensure_dirs()

    if not args.skip_build:
        build_fortran()

    output_dir = REPO_ROOT / "output"

    # Determine which figures to generate
    if args.only:
        to_run = {}
        for name in args.only:
            if name in ALL_FIGURES:
                to_run[name] = ALL_FIGURES[name]
            else:
                print(f"WARNING: unknown figure '{name}'. "
                      f"Available: {', '.join(ALL_FIGURES.keys())}")
    else:
        to_run = ALL_FIGURES

    n_total = len(to_run)
    n_ok = 0
    n_fail = 0

    print(f"\n{'=' * 60}")
    print(f"Generating {n_total} figure(s)")
    print(f"{'=' * 60}\n")

    for name, func in to_run.items():
        try:
            func(output_dir)
            n_ok += 1
        except Exception as exc:
            print(f"  ERROR in {name}: {exc}")
            n_fail += 1
        print()

    print(f"{'=' * 60}")
    print(f"Done: {n_ok} succeeded, {n_fail} failed out of {n_total}")
    print(f"Output: {FIGURE_DIR}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
