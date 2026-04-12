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


def parse_parts(output_dir: Path, k_index: int = 0) -> np.ndarray:
    """
    Parse output/parts.dat in multi-block gnuplot format.

    Each block starts with '# k = <value>' followed by n_eigenvalues rows of
    8 band-character weights.  Blocks are separated by blank lines.

    Parameters
    ----------
    output_dir : Path
    k_index : int
        Which k-point block to return (0 = first k-point, i.e. k=0).

    Returns
    -------
    parts : 2-D array (n_eigenvalues, 8)
    """
    path = output_dir / "parts.dat"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")

    blocks: list[np.ndarray] = []
    current: list[list[float]] = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                # New k-point block header
                if current:
                    blocks.append(np.array(current))
                    current = []
                continue
            if line == "":
                if current:
                    blocks.append(np.array(current))
                    current = []
                continue
            vals = [float(x) for x in line.split()]
            current.append(vals)

    if current:
        blocks.append(np.array(current))

    if not blocks:
        raise ValueError(f"No data blocks found in {path}")

    if k_index >= len(blocks):
        raise IndexError(
            f"Requested k_index={k_index} but only {len(blocks)} blocks in {path}"
        )

    return blocks[k_index]


def parse_parts_all_k(output_dir: Path) -> Tuple[np.ndarray, List[float]]:
    """
    Parse output/parts.dat returning ALL k-point blocks.

    Returns
    -------
    k_values : list of float
        k-magnitude for each block.
    all_parts : list of 2-D arrays (n_eigenvalues, 8)
        Band character for each k-point.
    """
    path = output_dir / "parts.dat"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")

    k_values: list[float] = []
    all_parts: list[np.ndarray] = []
    current: list[list[float]] = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("# k ="):
                # Save previous block
                if current:
                    all_parts.append(np.array(current))
                    current = []
                # Parse k value
                k_val = float(line.split("=")[1].strip())
                k_values.append(k_val)
                continue
            if line.startswith("#"):
                continue
            if line == "":
                if current:
                    all_parts.append(np.array(current))
                    current = []
                continue
            vals = [float(x) for x in line.split()]
            current.append(vals)

    if current:
        all_parts.append(np.array(current))

    return k_values, all_parts


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


def parse_optical_transitions(output_dir: Path) -> Tuple[np.ndarray, ...]:
    """Parse output/optical_transitions.dat.

    Handles Fortran G-edit descriptor output where very small exponents
    may omit the 'E' character (e.g. ``0.300918-110`` means ``0.300918E-110``).

    Returns
    -------
    cb_idx, vb_idx, energy, px, py, pz, f_osc : 1-D arrays
    """
    import re
    path = output_dir / "optical_transitions.dat"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")
    rows: list[list[float]] = []
    _fort_fix = re.compile(r"(\d\.\d+E[+-]?\d+)|(\d\.\d+)([+-]\d{2,3})(?=\s|$)")
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.strip() == "":
                continue
            # Fix missing-E in Fortran G-edit output: e.g. 0.300918-110 -> 0.300918E-110
            line = re.sub(r"(\d\.\d+)([+-]\d{2,3})(?=\s|$)", r"\1E\2", line)
            tokens = line.split()
            if len(tokens) >= 7:
                rows.append([float(t) for t in tokens[:7]])
    data = np.array(rows)
    return (data[:, 0].astype(int), data[:, 1].astype(int),
            data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6])


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

    fig, (ax_full, ax_zoom) = plt.subplots(1, 2, figsize=(9.5, 4.5),
                                            gridspec_kw={"width_ratios": [1.2, 1]})

    n_bands = eig.shape[0]
    for i in range(n_bands):
        ax_full.plot(k_vals, eig[i], color=BAND_COLORS[i], linewidth=1.2)
        ax_zoom.plot(k_vals, eig[i], color=BAND_COLORS[i], linewidth=1.2)

    # Full view
    ax_full.set_xlabel(r"$k$ (1/A)")
    ax_full.set_ylabel(r"$E$ (eV)")
    ax_full.set_title("Bulk GaAs 8-band dispersion")
    ax_full.axhline(0, color="grey", linewidth=0.4, linestyle="--")

    # Gamma-point energy annotations
    e0 = eig[:, 0]  # eigenvalues at k=0
    # Label distinct energies: SO, HH/LH, CB
    labels_placed = set()
    for i, e in enumerate(e0):
        rounded = round(e, 3)
        if rounded not in labels_placed:
            labels_placed.add(rounded)
            ax_full.annotate(f"{e:.3f} eV", xy=(k_vals[-1], e),
                           fontsize=7, color=BAND_COLORS[i], va="center")

    # Zoom near gap
    ax_zoom.set_xlim(-0.002, 0.06)
    ax_zoom.set_ylim(-0.5, 1.8)
    ax_zoom.set_xlabel(r"$k$ (1/A)")
    ax_zoom.set_title(r"Zoom near $\Gamma$")
    ax_zoom.axhline(0, color="grey", linewidth=0.4, linestyle="--")

    # Annotate band gap and SO splitting
    try:
        e_vb = max(e for e in e0 if e <= 0.01)
    except ValueError:
        e_vb = 0.0
    try:
        e_so = min(e for e in e0 if e <= 0.01)
    except ValueError:
        e_so = e_vb
    try:
        e_cb = min(e for e in e0 if e > 0.01)
    except ValueError:
        e_cb = max(e0) + 1.0

    gap = e_cb - e_vb
    so_split = e_vb - e_so

    ax_zoom.annotate(f"$E_g$ = {gap:.3f} eV", xy=(0.03, (e_vb + e_cb) / 2),
                    fontsize=8, color="grey", ha="left")
    ax_zoom.annotate(f"$\\Delta_{{SO}}$ = {so_split:.3f} eV", xy=(0.03, (e_so + e_vb) / 2),
                    fontsize=8, color="grey", ha="left")

    # Band labels on zoom
    ax_zoom.annotate("CB", xy=(0.04, e_cb + 0.03), fontsize=8, color="#17becf")
    ax_zoom.annotate("HH/LH", xy=(0.04, e_vb - 0.05), fontsize=8, color="#d62728")
    ax_zoom.annotate("SO", xy=(0.04, e_so - 0.05), fontsize=8, color="#ff7f0e")

    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "bulk_gaas_bands.png")
    plt.close(fig)
    print("  -> docs/figures/bulk_gaas_bands.png")


def fig_bulk_gaas_parts(output_dir: Path) -> None:
    """bulk_gaas_parts.png: band character at Gamma from bulk."""
    print("[figure] bulk_gaas_parts")
    # Always re-run bulk to ensure correct output (shared output dir gets overwritten)
    cfg = CONFIG_DIR / "bulk_gaas_kx.cfg"
    run_executable(EXE_BAND, cfg, REPO_ROOT, label="bulk_gaas_kx")
    parts = parse_parts(output_dir)

    if parts.size == 0:
        print("  WARNING: no parts data, skipping.")
        return

    n_ev = parts.shape[0]
    # Normalize: each row should sum to 1 (band character fraction)
    row_sums = parts.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0  # avoid division by zero
    parts = parts / row_sums

    fig, ax = plt.subplots(figsize=(6, 4))
    bar_width = 0.8
    x = np.arange(n_ev)
    bottom = np.zeros(n_ev)
    for b in range(min(8, parts.shape[1])):
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
    ax.set_ylim(0, 1.05)
    ax.legend(ncol=2, fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "bulk_gaas_parts.png")
    plt.close(fig)
    print("  -> docs/figures/bulk_gaas_parts.png")


def fig_bulk_gaas_parts_vs_k(output_dir: Path) -> None:
    """bulk_gaas_parts_vs_k.png: band character evolution along k."""
    print("[figure] bulk_gaas_parts_vs_k")
    cfg = CONFIG_DIR / "bulk_gaas_kxky.cfg"
    run_executable(EXE_BAND, cfg, REPO_ROOT, label="bulk_gaas_kxky")
    k_values, all_parts = parse_parts_all_k(output_dir)

    if len(all_parts) == 0:
        print("  WARNING: no parts data, skipping.")
        return

    n_k = len(k_values)
    n_ev = all_parts[0].shape[0]
    n_band = min(8, all_parts[0].shape[1])

    fig, axes = plt.subplots(2, 4, figsize=(12, 6), sharex=True)
    axes = axes.flatten()

    for ev in range(min(n_ev, 8)):
        ax = axes[ev]
        bottom = np.zeros(n_k)
        for b in range(n_band):
            vals = np.array([all_parts[k][ev, b] for k in range(n_k)])
            # Normalize
            row_sum = vals.sum()
            if row_sum > 0:
                vals = vals / np.sum([all_parts[k][ev, :].sum() for k in range(n_k)]) * n_k
                vals_raw = np.array([all_parts[k][ev, b] for k in range(n_k)])
                row_sums = np.array([all_parts[k][ev, :].sum() for k in range(n_k)])
                row_sums[row_sums == 0] = 1.0
                vals = vals_raw / row_sums
            ax.fill_between(k_values, bottom, bottom + vals,
                          color=BAND_COLORS[b], linewidth=0, alpha=0.9)
            bottom += vals
        ax.set_ylim(0, 1.05)
        ax.set_title(f"State {ev+1}", fontsize=9)
        ax.set_ylabel("Character", fontsize=7)
        if ev >= 4:
            ax.set_xlabel(r"$|k|$ (1/A)", fontsize=7)

    # Add legend to first panel
    for b in range(n_band):
        axes[0].plot([], [], color=BAND_COLORS[b], linewidth=6, label=BAND_NAMES[b])
    axes[0].legend(fontsize=5, ncol=2, loc="upper right")

    fig.suptitle("Bulk GaAs: band character evolution from $\\Gamma$", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "bulk_gaas_parts_vs_k.png")
    plt.close(fig)
    print("  -> docs/figures/bulk_gaas_parts_vs_k.png")


def fig_bulk_gaas_bands_110(output_dir: Path) -> None:
    """bulk_gaas_bands_110.png: E(k) along [110] (kxky)."""
    print("[figure] bulk_gaas_bands_110")
    cfg = CONFIG_DIR / "bulk_gaas_kxky.cfg"
    run_executable(EXE_BAND, cfg, REPO_ROOT, label="bulk_gaas_kxky")
    k_vals, eig = parse_eigenvalues(output_dir)

    n_bands = eig.shape[0]
    fig, ax = plt.subplots(figsize=(6, 5))
    for i in range(n_bands):
        ax.plot(k_vals, eig[i], color=BAND_COLORS[i], linewidth=1.2)

    ax.set_xlabel(r"$k_{[110]}$ (1/A)")
    ax.set_ylabel(r"$E$ (eV)")
    ax.set_title("Bulk GaAs $E(k)$ along $[110]$")
    ax.axhline(0, color="grey", linewidth=0.4, linestyle="--")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "bulk_gaas_bands_110.png")
    plt.close(fig)
    print("  -> docs/figures/bulk_gaas_bands_110.png")


def fig_bulk_gaas_warping(output_dir: Path) -> None:
    """bulk_gaas_warping.png: valence band E(k) comparing [100] vs [110]."""
    print("[figure] bulk_gaas_warping")
    # [100] direction
    cfg_100 = CONFIG_DIR / "bulk_gaas_kx.cfg"
    run_executable(EXE_BAND, cfg_100, REPO_ROOT, label="bulk_gaas_kx_100")
    k_100, eig_100 = parse_eigenvalues(output_dir)
    # Save before [110] run overwrites
    eig_100_saved = eig_100.copy()
    k_100_saved = k_100.copy()

    # [110] direction (kxky)
    cfg_110 = CONFIG_DIR / "bulk_gaas_kxky.cfg"
    run_executable(EXE_BAND, cfg_110, REPO_ROOT, label="bulk_gaas_kxky_110")
    k_110, eig_110 = parse_eigenvalues(output_dir)

    fig, (ax_full, ax_vb) = plt.subplots(1, 2, figsize=(10, 5))

    # Full bands
    for i in range(eig_100_saved.shape[0]):
        ax_full.plot(k_100_saved, eig_100_saved[i], color=BAND_COLORS[i],
                    linewidth=1.0, linestyle="-")
        ax_full.plot(k_110, eig_110[i], color=BAND_COLORS[i],
                    linewidth=1.0, linestyle="--")

    ax_full.set_xlabel(r"$|k|$ (1/A)")
    ax_full.set_ylabel(r"$E$ (eV)")
    ax_full.set_title("Bulk GaAs: [100] (solid) vs [110] (dashed)")
    ax_full.axhline(0, color="grey", linewidth=0.4, linestyle=":")

    # Zoom on valence bands to show warping
    for i in range(min(6, eig_100_saved.shape[0])):
        ax_vb.plot(k_100_saved, eig_100_saved[i], color=BAND_COLORS[i],
                  linewidth=1.2, linestyle="-", label="[100]" if i == 0 else "")
        ax_vb.plot(k_110, eig_110[i], color=BAND_COLORS[i],
                  linewidth=1.2, linestyle="--", label="[110]" if i == 0 else "")

    ax_vb.set_xlabel(r"$|k|$ (1/A)")
    ax_vb.set_ylabel(r"$E$ (eV)")
    ax_vb.set_title("Valence band warping (zoom)")
    ax_vb.set_ylim(-0.4, 0.05)
    ax_vb.axhline(0, color="grey", linewidth=0.4, linestyle=":")
    ax_vb.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "bulk_gaas_warping.png")
    plt.close(fig)
    print("  -> docs/figures/bulk_gaas_warping.png")


def fig_bulk_gaas_strained_bands(output_dir: Path) -> None:
    """bulk_gaas_strained_bands.png: strained GaAs on InP substrate."""
    print("[figure] bulk_gaas_strained_bands")
    cfg = CONFIG_DIR / "bulk_gaas_strained.cfg"
    run_executable(EXE_BAND, cfg, REPO_ROOT, label="bulk_gaas_strained")
    k_vals, eig = parse_eigenvalues(output_dir)

    n_bands = eig.shape[0]
    fig, ax = plt.subplots(figsize=(6, 5))
    for i in range(n_bands):
        ax.plot(k_vals, eig[i], color=BAND_COLORS[i], linewidth=1.2)

    ax.set_xlabel(r"$k_{[110]}$ (1/A)")
    ax.set_ylabel(r"$E$ (eV)")
    ax.set_title(r"GaAs strained on InP substrate ($a_{sub}=5.869$ A)")
    ax.axhline(0, color="grey", linewidth=0.4, linestyle="--")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "bulk_gaas_strained_bands.png")
    plt.close(fig)
    print("  -> docs/figures/bulk_gaas_strained_bands.png")


def fig_bulk_gaas_strain_comparison(output_dir: Path) -> None:
    """bulk_gaas_strain_comparison.png: unstrained vs strained GaAs."""
    print("[figure] bulk_gaas_strain_comparison")

    # Unstrained
    cfg_free = CONFIG_DIR / "bulk_gaas_kxky.cfg"
    run_executable(EXE_BAND, cfg_free, REPO_ROOT, label="bulk_gaas_unstrained")
    k_free, eig_free = parse_eigenvalues(output_dir)
    eig_free_saved = eig_free.copy()
    k_free_saved = k_free.copy()

    # Strained
    cfg_strain = CONFIG_DIR / "bulk_gaas_strained.cfg"
    run_executable(EXE_BAND, cfg_strain, REPO_ROOT, label="bulk_gaas_strained_cmp")
    k_strain, eig_strain = parse_eigenvalues(output_dir)

    fig, (ax_full, ax_vb) = plt.subplots(1, 2, figsize=(10, 5))

    for i in range(eig_free_saved.shape[0]):
        ax_full.plot(k_free_saved, eig_free_saved[i], color=BAND_COLORS[i],
                    linewidth=1.0, linestyle="--")
        ax_full.plot(k_strain, eig_strain[i], color=BAND_COLORS[i],
                    linewidth=1.2, linestyle="-")

    ax_full.set_xlabel(r"$k$ (1/A)")
    ax_full.set_ylabel(r"$E$ (eV)")
    ax_full.set_title("GaAs: unstrained (dashed) vs strained on InP (solid)")
    ax_full.axhline(0, color="grey", linewidth=0.4, linestyle=":")

    # Zoom on VB to show HH-LH splitting
    for i in range(min(6, eig_free_saved.shape[0])):
        ax_vb.plot(k_free_saved, eig_free_saved[i], color=BAND_COLORS[i],
                  linewidth=1.0, linestyle="--")
        ax_vb.plot(k_strain, eig_strain[i], color=BAND_COLORS[i],
                  linewidth=1.2, linestyle="-")

    ax_vb.set_xlabel(r"$k$ (1/A)")
    ax_vb.set_ylabel(r"$E$ (eV)")
    ax_vb.set_title("VB: HH/LH splitting from strain")
    ax_vb.set_ylim(-0.5, 0.05)
    ax_vb.axhline(0, color="grey", linewidth=0.4, linestyle=":")

    # Add legend
    ax_vb.plot([], [], "k--", linewidth=1.0, label="Unstrained")
    ax_vb.plot([], [], "k-", linewidth=1.2, label="Strained (InP)")
    ax_vb.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "bulk_gaas_strain_comparison.png")
    plt.close(fig)
    print("  -> docs/figures/bulk_gaas_strain_comparison.png")


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
    # Always re-run the QW to ensure fresh data
    cfg = CONFIG_DIR / "qw_alsbw_gasbw_inasw.cfg"
    run_executable(EXE_BAND, cfg, REPO_ROOT, label="qw_alsbw_gasbw_inasw")
    n_z = 101  # FDstep from config
    try:
        z, wf = parse_eigenfunctions_qw(output_dir, k_idx=1, n_ev=4, n_z=n_z)
    except FileNotFoundError:
        print("  WARNING: no eigenfunction data found, skipping.")
        return

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
        ax.set_xlabel(r"$z$ (\u00C5)")
        ax.set_title(f"State {idx + 1}")
        # Zoom to region where wavefunction is non-negligible
        nonzero = psi2_total > 0.01 * psi2_total.max()
        if np.any(nonzero):
            z_lo = z[nonzero][0] - 20
            z_hi = z[nonzero][-1] + 20
            ax.set_xlim(z_lo, z_hi)
    axes[0].set_ylabel(r"$|\psi(z)|^2$")
    fig.suptitle(r"QW probability density (AlSbW/GaSbW/InAsW, $k_{\parallel}=0$)", fontsize=12)
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

    # Normalize: each row should sum to 1 (band character fraction)
    row_sums = parts.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    parts = parts / row_sums

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
    ax.set_ylim(0, 1.05)
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
    bars_x = ax.bar(x_pos - width, [gx_vals[l] for l in labels], width, label=r"$g_x$", color="#1f77b4")
    bars_y = ax.bar(x_pos, [gy_vals[l] for l in labels], width, label=r"$g_y$", color="#ff7f0e")
    bars_z = ax.bar(x_pos + width, [gz_vals[l] for l in labels], width, label=r"$g_z$", color="#2ca02c")

    # Add value labels on bars — scale offset to data range
    all_heights = [abs(b.get_height()) for bars in [bars_x, bars_y, bars_z] for b in bars]
    max_h = max(all_heights) if all_heights else 1.0
    label_offset = 0.04 * max_h
    for bars in [bars_x, bars_y, bars_z]:
        for bar in bars:
            height = bar.get_height()
            if abs(height) > 0.01 * max_h:
                va = "bottom" if height > 0 else "top"
                offset = label_offset if height > 0 else -label_offset
                ax.text(bar.get_x() + bar.get_width() / 2, height + offset,
                       f"{height:.2f}", ha="center", va=va, fontsize=6)
            else:
                ax.text(bar.get_x() + bar.get_width() / 2, label_offset,
                       f"{height:.2f}", ha="center", va="bottom", fontsize=5,
                       color="grey")

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
    """sc_potential.png: band-edge profile with self-consistent band bending.

    Runs the SC GaAs/AlAs QW config.  If the SC loop produced
    ``sc_potential_profile.dat`` (with the converged band-bent edges) it
    is used; otherwise the flat-band ``potential_profile.dat`` is shown.
    """
    print("[figure] sc_potential")

    # Remove any stale sc_potential_profile.dat that may have been produced
    # by a prior wire SC run (5-column 2D format).  The QW run below will
    # produce a fresh potential_profile.dat with 4-column 1D format.
    sc_path = output_dir / "sc_potential_profile.dat"
    if sc_path.exists():
        sc_path.unlink()
        print("  Removed stale sc_potential_profile.dat")

    # Run SC config to get band bending from self-consistent calculation
    cfg = CONFIG_DIR / "sc_gaas_alas_qw.cfg"
    result = run_executable(EXE_BAND, cfg, REPO_ROOT, label="sc_gaas_alas_qw", timeout=600)
    if result.returncode != 0:
        print("  WARNING: SC run failed, skipping sc_potential figure.")
        return

    # Try the SC-converged profile first; fall back to the flat-band profile
    try:
        z, EV, EV_SO, EC = parse_sc_potential(output_dir)
        title = "SC band-edge profile (GaAs/AlAs QW)"
    except FileNotFoundError:
        try:
            z, EV, EV_SO, EC = parse_potential_profile(output_dir)
            title = "Band-edge profile (GaAs/AlAs QW, flat-band)"
        except FileNotFoundError:
            print("  WARNING: no potential profile data found, skipping.")
            return

    # Read flat-band profile for comparison
    try:
        _, EV_fb, _, EC_fb = parse_potential_profile(output_dir)
        has_flatband = True
    except FileNotFoundError:
        has_flatband = False

    fig, (ax_full, ax_well) = plt.subplots(1, 2, figsize=(9.5, 4.5),
                                            gridspec_kw={"width_ratios": [1.2, 1]})

    # Full profile
    ax_full.plot(z, EV, color="#d62728", linewidth=1.5, label="$E_V$")
    ax_full.plot(z, EV_SO, color="#ff7f0e", linewidth=1.5, label="$E_{\\Delta SO}$")
    ax_full.plot(z, EC, color="#17becf", linewidth=1.5, label="$E_C$")
    ax_full.fill_between(z, EV.min() - 0.1, EV, alpha=0.06, color="#d62728")
    ax_full.fill_between(z, EC, EC.max() + 0.1, alpha=0.06, color="#17becf")
    ax_full.set_xlabel(r"$z$ (\u00C5)")
    ax_full.set_ylabel("Energy (eV)")
    ax_full.set_title("GaAs/AlAs QW SC band edges")
    ax_full.legend(loc="best", fontsize=8)

    # Zoom on well region — show SC bending vs flat-band
    well_mask = (z >= -60) & (z <= 60)
    z_well = z[well_mask]
    ax_well.plot(z_well, EC[well_mask], color="#17becf", linewidth=1.5, label="$E_C$ (SC)")
    ax_well.plot(z_well, EV[well_mask], color="#d62728", linewidth=1.5, label="$E_V$ (SC)")
    if has_flatband:
        ax_well.plot(z_well, EC_fb[well_mask], "--", color="#17becf", linewidth=0.8,
                    alpha=0.6, label="$E_C$ (flat)")
        ax_well.plot(z_well, EV_fb[well_mask], "--", color="#d62728", linewidth=0.8,
                    alpha=0.6, label="$E_V$ (flat)")
    # Annotate band bending magnitude
    well_inner = (z >= -5) & (z <= 5)
    if np.any(well_inner):
        ec_sc_mid = np.mean(EC[well_inner])
        if has_flatband:
            ec_fb_mid = np.mean(EC_fb[well_inner])
            bending = ec_sc_mid - ec_fb_mid
            ax_well.annotate(f"SC shift: {bending*1000:.1f} meV",
                           xy=(0, ec_sc_mid), fontsize=8, color="grey", ha="center",
                           va="bottom" if bending > 0 else "top")
    ax_well.axvline(-50, color="grey", linewidth=0.5, linestyle=":")
    ax_well.axvline(50, color="grey", linewidth=0.5, linestyle=":")
    ax_well.set_xlabel(r"$z$ (\u00C5)")
    ax_well.set_ylabel("Energy (eV)")
    ax_well.set_title("Zoom: well band bending")
    ax_well.legend(loc="best", fontsize=7)

    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "sc_potential.png")
    plt.close(fig)
    print("  -> docs/figures/sc_potential.png")


def fig_sc_charge_density(output_dir: Path) -> None:
    """sc_charge_density.png: n(z) charge density from SC run."""
    print("[figure] sc_charge_density")
    # Clean stale 2D-format sc_charge.dat from wire runs
    sc_charge = output_dir / "sc_charge.dat"
    if sc_charge.exists():
        sc_charge.unlink()
    sc_path = output_dir / "sc_potential_profile.dat"
    if sc_path.exists():
        sc_path.unlink()

    # Need fresh SC run (shared output dir gets overwritten by other configs)
    cfg = CONFIG_DIR / "sc_gaas_alas_qw.cfg"
    result = run_executable(EXE_BAND, cfg, REPO_ROOT, label="sc_charge", timeout=600)
    if result.returncode != 0:
        print("  WARNING: SC run failed, skipping sc_charge_density figure.")
        return

    try:
        z, n_e, n_h = parse_sc_charge(output_dir, is_2d=False)
    except FileNotFoundError:
        print("  WARNING: sc_charge.dat not found after SC run, skipping.")
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
    """wire_subbands.png: E(k_z) dispersion for GaAs rectangular wire near the gap."""
    print("[figure] wire_subbands")
    cfg = CONFIG_DIR / "wire_gaas_rectangle.cfg"
    run_executable(EXE_BAND, cfg, REPO_ROOT, label="wire_gaas", timeout=300)
    try:
        k_vals, eig = parse_eigenvalues(output_dir)
    except FileNotFoundError:
        print("  WARNING: eigenvalues.dat not found for wire, skipping.")
        return

    n_bands = eig.shape[0]
    e0 = eig[:, 0]  # eigenvalues at k=0, sorted ascending

    # Classify VB/CB by finding the largest gap in the eigenvalue spectrum.
    # The band gap should be the largest gap between consecutive eigenvalues
    # near the expected gap region (between -0.5 and 2.0 eV for GaAs).
    gaps = np.diff(e0)
    # Only consider gaps in the relevant energy range
    gap_mask = (e0[:-1] > -1.0) & (e0[:-1] < 2.0)
    if np.any(gap_mask):
        gap_idx = np.argmax(gaps * gap_mask)  # largest gap in range
        n_vb = gap_idx + 1
        n_cb = n_bands - n_vb
    else:
        # Fallback: use config values
        n_vb = min(16, n_bands)
        n_cb = n_bands - n_vb

    e_vb_max = e0[n_vb - 1] if n_vb > 0 else e0[0]
    e_cb_min = e0[n_vb] if n_cb > 0 else e0[-1]

    # Only plot a subset of bands near the gap (up to 10 VB + 10 CB)
    max_plot_vb = min(n_vb, 10)
    max_plot_cb = min(n_cb, 10)
    plot_vb_start = n_vb - max_plot_vb  # highest VB states

    print(f"  Classification: {n_vb} VB + {n_cb} CB ({n_bands} total)")
    print(f"  VB max={e_vb_max:.4f}, CB min={e_cb_min:.4f}, "
          f"gap={e_cb_min - e_vb_max:.4f} eV")
    print(f"  Plotting {max_plot_vb} VB + {max_plot_cb} CB near the gap")

    fig, (ax_vb, ax_cb) = plt.subplots(1, 2, figsize=(8, 4.5), sharey=False)

    # VB subbands near the gap
    if max_plot_vb > 0:
        vb_colors = plt.cm.Reds_r(np.linspace(0.3, 0.8, max_plot_vb))
        for i, band_i in enumerate(range(plot_vb_start, n_vb)):
            ax_vb.plot(k_vals, eig[band_i], color=vb_colors[i], linewidth=0.8)

    ax_vb.axhline(e_vb_max, color="#d62728", linewidth=0.8, linestyle=":",
               label=f"VB top = {e_vb_max:.3f} eV")
    ax_vb.set_xlabel(r"$k_z$ (1/\u00C5)")
    ax_vb.set_ylabel(r"$E$ (eV)")
    ax_vb.set_title(f"VB subbands ({max_plot_vb} shown)")
    ax_vb.legend(loc="best", fontsize=7)

    # CB subbands near the gap
    if max_plot_cb > 0:
        cb_colors = plt.cm.Blues_r(np.linspace(0.3, 0.8, max_plot_cb))
        for i in range(max_plot_cb):
            band_i = n_vb + i
            ax_cb.plot(k_vals, eig[band_i], color=cb_colors[i], linewidth=0.8)

    ax_cb.axhline(e_cb_min, color="#17becf", linewidth=0.8, linestyle="--",
               label=f"CB bottom = {e_cb_min:.3f} eV")
    ax_cb.set_xlabel(r"$k_z$ (1/\u00C5)")
    ax_cb.set_ylabel(r"$E$ (eV)")
    ax_cb.set_title(f"CB subbands ({max_plot_cb} shown)")
    ax_cb.legend(loc="best", fontsize=7)

    wire_w = float(
        next(
            (l.split(":")[1].strip() for l in open(cfg) if l.strip().startswith("wire_width")),
            "63.0",
        )
    )
    wire_h = float(
        next(
            (l.split(":")[1].strip() for l in open(cfg) if l.strip().startswith("wire_height")),
            "63.0",
        )
    )
    fig.suptitle(f"GaAs wire ({wire_w:.0f}×{wire_h:.0f} Å, gap={e_cb_min-e_vb_max:.3f} eV)",
                fontsize=10)
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
    """wire_density_2d.png: 2D |psi(x,y)|^2 for CB-ground and VB-edge states."""
    print("[figure] wire_density_2d")

    # Read eigenvalues to find VB max and CB min
    try:
        _, eig = parse_eigenvalues(output_dir)
    except FileNotFoundError:
        cfg = CONFIG_DIR / "wire_gaas_rectangle.cfg"
        run_executable(EXE_BAND, cfg, REPO_ROOT, label="wire_gaas", timeout=300)
        _, eig = parse_eigenvalues(output_dir)

    e0 = eig[:, 0]
    n_bands = eig.shape[0]

    # Classify VB/CB by finding the largest gap in the eigenvalue spectrum
    # (same approach as fig_wire_subbands). The dense LAPACK solver in range
    # mode returns ALL eigenvalues in [emin, emax], not just numcb+numvb.
    gaps = np.diff(e0)
    gap_mask = (e0[:-1] > -1.0) & (e0[:-1] < 2.0)
    if np.any(gap_mask):
        gap_idx = np.argmax(gaps * gap_mask)
        n_vb = gap_idx + 1
        n_cb = n_bands - n_vb
    else:
        n_vb = min(16, n_bands)
        n_cb = n_bands - n_vb

    e_vb_max = e0[n_vb - 1] if n_vb > 0 else e0[0]
    e_cb_min = e0[n_vb] if n_cb > 0 else e0[-1]

    # Plot two panels: VB-edge state (highest VB) and CB-ground state (lowest CB)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

    for panel, (target_eV, label) in enumerate([
        (e_vb_max, r"VB-edge $|\psi|^2$"),
        (e_cb_min, r"CB-ground $|\psi|^2$")
    ]):
        ax = axes[panel]
        ev_idx = _find_ev_idx_nearest_to_energy(output_dir, target_eV)
        print(f"  Panel {panel+1}: ev_idx={ev_idx} (target E={target_eV:.4f} eV)")

        try:
            x, y, psi2 = parse_wire_eigenfunction_2d(output_dir, k_idx=1, ev_idx=ev_idx)
        except FileNotFoundError:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
            continue

        if x.size == 0:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
            continue

        state_energy = ""
        ev_file = output_dir / f"eigenfunctions_k_00001_ev_{ev_idx:05d}.dat"
        if ev_file.exists():
            with open(ev_file) as f:
                for line in f:
                    if line.startswith("# E ="):
                        state_energy = line.replace("# E =", "").strip()
                        break

        im = ax.pcolormesh(x, y, psi2, shading="auto", cmap="viridis")
        ax.set_xlabel(r"$x$ (\u00C5)")
        ax.set_ylabel(r"$y$ (\u00C5)")
        title = label
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


def fig_bulk_inas_bands(output_dir: Path) -> None:
    """bulk_inas_bands.png: 8-band E(k) dispersion from bulk InAs."""
    print("[figure] bulk_inas_bands")
    cfg = CONFIG_DIR / "bulk_inas_kx.cfg"
    result = run_executable(EXE_BAND, cfg, REPO_ROOT, label="bulk_inas_kx")
    if result.returncode != 0:
        print("  WARNING: bulk_inas_kx run failed, skipping.")
        return
    try:
        k_vals, eig = parse_eigenvalues(output_dir)
    except FileNotFoundError:
        print("  WARNING: eigenvalues.dat not found for bulk InAs, skipping.")
        return

    fig, (ax_full, ax_zoom) = plt.subplots(1, 2, figsize=(9.5, 4.5),
                                            gridspec_kw={"width_ratios": [1.2, 1]})

    n_bands = eig.shape[0]
    for i in range(n_bands):
        ax_full.plot(k_vals, eig[i], color=BAND_COLORS[i], linewidth=1.2)
        ax_zoom.plot(k_vals, eig[i], color=BAND_COLORS[i], linewidth=1.2)

    # Full view
    ax_full.set_xlabel(r"$k$ (1/A)")
    ax_full.set_ylabel(r"$E$ (eV)")
    ax_full.set_title("Bulk InAs 8-band dispersion")
    ax_full.axhline(0, color="grey", linewidth=0.4, linestyle="--")

    # Gamma-point energy annotations
    e0 = eig[:, 0]
    labels_placed = set()
    for i, e in enumerate(e0):
        rounded = round(e, 3)
        if rounded not in labels_placed:
            labels_placed.add(rounded)
            ax_full.annotate(f"{e:.3f} eV", xy=(k_vals[-1], e),
                           fontsize=7, color=BAND_COLORS[i], va="center")

    # Zoom near gap
    ax_zoom.set_xlim(-0.002, 0.06)
    ax_zoom.set_ylim(-1.0, 1.2)
    ax_zoom.set_xlabel(r"$k$ (1/A)")
    ax_zoom.set_title(r"Zoom near $\Gamma$")
    ax_zoom.axhline(0, color="grey", linewidth=0.4, linestyle="--")

    # Annotate band gap and SO splitting
    try:
        e_vb = max(e for e in e0 if e <= 0.01)
    except ValueError:
        e_vb = 0.0
    try:
        e_so = min(e for e in e0 if e <= 0.01)
    except ValueError:
        e_so = e_vb
    try:
        e_cb = min(e for e in e0 if e > 0.01)
    except ValueError:
        e_cb = max(e0) + 1.0

    gap = e_cb - e_vb
    so_split = e_vb - e_so

    ax_zoom.annotate(f"$E_g$ = {gap:.3f} eV", xy=(0.03, (e_vb + e_cb) / 2),
                    fontsize=8, color="grey", ha="left")
    ax_zoom.annotate(f"$\\Delta_{{SO}}$ = {so_split:.3f} eV", xy=(0.03, (e_so + e_vb) / 2),
                    fontsize=8, color="grey", ha="left")

    # Band labels on zoom
    ax_zoom.annotate("CB", xy=(0.04, e_cb + 0.03), fontsize=8, color="#17becf")
    ax_zoom.annotate("HH/LH", xy=(0.04, e_vb - 0.05), fontsize=8, color="#d62728")
    ax_zoom.annotate("SO", xy=(0.04, e_so - 0.05), fontsize=8, color="#ff7f0e")

    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "bulk_inas_bands.png")
    plt.close(fig)
    print("  -> docs/figures/bulk_inas_bands.png")


def fig_qcse_stark_shift(output_dir: Path) -> None:
    """qcse_stark_shift.png: QW subbands at 0 field vs -70 kV/cm — side by side."""
    print("[figure] qcse_stark_shift")
    cfg_zero = CONFIG_DIR / "sc_qcse_gaas_algaas.cfg"
    cfg_field = CONFIG_DIR / "sc_qcse_gaas_algaas_ef.cfg"

    # Run zero-field config
    result = run_executable(EXE_BAND, cfg_zero, REPO_ROOT, label="qcse_zero_field",
                           timeout=600)
    if result.returncode != 0:
        print("  WARNING: zero-field QCSE run failed, skipping.")
        return
    try:
        z0, EV0, _, EC0 = parse_potential_profile(output_dir)
    except FileNotFoundError:
        print("  WARNING: potential_profile.dat not found for zero field, skipping.")
        return
    # Also get eigenvalues at k=0
    try:
        _, eig0 = parse_eigenvalues(output_dir)
        e0_vals = eig0[:, 0]
    except FileNotFoundError:
        e0_vals = None

    # Run with-field config
    result = run_executable(EXE_BAND, cfg_field, REPO_ROOT, label="qcse_efield",
                           timeout=600)
    if result.returncode != 0:
        print("  WARNING: electric-field QCSE run failed, skipping.")
        return
    try:
        z_ef, EV_ef, _, EC_ef = parse_potential_profile(output_dir)
    except FileNotFoundError:
        print("  WARNING: potential_profile.dat not found for e-field, skipping.")
        return
    try:
        _, eig_ef = parse_eigenvalues(output_dir)
        ef_vals = eig_ef[:, 0]
    except FileNotFoundError:
        ef_vals = None

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(10, 5))

    # Left panel: band-edge profiles
    ax_left.plot(z0, EC0, color="#17becf", linewidth=1.5, label=r"$E_C$ (0 field)")
    ax_left.plot(z0, EV0, color="#d62728", linewidth=1.5, label=r"$E_V$ (0 field)")
    ax_left.plot(z_ef, EC_ef, "--", color="#17becf", linewidth=1.3,
                alpha=0.8, label=r"$E_C$ (-70 kV/cm)")
    ax_left.plot(z_ef, EV_ef, "--", color="#d62728", linewidth=1.3,
                alpha=0.8, label=r"$E_V$ (-70 kV/cm)")
    # Mark well region
    ax_left.axvline(-30, color="grey", linewidth=0.5, linestyle=":")
    ax_left.axvline(30, color="grey", linewidth=0.5, linestyle=":")
    ax_left.set_xlabel(r"$z$ (\u00C5)")
    ax_left.set_ylabel("Energy (eV)")
    ax_left.set_title("QCSE: Band-edge tilt under field")
    ax_left.legend(loc="best", fontsize=7)

    # Right panel: energy level diagram with horizontal lines
    if e0_vals is not None and ef_vals is not None:
        n_show = min(len(e0_vals), len(ef_vals), 8)
        col_zero = 1
        col_field = 2
        line_len = 0.3
        for i in range(n_show):
            # 0-field level
            ax_right.hlines(e0_vals[i], col_zero - line_len / 2,
                           col_zero + line_len / 2,
                           colors="#1f77b4", linewidth=2.0)
            # With-field level
            ax_right.hlines(ef_vals[i], col_field - line_len / 2,
                           col_field + line_len / 2,
                           colors="#ff7f0e", linewidth=2.0)
            # Annotate Stark shift (delta_E) between the pair
            delta_e = ef_vals[i] - e0_vals[i]
            if abs(delta_e) > 1e-6:
                mid_e = (e0_vals[i] + ef_vals[i]) / 2
                ax_right.annotate(f"{delta_e*1000:+.1f}",
                                xy=((col_zero + col_field) / 2, mid_e),
                                fontsize=6, color="grey", ha="center",
                                va="center")
                # Dotted connector between the two levels
                ax_right.plot([col_zero + line_len / 2, col_field - line_len / 2],
                             [e0_vals[i], ef_vals[i]],
                             ":", color="grey", linewidth=0.6)

        # Column labels
        ax_right.set_xticks([col_zero, col_field])
        ax_right.set_xticklabels(["0 field", "-70 kV/cm"])
        ax_right.set_ylabel("Energy (eV)")
        ax_right.set_title("Subband Stark shift (meV)")

        # Manual legend via proxy artists
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], color="#1f77b4", linewidth=2, label="0 field"),
            Line2D([0], [0], color="#ff7f0e", linewidth=2, label="-70 kV/cm"),
        ]
        ax_right.legend(handles=legend_elements, loc="best", fontsize=8)
    else:
        ax_right.text(0.5, 0.5, "No eigenvalue data", transform=ax_right.transAxes,
                     ha="center", va="center")

    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "qcse_stark_shift.png")
    plt.close(fig)
    print("  -> docs/figures/qcse_stark_shift.png")


def fig_qw_gaas_algaas_subbands(output_dir: Path) -> None:
    """qw_gaas_algaas_subbands.png: Type-I GaAs/AlGaAs QW E(k_parallel) subbands."""
    print("[figure] qw_gaas_algaas_subbands")
    cfg = REPO_ROOT / "docs" / "benchmarks" / "qw_gaas_algaas.cfg"
    result = run_executable(EXE_BAND, cfg, REPO_ROOT, label="qw_gaas_algaas")
    if result.returncode != 0:
        print("  WARNING: qw_gaas_algaas run failed, skipping.")
        return
    try:
        k_vals, eig = parse_eigenvalues(output_dir)
    except FileNotFoundError:
        print("  WARNING: eigenvalues.dat not found for GaAs/AlGaAs QW, skipping.")
        return

    fig, ax = plt.subplots(figsize=(6, 5))
    n_bands = eig.shape[0]
    for i in range(n_bands):
        energy_mid = np.mean(eig[i])
        if energy_mid > 0.5:
            color = "#17becf"
            alpha = 0.85
        else:
            color = "#d62728"
            alpha = 0.85
        ax.plot(k_vals, eig[i], color=color, linewidth=0.9, alpha=alpha)

    ax.set_xlabel(r"$k_{\parallel}$ (1/A)")
    ax.set_ylabel(r"$E$ (eV)")
    ax.set_title("GaAs/Al$_{0.3}$Ga$_{0.7}$As QW subbands")
    ax.axhline(0, color="grey", linewidth=0.4, linestyle="--")

    # Annotate VB top and CB bottom
    e0 = eig[:, 0]
    try:
        e_vb = max(e for e in e0 if e <= 0.01)
    except ValueError:
        e_vb = 0.0
    try:
        e_cb = min(e for e in e0 if e > 0.01)
    except ValueError:
        e_cb = max(e0) + 1.0
    gap = e_cb - e_vb
    ax.annotate(f"Gap = {gap*1000:.1f} meV", xy=(k_vals[-1] * 0.6, (e_vb + e_cb) / 2),
               fontsize=8, color="grey", ha="center")

    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "qw_gaas_algaas_subbands.png")
    plt.close(fig)
    print("  -> docs/figures/qw_gaas_algaas_subbands.png")


def fig_convergence_grid_spacing(output_dir: Path) -> None:
    """convergence_grid_spacing.png: eigenvalue error vs dz at fixed FD order.

    Runs the same QW config with varying FDstep (51, 101, 201, 401, 801)
    and plots error against the finest grid (801) as reference.
    """
    print("[figure] convergence_grid_spacing")
    base_cfg = CONFIG_DIR / "qw_alsbw_gasbw_inasw.cfg"
    if not base_cfg.exists():
        print("  WARNING: qw_alsbw_gasbw_inasw.cfg not found, skipping.")
        return

    base_text = base_cfg.read_text()

    fd_steps = [51, 101, 201, 401, 801]
    energies_at_gamma: Dict[int, np.ndarray] = {}
    dz_values: Dict[int, float] = {}

    # Total z-range from config: -150 to 150 = 300 A
    # dz = 300 / (FDstep - 1)
    total_range = 300.0

    for step in fd_steps:
        lines = base_text.splitlines()
        modified_lines = []
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("FDstep"):
                modified_lines.append(f"FDstep: {step}")
            elif stripped.startswith("waveVectorMax"):
                modified_lines.append("waveVectorMax: 0.0")
            elif stripped.startswith("waveVectorStep"):
                modified_lines.append("waveVectorStep: 1")
            elif stripped.startswith("waveVector:"):
                modified_lines.append("waveVector: k0")
            elif stripped.startswith("numcb"):
                modified_lines.append("numcb: 4")
            elif stripped.startswith("numvb"):
                modified_lines.append("numvb: 4")
            else:
                modified_lines.append(line)
        modified_text = "\n".join(modified_lines) + "\n"

        tmp_cfg = REPO_ROOT / "input.cfg"
        tmp_cfg.write_text(modified_text)

        print(f"  [FDstep={step}] Running QW ...")
        result = subprocess.run(
            [str(EXE_BAND)],
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
            timeout=300,
        )
        if result.returncode != 0:
            print(f"  [FDstep={step}] FAILED, skipping.")
            continue

        try:
            _, eig = parse_eigenvalues(output_dir)
            energies_at_gamma[step] = eig[:, 0]
            dz_values[step] = total_range / (step - 1)
        except FileNotFoundError:
            print(f"  [FDstep={step}] eigenvalues.dat not found, skipping.")
            continue

    if len(energies_at_gamma) < 2:
        print("  WARNING: not enough grid-spacing results for convergence plot, skipping.")
        return

    fig, ax = plt.subplots(figsize=(6, 4.5))
    steps_sorted = sorted(energies_at_gamma.keys())
    ref_step = max(steps_sorted)
    ref_e = energies_at_gamma[ref_step]

    for band_idx in range(min(ref_e.shape[0], 8)):
        diffs = []
        dz_plot = []
        for step in steps_sorted:
            if step == ref_step:
                continue
            e = energies_at_gamma[step]
            if band_idx < e.shape[0]:
                diffs.append(abs(e[band_idx] - ref_e[band_idx]))
                dz_plot.append(dz_values[step])
        if len(dz_plot) >= 2:
            ax.loglog(dz_plot, diffs, "o-", linewidth=1.2, markersize=5,
                     label=f"Band {band_idx + 1}", color=BAND_COLORS[band_idx])

    # Reference slope line (2nd-order convergence: error ~ dz^2)
    # Collect dz values cleanly from non-reference steps
    all_dz = sorted(dz_values[s] for s in steps_sorted if s != ref_step)
    if len(all_dz) >= 2:
        dz_arr = np.array(all_dz)
        # Scale the reference line to pass through the mid-point of first band data
        ref_diffs = []
        ref_dzs = []
        for step in steps_sorted:
            if step == ref_step:
                continue
            e = energies_at_gamma[step]
            if 0 < e.shape[0]:
                ref_diffs.append(abs(e[0] - ref_e[0]))
                ref_dzs.append(dz_values[step])
        if ref_diffs:
            mid_idx = len(ref_dzs) // 2
            scale = ref_diffs[mid_idx] / ref_dzs[mid_idx] ** 2
            ax.loglog(dz_arr, scale * dz_arr ** 2, "k--", linewidth=0.8, alpha=0.5,
                     label=r"$\propto \Delta z^2$")

    ax.set_xlabel(r"Grid spacing $\Delta z$ (\u00C5)")
    ax.set_ylabel(r"$|E_n - E_n^{\rm ref}|$ (eV)")
    ax.set_title(f"Grid convergence (ref: FDstep={ref_step})")
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "convergence_grid_spacing.png")
    plt.close(fig)
    print("  -> docs/figures/convergence_grid_spacing.png")


def fig_benchmark_inasw_gasbw_broken_gap(output_dir: Path) -> None:
    """benchmark_inasw_gasbw_broken_gap.png: broken-gap QW E(k) subbands."""
    print("[figure] benchmark_inasw_gasbw_broken_gap")
    cfg = CONFIG_DIR / "qw_inasw_gasbw_broken_gap.cfg"
    result = run_executable(EXE_BAND, cfg, REPO_ROOT, label="qw_inasw_gasbw_broken_gap",
                  timeout=300)
    if result.returncode != 0:
        print("  WARNING: broken-gap run failed, skipping.")
        return
    try:
        k_vals, eig = parse_eigenvalues(output_dir)
    except FileNotFoundError:
        print("  WARNING: eigenvalues.dat not found for broken-gap QW, skipping.")
        return

    fig, ax = plt.subplots(figsize=(6, 5))
    n_bands = eig.shape[0]

    # For broken-gap systems, VB and CB overlap — use energy-based coloring
    for i in range(n_bands):
        energy_mid = np.mean(eig[i])
        if energy_mid > 0.2:
            color = "#17becf"  # CB-like
        elif energy_mid > -0.2:
            color = "#2ca02c"  # near-gap / mixed
        else:
            color="#d62728"  # VB-like
        ax.plot(k_vals, eig[i], color=color, linewidth=0.9, alpha=0.85)

    ax.set_xlabel(r"$k_{\parallel}$ (1/A)")
    ax.set_ylabel(r"$E$ (eV)")
    ax.set_title("Broken-gap InAsW/GaSbW QW subbands")
    ax.axhline(0, color="grey", linewidth=0.4, linestyle="--")

    # Annotate the overlap region
    e0 = eig[:, 0]
    e_cb_min = min(e for e in e0 if e > 0.01) if any(e > 0.01 for e in e0) else 0
    e_vb_max = max(e for e in e0 if e <= 0.01) if any(e <= 0.01 for e in e0) else 0
    if e_vb_max > e_cb_min:
        overlap = e_vb_max - e_cb_min
        ax.annotate(f"Overlap: {overlap*1000:.1f} meV",
                   xy=(k_vals[len(k_vals)//2], (e_vb_max + e_cb_min) / 2),
                   fontsize=8, color="grey", ha="center")
    else:
        gap = e_cb_min - e_vb_max
        ax.annotate(f"Gap: {gap*1000:.1f} meV",
                   xy=(k_vals[len(k_vals)//2], (e_vb_max + e_cb_min) / 2),
                   fontsize=8, color="grey", ha="center")

    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "benchmark_inasw_gasbw_broken_gap.png")
    plt.close(fig)
    print("  -> docs/figures/benchmark_inasw_gasbw_broken_gap.png")


# ---------------------------------------------------------------------------
# New QW dispersion, optics, and potential-profile figures
# ---------------------------------------------------------------------------


def fig_qw_dispersion_gaas_algaas(output_dir: Path) -> None:
    """qw_dispersion_gaas_algaas.png: GaAs/AlGaAs QW E(k_parallel) with band-type coloring."""
    print("[figure] qw_dispersion_gaas_algaas")
    cfg = CONFIG_DIR / "qw_gaas_algaas_kpar.cfg"
    result = run_executable(EXE_BAND, cfg, REPO_ROOT, label="qw_gaas_algaas_kpar",
                           timeout=600)
    if result.returncode != 0:
        print("  WARNING: qw_gaas_algaas_kpar run failed, skipping.")
        return
    try:
        k_vals, eig = parse_eigenvalues(output_dir)
    except FileNotFoundError:
        print("  WARNING: eigenvalues.dat not found, skipping.")
        return

    fig, ax = plt.subplots(figsize=(6, 5))
    n_bands = eig.shape[0]
    for i in range(n_bands):
        energy_mid = np.mean(eig[i])
        if energy_mid > 0.5:
            color = "#17becf"  # CB - cyan
        elif energy_mid <= 0.5 and energy_mid > 0:
            color = "#d62728"  # VB - red
        else:
            color = "#ff7f0e"  # SO - orange
        ax.plot(k_vals, eig[i], color=color, linewidth=0.9, alpha=0.85)

    # HH/LH labels near VB top at k=0
    e0 = eig[:, 0]
    vb_top = [e for e in e0 if -0.5 < e <= 0.5]
    if vb_top:
        e_vb_max = max(vb_top)
        e_vb_min = min(vb_top)
        ax.annotate("HH", xy=(k_vals[1], e_vb_max), fontsize=8, color="#d62728",
                    va="bottom")
        if abs(e_vb_max - e_vb_min) > 0.01:
            ax.annotate("LH", xy=(k_vals[1], e_vb_min), fontsize=8, color="#d62728",
                        va="top")

    ax.set_xlabel(r"$k_{\parallel}$ (1/A)")
    ax.set_ylabel(r"$E$ (eV)")
    ax.set_title(r"GaAs/Al$_{0.3}$Ga$_{0.7}$As QW Dispersion")
    ax.axhline(0, color="grey", linewidth=0.4, linestyle="--")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "qw_dispersion_gaas_algaas.png")
    plt.close(fig)
    print("  -> docs/figures/qw_dispersion_gaas_algaas.png")


def fig_qw_dispersion_broken_gap(output_dir: Path) -> None:
    """qw_dispersion_broken_gap.png: InAsW/GaSbW broken-gap QW with anticrossing."""
    print("[figure] qw_dispersion_broken_gap")
    cfg = CONFIG_DIR / "qw_inas_gasb_broken_gap_kpar.cfg"
    result = run_executable(EXE_BAND, cfg, REPO_ROOT,
                           label="qw_inas_gasb_broken_gap_kpar", timeout=600)
    if result.returncode != 0:
        print("  WARNING: broken-gap kpar run failed, skipping.")
        return
    try:
        k_vals, eig = parse_eigenvalues(output_dir)
    except FileNotFoundError:
        print("  WARNING: eigenvalues.dat not found, skipping.")
        return

    fig, ax = plt.subplots(figsize=(6, 5))
    n_bands = eig.shape[0]
    for i in range(n_bands):
        energy_mid = np.mean(eig[i])
        if energy_mid > 0.1:
            color = "#17becf"  # CB-like - cyan
        elif energy_mid >= -0.1:
            color = "#2ca02c"  # mixed / near-gap - green
        else:
            color = "#d62728"  # VB-like - red
        ax.plot(k_vals, eig[i], color=color, linewidth=0.9, alpha=0.85)

    # Find anticrossing: minimum gap between any two bands across k-points
    min_gap_overall = np.inf
    anticrossing_k = k_vals[0]
    n_k = len(k_vals)
    for ki in range(n_k):
        energies_at_k = sorted(eig[:, ki])
        for j in range(len(energies_at_k) - 1):
            gap_here = energies_at_k[j + 1] - energies_at_k[j]
            if gap_here < min_gap_overall:
                min_gap_overall = gap_here
                anticrossing_k = k_vals[ki]

    ax.axvline(anticrossing_k, color="grey", linewidth=0.8, linestyle="--", alpha=0.7)
    ax.annotate(f"Anticrossing at k ~ {anticrossing_k:.2f} 1/A",
                xy=(anticrossing_k, ax.get_ylim()[1] * 0.85 if ax.get_ylim()[1] > 0 else 0.5),
                fontsize=8, color="grey", ha="left",
                xytext=(anticrossing_k + 0.01, 0.0))

    ax.set_xlabel(r"$k_{\parallel}$ (1/A)")
    ax.set_ylabel(r"$E$ (eV)")
    ax.set_title("InAsW/GaSbW Broken-Gap QW Dispersion")
    ax.axhline(0, color="grey", linewidth=0.4, linestyle="--")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "qw_dispersion_broken_gap.png")
    plt.close(fig)
    print("  -> docs/figures/qw_dispersion_broken_gap.png")


def fig_qw_optical_matrix_elements(output_dir: Path) -> None:
    """qw_optical_matrix_elements.png: grouped bar chart of optical matrix elements."""
    print("[figure] qw_optical_matrix_elements")
    cfg = CONFIG_DIR / "qw_gaas_algaas_optics.cfg"
    result = run_executable(EXE_GFACTOR, cfg, REPO_ROOT,
                           label="qw_gaas_algaas_optics", timeout=1200)
    if result.returncode != 0:
        print("  WARNING: qw_gaas_algaas_optics run failed, skipping.")
        return
    try:
        cb_idx, vb_idx, energy, px, py, pz, f_osc = parse_optical_transitions(output_dir)
    except FileNotFoundError:
        print("  WARNING: optical_transitions.dat not found, skipping.")
        return

    # Filter to positive-energy transitions only
    pos_mask = energy > 0
    if not np.any(pos_mask):
        print("  WARNING: no positive-energy transitions found, skipping.")
        return
    cb_idx = cb_idx[pos_mask]
    vb_idx = vb_idx[pos_mask]
    energy = energy[pos_mask]
    px = px[pos_mask]
    py = py[pos_mask]
    pz = pz[pos_mask]
    f_osc = f_osc[pos_mask]

    # Sort by oscillator strength descending, keep top entries
    sort_order = np.argsort(-f_osc)
    max_show = min(15, len(sort_order))
    sort_order = sort_order[:max_show]

    cb_idx = cb_idx[sort_order]
    vb_idx = vb_idx[sort_order]
    energy = energy[sort_order]
    px = px[sort_order]
    py = py[sort_order]
    pz = pz[sort_order]

    n_trans = len(cb_idx)
    labels = [f"CB{cb_idx[i]}-VB{vb_idx[i]}" for i in range(n_trans)]

    fig, ax = plt.subplots(figsize=(max(8, n_trans * 0.7), 5))
    x = np.arange(n_trans)
    width = 0.25

    ax.bar(x - width, px, width, label=r"$|p_x|^2$", color="#1f77b4", linewidth=0)
    ax.bar(x, py, width, label=r"$|p_y|^2$", color="#2ca02c", linewidth=0)
    ax.bar(x + width, pz, width, label=r"$|p_z|^2$", color="#d62728", linewidth=0)

    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Matrix element (arb. units)")
    ax.set_title(r"GaAs/AlGaAs QW Optical Matrix Elements at $k=0$")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "qw_optical_matrix_elements.png")
    plt.close(fig)
    print("  -> docs/figures/qw_optical_matrix_elements.png")


def fig_qw_potential_profile_gaas(output_dir: Path) -> None:
    """qw_potential_profile_gaas.png: GaAs/AlGaAs QW band-edge profile."""
    print("[figure] qw_potential_profile_gaas")
    cfg = CONFIG_DIR / "qw_gaas_algaas_kpar.cfg"
    result = run_executable(EXE_BAND, cfg, REPO_ROOT, label="qw_gaas_algaas_kpar_pp",
                           timeout=600)
    if result.returncode != 0:
        print("  WARNING: qw_gaas_algaas_kpar run failed, skipping.")
        return
    try:
        z, EV, EV_SO, EC = parse_potential_profile(output_dir)
    except FileNotFoundError:
        print("  WARNING: potential_profile.dat not found, skipping.")
        return

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(z, EC, color="#17becf", linewidth=1.5, label=r"$E_C$")
    ax.plot(z, EV, color="#d62728", linewidth=1.5, label=r"$E_V$")
    ax.plot(z, EV_SO, color="#ff7f0e", linewidth=1.5, label=r"$E_{\Delta SO}$")

    # Dashed vertical lines at material boundaries
    ax.axvline(-50, color="grey", linewidth=0.8, linestyle="--")
    ax.axvline(50, color="grey", linewidth=0.8, linestyle="--")

    # Shade the well region
    well_mask = (z >= -50) & (z <= 50)
    if np.any(well_mask):
        ax.fill_between(z[well_mask], EV[well_mask].min() - 0.05,
                        EC[well_mask].max() + 0.05, alpha=0.08, color="#17becf")

    ax.set_xlabel(r"$z$ (\u00C5)")
    ax.set_ylabel("Energy (eV)")
    ax.set_title(r"GaAs/Al$_{0.3}$Ga$_{0.7}$As QW Band Edge Profile")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "qw_potential_profile_gaas.png")
    plt.close(fig)
    print("  -> docs/figures/qw_potential_profile_gaas.png")


# ===========================================================================
# Main
# ===========================================================================

ALL_FIGURES = {
    "bulk_gaas_bands": fig_bulk_gaas_bands,
    "bulk_gaas_parts": fig_bulk_gaas_parts,
    "bulk_gaas_parts_vs_k": fig_bulk_gaas_parts_vs_k,
    "bulk_gaas_bands_110": fig_bulk_gaas_bands_110,
    "bulk_gaas_warping": fig_bulk_gaas_warping,
    "bulk_gaas_strained_bands": fig_bulk_gaas_strained_bands,
    "bulk_gaas_strain_comparison": fig_bulk_gaas_strain_comparison,
    "bulk_inas_bands": fig_bulk_inas_bands,
    "qw_alsbw_gasbw_inasw_bands": fig_qw_alsbw_gasbw_inasw_bands,
    "qw_potential_profile": fig_qw_potential_profile,
    "qw_wavefunctions": fig_qw_wavefunctions,
    "qw_parts": fig_qw_parts,
    "qw_gaas_algaas_subbands": fig_qw_gaas_algaas_subbands,
    "qcse_stark_shift": fig_qcse_stark_shift,
    "gfactor_components": fig_gfactor_components,
    "gfactor_zeeman": fig_gfactor_zeeman,
    "sc_potential": fig_sc_potential,
    "sc_charge_density": fig_sc_charge_density,
    "sc_convergence": fig_sc_convergence,
    "wire_subbands": fig_wire_subbands,
    "wire_density_2d": fig_wire_density_2d,
    "convergence_fd_order": fig_convergence_fd_order,
    "convergence_grid_spacing": fig_convergence_grid_spacing,
    "timing_dense_vs_sparse": fig_timing_dense_vs_sparse,
    "benchmark_inasw_gasbw_broken_gap": fig_benchmark_inasw_gasbw_broken_gap,
    "qw_dispersion_gaas_algaas": fig_qw_dispersion_gaas_algaas,
    "qw_dispersion_broken_gap": fig_qw_dispersion_broken_gap,
    "qw_optical_matrix_elements": fig_qw_optical_matrix_elements,
    "qw_potential_profile_gaas": fig_qw_potential_profile_gaas,
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
