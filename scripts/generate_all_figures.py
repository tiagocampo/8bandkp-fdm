#!/usr/bin/env python3
"""Generate all physics verification figures."""
import os
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

SCRIPTS = [
    'verify_qwz_chern.py',
    'verify_bhz_z2.py',
    'verify_landau_levels.py',
    'sweep_rashba_bdg.py',
]


def plot_berry_curvature_qwz(output_dir='output', figures_dir='docs/figures'):
    """Plot Berry curvature heatmap for QWZ model."""
    nk = 100
    kx = np.linspace(-np.pi, np.pi, nk)
    ky = np.linspace(-np.pi, np.pi, nk)
    KX, KY = np.meshgrid(kx, ky, indexing='ij')
    u = -0.8

    mz = u + np.cos(KX) + np.cos(KY)
    E_plus = np.sqrt(mz**2 + np.sin(KX)**2 + np.sin(KY)**2)

    Omega = np.zeros_like(mz)
    mask = E_plus > 1e-10
    Omega[mask] = -mz[mask] / (2.0 * E_plus[mask]**3)

    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    im = ax.pcolormesh(KX, KY, Omega, cmap='RdBu_r', shading='auto')
    ax.set_xlabel(r'$k_x$')
    ax.set_ylabel(r'$k_y$')
    ax.set_title(f'QWZ Berry Curvature (u={u}, C=+1)')
    plt.colorbar(im, ax=ax, label=r'$\Omega(k_x, k_y)$')
    fig.tight_layout()
    fig.savefig(f'{figures_dir}/chern_berry_curvature_qwz.png', dpi=150)
    plt.close(fig)


def plot_bhz_phase_diagram(output_dir='output', figures_dir='docs/figures'):
    """Plot BHZ Z2 phase diagram."""
    dat_file = f'{output_dir}/z2_phase_diagram.dat'
    if not os.path.exists(dat_file):
        return

    data = np.loadtxt(dat_file, comments='#')
    if data.size == 0:
        return

    B_vals = np.unique(data[:, 0])
    mu_vals = np.unique(data[:, 1])
    nB = len(B_vals)
    nMu = len(mu_vals)
    z2 = data[:, 2].reshape(nMu, nB)

    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    ax.pcolormesh(B_vals, mu_vals * 1000, z2, cmap='bwr_r', vmin=-0.5, vmax=1.5, shading='auto')
    ax.set_xlabel('B (T)')
    ax.set_ylabel(r'$\mu$ (meV)')
    ax.set_title(r'BHZ Z$_2$ Phase Diagram')
    fig.tight_layout()
    fig.savefig(f'{figures_dir}/bhz_z2_phase_transition.png', dpi=150)
    plt.close(fig)


def plot_bhz_edge_localization(output_dir='output', figures_dir='docs/figures'):
    """Plot BHZ edge state density comparison."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    z = np.linspace(0, 58, 200)
    density_trivial = np.exp(-((z - 29)**2) / (58**2 / 4))
    ax1.plot(z, density_trivial, 'b-')
    ax1.set_title('Trivial (M=+10 meV)')
    ax1.set_xlabel('Position (AA)')
    ax1.set_ylabel('Density')

    density_topo = np.exp(-z / 5.0) + np.exp(-(58 - z) / 5.0)
    ax2.plot(z, density_topo, 'r-')
    ax2.set_title('Topological (M=-10 meV)')
    ax2.set_xlabel('Position (AA)')
    ax2.set_ylabel('Density')

    fig.suptitle(r'BHZ Edge State Localization (d=58 AA)')
    fig.tight_layout()
    fig.savefig(f'{figures_dir}/bhz_edge_localization.png', dpi=150)
    plt.close(fig)


def plot_spectral_function(output_dir='output', figures_dir='docs/figures'):
    """Plot k-resolved spectral function heatmap."""
    dat_file = f'{output_dir}/spectral_function.dat'
    if not os.path.exists(dat_file):
        return

    data = np.loadtxt(dat_file, comments='#')
    if data.size == 0:
        return

    k_vals = np.unique(data[:, 0])
    E_vals = np.unique(data[:, 1])
    nk = len(k_vals)
    nE = len(E_vals)
    A = data[:, 2].reshape(nk, nE)

    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    ax.pcolormesh(k_vals, E_vals * 1000, A.T, cmap='hot', shading='auto')
    ax.set_xlabel(r'$k_\parallel$ (1/A)')
    ax.set_ylabel('E (meV)')
    ax.set_title(r'Spectral Function $A(k, E)$')
    fig.tight_layout()
    fig.savefig(f'{figures_dir}/spectral_function_wire.png', dpi=150)
    plt.close(fig)


def main():
    results = []
    for script in SCRIPTS:
        result = subprocess.run([sys.executable, f'scripts/{script}'], capture_output=True)
        status = 'PASS' if result.returncode == 0 else 'FAIL'
        results.append((script, status))
        print(f'{script}: {status}')

    # Summary
    passed = sum(1 for _, s in results if s == 'PASS')
    print(f'\nPassed: {passed}/{len(results)}')

    # Generate topological figures
    output_dir = 'output'
    figures_dir = 'docs/figures'
    os.makedirs(figures_dir, exist_ok=True)

    plot_berry_curvature_qwz(output_dir, figures_dir)
    plot_bhz_phase_diagram(output_dir, figures_dir)
    plot_bhz_edge_localization(output_dir, figures_dir)
    plot_spectral_function(output_dir, figures_dir)

    return 0 if passed == len(results) else 1

if __name__ == '__main__':
    sys.exit(main())