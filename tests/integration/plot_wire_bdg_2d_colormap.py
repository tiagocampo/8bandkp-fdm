#!/usr/bin/env python3
"""Issue 11 / Plot 6: Wire BdG 2D minigap colormap.

Reads `output/z2_phase_diagram.dat` (produced by the 2D regression test
`tests/regression/test_wire_bdg_topological_2d.sh`) and plots the
(B, mu) minigap as a 2D colormap. Also overlays the 1D B_crit(wire) =
2.8 T horizontal marker from AE3.

Args (optional): <z2_phase_diagram_path> <output_png>
"""
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # headless
import matplotlib.pyplot as plt
import numpy as np


def read_phase(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 4:
                rows.append([float(x) for x in parts[:4]])
    return rows


def main():
    phase_path = Path(sys.argv[1] if len(sys.argv) > 1 else "output/z2_phase_diagram.dat")
    out_path = Path(sys.argv[2] if len(sys.argv) > 2 else "output/wire_bdg_2d_colormap.png")

    rows = read_phase(phase_path)
    if not rows:
        print(f"FAIL: no rows in {phase_path}")
        sys.exit(1)
    rows = np.array(rows)
    Bs = np.unique(rows[:, 0])
    mus = np.unique(rows[:, 1])
    # Build gap grid (B along rows, mu along columns)
    gap_grid = np.zeros((len(Bs), len(mus)))
    z2_grid = np.zeros((len(Bs), len(mus)), dtype=int)
    for r in rows:
        iB = int(np.where(Bs == r[0])[0][0])
        imu = int(np.where(mus == r[1])[0][0])
        gap_grid[iB, imu] = r[3]
        z2_grid[iB, imu] = int(r[2])

    fig, ax = plt.subplots(figsize=(8, 5))
    # Plot minigap in meV
    im = ax.pcolormesh(mus, Bs, gap_grid * 1000.0, shading='auto', cmap='viridis')
    ax.axhline(2.8, color='red', linestyle='--', linewidth=1.5,
               label='B_crit (wire 1D) = 2.8 T')
    # Overlay Z2=1 (topological) markers
    topo_mask = z2_grid == 1
    if topo_mask.any():
        topo_Bs = Bs[np.where(topo_mask.any(axis=1))[0]]
        for B in topo_Bs:
            ax.scatter(mus, np.full_like(mus, B), marker='o', s=20,
                       facecolors='none', edgecolors='red', linewidths=1.0,
                       alpha=0.6)
    ax.set_xlabel('μ (eV)')
    ax.set_ylabel('B (T)')
    ax.set_title('Wire BdG 2D minigap colormap '
                 '(InAs/GaAs, μ=0.659-0.661 eV, Issue 07 / U10)')
    fig.colorbar(im, ax=ax, label='minigap (meV)')
    ax.legend(loc='upper right')
    fig.savefig(out_path, dpi=120, bbox_inches='tight')
    plt.close(fig)
    print(f"wrote {out_path}")
    print(f"  {len(Bs)} B-points x {len(mus)} mu-points; "
          f"gap range [{gap_grid.min()*1000:.3f}, {gap_grid.max()*1000:.3f}] meV")


if __name__ == "__main__":
    main()