#!/usr/bin/env python3
"""Sweep B field for Rashba wire + s-wave pairing to find Majorana transition.

Analytical B_crit = sqrt(mu^2 + Delta^2) / (g * mu_B)
where mu_B = 0.05788 meV/T
"""

import subprocess
from pathlib import Path
import re
import sys

REPO = Path(__file__).resolve().parent.parent
EXE = REPO / "build" / "src" / "topologicalAnalysis"
OUT_DIR = Path("docs/lecture/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def run_bdg(B, workdir):
    """Run BdG mode with given B field, return min_gap in meV."""
    cfg = workdir / "input.cfg"
    output = workdir / "output"
    output.mkdir(exist_ok=True)

    # Field order must match input_parser.f90 topology_block sequential reads.
    # See topology_rashba_phase.cfg for the canonical order.
    config_text = f"""\
waveVector: kz
waveVectorMax: 0.1
waveVectorStep: 11
confinement:  2
FDstep: 1
FDorder: 2
numLayers:  1
wire_nx: 11
wire_ny: 11
wire_dx: 3.0
wire_dy: 3.0
wire_shape: rectangle
wire_width: 50.0
wire_height: 50.0
numRegions: 1
region: GaAs  0.0  100.0
numcb: 4
numvb: 4
ExternalField: 0  EF
EFParams: 0.0005
whichBand: 0
bandIdx: 1
SC: 0
topology: T
mode:  bdg
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
ldos_num_E: 200
bdg: T
mu: 0.0005
delta_0: 0.0003
g_factor: 2.0
B_vec: 0.0 0.0 {B:.1f}
gauge: landau
"""

    cfg.write_text(config_text)

    result = subprocess.run(
        [str(EXE)],
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=120
    )
    if result.returncode != 0:
        print(f"  WARNING: B={B}T returned exit code {result.returncode}", file=sys.stderr)
        if result.stderr:
            print(f"  stderr: {result.stderr[:200]}", file=sys.stderr)
        return None

    topo_file = output / "topology_result.dat"
    if not topo_file.exists():
        print(f"  WARNING: B={B}T - topology_result.dat not produced", file=sys.stderr)
        return None

    content = topo_file.read_text()
    match = re.search(r"# Min gap \(eV\):\s*([0-9.e+-]+)", content)
    if match:
        min_gap_eV = float(match.group(1))
        return min_gap_eV * 1000  # convert to meV
    return None


def main():
    try:
        import matplotlib.pyplot as plt
        import tempfile
    except ImportError:
        print("matplotlib not available")
        return

    B_vals = [b / 10 for b in range(0, 105, 5)]  # 0 to 10 T in 0.5 T steps
    gaps = []

    print("Sweeping B field for Rashba BdG Majorana transition...")
    for B in B_vals:
        workdir = Path(tempfile.mkdtemp(prefix=f"bdg_sweep_{B}_"))
        g = run_bdg(B, workdir)
        gaps.append(g)
        status = f"{g:.4f} meV" if g is not None else "FAILED"
        print(f"  B={B:.1f} T: min_gap={status}")

    # Analytical B_crit
    mu = 0.5  # meV
    Delta = 0.3  # meV
    g_factor = 2.0
    mu_B = 0.05788  # meV/T
    B_crit = (mu**2 + Delta**2)**0.5 / (g_factor * mu_B)
    print(f"\nAnalytical B_crit = {B_crit:.1f} T")

    # Plot
    valid_B = [B for B, g in zip(B_vals, gaps) if g is not None]
    valid_g = [g for g in gaps if g is not None]

    fig, ax = plt.subplots(figsize=(7, 5))
    if valid_B:
        ax.plot(valid_B, valid_g, 'o-', label='Computed min gap', markersize=6)
    ax.axvline(x=B_crit, color='red', linestyle='--', label=f'B_crit={B_crit:.1f} T')
    ax.set_xlabel('Magnetic field B (T)')
    ax.set_ylabel('Min gap (meV)')
    ax.set_title('BdG Majorana Phase Transition (Rashba wire + s-wave pairing)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'rashba_majorana_phase_diagram.png', dpi=150)
    print(f"\nSaved rashba_majorana_phase_diagram.png")

    # Also save data
    data_file = OUT_DIR / 'rashba_majorana_phase_diagram.txt'
    with open(data_file, 'w') as f:
        f.write("# B (T)\tmin_gap (meV)\n")
        for B, g in zip(B_vals, gaps):
            f.write(f"{B:.1f}\t{g if g is not None else 'FAILED'}\n")
    print(f"Saved rashba_majorana_phase_diagram.txt")


if __name__ == "__main__":
    main()