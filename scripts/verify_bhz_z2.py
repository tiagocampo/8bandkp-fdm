#!/usr/bin/env python3
"""
Verify BHZ Z2 invariant vs wire width.
Sweeps wire_width from 40 to 100 Angstrom and plots Z2 transition.

Expected behavior:
- wire_width < 70A: M = +10 meV (trivial) -> Z2 = 0
- wire_width >= 70A: M = -10 meV (topological) -> Z2 = 1

Uses topologicalAnalysis executable with mode=qshe, compute_z2=T.
"""

import argparse
import subprocess
import os
import re
import tempfile
import shutil
import sys
from pathlib import Path

# Configuration
REPO = Path(__file__).resolve().parent.parent
CONFIG_TEMPLATE = REPO / "tests/regression/configs/topology_bhz_z2_trivial.cfg"
OUTPUT_DIR = REPO / "docs" / "lecture" / "figures"

# Sweep parameters
WIDTH_MIN = 40
WIDTH_MAX = 100
WIDTH_STEP = 2


def parse_config_template(template_path: Path):
    """Read the template config file."""
    with open(template_path) as f:
        return f.read()


def run_topological_analysis(exe: Path, width_angstrom, config_template_path: Path):
    """Run topologicalAnalysis for a given wire width."""
    # Create temporary working directory
    workdir = tempfile.mkdtemp(prefix="bhz_z2_sweep_")

    try:
        # Create input config with modified wire_width
        config = parse_config_template(config_template_path)
        config = config.replace("wire_width: 58.0", f"wire_width: {width_angstrom:.1f}")

        input_cfg = Path(workdir) / "input.cfg"
        with open(input_cfg, "w") as f:
            f.write(config)

        output_dir = Path(workdir) / "output"
        output_dir.mkdir()

        # Run topologicalAnalysis
        result = subprocess.run(
            [str(exe.absolute())],
            cwd=workdir,
            capture_output=True,
            text=True,
            timeout=120
        )

        if result.returncode != 0:
            print(f"  Warning: executable returned {result.returncode}")
            return None, result.stdout + result.stderr

        # Read topology_result.dat
        topo_file = output_dir / "topology_result.dat"
        if not topo_file.exists():
            return None, "topology_result.dat not found"

        with open(topo_file) as f:
            content = f.read()

        # Parse Z2 invariant from line like "# Z2 invariant: 0"
        match = re.search(r'Z2 invariant:\s*(\d+)', content)
        if match:
            z2 = int(match.group(1))
            return z2, content
        else:
            return None, f"Z2 not found in output. Content:\n{content[:500]}"

    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser(description='Verify BHZ Z2 invariant computation')
    parser.add_argument('--exe', type=str, default='build/src/topologicalAnalysis',
                        help='Executable path')
    parser.add_argument('--config', type=str, help='Config template file path')
    parser.add_argument('--tolerance', type=float, default=0.1,
                        help='Acceptable deviation %%')
    args = parser.parse_args()
    exe = REPO / args.exe if not Path(args.exe).is_absolute() else Path(args.exe)
    config_template = Path(args.config) if args.config else CONFIG_TEMPLATE

    print("BHZ Z2 Verification Sweep")
    print("=" * 50)
    print(f"Sweeping wire_width from {WIDTH_MIN} to {WIDTH_MAX} Angstrom")
    print(f"Expected: Z2=0 for width < 70A, Z2=1 for width >= 70A")
    print()

    # Check executable exists
    if not exe.exists():
        print(f"Error: {exe} not found.")
        print("Run 'make all' first to build the executables.")
        return 1

    # Sweep widths
    widths = list(range(WIDTH_MIN, WIDTH_MAX + 1, WIDTH_STEP))
    results = []

    for w in widths:
        print(f"  Testing width={w}A... ", end="", flush=True)
        z2, info = run_topological_analysis(exe, w, config_template)
        if z2 is not None:
            print(f"Z2={z2}")
            results.append((w, z2))
        else:
            print(f"FAILED")
            results.append((w, None))
            print(f"    Error: {info[:200]}")

    print()
    print("Results:")
    print("-" * 30)
    for w, z2 in results:
        status = str(z2) if z2 is not None else "ERROR"
        marker = " <-- topological" if z2 == 1 else (" <-- trivial" if z2 == 0 else "")
        print(f"  width={w:3d}A: Z2={status}{marker}")

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Generate text summary file
    summary_file = OUTPUT_DIR / "verify_bhz_z2_sweep.txt"
    with open(summary_file, "w") as f:
        f.write("BHZ Z2 Verification Sweep Results\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Wire width sweep: {WIDTH_MIN} to {WIDTH_MAX} A, step {WIDTH_STEP}\n")
        f.write(f"Expected transition at 70A (M changes from +10 to -10 meV)\n\n")
        f.write("Width (A)  Z2  Phase\n")
        f.write("-" * 25 + "\n")
        for w, z2 in results:
            phase = "topological" if z2 == 1 else ("trivial" if z2 == 0 else "error")
            f.write(f"  {w:3d}      {z2}   {phase}\n")

    print(f"\nSummary written to {summary_file}")

    # Text-based summary
    print()
    print("Z2 vs Wire Width:")
    print("-" * 50)
    for w, z2 in results:
        if z2 is not None:
            phase = "TRIVIAL (Z2=0)" if z2 == 0 else "TOPOLOGICAL (Z2=1)"
            marker = " ***" if (w >= 70 and z2 == 1) or (w < 70 and z2 == 0) else ""
        else:
            phase = "ERROR"
            marker = ""
        print(f"  {w:3d}A: {phase}{marker}")

    # Python plot with matplotlib
    try:
        import matplotlib.pyplot as plt
        import numpy as np

        valid_results = [(w, z2) for w, z2 in results if z2 is not None]
        if valid_results:
            widths_valid, z2_valid = zip(*valid_results)

            fig, ax = plt.subplots(figsize=(10, 6))

            # Step plot to show Z2 transition
            ax.step(widths_valid, z2_valid, where='post', linewidth=2,
                    marker='o', markersize=8, color='navy')

            # Fill areas
            ax.fill_between(widths_valid, 0, z2_valid, step='post',
                           alpha=0.2, color='blue')
            ax.fill_between(widths_valid, z2_valid, 1.5, step='post',
                           alpha=0.2, color='orange')

            # Mark transition region
            ax.axvline(x=70, color='red', linestyle='--', linewidth=1.5,
                      label='Critical width (70A)')

            ax.set_xlabel('Wire Width (Angstrom)', fontsize=12)
            ax.set_ylabel('Z2 Invariant', fontsize=12)
            ax.set_title('BHZ Model: Z2 Invariant vs Wire Width', fontsize=14)
            ax.set_ylim(-0.1, 1.5)
            ax.set_xlim(WIDTH_MIN - 2, WIDTH_MAX + 2)
            ax.set_yticks([0, 1])
            ax.set_yticklabels(['0 (Trivial)', '1 (Topological)'])
            ax.legend(loc='center right')
            ax.grid(True, alpha=0.3)

            # Add annotations
            ax.annotate('Trivial phase\n(M = +10 meV)', xy=(55, 0.15),
                       fontsize=10, ha='center', color='navy')
            ax.annotate('Topological phase\n(M = -10 meV)', xy=(85, 0.85),
                       fontsize=10, ha='center', color='darkorange')

            plt.tight_layout()

            fig_file = OUTPUT_DIR / "verify_bhz_z2_width_sweep.png"
            plt.savefig(fig_file, dpi=150)
            print(f"\nFigure saved to {fig_file}")
            plt.close()

    except ImportError:
        print("\nmatplotlib not available, skipping figure generation")

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())