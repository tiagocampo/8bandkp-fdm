#!/bin/bash
# Issue 07 (U10): 2D wire BdG minigap colormap regression.
# B swept (5 values), mu fixed at conduction-band edge (0.659-0.661 eV).
# AC: colormap is non-flat (gap varies across the grid). Compares against
# golden reference at tests/regression/data/.
#
# Args: <topologicalAnalysis_exe> <config_file> <golden_phase_dat>
set -euo pipefail
EXE="$(realpath "$1")"
CONFIG="$(realpath "$2")"
GOLDEN="$(realpath "$3")"
HERE="$(cd "$(dirname "$0")" && pwd)"

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

# Use a single common warm-up of FEAST for the sweep so the result reproduces
# the saved golden within the file tolerance. Run the executable in the workdir.
WORKDIR_INPUT="$WORKDIR/input.toml"
cp "$CONFIG" "$WORKDIR_INPUT"

cd "$WORKDIR"
OMP_NUM_THREADS=4 "$EXE" > run.log 2>&1 || {
  echo "FAIL: topologicalAnalysis exited non-zero on 2D config"
  tail -20 run.log | sed 's/^/  /'
  exit 1
}

if [ ! -f output/z2_phase_diagram.dat ]; then
  echo "FAIL: output/z2_phase_diagram.dat not produced"
  exit 1
fi

python3 - "$GOLDEN" "output/z2_phase_diagram.dat" <<'PYEOF'
"""Verify the 2D colormap is non-flat and matches the golden reference.

Non-flat: gap values vary across (B, mu) grid points.
Tolerance: 1e-8 eV absolute (golden was generated from the same executable
on the same git revision, so the run is deterministic).
"""
import sys
import math
import re

GOLDEN = sys.argv[1]
ACTUAL = sys.argv[2]


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


golden = read_phase(GOLDEN)
actual = read_phase(ACTUAL)

if len(golden) != len(actual):
    print(f"FAIL: row count mismatch (golden={len(golden)}, actual={len(actual)})")
    sys.exit(1)

# Compare row by row.
max_gap_err = 0.0
for g, a in zip(golden, actual):
    if abs(g[0] - a[0]) > 1e-9 or abs(g[1] - a[1]) > 1e-9:
        print(f"FAIL: grid point mismatch B={g[0]} mu={g[1]} vs B={a[0]} mu={a[1]}")
        sys.exit(1)
    if abs(g[2] - a[2]) > 1e-9:
        print(f"FAIL: z2 mismatch at B={g[0]} mu={g[1]} (golden={int(g[2])}, actual={int(a[2])})")
        sys.exit(1)
    max_gap_err = max(max_gap_err, abs(g[3] - a[3]))

if max_gap_err > 1e-8:
    print(f"FAIL: max gap drift {max_gap_err:.3e} eV exceeds 1e-8 tolerance")
    sys.exit(1)

# Non-flat check: at least 2 distinct gap values across the colormap.
gaps = sorted({round(row[3], 12) for row in actual})
if len(gaps) < 2:
    print(f"FAIL: colormap is flat (all gaps={gaps})")
    sys.exit(1)

# Report span.
gap_min = min(r[3] for r in actual)
gap_max = max(r[3] for r in actual)
print(f"PASS: 2D colormap non-flat "
      f"({len(gaps)} distinct gap values, "
      f"min={gap_min*1000:.3f} meV, "
      f"max={gap_max*1000:.3f} meV, "
      f"span={max_gap_err:.3e} eV)")
PYEOF