#!/bin/bash
# U10 T2: smoke test that the new sibling fixture
# wire_inas_gaas_bdg_topological_phase.toml emits BOTH config-driven writer
# outputs (z2_phase_diagram_phase.dat + wire_slim_pfaffian_witness_phase.dat)
# in the lecture 13 acceptance-gate regex format. Mirrors
# test_wire_slim_pfaffian_witness.sh (canonical path).
#
# Args: <topologicalAnalysis_exe> <config_file>
set -euo pipefail
EXE="$(realpath "$1")"
CONFIG="$(realpath "$2")"
HERE="$(cd "$(dirname "$0")" && pwd)"

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

cp "$CONFIG" "$WORKDIR/input.toml"
cd "$WORKDIR"
OMP_NUM_THREADS=4 "$EXE" > run.log 2>&1 || {
    echo "FAIL: topologicalAnalysis exited non-zero on $CONFIG"
    tail -20 run.log | sed 's/^/  /'
    exit 1
}

python3 - "$WORKDIR/output/z2_phase_diagram_phase.dat" \
           "$WORKDIR/output/wire_slim_pfaffian_witness_phase.dat" <<'PYEOF'
"""Verify the new phase-fixture output files are present + well-formed.

Config-driven writer path (T2): both files must exist (overridden names)
and the slim Pfaffian witness file must match the lecture 13 regex
(B=val |Pf|=val). The phase diagram file is asserted to have at least
one (B, mu, z2, gap) data row in the canonical 4-column format.
"""
import sys
import re

PHASE_PATH = sys.argv[1]
WITNESS_PATH = sys.argv[2]

# --- Slim Pfaffian witness file ---
try:
    text = open(WITNESS_PATH).read()
except OSError as exc:
    print(f"FAIL: cannot open {WITNESS_PATH}: {exc}")
    sys.exit(1)

data_rows = [
    m for m in re.finditer(r"^B=([\d.eE+-]+)\s+\|Pf\|=([\d.eE+-]+)\s*$",
                            text, flags=re.MULTILINE)
]
if not data_rows:
    print(f"FAIL: {WITNESS_PATH} matched no B=val |Pf|=val rows (lecture 13 regex)")
    sys.exit(1)

bad = []
for m in data_rows:
    B = float(m.group(1))
    pf = float(m.group(2))
    if B < 0.0 or pf < 0.0 or not (pf < float('inf')):
        bad.append((B, pf))
if bad:
    print(f"FAIL: {len(bad)} rows with bad B / |Pf| values: {bad[:3]}")
    sys.exit(1)

# --- Phase diagram file ---
try:
    phase_text = open(PHASE_PATH).read()
except OSError as exc:
    print(f"FAIL: cannot open {PHASE_PATH}: {exc}")
    sys.exit(1)

phase_data_rows = [
    line for line in phase_text.splitlines()
    if line.strip() and not line.lstrip().startswith('#')
]
if not phase_data_rows:
    print(f"FAIL: {PHASE_PATH} has no data rows")
    sys.exit(1)

for line in phase_data_rows:
    parts = line.split()
    if len(parts) != 4:
        print(f"FAIL: {PHASE_PATH} row has {len(parts)} cols, expected 4: {line}")
        sys.exit(1)

print(f"PASS: phase fixture emits both files under config-driven names "
      f"(phase={len(phase_data_rows)} rows, witness={len(data_rows)} B rows, "
      f"all regex-matched and finite)")
PYEOF