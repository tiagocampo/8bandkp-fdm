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

# T3 (U10): wire_BdG path emits Pfaffian-native z2 ∈ {-1, +1, 0}.
# Pre-T3 the remap at main_topology.f90:1420-1423 inverted the sign so the
# column read {0, 1}; gate precondition gate_row_colormap_present
# (verify_majorana_polarization.py:74) needs z2 == -1 to match the wire
# sweep. Pin the native convention here so the integration test fails if
# the remap is reintroduced.
z2_col = [int(round(float(line.split()[2]))) for line in phase_data_rows]
allowed = {-1, 0, 1}
bad_z2 = [v for v in z2_col if v not in allowed]
if bad_z2:
    print(f"FAIL: {PHASE_PATH} z2 column has non-native values {set(bad_z2)}; "
          f"expected subset of {allowed} (Pfaffian-native convention)")
    sys.exit(1)
# Stronger pin: at least one row must read z2 == -1 (the topological cell
# at low B, pre-B_crit). Pre-T3 the remap emits {0, 1} so the value -1
# never appears; this assertion turns the TDD-red into a true fail.
if -1 not in z2_col:
    print(f"FAIL: {PHASE_PATH} z2 column never reads -1 (topological cell); "
          f"got values={sorted(set(z2_col))}; the wire path is still "
          f"emitting the inverted convention")
    sys.exit(1)
non_trivial = [v for v in z2_col if v != 0]
print(f"INFO: phase diagram z2 column native (n={len(z2_col)}, "
      f"non-closure={len(non_trivial)}, values={sorted(set(z2_col))})")

print(f"PASS: phase fixture emits both files under config-driven names "
      f"(phase={len(phase_data_rows)} rows, witness={len(data_rows)} B rows, "
      f"all regex-matched and finite, native z2 column)")
PYEOF
