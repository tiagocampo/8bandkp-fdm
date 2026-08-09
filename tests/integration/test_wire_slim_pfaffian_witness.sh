#!/bin/bash
# U10 T1b: smoke test that the per-B min-|Pf| proxy producer
# (write_wire_slim_pfaffian_witness in src/io/outputFunctions.f90) emits
# output/wire_slim_pfaffian_witness.dat with nB rows in the regex format
# the lecture 13 acceptance-gate reader parses (B=val |Pf|=val).
#
# Drives the existing wire_bdg 2D sweep config and parses the emitted file.
# Full non-flat assertion (T4) is a separate ticket.
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

python3 - "$WORKDIR/output/wire_slim_pfaffian_witness.dat" <<'PYEOF'
"""Verify the per-B min-|Pf| proxy file is present + well-formed.

Emerges from write_wire_slim_pfaffian_witness (outputFunctions.f90, U10 T1b)
on the wire_bdg sweep path. Header rows start with '#', data rows are
'B=<value> <space> |Pf|=<value>' (literal-token format per
lecture_13_topological.py:226 regex).
"""
import sys
import re

PATH = sys.argv[1]
try:
    text = open(PATH).read()
except OSError as exc:
    print(f"FAIL: cannot open {PATH}: {exc}")
    sys.exit(1)

# Skip '#' header lines; assert at least one data row matches the lecture
# 13 acceptance-gate reader regex.
data_rows = [
    m for m in re.finditer(r"^B=([\d.eE+-]+)\s+\|Pf\|=([\d.eE+-]+)\s*$",
                            text, flags=re.MULTILINE)
]
if not data_rows:
    print(f"FAIL: {PATH} matched no B=val |Pf|=val rows (lecture 13 regex)")
    sys.exit(1)

# Each B value must be a non-negative number; each |Pf| must be
# non-negative (the per-(B, mu) max-site |Pf| is always >= 0, and the
# per-B min over the mu-window inherits the floor).
bad = []
for m in data_rows:
    B = float(m.group(1))
    pf = float(m.group(2))
    if B < 0.0 or pf < 0.0 or not (pf < float('inf')):
        bad.append((B, pf))
if bad:
    print(f"FAIL: {len(bad)} rows with bad B / |Pf| values: {bad[:3]}")
    sys.exit(1)

print(f"PASS: wire_slim_pfaffian_witness.dat present "
      f"({len(data_rows)} B rows, regex-matched, all finite + non-negative)")
PYEOF
