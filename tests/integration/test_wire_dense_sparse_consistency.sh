#!/bin/bash
# Integration test: reduced wire repro should produce identical eigenvalue
# ordering for dense fallback and FEAST once branch tracking is deterministic.
# Args: <bandStructure_exe> <dense_config> <sparse_config>
set -euo pipefail

EXE="$1"
DENSE_CFG="$2"
SPARSE_CFG="$3"

DENSE_DIR=$(mktemp -d)
SPARSE_DIR=$(mktemp -d)
trap "rm -rf $DENSE_DIR $SPARSE_DIR" EXIT

/bin/cp "$DENSE_CFG" "$DENSE_DIR/input.cfg"
/bin/cp "$SPARSE_CFG" "$SPARSE_DIR/input.cfg"
mkdir -p "$DENSE_DIR/output" "$SPARSE_DIR/output"

(cd "$DENSE_DIR" && "$EXE" > test_output.log 2>&1)
(cd "$SPARSE_DIR" && "$EXE" > test_output.log 2>&1)

python3 - "$DENSE_DIR/output/eigenvalues.dat" "$SPARSE_DIR/output/eigenvalues.dat" <<'PY'
import sys
from pathlib import Path

def parse(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        rows.append([float(x) for x in s.split()])
    return rows

dense = parse(sys.argv[1])
sparse = parse(sys.argv[2])
if len(dense) != len(sparse):
    raise SystemExit(f"row count mismatch: dense={len(dense)} sparse={len(sparse)}")

tol = 1.0e-8
for i, (drow, srow) in enumerate(zip(dense, sparse), start=1):
    if len(drow) != len(srow):
        raise SystemExit(f"column count mismatch at row {i}: dense={len(drow)} sparse={len(srow)}")
    for j, (dv, sv) in enumerate(zip(drow, srow), start=1):
        if abs(dv - sv) > tol:
            raise SystemExit(
                f"eigenvalue mismatch at row {i}, col {j}: dense={dv:.12g} sparse={sv:.12g}"
            )

print("PASS: dense and FEAST wire sweeps match exactly")
PY
