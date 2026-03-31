#!/bin/bash
# Integration test: QCSE GaAs/AlGaAs QW with electric field (-70 kV/cm)
# Validates Stark shift: CB1 should red-shift relative to zero-field case
set -euo pipefail

EXE="$1"
CONFIG="$2"
REF_DIR="$3"
COMPARE="$4"

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

/bin/cp "$CONFIG" "$WORKDIR/input.cfg"
mkdir -p "$WORKDIR/output"

cd "$WORKDIR"
"$EXE" > test_output.log 2>&1
RC=$?
if [ $RC -ne 0 ]; then
    echo "FAIL: bandStructure returned exit code $RC"
    cat test_output.log
    exit 1
fi

if [ ! -f "$WORKDIR/output/eigenvalues.dat" ]; then
    echo "FAIL: eigenvalues.dat not produced"
    exit 1
fi

python3 "$COMPARE" "$REF_DIR/eigenvalues.dat" "$WORKDIR/output/eigenvalues.dat" --tolerance 1e-8
