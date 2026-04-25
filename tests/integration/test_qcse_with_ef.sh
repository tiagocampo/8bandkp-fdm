#!/bin/bash
# Integration test: QCSE GaAs/AlGaAs QW with electric field from config EFParams
# Validates both regression output equality and config-derived Stark shift semantics
set -euo pipefail

EXE="$1"
CONFIG="$2"
REF_DIR="$3"
COMPARE="$4"
VERIFY_STARK="$(dirname "$0")/verify_stark_shift.py"
ZERO_FIELD_REF="$(dirname "$REF_DIR")/sc_qcse_gaas_algaas/eigenvalues.dat"

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
python3 "$VERIFY_STARK" "$CONFIG" "$ZERO_FIELD_REF" "$WORKDIR/output/eigenvalues.dat"

# QW benchmark checks on the zero-field reference data
VERIFY_QW="$(dirname "$0")/verify_qw_benchmarks.py"
if [ -f "$VERIFY_QW" ] && [ -f "$ZERO_FIELD_REF" ]; then
    python3 "$VERIFY_QW" "$ZERO_FIELD_REF" "$CONFIG"
fi
