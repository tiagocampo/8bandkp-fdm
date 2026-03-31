#!/bin/bash
# Integration test: Modulation-doped AlAs/GaAs/AlAs QW self-consistent calculation
# Tests SC convergence with charge neutrality (fermi_mode=0) and spatial doping
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

# Check that SC loop actually converged
if ! grep -q "SC loop converged" test_output.log; then
    echo "FAIL: SC loop did not converge"
    cat test_output.log
    exit 1
fi

python3 "$COMPARE" "$REF_DIR/eigenvalues.dat" "$WORKDIR/output/eigenvalues.dat" --tolerance 1e-8
