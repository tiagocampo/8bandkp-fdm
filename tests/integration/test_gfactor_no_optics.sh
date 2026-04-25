#!/bin/bash
# Integration test: g-factor run with optics disabled must not emit optical transitions.
# Args: <gfactorCalculation_exe> <config_file>
set -euo pipefail

EXE="$1"
CONFIG="$2"

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

/bin/cp "$CONFIG" "$WORKDIR/input.cfg"
mkdir -p "$WORKDIR/output"

cd "$WORKDIR"
"$EXE" > test_output.log 2>&1
RC=$?
if [ $RC -ne 0 ]; then
    echo "FAIL: gfactorCalculation returned exit code $RC"
    cat test_output.log
    exit 1
fi

if [ ! -f "$WORKDIR/output/gfactor.dat" ]; then
    echo "FAIL: gfactor.dat not produced"
    exit 1
fi

if [ -f "$WORKDIR/output/optical_transitions.dat" ]; then
    echo "FAIL: optical_transitions.dat should not be produced when Optics is disabled"
    exit 1
fi

echo "PASS: optics-disabled gfactor run produced only gfactor outputs"
