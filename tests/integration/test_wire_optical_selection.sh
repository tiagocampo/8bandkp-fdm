#!/bin/bash
# Integration test: wire optical transitions should be built from the states
# nearest the actual band edge, not from the deepest valence states in the FEAST
# window.
set -euo pipefail

EXE="$1"
CONFIG="$2"
VERIFY="$3"

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

for required in output/gfactor.dat output/optical_transitions.dat; do
    if [ ! -f "$required" ]; then
        echo "FAIL: $required not produced"
        exit 1
    fi
done

python3 "$VERIFY" test_output.log output/optical_transitions.dat
