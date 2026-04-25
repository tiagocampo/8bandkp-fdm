#!/bin/bash
# Integration test: interband QW absorption should preserve the expected
# qualitative TE/TM ordering for the simple GaAs/AlGaAs benchmark.
set -euo pipefail

EXE="$1"
CONFIG="$2"
VERIFY="$3"

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

/bin/cp "$CONFIG" "$WORKDIR/input.cfg"

cd "$WORKDIR"
"$EXE" > test_output.log 2>&1
RC=$?
if [ $RC -ne 0 ]; then
    echo "FAIL: bandStructure returned exit code $RC"
    cat test_output.log
    exit 1
fi

for required in output/absorption_TE.dat output/absorption_TM.dat; do
    if [ ! -f "$required" ]; then
        echo "FAIL: $required not produced"
        exit 1
    fi
done

python3 "$VERIFY" output/absorption_TE.dat output/absorption_TM.dat
