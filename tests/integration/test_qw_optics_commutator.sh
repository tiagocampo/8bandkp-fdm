#!/bin/bash
# Integration test: QW optical properties regression (commutator-based velocity)
# Compares absorption_TE.dat and absorption_TM.dat against golden reference.
# Args: <opticalProperties_exe> <config_file> <ref_data_dir> <compare_script>
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
    echo "FAIL: opticalProperties returned exit code $RC"
    cat test_output.log
    exit 1
fi

for required in output/absorption_TE.dat output/absorption_TM.dat; do
    if [ ! -f "$required" ]; then
        echo "FAIL: $required not produced"
        exit 1
    fi
done

python3 "$COMPARE" "$REF_DIR/absorption_TE.dat" "$WORKDIR/output/absorption_TE.dat" --tolerance 1e-8
python3 "$COMPARE" "$REF_DIR/absorption_TM.dat" "$WORKDIR/output/absorption_TM.dat" --tolerance 1e-8
