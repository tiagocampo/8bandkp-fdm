#!/bin/bash
# Integration test: topological sweep mode.
# Args: <topologicalAnalysis_exe> <config_file> <expected_model>
set -euo pipefail

EXE="$(realpath "$1")"
CONFIG="$(realpath "$2")"
EXPECTED_MODEL="$3"

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

/bin/cp "$CONFIG" "$WORKDIR/input.cfg"
mkdir -p "$WORKDIR/output"

cd "$WORKDIR"
"$EXE" > test_output.log 2>&1

if ! grep -q "sweep_model=${EXPECTED_MODEL}" test_output.log; then
    echo "FAIL: expected sweep_model=${EXPECTED_MODEL}"
    cat test_output.log
    exit 1
fi

if [ ! -s "$WORKDIR/output/z2_phase_diagram.dat" ]; then
    echo "FAIL: z2_phase_diagram.dat not produced or empty"
    cat test_output.log
    exit 1
fi

if [ ! -s "$WORKDIR/output/z2_transitions.dat" ]; then
    echo "FAIL: z2_transitions.dat not produced or empty"
    cat test_output.log
    exit 1
fi

ROWS=$(awk '$1 !~ /^#/ && NF >= 4 {n++} END{print n+0}' "$WORKDIR/output/z2_phase_diagram.dat")
if [ "$ROWS" -ne 6 ]; then
    echo "FAIL: expected 6 sweep rows for nB=3,nMu=2, got $ROWS"
    cat "$WORKDIR/output/z2_phase_diagram.dat"
    exit 1
fi

awk '$1 !~ /^#/ && NF >= 4 {if ($4+0 < 0.0) bad=1} END{exit bad ? 1 : 0}' \
    "$WORKDIR/output/z2_phase_diagram.dat" || {
    echo "FAIL: sweep gap contains negative values"
    cat "$WORKDIR/output/z2_phase_diagram.dat"
    exit 1
}

echo "PASS: topology sweep regression (${EXPECTED_MODEL})"
