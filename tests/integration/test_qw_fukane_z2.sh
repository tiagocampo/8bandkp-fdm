#!/bin/bash
# COVERAGE: observable=z2_invariant geometry=QW
# Integration test: QW Fu-Kane Z2 invariant path
# Args: <topologicalAnalysis_exe> <config_file>
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
    echo "FAIL: topologicalAnalysis returned exit code $RC"
    cat test_output.log
    exit 1
fi

if [ ! -f "$WORKDIR/output/topology_result.dat" ]; then
    echo "FAIL: topology_result.dat not produced"
    cat test_output.log
    exit 1
fi

Z2=$(grep "Z2 invariant:" "$WORKDIR/output/topology_result.dat" | awk '{print $4}')
if [ "$Z2" != "0" ]; then
    echo "FAIL: Expected Z2 invariant 0 for symmetric GaAs QW, got '$Z2'"
    cat "$WORKDIR/output/topology_result.dat"
    exit 1
fi

MIN_GAP=$(grep "Min gap" "$WORKDIR/output/topology_result.dat" | awk '{print $5}')
awk -v gap="$MIN_GAP" 'BEGIN { target = 1.597960; tol = 0.005; exit !(gap > target - tol && gap < target + tol) }' || {
    echo "FAIL: Expected QW Fu-Kane min gap near 1.597960 eV, got '$MIN_GAP'"
    cat "$WORKDIR/output/topology_result.dat"
    exit 1
}

echo "PASS: QW Fu-Kane Z2 regression"
cat "$WORKDIR/output/topology_result.dat"
