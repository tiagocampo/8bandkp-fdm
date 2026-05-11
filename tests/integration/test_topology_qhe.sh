#!/bin/bash
# COVERAGE: observable=chern_number geometry=QW
# Integration test: topological analysis QHE mode (Chern number)
# Args: <topologicalAnalysis_exe> <config_file>
set -euo pipefail

EXE="$1"
CONFIG="$2"

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

# Setup
/bin/cp "$CONFIG" "$WORKDIR/input.cfg"
mkdir -p "$WORKDIR/output"

# Run
cd "$WORKDIR"
"$EXE" > test_output.log 2>&1
RC=$?
if [ $RC -ne 0 ]; then
    echo "FAIL: topologicalAnalysis returned exit code $RC"
    cat test_output.log
    exit 1
fi

# Check output files
if [ ! -f "$WORKDIR/output/topology_result.dat" ]; then
    echo "FAIL: topology_result.dat not produced"
    cat test_output.log
    exit 1
fi

echo "PASS: topologicalAnalysis QHE test"
echo "Results:"
cat "$WORKDIR/output/topology_result.dat"