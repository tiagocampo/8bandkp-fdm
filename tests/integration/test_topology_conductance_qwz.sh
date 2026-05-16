#!/bin/bash
# COVERAGE: observable=conductance geometry=QW
# Integration test: QWZ conductance mode
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

if ! grep -q "conductance_method: kubo_chern" test_output.log; then
    echo "FAIL: conductance_method was not parsed from config"
    cat test_output.log
    exit 1
fi
if ! grep -q "berry_nk:" test_output.log; then
    echo "FAIL: berry_nk was not parsed from config"
    cat test_output.log
    exit 1
fi
if ! grep -q "landauer_energy:" test_output.log; then
    echo "FAIL: landauer_energy was not parsed from config"
    cat test_output.log
    exit 1
fi

if [ ! -f "$WORKDIR/output/topology_result.dat" ]; then
    echo "FAIL: topology_result.dat not produced"
    cat test_output.log
    exit 1
fi

CONDUCTANCE_XY=$(grep "Conductance xy" "$WORKDIR/output/topology_result.dat" | awk '{print $5}')
awk -v sigma="$CONDUCTANCE_XY" 'BEGIN { target = 1.0; tol = 1.0e-6; exit !(sigma > target - tol && sigma < target + tol) }' || {
    echo "FAIL: Expected conductance_xy near 1, got '$CONDUCTANCE_XY'"
    cat "$WORKDIR/output/topology_result.dat"
    exit 1
}

echo "PASS: QWZ conductance regression"
cat "$WORKDIR/output/topology_result.dat"

awk '/conductance_method:|berry_nk:|landauer_energy:/ {next} {print} END {print "compute_spectral: F"}' \
    "$CONFIG" > "$WORKDIR/input.cfg"
rm -rf "$WORKDIR/output"
mkdir -p "$WORKDIR/output"
"$EXE" > test_default_method_output.log 2>&1
RC=$?
if [ $RC -ne 0 ]; then
    echo "FAIL: default conductance_method smoke returned exit code $RC"
    cat test_default_method_output.log
    exit 1
fi

if ! grep -q "compute_spectral: F" test_default_method_output.log; then
    echo "FAIL: compute_spectral was not preserved for subsequent parsing"
    cat test_default_method_output.log
    exit 1
fi

CONDUCTANCE_XY=$(grep "Conductance xy" "$WORKDIR/output/topology_result.dat" | awk '{print $5}')
awk -v sigma="$CONDUCTANCE_XY" 'BEGIN { target = 1.0; tol = 1.0e-6; exit !(sigma > target - tol && sigma < target + tol) }' || {
    echo "FAIL: Expected default-method conductance_xy near 1, got '$CONDUCTANCE_XY'"
    cat test_default_method_output.log
    cat "$WORKDIR/output/topology_result.dat"
    exit 1
}

echo "PASS: QWZ conductance default-method smoke"
