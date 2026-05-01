#!/bin/bash
# Integration test: QWZ Chern number regression
# Tests the QWZ model at three values of the mass parameter u:
#   u=-0.8 -> C=+1 (topological, for u < -2)
#   u=0.5  -> C=-1 (topological, for -2 < u < 0)
#   u=2.5  -> C=0  (trivial, for u > 0)
#
# Args: <topologicalAnalysis_exe> <config_u-0.8> <config_u0.5> <config_u2.5>
set -euo pipefail

EXE="$1"
CONFIG_U1="$2"   # u=-0.8
CONFIG_U2="$3"   # u=0.5
CONFIG_U3="$4"   # u=2.5

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

mkdir -p "$WORKDIR/output"

echo "=== QWZ Chern Number Regression Tests ==="
echo ""

# Test 1: u=-0.8 -> C=+1
echo "Test 1: u=-0.8 (expected C=+1)"
/bin/cp "$CONFIG_U1" "$WORKDIR/input.cfg"
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
    exit 1
fi

# Parse Chern number from line like "# Chern number: 0"
CHERN1=$(grep "Chern number:" "$WORKDIR/output/topology_result.dat" | awk '{print $4}')
echo "  Computed Chern number: C=$CHERN1"
if [ "$CHERN1" != "+1" ] && [ "$CHERN1" != "1" ]; then
    echo "FAIL: Expected C=+1 for u=-0.8, got C=$CHERN1"
    echo "  NOTE: QWZ model with u=-0.8 should be in topological phase (C=+1)"
    exit 1
fi
echo "  PASS"

# Test 2: u=0.5 -> C=-1
echo ""
echo "Test 2: u=0.5 (expected C=-1)"
/bin/cp "$CONFIG_U2" "$WORKDIR/input.cfg"
cd "$WORKDIR"
"$EXE" > test_output.log 2>&1
RC=$?
if [ $RC -ne 0 ]; then
    echo "FAIL: topologicalAnalysis returned exit code $RC"
    cat test_output.log
    exit 1
fi

CHERN2=$(grep "Chern number:" "$WORKDIR/output/topology_result.dat" | awk '{print $4}')
echo "  Computed Chern number: C=$CHERN2"
if [ "$CHERN2" != "-1" ]; then
    echo "FAIL: Expected C=-1 for u=0.5, got C=$CHERN2"
    echo "  NOTE: QWZ model with u=0.5 should be in topological phase (C=-1)"
    exit 1
fi
echo "  PASS"

# Test 3: u=2.5 -> C=0
echo ""
echo "Test 3: u=2.5 (expected C=0)"
/bin/cp "$CONFIG_U3" "$WORKDIR/input.cfg"
cd "$WORKDIR"
"$EXE" > test_output.log 2>&1
RC=$?
if [ $RC -ne 0 ]; then
    echo "FAIL: topologicalAnalysis returned exit code $RC"
    cat test_output.log
    exit 1
fi

CHERN3=$(grep "Chern number:" "$WORKDIR/output/topology_result.dat" | awk '{print $4}')
echo "  Computed Chern number: C=$CHERN3"
if [ "$CHERN3" != "0" ]; then
    echo "FAIL: Expected C=0 for u=2.5, got C=$CHERN3"
    echo "  NOTE: QWZ model with u=2.5 should be in trivial phase (C=0)"
    exit 1
fi
echo "  PASS"

echo ""
echo "=== All QWZ Chern number tests passed ==="
echo "Summary: u=-0.8 -> C=+1, u=0.5 -> C=-1, u=2.5 -> C=0"