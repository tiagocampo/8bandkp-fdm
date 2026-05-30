#!/bin/bash
# COVERAGE: observable=z2_invariant geometry=QW
# Integration test: BHZ Z2 invariant regression
# Tests the BHZ wire model for trivial (Z2=0) and topological (Z2=1) phases
#
# Args: <topologicalAnalysis_exe> <config_trivial> <config_topological>
set -euo pipefail

EXE="$1"
CONFIG_TRIVIAL="$2"   # d=58A, M=+10 meV -> Z2=0 (trivial)
CONFIG_TOPOL="$3"     # d=70A, M=-10 meV -> Z2=1 (topological)

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

mkdir -p "$WORKDIR/output"

echo "=== BHZ Z2 Invariant Regression Tests ==="
echo ""

# Test 1: Trivial case (Z2=0)
echo "Test 1: BHZ trivial (d=58A, M=+10 meV) -> Z2=0"
/bin/cp "$CONFIG_TRIVIAL" "$WORKDIR/input.toml"
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

# Parse Z2 from line like "# Z2 invariant: 0"
Z2_1=$(grep "Z2 invariant:" "$WORKDIR/output/topology_result.dat" | awk '{print $4}')
echo "  Computed Z2 invariant: Z2=$Z2_1"
if [ "$Z2_1" != "0" ]; then
    echo "FAIL: Expected Z2=0 for trivial BHZ, got Z2=$Z2_1"
    exit 1
fi
echo "  PASS"

# Test 2: Topological case (Z2=1)
echo ""
echo "Test 2: BHZ topological (d=70A, M=-10 meV) -> Z2=1"
/bin/cp "$CONFIG_TOPOL" "$WORKDIR/input.toml"
cd "$WORKDIR"
"$EXE" > test_output.log 2>&1
RC=$?
if [ $RC -ne 0 ]; then
    echo "FAIL: topologicalAnalysis returned exit code $RC"
    cat test_output.log
    exit 1
fi

Z2_2=$(grep "Z2 invariant:" "$WORKDIR/output/topology_result.dat" | awk '{print $4}')
echo "  Computed Z2 invariant: Z2=$Z2_2"
if [ "$Z2_2" != "1" ]; then
    echo "FAIL: Expected Z2=1 for topological BHZ, got Z2=$Z2_2"
    exit 1
fi
echo "  PASS"

echo ""
echo "=== All BHZ Z2 tests passed ==="
echo "Summary: trivial d=58A -> Z2=0, topological d=70A -> Z2=1"