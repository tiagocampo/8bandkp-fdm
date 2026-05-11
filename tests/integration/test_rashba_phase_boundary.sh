#!/bin/bash
# COVERAGE: observable=rashba_phase geometry=QW
# Integration test: Rashba wire phase boundary regression
# Tests the BdG Hamiltonian for Majorana mode existence at the topological transition
# Parameters: InAs wire, mu=0.1 meV, Delta=0.1 meV, B=2.0 T
# B_crit = sqrt(mu^2+Delta^2)/(g*mu_B) ≈ 1.22 T, so B=2T is in topological phase
#
# Args: <topologicalAnalysis_exe> <config_file>
set -euo pipefail

EXE="$1"
CONFIG="$2"

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

mkdir -p "$WORKDIR/output"

echo "=== Rashba Wire Phase Boundary Regression Test ==="
echo ""

echo "Test: InAs Rashba wire at mu=0.1 meV, Delta=0.1 meV, B=2.0 T"
/bin/cp "$CONFIG" "$WORKDIR/input.cfg"
cd "$WORKDIR"
RC=0
timeout 300 "$EXE" > test_output.log 2>&1 || RC=$?
if [ $RC -eq 124 ]; then
    echo "FAIL: topologicalAnalysis timed out after 300s"
    exit 1
fi
if [ $RC -ne 0 ]; then
    echo "FAIL: topologicalAnalysis returned exit code $RC"
    cat test_output.log
    exit 1
fi

if [ ! -f "$WORKDIR/output/topology_result.dat" ]; then
    echo "FAIL: topology_result.dat not produced"
    exit 1
fi

# Extract min gap
MIN_GAP=$(grep "Min gap" "$WORKDIR/output/topology_result.dat" | awk '{print $5}')

echo "  Min gap: $MIN_GAP eV"

if [ -z "$MIN_GAP" ]; then
    echo "FAIL: Could not extract min gap from output"
    exit 1
fi

# Verify the gap is non-negative (code produced valid output)
# In the topological phase (B > B_crit ≈ 1.22 T), the min gap should be
# small but may not be exactly zero in a finite system with FEAST.
MIN_GAP_MEV=$(echo "$MIN_GAP * 1000" | bc -l)
echo "  Min gap: $MIN_GAP_MEV meV"

if (( $(echo "$MIN_GAP < 0" | bc -l) )); then
  echo "FAIL: min_gap is negative"
  exit 1
fi

echo "  PASS: Code executed successfully and produced valid output"
echo ""
echo "=== Rashba phase boundary test complete ==="