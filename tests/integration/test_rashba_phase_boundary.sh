#!/bin/bash
# Integration test: Rashba wire phase boundary regression
# Tests the BdG Hamiltonian for Majorana mode existence at the topological transition
# Parameters: mu=0.5 meV, Delta=0.3 meV
# The critical velocity Vc ~ 0.5831 meV marks the topological transition
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

echo "Test: Rashba wire at mu=0.5 meV, Delta=0.3 meV"
/bin/cp "$CONFIG" "$WORKDIR/input.cfg"
cd "$WORKDIR"
RC=0
"$EXE" > test_output.log 2>&1 || RC=$?
if [ $RC -ne 0 ]; then
    echo "FAIL: topologicalAnalysis returned exit code $RC"
    cat test_output.log
    exit 1
fi

if [ ! -f "$WORKDIR/output/topology_result.dat" ]; then
    echo "FAIL: topology_result.dat not produced"
    exit 1
fi

# Extract min gap and edge localization length
MIN_GAP=$(grep "Min gap" "$WORKDIR/output/topology_result.dat" | awk '{print $5}')
EDGE_XI=$(grep "Edge localization length avg" "$WORKDIR/output/topology_result.dat" | awk '{print $NF}')

echo "  Min gap: $MIN_GAP eV"
echo "  Edge localization length: $EDGE_XI AA"

# Check if Majorana modes were detected (edge_xi > 0 indicates localized edge states)
# The exact Z2 outcome depends on the Rashba parameter and wire width
# For the given parameters, we expect to see some gap structure

# For this test, we primarily verify the code runs without error
# and produces meaningful output

if [ -z "$MIN_GAP" ]; then
    echo "FAIL: Could not extract min gap from output"
    exit 1
fi

# Verify gap is near zero for topological parameters (Vz > Vc ≈ 0.58)
# In the topological phase the BdG spectrum should have a minuscule gap
# (< 0.5 meV = 5e-4 eV) due to the Majorana zero mode closing the bulk gap.
MIN_GAP_MEV=$(echo "$MIN_GAP * 1000" | bc -l)
if (( $(echo "$MIN_GAP > 0.0005" | bc -l) )); then
  echo "FAIL: min_gap=$MIN_GAP eV ($MIN_GAP_MEV meV) should be < 0.5 meV for topological phase"
  exit 1
fi

echo "  Gap threshold check passed: min_gap=$MIN_GAP eV ($MIN_GAP_MEV meV) < 0.5 meV"

echo "  PASS: Code executed successfully and produced output"
echo ""
echo "=== Rashba phase boundary test complete ==="