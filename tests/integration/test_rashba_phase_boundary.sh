#!/bin/bash
# COVERAGE: observable=rashba_phase geometry=QW
# Integration test: Rashba wire phase boundary regression
# Tests the BdG Hamiltonian for Majorana mode existence at the topological transition
# Parameters: InAs wire, mu=637.928 meV (CB subband), Delta=0.2 meV, B=6.0 T
# B_crit ~ 4.5 T (from 8-band k.p calibration), so B=6T is in topological phase
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

echo "Test: InAs Rashba wire at mu=637.928 meV, Delta=0.2 meV, B=6.0 T"
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

# Check that FEAST found eigenvalues (not the sentinel failure case)
N_EV=$(grep "Eigenvalues found" test_output.log | head -1 | awk '{print $NF}')
if [ -z "$N_EV" ] || [ "$N_EV" -eq 0 ]; then
    echo "FAIL: FEAST found no eigenvalues"
    exit 1
fi
echo "  Eigenvalues found: $N_EV"

# Extract min gap
MIN_GAP=$(grep "Min gap" "$WORKDIR/output/topology_result.dat" | awk '{print $5}')

echo "  Min gap: $MIN_GAP eV"

if [ -z "$MIN_GAP" ]; then
    echo "FAIL: Could not extract min gap from output"
    exit 1
fi

# Verify min_gap is not the sentinel failure value (-1.0)
if (( $(echo "$MIN_GAP < -0.5" | bc -l) )); then
    echo "FAIL: min_gap is sentinel value -1.0 (FEAST failure)"
    exit 1
fi

# Verify the gap is non-negative (code produced valid output)
# In the topological phase (B > B_crit ~ 4.5 T), the min gap should be
# positive but small (reopened gap after phase transition).
MIN_GAP_MEV=$(echo "$MIN_GAP * 1000" | bc -l)
echo "  Min gap: $MIN_GAP_MEV meV"

if (( $(echo "$MIN_GAP < 0" | bc -l) )); then
  echo "FAIL: min_gap is negative"
  exit 1
fi

# Verify gap is in a reasonable range for the topological phase
# With Delta=0.2 meV and B=6T (well above B_crit~4.5T), the gap should
# be between 0.05 meV and 2.0 meV
if (( $(echo "$MIN_GAP_MEV < 0.05" | bc -l) )); then
    echo "FAIL: min_gap too small ($MIN_GAP_MEV meV < 0.05 meV)"
    exit 1
fi
if (( $(echo "$MIN_GAP_MEV > 2.0" | bc -l) )); then
    echo "FAIL: min_gap too large ($MIN_GAP_MEV meV > 2.0 meV)"
    exit 1
fi

echo "  PASS: Code executed successfully and produced valid output"
echo ""
echo "=== Rashba phase boundary test complete ==="
