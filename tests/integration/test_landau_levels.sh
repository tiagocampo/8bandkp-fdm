#!/bin/bash
# COVERAGE: observable=landau_levels geometry=bulk material=InAs
# Integration test: Zeeman splitting for InAs at B=5T
# Verifies conduction band eigenvalues show correct Zeeman splitting
# Delta_E = 2 * g_factor * mu_B * B where mu_B = 5.788e-5 eV/T
# Args: <bandStructure_exe> <config_file>
set -euo pipefail

EXE="$1"
CONFIG="$2"

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

# Setup
/bin/cp "$CONFIG" "$WORKDIR/input.toml"
mkdir -p "$WORKDIR/output"

# Run
cd "$WORKDIR"
RC=0
"$EXE" > test_output.log 2>&1 || RC=$?
if [ $RC -ne 0 ]; then
    echo "FAIL: bandStructure returned exit code $RC"
    cat test_output.log
    exit 1
fi

# Parse eigenvalues from eigenvalues.dat
if [ ! -f "$WORKDIR/output/eigenvalues.dat" ]; then
    echo "FAIL: eigenvalues.dat not produced"
    cat test_output.log
    exit 1
fi

# Extract conduction band eigenvalues (columns 8 and 9) in eV
E_CB1=$(awk 'NR==2 {print $8}' "$WORKDIR/output/eigenvalues.dat")
E_CB2=$(awk 'NR==2 {print $9}' "$WORKDIR/output/eigenvalues.dat")

if [ -z "$E_CB1" ] || [ -z "$E_CB2" ]; then
    echo "FAIL: could not parse eigenvalues"
    cat "$WORKDIR/output/eigenvalues.dat"
    exit 1
fi

echo "Parsed eigenvalues: E_CB1 = $E_CB1 eV, E_CB2 = $E_CB2 eV"

# Expected: Zeeman splitting Delta = 2 * g * mu_B * B
# With g_factor=2.0 (default), mu_B=5.788e-5 eV/T, B=5T:
# Delta = 2 * 2.0 * 5.788e-5 * 5.0 = 1.1576e-3 eV
# Eg(InAs) = 0.417 eV
# E_CB1 = 0.417 - Delta/2 = 0.416421 eV
# E_CB2 = 0.417 + Delta/2 = 0.417579 eV
DELTA_EXPECTED=0.001158
EG_INAS=0.417
ECB1_EXPECTED=$(python3 -c "print($EG_INAS - $DELTA_EXPECTED/2)")
ECB2_EXPECTED=$(python3 -c "print($EG_INAS + $DELTA_EXPECTED/2)")

echo "Expected: E_CB1 = $ECB1_EXPECTED eV, E_CB2 = $ECB2_EXPECTED eV"

TOLERANCE=0.001

PASS=true

DIFF1=$(python3 -c "print(abs($E_CB1 - $ECB1_EXPECTED))")
DIFF2=$(python3 -c "print(abs($E_CB2 - $ECB2_EXPECTED))")

if (( $(echo "$DIFF1 > $TOLERANCE" | bc -l) )); then
    echo "FAIL: E_CB1 = $E_CB1 eV differs from expected $ECB1_EXPECTED eV by $DIFF1 eV"
    PASS=false
fi

if (( $(echo "$DIFF2 > $TOLERANCE" | bc -l) )); then
    echo "FAIL: E_CB2 = $E_CB2 eV differs from expected $ECB2_EXPECTED eV by $DIFF2 eV"
    PASS=false
fi

if [ "$PASS" = true ]; then
    echo "PASS: Zeeman splitting eigenvalues within tolerance"
    echo "E_CB1 = $E_CB1 eV (expected $ECB1_EXPECTED eV)"
    echo "E_CB2 = $E_CB2 eV (expected $ECB2_EXPECTED eV)"
    exit 0
else
    exit 1
fi
