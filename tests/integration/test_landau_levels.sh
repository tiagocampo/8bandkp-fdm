#!/bin/bash
# Integration test: Landau level eigenvalues for InAs at B=5T
# Verifies E_0 ≈ 11.13 meV and E_1 ≈ 33.39 meV (analytical: E_n = hbar*omega_c*(n+1/2))
# Args: <bandStructure_exe> <config_file>
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
    echo "FAIL: bandStructure returned exit code $RC"
    cat test_output.log
    exit 1
fi

# Parse eigenvalues from eigenvalues.dat
# Format: "#k, values" followed by lines: "k_value e1 e2 ... e8"
# Bands 7 and 8 (0-indexed: 6, 7) are conduction bands
if [ ! -f "$WORKDIR/output/eigenvalues.dat" ]; then
    echo "FAIL: eigenvalues.dat not produced"
    cat test_output.log
    exit 1
fi

# Extract first two conduction band eigenvalues (bands 7 and 8)
E0=$(awk 'NR==2 {print $8}' "$WORKDIR/output/eigenvalues.dat")
E1=$(awk 'NR==2 {print $9}' "$WORKDIR/output/eigenvalues.dat")

if [ -z "$E0" ] || [ -z "$E1" ]; then
    echo "FAIL: could not parse eigenvalues"
    cat "$WORKDIR/output/eigenvalues.dat"
    exit 1
fi

echo "Parsed eigenvalues: E_0 = $E0 meV, E_1 = $E1 meV"

# Analytical values: hbar*omega_c = 22.26 meV for InAs at B=5T
# E_0 = 1/2 * hbar*omega_c = 11.13 meV
# E_1 = 3/2 * hbar*omega_c = 33.39 meV
TOLERANCE=1.0

E0_REF=11.13
E1_REF=33.39

PASS=true

DIFF0=$(python3 -c "print(abs($E0 - $E0_REF))")
DIFF1=$(python3 -c "print(abs($E1 - $E1_REF))")

if (( $(echo "$DIFF0 > $TOLERANCE" | bc -l) )); then
    echo "FAIL: E_0 = $E0 meV differs from expected $E0_REF meV by $DIFF0 meV (tolerance = $TOLERANCE meV)"
    PASS=false
fi

if (( $(echo "$DIFF1 > $TOLERANCE" | bc -l) )); then
    echo "FAIL: E_1 = $E1 meV differs from expected $E1_REF meV by $DIFF1 meV (tolerance = $TOLERANCE meV)"
    PASS=false
fi

if [ "$PASS" = true ]; then
    echo "PASS: Landau level eigenvalues within tolerance"
    echo "E_0 = $E0 meV (expected $E0_REF meV)"
    echo "E_1 = $E1 meV (expected $E1_REF meV)"
    exit 0
else
    exit 1
fi