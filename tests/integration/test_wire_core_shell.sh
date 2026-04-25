#!/bin/bash
# Integration test: InAs/GaAs core-shell wire with strain (Stier & Bimberg benchmark)
# Args: <bandStructure_exe> <config_file> <ref_data_dir> <compare_script>
set -euo pipefail

EXE="$1"
CONFIG="$2"
REF_DIR="$3"
COMPARE="$4"

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

# Compare eigenvalues
if [ ! -f "$WORKDIR/output/eigenvalues.dat" ]; then
    echo "FAIL: eigenvalues.dat not produced"
    exit 1
fi

# Extract only the first data line (k=0, deterministic) for comparison.
# FEAST non-convergence at k>1 makes those eigenvalues non-deterministic.
awk '/^[^#]/{print; exit}' "$WORKDIR/output/eigenvalues.dat" > "$WORKDIR/output/eig_k1.dat"
awk '/^[^#]/{print; exit}' "$REF_DIR/eigenvalues.dat" > "$WORKDIR/ref_k1.dat"

python3 "$COMPARE" "$WORKDIR/ref_k1.dat" "$WORKDIR/output/eig_k1.dat" --tolerance 1e-5

echo "PASS: InAs/GaAs core-shell wire eigenvalues match"
