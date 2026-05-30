#!/bin/bash
# COVERAGE: observable=SC_convergence geometry=wire material=GaAs
# Integration test: SC wire self-consistent convergence regression
# Args: <bandStructure_exe> <config_file> <ref_data_dir> <compare_script>
set -euo pipefail

EXE="$1"
CONFIG="$2"
REF_DIR="$3"
COMPARE="$4"

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

# Setup
/bin/cp "$CONFIG" "$WORKDIR/input.toml"
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

# Check SC convergence marker
if ! grep -q "Wire SC loop complete" test_output.log; then
    echo "FAIL: SC loop did not complete (missing 'Wire SC loop complete')"
    cat test_output.log
    exit 1
fi
echo "PASS: SC loop completed"

# Compare eigenvalues
if [ ! -f "$WORKDIR/output/eigenvalues.dat" ]; then
    echo "FAIL: eigenvalues.dat not produced"
    exit 1
fi

# Extract only the first data line (k=1, deterministic) for comparison.
# Use awk to avoid SIGPIPE with pipefail from grep|head.
awk '/^[^#]/{print; exit}' "$WORKDIR/output/eigenvalues.dat" > "$WORKDIR/output/eig_k1.dat"
awk '/^[^#]/{print; exit}' "$REF_DIR/eigenvalues.dat" > "$WORKDIR/ref_k1.dat"

python3 "$COMPARE" "$WORKDIR/ref_k1.dat" "$WORKDIR/output/eig_k1.dat" --tolerance 1e-8

echo "PASS: SC wire eigenvalues match golden data"
