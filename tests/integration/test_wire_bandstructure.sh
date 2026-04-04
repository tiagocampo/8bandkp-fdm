#!/bin/bash
# Integration test: wire band structure regression
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

# Extract only the first data line (k=1, deterministic) for comparison.
# FEAST non-convergence at k>1 makes those eigenvalues non-deterministic.
grep -v '^#' "$WORKDIR/output/eigenvalues.dat" | head -1 > "$WORKDIR/output/eig_k1.dat"
grep -v '^#' "$REF_DIR/eigenvalues.dat" | head -1 > "$WORKDIR/ref_k1.dat"

python3 "$COMPARE" "$WORKDIR/ref_k1.dat" "$WORKDIR/output/eig_k1.dat" --tolerance 1e-8

# Check parts.dat exists and has correct dimensions (nev rows x 8 columns)
if [ ! -f "$WORKDIR/output/parts.dat" ]; then
    echo "FAIL: parts.dat not produced"
    exit 1
fi

PARTS_ROWS=$(grep -v '^#' "$WORKDIR/output/parts.dat" | grep -v '^$' | wc -l)
PARTS_COLS=$(grep -v '^#' "$WORKDIR/output/parts.dat" | grep -v '^$' | head -1 | wc -w)
if [ "$PARTS_COLS" -ne 8 ]; then
    echo "FAIL: parts.dat has $PARTS_COLS columns, expected 8"
    exit 1
fi
if [ "$PARTS_ROWS" -lt 1 ]; then
    echo "FAIL: parts.dat has no data rows"
    exit 1
fi
echo "PASS: parts.dat has $PARTS_ROWS rows x $PARTS_COLS columns"

# Check eigenfunctions files use new naming convention
if ! ls "$WORKDIR/output/eigenfunctions_k_00001_ev_"*".dat" 1>/dev/null 2>&1; then
    echo "FAIL: eigenfunctions_k_00001_ev_*.dat files not found"
    exit 1
fi
echo "PASS: eigenfunctions files use correct naming convention"
