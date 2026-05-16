#!/bin/bash
# COVERAGE: observable=landau_levels geometry=bulk material=GaAs
# COVERAGE: observable=landau_levels geometry=bulk material=InAs
# Integration test: Landau mode (confinement=3) bulk band structure regression
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

python3 "$COMPARE" "$REF_DIR/eigenvalues.dat" "$WORKDIR/output/eigenvalues.dat" --tolerance 1e-8

# Check for landau_fan.dat if the config uses landau_sweep: B
if grep -q "landau_sweep: B" "$CONFIG" 2>/dev/null || grep -q "landau_sweep:B" "$CONFIG" 2>/dev/null; then
    if [ ! -f "$WORKDIR/output/landau_fan.dat" ]; then
        echo "FAIL: landau_fan.dat not produced for B-sweep config"
        exit 1
    fi
    echo "PASS: landau_fan.dat produced"
fi
