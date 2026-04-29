#!/bin/bash
# Integration test: optional optics/exciton/scattering sections must not
# consume subsequent FEAST or strain settings in wire configs.
#
# Args:
#   <bandStructure_exe> <wire_dense_cfg> <wire_strain_cfg>
set -euo pipefail

EXE=$(realpath "$1")
DENSE_CFG=$(realpath "$2")
STRAIN_CFG=$(realpath "$3")

run_case() {
    local cfg="$1"
    local workdir
    workdir=$(mktemp -d)
    /bin/cp "$cfg" "$workdir/input.cfg"
    sed -i 's/^waveVectorStep:.*/waveVectorStep: 1/' "$workdir/input.cfg"
    sed -i 's/^wire_nx:.*/wire_nx: 7/' "$workdir/input.cfg"
    sed -i 's/^wire_ny:.*/wire_ny: 7/' "$workdir/input.cfg"
    mkdir -p "$workdir/output"
    (
        cd "$workdir"
        "$EXE" > test_output.log 2>&1
    )
    echo "$workdir"
}

WORKDIR_DENSE=$(run_case "$DENSE_CFG")
trap 'rm -rf "$WORKDIR_DENSE" "${WORKDIR_STRAIN:-}"' EXIT

if ! grep -q "feast_emin:" "$WORKDIR_DENSE/test_output.log"; then
    echo "FAIL: feast_emin was not parsed from wire dense config"
    cat "$WORKDIR_DENSE/test_output.log"
    exit 1
fi

if ! grep -q "feast_emax:" "$WORKDIR_DENSE/test_output.log"; then
    echo "FAIL: feast_emax was not parsed from wire dense config"
    cat "$WORKDIR_DENSE/test_output.log"
    exit 1
fi

if ! grep -q "feast_m0:" "$WORKDIR_DENSE/test_output.log"; then
    echo "FAIL: feast_m0 was not parsed from wire dense config"
    cat "$WORKDIR_DENSE/test_output.log"
    exit 1
fi

# Negative feast_m0 selects dense LAPACK path (not FEAST).
# Verify by checking for the "Manual energy window" print unique to dense path.
if ! grep -q "Manual energy window" "$WORKDIR_DENSE/test_output.log"; then
    echo "FAIL: negative feast_m0 did not select dense LAPACK in wire config"
    cat "$WORKDIR_DENSE/test_output.log"
    exit 1
fi

WORKDIR_STRAIN=$(run_case "$STRAIN_CFG")

if ! grep -q "strain:" "$WORKDIR_STRAIN/test_output.log"; then
    echo "FAIL: strain flag was not parsed from wire strain config"
    cat "$WORKDIR_STRAIN/test_output.log"
    exit 1
fi

if ! grep -q "strain_ref:" "$WORKDIR_STRAIN/test_output.log"; then
    echo "FAIL: strain reference was not parsed from wire strain config"
    cat "$WORKDIR_STRAIN/test_output.log"
    exit 1
fi

if ! grep -q "strain_solver:" "$WORKDIR_STRAIN/test_output.log"; then
    echo "FAIL: strain solver was not parsed from wire strain config"
    cat "$WORKDIR_STRAIN/test_output.log"
    exit 1
fi

if [ ! -f "$WORKDIR_STRAIN/output/strain.dat" ]; then
    echo "FAIL: strain.dat was not produced for strained wire config"
    cat "$WORKDIR_STRAIN/test_output.log"
    exit 1
fi

echo "PASS: wire optional-section parsing preserves FEAST and strain settings"
