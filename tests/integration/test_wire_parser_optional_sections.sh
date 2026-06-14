#!/bin/bash
# COVERAGE: observable=CB_ground_state geometry=wire material=GaAs
# Integration test: optional optics/exciton/scattering sections must not
# consume subsequent solver or strain settings in wire configs.
#
# Verifies that [solver] emin/emax are used (via "Manual energy window" print)
# and that strain computation runs (strain.dat produced).
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
    /bin/cp "$cfg" "$workdir/input.toml"
    sed -i 's/^nsteps = .*/nsteps = 2/' "$workdir/input.toml"
    sed -i 's/^nx = .*/nx = 7/' "$workdir/input.toml"
    sed -i 's/^ny = .*/ny = 7/' "$workdir/input.toml"
    mkdir -p "$workdir/output"
    (
        cd "$workdir"
        "$EXE" > test_output.log 2>&1
    )
    echo "$workdir"
}

WORKDIR_DENSE=$(run_case "$DENSE_CFG")
trap 'rm -rf "$WORKDIR_DENSE" "${WORKDIR_STRAIN:-}"' EXIT

# Verify solver energy window was parsed by checking the "Manual energy window"
# print which includes the solver emin/emax values.
if ! grep -q "Manual energy window" "$WORKDIR_DENSE/test_output.log"; then
    echo "FAIL: solver emin/emax not used (no Manual energy window print)"
    cat "$WORKDIR_DENSE/test_output.log"
    exit 1
fi

# Verify the energy window contains the expected solver values [-1.5, 2.0]
if ! grep -q "\[  -1.5000000000000000      ,   2.0000000000000000" "$WORKDIR_DENSE/test_output.log"; then
    echo "FAIL: solver emin/emax values incorrect in energy window"
    cat "$WORKDIR_DENSE/test_output.log"
    exit 1
fi

# method = "DENSE" selects dense LAPACK path (not FEAST), which prints
# "Dense eigensolver" or runs through the dense path without FEAST.
# The "Manual energy window" print itself confirms dense path was taken.

WORKDIR_STRAIN=$(run_case "$STRAIN_CFG")

# Verify strain computation ran by checking for strain.dat output
if [ ! -f "$WORKDIR_STRAIN/output/strain.dat" ]; then
    echo "FAIL: strain.dat was not produced for strained wire config"
    cat "$WORKDIR_STRAIN/test_output.log"
    exit 1
fi

# Verify strain output is non-empty (has data rows beyond header)
STRAIN_LINES=$(grep -c -v '^#' "$WORKDIR_STRAIN/output/strain.dat" 2>/dev/null || echo 0)
if [ "$STRAIN_LINES" -lt 1 ]; then
    echo "FAIL: strain.dat has no data rows"
    cat "$WORKDIR_STRAIN/test_output.log"
    exit 1
fi

echo "PASS: wire optional-section parsing preserves solver and strain settings"
