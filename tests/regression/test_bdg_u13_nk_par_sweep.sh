#!/bin/bash
# U13 T3: nk_par > 1 smoke. Runs a committed wire-bdg fixture with
# nk_par = 4 over a 4-point k_par grid near Gamma in the wire free-z
# direction. Verifies the executable exits 0 and the z2_phase_diagram is
# produced (NOT compared bit-for-bit — this is a smoke test for the new
# dispatch branch only).
#
# Args: <topologicalAnalysis_exe> <config_file>
set -euo pipefail
EXE="$(realpath "$1")"
CONFIG="$(realpath "$2")"

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

/bin/cp "$CONFIG" "$WORKDIR/input.toml"

cd "$WORKDIR"
OMP_NUM_THREADS=4 "$EXE" > run.log 2>&1 || {
  echo "FAIL: topologicalAnalysis exited non-zero on nk_par=4"
  tail -20 run.log | sed 's/^/  /'
  exit 1
}

if [ ! -f output/z2_phase_diagram.dat ]; then
  echo "FAIL: output/z2_phase_diagram.dat not produced (nk_par=4 dispatch failed)"
  exit 1
fi

echo "PASS: nk_par=4 dispatch ran and produced output/z2_phase_diagram.dat"
