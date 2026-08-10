#!/bin/bash
# U13 T3: nk_par > 1 smoke. Runs the canonical wire-bdg fixture with
# nk_par = 4 over a 4-point k_par grid spanning [-0.05, +0.05] 1/A in
# the wire free-z direction. Verifies the executable exits 0 and the
# z2_phase_diagram is produced (NOT compared bit-for-bit — this is a
# smoke test for the new dispatch branch only).
#
# Args: <topologicalAnalysis_exe> <config_file>
set -euo pipefail
EXE="$(realpath "$1")"
CONFIG="$(realpath "$2")"
HERE="$(cd "$(dirname "$0")" && pwd)"

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

cp "$CONFIG" "$WORKDIR/input.toml"
# Inject nk_par + k_par_min/k_par_max into [bdg] after the delta_0 line.
# Keep k_par_min/max near Γ (±0.001 1/A) so the broadened Bloch stack's
# eigenvalues stay inside the canonical fixture's tight FEAST window
# (±0.005 eV) and FEST converges quickly. The nk_par=4 path is exercised
# end-to-end; full physics-correct sweep coverage is T4's destination.
sed -i '/^delta_0 /a nk_par = 4\nk_par_min = -0.001\nk_par_max = 0.001' "$WORKDIR/input.toml"
# Shrink the sweep grid to nB=1, nMu=1 (single (B, mu) cell). This is a
# smoke test for dispatch wiring only, not sweep coverage; the equiv
# test (regression_bdg_u13_nk_par_equiv) covers the full grid at nk_par=1.
sed -i 's/^gap_sweep_nB = 5/gap_sweep_nB = 1/' "$WORKDIR/input.toml"
sed -i 's/^gap_sweep_nMu = 2/gap_sweep_nMu = 1/' "$WORKDIR/input.toml"

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