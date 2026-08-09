#!/bin/bash
# U13 T3: nk_par=1 vs absent regression guard. Runs the canonical wire-bdg
# 2D fixture twice: once with [bdg].nk_par ABSENT (U10 default path),
# once with [bdg].nk_par = 1 (T3 explicit knob). Diffs every emitted
# output file byte-for-byte. Guard FAILS if dispatch diverges.
#
# Args: <topologicalAnalysis_exe> <config_file>
set -euo pipefail
EXE="$(realpath "$1")"
CONFIG="$(realpath "$2")"
HERE="$(cd "$(dirname "$0")" && pwd)"

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

# --- Run 1: nk_par absent (U10 baseline) ---
RUN1="$WORKDIR/run1"
mkdir -p "$RUN1"
cp "$CONFIG" "$RUN1/input.toml"
# Strip any pre-existing nk_par line (the canonical 2D fixture has none).
sed -i '/^nk_par /d' "$RUN1/input.toml"
cd "$RUN1"
OMP_NUM_THREADS=4 "$EXE" > run.log 2>&1 || {
  echo "FAIL: topologicalAnalysis exited non-zero on Run 1 (nk_par absent)"
  tail -20 run.log | sed 's/^/  /'
  exit 1
}

# --- Run 2: nk_par = 1 (T3 explicit) ---
RUN2="$WORKDIR/run2"
mkdir -p "$RUN2"
cp "$CONFIG" "$RUN2/input.toml"
# Inject nk_par = 1 into the [bdg] section (after the delta_0 line).
sed -i '/^delta_0 /a nk_par = 1' "$RUN2/input.toml"
cd "$RUN2"
OMP_NUM_THREADS=4 "$EXE" > run.log 2>&1 || {
  echo "FAIL: topologicalAnalysis exited non-zero on Run 2 (nk_par = 1)"
  tail -20 run.log | sed 's/^/  /'
  exit 1
}

# --- Diff all emitted output files ---
OUTPUT_FILES=(
  z2_phase_diagram.dat
  wire_slim_pfaffian_witness.dat
)
DIFF_FOUND=0
for f in "${OUTPUT_FILES[@]}"; do
  if [ ! -f "$RUN1/output/$f" ]; then
    echo "FAIL: Run 1 missing output/$f"
    exit 1
  fi
  if [ ! -f "$RUN2/output/$f" ]; then
    echo "FAIL: Run 2 missing output/$f"
    exit 1
  fi
  if ! diff -q "$RUN1/output/$f" "$RUN2/output/$f" >/dev/null; then
    echo "FAIL: output/$f diverges between nk_par absent and nk_par=1"
    diff "$RUN1/output/$f" "$RUN2/output/$f" | head -20 | sed 's/^/  /'
    DIFF_FOUND=1
  fi
done

if [ $DIFF_FOUND -ne 0 ]; then
  exit 1
fi
echo "PASS: nk_par absent and nk_par=1 produce byte-identical outputs"