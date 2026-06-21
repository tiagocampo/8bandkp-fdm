#!/bin/bash
# Integration test: wire BdG strain-omission regression (#04).
# Runs the Python verifier which toggles [strain] on/off on a strained
# InAs/GaAs wire BdG config and asserts the BdG spectrum shifts.
# Args: <topologicalAnalysis_exe> <config_file>
set -euo pipefail

EXE="$(realpath "$1")"
CONFIG="$(realpath "$2")"
HERE="$(cd "$(dirname "$0")" && pwd)"

python3 "$HERE/verify_wire_bdg_strain_shift.py" "$EXE" "$CONFIG"
