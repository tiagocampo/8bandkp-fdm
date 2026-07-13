#!/bin/bash
# Integration test: BdG LDOS + A(k,E) + Nambu-resolved LDOS (Issue 06 / U9).
# Thin wrapper around verify_bdg_spectral.py (no COVERAGE annotations per
# tests/integration/AGENTS.md convention -- annotations live in the verifier).
set -euo pipefail

EXE="$(realpath "$1")"
CONFIG="$(realpath "$2")"

python3 "$(dirname "$0")/verify_bdg_spectral.py" "$EXE" "$CONFIG"