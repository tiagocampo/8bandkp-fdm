#!/bin/bash
# Issue 05 / Unit U7: Dense-QW BdG rung end-to-end certification.
# Thin wrapper around verify_dense_qw_bdg_rung.py (no COVERAGE annotations
# per tests/integration/AGENTS.md convention -- annotations live in the
# verifier).
set -euo pipefail
EXE="$(realpath "$1")"
CONFIG="$(realpath "$2")"
python3 "$(dirname "$0")/verify_dense_qw_bdg_rung.py" "$EXE" "$CONFIG"