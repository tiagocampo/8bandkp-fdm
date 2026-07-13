#!/bin/bash
# Review-finding #4 regression gate (thin wrapper).
# Routes ctest through verify_qw_bandstructure_dense_feast.py, which runs
# bandStructure twice on the SAME QW config (method=DENSE + method=FEAST)
# and asserts the gap-straddling eigenvalue bands agree within tolerance.
# Per tests/integration/AGENTS.md convention, COVERAGE annotations live in
# the verifier Python file (not in this wrapper).
# Args: <build_dir> <source_dir>
set -euo pipefail
BUILD_DIR="$(realpath "$1")"
SOURCE_DIR="$(realpath "$2")"
python3 "$(dirname "$0")/verify_qw_bandstructure_dense_feast.py" "$BUILD_DIR" "$SOURCE_DIR"
