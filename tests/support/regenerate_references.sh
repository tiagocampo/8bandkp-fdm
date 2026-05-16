#!/bin/bash
# Regenerates krylov_reference_data.f90 from the generator.
# Usage: cmake --build build --target regenerate_krylov_references

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/../../build/tests/support"
GENERATOR="${BUILD_DIR}/generate_krylov_reference"
OUTPUT="${SCRIPT_DIR}/krylov_reference_data.f90"

if [ ! -x "$GENERATOR" ]; then
    echo "Error: generator not found at $GENERATOR" >&2
    echo "Build it first: cmake --build build --target generate_krylov_reference" >&2
    exit 1
fi

echo "Regenerating Krylov reference data..."
ERRFILE=$(mktemp)
"$GENERATOR" 2>"$ERRFILE" | sed -n '/^module krylov/,$p' > "$OUTPUT"
if [ $? -ne 0 ]; then
    echo "Error: generator failed. stderr:" >&2
    cat "$ERRFILE" >&2
    rm -f "$ERRFILE"
    exit 1
fi
rm -f "$ERRFILE"
echo "Done: $(wc -l < "$OUTPUT") lines written to $OUTPUT"
