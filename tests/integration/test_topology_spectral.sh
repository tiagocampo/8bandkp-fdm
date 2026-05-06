#!/bin/bash
# Integration test: topological spectral mode
# Args: <topologicalAnalysis_exe> <config_file> <require_nonzero>
set -euo pipefail

EXE="$(realpath "$1")"
CONFIG="$(realpath "$2")"
REQUIRE_NONZERO="${3:-1}"

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

/bin/cp "$CONFIG" "$WORKDIR/input.cfg"
mkdir -p "$WORKDIR/output"

cd "$WORKDIR"
"$EXE" > test_output.log 2>&1

if [ ! -s "$WORKDIR/output/spectral_function.dat" ]; then
    echo "FAIL: spectral_function.dat not produced or empty"
    cat test_output.log
    exit 1
fi

MAX_A=$(awk 'BEGIN{m=0} $1 !~ /^#/ && NF >= 3 {if ($3+0 > m) m=$3+0} END{printf "%.12e\n", m}' \
    "$WORKDIR/output/spectral_function.dat")

if [ "$REQUIRE_NONZERO" = "1" ]; then
    awk -v max_a="$MAX_A" 'BEGIN{exit !(max_a > 0.0)}' || {
        echo "FAIL: expected nonzero spectral weight, max A=$MAX_A"
        cat test_output.log
        exit 1
    }
fi

echo "PASS: topology spectral regression"
echo "max A=$MAX_A"
