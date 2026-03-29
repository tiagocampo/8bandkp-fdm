#!/bin/bash
# Integration test: g-factor calculation regression
# Args: <gfactorCalculation_exe> <config_file> <ref_data_dir> <compare_script>
set -euo pipefail

EXE="$1"
CONFIG="$2"
REF_DIR="$3"
COMPARE="$4"

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

/bin/cp "$CONFIG" "$WORKDIR/input.cfg"
mkdir -p "$WORKDIR/output"

cd "$WORKDIR"
"$EXE" > test_output.log 2>&1
RC=$?
if [ $RC -ne 0 ]; then
    echo "FAIL: gfactorCalculation returned exit code $RC"
    cat test_output.log
    exit 1
fi

# Compare g-factor values extracted from stdout
# Reference output has lines like:
#   gx
#   (-16.225364267971994,-2.73e-16)
#    -16.225364267971987
# We compare the real g-values (the third line after gx/gy/gz)

extract_gvalue() {
    grep -A2 "^[[:space:]]*$1" "$2" | tail -1 | awk '{print $1}'
}

REF_FILE="$REF_DIR/output.txt"
TEST_FILE="$WORKDIR/test_output.log"

if [ ! -f "$REF_FILE" ]; then
    echo "FAIL: reference output.txt not found at $REF_FILE"
    exit 1
fi

TOL=1e-6
PASS=true

for comp in gx gy gz; do
    REF_VAL=$(extract_gvalue "$comp" "$REF_FILE")
    TEST_VAL=$(extract_gvalue "$comp" "$TEST_FILE")

    if [ -z "$REF_VAL" ] || [ -z "$TEST_VAL" ]; then
        echo "FAIL: could not extract $comp from output"
        PASS=false
        continue
    fi

    python3 -c "
import sys
r, t = float('$REF_VAL'), float('$TEST_VAL')
tol = float('$TOL')
if r == 0:
    err = abs(t)
else:
    err = abs((t - r) / r)
if err > tol:
    print(f'FAIL: $comp ref={r:.10e} test={t:.10e} rel_err={err:.2e}')
    sys.exit(1)
else:
    print(f'PASS: $comp ref={r:.10e} test={t:.10e} rel_err={err:.2e}')
" || PASS=false
done

if $PASS; then
    echo "All g-factor values match"
    exit 0
else
    exit 1
fi
