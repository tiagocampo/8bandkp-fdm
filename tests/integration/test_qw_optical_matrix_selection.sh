#!/bin/bash
# Integration test: QW k=0 optical matrix elements should preserve the
# expected TE/TM selection rules for the simple GaAs/AlGaAs benchmark.
# Args: <gfactorCalculation_exe> <config_file>
set -euo pipefail

EXE="$1"
CONFIG="$2"

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

/bin/cp "$CONFIG" "$WORKDIR/input.cfg"

cd "$WORKDIR"
"$EXE" > test_output.log 2>&1

if [ ! -f output/optical_transitions.dat ]; then
    echo "FAIL: output/optical_transitions.dat not produced"
    cat test_output.log
    exit 1
fi

python3 - output/optical_transitions.dat <<'PY'
import sys
import numpy as np

data = np.loadtxt(sys.argv[1], comments="#")
if data.ndim == 1:
    data = data.reshape(1, -1)

te = data[:, 3] + data[:, 4]
tm = data[:, 5]

hh_like = data[te / np.maximum(tm, 1e-300) > 100.0]
lh_like = data[tm / np.maximum(te, 1e-300) > 100.0]

if hh_like.size == 0:
    raise SystemExit("FAIL: no strongly TE-dominant QW transition found")
if lh_like.size == 0:
    raise SystemExit("FAIL: no strongly TM-dominant QW transition found")

lowest_te = hh_like[np.argmin(hh_like[:, 2])]
lowest_tm = lh_like[np.argmin(lh_like[:, 2])]

print(
    "TE-dominant transition:",
    f"CB{int(lowest_te[0])}->VB{int(lowest_te[1])}",
    f"dE={lowest_te[2]:.6f}",
    f"TE={lowest_te[3] + lowest_te[4]:.6e}",
    f"TM={lowest_te[5]:.6e}",
)
print(
    "TM-dominant transition:",
    f"CB{int(lowest_tm[0])}->VB{int(lowest_tm[1])}",
    f"dE={lowest_tm[2]:.6f}",
    f"TE={lowest_tm[3] + lowest_tm[4]:.6e}",
    f"TM={lowest_tm[5]:.6e}",
)

if not lowest_te[2] < lowest_tm[2]:
    raise SystemExit("FAIL: lowest TE-dominant transition is not below the lowest TM-dominant transition")

print("PASS: QW k=0 optical matrix elements preserve TE/TM ordering")
PY
