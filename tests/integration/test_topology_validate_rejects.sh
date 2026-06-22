#!/bin/bash
# U8: BdG solver-window validate guard. A BdG config with a Gershgorin-scale
# [solver] window must error stop at validate_semantic.
# Args: <topologicalAnalysis_exe>
set -uo pipefail
EXE="$1"
WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT
PASS=0
FAIL=0

run_test() {
  local name="$1"
  local pat="$2"
  "$EXE" > "$WORKDIR/test_output.log" 2>&1 && RC=0 || RC=$?
  if [ "$RC" -ne 0 ] && grep -qi "$pat" "$WORKDIR/test_output.log"; then
    PASS=$((PASS + 1))
  else
    FAIL=$((FAIL + 1))
    echo "FAIL: $name — exit=$RC, pattern='$pat'"
    grep -i 'error\|STOP' "$WORKDIR/test_output.log" | head -3 | sed 's/^/  /'
  fi
}

cd "$WORKDIR"

# T3: BdG mode + Gershgorin-scale solver window -> reject
cat > input.toml << 'EOF'
confinement = "wire"
FDorder = 2
fd_step = 1
[bands]
num_cb = 4
num_vb = 8
[wire]
nx = 9
ny = 9
dx = 5.0
dy = 5.0
[wire.geometry]
shape = "rectangle"
width = 45.0
height = 45.0
[[region]]
material = "InAs"
inner = 0.0
outer = 45.0
[bdg]
mu = 0.66
delta_0 = 0.0002
[topology]
mode = "bdg"
[solver]
method = "FEAST"
mode = "ENERGY"
emin = -70.0
emax = 70.0
EOF
run_test "T3_bdg_wide_window" "BdG solver window"

if [ "$FAIL" -ne 0 ]; then
  echo "FAIL: topology validate rejects ($FAIL failed)"
  exit 1
fi
echo "PASS: topology validate rejects ($PASS case(s))"
exit 0
