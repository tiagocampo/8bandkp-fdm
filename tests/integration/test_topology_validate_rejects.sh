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

# Shared wire baseline for T3/T4/T5 (conf, FD, grid, geometry, single InAs region).
write_base() {
  cat <<'EOF'
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
EOF
}

# T3: BdG mode + Gershgorin-scale solver window -> reject
{ write_base; cat <<'EOF'
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
} > input.toml
run_test "T3_bdg_wide_window" "BdG solver window"

# T4: BdG mode + axial B_vec (Peierls would silently early-return) -> reject
{ write_base; cat <<'EOF'
[bdg]
mu = 0.66
delta_0 = 0.0002
B_vec = [0.0, 0.0, 1.0]
g_factor = 15.0
[topology]
mode = "bdg"
EOF
} > input.toml
run_test "T4_bdg_axial_B" "Bx nonzero"

# T5: sweep mode + wire_bdg + Gershgorin-scale solver window -> reject
{ write_base; cat <<'EOF'
[bdg]
mu = 0.66
delta_0 = 0.0002
B_vec = [1.0, 0.0, 0.0]
g_factor = 15.0
[topology]
mode = "sweep"
sweep_model = "wire_bdg"
gap_sweep_B_min = 0.0
gap_sweep_B_max = 5.0
gap_sweep_nB = 5
gap_sweep_mu_min = 0.65
gap_sweep_mu_max = 0.67
gap_sweep_nMu = 3
compute_gap_sweep = true
[solver]
method = "FEAST"
mode = "ENERGY"
emin = -70.0
emax = 70.0
EOF
} > input.toml
run_test "T5_sweep_wire_bdg_wide_window" "wire_bdg sweep"

if [ "$FAIL" -ne 0 ]; then
  echo "FAIL: topology validate rejects ($FAIL failed)"
  exit 1
fi
echo "PASS: topology validate rejects ($PASS case(s))"
exit 0
