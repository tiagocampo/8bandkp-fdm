#!/bin/bash
# COVERAGE: observable=validation geometry=none material=none
# Integration test: verify validate() and validate_semantic() reject bad configs
#
# Each test case writes a minimal TOML config with one invalid field,
# runs bandStructure, and verifies non-zero exit + expected error pattern.
#
# This covers error-stop branches that pFUnit 4.x cannot test
# (no expectErrorStop support — see ADR 0002).
set -uo pipefail

EXE="$1"
WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT
PASS=0
FAIL=0

run_test() {
  local name="$1"
  local expected_pattern="$2"

  "$EXE" > "$WORKDIR/test_output.log" 2>&1 && RC=0 || RC=$?

  if [ "$RC" -ne 0 ] && grep -qi "$expected_pattern" "$WORKDIR/test_output.log"; then
    PASS=$((PASS + 1))
  else
    FAIL=$((FAIL + 1))
    echo "FAIL: $name — exit=$RC, pattern='$expected_pattern'"
    if [ -f "$WORKDIR/test_output.log" ]; then
      grep -i 'error\|STOP' "$WORKDIR/test_output.log" | head -3 | sed 's/^/  /'
    fi
  fi
}

# All tests run in WORKDIR — bandStructure reads input.toml from cwd
cd "$WORKDIR"

# V1: invalid confinement string (caught by parser)
cat > input.toml << 'EOF'
confinement = "invalid"
FDorder = 2
[[material]]
name = "GaAs"
EOF
run_test "V1_bad_confinement" "confinement"

# V2: num_cb = 0 (validate rejects)
cat > input.toml << 'EOF'
confinement = "bulk"
FDorder = 2
[bands]
num_cb = 0
num_vb = 6
[[material]]
name = "GaAs"
EOF
run_test "V2_num_cb_zero" "num_cb must be >= 1"

# V3: num_vb = 0 (validate rejects)
cat > input.toml << 'EOF'
confinement = "bulk"
FDorder = 2
[bands]
num_cb = 2
num_vb = 0
[[material]]
name = "GaAs"
EOF
run_test "V3_num_vb_zero" "num_vb must be >= 1"

# V4: FDorder = 3 (not in {2,4,6,8,10})
cat > input.toml << 'EOF'
confinement = "bulk"
FDorder = 3
[bands]
num_cb = 2
num_vb = 6
[[material]]
name = "GaAs"
EOF
run_test "V4_bad_FDorder" "FDorder must be"

# V5: QW with fd_step = 2 (below minimum 3)
cat > input.toml << 'EOF'
confinement = "qw"
FDorder = 2
fd_step = 2
[bands]
num_cb = 2
num_vb = 6
[[material]]
name = "GaAs"
z_min = 0.0
z_max = 100.0
EOF
run_test "V5_qw_fd_step_2" "fd_step"

# V6: which_band out of range (semantic check)
cat > input.toml << 'EOF'
confinement = "bulk"
FDorder = 2
which_band = 5
[bands]
num_cb = 2
num_vb = 6
[[material]]
name = "GaAs"
EOF
run_test "V6_bad_which_band" "which_band"

# V7: no [[material]] entries
cat > input.toml << 'EOF'
confinement = "bulk"
FDorder = 2
[bands]
num_cb = 2
num_vb = 6
EOF
run_test "V7_no_material" "material"

# V8: wave_vector mode invalid
cat > input.toml << 'EOF'
confinement = "bulk"
FDorder = 2
[bands]
num_cb = 2
num_vb = 6
[wave_vector]
mode = "bad_mode"
max = 0.1
nsteps = 10
[[material]]
name = "GaAs"
EOF
run_test "V8_bad_wv_mode" "wave_vector mode"

echo ""
echo "========================================"
echo " Validation rejection: $PASS passed, $FAIL failed"
echo "========================================"

if [ "$FAIL" -gt 0 ]; then
  exit 1
fi
