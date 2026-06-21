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

GF_EXE="$(dirname "$EXE")/gfactorCalculation"

# V9: bulk gfactor with band_idx out of range
cat > input.toml << 'EOF'
confinement = "bulk"
FDorder = 2
which_band = 0
band_idx = 2
[bands]
num_cb = 2
num_vb = 6
[[material]]
name = "GaAs"
EOF
"$GF_EXE" > "$WORKDIR/test_output.log" 2>&1 && RC=0 || RC=$?
if [ "$RC" -ne 0 ] && grep -qi "bandIdx" "$WORKDIR/test_output.log"; then
  PASS=$((PASS + 1))
else
  FAIL=$((FAIL + 1))
  echo "FAIL: V9_bulk_band_idx_out_of_range — exit=$RC, pattern='bandIdx'"
  grep -i 'error\|STOP' "$WORKDIR/test_output.log" | head -3 | sed 's/^/  /'
fi

# V10: unknown fermi_mode silently falls back — now errors
cat > input.toml << 'EOF'
confinement = "bulk"
FDorder = 2
[bands]
num_cb = 2
num_vb = 6
[sc]
fermi_mode = "typo"
max_iterations = 100
tolerance = 1e-6
[[material]]
name = "GaAs"
EOF
run_test "V10_bad_fermi_mode" "fermi_mode"

# V11: [solver] with invalid method
cat > input.toml << 'EOF'
confinement = "bulk"
FDorder = 2
[bands]
num_cb = 2
num_vb = 6
[solver]
method = "INVALID"
[[material]]
name = "GaAs"
EOF
run_test "V11_bad_solver_method" "solver%method"

# V12: [solver] with invalid mode
cat > input.toml << 'EOF'
confinement = "bulk"
FDorder = 2
[bands]
num_cb = 2
num_vb = 6
[solver]
mode = "INVALID"
[[material]]
name = "GaAs"
EOF
run_test "V12_bad_solver_mode" "solver%mode"

# V13: FEAST solver + INDEX mode (wire — eigensolver_config_validate rejects)
cat > input.toml << 'EOF'
confinement = "wire"
FDorder = 2
fd_step = 1
[wave_vector]
mode = "kz"
max = 0.1
nsteps = 21
[bands]
num_cb = 8
num_vb = 16
[wire]
nx = 21
ny = 21
dx = 3.0
dy = 3.0
[wire.geometry]
shape = "rectangle"
width = 63.0
height = 63.0
[[region]]
material = "GaAs"
inner = 0.0
outer = 100.0
[solver]
method = "FEAST"
mode = "INDEX"
EOF
run_test "V13_feast_index_mode" "FEAST.*INDEX"

# V14: legacy [feast] section must be rejected with a migration message
cat > input.toml << 'EOF'
confinement = "bulk"
FDorder = 2
[[material]]
name = "GaAs"
[feast]
emin = -1.5
emax = 2.0
EOF
run_test "V14_legacy_feast_section" "feast.*removed"

# V15: FEAST + INDEX combination rejected at input validation (defs), not eigensolver
cat > input.toml << 'EOF'
confinement = "wire"
FDorder = 2
[wire]
nx = 5
ny = 5
dx = 3.0
dy = 3.0
[wire.geometry]
shape = "rectangle"
width = 12.0
height = 12.0
[[region]]
material = "GaAs"
inner = 0.0
outer = 100.0
[solver]
method = "FEAST"
mode = "INDEX"
EOF
run_test "V15_feast_index_rejected" "FEAST.*INDEX"

# V16: bulk band count must equal 8 (bulk is fully diagonalized; window < 8
# would write out of bounds). A bulk evnum of 7 (num_cb=2 + num_vb=5) is rejected.
cat > input.toml << 'EOF'
confinement = "bulk"
FDorder = 2
[bands]
num_cb = 2
num_vb = 5
[[material]]
name = "GaAs"
EOF
run_test "V16_bulk_band_count_ne8" "bulk band count"

# Vnew1: partial energy window — emin set, emax left at 0 (auto sentinel).
# validate() must reject: a partial window is ambiguous (0 means auto).
cat > input.toml << 'EOF'
confinement = "bulk"
FDorder = 2
[bands]
num_cb = 2
num_vb = 6
[solver]
method = "DENSE"
emin = -2.0
[[material]]
name = "GaAs"
EOF
run_test "Vnew1_partial_window_emin_only" "both solver emin and emax"

# Vnew2: bulk + FEAST. Bulk is always 8x8 dense by nature; FEAST is
# permanently unsupported for bulk (PRD invariant 1). AUTO->DENSE is fine,
# but an explicit method=FEAST must be rejected.
cat > input.toml << 'EOF'
confinement = "bulk"
FDorder = 2
[bands]
num_cb = 2
num_vb = 6
[solver]
method = "FEAST"
[[material]]
name = "GaAs"
EOF
run_test "Vnew2_bulk_feast" "bulk is always"

echo ""
echo "========================================"
echo " Validation rejection: $PASS passed, $FAIL failed"
echo "========================================"

if [ "$FAIL" -gt 0 ]; then
  exit 1
fi
