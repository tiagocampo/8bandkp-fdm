#!/bin/bash
# Issue 04 Fix Round 1 (Important 1): error-stop path for majorana_polarization.
#
# Per CLAUDE.md "No silent corrections.", the defensive zero-fill on size
# mismatch was replaced with `error stop`. pFUnit 4.x cannot catch
# Fortran `error stop` (per ADR 0002 + test_parameters.pf precedent), so
# this is verified by spawning a standalone Fortran driver compiled at
# build time and checking the exit code + error message.
#
# Args: <driver_exe>
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
    echo "PASS: $name (exit=$RC)"
  else
    FAIL=$((FAIL + 1))
    echo "FAIL: $name — exit=$RC, pattern='$expected_pattern'"
    if [ -f "$WORKDIR/test_output.log" ]; then
      grep -i 'error\|STOP' "$WORKDIR/test_output.log" | head -3 | sed 's/^/  /'
    fi
  fi
}

# T1: undersized evec_bdg (half_n instead of 2*half_n) must error stop.
run_test "size_mismatch_triggers_error_stop" "eigenvector size mismatch"

if [ "$FAIL" -ne 0 ]; then
  exit 1
fi
echo "All $PASS majorana_polarization error-stop tests passed."
exit 0
