# Issue 03: Parser cleanup — stop 1→error stop, check_optional_stat, stack-alloc table

**Priority:** Important (I1 + I2 + I4)
**Commit:** `fix(parser): replace stop 1 with error stop, add check_optional_stat, stack-alloc table`
**Files:** `src/io/input_parser.f90`, `src/physics/hamiltonian_wire.f90`

## Problem

**I1:** 20 occurrences of `print *, 'Error: ...'; stop 1` in `input_parser.f90` (lines 40, 106, 155, 161, 175, 196, 202, 217, 264, 288, 306, 325, 339, 362, 379, 393, 914, 931, 948, 962). Per CLAUDE.md, `stop 1` without message is deprecated.

**I2:** Polygon vertex coordinates (input_parser.f90:297-300) and b_field components (input_parser.f90:450-455) extract `stat` but silently fall back to 0.0 on type mismatch. Other numeric fields use `check_optional_stat` which errors on type mismatch.

**I4:** `hamiltonian_wire.f90` lines 1019 and 1087 use `type(kp_entry), allocatable :: table(:)` causing heap alloc/dealloc per call in the kz-sweep hot path. `get_kp_block_table()` always returns 52 entries — automatic array suffices.

## Fix

- [ ] **I1:** Replace all 20 `print *, 'Error: ...'; stop 1` with `error stop 'Error: ...'` (single line each)

- [ ] **I2 polygon:** After each `get_value(vertex, 'x', ...)` and `get_value(vertex, 'y', ...)` in `parse_wire`, add:
  ```fortran
  call check_optional_stat(stat, 'x', 'wire.geometry.polygon')
  call check_optional_stat(stat, 'y', 'wire.geometry.polygon')
  ```

- [ ] **I2 b_field:** Replace the 3 `if (stat /= 0) cfg%b_field%components(N) = 0.0_dp` lines with:
  ```fortran
  call get_value(bf_tbl, 'Bx', cfg%b_field%components(1), 0.0_dp, stat=stat)
  call check_optional_stat(stat, 'Bx', 'b_field')
  ```
  (Same for By, Bz)

- [ ] **I4:** In `hamiltonian_wire.f90:1019` and `hamiltonian_wire.f90:1087`, change:
  ```fortran
  type(kp_entry), allocatable :: table(:)
  ```
  to:
  ```fortran
  type(kp_entry) :: table(52)
  ```

## Verification

- [ ] All 107 existing tests pass (`ctest --test-dir build`)
- [ ] No behavioral change for valid configs — parser errors only fire on invalid input
