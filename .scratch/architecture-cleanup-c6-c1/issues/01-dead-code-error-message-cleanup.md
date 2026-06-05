# Dead Code Removal and Error Message Cleanup (C6)

**Type:** AFK
**Blocked by:** #00 (branch must exist)

## What to Build

Remove dead code and fix misleading error messages across 5 files. All changes are mechanical — deletions and string replacements with zero behavioral impact. Commit as `refactor: dead code removal and error message cleanup`.

Seven items:

1. **Delete `use confinement_init`** from `hamiltonianConstructor.f90` — the import is only referenced in a comment on line 295. Remove the import; leave the comment.

2. **Delete `add_zeeman_coo` subroutine** from `magnetic_field.f90` (lines 17–49) and remove from public exports (line 8). Zero callers. The `stop 1` inside it becomes moot.

3. **Fix 4 misleading "FEAST convergence" messages in zheevx paths** in `main.f90`:
   - Line 645 (QW k-sweep): `'FEAST convergence failed at k-point 1'` → `'zheevx diagonalization failed'`
   - Line 709 (Landau k-sweep): `'FEAST convergence failed at k-point 1'` → `'zheevx diagonalization failed'`
   - Line 786 (Landau B-sweep): `'FEAST convergence failed at k-point 1'` → `'zheevx diagonalization failed in B-sweep'`

4. **Fix hardcoded "k-point 1" on line 327** in `main.f90`: inside the wire k-sweep loop, the `error stop` says "k-point 1" for all k values. Change to `'FEAST convergence failed'` — the actual k value is already printed on line 325. Fortran `error stop` only accepts character literals, so the detail must live in the `print` above.

5. **`stop "error diag"` → `error stop "error diag"`** on `main.f90` line 586. This is a fatal diagonalization error exit.

6. **`stop 'evnum...'` → `error stop 'evnum...'`** on `main_gfactor.f90` line 272. Fatal error exit.

7. **Delete `setup_alloc_sweep` subroutine** from `simulation_setup.f90` (lines 494–509) and remove from public exports (line 26). Zero callers outside its own module.

**Do NOT change line 380** (`stop` with comment "wire mode complete") — this is a normal successful exit, not an error. The early-return pattern for `program` units requires `stop`.

## Acceptance Criteria

- [ ] All 7 items completed
- [ ] `cmake --build build` succeeds
- [ ] `ctest --test-dir build -j4 --output-on-failure` passes all tests
- [ ] `grep -rn 'add_zeeman_coo' src/` returns no results
- [ ] `grep -rn 'setup_alloc_sweep' src/` returns no results
- [ ] `grep -rn "FEAST convergence failed at k-point 1" src/` returns no results
- [ ] Line 380 in `main.f90` is still `stop  ! wire mode complete` (unchanged)
- [ ] Single commit with message `refactor: dead code removal and error message cleanup`

## Blocked by

- #00 (branch must exist)

## User Stories Covered

- #1–4: Correct error messages for diagonalization failures
- #5–7: Dead code removed from modules
- #8–9: Deprecated stop replaced with error stop
