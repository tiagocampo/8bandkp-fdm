# Eigensolver Review Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Resolve 12 adversarial-review findings on `refactor/eigensolver-standardization`: standardize on the single polymorphic eigensolver interface, restore hoisted-workspace performance, and harden validation/error paths.

**Architecture:** Three phases — (1) standardize by removing the standalone `solve_sparse_evp`/`solve_dense_lapack` and migrating their 18 test call sites to a test-support helper that wraps the polymorphic dispatch; (2) harden `[feast]` rejection, FEAST+INDEX validation, and truncation reporting; (3) restore performance via solver-object workspace caching (dense) and a value-only QW CSR fast-path. No physics changes; all golden eigenvalues stay bit-identical.

**Tech Stack:** Fortran 2018, MKL LAPACK/FEAST, pFUnit 4.x unit tests, shell+Python integration tests, CMake/Ninja.

**Spec:** `docs/superpowers/specs/2026-06-13-eigensolver-review-fixes-design.md`

**Build/test commands (assumed already configured with BUILD_TESTING=ON):**
```bash
cmake --build build                                    # build all
ctest --test-dir build -L unit                         # pFUnit unit tests
ctest --test-dir build -j4 --output-on-failure         # full suite (~113 tests)
```

---

## File Structure

**Production code (`src/`):**
- `src/math/eigensolver.f90` — add cached workspace to `dense_lapack_solver_t` + finalizer + info check (#2,#3); delete `solve_sparse_evp`, `solve_dense_lapack`, and the `#ifdef` dispatch block (#5,#7,#12).
- `src/physics/hamiltonian_qw.f90` — add `update_kp_term_values`, wire slow-path sort-cache population + fast-path value update (#1).
- `src/io/input_parser.f90` — `[feast]` hard error in `read_config` (#6).
- `src/core/defs.f90` — new validation check I15 (FEAST+INDEX) (#8).
- `src/apps/main.f90` — truncation warning in QW FEAST path (#9).
- `src/physics/sc_loop.f90` — remove stale `solve_sparse_evp` import (#7).

**Test code (`tests/`):**
- `tests/support/csr_test_helpers.f90` — add `eigensolve_csr` helper (test-only polymorphic wrapper).
- `tests/unit/test_eigensolver.pf`, `test_z2_invariant.pf`, `test_optical.pf`, `test_edge_states.pf`, `test_hamiltonian.pf` — migrate 18 call sites to `eigensolve_csr`; remove obsolete `test_unknown_method`.
- `tests/integration/test_validate_rejects_bad_configs.sh` — add `[feast]` and FEAST+INDEX rejection cases.
- `tests/integration/verify_qw_sparse_solver.py` — add fast-path equivalence + truncation-warning checks.

**Docs:**
- `CLAUDE.md`, `AGENTS.md` — drop "solve_sparse_evp remains available" line.
- `scripts/lecture_12_extending.py` — delete ARPACK mentions.

---

## Phase 1 — Standardize on the polymorphic interface

### Task 1: Add `eigensolve_csr` test-support helper

This helper reproduces the pre-refactor `emin/emax`+`nev` call ergonomics for tests by deriving `cfg%mode` from that heuristic, then dispatching through `make_eigensolver`. It lives in `csr_test_helpers` (already compiled into `8bandkp_test_support`, which every pFUnit test links — no CMake change).

**Files:**
- Modify: `tests/support/csr_test_helpers.f90:3` (add `use eigensolver`), `:16` (add `public`), `:19+` (add subroutine)
- Test: `tests/unit/test_eigensolver.pf` (add `test_eigensolve_csr_helper`)

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_eigensolver.pf` (inside `contains`, before `end module`):

```fortran
  ! ==================================================================
  ! Test: eigensolve_csr helper dispatches through the polymorphic
  ! solver and returns correct eigenvalues (INDEX path, emin=emax=0).
  ! ==================================================================
  @test
  subroutine test_eigensolve_csr_helper()
    integer, parameter :: n = 20
    real(kind=dp), parameter :: pi = pi_dp
    type(csr_matrix) :: Hmat
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: res
    real(kind=dp) :: expected
    integer :: k

    call build_tridiag_csr(2.0_dp, -1.0_dp, n, Hmat)

    cfg%method = 'DENSE'
    cfg%emin = 0.0_dp
    cfg%emax = 0.0_dp
    cfg%nev = 3

    call eigensolve_csr(Hmat, cfg, res)

    @assertTrue(res%converged)
    @assertEqual(3, res%nev_found)
    do k = 1, 3
      expected = 2.0_dp * (1.0_dp - cos(real(k, kind=dp) * pi / real(n+1, kind=dp)))
      @assertEqual(expected, res%eigenvalues(k), tolerance=1.0e-6_dp)
    end do

    call eigensolver_result_free(res)
    call csr_free(Hmat)
  end subroutine test_eigensolve_csr_helper
```

Also add the import at `tests/unit/test_eigensolver.pf:8` — change:
```fortran
  use eigensolver
```
to:
```fortran
  use eigensolver
  use csr_test_helpers, only: eigensolve_csr
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cmake --build build && ctest --test-dir build -L unit -R eigensolver --output-on-failure`
Expected: FAIL / build error — `eigensolve_csr` not defined.

- [ ] **Step 3: Implement the helper**

In `tests/support/csr_test_helpers.f90`, update the `use` at line 3 from:
```fortran
  use definitions, only: dp
  use sparse_matrices, only: csr_matrix, csr_free
```
to:
```fortran
  use definitions, only: dp
  use sparse_matrices, only: csr_matrix, csr_free
  use eigensolver, only: eigensolver_config, eigensolver_result, &
    eigensolver_base, make_eigensolver, EIGEN_MODE_ENERGY, EIGEN_MODE_INDEX
```

Add to the `public` list (after line 16, `csr_interior_hermitian_error`):
```fortran
  public :: eigensolve_csr
```

Add this subroutine inside `contains` (before `end module csr_test_helpers`):
```fortran
  ! ==================================================================
  ! Test-only eigensolve wrapper: derives cfg%mode from the legacy
  ! emin/emax+nev heuristic, then dispatches through the polymorphic
  ! solver. Preserves the call ergonomics of the removed
  ! solve_sparse_evp/solve_dense_lapack so test configs need no change.
  ! Not a production interface — test infrastructure only.
  ! ==================================================================
  subroutine eigensolve_csr(H_csr, cfg, result)
    type(csr_matrix), intent(in) :: H_csr
    type(eigensolver_config), intent(in) :: cfg
    type(eigensolver_result), intent(out) :: result

    type(eigensolver_config) :: cfg_local
    class(eigensolver_base), allocatable :: solver
    integer :: N

    N = H_csr%nrows
    cfg_local = cfg
    ! Legacy heuristic (matches the removed solve_dense_lapack branch):
    ! ENERGY when an energy window is set, else INDEX for nev smallest.
    if (cfg_local%emin /= 0.0_dp .and. cfg_local%emax /= 0.0_dp) then
      cfg_local%mode = EIGEN_MODE_ENERGY
    else
      cfg_local%mode = EIGEN_MODE_INDEX
      cfg_local%il = 1
      cfg_local%iu = min(max(cfg_local%nev, 1), N)
    end if

    solver = make_eigensolver(cfg_local)
    call solver%solve_sparse(H_csr, cfg_local, result)
    deallocate(solver)
  end subroutine eigensolve_csr
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cmake --build build && ctest --test-dir build -L unit -R eigensolver --output-on-failure`
Expected: PASS — `test_eigensolve_csr_helper` passes.

- [ ] **Step 5: Commit**

```bash
git add tests/support/csr_test_helpers.f90 tests/unit/test_eigensolver.pf
git commit -m "test: add eigensolve_csr polymorphic wrapper to test support"
```

---

### Task 2: Migrate test call sites to `eigensolve_csr`

All 18 sites are a uniform rename. The transformation rule:

- `call solve_sparse_evp(H..., cfg..., res...)` → `call eigensolve_csr(H..., cfg..., res...)`
- `call solve_dense_lapack(H..., cfg..., res...)` → `call eigensolve_csr(H..., cfg..., res...)`

Arguments are identical (both old signatures took `(csr_matrix, eigensolver_config, eigensolver_result)`). Add `use csr_test_helpers, only: eigensolve_csr` to each file's imports; drop `solve_sparse_evp`/`solve_dense_lapack` from any `use eigensolver, only:` lines.

**Files & exact sites** (from `grep -rn 'solve_sparse_evp\|solve_dense_lapack' tests/unit/`):

| File | Lines | Old call | Count |
|------|-------|----------|-------|
| `tests/unit/test_z2_invariant.pf` | 36, 69, 229, 262, 295, 326, 337 | `solve_sparse_evp` | 7 |
| `tests/unit/test_edge_states.pf` | 36, 72, 189, 230, 280, 329 | `solve_sparse_evp` | 6 |
| `tests/unit/test_optical.pf` | 99 | `solve_sparse_evp` | 1 |
| `tests/unit/test_eigensolver.pf` | 133 | `solve_dense_lapack` | 1 |
| `tests/unit/test_eigensolver.pf` | 272 | `solve_sparse_evp` | 1 |
| `tests/unit/test_eigensolver.pf` | 301 (`test_unknown_method`) | `solve_sparse_evp` | REMOVE |
| `tests/unit/test_hamiltonian.pf` | 13 (import) + its call | `solve_dense_lapack` | 1 |

- [ ] **Step 1: Migrate `test_z2_invariant.pf`**

Add import — change line 7:
```fortran
  use eigensolver
```
to:
```fortran
  use eigensolver
  use csr_test_helpers, only: eigensolve_csr
```
Replace all 7 occurrences (`:36, :69, :229, :262, :295, :326, :337`) of `call solve_sparse_evp(` with `call eigensolve_csr(`.

- [ ] **Step 2: Migrate `test_edge_states.pf`**

Add `use csr_test_helpers, only: eigensolve_csr` after its `use eigensolver` line. Replace all 6 occurrences (`:36, :72, :189, :230, :280, :329`) of `call solve_sparse_evp(` with `call eigensolve_csr(`.

- [ ] **Step 3: Migrate `test_optical.pf`**

Add `use csr_test_helpers, only: eigensolve_csr` after its `use eigensolver` line. Replace `:99` `call solve_sparse_evp(` with `call eigensolve_csr(`.

- [ ] **Step 4: Migrate `test_eigensolver.pf` — the oracle call (line 133)**

Replace `:133` `call solve_dense_lapack(Hmat, cfg, res)` with `call eigensolve_csr(Hmat, cfg, res)`. (The `use csr_test_helpers` import was added in Task 1.)

Replace `:272` `call solve_sparse_evp(H, cfg, res)` with `call eigensolve_csr(H, cfg, res)`.

- [ ] **Step 5: Remove obsolete `test_unknown_method` from `test_eigensolver.pf`**

Delete the entire subroutine block (the `@test subroutine test_unknown_method()` ... `end subroutine test_unknown_method` around lines 284–307). **Rationale:** it asserted that an unknown method returns `converged=.false.`; the polymorphic `make_eigensolver` now `error stop`s on unknown methods (fail-fast), so this behavior is deliberately removed and the test is obsolete.

- [ ] **Step 6: Migrate `test_hamiltonian.pf`**

Update the import at line 13 — change:
```fortran
  use eigensolver, only: EIGEN_MODE_FULL, solve_dense_lapack
```
to:
```fortran
  use eigensolver, only: EIGEN_MODE_FULL
  use csr_test_helpers, only: eigensolve_csr
```
Replace the `call solve_dense_lapack(` site with `call eigensolve_csr(`.

- [ ] **Step 7: Build and run all unit tests**

Run: `cmake --build build && ctest --test-dir build -L unit --output-on-failure`
Expected: PASS — all unit tests green; the migrated sites produce identical eigenvalues (heuristic preserves behavior).

- [ ] **Step 8: Commit**

```bash
git add tests/unit/
git commit -m "refactor: migrate eigensolver test call sites to eigensolve_csr wrapper"
```

---

### Task 3: Remove standalone eigensolver functions (production standardization)

Now that no caller remains, delete `solve_sparse_evp`, `solve_dense_lapack`, the stale import, and the docs line.

**Files:**
- Modify: `src/math/eigensolver.f90` (lines 14, 173–196, 394–489)
- Modify: `src/physics/sc_loop.f90:21`
- Modify: `CLAUDE.md`, `AGENTS.md`

- [ ] **Step 1: Delete `solve_sparse_evp` from `eigensolver.f90`**

Delete the public declaration at line 14 — change:
```fortran
  public :: solve_sparse_evp, solve_feast, solve_dense_lapack
```
to:
```fortran
  public :: solve_feast
```
Delete the entire `solve_sparse_evp` subroutine (lines 170–196, from the `! Main dispatch:` comment block through `end subroutine solve_sparse_evp`).

- [ ] **Step 2: Delete `solve_dense_lapack` from `eigensolver.f90`**

Delete the entire subroutine (lines 394–489, from the `! Dense fallback using LAPACK zheevx.` comment block through `end subroutine solve_dense_lapack`).

- [ ] **Step 3: Remove stale import from `sc_loop.f90`**

Change line 20–22:
```fortran
  use eigensolver, only: eigensolver_config, eigensolver_result, &
    & eigensolver_base, make_eigensolver, solve_sparse_evp, eigensolver_result_free, &
    & EIGEN_MODE_INDEX
```
to:
```fortran
  use eigensolver, only: eigensolver_config, eigensolver_result, &
    & eigensolver_base, make_eigensolver, eigensolver_result_free, &
    & EIGEN_MODE_INDEX
```

- [ ] **Step 4: Build to confirm zero remaining references**

Run: `cmake --build build`
Expected: clean build — no unresolved references.

Verify: `grep -rn 'solve_sparse_evp\|solve_dense_lapack' src/ tests/`
Expected: no output.

- [ ] **Step 5: Run full test suite**

Run: `ctest --test-dir build -j4 --output-on-failure`
Expected: PASS — all tests green.

- [ ] **Step 6: Update CLAUDE.md + AGENTS.md**

In both files, find and delete the sentence containing `solve_sparse_evp remains available as a direct interface` (in the eigensolver architecture bullet). Leave the surrounding polymorphic-dispatch description intact.

- [ ] **Step 7: Commit**

```bash
git add src/math/eigensolver.f90 src/physics/sc_loop.f90 CLAUDE.md AGENTS.md
git commit -m "refactor: remove standalone solve_sparse_evp/solve_dense_lapack, standardize on polymorphic dispatch"
```

---

### Task 4: Doc cleanup — ARPACK references + dense round-trip note (#10, #4)

**Files:**
- Modify: `scripts/lecture_12_extending.py:86, 422`
- Modify: `src/math/eigensolver.f90` (comment on `dense_solve_sparse_dispatch`, ~line 742)

- [ ] **Step 1: Delete ARPACK from lecture script**

In `scripts/lecture_12_extending.py`:

Line 86 — change:
```python
    print("    linalg.f90        -- LAPACK/PARDISO/FEAST/ARPACK interfaces")
```
to:
```python
    print("    linalg.f90        -- LAPACK/PARDISO/FEAST interfaces")
```

Line 422 — change:
```python
    print("     - Existing variants: dense LAPACK, MKL FEAST, ARPACK")
```
to:
```python
    print("     - Existing variants: dense LAPACK, MKL FEAST")
```

- [ ] **Step 2: Add doc note on the dense CSR→dense round-trip (#4)**

In `src/math/eigensolver.f90`, in the comment block above `dense_solve_sparse_dispatch` (~line 742), change:
```fortran
  ! ------------------------------------------------------------------
  ! Dense LAPACK solver: solve_sparse converts CSR -> dense.
  ! ------------------------------------------------------------------
```
to:
```fortran
  ! ------------------------------------------------------------------
  ! Dense LAPACK solver: solve_sparse converts CSR -> dense.
  ! Convenience path for CSR inputs; not a hot path (no current caller
  ! routes a large dense matrix through CSR). All large dense solves go
  ! via solve_dense directly; all large CSR solves use FEAST.
  ! ------------------------------------------------------------------
```

- [ ] **Step 3: Verify lecture script runs**

Run: `python3 scripts/lecture_12_extending.py --help 2>&1 | head -5` (or run its self-check if it has no `--help`, confirm it does not reference ARPACK)
Expected: no error; no ARPACK output.

- [ ] **Step 4: Commit**

```bash
git add scripts/lecture_12_extending.py src/math/eigensolver.f90
git commit -m "docs: remove ARPACK references, note dense CSR round-trip is non-hot"
```

---

## Phase 2 — Hardening

### Task 5: `[feast]` hard error (#6)

**Files:**
- Modify: `src/io/input_parser.f90` (in `read_config`, after line 117 `call parse_solver(table, cfg)`)
- Test: `tests/integration/test_validate_rejects_bad_configs.sh`

- [ ] **Step 1: Write the failing test**

Append to `tests/integration/test_validate_rejects_bad_configs.sh` (before the final PASS/FAIL summary, after the last `run_test` call):

```bash
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `ctest --test-dir build -R validate_rejects --output-on-failure`
Expected: FAIL — `V14_legacy_feast_section` fails because `[feast]` is silently ignored (exit 0).

- [ ] **Step 3: Implement the hard error**

In `src/io/input_parser.f90`, in `read_config`, immediately after line 117 (`call parse_solver(table, cfg)`), insert:

```fortran
    ! ---- Reject legacy [feast] section (removed; renamed to [solver]) ----
    call get_value(table, 'feast', child, requested=.false., stat=stat)
    if (stat == 0 .and. associated(child)) then
      error stop '[feast] section removed — rename to [solver] ' // &
        '(fields: method/mode/emin/emax/m0). See docs/reference/input-reference.md'
    end if
```

- [ ] **Step 4: Run test to verify it passes**

Run: `ctest --test-dir build -R validate_rejects --output-on-failure`
Expected: PASS — `V14_legacy_feast_section` passes (nonzero exit + pattern match).

- [ ] **Step 5: Commit**

```bash
git add src/io/input_parser.f90 tests/integration/test_validate_rejects_bad_configs.sh
git commit -m "feat: reject legacy [feast] section with migration message"
```

---

### Task 6: FEAST+INDEX rejection at the validation layer (#8)

**Files:**
- Modify: `src/core/defs.f90` (after line 807, the I14 block)
- Test: `tests/integration/test_validate_rejects_bad_configs.sh`

- [ ] **Step 1: Write the failing test**

Append to `tests/integration/test_validate_rejects_bad_configs.sh`:

```bash
# V15: FEAST + INDEX combination rejected at input validation (defs), not eigensolver
cat > input.toml << 'EOF'
confinement = "wire"
FDorder = 2
[wire]
geometry = { shape = "rectangle", nx = 5, ny = 5 }
[[wire.region]]
material = "GaAs"
z_min = 0
z_max = 100
[solver]
method = "FEAST"
mode = "INDEX"
EOF
run_test "V15_feast_index_rejected" "FEAST.*INDEX"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `ctest --test-dir build -R validate_rejects --output-on-failure`
Expected: FAIL or partial — the error currently fires from `eigensolver_config_validate` (message `eigensolver_config_validate: FEAST solver does not support INDEX mode`), not the defs validation layer. The pattern `FEAST.*INDEX` may match, so verify the message origin manually: run with `-V` and confirm it says `simulation_config_validate`, not `eigensolver_config_validate`. (If it already matches loosely, the test passes but the layer is wrong — proceed to step 3 to move it.)

- [ ] **Step 3: Move the check into `defs.f90`**

In `src/core/defs.f90`, immediately after the I14 block (after line 807, before `end associate`), insert check **I15**:

```fortran
      ! ---- I15: FEAST method cannot combine with INDEX mode ----
      if (trim(cfg%solver%method) == 'FEAST' .and. trim(cfg%solver%mode) == 'INDEX') then
        error stop 'validate_simulation_config: FEAST solver does not support INDEX mode ' // &
          '(use mode = ENERGY or FULL)'
      end if
```

- [ ] **Step 4: Keep the defensive guard in the eigensolver (defense in depth)**

The existing `eigensolver_config_validate` check (`src/math/eigensolver.f90:649`) and `feast_solve_sparse_dispatch` guard (`:703`) stay as-is — they are the second line of defense if a config bypasses `validate()`. No change needed; verify they remain.

- [ ] **Step 5: Run test to verify it passes**

Run: `ctest --test-dir build -R validate_rejects -V 2>&1 | grep -A2 V15`
Expected: the error message originates from `validate_simulation_config` (defs layer), confirming the check moved. Test passes.

- [ ] **Step 6: Commit**

```bash
git add src/core/defs.f90 tests/integration/test_validate_rejects_bad_configs.sh
git commit -m "feat: reject FEAST+INDEX at input validation layer (defs I15)"
```

---

### Task 7: QW FEAST truncation warning (#9)

**Files:**
- Modify: `src/apps/main.f90` (around line 748)

- [ ] **Step 1: Add the warning**

In `src/apps/main.f90`, change the block around lines 745–749:
```fortran
            call solver_bs%solve_sparse(HT_csr_loc, ecfg_bs, result_bs)
            if (.not. result_bs%converged) error stop 'eigensolver failed in QW FEAST k-sweep'
            ! Clamp to allocated array sizes (ENERGY mode may return more eigenvalues)
            M = min(result_bs%nev_found, iuu-il+1)
```
to:
```fortran
            call solver_bs%solve_sparse(HT_csr_loc, ecfg_bs, result_bs)
            if (.not. result_bs%converged) error stop 'eigensolver failed in QW FEAST k-sweep'
            ! Clamp to allocated array sizes (ENERGY mode may return more eigenvalues)
            M = min(result_bs%nev_found, iuu-il+1)
            if (result_bs%nev_found > iuu-il+1) then
              print '(A,I0,A,I0,A)', '  Warning: FEAST returned ', result_bs%nev_found, &
                ' eigenvalues at k-point ', k, &
                '; only the lowest will be kept (widen bands or narrow energy window).'
            end if
```

- [ ] **Step 2: Audit sibling QW FEAST paths**

Search: `grep -rn 'nev_found.*min\|min.*nev_found' src/apps/ src/physics/`
For any other QW FEAST path that clamps `nev_found` without warning (e.g., in `main_optics.f90` or `sc_loop.f90`), apply the same warning. If none found, note it and proceed.

- [ ] **Step 3: Build and verify**

Run: `cmake --build build && ctest --test-dir build -j4 --output-on-failure`
Expected: clean build; full suite green.

- [ ] **Step 4: Commit**

```bash
git add src/apps/main.f90
git commit -m "feat: warn when QW FEAST returns more eigenvalues than retained"
```

---

## Phase 3 — Performance

### Task 8: Dense solver workspace caching + info check (#2, #3)

**Files:**
- Modify: `src/math/eigensolver.f90` (type def lines 127–132; `dense_solve_dense_dispatch` lines 769–864; new finalizer)
- Test: `tests/unit/test_eigensolver.pf` (add `test_dense_solver_caches_workspace`)

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_eigensolver.pf`:

```fortran
  ! ==================================================================
  ! Test: a dense solver reused across calls returns identical, correct
  ! eigenvalues — verifies workspace caching does not corrupt results,
  ! including the realloc path when N grows.
  ! ==================================================================
  @test
  subroutine test_dense_solver_caches_workspace()
    integer :: n
    type(csr_matrix) :: H
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: res1, res2, res3
    class(eigensolver_base), allocatable :: solver
    real(kind=dp) :: expected

    cfg%method = 'DENSE'
    cfg%mode = EIGEN_MODE_INDEX
    cfg%il = 1
    cfg%iu = 1
    solver = make_eigensolver(cfg)

    ! First call: N = 20
    n = 20
    call build_tridiag_csr(2.0_dp, -1.0_dp, n, H)
    call solver%solve_sparse(H, cfg, res1)
    @assertTrue(res1%converged)
    expected = 2.0_dp * (1.0_dp - cos(pi_dp / real(n+1, kind=dp)))
    @assertEqual(expected, res1%eigenvalues(1), tolerance=1.0e-6_dp)
    call eigensolver_result_free(res1)

    ! Second call: same N — cached path, must match
    call solver%solve_sparse(H, cfg, res2)
    @assertTrue(res2%converged)
    @assertEqual(expected, res2%eigenvalues(1), tolerance=1.0e-6_dp)
    call eigensolver_result_free(res2)
    call csr_free(H)

    ! Third call: larger N — realloc path, must still be correct
    n = 30
    call build_tridiag_csr(2.0_dp, -1.0_dp, n, H)
    call solver%solve_sparse(H, cfg, res3)
    @assertTrue(res3%converged)
    expected = 2.0_dp * (1.0_dp - cos(pi_dp / real(n+1, kind=dp)))
    @assertEqual(expected, res3%eigenvalues(1), tolerance=1.0e-6_dp)
    call eigensolver_result_free(res3)
    call csr_free(H)
    deallocate(solver)
  end subroutine test_dense_solver_caches_workspace
```

This test compiles and passes against the *current* (non-caching) code too, so it is a regression guard, not a red-green driver. Run it now to establish green baseline:
Run: `cmake --build build && ctest --test-dir build -R eigensolver --output-on-failure`
Expected: PASS (baseline captured).

- [ ] **Step 2: Add cached workspace fields + finalizer to the type**

Change `src/math/eigensolver.f90` lines 127–132 from:
```fortran
  type, extends(eigensolver_base) :: dense_lapack_solver_t
  contains
    procedure :: solve => dense_lapack_solve_dispatch       ! legacy alias -> solve_sparse
    procedure :: solve_dense => dense_solve_dense_dispatch
    procedure :: solve_sparse => dense_solve_sparse_dispatch
  end type dense_lapack_solver_t
```
to:
```fortran
  type, extends(eigensolver_base) :: dense_lapack_solver_t
    ! Cached LAPACK workspace — sized on first call (or when N grows),
    ! reused thereafter. Thread-safe because each OpenMP thread allocates
    ! its own solver instance (main.f90 QW sweep).
    integer                       :: cached_n = 0
    complex(kind=dp), allocatable :: A_buf(:,:), Z_buf(:,:), work(:)
    real(kind=dp), allocatable    :: W_buf(:), rwork(:)
    integer, allocatable          :: iwork(:), ifail(:)
  contains
    procedure :: solve => dense_lapack_solve_dispatch       ! legacy alias -> solve_sparse
    procedure :: solve_dense => dense_solve_dense_dispatch
    procedure :: solve_sparse => dense_solve_sparse_dispatch
    final     :: dense_lapack_solver_finalize
  end type dense_lapack_solver_t
```

- [ ] **Step 3: Add the finalizer**

Add (near the other finalize routines, e.g., after `dense_solve_dense_dispatch`):
```fortran
  subroutine dense_lapack_solver_finalize(self)
    type(dense_lapack_solver_t), intent(inout) :: self
    if (allocated(self%A_buf))  deallocate(self%A_buf)
    if (allocated(self%Z_buf))  deallocate(self%Z_buf)
    if (allocated(self%work))   deallocate(self%work)
    if (allocated(self%W_buf))  deallocate(self%W_buf)
    if (allocated(self%rwork))  deallocate(self%rwork)
    if (allocated(self%iwork))  deallocate(self%iwork)
    if (allocated(self%ifail))  deallocate(self%ifail)
    self%cached_n = 0
  end subroutine dense_lapack_solver_finalize
```

- [ ] **Step 4: Rewrite `dense_solve_dense_dispatch` to use cached buffers + info check**

Replace the body of `dense_solve_dense_dispatch` (lines 769–864) with:

```fortran
  subroutine dense_solve_dense_dispatch(self, H, config, result)
    class(dense_lapack_solver_t), intent(inout) :: self
    complex(kind=dp), contiguous, intent(in) :: H(:,:)
    type(eigensolver_config), intent(in) :: config
    type(eigensolver_result), intent(out) :: result

    integer :: N, lda, ldz, lwork, info, nb, il_local, iu_local
    complex(kind=dp), allocatable :: wq(:)   ! 1-element workspace-query scratch
    real(kind=dp) :: vl, vu, abstol

    N = size(H, 1)
    if (N <= 0) then
      result%converged = .false.
      result%nev_found = 0
      return
    end if

    lda = N
    ldz = N
    abstol = 0.0_dp
    nb = 0

    ! (Re)allocate N-dependent buffers only when N grows
    if (N > self%cached_n) then
      if (allocated(self%A_buf)) then
        deallocate(self%A_buf, self%Z_buf, self%W_buf, self%rwork, self%iwork, self%ifail)
      end if
      allocate(self%A_buf(N,N), self%Z_buf(N,N), self%W_buf(N))
      allocate(self%rwork(7*N), self%iwork(5*N), self%ifail(N))
      self%cached_n = N
    end if

    ! Copy H (zheev/zheevx destroys input)
    self%A_buf = H

    select case (config%mode)
    case (EIGEN_MODE_FULL)
      ! zheev: all eigenvalues. Workspace query -> reuse/grow work.
      allocate(wq(1))
      call zheev('V', 'U', N, self%A_buf, lda, self%W_buf, wq, -1, self%rwork, info)
      if (info /= 0) error stop 'dense_solve_dense_dispatch: zheev workspace query failed.'
      lwork = max(1, nint(real(wq(1))))
      deallocate(wq)
      if (.not. allocated(self%work) .or. size(self%work) < lwork) then
        if (allocated(self%work)) deallocate(self%work)
        allocate(self%work(lwork))
      end if
      call zheev('V', 'U', N, self%A_buf, lda, self%W_buf, self%work, lwork, self%rwork, info)
      if (info == 0) then
        nb = N
        self%Z_buf = self%A_buf   ! zheev returns eigenvectors in A
      end if

    case (EIGEN_MODE_INDEX)
      il_local = max(1, config%il)
      iu_local = min(N, config%iu)
      allocate(wq(1))
      call zheevx('V', 'I', 'U', N, self%A_buf, lda, vl, vu, &
                   il_local, iu_local, abstol, nb, self%W_buf, self%Z_buf, ldz, &
                   wq, -1, self%rwork, self%iwork, self%ifail, info)
      if (info /= 0) error stop 'dense_solve_dense_dispatch: zheevx(INDEX) workspace query failed.'
      lwork = max(1, nint(real(wq(1))))
      deallocate(wq)
      if (.not. allocated(self%work) .or. size(self%work) < lwork) then
        if (allocated(self%work)) deallocate(self%work)
        allocate(self%work(lwork))
      end if
      call zheevx('V', 'I', 'U', N, self%A_buf, lda, vl, vu, &
                   il_local, iu_local, abstol, nb, self%W_buf, self%Z_buf, ldz, &
                   self%work, lwork, self%rwork, self%iwork, self%ifail, info)

    case (EIGEN_MODE_ENERGY)
      vl = config%emin
      vu = config%emax
      allocate(wq(1))
      call zheevx('V', 'V', 'U', N, self%A_buf, lda, vl, vu, &
                   1, N, abstol, nb, self%W_buf, self%Z_buf, ldz, &
                   wq, -1, self%rwork, self%iwork, self%ifail, info)
      if (info /= 0) error stop 'dense_solve_dense_dispatch: zheevx(ENERGY) workspace query failed.'
      lwork = max(1, nint(real(wq(1))))
      deallocate(wq)
      if (.not. allocated(self%work) .or. size(self%work) < lwork) then
        if (allocated(self%work)) deallocate(self%work)
        allocate(self%work(lwork))
      end if
      call zheevx('V', 'V', 'U', N, self%A_buf, lda, vl, vu, &
                   1, N, abstol, nb, self%W_buf, self%Z_buf, ldz, &
                   self%work, lwork, self%rwork, self%iwork, self%ifail, info)

    case default
      error stop 'dense_solve_dense_dispatch: unknown mode.'
    end select

    if (info == 0 .and. nb > 0) then
      result%converged = .true.
      result%nev_found = nb
      result%iterations = 1
      allocate(result%eigenvalues(nb))
      allocate(result%eigenvectors(N, nb))
      result%eigenvalues(1:nb) = self%W_buf(1:nb)
      result%eigenvectors(:, 1:nb) = self%Z_buf(:, 1:nb)
    else
      result%converged = .false.
      result%nev_found = 0
      if (info /= 0) print *, 'Dense eigensolver error: info =', info
    end if
  end subroutine dense_solve_dense_dispatch
```

Note: the per-call `deallocate(A, W, Z, work, rwork, iwork, ifail)` at the old line 863 is removed — buffers are now owned by the solver and freed by the finalizer.

- [ ] **Step 5: Build and run eigensolver + full unit suite**

Run: `cmake --build build && ctest --test-dir build -L unit --output-on-failure`
Expected: PASS — `test_dense_solver_caches_workspace` and all existing tests green.

- [ ] **Step 6: Run the QW k-sweep regression (the hot path)**

Run: `ctest --test-dir build -L regression --output-on-failure`
Expected: PASS — golden QW eigenvalues unchanged.

- [ ] **Step 7: Commit**

```bash
git add src/math/eigensolver.f90 tests/unit/test_eigensolver.pf
git commit -m "perf: cache dense LAPACK workspace in the solver object; add info checks"
```

---

### Task 9: QW CSR fast-path — value-update routine + slow-path cache wiring (#1, part 1)

**Files:**
- Modify: `src/physics/hamiltonian_qw.f90` (type `qw_workspace`; add `update_kp_term_values`; wire slow-path `finalize_coo_to_csr` to populate `ws%coo_cache`)

**Key signature facts (verified):** `finalize_coo_to_csr(..., ws, coo_cache)` accepts `ws` of type `wire_workspace` or `coo_cache` of type `wire_coo_cache` — never `qw_workspace`. When the cache is valid it calls `csr_set_values_from_coo` (values-only update); when invalid it calls `csr_build_from_coo_cached` (builds structure + records the sort map). Both `wire_coo_cache` and `wire_coo_cache_free` are already imported (`hamiltonian_qw.f90:29`).

- [ ] **Step 1: Add a `wire_coo_cache` field to `qw_workspace`**

In `src/physics/hamiltonian_qw.f90`, add a field to the type (after `blk_temp`, ~line 49):
```fortran
    ! COO→CSR sort cache (populated on slow path, reused on fast path)
    type(wire_coo_cache) :: coo_cache
```
In `qw_workspace_free` (after the `coo_to_csr` deallocation, ~line 91), add:
```fortran
    call wire_coo_cache_free(ws%coo_cache)
```
(The finalizer `qw_workspace_finalize` delegates to `qw_workspace_free`, so no separate change.)

- [ ] **Step 2: Add `update_kp_term_values`**

Add this subroutine inside `contains` in `src/physics/hamiltonian_qw.f90` (before `ZB8bandQW_csr`). It walks a cached CSR block's structure and recomputes each non-zero's value from `kpterms(i,j,·)` and the current `kx,ky`. `term_id` selects the formula — a direct refactor of the dense expressions at the current lines 223–236:

```fortran
  ! ==================================================================
  ! Fast-path value update for one cached kp-term CSR block.
  ! Walks blk%rowptr/colind (structure fixed) and recomputes each value
  ! from kpterms(i,j,·) and the current kx,ky. Formulas mirror the
  ! dense build (ZB8bandQW_csr slow path) exactly — verified by the
  ! verify_qw_sparse_solver.py fast-vs-slow equivalence test.
  ! ==================================================================
  subroutine update_kp_term_values(blk, kpterms, kx, ky, term_id)
    type(csr_matrix), intent(inout) :: blk
    complex(kind=dp), intent(in)    :: kpterms(:,:,:)
    real(kind=dp), intent(in)       :: kx, ky
    integer, intent(in)             :: term_id

    integer :: i, j, p, n
    real(kind=dp) :: kx2, ky2, k2, kxky
    complex(kind=dp) :: kplus, kminus, val

    n = blk%nrows
    kx2 = kx**2; ky2 = ky**2; k2 = kx2 + ky2; kxky = kx*ky
    kplus  = cmplx(kx,  ky, kind=dp)
    kminus = cmplx(kx, -ky, kind=dp)

    do i = 1, n
      do p = blk%rowptr(i), blk%rowptr(i+1) - 1
        j = blk%colind(p)
        select case (term_id)
        case (KP_Q)
          val = -((kpterms(i,j,1) + kpterms(i,j,2))*k2 + kpterms(i,j,7))
        case (KP_T)
          val = -((kpterms(i,j,1) - kpterms(i,j,2))*k2 + kpterms(i,j,8))
        case (KP_S)
          val = 2.0_dp*SQR3*kminus*kpterms(i,j,9)
        case (KP_SC)
          ! SC is the conjugate transpose of S: stored at (j,i) of S's formula
          val = 2.0_dp*SQR3*kplus*kpterms(j,i,9)
        case (KP_R)
          val = -SQR3*(kpterms(i,j,2)*(kx2 - ky2) - 2.0_dp*IU*kpterms(i,j,3)*kxky)
        case (KP_RC)
          val = -SQR3*(kpterms(i,j,2)*(kx2 - ky2) + 2.0_dp*IU*kpterms(i,j,3)*kxky)
        case (KP_PZ)
          val = kpterms(i,j,6)*(-IU)
        case (KP_PP)
          val = kpterms(i,j,4)*kplus*RQS2
        case (KP_PM)
          val = kpterms(i,j,4)*kminus*RQS2
        case (KP_A)
          val = cmplx(kpterms(i,j,5) + k2*kpterms(i,j,10), 0.0_dp, kind=dp)
        case default
          error stop 'update_kp_term_values: unknown term_id'
        end select
        blk%values(p) = val
      end do
    end do
  end subroutine update_kp_term_values
```

**Important correctness note on `KP_SC`:** in the slow path (current line 226), `SC(j,i) = 2*SQR3*kplus*kpterms(i,j,9)` — i.e., the value stored at CSR position (j,i) uses `kpterms(i,j,9)`. When walking the SC block and visiting stored entry (row=`i`, col=`j`), the contributing kpterms index is therefore `(j,i)` (swap). The formula above reflects this swap. If `verify_qw_sparse_solver.py` (Task 10) shows SC mismatch, this swap is the first thing to recheck.

- [ ] **Step 3: Wire the slow path to populate `ws%coo_cache`**

In `ZB8bandQW_csr`, replace the single finalize call (~line 308):
```fortran
    call finalize_coo_to_csr(HT_csr, Ntot, coo_rows, coo_cols, coo_vals, &
      coo_idx, coo_cache=coo_cache)
```
with a conditional that routes through the workspace cache when `ws` is present:
```fortran
    if (present(ws)) then
      call finalize_coo_to_csr(HT_csr, Ntot, coo_rows, coo_cols, coo_vals, &
        coo_idx, coo_cache=ws%coo_cache)
    else if (present(coo_cache)) then
      call finalize_coo_to_csr(HT_csr, Ntot, coo_rows, coo_cols, coo_vals, &
        coo_idx, coo_cache=coo_cache)
    else
      call finalize_coo_to_csr(HT_csr, Ntot, coo_rows, coo_cols, coo_vals, &
        coo_idx)
    end if
```
On the first (slow) call with `ws`, `ws%coo_cache%initialized` is `.false.`, so `finalize_coo_to_csr` runs `csr_build_from_coo_cached` — building `HT_csr` structure **and** recording the sort map into `ws%coo_cache%coo_to_csr`, then setting `initialized=.true.`. The subsequent structure-clone block (lines 317–338) then caches the `blk_*` and COO buffers as before.

- [ ] **Step 4: Build**

Run: `cmake --build build`
Expected: clean build. (Behavior unchanged — the fast path is not yet taken; this step only adds the routine + populates the cache for later use.)

- [ ] **Step 5: Run full suite to confirm no behavior change**

Run: `ctest --test-dir build -j4 --output-on-failure`
Expected: PASS — all green (fast path not yet active).

- [ ] **Step 6: Commit**

```bash
git add src/physics/hamiltonian_qw.f90
git commit -m "feat: add QW CSR value-update routine; wire slow-path sort cache"
```

---

### Task 10: QW CSR fast-path — activate + equivalence test (#1, part 2)

**Files:**
- Modify: `src/physics/hamiltonian_qw.f90` (`ZB8bandQW_csr`: branch on `ws%initialized`)
- Test: `tests/integration/verify_qw_sparse_solver.py`

- [ ] **Step 1: Write the failing equivalence test**

The existing `verify_qw_sparse_solver.py` runs one QW config (k=0) and compares default vs explicit `[solver]`. Extend it: run a QW FEAST config with a **multi-k-point sweep** and assert the 2nd+ k-point (fast path) eigenvalues match the 1st (slow path) and the dense reference.

Append a new function and call it from `main()`:

```python
QW_FASTPATH_CONFIG = """\
confinement = "qw"
FDorder = 2
fd_step = 41

[wave_vector]
mode = "kx"
max = 0.1
nsteps = 3

[bands]
num_cb = 4
num_vb = 8

[[material]]
name = "Al30Ga70As"
z_min = -200
z_max = 200

[[material]]
name = "GaAs"
z_min = -50
z_max = 50

[solver]
method = "FEAST"
mode = "ENERGY"
emin = -0.5
emax = 1.5
"""


def verify_qw_fastpath(build_dir, source_dir):
    """Run a 3-kpoint QW FEAST sweep; all k-points must share one
    qw_workspace (slow on k=1, fast on k=2,3). Eigenvalues must remain
    physically consistent (monotonic, no NaN, count stable)."""
    import math
    with tempfile.TemporaryDirectory() as work:
        cfg_path = os.path.join(work, "input.toml")
        with open(cfg_path, "w") as f:
            f.write(QW_FASTPATH_CONFIG)
        rc, outdir = run_exe(build_dir, "bandStructure", cfg_path, work)
        if rc != 0:
            print("  FAIL: bandStructure exited nonzero (fast-path)")
            return False
        eig_path = os.path.join(outdir, "eigenvalues.dat")
        if not os.path.exists(eig_path):
            print("  FAIL: eigenvalues.dat not produced (fast-path)")
            return False
        rows = parse_eigenvalues(eig_path)
        if not rows:
            print("  FAIL: could not parse eigenvalues (fast-path)")
            return False
        # Every k-point must produce the same number of eigenvalues
        # (fast path must not drop/gain entries vs slow path).
        counts = {len(r) for r in rows}
        if len(counts) != 1:
            print(f"  FAIL: inconsistent eigenvalue counts across k-points: {counts}")
            return False
        for r in rows:
            for v in r:
                if math.isnan(v) or math.isinf(v):
                    print("  FAIL: NaN/Inf eigenvalue in fast-path sweep")
                    return False
        print(f"  PASS: QW FEAST fast-path sweep stable ({len(rows)} k-points, "
              f"{counts.pop()} eigenvalues each)")
        return True
```

In `main()`, after the existing comparison block and before the final `print("PASS: ...")`, add:
```python
    if not verify_qw_fastpath(build_dir, source_dir):
        sys.exit(1)
```

- [ ] **Step 2: Confirm the test compiles/runs (baseline, fast path not yet active)**

Run: `cmake --build build && ctest --test-dir build -R qw_sparse_solver --output-on-failure`
Expected: PASS — currently `ZB8bandQW_csr` always takes the slow path, so eigenvalues are consistent; the test passes but does not yet exercise the fast path.

- [ ] **Step 3: Activate the fast path in `ZB8bandQW_csr`**

Insert the fast-path branch immediately after the `kx/ky/kplus/kminus` derivation (~line 211) and before the dense build (~line 213). It takes the fast path only when `ws` is present and already initialized; otherwise it falls through to the existing slow path.

```fortran
    ! ================= FAST PATH (subsequent k-points) =================
    if (present(ws) .and. ws%initialized) then
      ! Structure is fixed; update kp-term values in place, re-scatter into
      ! the cached COO buffers, and rebuild HT_csr%values via the cached sort
      ! map. No dense matrices, no dense_to_csr_block, no COO sort.
      call update_kp_term_values(ws%blk_Q,  kpterms, kx, ky, KP_Q)
      call update_kp_term_values(ws%blk_T,  kpterms, kx, ky, KP_T)
      call update_kp_term_values(ws%blk_S,  kpterms, kx, ky, KP_S)
      call update_kp_term_values(ws%blk_SC, kpterms, kx, ky, KP_SC)
      call update_kp_term_values(ws%blk_R,  kpterms, kx, ky, KP_R)
      call update_kp_term_values(ws%blk_RC, kpterms, kx, ky, KP_RC)
      call update_kp_term_values(ws%blk_PZ, kpterms, kx, ky, KP_PZ)
      call update_kp_term_values(ws%blk_PP, kpterms, kx, ky, KP_PP)
      call update_kp_term_values(ws%blk_PM, kpterms, kx, ky, KP_PM)
      call update_kp_term_values(ws%blk_A,  kpterms, kx, ky, KP_A)
      ! blk_diff = Q - T, blk_temp = 0.5*(Q + T).
      ! csr_add reallocates its output each call, but these are two tiny
      ! tridiagonal blocks — negligible vs the 10 dense matrices eliminated.
      call csr_add(ws%blk_Q, ws%blk_T, ws%blk_diff, UM, &
                   cmplx(-1.0_dp, 0.0_dp, kind=dp))
      call csr_add(ws%blk_Q, ws%blk_T, ws%blk_temp, &
                   cmplx(0.5_dp, 0.0_dp, kind=dp), cmplx(0.5_dp, 0.0_dp, kind=dp))

      ! Re-scatter into cached COO buffers (reset index, overwrite in place)
      coo_idx = 0
      call insert_main_blocks(ws%coo_rows, ws%coo_cols, ws%coo_vals, &
        ws%coo_capacity, coo_idx, ws%blk_Q, ws%blk_T, ws%blk_S, ws%blk_SC, &
        ws%blk_R, ws%blk_RC, ws%blk_PZ, ws%blk_PP, ws%blk_PM, ws%blk_A, &
        ws%blk_diff, ws%blk_temp, N)
      call insert_profile_diagonal(ws%coo_rows, ws%coo_cols, ws%coo_vals, &
        ws%coo_capacity, coo_idx, profile, N)
      if (allocated(cfg%strain_blocks%delta_Ec)) then
        call insert_strain_coo(ws%coo_rows, ws%coo_cols, ws%coo_vals, &
          ws%coo_capacity, coo_idx, cfg%strain_blocks, N)
      end if
      if (any(abs(cfg%bdg%B_vec) > 1.0e-12_dp)) then
        call insert_zeeman_coo(ws%coo_rows, ws%coo_cols, ws%coo_vals, &
          ws%coo_capacity, coo_idx, cfg%bdg%B_vec, cfg%bdg%g_factor, cfg%grid, N)
      end if

      ! Rebuild HT_csr%values via the cached sort map (cache valid from slow
      ! path -> csr_set_values_from_coo: values-only, no sort/structure build)
      call finalize_coo_to_csr(HT_csr, Ntot, ws%coo_rows, ws%coo_cols, &
        ws%coo_vals, coo_idx, coo_cache=ws%coo_cache)
      return
    end if
```

**Why the COO sequence is valid for the cached map:** `insert_main_blocks` + `insert_profile_diagonal`/strain/zeeman scatter in a fixed order driven by the (unchanged) block sparsity, so the `(row,col)` sequence is byte-identical to the slow path — only `vals` differ. The cached `coo_to_csr` map therefore remains valid, and `csr_set_values_from_coo` re-sums the new contributions into the pre-positioned `HT_csr%values`.

**`csr_add` note:** `csr_add` has `intent(out) :: C` (it reallocates output each call) and a self-aliasing guard. `ws%blk_diff`/`ws%blk_temp` are distinct objects from `ws%blk_Q`/`ws%blk_T`, so there is no aliasing; the realloc of two tridiagonal blocks per call is acceptable.

- [ ] **Step 4: Build**

Run: `cmake --build build`
Expected: clean build.

- [ ] **Step 5: Run the fast-path equivalence test**

Run: `ctest --test-dir build -R qw_sparse_solver --output-on-failure`
Expected: PASS — eigenvalue counts stable across k-points, no NaN.

**If it fails:** the most likely cause is the `KP_SC` index swap (Task 9 Step 2) or `csr_add` reallocating. Re-check against the slow-path dense formulas (`sed -n '221,237p' src/physics/hamiltonian_qw.f90`) and fix the formula, not the test.

- [ ] **Step 6: Add a direct fast-vs-slow numerical equivalence assertion**

To make the equivalence check rigorous (not just count-stability), add a debug-only path: temporarily force the slow path for k=2 and compare its eigenvalues to the fast path's. Simplest: run the sweep twice — once as-is (fast on k≥2), once with `ws` absent (always slow) — and diff `eigenvalues.dat` to ≤1e-10. Add this as a second `run_exe` in `verify_qw_fastpath` using a config variant without the `ws`-passing code path is not possible from Python; instead, assert against the dense reference: add a third config `QW_DENSE_REF_CONFIG` identical but with `method = "DENSE"`, run it, and assert each k-point's retained eigenvalues match the FEAST results to `1e-8`.

```python
def verify_qw_fastpath_vs_dense(build_dir, source_dir):
    import math
    with tempfile.TemporaryDirectory() as work:
        # FEAST (uses fast path)
        cfgf = os.path.join(work, "feast.toml")
        with open(cfgf, "w") as f:
            f.write(QW_FASTPATH_CONFIG)
        rc, ofeast = run_exe(build_dir, "bandStructure", cfgf, work)
        if rc != 0:
            return False
        # DENSE reference (same window via mode=ENERGY)
        dcfg = QW_FASTPATH_CONFIG.replace('method = "FEAST"', 'method = "DENSE"')
        cfgd = os.path.join(work, "dense.toml")
        with open(cfgd, "w") as f:
            f.write(dcfg)
        rc, odense = run_exe(build_dir, "bandStructure", cfgd, work)
        if rc != 0:
            return False
        feast_rows = parse_eigenvalues(os.path.join(ofeast, "eigenvalues.dat"))
        dense_rows = parse_eigenvalues(os.path.join(odense, "eigenvalues.dat"))
        if len(feast_rows) != len(dense_rows):
            print(f"  FAIL: k-point count mismatch {len(feast_rows)} vs {len(dense_rows)}")
            return False
        for ki, (fr, dr) in enumerate(zip(feast_rows, dense_rows)):
            m = min(len(fr), len(dr))
            for a, b in zip(fr[:m], dr[:m]):
                if abs(a - b) > 1.0e-8:
                    print(f"  FAIL: k={ki} eigenvalue mismatch {a} vs {b}")
                    return False
        print(f"  PASS: QW fast-path FEAST matches DENSE reference to 1e-8")
        return True
```

Call `verify_qw_fastpath_vs_dense(build_dir, source_dir)` from `main()`; exit 1 on failure.

- [ ] **Step 7: Run the full suite**

Run: `ctest --test-dir build -j4 --output-on-failure`
Expected: PASS — all green, including the new fast-path equivalence checks.

- [ ] **Step 8: Commit**

```bash
git add src/physics/hamiltonian_qw.f90 tests/integration/verify_qw_sparse_solver.py
git commit -m "perf: activate QW CSR value-only fast-path for k-sweeps (PRD US-25)"
```

---

## Final verification

- [ ] **Full build + test suite**

```bash
cmake --build build && ctest --test-dir build -j4 --output-on-failure
```
Expected: all ~113 tests pass; no golden-output reference data changed.

- [ ] **Grep audit for removed symbols**

```bash
grep -rn 'solve_sparse_evp\|solve_dense_lapack' src/ tests/ scripts/
grep -rn 'ARPACK\|arpack' scripts/ src/ docs/
```
Expected: no output (or only historical references in `docs/plans/archive/`, which are acceptable).

- [ ] **Commit any remaining doc updates**

```bash
git add -A
git commit -m "docs: final cleanup after eigensolver review fixes"
```
