# PR #39 Review Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Apply the corrective fixes from the PR #39 max-effort code review (Groups A+B+safe-C) so the opt-in FEAST paths are crash-free, never silently truncate a spectrum, and honor the bulk-never-FEAST and both-or-neither-window invariants — with golden AUTO/DENSE outputs unchanged.

**Architecture:** Seven localized fixes across `eigensolver.f90`, `defs.f90`, `hamiltonianConstructor.f90`, `main.f90`, `main_optics.f90`, `wire_setup.f90`. The largest is a shared `reconcile_band_slice` offset-helper that replaces three divergent hand-rolled FEAST→INDEX reconciliation blocks (the source of the OOB + wrong-bands bugs). Two `validate()` additions enforce spec invariants; one guard fix in `solve_feast` makes `converged`/`nev_found` consistent; three small cleanups remove dead code and latent crashes.

**Tech Stack:** Fortran 2018 (`-std=f2018`), gfortran, MKL (FEAST/PARDISO), pFUnit 4.x unit tests, ctest + Python verify scripts, OpenMP. Build: `cmake --build build`. Test: `OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -j4 --output-on-failure`.

**Spec:** `docs/superpowers/specs/2026-06-21-pr39-review-fixes-design.md` (commit `400e50d`).

**Plan-level refinement from spec:** the reconciliation helper returns the integer slice offset (`reconcile_band_slice(nev_found, il, iu, idx_lo)`), not array copies. Passing `eig(:,k)` / `eigv(:,:,k)` array-sections into `contiguous` dummies would create stack temporaries inside the OpenMP region — the known -O3+OpenMP corruption mode (`docs/solutions` / `project_segfault_fix`). Centralizing the offset decision (the divergent/buggy part) is the DRY win; the mechanical copy stays inline at each site.

**Branch:** `refactor/architecture-deepening` (PR #39's branch — fixes land on top of the 10-issue work, as follow-up commits).

---

## File Structure

| File | Responsibility | Tasks |
|------|----------------|-------|
| `src/core/defs.f90` | `validate()` — add A4 partial-window rejection + A3 bulk+FEAST rejection | 1, 2 |
| `src/math/eigensolver.f90` | `solve_feast` info=3 guard (A2); `reconcile_band_slice` helper (A1a); cache guard (B3); delete `asw_evals`+`asw_apply_margin` (C1) | 3, 4, 9, 10 |
| `src/apps/main_optics.f90` | QW optics k-sweep → use helper (A1b) | 5 |
| `src/apps/main.f90` | Landau k-sweep → helper (A1c); QW band-structure CSR sweep → helper (A1d) | 6, 7 |
| `src/physics/hamiltonianConstructor.f90` | `kp_scalar_block` error-stop (B2) | 8 |
| `src/physics/wire_setup.f90` | delete `wire_setup_adopt_precomputed` (C1) | 10 |
| `tests/unit/test_eigensolver.pf` | A2 info=3 test; A1a helper tests | 3, 4 |
| `tests/integration/test_validate_rejects_bad_configs.sh` | A3, A4 rejection cases | 1, 2 |
| `tests/integration/verify_qw_bandstructure_dense_feast.py` + `.sh` wrapper + `tests/CMakeLists.txt` | A1d DENSE-vs-FEAST QW band-structure regression | 7 |

**Conventions for every task:** after editing Fortran, rebuild with `cmake --build build` (mind the stale-`.mod` gotcha — `rm -f *.mod` only if you see type-mismatch errors). Commit each task with a Conventional Commit message; **no** Co-Authored-By / Generated-with trailer (attribution is blank globally).

---

### Task 1: A4 — reject partial `[solver]` energy windows in `validate()`

**Files:**
- Modify: `src/core/defs.f90` (inside `validate_simulation_config`, next to check I15 ~line 815-822)
- Modify: `tests/integration/test_validate_rejects_bad_configs.sh` (add a case)

- [ ] **Step 1: Write the failing rejection test**

Append this case inside `test_validate_rejects_bad_configs.sh` (after the last existing `run_test` block, before the final pass/fail summary):

```bash
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
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `ctest --test-dir build -R test_validate_rejects --output-on-failure`
Expected: FAIL — `Vnew1_partial_window_emin_only` reports exit=0 (bandStructure currently accepts the partial window) or the pattern isn't found.

- [ ] **Step 3: Implement the guard in `validate()`**

In `src/core/defs.f90`, inside `validate_simulation_config`, immediately before the existing I15 block (`! ---- I15: FEAST method cannot combine with INDEX mode ----`), add:

```fortran
      ! ---- I14b: a partial energy window is ambiguous ----
      ! 0 is the auto sentinel for BOTH emin and emax. Setting exactly one
      ! of them is ambiguous (is the 0 a literal bound or "auto"?), and the
      ! window authority's OR-override would otherwise produce a degenerate
      ! [emin, 0] window. Require both-or-neither. (Review finding #7.)
      if ((cfg%solver%emin == 0.0_dp) .neqv. (cfg%solver%emax == 0.0_dp)) then
        error stop 'validate_simulation_config: set both solver emin and emax, ' // &
          'or neither (0 = auto). A partial energy window is ambiguous.'
      end if
```

- [ ] **Step 4: Rebuild and run the rejection test to verify it passes**

Run: `cmake --build build && ctest --test-dir build -R test_validate_rejects --output-on-failure`
Expected: PASS — `Vnew1_partial_window_emin_only` matches "both solver emin and emax".

- [ ] **Step 5: Commit**

```bash
git add src/core/defs.f90 tests/integration/test_validate_rejects_bad_configs.sh
git commit -m "fix(validate): reject partial [solver] energy windows (#7)"
```

---

### Task 2: A3 — permanently reject `bulk + FEAST` in `validate()`

**Files:**
- Modify: `src/core/defs.f90` (inside `validate_simulation_config`, next to the I15 block)
- Modify: `tests/integration/test_validate_rejects_bad_configs.sh`

- [ ] **Step 1: Write the failing rejection test**

Append to `test_validate_rejects_bad_configs.sh`:

```bash
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
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `ctest --test-dir build -R test_validate_rejects --output-on-failure`
Expected: FAIL — `Vnew2_bulk_feast` exit=0 (currently accepted) or pattern absent.

- [ ] **Step 3: Implement the guard**

In `src/core/defs.f90`, inside `validate_simulation_config`, right after the I14b block added in Task 1 (still before I15), add:

```fortran
      ! ---- I14c: bulk never offers FEAST (permanent rejection) ----
      ! Bulk is always 8x8 dense by nature. AUTO resolves to DENSE for bulk,
      ! but an explicit method=FEAST must be rejected outright (PRD inv 1).
      ! Contrast FEAST+INDEX, rejected structurally at I15 below. (Review #5.)
      if (trim(cfg%confinement) == 'bulk' .and. trim(cfg%solver%method) == 'FEAST') then
        error stop 'validate_simulation_config: bulk is always 8x8 dense; FEAST ' // &
          'is not supported for bulk (use method = AUTO or DENSE).'
      end if
```

- [ ] **Step 4: Rebuild and run the rejection test**

Run: `cmake --build build && ctest --test-dir build -R test_validate_rejects --output-on-failure`
Expected: PASS for `Vnew2_bulk_feast` (and Task 1's case still passes).

- [ ] **Step 5: Commit**

```bash
git add src/core/defs.f90 tests/integration/test_validate_rejects_bad_configs.sh
git commit -m "fix(validate): permanently reject bulk + FEAST (#5)"
```

---

### Task 3: A2 — `solve_feast` must not populate the result on `info=3`

**Files:**
- Modify: `src/math/eigensolver.f90` (the result-population guard in `solve_feast`, ~line 424)
- Modify: `tests/unit/test_eigensolver.pf` (new `@test`)

- [ ] **Step 1: Write the failing unit test**

Append to `tests/unit/test_eigensolver.pf`, inside the `contains` block (after `test_feast_empty_matrix`):

```fortran
  ! ==================================================================
  ! Test: a FEAST solve whose energy window holds more eigenvalues than
  ! the subspace can hold (info=3, retries exhausted) must report
  ! converged=.false. AND nev_found=0 — never a truncated spectrum that a
  ! caller might silently use. (Review finding #3.)
  ! ==================================================================
  @test
  subroutine test_feast_subspace_saturated_yields_zero()
    integer, parameter :: n = 20
    real(kind=dp) :: diag(n)
    type(csr_matrix) :: H
    type(eigensolver_config) :: cfg
    type(eigensolver_result) :: res
    integer :: i

    ! Diagonal eigenvalues 1..20. Window [0.5, 20.5] holds all 20.
    do i = 1, n
      diag(i) = real(i, kind=dp)
    end do
    call build_diagonal_csr(diag, n, H)

    cfg%method = 'FEAST'
    cfg%emin  = 0.5_dp
    cfg%emax  = 20.5_dp
    cfg%nev   = 1
    cfg%m0    = 2           ! starts at 2; retries double to 4,8,16 — all < 20

    call solve_feast(H, cfg, res)

    @assertFalse(res%converged)
    @assertEqual(0, res%nev_found)   ! <-- the #3 fix: no truncated spectrum

    call eigensolver_result_free(res)
    call csr_free(H)
  end subroutine test_feast_subspace_saturated_yields_zero
```

- [ ] **Step 2: Build and run the test to verify it fails**

Run: `cmake --build build && ctest --test-dir build -R test_eigensolver --output-on-failure`
Expected: FAIL — `@assertEqual(0, res%nev_found)` fails because the current `info >= 0` guard populates `nev_found = 16` (M=16 after the last retry) even though `converged = .false.`.

- [ ] **Step 3: Tighten the guard in `solve_feast`**

In `src/math/eigensolver.f90`, in subroutine `solve_feast`, change the result-population guard (~line 424):

```fortran
    ! was: if (M > 0 .and. M <= M0 .and. info >= 0) then
    if (M > 0 .and. M <= M0 .and. (info == 0 .or. info == 2)) then
```

(Only genuine convergence — FEAST `info` 0 or 2 — populates the result. `info = 3` means the subspace was too small even after retries, so the returned M is incomplete; report `nev_found = 0` instead of a truncated spectrum callers might silently use. `result%converged` already reflects exactly `(info == 0 .or. info == 2)`, so the two are now consistent.)

- [ ] **Step 4: Rebuild and run the test to verify it passes**

Run: `cmake --build build && ctest --test-dir build -R test_eigensolver --output-on-failure`
Expected: PASS — `test_feast_subspace_saturated_yields_zero` passes, and all existing `test_feast_*` tests still pass (they are converged cases, unaffected by the tighter guard).

- [ ] **Step 5: Commit**

```bash
git add src/math/eigensolver.f90 tests/unit/test_eigensolver.pf
git commit -m "fix(eigensolver): solve_feast returns nev_found=0 on info=3 (#3)"
```

---

### Task 4: A1a — add the `reconcile_band_slice` offset helper

**Files:**
- Modify: `src/math/eigensolver.f90` (new public subroutine + export)
- Modify: `tests/unit/test_eigensolver.pf` (two `@test`s)

- [ ] **Step 1: Write the failing unit tests**

Append to `tests/unit/test_eigensolver.pf`:

```fortran
  ! ==================================================================
  ! Tests for reconcile_band_slice: the single source of the FEAST-full /
  ! DENSE-INDEX reconciliation offset. (Review findings #1, #2, #4, #10.)
  ! ==================================================================
  @test
  subroutine test_reconcile_band_slice_full_spectrum()
    integer :: idx_lo
    ! FEAST+ENERGY returned the full in-window spectrum (nev_found=20 >= iu=18):
    ! extract the global [il=15, iu=18] slice.
    call reconcile_band_slice(20, 15, 18, idx_lo)
    @assertEqual(15, idx_lo)
  end subroutine test_reconcile_band_slice_full_spectrum

  @test
  subroutine test_reconcile_band_slice_pre_sliced()
    integer :: idx_lo
    ! DENSE+INDEX already returned exactly [il, iu] (nev_found == nev_target = 4):
    ! take it verbatim from offset 1.
    call reconcile_band_slice(4, 15, 18, idx_lo)
    @assertEqual(1, idx_lo)
  end subroutine test_reconcile_band_slice_pre_sliced
```

- [ ] **Step 2: Build and run to verify they fail**

Run: `cmake --build build && ctest --test-dir build -R test_eigensolver --output-on-failure`
Expected: FAIL — compile error: `reconcile_band_slice` has no implicit interface / not defined.

- [ ] **Step 3: Implement the helper**

In `src/math/eigensolver.f90`:

(a) Add to the public list near the top of the module (next to the other `public ::` lines, ~line 12-23):

```fortran
  public :: reconcile_band_slice
```

(b) Add the subroutine in the `contains` block (e.g. right after `eigensolver_config_validate`):

```fortran
  ! ==================================================================
  ! Reconcile a solver result to the DENSE+INDEX [il, iu] band window and
  ! return the offset into result%eigenvalues/eigenvectors to copy from.
  ! Single source of the reconciliation decision (was copy-pasted — and had
  ! diverged into an OOB + a wrong-bands bug — across the three k-sweeps).
  !
  ! Two result shapes are handled:
  !   - FEAST+ENERGY full in-window spectrum (nev_found >= iu): the global
  !     [il, iu] slice lives at offset il.
  !   - DENSE+INDEX pre-sliced result (nev_found == iu-il+1): already the
  !     requested slice; offset 1.
  !   - partial (nev_target < nev_found < iu): the FEAST window truncated
  !     the spectrum — refuse rather than read out of bounds. (pFUnit 4.x
  !     cannot assert error-stop; this branch is defense-in-depth, verified
  !     by the call sites' nev_found < nev_target guards that fire first.)
  ! ==================================================================
  subroutine reconcile_band_slice(nev_found, il, iu, idx_lo)
    integer, intent(in)  :: nev_found, il, iu
    integer, intent(out) :: idx_lo
    integer :: nev_target
    character(len=160) :: msg

    nev_target = iu - il + 1

    if (nev_found >= iu) then
      idx_lo = il
    else if (nev_found == nev_target) then
      idx_lo = 1
    else
      write(msg, '(A,I0,A,I0,A,I0,A,I0,A)') &
        'solver returned ', nev_found, ' eigenvalues; need >= ', iu, &
        ' to extract global bands [', il, ',', iu, ']. Widen the FEAST energy window.'
      error stop 'reconcile_band_slice: FEAST window truncated — ' // trim(msg)
    end if
  end subroutine reconcile_band_slice
```

- [ ] **Step 4: Rebuild and run the tests to verify they pass**

Run: `cmake --build build && ctest --test-dir build -R test_eigensolver --output-on-failure`
Expected: PASS — both new tests pass; existing tests unaffected.

- [ ] **Step 5: Commit**

```bash
git add src/math/eigensolver.f90 tests/unit/test_eigensolver.pf
git commit -m "refactor(eigensolver): add reconcile_band_slice offset helper (#1/#2/#4/#10)"
```

---

### Task 5: A1b — route the QW optics k-sweep through the helper

**Files:**
- Modify: `src/apps/main_optics.f90` (the QW k_par sweep OpenMP block, ~lines 383-407)

- [ ] **Step 1: Read the current block**

Run: `sed -n '380,410p' src/apps/main_optics.f90`
Confirm the block looks like:
```fortran
          M_loc = qw_result%nev_found
          if (M_loc < nev_target) then
            ... error stop 'insufficient eigenvalues for QW optics accumulation'
          end if
          if (M_loc == nev_target) then
            idx_lo = 1
          else
            idx_lo = il
          end if
          eig(1:nev_target, k) = qw_result%eigenvalues(idx_lo:idx_lo+nev_target-1)
          eigv(:, 1:nev_target, k) = qw_result%eigenvectors(:, idx_lo:idx_lo+nev_target-1)
```

- [ ] **Step 2: Replace the offset decision with the helper call**

Replace the `if (M_loc == nev_target) ... else ... end if` block and the two copy lines with:

```fortran
          call reconcile_band_slice(qw_result%nev_found, il, iuu, idx_lo)
          eig(1:nev_target, k) = qw_result%eigenvalues(idx_lo:idx_lo+nev_target-1)
          eigv(:, 1:nev_target, k) = qw_result%eigenvectors(:, idx_lo:idx_lo+nev_target-1)
```

Keep the preceding `if (qw_result%nev_found < nev_target) error stop` guard (it fires first). Keep `idx_lo` in the OpenMP `private(...)` clause (it is still used). Drop the now-unused `M_loc` assignment only if no later line references it — grep first: `grep -n M_loc src/apps/main_optics.f90`. If `M_loc` is now unused, remove its declaration/assignment.

- [ ] **Step 3: Build and run the QW-optics FEAST regression**

Run: `cmake --build build && ctest --test-dir build -R "qw_optics" --output-on-failure`
Expected: PASS — `verify_qw_optics_dense_feast` still passes (the AUTO-window path returns the full spectrum, so the helper takes offset `il` exactly as before; behavior unchanged).

- [ ] **Step 4: Commit**

```bash
git add src/apps/main_optics.f90
git commit -m "refactor(optics): route QW k-sweep through reconcile_band_slice"
```

---

### Task 6: A1c — route the Landau k-sweep through the helper

**Files:**
- Modify: `src/apps/main.f90` (the Landau k-sweep OpenMP block, ~lines 529-556)

- [ ] **Step 1: Read the current block**

Run: `sed -n '528,558p' src/apps/main.f90`
Confirm the `M = landau_result%nev_found` → `if (M < nev_target) error stop` → `if (M == nev_target) idx_lo=1 else idx_lo=il` → copy block.

- [ ] **Step 2: Replace the offset decision with the helper call**

Replace the `if (M == nev_target) ... else ... end if` block + copy lines with:

```fortran
          call reconcile_band_slice(landau_result%nev_found, il, iuu, idx_lo)
          eig(1:nev_target, k) = landau_result%eigenvalues(idx_lo:idx_lo+nev_target-1)
          eigv(:, 1:nev_target, k) = landau_result%eigenvectors(:, idx_lo:idx_lo+nev_target-1)
```

Keep the preceding `if (landau_result%nev_found < nev_target) error stop` guard. Keep `idx_lo` in the `private(...)` clause.

- [ ] **Step 3: Build and run the Landau FEAST regression**

Run: `cmake --build build && ctest --test-dir build -R "landau" --output-on-failure`
Expected: PASS — `verify_landau_dense_feast` still passes (AUTO path unchanged; DENSE+INDEX gives `nev_found == nev_target` → offset 1).

- [ ] **Step 4: Commit**

```bash
git add src/apps/main.f90
git commit -m "refactor(main): route Landau k-sweep through reconcile_band_slice"
```

---

### Task 7: A1d — route the QW band-structure CSR sweep through the helper + lock with a verify test

**Files:**
- Modify: `src/apps/main.f90` (the QW FEAST CSR k-sweep block, ~lines 854-880)
- Create: `tests/integration/verify_qw_bandstructure_dense_feast.py`
- Create: `tests/integration/test_qw_bandstructure_dense_feast.sh`
- Modify: `tests/CMakeLists.txt`

**Note:** this site currently uses `M = min(nev_found, nev_target); eig(1:M) = eigenvalues(1:M)` (lowest-M) — the divergent/incorrect path. Switching to the helper changes opt-in-FEAST output from "deepest valence" to the correct gap-straddling `[il, iuu]` slice. AUTO→DENSE is untouched, so golden outputs don't move. The new verify test locks the corrected FEAST path.

- [ ] **Step 1: Write the failing verify test (mirrors `verify_qw_optics_dense_feast.py`)**

Create `tests/integration/verify_qw_bandstructure_dense_feast.py`. Copy the structure of `verify_qw_optics_dense_feast.py` (same `run_exe` / parse helpers from `star_helpers`, same DENSE-vs-FEAST two-run pattern), but compare **band-structure eigenvalues** (`output/eigenvalues.dat`) instead of optical spectra. The differing core logic:

```python
"""
Verify QW band structure: FEAST+ENERGY (sweep-envelope window) vs DENSE+INDEX.
FEAST-enable regression gate for the QW band-structure path (review #4).
Runs bandStructure twice (method=DENSE, method=FEAST) on a QW config and
asserts the gap-straddling eigenvalue bands agree within tolerance.
"""
# COVERAGE: observable=band_structure geometry=qw material=GaAs/AlGaAs backend=FEAST-vs-DENSE

import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from star_helpers import run_exe, TOL

EIG_TOL = 1.0e-8  # eigenvalues: FEAST vs LAPACK agree to ~machine precision

def parse_eigenvalues(path):
    """Parse output/eigenvalues.dat -> (k_points, evals[nk, n_bands]). Skip '#' header lines."""
    rows = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            rows.append([float(x) for x in s.split()])
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1:]

def main():
    cfg = os.path.join(os.path.dirname(__file__), '..', 'regression', 'configs',
                       'qw_gaas_algaas_bandstructure.toml')  # see Step 5
    # ... two runs: write cfg with [solver] method=DENSE, then method=FEAST,
    #     each via run_exe('bandStructure') into a separate cwd ...
    k_d, ev_d = parse_eigenvalues(os.path.join(out_dense, 'eigenvalues.dat'))
    k_f, ev_f = parse_eigenvalues(os.path.join(out_feast, 'eigenvalues.dat'))
    assert k_d.shape == k_f.shape, f"k-grid mismatch {k_d.shape} vs {k_f.shape}"
    # Compare the gap-straddling band window (same column count both runs)
    n = min(ev_d.shape[1], ev_f.shape[1])
    worst = float(np.max(np.abs(ev_d[:, :n] - ev_f[:, :n])))
    assert worst < EIG_TOL, f"QW band-structure DENSE vs FEAST worst diff {worst:.3e} > {EIG_TOL:.0e}"
    print(f"PASS: QW band-structure DENSE vs FEAST worst eigenvalue diff {worst:.3e}")
    # DISPATCH GATE: the FEAST run's stdout must contain the FEAST backend tag
    # (proves method=FEAST was honored, not silently overridden to DENSE).
    assert 'FEAST' in feast_stdout, "FEAST run did not report FEAST backend"
    print("PASS: dispatch gate (FEAST backend honored)")

if __name__ == '__main__':
    main()
```

Fill in the two-run scaffolding by copying the equivalent section from `verify_qw_optics_dense_feast.py` (the `run_exe` calls into per-backend `cwd`s, writing `[solver] method = ...` into each `input.toml`).

Create the shell wrapper `tests/integration/test_qw_bandstructure_dense_feast.sh`:
```bash
#!/bin/bash
# COVERAGE: observable=band_structure geometry=qw material=GaAs/AlGaAs
set -euo pipefail
EXE="$1"
cd "$(dirname "$0")"
python3 verify_qw_bandstructure_dense_feast.py "$EXE"
```
`chmod +x tests/integration/test_qw_bandstructure_dense_feast.sh`.

- [ ] **Step 2: Create the QW band-structure config**

Create `tests/regression/configs/qw_gaas_algaas_bandstructure.toml` — a standard GaAs/AlGaAs QW with an in-plane k-sweep. Copy an existing QW config from `tests/regression/configs/` (e.g. one used by `verify_8band_rung3_qw.py`), and ensure it has `[wave_vector]` (mode `kxky`, `nsteps` odd ≥ 5), `[bands]` (num_cb/num_vb), `[[material]]` GaAs well + AlGaAs barriers, and **no** `[solver]` section (the verify script injects `method`). Keep `fd_step` modest (e.g. 60) so FEAST is fast.

- [ ] **Step 3: Register the test in CMake**

In `tests/CMakeLists.txt`, add (mirroring the optics entry — find it with `grep -n qw_optics tests/CMakeLists.txt` and copy the block):

```cmake
add_test(NAME verify_qw_bandstructure_dense_feast
         COMMAND ${bash} ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_qw_bandstructure_dense_feast.sh
                 $<TARGET_FILE:bandStructure>)
set_tests_properties(verify_qw_bandstructure_dense_feast PROPERTIES LABELS "regression")
```

- [ ] **Step 4: Build and run the verify test to verify it fails**

Run: `cmake --build build && ctest --test-dir build -R verify_qw_bandstructure_dense_feast --output-on-failure`
Expected: FAIL — the QW FEAST CSR path still uses lowest-M (`eig(1:M)=eigenvalues(1:M)`), so FEAST returns the deepest valence bands while DENSE returns the `[il,iuu]` gap-straddling bands → worst diff >> tolerance.

- [ ] **Step 5: Route the sweep through the helper**

In `src/apps/main.f90`, the QW FEAST CSR k-sweep block (~lines 854-880). Replace:
```fortran
            M = min(result_bs%nev_found, iuu-il+1)
            if (result_bs%nev_found > iuu-il+1) then
              ... 'Warning: ... returned ... only the lowest will be kept ...'
            end if
            if (result_bs%nev_found < (iuu - il + 1)) then
              ... 'Warning: ... returned only ... missing bands zero-filled ...'
            end if
            eig(1:M, k) = result_bs%eigenvalues(1:M)
            eigv(:, 1:M, k) = result_bs%eigenvectors(:, 1:M)
```
with:
```fortran
            nev_target = iuu - il + 1
            if (result_bs%nev_found < nev_target) then
              print '(A,A,A,I0,A,I0,A)', ' ERROR: ', solver_bs%backend_name(), &
                ' returned only ', result_bs%nev_found, ' eigenvalues at k-point ', k, &
                '; widen the energy window'
              error stop 'insufficient eigenvalues for QW band-structure k-sweep'
            end if
            call reconcile_band_slice(result_bs%nev_found, il, iuu, idx_lo)
            eig(1:nev_target, k) = result_bs%eigenvalues(idx_lo:idx_lo+nev_target-1)
            eigv(:, 1:nev_target, k) = result_bs%eigenvectors(:, idx_lo:idx_lo+nev_target-1)
```
(Ensure `idx_lo` and `nev_target` are declared in this scope — `idx_lo` likely already is; add `integer :: nev_target, idx_lo` if not. `M` may become unused — remove it if so.)

- [ ] **Step 6: Build and run the verify test to verify it passes**

Run: `cmake --build build && ctest --test-dir build -R verify_qw_bandstructure_dense_feast --output-on-failure`
Expected: PASS — DENSE and FEAST now return the same gap-straddling bands.

- [ ] **Step 7: Commit**

```bash
git add src/apps/main.f90 tests/integration/verify_qw_bandstructure_dense_feast.py \
        tests/integration/test_qw_bandstructure_dense_feast.sh \
        tests/regression/configs/qw_gaas_algaas_bandstructure.toml tests/CMakeLists.txt
git commit -m "fix(main): QW band-structure FEAST returns gap-straddling bands via helper (#4)"
```

---

### Task 8: B2 — `kp_scalar_block` error-stops on an unknown tag

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90` (`kp_scalar_block`, ~lines 777-796)

- [ ] **Step 1: Read the current function**

Run: `sed -n '777,797p' src/physics/hamiltonianConstructor.f90`
Confirm the `pure function kp_scalar_block(...)` with `case default; val = cmplx(0.0_dp, 0.0_dp, kind=dp)`.

- [ ] **Step 2: Make the failure altitude match `kp_block_dense`**

Change the function — drop `pure` (required: `error stop` is not allowed in a `pure` function) and make the `case default` an `error stop`, matching `kp_block_dense` (~line 726-727):

```fortran
  function kp_scalar_block(tag, Q, T, S, SC, R, RC, PP, PM, PZ, A) result(val)
    integer, intent(in) :: tag
    complex(kind=dp), intent(in) :: Q, T, S, SC, R, RC, PP, PM, PZ, A
    complex(kind=dp) :: val

    select case (tag)
    case (KP_Q);  val = Q
    case (KP_T);  val = T
    case (KP_S);  val = S
    case (KP_SC); val = SC
    case (KP_R);  val = R
    case (KP_RC); val = RC
    case (KP_PP); val = PP
    case (KP_PM); val = PM
    case (KP_PZ); val = PZ
    case (KP_A);  val = A
    case default
      error stop 'kp_scalar_block: unknown base kp_term tag'
    end select
  end function kp_scalar_block
```

- [ ] **Step 3: Build and run the bulk/kp-block regression**

Run: `cmake --build build && ctest --test-dir build -R "bulk|confinement_assembly|eigensolver" --output-on-failure`
Expected: PASS — no current code path passes an invalid tag, so behavior is unchanged for valid inputs; the change only hardens the failure mode.

- [ ] **Step 4: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90
git commit -m "fix(blocks): kp_scalar_block error-stops on unknown tag (#13)"
```

---

### Task 9: B3 — dense-solver cache must be correct under shrinking N

**Files:**
- Modify: `src/math/eigensolver.f90` (`dense_solve_dense_dispatch`, ~line 846)

- [ ] **Step 1: Read the current cache logic**

Run: `sed -n '844,854p' src/math/eigensolver.f90`
Confirm: `if (N > self%cached_n) then ... deallocate/reallocate ... self%cached_n = N`, followed by `self%A_buf = H`.

- [ ] **Step 2: Reallocate on any N change, not only growth**

Change the condition so the buffers are exact-sized for the current N (the hot same-N path still skips the reallocation):

```fortran
    ! (Re)allocate N-dependent buffers whenever N changes (grow OR shrink).
    ! The previous grow-only policy left oversized buffers and then did a
    ! non-conformable A_buf = H assignment if a solver was ever reused at a
    ! smaller N. Same-N (the hot sweep path) still reuses the buffers.
    if (N /= self%cached_n) then
      if (allocated(self%A_buf)) then
        deallocate(self%A_buf, self%Z_buf, self%W_buf, self%rwork, self%iwork, self%ifail)
      end if
      allocate(self%A_buf(N,N), self%Z_buf(N,N), self%W_buf(N))
      allocate(self%rwork(7*N), self%iwork(5*N), self%ifail(N))
      self%cached_n = N
    end if
```

- [ ] **Step 3: Build and run the dense-solver regression**

Run: `cmake --build build && ctest --test-dir build -R "eigensolver|landau|qw_optics" --output-on-failure`
Expected: PASS — same-N reuse is unchanged; the fix only adds correctness for the (currently unreached) shrinking-N case.

- [ ] **Step 4: Commit**

```bash
git add src/math/eigensolver.f90
git commit -m "fix(eigensolver): reallocate dense cache on any N change (#8)"
```

---

### Task 10: C1 — remove dead code (`adopt_precomputed`, `asw_evals`, `asw_apply_margin`)

**Files:**
- Modify: `src/physics/wire_setup.f90` (delete `wire_setup_adopt_precomputed` + export)
- Modify: `src/math/eigensolver.f90` (delete `asw_evals`, `asw_apply_margin`, the interface entry)
- Verify no caller exists before each deletion.

- [ ] **Step 1: Confirm zero callers (safety gate)**

Run both; each must return only definition/export lines, NO call site:
```bash
grep -rn "wire_setup_adopt_precomputed" src/ tests/
grep -rn "asw_evals" src/ tests/
```
If any caller appears, STOP and surface it — do not delete.

- [ ] **Step 2: Delete `wire_setup_adopt_precomputed`**

In `src/physics/wire_setup.f90`: remove the `public :: wire_setup_adopt_precomputed` from the public list (line ~50), and delete the entire `subroutine wire_setup_adopt_precomputed(...)` body (~lines 118-159). The header comment block describing the two entry points (lines ~21-28) should be trimmed to mention only `wire_setup_init`.

- [ ] **Step 3: Delete `asw_evals` + `asw_apply_margin` (also removes the duplicated margin constants — #11)**

In `src/math/eigensolver.f90`:
- Remove `module procedure :: asw_evals` from the `apply_solver_window` interface block (~line 58). The interface keeps `asw_envelope` and `asw_single`.
- Delete the `asw_evals` subroutine (~lines 574-593) and the private `asw_apply_margin` helper (~lines 503-511). The `margin_frac` / `margin_floor` constants lived only in `asw_apply_margin`, so they vanish with it (resolves #11 by deletion — DRY, single remaining user is `auto_compute_energy_window`).
- Trim the "Variant 3 — spectral (eigenvalue-array) source" comment block (~lines 566-593) since the variant is gone.

- [ ] **Step 4: Build and run the full unit + dispatch regression**

Run: `cmake --build build && ctest --test-dir build -L unit --output-on-failure && ctest --test-dir build -R "eigensolver|solver_derivation" --output-on-failure`
Expected: PASS — compiles (no dangling references); `test_eigensolver` / `test_solver_derivation` still pass (they never used the removed APIs).

- [ ] **Step 5: Commit**

```bash
git add src/physics/wire_setup.f90 src/math/eigensolver.f90
git commit -m "refactor: remove dead code (adopt_precomputed, asw_evals) + dup margin (#14/#11)"
```

---

### Task 11: Full regression gate

- [ ] **Step 1: Configure with tests (if not already)**

```bash
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl \
    -DBUILD_TESTING=ON -DPFUNIT_DIR=$HOME/.local/pfunit/PFUNIT-4.16/cmake
cmake --build build
```

- [ ] **Step 2: Run the full suite, OMP-capped**

Run: `OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -j4 --output-on-failure`
Expected: ALL PASS — the prior 123 tests plus the new ones (`Vnew1`, `Vnew2`, `test_feast_subspace_saturated_yields_zero`, `test_reconcile_band_slice_*`, `verify_qw_bandstructure_dense_feast`), 0 failures.

- [ ] **Step 3: Confirm golden regression outputs are unchanged**

Run: `git diff --stat tests/regression/data/`
Expected: empty — no golden output moved (all AUTO/DENSE default paths are untouched).

- [ ] **Step 4: Commit any CMake/test-registration tidy, then report**

```bash
git add -A
git status   # confirm only intended test/config files
git commit -m "test: register verify_qw_bandstructure_dense_feast gate"  # if not already committed in Task 7
```

---

## Self-Review (completed)

**1. Spec coverage** — A1→Tasks 4-7; A2→Task 3; A3→Task 2; A4→Task 1; B2→Task 8; B3→Task 9; C1(#14+#11)→Task 10; full gate→Task 11. Deferred #9/#12/#15 are in BACKLOG, intentionally unimplemented. ✓
**2. Placeholder scan** — Every code step shows the actual Fortran/Python/bash. The verify-script two-run scaffolding intentionally reuses `verify_qw_optics_dense_feast.py`'s `run_exe` section (named explicitly, with the differing parse/compare logic shown in full). No "TBD"/"add error handling". ✓
**3. Type consistency** — `reconcile_band_slice(nev_found, il, iu, idx_lo)` signature is identical in Task 4 (definition) and Tasks 5/6/7 (call sites). `nev_target`, `idx_lo`, `iuu`/`il` names match across sites. ✓
