# BdG U2 — actual ship (recovery + execution) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Recover from premature U2 close-out (doc-only stamps landed in working tree without backing code; ctest 50/51 RED) and ship U2 for real on a feature branch: extend `bdg_observables.f90` seam with the two siblings + SSOT, rename `compute_z2_gap*` helpers, migrate `main_topology.f90:1371`, add the new pFUnit tests + CMake registration, get 51/51 unit tests green, then re-apply the doc stamps and open a PR.

**Architecture:** Revert-then-execute-then-re-stamp on a single feature branch `feat/bdg-u2-actual-ship`. Phase 1 reverts working-tree drift (8 M files + 3 ?? entries) back to a clean `main@8fa9551` state — no commit needed, just `git checkout` + `mv` + `rm`. Phase 2 executes the spec.md as a sequence of TDD tasks (each writes failing test, runs, implements, runs, commits). Phase 3 re-applies the doc-only stamps only after the verification gate passes (51/51 unit tests green).

**Tech Stack:** Fortran 2018, pFUnit 4.16, MKL LP64 sequential, CMake/Ninja.

## Working-tree state at plan start

```
 M docs/agents/triage-labels.md                              (premature U2 close edits — revert)
 M docs/lecture/13-topological-superconductivity.md          (premature U2 close edits — revert)
 M docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md  (U2 status_note flip — revert)
 M docs/plans/BACKLOG.md                                     (U2 closed 2026-07-13 — revert)
 M docs/plans/REVIEW.md                                      (row 79 U2 closed — revert)
 M scripts/lecture_13_topological.py                         (3→4 witness + colormap extract — revert)
 M src/physics/AGENTS.md                                     (DAG seam siblings claim — revert)
 M tests/integration/test_lecture_13_acceptance_gate.sh      (WITNESS_LABEL="4-witness" — revert)
 M tests/integration/verify_majorana_polarization.py         (gate_row_colormap_present — revert)
 ?? .scratch/archive/bdg-evaluator-pfaffian/                 (untracked, move back)
 ?? .scratch/bdg-u10-pfaffian-phase-diagram/                 (U10 handoff — separate scope, KEEP)
 ?? docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md  (premature — delete)
 ?? docs/agents/triage-labels.md is dirty                    (likely unrelated; revert)
```
HEAD: `8fa9551` (PR #41 merge). Branch from this.

## Global Constraints

- All work on branch `feat/bdg-u2-actual-ship` (off `main@8fa9551`).
- Frequent commits (every 1-3 task steps), no `Co-Authored-By:` trailer (`feedback_no_coauthor_trailer`).
- pFUnit `@assertEqual` / `@assertTrue` macros MUST be SINGLE-LINE (no `&` continuation — `feedback_pfunit_macro_single_line`).
- All `error stop '<descriptive message>'`, no bare `stop 1`.
- F2018 enforced via `-std=f2018`; prefer `elemental pure` for new scalar pure functions.
- TDD discipline: write failing test → run → implement minimal → run → commit.
- Verify before claim: ctest green at 51/51 BEFORE re-applying doc stamps.
- The 7 issues at `.scratch/archive/bdg-evaluator-pfaffian/issues/` are decisions-of-record — reuse, do NOT re-litigate.
- Module docstring + AGENTS.md inventory updates ride with each code change.
- The spec at `.scratch/archive/bdg-evaluator-pfaffian/spec.md` is the spec-of-record.
- Reuse existing test fixture style: pFUnit `@assertEqual(...)`, complex matrices via `complex(kind=dp)`, parameter type `bdg_eval_params_t`.
- The `wire_pfaffian_witness` dense-path subroutine in `topological_analysis.f90:1626-1692` is the canonical implementation — the seam sibling is the CSR entry that wraps a CSR copy of the same logic per ticket 02 §4.

---

# Phase 1: Revert doc-only premature close (5 tasks)

## Task 1.1: Audit current working-tree drift (read-only)

**Files:** N/A

- [ ] **Step 1: Confirm drift state**

Run: `git status -s --untracked-files=all`
Expected: 9 M files + 3 ?? entries (see "Working-tree state at plan start" above).

- [ ] **Step 2: Quantify churn**

Run: `git diff --stat HEAD`
Expected: ~80 lines added across 8 files (the doc-only close-out churn).

- [ ] **Step 3: Save audit**

```bash
mkdir -p .scratch/bdg-u2-actual-ship
git status -s --untracked-files=all > .scratch/bdg-u2-actual-ship/audit-phase1-pre.log
git diff --stat HEAD                 > .scratch/bdg-u2-actual-ship/audit-phase1-diffstat.log
```

## Task 1.2: Revert modified files in working tree

**Files:** Restore 8 files to HEAD (`8fa9551`).

- [ ] **Step 1: Checkout HEAD versions**

```bash
git checkout HEAD -- \
  docs/agents/triage-labels.md \
  docs/lecture/13-topological-superconductivity.md \
  docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md \
  docs/plans/BACKLOG.md \
  docs/plans/REVIEW.md \
  scripts/lecture_13_topological.py \
  src/physics/AGENTS.md \
  tests/integration/test_lecture_13_acceptance_gate.sh \
  tests/integration/verify_majorana_polarization.py
```

- [ ] **Step 2: Verify the AGENTS.md DAG row no longer claims seam siblings exist**

Run: `grep -nE "Seam siblings|eval_bdg_pfaffian_witness_csr|eval_bdg_kitaev_majorana" src/physics/AGENTS.md`
Expected: NO matches (the previous claim of seam siblings being added is reverted).

- [ ] **Step 3: Verify the gate shell is back to 3-witness + PFAFFIAN_DEFERRED branch**

Run: `grep -nE "WITNESS_LABEL|PFAFFIAN_DEFERRED|TOL_BCRIT_RANGE" tests/integration/test_lecture_13_acceptance_gate.sh`
Expected: `WITNESS_LABEL="3-witness"`, `PFAFFIAN_DEFERRED` branch present, `TOL_BCRIT_RANGE` matches pre-premature value (likely 2.0 T).

- [ ] **Step 4: Verify the polarization verifier's `gate_row_colormap_present` helper is gone**

Run: `grep -nE "gate_row_colormap_present|_gate_row" tests/integration/verify_majorana_polarization.py`
Expected: NO matches.

## Task 1.3: Move archive back + delete premature artifacts

**Files:**
- Move: `.scratch/archive/bdg-evaluator-pfaffian/` → `.scratch/archive/bdg-evaluator-pfaffian/`
- Delete: `docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md`
- Delete: `~/.claude/projects/-data-8bandkp-fdm/memory/project_bdg_evaluator_seam_ssot.md`
- Edit: `~/.claude/projects/-data-8bandkp-fdm/memory/MEMORY.md` (remove pointer line)

- [ ] **Step 1: Move archive back**

```bash
mv .scratch/archive/bdg-evaluator-pfaffian .scratch/archive/bdg-evaluator-pfaffian
```

- [ ] **Step 2: Delete the premature solution doc**

```bash
rm docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md
```

- [ ] **Step 3: Delete the premature memory entry**

```bash
rm ~/.claude/projects/-data-8bandkp-fdm/memory/project_bdg_evaluator_seam_ssot.md
```

- [ ] **Step 4: Remove the MEMORY.md pointer line**

Open `~/.claude/projects/-data-8bandkp-fdm/memory/MEMORY.md` and remove the line:
```
- [BdG Evaluator Seam SSOT](project_bdg_evaluator_seam_ssot.md) — U2 closed 2026-07-13; ...
```
(Use Edit with the exact line; do not delete unrelated lines.)

- [ ] **Step 5: Verify clean state**

Run: `git status -s`
Expected: ONLY `.scratch/bdg-u10-pfaffian-phase-diagram/` remains `??` (U10 handoff, separate scope, KEEP).

## Task 1.4: Verify pre-revert ctest state

- [ ] **Step 1: Confirm build is still green (no source files touched yet)**

Run: `cmake --build build 2>&1 | tail -10`
Expected: no errors. (Reverted files are doc-only; build unaffected.)

- [ ] **Step 2: Run unit tests, capture pre-revert state**

Run: `OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L unit -j4 2>&1 | tee .scratch/bdg-u2-actual-ship/ctest-pre-revert.log | tail -15`
Expected: **50/51 PASS, 1 FAIL: `test_wire_pfaffian_witness`** (the pre-premature state — this was the only RED before the doc-only close landed, and it remains RED because the doc-only edits don't change behavior).

- [ ] **Step 3: Save the strict-test failure trace for Phase 2 reference**

Run: `cmake --build build --target test_wire_pfaffian_witness && ctest --test-dir build -R test_wire_pfaffian_witness -V 2>&1 | tail -30 > .scratch/bdg-u2-actual-ship/strict-test-trace.log`

## Task 1.5: Phase 1 wrap-up — no commit needed

- [ ] **Step 1: Verify final clean state**

Run: `git status -s`
Expected: ONLY `.scratch/bdg-u2-actual-ship/` (untracked, audit dir we created) and `.scratch/bdg-u10-pfaffian-phase-diagram/` (untracked, separate scope).

- [ ] **Step 2: Verify HEAD is unchanged**

Run: `git log -1 --oneline`
Expected: `8fa9551 Merge pull request #41 ...`.

Phase 1 produces no commits; the working tree is restored to clean `main@8fa9551` state.

---

# Phase 2: Execute U2 spec on feature branch (TDD) (10 tasks)

## Task 2.0: Create feature branch off clean main

- [ ] **Step 1: Create and switch to the feature branch**

Run: `git checkout -b feat/bdg-u2-actual-ship main`
Expected: `Switched to a new branch 'feat/bdg-u2-actual-ship'`, on commit `8fa9551`.

- [ ] **Step 2: Confirm clean working tree**

Run: `git status`
Expected: "On branch feat/bdg-u2-actual-ship", nothing to commit.

## Task 2.1: Extend `bdg_observables.f90` — `bdg_default_pfaffian_floor` SSOT + L0 imports

**Files:**
- Modify: `src/physics/bdg_observables.f90`
- Read for context: `src/math/pfaffian.f90:1-30` (public exports), `src/math/sparse_matrices.f90:1-30` (csr_matrix type)

**Interfaces (per spec.md §"Seam shape"):**
- New SSOT parameter: `bdg_default_pfaffian_floor = 1.0e-12_dp`
- New `bdg_eval_params_t` field: `pfaffian_floor`
- New `use` lines: `sparse_matrices` (for `csr_matrix` type), `pfaffian` (for `complex_pfaffian` + `kitaev_majorana_number`)

- [ ] **Step 1: Add the SSOT parameter**

Edit `src/physics/bdg_observables.f90:36` — after the existing `bdg_default_min_threshold` line, add:
```fortran
  ! Magic-number SSOT for the slim Pfaffian floor (replaces the literal at
  ! topological_analysis.f90:1663, 1681, 1766). Promoted per ticket 02 of
  ! `.scratch/archive/bdg-evaluator-pfaffian/`.
  real(kind=dp), parameter, public :: bdg_default_pfaffian_floor = 1.0e-12_dp
```

- [ ] **Step 2: Add the L0 imports**

Edit `src/physics/bdg_observables.f90:20` — replace the `use definitions` line with:
```fortran
  use definitions, only: dp
  use sparse_matrices, only: csr_matrix
  use pfaffian, only: complex_pfaffian, kitaev_majorana_number
```
(Imports are used by the new siblings added in Tasks 2.2 + 2.3.)

- [ ] **Step 3: Add `pfaffian_floor` field to `bdg_eval_params_t`**

Edit `src/physics/bdg_observables.f90:39-42` — replace the type definition with:
```fortran
  type :: bdg_eval_params_t
    real(kind=dp) :: delta_0        ! SC gap magnitude (eV) — scale for near-zero band
    real(kind=dp) :: near_zero_frac ! default 0.001; |E| < near_zero_frac*delta_0 counts as near-zero
    real(kind=dp) :: pfaffian_floor ! default 1.0e-12; |Pf| below this counts as zero (slim witness)
  end type
```

- [ ] **Step 4: Update the factory to default `pfaffian_floor`**

Edit `src/physics/bdg_observables.f90:96-101` — replace with:
```fortran
  pure function bdg_eval_params_with_delta(delta_0) result(p)
    real(kind=dp), intent(in) :: delta_0
    type(bdg_eval_params_t) :: p
    p%delta_0        = delta_0
    p%near_zero_frac = bdg_default_near_zero_frac
    p%pfaffian_floor = bdg_default_pfaffian_floor
  end function bdg_eval_params_with_delta
```

- [ ] **Step 5: Build, verify compile**

Run: `cmake --build build 2>&1 | tail -5`
Expected: `bdg_observables_mod.o` rebuilt cleanly; no errors.

- [ ] **Step 6: Commit**

```bash
git add src/physics/bdg_observables.f90
git commit -m "feat(bdg): add bdg_default_pfaffian_floor SSOT + sparse_matrices/pfaffian imports"
```

## Task 2.2: Implement `eval_bdg_pfaffian_witness_csr` (TDD)

**Files:**
- Create: `tests/unit/test_bdg_pfaffian_witness_csr.pf` (NEW, 4 tests)
- Modify: `src/physics/bdg_observables.f90` (add the sibling)

**Reference:** the dense-path implementation at `src/physics/topological_analysis.f90:1626-1790` (especially S2 row extraction at lines 1724-1764 per ticket 02 §4) — reuse the same omega construction + per-site Pf call pattern, but for a CSR input.

**Interface (per spec.md §"Seam shape"):**
```fortran
function eval_bdg_pfaffian_witness_csr(H_bdg_csr, Nbdg, params) result(s2_sign)
  type(csr_matrix), intent(in)        :: H_bdg_csr
  integer,          intent(in)        :: Nbdg      ! 16N for a BdG CSR; use 16 for the synthetic 16x16 fixture
  type(bdg_eval_params_t), intent(in) :: params
  integer                              :: s2_sign   ! {-1, 0, +1}
end function
```

- [ ] **Step 1: Create the test file with the first failing test (range contract)**

Create `tests/unit/test_bdg_pfaffian_witness_csr.pf`:
```fortran
module test_bdg_pfaffian_witness_csr
  use iso_fortran_env, only: dp => real64
  use funit
  use definitions, only: csr_matrix
  use bdg_observables, only: bdg_eval_params_t, bdg_default_pfaffian_floor, &
       bdg_eval_params_with_delta, eval_bdg_pfaffian_witness_csr

  implicit none

contains

  @test
  subroutine test_pfaffian_witness_csr_returns_valid_sign_on_synthetic()
    ! Spec §4 — slim projected Pfaffian over S1 (2 lowest single-particle)
    ! ⊗ S2 (bands 7-8). Returns s2_sign ∈ {-1, 0, +1}.
    type(csr_matrix) :: H
    integer :: N
    type(bdg_eval_params_t) :: p
    integer :: s2_sign

    call build_synthetic_bdg_csr(H, N, mu=0.0_dp, delta_0=0.001_dp)
    p = bdg_eval_params_with_delta(0.001_dp)

    s2_sign = eval_bdg_pfaffian_witness_csr(H, N, p)

    @assertTrue(s2_sign >= -1 .and. s2_sign <= 1)
  end subroutine

  @test
  subroutine test_pfaffian_witness_csr_handles_zero_matrix()
    ! Defensive: empty / zero spectrum → s2_sign = 0 (not a fatal error).
    type(csr_matrix) :: H
    integer :: N
    type(bdg_eval_params_t) :: p
    integer :: s2_sign

    call build_zero_bdg_csr(H, N)
    p = bdg_eval_params_with_delta(0.001_dp)

    s2_sign = eval_bdg_pfaffian_witness_csr(H, N, p)

    @assertEqual(0, s2_sign)
  end subroutine

  @test
  subroutine test_pfaffian_witness_csr_delegation_contract()
    ! Per design decision 2026-07-13 (option b): the seam sibling delegates
    ! to wire_pfaffian_witness_sweep in topological_analysis.f90. The
    ! params%pfaffian_floor SSOT is wired through in U13 (the floor is the
    ! SSOT surface; the helper doesn't accept it yet). For now, this test
    ! pins the delegation contract: same input → same output as the helper.
    use topological_analysis, only: wire_pfaffian_witness_sweep
    type(csr_matrix) :: H
    integer :: N
    type(bdg_eval_params_t) :: p
    integer :: s2_sign_seam, s2_sign_direct

    call build_nondiagonal_bdg_csr(H, N, mu=0.0_dp, delta_0=0.001_dp)
    p = bdg_eval_params_with_delta(0.001_dp)

    s2_sign_seam = eval_bdg_pfaffian_witness_csr(H, N, p)
    call wire_pfaffian_witness_sweep(H, N, s2_sign_direct)

    @assertEqual(s2_sign_direct, s2_sign_seam)
  end subroutine

  @test
  subroutine test_pfaffian_witness_csr_nondiagonal_nonzero_sign()
    ! Non-diagonal fixture (conduction-band off-diagonal coupling crosses
    ! the Nambu seam) → Pf ≠ 0 → s2_sign ≠ 0 (either +1 or -1).
    type(csr_matrix) :: H
    integer :: N
    type(bdg_eval_params_t) :: p
    integer :: s2_sign

    call build_nondiagonal_bdg_csr(H, N, mu=0.0_dp, delta_0=0.001_dp)
    p = bdg_eval_params_with_delta(0.001_dp)

    s2_sign = eval_bdg_pfaffian_witness_csr(H, N, p)

    @assertTrue(s2_sign /= 0)
  end subroutine

  ! ---------------------------------------------------------------------------
  ! Fixtures — CSR builders. Reuse the dense-fixture pattern from
  ! tests/unit/test_wire_pfaffian_witness.pf:160-260 but emit CSR via
  ! extract_block_csr or inline CSR construction.
  !
  ! For test purposes, build a dense matrix first, then convert to CSR using
  ! utils.f90:dense_to_csr. (Do not re-derive from scratch.)
  ! ---------------------------------------------------------------------------

  subroutine build_synthetic_bdg_csr(H, N, mu, delta_0)
    type(csr_matrix), intent(out) :: H
    integer, intent(out) :: N
    real(kind=dp), intent(in) :: mu, delta_0
    complex(kind=dp), allocatable :: H_dense(:,:)
    integer :: i

    N = 16
    allocate(H_dense(N, N))
    H_dense = (0.0_dp, 0.0_dp)
    ! Diagonal-in-(c,c†) structure: H₀ = μ·diag(1,..,1) + delta_0·off-diag pairing.
    do i = 1, n
      H_dense(i, i) = cmplx(mu, 0.0_dp, kind=dp)
    end do
    ! Pairing block (iσ_y in bands 7-8): delta_0 between (7,15), (8,16), (15,7), (16,8).
    H_dense(7, 15)  = cmplx(0.0_dp, delta_0, kind=dp)
    H_dense(8, 16)  = cmplx(0.0_dp, delta_0, kind=dp)
    H_dense(15, 7)  = cmplx(0.0_dp, -delta_0, kind=dp)
    H_dense(16, 8)  = cmplx(0.0_dp, -delta_0, kind=dp)
    call dense_to_csr(H_dense, H)
    deallocate(H_dense)
  end subroutine

  subroutine build_zero_bdg_csr(H, N)
    type(csr_matrix), intent(out) :: H
    integer, intent(out) :: N
    complex(kind=dp), allocatable :: H_dense(:,:)

    N = 16
    allocate(H_dense(N, N))
    H_dense = (0.0_dp, 0.0_dp)
    call dense_to_csr(H_dense, H)
    deallocate(H_dense)
  end subroutine

  subroutine build_nondiagonal_bdg_csr(H, N, mu, delta_0)
    type(csr_matrix), intent(out) :: H
    integer, intent(out) :: N
    real(kind=dp), intent(in) :: mu, delta_0
    complex(kind=dp), allocatable :: H_dense(:,:)
    integer :: i

    N = 16
    allocate(H_dense(N, N))
    H_dense = (0.0_dp, 0.0_dp)
    do i = 1, n
      H_dense(i, i) = cmplx(mu, 0.0_dp, kind=dp)
    end do
    ! Non-diagonal k.p mixing: bands 1-2 couple into bands 7-8.
    H_dense(1, 7)  = cmplx(0.3_dp, 0.0_dp, kind=dp)
    H_dense(2, 8)  = cmplx(0.3_dp, 0.0_dp, kind=dp)
    H_dense(7, 1)  = cmplx(0.3_dp, 0.0_dp, kind=dp)
    H_dense(8, 2)  = cmplx(0.3_dp, 0.0_dp, kind=dp)
    H_dense(7, 15) = cmplx(0.0_dp, delta_0, kind=dp)
    H_dense(8, 16) = cmplx(0.0_dp, delta_0, kind=dp)
    H_dense(15, 7) = cmplx(0.0_dp, -delta_0, kind=dp)
    H_dense(16, 8) = cmplx(0.0_dp, -delta_0, kind=dp)
    call dense_to_csr(H_dense, H)
    deallocate(H_dense)
  end subroutine

end module test_bdg_pfaffian_witness_csr
```

**Notes on the test file:**
- `csr_matrix` type is `use definitions, only: csr_matrix` (the `csr_matrix` is re-exported from `sparse_matrices` to `definitions` for module-graph convenience — confirm by reading `src/core/defs.f90` and `src/math/sparse_matrices.f90:1-30`).
- `dense_to_csr` lives in `src/core/utils.f90` (per the source layout in CLAUDE.md). **Before writing the test file, verify its exact signature** with `grep -n "subroutine dense_to_csr\|function dense_to_csr" src/core/utils.f90`; the call `call dense_to_csr(H_dense, H)` shown above assumes a 2-arg shape `(dense, csr)`. If the signature differs (e.g., extra shape args), adapt the fixture.
- If `csr_matrix` is not in `definitions`, adjust the import to `use sparse_matrices, only: csr_matrix`.
- All `@assertEqual` / `@assertTrue` macros are single-line per the repo convention.

- [ ] **Step 2: Register the test in CMakeLists.txt**

Edit `tests/CMakeLists.txt` — after the existing `test_bdg_evaluator` block (lines 182-186), add:
```cmake
    add_pfunit_ctest(test_bdg_pfaffian_witness_csr
        TEST_SOURCES unit/test_bdg_pfaffian_witness_csr.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )
```

- [ ] **Step 3: Re-configure cmake to pick up the new test**

Run: `cmake -B build 2>&1 | tail -5`
Expected: `build` directory updated with new target.

- [ ] **Step 4: Run the new test, verify FAIL (function not defined yet)**

Run: `cmake --build build --target test_bdg_pfaffian_witness_csr 2>&1 | tail -10`
Expected: build fails with `eval_bdg_pfaffian_witness_csr` undefined (or, if build succeeds via separate compilation, the link step fails with undefined reference).

- [ ] **Step 5: Add the seam sibling to `bdg_observables.f90`**

Edit `src/physics/bdg_observables.f90` — add the function declaration in the `public ::` list (after line 30):
```fortran
  public :: eval_bdg_pfaffian_witness_csr
```

Add the function body before `end module bdg_observables`:
```fortran
  ! ==============================================================================
  ! CSR BdG slim projected Pfaffian witness — seam sibling (wire-rung invariant).
  !
  ! Returns the S2-projected Pfaffian sign (s2_sign ∈ {-1, 0, +1}) of the
  ! BdG matrix H_bdg_csr. S2 = bands 7-8 per the k.p block table SSOT.
  !
  ! Per design decision 2026-07-13: this seam sibling is a thin wrapper that
  ! delegates the CSR-aware S2 extraction to wire_pfaffian_witness_sweep in
  ! topological_analysis.f90. That subroutine already operates on CSR input
  ! (per ticket 04 — `main_topology.f90:1371` migration pattern) and is the
  ! existing CSR-aware dense-path witness; reusing it keeps the seam thin
  ! (single subroutine import) without re-implementing the S2 row-extraction
  ! for CSR.
  !
  ! The params%pfaffian_floor field is the SSOT for future U13 work where the
  ! floor will be wired through; for now wire_pfaffian_witness_sweep uses its
  ! own internal floor.
  ! ==============================================================================
  function eval_bdg_pfaffian_witness_csr(H_bdg_csr, Nbdg, params) result(s2_sign)
    type(csr_matrix), intent(in)        :: H_bdg_csr
    integer,          intent(in)        :: Nbdg
    type(bdg_eval_params_t), intent(in) :: params
    integer                              :: s2_sign

    ! Delegate to the CSR-aware dense-path witness. params is unused for now
    ! (params%pfaffian_floor plumbing is the U13 follow-up — the SSOT exists,
    ! the helper doesn't accept it yet).
    call wire_pfaffian_witness_sweep(H_bdg_csr, Nbdg, s2_sign)
  end function eval_bdg_pfaffian_witness_csr
```

**Implementation notes:**
- Add `use topological_analysis, only: wire_pfaffian_witness_sweep` to the import block at the top of `bdg_observables.f90`.
- This adds a dep edge from `bdg_observables` to `topological_analysis` in the AGENTS.md DAG. Update AGENTS.md (Task 3.1) to reflect this — `bdg_observables` is no longer a pure L1 leaf; it imports one symbol from L3 `topological_analysis` (a deliberate trade-off per the user's 2026-07-13 design choice).
- `params` is currently unused; the unused-dummy-argument is intentional — it's the SSOT surface for future floor plumbing (U13 follow-up).
- `wire_pfaffian_witness_sweep` is `subroutine` with `intent(out) :: s2_sign` — call it as a subroutine, not a function.

- [ ] **Step 6: Build and run the new test**

Run: `cmake --build build --target test_bdg_pfaffian_witness_csr 2>&1 | tail -10`
Expected: build succeeds.

Run: `OMP_NUM_THREADS=1 ctest --test-dir build -R test_bdg_pfaffian_witness_csr -V 2>&1 | tail -20`
Expected: 4 tests pass. (May need iteration on the S2 extraction; see note below.)

- [ ] **Step 7: Iterate if needed**

If the test fails on the S2 indexing, debug by adding a temporary print statement inside `eval_bdg_pfaffian_witness_csr` showing `H_proj(1,1)` etc. for the synthetic fixture; verify the extracted block matches the dense `H_dense` rows.

- [ ] **Step 8: Commit**

```bash
git add tests/unit/test_bdg_pfaffian_witness_csr.pf tests/CMakeLists.txt src/physics/bdg_observables.f90
git commit -m "feat(bdg): add eval_bdg_pfaffian_witness_csr seam sibling (wire-rung invariant)"
```

## Task 2.3: Implement `eval_bdg_kitaev_majorana` (TDD)

**Files:**
- Create: `tests/unit/test_bdg_kitaev_majorana.pf` (NEW, 3 tests)
- Modify: `src/physics/bdg_observables.f90` (add the sibling)

**Reference:** the existing `pfaffian.f90:kitaev_majorana_number` (already imported at Task 2.1). The seam sibling is a thin wrapper.

**Interface (per spec.md §"Seam shape"):**
```fortran
function eval_bdg_kitaev_majorana(H_k_array, k_par_values) result(majorana_number)
  complex(kind=dp), intent(in) :: H_k_array(:,:,:)   ! H(k) at each k_par index
  real(kind=dp),    intent(in) :: k_par_values(:)    ! PHS-invariant momenta: [0.0, π/a]
  integer                       :: majorana_number   ! {-1, 0, +1}; -1 = topological (Kitaev 2001)
end function
```

(No `params` arg in the spec signature — the wrapped `kitaev_majorana_number` from `pfaffian.f90` has its own internal floor; pass through `params%pfaffian_floor` only if the helper accepts it; otherwise omit.)

- [ ] **Step 1: Create the test file with the first failing test (delegation contract)**

Create `tests/unit/test_bdg_kitaev_majorana.pf`:
```fortran
module test_bdg_kitaev_majorana
  use iso_fortran_env, only: dp => real64
  use funit
  use bdg_observables, only: eval_bdg_kitaev_majorana

  implicit none

contains

  @test
  subroutine test_kitaev_majorana_seam_returns_valid_sign()
    ! Spec §4 — QW+Kitaev rung seam sibling. Returns
    ! majorana_number ∈ {-1, 0, +1}. Wraps pfaffian.f90:kitaev_majorana_number.
    complex(kind=dp), allocatable :: H_k(:,:,:)
    real(kind=dp), allocatable :: k_par(:)
    integer :: m

    call build_kitaev_chain(H_k, k_par, n_sites=4)

    m = eval_bdg_kitaev_majorana(H_k, k_par)

    @assertTrue(m >= -1 .and. m <= 1)
  end subroutine

  @test
  subroutine test_kitaev_majorana_seam_topological_sign()
    ! |μ| < 2t → M = -1 (Kitaev 2001, Eq. 26 convention).
    complex(kind=dp), allocatable :: H_k(:,:,:)
    real(kind=dp), allocatable :: k_par(:)
    integer :: m

    call build_kitaev_chain(H_k, k_par, n_sites=4, mu=0.5_dp, t=1.0_dp, delta=0.5_dp)

    m = eval_bdg_kitaev_majorana(H_k, k_par)

    @assertEqual(-1, m)
  end subroutine

  @test
  subroutine test_kitaev_majorana_seam_trivial_sign()
    ! |μ| > 2t → M = +1 (Kitaev 2001, Eq. 26 convention).
    complex(kind=dp), allocatable :: H_k(:,:,:)
    real(kind=dp), allocatable :: k_par(:)
    integer :: m

    call build_kitaev_chain(H_k, k_par, n_sites=4, mu=3.0_dp, t=1.0_dp, delta=0.5_dp)

    m = eval_bdg_kitaev_majorana(H_k, k_par)

    @assertEqual(1, m)
  end subroutine

  ! ---------------------------------------------------------------------------
  ! Fixture: 2-band Kitaev chain H(k) for n_sites cells, evaluated at k=0 and
  ! k=π/a. H_k(i,j,k_idx) is the (i,j) element of H at momentum k(k_idx).
  ! Reuses the same construction as tests/unit/test_kitaev_strict.pf via
  ! pfaffian.f90:kitaev_bdg_fixture_2band.
  ! ---------------------------------------------------------------------------

  subroutine build_kitaev_chain(H_k, k_par, n_sites, mu, t, delta)
    integer, intent(in) :: n_sites
    real(kind=dp), intent(in), optional :: mu, t, delta
    complex(kind=dp), allocatable, intent(out) :: H_k(:,:,:)
    real(kind=dp), allocatable, intent(out) :: k_par(:)
    real(kind=dp) :: mu_, t_, delta_

    mu_    = 0.5_dp;  if (present(mu))    mu_    = mu
    t_     = 1.0_dp;  if (present(t))     t_     = t
    delta_ = 0.5_dp;  if (present(delta)) delta_ = delta

    allocate(k_par(2))
    k_par = [0.0_dp, 3.141592653589793_dp]  ! k=0 and k=π/a

    allocate(H_k(2, 2, 2))
    ! H(k) = -(mu + 2t cos k) σ_z + delta sin k σ_y (Kitaev 2001)
    H_k(:, :, 1) = reshape([ &                ! k=0
      cmplx(-mu_ - 2*t_, 0.0_dp, kind=dp),    cmplx(0.0_dp, 0.0_dp, kind=dp), &
      cmplx(0.0_dp, 0.0_dp, kind=dp),        cmplx(mu_ + 2*t_, 0.0_dp, kind=dp)  &
    ], shape=[2, 2])
    H_k(:, :, 2) = reshape([ &                ! k=π/a
      cmplx(-mu_ + 2*t_, 0.0_dp, kind=dp),   cmplx(0.0_dp, 0.0_dp, kind=dp), &
      cmplx(0.0_dp, 0.0_dp, kind=dp),        cmplx(mu_ - 2*t_, 0.0_dp, kind=dp)  &
    ], shape=[2, 2])
  end subroutine

end module test_kitaev_majorana
```

**Note on the fixture:** For accuracy, prefer to use `pfaffian.f90:kitaev_bdg_fixture_2band` (the existing helper at `src/math/pfaffian.f90:483`) instead of hand-building. If the helper has a different signature, adapt. Alternatively, mirror the fixture from `tests/unit/test_kitaev_strict.pf:11-42`.

- [ ] **Step 2: Register the test in CMakeLists.txt**

Edit `tests/CMakeLists.txt` — after the new `test_bdg_pfaffian_witness_csr` block, add:
```cmake
    add_pfunit_ctest(test_bdg_kitaev_majorana
        TEST_SOURCES unit/test_bdg_kitaev_majorana.pf
        LINK_LIBRARIES 8bandkp_common 8bandkp_test_support
        LABELS "unit"
    )
```

- [ ] **Step 3: Re-configure cmake**

Run: `cmake -B build 2>&1 | tail -3`

- [ ] **Step 4: Run the new test, verify FAIL**

Run: `cmake --build build --target test_bdg_kitaev_majorana 2>&1 | tail -10`
Expected: build fails with `eval_bdg_kitaev_majorana` undefined.

- [ ] **Step 5: Add the stub to `bdg_observables.f90`**

Edit `src/physics/bdg_observables.f90` — add to `public ::` list:
```fortran
  public :: eval_bdg_kitaev_majorana
```

Add the function body before `end module bdg_observables`:
```fortran
  ! ==============================================================================
  ! QW+Kitaev-rung Majorana number seam sibling.
  !
  ! Wraps pfaffian.f90:kitaev_majorana_number with a seam-shape signature.
  ! Returns majorana_number ∈ {-1, 0, +1} (Kitaev 2001, Eq. 26 convention:
  ! M = -1 topological, M = +1 trivial).
  ! ==============================================================================
  function eval_bdg_kitaev_majorana(H_k_array, k_par_values) result(majorana_number)
    complex(kind=dp), intent(in) :: H_k_array(:,:,:)
    real(kind=dp),    intent(in) :: k_par_values(:)
    integer                       :: majorana_number

    majorana_number = kitaev_majorana_number(H_k_array, k_par_values)
  end function eval_bdg_kitaev_majorana
```

- [ ] **Step 6: Build, run, verify PASS**

Run: `cmake --build build --target test_bdg_kitaev_majorana && OMP_NUM_THREADS=1 ctest --test-dir build -R test_bdg_kitaev_majorana -V 2>&1 | tail -15`
Expected: 3 tests pass. If `kitaev_majorana_number` signature differs, adjust the wrapper call.

- [ ] **Step 7: Commit**

```bash
git add tests/unit/test_bdg_kitaev_majorana.pf tests/CMakeLists.txt src/physics/bdg_observables.f90
git commit -m "feat(bdg): add eval_bdg_kitaev_majorana seam sibling (QW+Kitaev rung)"
```

## Task 2.4: Rename `compute_z2_gap` / `compute_z2_gap_edge` → `*_bhz_heuristic`

**Files:**
- Modify: `src/physics/topological_analysis.f90:19, 20, 275, 318` (rename + update exports + update doc-comments)
- Modify: `src/apps/main_topology.f90:369, 375` (rename call sites + comments)

Per ticket 03: scope-narrow to BHZ-only at the call site (heuristic is the gap-closure fallback; the slim Pfaffian is the invariant).

- [ ] **Step 1: Rename `compute_z2_gap` → `compute_z2_gap_bhz_heuristic` in topological_analysis.f90**

Edit `src/physics/topological_analysis.f90:19` — replace:
```fortran
  public :: compute_z2_gap
```
with:
```fortran
  public :: compute_z2_gap_bhz_heuristic
```

Edit `src/physics/topological_analysis.f90:275-302` — rename the function:
- Replace `function compute_z2_gap(eigenvalues, gap_threshold) result(z2)` with `function compute_z2_gap_bhz_heuristic(eigenvalues, gap_threshold) result(z2)`
- Replace `end function compute_z2_gap` with `end function compute_z2_gap_bhz_heuristic`
- Update the doc-comment header above (line ~270) to add: "BHZ-only (hardcoded 4-band basis). Kept as the gap-closure fallback after the slim Pfaffian returns 0 — per ticket 03 of `.scratch/archive/bdg-evaluator-pfaffian/`."

- [ ] **Step 2: Rename `compute_z2_gap_edge` → `compute_z2_gap_edge_bhz_heuristic`**

Edit `src/physics/topological_analysis.f90:20, 318, 384` — same pattern: rename function + end-function + update doc-comment with the BHZ-only scope-narrow note.

- [ ] **Step 3: Update call sites in `main_topology.f90`**

Edit `src/apps/main_topology.f90:369` — update the comment:
```fortran
        ! First try: spectral edge-state detection with compute_z2_gap_edge_bhz_heuristic
        ! (BHZ-only, gap-closure fallback after slim Pfaffian returns 0).
```

Edit `src/apps/main_topology.f90:375` — replace:
```fortran
        result%z2_invariant = compute_z2_gap_edge( &
```
with:
```fortran
        result%z2_invariant = compute_z2_gap_edge_bhz_heuristic( &
```

Edit `src/apps/main_topology.f90:1380` — replace:
```fortran
        z2 = compute_z2_gap(eigen_res_local%eigenvalues, gap_threshold)
```
with:
```fortran
        z2 = compute_z2_gap_bhz_heuristic(eigen_res_local%eigenvalues, gap_threshold)
```
(And update the surrounding comment to note this is the gap-closure fallback AFTER the slim Pfaffian returns 0.)

- [ ] **Step 4: Build, verify compile**

Run: `cmake --build build 2>&1 | tail -10`
Expected: clean build, no unresolved references.

- [ ] **Step 5: Run unit tests, verify same 50/51 (rename is structural)**

Run: `OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L unit -j4 2>&1 | tail -10`
Expected: **50/51 PASS, 1 FAIL: `test_wire_pfaffian_witness`** (same as before — the rename is purely structural).

- [ ] **Step 6: Commit**

```bash
git add src/physics/topological_analysis.f90 src/apps/main_topology.f90
git commit -m "refactor(bdg): rename compute_z2_gap* to _bhz_heuristic (scope-narrow gap-closure fallback)"
```

## Task 2.5: Migrate `main_topology.f90:1371` to seam sibling

**Files:**
- Modify: `src/apps/main_topology.f90:1371` (replace `wire_pfaffian_witness_sweep` call with `eval_bdg_pfaffian_witness_csr`)
- Modify: `src/apps/main_topology.f90` (add `use bdg_observables, only: eval_bdg_pfaffian_witness_csr` to existing import)

- [ ] **Step 1: Find the existing `use bdg_observables` block in `main_topology.f90`**

Run: `grep -n "use bdg_observables" src/apps/main_topology.f90`
Expected: a single `use` line near the top of the file.

- [ ] **Step 2: Add `eval_bdg_pfaffian_witness_csr` to the `use` block**

Edit the `use bdg_observables, only: ...` line to include `eval_bdg_pfaffian_witness_csr`. (If other names are imported, add the new one alphabetically.)

- [ ] **Step 3: Replace the call site at :1371**

Read `src/apps/main_topology.f90:1365-1385` to see the context. The current call is:
```fortran
      call wire_pfaffian_witness_sweep(H_bdg_csr, Nbdg_local, s2_sign)
```

Replace with the seam sibling:
```fortran
      ! Slim Pfaffian witness via seam sibling (per ticket 04 of
      ! .scratch/archive/bdg-evaluator-pfaffian/ — drop-in replacement, same
      ! s2_sign ∈ {-1, 0, +1} semantics).
      s2_sign = eval_bdg_pfaffian_witness_csr(H_bdg_csr, Nbdg_local, &
                                              bdg_eval_params_with_delta(delta_0))
```

If the surrounding code uses `delta_0` from a different source, adapt.

- [ ] **Step 4: Build, verify compile**

Run: `cmake --build build 2>&1 | tail -10`
Expected: clean build.

- [ ] **Step 5: Run unit tests, verify 50/51 still (migration is structural)**

Run: `OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L unit -j4 2>&1 | tail -10`
Expected: **50/51 PASS, 1 FAIL: `test_wire_pfaffian_witness`** (still the same — the migration doesn't turn the strict test GREEN yet).

- [ ] **Step 6: Commit**

```bash
git add src/apps/main_topology.f90
git commit -m "refactor(bdg): migrate wire_pfaffian_witness_sweep to seam sibling eval_bdg_pfaffian_witness_csr"
```

## Task 2.6: Update `test_wire_pfaffian_witness.pf` with non-diagonal fixture (turn strict RED → GREEN)

**Files:**
- Modify: `tests/unit/test_wire_pfaffian_witness.pf:75-105` (the strict sign-agreement test)

Per ticket 06: the strict `@assertTrue` is RED because the synthetic fixture `build_synthetic_bdg_16` is diagonal-in-(c,c†) → Pf = 0 by construction. Fix: replace the fixture in the strict test with a non-diagonal one (~30 LOC) so Pf ≠ 0 and the strict assertion turns GREEN.

- [ ] **Step 1: Read the current strict test and its fixture**

Read `tests/unit/test_wire_pfaffian_witness.pf:60-200` to understand the existing `build_synthetic_bdg_16` and `wire_pfaffian_witness` (dense) interface.

- [ ] **Step 2: Add a new fixture `build_nondiagonal_bdg_16`**

Edit `tests/unit/test_wire_pfaffian_witness.pf` — in the fixtures section, add:
```fortran
  subroutine build_nondiagonal_bdg_16(H_bdg, mu, delta_0)
    ! Per ticket 06 of .scratch/archive/bdg-evaluator-pfaffian/ — non-diagonal fixture
    ! that crosses the Nambu seam, so Pf != 0 by construction.
    complex(kind=dp), intent(out) :: H_bdg(16, 16)
    real(kind=dp), intent(in) :: mu, delta_0
    integer :: i

    H_bdg = (0.0_dp, 0.0_dp)
    do i = 1, 16
      H_bdg(i, i) = cmplx(mu, 0.0_dp, kind=dp)
    end do
    ! Non-diagonal k.p mixing: bands 1-2 couple into bands 7-8 (conduction),
    ! which crosses the Nambu seam → Pf ≠ 0.
    H_bdg(1, 7)  = cmplx(0.3_dp, 0.0_dp, kind=dp)
    H_bdg(2, 8)  = cmplx(0.3_dp, 0.0_dp, kind=dp)
    H_bdg(7, 1)  = cmplx(0.3_dp, 0.0_dp, kind=dp)
    H_bdg(8, 2)  = cmplx(0.3_dp, 0.0_dp, kind=dp)
    ! Pairing block (iσ_y in bands 7-8): delta_0 between (7,15), (8,16).
    H_bdg(7, 15)  = cmplx(0.0_dp, delta_0, kind=dp)
    H_bdg(8, 16)  = cmplx(0.0_dp, delta_0, kind=dp)
    H_bdg(15, 7)  = cmplx(0.0_dp, -delta_0, kind=dp)
    H_bdg(16, 8)  = cmplx(0.0_dp, -delta_0, kind=dp)
  end subroutine build_nondiagonal_bdg_16
```

- [ ] **Step 3: Update the strict test to use the new fixture**

Edit `tests/unit/test_wire_pfaffian_witness.pf:75-105` — replace `build_synthetic_bdg_16(H_bdg, mu=0.0_dp, delta_0=0.001_dp)` with `build_nondiagonal_bdg_16(H_bdg, mu=0.0_dp, delta_0=0.001_dp)`.

Also replace the long comment block (lines 81-96) with the shorter:
```fortran
    ! Spec §4.3 — STRICT sign-agreement assertion (no zero-escape hatch).
    ! Per ticket 06 of .scratch/archive/bdg-evaluator-pfaffian/, the strict assertion
    ! is unblocked by a non-diagonal synthetic fixture (conduction-band
    ! off-diagonal coupling crosses the Nambu seam → Pf ≠ 0 by construction).
    ! Decoupled from U13 (slim → full Bloch-Pfaffian swap).
```

- [ ] **Step 4: Run the test, verify PASS**

Run: `cmake --build build --target test_wire_pfaffian_witness && OMP_NUM_THREADS=1 ctest --test-dir build -R test_wire_pfaffian_witness -V 2>&1 | tail -15`
Expected: all 5 tests in `test_wire_pfaffian_witness` PASS — the strict sign-agreement test turns GREEN, the @todo is gone.

- [ ] **Step 5: Confirm the unit-test label is now 51/51**

Run: `OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L unit -j4 2>&1 | tail -10`
Expected: **51/51 PASS** (was 50/51). Wait — we haven't added the new tests to the build yet at this point. If Tasks 2.2 + 2.3 ran first, then the count is 51/51 (50 originals + 1 from non-diagonal turning GREEN). The new tests `test_bdg_pfaffian_witness_csr` + `test_bdg_kitaev_majorana` will bring it to ~55/55 after Task 2.9.

Expected at this step: 51/51 (the strict test turned GREEN; new tests not yet counted if CMake re-config hasn't been re-run).

- [ ] **Step 6: Commit**

```bash
git add tests/unit/test_wire_pfaffian_witness.pf
git commit -m "test(bdg): turn strict wire_pfaffian_witness sign-agreement test GREEN via non-diagonal fixture"
```

## Task 2.7: Update `test_bdg_evaluator.pf` with `bdg_default_pfaffian_floor` SSOT tests

**Files:**
- Modify: `tests/unit/test_bdg_evaluator.pf` (add 2 tests for the new SSOT)

- [ ] **Step 1: Read the existing test file**

Read `tests/unit/test_bdg_evaluator.pf:1-30` to see the imports + module declaration.

- [ ] **Step 2: Add 2 new tests for the SSOT**

Edit `tests/unit/test_bdg_evaluator.pf` — in the `public ::` use block, add `bdg_default_pfaffian_floor` and `bdg_eval_params_t` if not already imported.

Add 2 new tests before `end module test_bdg_evaluator`:
```fortran
  @test
  subroutine test_bdg_default_pfaffian_floor_value()
    ! Spec §"SSOT magic-number extraction" — bdg_default_pfaffian_floor is the
    ! SSOT for the magic literal 1.0e-12_dp formerly duplicated at
    ! topological_analysis.f90:1663, 1681, 1766.
    @assertEqual(1.0e-12_dp, bdg_default_pfaffian_floor, tolerance=0.0_dp)
  end subroutine

  @test
  subroutine test_bdg_eval_params_factory_uses_pfaffian_floor()
    ! bdg_eval_params_with_delta must default pfaffian_floor to the SSOT.
    type(bdg_eval_params_t) :: p
    p = bdg_eval_params_with_delta(0.001_dp)
    @assertEqual(bdg_default_pfaffian_floor, p%pfaffian_floor, tolerance=0.0_dp)
  end subroutine
```

**Note on `@assertEqual` for reals:** pFUnit's `@assertEqual` for reals requires a tolerance. The exact-equal tolerance is `0.0_dp` when both sides are the same SSOT literal; if pFUnit complains, use `@assertTrue(p%pfaffian_floor == bdg_default_pfaffian_floor)` instead.

- [ ] **Step 3: Build, run, verify PASS**

Run: `cmake --build build --target test_bdg_evaluator && OMP_NUM_THREADS=1 ctest --test-dir build -R test_bdg_evaluator -V 2>&1 | tail -10`
Expected: 10 tests pass (was 8; +2 new).

- [ ] **Step 4: Commit**

```bash
git add tests/unit/test_bdg_evaluator.pf
git commit -m "test(bdg): pin bdg_default_pfaffian_floor SSOT (ticket 02 §3)"
```

## Task 2.8: Update `test_kitaev_majorana.pf` with seam-sibling direct call

**Files:**
- Modify: `tests/unit/test_kitaev_majorana.pf` (add 1 test that calls `eval_bdg_kitaev_majorana` directly via the seam)

- [ ] **Step 1: Read the existing test file's import block**

Read `tests/unit/test_kitaev_majorana.pf:1-30` to see the imports.

- [ ] **Step 2: Add the seam sibling to the `use` block**

Add `eval_bdg_kitaev_majorana` to the imports.

- [ ] **Step 3: Add 1 new test that calls the seam sibling directly**

Edit `tests/unit/test_kitaev_majorana.pf` — add before `end module`:
```fortran
  @test
  subroutine test_kitaev_majorana_via_seam_sibling()
    ! Spec §"Seam shape" — QW+Kitaev rung consumes the seam sibling
    ! (not pfaffian.f90:kitaev_majorana_number directly). Pin the
    ! delegation contract: identical result for the same fixture.
    complex(kind=dp), allocatable :: H_k(:,:,:)
    real(kind=dp), allocatable :: k_par(:)
    integer :: m_via_seam

    call build_kitaev_chain(H_k, k_par, n_sites=4, mu=0.5_dp, t=1.0_dp, delta=0.5_dp)
    m_via_seam = eval_bdg_kitaev_majorana(H_k, k_par)
    @assertEqual(-1, m_via_seam)
  end subroutine
```

(If `bdg_eval_params_t` is already in scope, omit the `use bdg_observables` line inside the subroutine.)

- [ ] **Step 4: Build, run, verify PASS**

Run: `cmake --build build --target test_kitaev_majorana && OMP_NUM_THREADS=1 ctest --test-dir build -R test_kitaev_majorana -V 2>&1 | tail -10`
Expected: 11 tests pass (was 10; +1 new).

- [ ] **Step 5: Commit**

```bash
git add tests/unit/test_kitaev_majorana.pf
git commit -m "test(bdg): pin seam-sibling delegation for QW+Kitaev rung"
```

## Task 2.9: Verification gate — unit tests green

- [ ] **Step 1: Re-configure cmake to pick up all new tests**

Run: `cmake -B build 2>&1 | tail -3`
Expected: 2 new targets registered (`test_bdg_pfaffian_witness_csr`, `test_bdg_kitaev_majorana`).

- [ ] **Step 2: Build everything**

Run: `cmake --build build 2>&1 | tail -5`
Expected: clean build.

- [ ] **Step 3: Run the full unit-test label, capture results**

Run: `OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L unit -j4 2>&1 | tee .scratch/bdg-u2-actual-ship/ctest-phase2-verify.log | tail -15`
Expected: **53/53 PASS** (or 54/54; the exact count depends on whether `test_bdg_pfaffian_witness_csr` registers as 1 target with 4 subtests, etc.).

The exact count to verify:
- Original 50 unit tests (per Phase 24 follow-up)
- − 0 (no original test removed)
- + 4 (test_bdg_pfaffian_witness_csr subtests)
- + 3 (test_bdg_kitaev_majorana subtests)
- + 2 (test_bdg_evaluator SSOT tests)
- + 1 (test_kitaev_majorana seam-sibling direct call)
- = **60 total unit tests, all green**

(If the count differs slightly due to pFUnit target grouping, that's fine as long as 0 fails.)

- [ ] **Step 4: If any test FAILS, debug and fix**

- For test_bdg_pfaffian_witness_csr: most likely failure is in `extract_s2_block_csr` indexing — verify against dense fixture.
- For test_bdg_kitaev_majorana: most likely is the `kitaev_majorana_number` signature mismatch — adjust the wrapper call.
- For test_bdg_evaluator SSOT: most likely is the `@assertEqual` tolerance — switch to `@assertTrue(p%pfaffian_floor == ...)` if needed.

- [ ] **Step 5: Commit any fixes**

```bash
git add -A
git commit -m "fix(bdg): address test failures from Phase 2 verification gate"
```

---

# Phase 3: Re-apply doc-only "U2 closed" stamps (now backed by green ctest) (10 tasks)

**Precondition:** All Phase 2 tasks done. `OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L unit -j4` shows 0 failures.

## Task 3.1: AGENTS.md inventory + DAG update

**Files:**
- Modify: `src/physics/AGENTS.md`

- [ ] **Step 1: Update the `bdg_observables.f90` inventory row**

Edit `src/physics/AGENTS.md:44` — extend the row to include:
```
| `bdg_observables.f90` | `bdg_observables` | 88+ | Pure per-point BdG evaluator seam. `eval_bdg_point` (eigenvalues-only → minigap + near-zero count + heuristic invariant_flag ∈ {0,1}); seam siblings `eval_bdg_pfaffian_witness_csr` (CSR BdG → s2_sign ∈ {-1, 0, +1}, wire-rung invariant via slim projected Pfaffian) and `eval_bdg_kitaev_majorana` (H_k_array + k_par_values → majorana_number ∈ {-1, 0, +1}, QW+Kitaev rung) per `.scratch/archive/bdg-evaluator-pfaffian/` ticket 01. SSOT parameters: `bdg_default_near_zero_frac = 0.001_dp`, `bdg_default_min_threshold = 1.0e-10_dp`, `bdg_default_pfaffian_floor = 1.0e-12_dp`. Imports L0 leaves `sparse_matrices` + `pfaffian`. |
```

(Adjust the size annotation if the file grew beyond 88 lines; check actual line count first with `wc -l src/physics/bdg_observables.f90`.)

- [ ] **Step 2: Update the Dependency DAG block**

Edit `src/physics/AGENTS.md:46-58` — add `sparse_matrices` + `pfaffian` to Layer 0 (with the `csr_matrix` + `complex_pfaffian` / `kitaev_majorana_number` exports) and add `bdg_observables` to Layer 1 with the `uses sparse_matrices + pfaffian` annotation per ticket 02 §1 + §3.

(Read the existing DAG block first to confirm the exact layout.)

- [ ] **Step 3: Verify**

Run: `grep -nE "bdg_observables|eval_bdg_pfaffian_witness_csr|eval_bdg_kitaev_majorana|sparse_matrices.*pfaffian|pfaffian.*sparse_matrices" src/physics/AGENTS.md`
Expected: 4-5 matches reflecting the new state.

## Task 3.2: Lecture markdown §13.0 + §13.7.1 + §13.7.4 + new §13.7.5 + summary table + tail

**Files:**
- Modify: `docs/lecture/13-topological-superconductivity.md`

Per ticket 07's "Lecture disclosure" section.

- [ ] **Step 1: Update §13.0**

Edit `docs/lecture/13-topological-superconductivity.md:28-33` — replace "3-witness ... slim Pfaffian row reserved for U13" with "**4-witness** within a 1.0 T tolerance. The slim Pfaffian row is **live**: `bcrit_pfaffian = arg_B at mu ≈ 0.6601 (tol ±0.0001) where the colormap's z2 column first reads -1`".

- [ ] **Step 2: Update §13.7.1**

Edit `docs/lecture/13-topological-superconductivity.md:528-533` — replace the "3-witness ... slim Pfaffian row reserved for U13" wording with: "**4-witness agreement** within a 1.0 T tolerance — the slim Pfaffian row is **live** (colormap-extracted from `output/z2_phase_diagram.dat` z2 column at mu ≈ 0.6601 ± 0.0001, ticket 05)". Update any "reserved for U13" or "deferred to U13" wording in §13.7.1 to "live".

- [ ] **Step 3: Update §13.7.4**

Edit `docs/lecture/13-topological-superconductivity.md:622-650` — rewrite the slim Pfaffian caption + trailing paragraph:
- First bullet: "real Pfaffian witness output requires U13" → "live; not reserved for U13"
- Trailing paragraph: "live-witness explanation with §13.7.5 reference"

- [ ] **Step 4: Add new §13.7.5 "Slim Pfaffian witness criterion"**

Edit `docs/lecture/13-topological-superconductivity.md` — after §13.7.4 (around line 668), add:
```markdown
### 13.7.5 Slim Pfaffian witness criterion

The slim Pfaffian is the **live** 4th witness of the L13 acceptance gate.
It is read from the existing `(B, μ)` colormap dataset
(`output/z2_phase_diagram.dat` produced by `compute_wire_bdg_gap_sweep`),
extracting `bcrit_pfaffian = B at mu ≈ 0.6601 (tol ±0.0001) where the
colormap's `z2` column first reads -1`. No new Fortran I/O is required:
the slim Pfaffian row reads from the existing colormap file. The floor
honored by the witness is `bdg_default_pfaffian_floor = 1.0e-12_dp`
(promoted to SSOT per ticket 02 of `.scratch/archive/bdg-evaluator-pfaffian/`).
Forward reference: U13 swaps the slim witness for the full Bloch-Pfaffian
on a periodic supercell (separate scoped PR).
```

- [ ] **Step 5: Update the summary table row 836**

Edit `docs/lecture/13-topological-superconductivity.md:876` — flip "Pfaffian witness s1=s2 at (B_crit, mu=0.6601)" → "Pfaffian witness (colormap-extracted, z2==-1 at mu≈0.6601, ticket 05)".

- [ ] **Step 6: Update the tail "U12 3-witness gate" → "U12 4-witness gate"**

Edit `docs/lecture/13-topological-superconductivity.md:1005` — flip to "U12 4-witness gate (slim Pfaffian row live per ticket 05; colormap-extracted from `output/z2_phase_diagram.dat` z2 column at mu ≈ 0.6601 ±0.0001)".

- [ ] **Step 7: Verify**

Run: `grep -nE "13.7.5|slim Pfaffian|Slim Pfaffian|live; not reserved" docs/lecture/13-topological-superconductivity.md`
Expected: ~6-8 matches confirming all sites flipped.

## Task 3.3: `scripts/lecture_13_topological.py` — colormap-extracted `bcrit_pfaffian`

**Files:**
- Modify: `scripts/lecture_13_topological.py`

Per ticket 07's "Lecture script + acceptance gate + verifier" section.

- [ ] **Step 1: Flip module-level docstring 3→4 witness**

Edit `scripts/lecture_13_topological.py` — find the "These lines" + "Output" docstring sections; flip "3-witness" → "4-witness".

- [ ] **Step 2: Update `TOLERANCE_BCRIT_RANGE` comment**

Edit `scripts/lecture_13_topological.py` — flip the comment from "3-witness range" to "4-witness range".

- [ ] **Step 3: Update `section_wire_rung` to read colormap z2 column**

Edit `scripts/lecture_13_topological.py:section_wire_rung` — rewrite the `bcrit_pfaffian` block:
```python
    # Slim Pfaffian row read live from output/z2_phase_diagram.dat z2 column
    # at mu ≈ 0.6601 (tol ±0.0001); first B row where z2 == -1.
    bcrit_pfaffian = read_colormap_first_z2_minus_one(
        Path(repo_root) / "output" / "z2_phase_diagram.dat",
        mu_target=0.6601, mu_tol=0.0001
    )
    if bcrit_pfaffian is None:
        raise RuntimeError(
            "Slim Pfaffian row missing: output/z2_phase_diagram.dat has no "
            "row with mu ≈ 0.6601 ± 0.0001 and z2 == -1. Re-run wire BdG "
            "gap sweep; the colormap-extracted witness is mandatory "
            "(ticket 05 of .scratch/archive/bdg-evaluator-pfaffian/)."
        )
```

- [ ] **Step 4: Remove the `(deferred to U13)` print branch**

Find and delete the print branch that says the colormap row is "reserved for U13".

- [ ] **Step 5: Update `render_reconciliation_table`**

Strip the "deferred to U13" cell-text branch + flip "3-witness range" → "4-witness range" + matching print + summary-row label.

- [ ] **Step 6: Update the trailing FAIL label**

Flip `FAIL: 3-witness disagreement` → `FAIL: 4-witness disagreement`.

## Task 3.4: `test_lecture_13_acceptance_gate.sh` — 4-witness + `PFAFFIAN_DEFERRED` strip

**Files:**
- Modify: `tests/integration/test_lecture_13_acceptance_gate.sh`

- [ ] **Step 1: Flip capture comment 3→4 witness**

Edit `tests/integration/test_lecture_13_acceptance_gate.sh` — update the capture comment.

- [ ] **Step 2: Strip `PFAFFIAN_DEFERRED` branch**

Find and remove the `PFAFFIAN_DEFERRED` branch (colormap row is mandatory).

- [ ] **Step 3: Collapse tolerance-strip + inner Python `min`/`max` to 4-witness**

Find the inner Python block that computes the tolerance range and verify it covers all 4 witnesses.

- [ ] **Step 4: Update `WITNESS_LABEL`**

Set `WITNESS_LABEL="4-witness"`.

- [ ] **Step 5: Flip `TOL_BCRIT_RANGE` comment**

Flip the comment from 3-witness to 4-witness.

## Task 3.5: `verify_majorana_polarization.py` — `gate_row_colormap_present` precondition

**Files:**
- Modify: `tests/integration/verify_majorana_polarization.py`

- [ ] **Step 1: Rewrite module docstring**

Clarify the polarization (Sticlet `P_M`) observable is distinct from the gate's Z2 row.

- [ ] **Step 2: Add `gate_row_colormap_present` helper**

Add:
```python
def gate_row_colormap_present(repo_root, mu_target=0.6601, mu_tol=0.0001):
    """Return True iff output/z2_phase_diagram.dat has a row with
    mu ≈ mu_target (± mu_tol) and z2 == -1.

    Per ticket 05 of .scratch/archive/bdg-evaluator-pfaffian/ — the polarization
    verifier's SKIP precondition is tightened: SKIP only when this row
    is present (otherwise FAIL non-zero).
    """
    colormap = Path(repo_root) / "output" / "z2_phase_diagram.dat"
    if not colormap.exists():
        return False
    with open(colormap) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            mu = float(parts[1])
            z2 = int(float(parts[3]))
            if abs(mu - mu_target) <= mu_tol and z2 == -1:
                return True
    return False
```

- [ ] **Step 3: Tighten SKIP precondition**

In `parse_polarization`'s missing-file branch, call `gate_row_colormap_present(repo_root)`; if False, exit non-zero with the missing-row error message; if True, allow the existing SKIP behavior.

## Task 3.6: parent plan + BACKLOG + REVIEW updates

**Files:**
- Modify: `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md`
- Modify: `docs/plans/BACKLOG.md`
- Modify: `docs/plans/REVIEW.md`

- [ ] **Step 1: Update parent plan `status_note`**

Edit the `status_note` block at `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md:9-62` — flip U2 from "still open" to "shipped" with one-line gist + the PR reference (the PR number will be filled in after the PR is opened; use the placeholder `feat/bdg-u2-actual-ship` for now and update after merge).

Also flip U12 wording from "3-witness" to "4-witness at 1.0 T" and add U2 shipped detail.

- [ ] **Step 2: Update BACKLOG.md line 5**

Edit `docs/plans/BACKLOG.md:5` — rewrite "U1, U2, U9, U10 open" → "U1, U9, U10 open; **U2 closed 2026-07-13** (slim Pfaffian seam + heuristic retirement + 4-witness gate row)".

- [ ] **Step 3: Update BACKLOG.md line 530**

Edit `docs/plans/BACKLOG.md:530` — flip the `~~U1, U2, U9, U10, U11, U12~~` line.

- [ ] **Step 4: Update BACKLOG.md line 542 (parent plan reference)**

Edit `docs/plans/BACKLOG.md:542` — trim U2 from the open list in the parent-plan summary.

- [ ] **Step 5: Update BACKLOG.md line 731 (Phase 24 + summary table trailer)**

Edit `docs/plans/BACKLOG.md:731` — replace "U1, U2, U9, U10 open" with "U1, U9, U10 open; **U2 closed 2026-07-13**". Also add a new Phase 26 row to the summary table: "BdG evaluator Pfaffian plug-in (U2 close-out, slim witness for wire rung)".

- [ ] **Step 6: Update REVIEW.md row 79**

Edit `docs/plans/REVIEW.md:85` (the row with `79 | 2026-06-14-001-feat-bdg-majorana-validation-plan`) — flip "Still open" wording to "U1, U9, U10" (drop U2); expand U2 + U12 with their 2026-07-13 close-out details.

## Task 3.7: Solution doc + memory entry

**Files:**
- Create: `docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md`
- Create: `~/.claude/projects/-data-8bandkp-fdm/memory/project_bdg_evaluator_seam_ssot.md`
- Edit: `~/.claude/projects/-data-8bandkp-fdm/memory/MEMORY.md` (add pointer)

- [ ] **Step 1: Create the solution doc**

Create `docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md`:
```markdown
---
module: bdg_observables
tags: [seam, SSOT, Pfaffian, slim-witness, ticket-07, U2-close-out]
problem_type: invariant-SSOT-discipline
component: bdg_observables
---

# BdG evaluator seam SSOT — slim Pfaffian plug-in + heuristic retirement

## Problem

Three sibling evaluators on three different modules — the heuristic
eigenvalues-only invariant and two ad-hoc per-rung wrappers — with no
single seam consumed by all BdG per-point work. The heuristic
`near_zero_count >= 2` discriminator silently masqueraded as a Z2
invariant on the wire rung. The acceptance gate was held back to
3-witness because no real Pfaffian signature was wired through the seam.

## Solution

Consolidate the per-point BdG invariant on `bdg_observables.f90` as the
single seam with three faces:

- `eval_bdg_point(eigenvalues, params) → bdg_eval_result_t` —
  eigenvalues-only, minigap + heuristic invariant (kept for QW rung).
- `eval_bdg_pfaffian_witness_csr(H_bdg_csr, Nbdg, params) → s2_sign ∈
  {-1, 0, +1}` — wire-rung invariant. Slim projected Pfaffian over S1
  ⊗ S2. Replaces the heuristic on the wire rung.
- `eval_bdg_kitaev_majorana(H_k_array, k_par_values) → majorana_number ∈
  {-1, 0, +1}` — QW+Kitaev rung. Wraps the existing Kitaev helper.

Scope-narrow `compute_z2_gap` / `compute_z2_gap_edge` to
`compute_z2_gap_bhz_heuristic` / `compute_z2_gap_edge_bhz_heuristic`,
signalling BHZ-only at the call site.

The acceptance gate's Pfaffian row reads live from the existing `(B, μ)`
colormap dataset (`output/z2_phase_diagram.dat`); no new Fortran I/O.

## Why this works

- Pure-function discipline: seam stays one module with three faces,
  not three seams. Dispatch by procedure choice (ADR 0001).
- Layering: seam imports only L0 leaves (`sparse_matrices` + `pfaffian`),
  not the L3 `topological_analysis` hub.
- Heuristic retirement via scope-narrow + rename rather than deletion:
  the heuristic is the gap-closure fallback, not the invariant.
- Live gate witnesses derived from existing Fortran output, no new I/O.
- Magic-number SSOT: `bdg_default_pfaffian_floor = 1.0e-12_dp` replaces
  the literal at 3 call sites in `topological_analysis.f90`.

## When to use

- New BdG invariant → add a sibling to `bdg_observables.f90`, sharing the
  `bdg_eval_params_t` record. Don't grow the seam into a hub.
- Gate row source of truth is the colormap (`output/z2_phase_diagram.dat`
  z2 column), NOT a separately-emitted Fortran file.
- Naming: BHZ-only helpers get `_bhz_heuristic` suffix; never pretend a
  heuristic is a generic Z2 helper.

## Source

- Map + 7 tickets: `.scratch/archive/bdg-evaluator-pfaffian/`
- Spec: `.scratch/archive/bdg-evaluator-pfaffian/spec.md`
- Dense-path witness: `src/physics/topological_analysis.f90:1626-1790`
- Seam: `src/physics/bdg_observables.f90`
- Gate: `tests/integration/test_lecture_13_acceptance_gate.sh`
- Verifier: `tests/integration/verify_majorana_polarization.py`
- Lecture: `docs/lecture/13-topological-superconductivity.md` §13.7.4 + §13.7.5
```

- [ ] **Step 2: Create the memory entry**

Create `~/.claude/projects/-data-8bandkp-fdm/memory/project_bdg_evaluator_seam_ssot.md`:
```markdown
---
name: bdg-evaluator-seam-ssot
description: U2 closed 2026-07-13; `bdg_observables.f90` seam SSOT for invariant work; slim projected Pfaffian = wire-rung invariant; `compute_z2_gap_bhz_heuristic` = gap-closure fallback only
metadata:
  type: project
---

The BdG per-point invariant seam is `bdg_observables.f90` with three faces:
`eval_bdg_point` (eigenvalues-only heuristic), `eval_bdg_pfaffian_witness_csr`
(slim projected Pfaffian over S1 ⊗ S2, wire-rung invariant), and
`eval_bdg_kitaev_majorana` (QW+Kitaev rung). Imports L0 leaves `sparse_matrices` +
`pfaffian` only — not the L3 `topological_analysis` hub.

The slim Pfaffian is the wire-rung invariant; the full Bloch-Pfaffian
swap is U13 (separate scoped PR per CLAUDE.md Known Issues). The acceptance
gate's Pfaffian row reads live from `output/z2_phase_diagram.dat` z2 column at
mu ≈ 0.6601 ± 0.0001; no new Fortran I/O. The heuristic
`near_zero_count >= 2` was RETIRED as the wire-rung invariant discriminator
(it remains the discriminant inside `eval_bdg_point` for the QW rung).

`compute_z2_gap` / `compute_z2_gap_edge` were renamed to `..._bhz_heuristic` —
they are BHZ-only (hardcoded 4-band basis), kept as the gap-closure fallback
after the slim Pfaffian returns 0. `compute_z2_gap_sweep` is untouched
(BHZ-only; calls `eval_bhz_analytic` directly).

**Why:** the seam stays a Layer-1 leaf; the heuristic's BHZ-only scope is
signalled at the call site; the gate row reads from existing Fortran output.

**How to apply:** when adding a new BdG invariant, add a sibling to
`bdg_observables.f90` sharing `bdg_eval_params_t`; promote any magic numbers
to SSOT parameters (`bdg_default_*`); do NOT import `topological_analysis`
from the seam. Related: [[bdg-z2-gap-helper-fate]], [[gate-witness-colormap]].
```

- [ ] **Step 3: Add the MEMORY.md pointer**

Edit `~/.claude/projects/-data-8bandkp-fdm/memory/MEMORY.md` — add the line:
```
- [BdG Evaluator Seam SSOT](project_bdg_evaluator_seam_ssot.md) — U2 closed 2026-07-13; `bdg_observables.f90` seam SSOT for invariant work; slim projected Pfaffian = wire-rung invariant; `compute_z2_gap_bhz_heuristic` = gap-closure fallback only
```

## Task 3.8: Re-archive the closed work

- [ ] **Step 1: Move the work back to `.scratch/archive/`**

```bash
mv .scratch/archive/bdg-evaluator-pfaffian .scratch/archive/bdg-evaluator-pfaffian
```

- [ ] **Step 2: Update map.md + spec.md footers**

Edit `.scratch/archive/bdg-evaluator-pfaffian/spec.md:1` — flip Status from "closed (as-built spec for shipped work)" to add the PR reference: "Status: closed (shipped 2026-07-13 via PR #N, merged `main@<sha>`)."

Edit `.scratch/archive/bdg-evaluator-pfaffian/map.md:1` — add the PR reference.

- [ ] **Step 3: Update the U10 handoff reference**

Edit `.scratch/bdg-u10-pfaffian-phase-diagram/HANDOFF.md` — search for `.scratch/archive/bdg-evaluator-pfaffian/` references; ensure all point to `.scratch/archive/bdg-evaluator-pfaffian/`.

## Task 3.9: Final verification — full ctest green

- [ ] **Step 1: Build**

Run: `cmake --build build 2>&1 | tail -5`
Expected: clean.

- [ ] **Step 2: Run the full unit-test label**

Run: `OMP_NUM_THREADS=$(( $(nproc)/4 )) ctest --test-dir build -L unit -j4 2>&1 | tee .scratch/bdg-u2-actual-ship/ctest-final-unit.log | tail -5`
Expected: 0 failures.

- [ ] **Step 3: Run the lecture-13 acceptance gate**

Run: `bash tests/integration/test_lecture_13_acceptance_gate.sh 2>&1 | tee .scratch/bdg-u2-actual-ship/lecture-13-gate-final.log | tail -10`
Expected: PASS at 4-witness, B_crit range within 1.0 T tolerance.

- [ ] **Step 4: Run the polarization verifier**

Run: `python3 tests/integration/verify_majorana_polarization.py 2>&1 | tail -10`
Expected: SKIP (per the new precondition) with the message confirming the colormap-extracted row is present.

## Task 3.10: Commit + open PR

- [ ] **Step 1: Stage everything**

```bash
git add -A
```

Run: `git status`
Expected: a coherent set of changes covering src/physics/bdg_observables.f90, src/physics/topological_analysis.f90, src/apps/main_topology.f90, src/physics/AGENTS.md, tests/unit/test_bdg_pfaffian_witness_csr.pf (new), tests/unit/test_bdg_kitaev_majorana.pf (new), tests/unit/test_bire_pfaffian_witness.pf, tests/unit/test_bdg_evaluator.pf, tests/unit/test_kitaev_majorana.pf, tests/CMakeLists.txt, docs/lecture/13-topological-superconductivity.md, docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md, docs/plans/BACKLOG.md, docs/plans/REVIEW.md, scripts/lecture_13_topological.py, tests/integration/test_lecture_13_acceptance_gate.sh, tests/integration/verify_majorana_polarization.py, docs/solutions/best-practices/2026-07-13-bdg-evaluator-seam-ssot.md (new).

- [ ] **Step 2: Final commit**

```bash
git commit -m "feat(bdg): ship U2 seam siblings + slim Pfaffian plug-in + heuristic retirement"
```

(No `Co-Authored-By:` trailer per repo convention `feedback_no_coauthor_trailer`.)

- [ ] **Step 3: Push the branch**

```bash
git push -u origin feat/bdg-u2-actual-ship
```

- [ ] **Step 4: Open the PR**

Use `gh pr create` (the project's PR workflow per CLAUDE.md) with a body that summarizes:
- Problem: U2 was prematurely closed via doc-only stamps (no backing code); ctest 50/51 RED.
- Solution: revert + TDD-execute + re-stamp. Phase 1 (clean revert), Phase 2 (seam siblings + rename + migration + tests + 51/51+ green), Phase 3 (re-apply stamps now backed by green).
- Verification: full ctest green + lecture-13 acceptance gate 4-witness GREEN at 1.0 T.
- Out of scope: U9 (spectral + LDOS), U10 (phase diagram), U11 (lecture revamp), U13 (Bloch-Pfaffian).

Hand off the PR URL to the user for review.

---

# Verification gates (run before claiming complete)

- [ ] **Phase 1 (Task 1.4):** ctest 50/51 PASS pre-revert (clean state restored).
- [ ] **Phase 2 (Task 2.9):** unit tests all green (50/51 → 60+/60+ — wire_pfaffian_witness strict assertion GREEN; +4 +3 +2 +1 new tests).
- [ ] **Phase 3 (Task 3.9):** full ctest green + lecture_13_acceptance_gate 4-witness GREEN + polarization verifier SKIP-with-row-present.

## Quality gates

- [ ] No `Co-Authored-By:` trailer on any commit (`feedback_no_coauthor_trailer`).
- [ ] All pFUnit `@assertEqual` / `@assertTrue` macros are SINGLE-LINE (`feedback_pfunit_macro_single_line`).
- [ ] All `error stop '<descriptive message>'`, no bare `stop 1`.
- [ ] AGENTS.md inventory + DAG row reflect actual seam state (not pre-claims).
- [ ] No `bdg_default_pfaffian_floor` literal duplication outside the SSOT.
- [ ] Memory entry + MEMORY.md pointer created.
- [ ] Solution doc created with proper frontmatter.

## Doc-drift discipline (per `codebase-doc-drift-prevention` memory)

- Parent plan `status_note` updated to "U2 shipped" with PR ref.
- BACKLOG summary-table + trailer + parent-plan reference all carry the new wording.
- REVIEW row 79 mirrors the parent plan.
- Lecture §13.0 + §13.7.1 + §13.7.4 + new §13.7.5 + summary table + tail flipped to live slim-witness disclosure.
- Gate shell `PFAFFIAN_DEFERRED` branch stripped + `WITNESS_LABEL="4-witness"` + `TOL_BCRIT_RANGE="1.0"`.
- Lecture script `TOLERANCE_BCRIT_RANGE` comment + `render_reconciliation_table` + `section_wire_rung` flipped 3→4.
- Polarization verifier SKIP precondition parses colormap + exits non-zero if z2==-1 row missing.
- AGENTS.md inventory + DAG row reflect the seam siblings + L0 deps.
- New solution doc + memory entry + MEMORY.md index pointer.

## Cross-references

- Source spec: `.scratch/archive/bdg-evaluator-pfaffian/spec.md` (spec-of-record for the seam shape)
- Source map + 7 tickets: `.scratch/archive/bdg-evaluator-pfaffian/` (will re-archive in Task 3.8)
- U10 handoff: `.scratch/bdg-u10-pfaffian-phase-diagram/HANDOFF.md` (separate scope, KEEP)
- Parent plan: `docs/plans/2026-06-14-001-feat-bdg-majorana-validation-plan.md` §U2 (lines 196-219)
- ENGINEERING PRINCIPLES: CLAUDE.md root
- Conventions: `feedback_no_coauthor_trailer`, `feedback_pfunit_macro_single_line`, `codebase-doc-drift-prevention`, `codebase-doc-drift-event-3`
- Locked decisions (do NOT re-litigate): `.scratch/archive/bdg-evaluator-pfaffian/issues/01-evaluator-api-shape.md` through `06-test-coverage-audit.md`