---
title: "CSR test infrastructure: stop 1, OOB access, tautological assertions"
date: 2026-05-09
category: docs/solutions/test-failures
module: tests/support
problem_type: test_failure
component: testing_framework
symptoms:
  - assert_csr_structural_invariants aborted entire test process via stop 1 instead of reporting failures through pFUnit
  - csr_hermitian_error segfaulted on non-square matrices due to OOB dense array indexing
  - kpterms_2d cleanup freed only 15 of 17 elements across 4 test routines
  - R12 verification rung contained a tautological overlap check that always passed
  - No test coverage for csr_add with disjoint (non-overlapping) sparsity patterns
root_cause: logic_error
resolution_type: test_fix
severity: medium
related_components:
  - csr_test_helpers
  - krylov_helpers
  - test_hamiltonian_2d
  - verify_8band_rung3_qw
tags: [csr, pfunit, bounds-checking, test-infrastructure, verification-ladder, fortran]
---

# CSR test infrastructure: stop 1, OOB access, tautological assertions

## Problem

Six code review findings in CSR structural testing infrastructure: a test helper using `stop 1` instead of returning pass/fail, an out-of-bounds access in a Hermitian-check helper for non-square matrices, a missing dimension guard in Krylov comparison, a hardcoded loop bound in 2D kpterms cleanup (freeing only 15 of 17 elements), a tautological material overlap assertion in the verification ladder, and a missing disjoint-sparsity test case.

These issues were discovered during a multi-agent code review of the CSR structural testing plan execution. The test infrastructure was originally built as prevention for topological index-logic bugs documented in [topological-magnetic-index-logic-errors](../logic-errors/topological-magnetic-index-logic-errors-2026-05-08.md) (session history).

## Symptoms

- `assert_csr_structural_invariants` calling `stop 1` killed the entire test process instead of reporting a single assertion failure through pFUnit, preventing identification of which specific invariant was violated
- `csr_hermitian_error` crashed with a segfault when given a non-square CSR matrix because `dense(j,i)` used `j` up to `ncols` as a row index in a `dense(nrows, ncols)` array
- `krylov_compare` silently produced wrong results or crashed on dimension mismatch between reference and computed vectors
- `kpterms_2d` cleanup loop `do ij = 1, 15` silently leaked 2 of 17 CSR matrices in 4 test subroutines
- R12 verification rung compared `material_overlap_mev` (computed from hardcoded GaSbW/InAsW parameters) against `expected_overlap_mev = 142.0` (same hardcoded constants) -- the check could never fail

## What Didn't Work

- The `stop 1` approach seemed acceptable because it did halt on invariant violations, but it prevented pFUnit from generating proper test reports and masked which specific check failed. Individual `check_*` functions already returned logical -- only the top-level aggregator used `stop 1` (session history: brainstorm phase identified this gap).
- The tautological overlap check appeared valid because it passed consistently, but the pass was guaranteed by construction. A useful regression test must derive expected and actual values from independent sources.
- Hardcoded `15` in cleanup loops worked for the current array size but broke silently when `confinementInitialization_2d` was extended to produce 17 kpterms.

## Solution

**Fix 1: Convert subroutine to logical function.** All 28 callers updated.

```fortran
! Before
subroutine assert_csr_structural_invariants(mat, require_diagonal)
  if (.not. check_rowptr_monotonic(mat)) stop 1
  ...

! After
function assert_csr_structural_invariants(mat, require_diagonal) result(ok)
  logical :: ok
  ok = .true.
  if (.not. check_rowptr_monotonic(mat)) then; ok = .false.; return; end if
  ...
```

Callers changed from `call assert_csr_structural_invariants(A, .true.)` to `@assertTrue(assert_csr_structural_invariants(A, .true.), message="...")`.

**Fix 2: Clamp Hermitian error loops.**

```fortran
! Before
do j = 1, mat%ncols
  do i = 1, mat%nrows
    max_err = max(max_err, abs(dense(i,j) - conjg(dense(j,i))))

! After
do j = 1, min(mat%nrows, mat%ncols)
  do i = 1, min(mat%nrows, mat%ncols)
    max_err = max(max_err, abs(dense(i,j) - conjg(dense(j,i))))
```

**Fix 3: Dimension guard in krylov_compare.**

```fortran
if (size(reference, 1) /= n .or. size(reference, 2) /= k) then
  pass = .false.
  if (present(failing_iteration)) failing_iteration = 0
  if (present(max_divergence)) max_divergence = huge(1.0_dp)
  return
end if
```

**Fix 4: Use `size()` instead of magic number.** Applied at 4 sites in `test_hamiltonian_2d.pf`.

```fortran
! Before
do ij = 1, 15
  call csr_free(kpterms_2d(ij))

! After
do ij = 1, size(kpterms_2d)
  call csr_free(kpterms_2d(ij))
```

**Fix 5: Removed tautological R12 overlap check.** The self-referential 17-line block was deleted. The meaningful physics check (VB state above InAs CB edge) was preserved.

**Fix 6: Added csr_add disjoint sparsity test.** New test covering union of two CSR matrices with non-overlapping patterns (A = diagonal, B = off-diagonal).

## Why This Works

- **Fix 1:** pFUnit's `@assertTrue` captures failure context (source line, custom message) and continues running other tests. The `stop 1` bypassed the entire test framework. The early-return pattern reports the first violation rather than cascading into subsequent checks on an already-invalid matrix.
- **Fix 2:** For Hermiticity checking, only the square submatrix intersection matters. The dense array has shape `(nrows, ncols)`, so accessing `dense(j,i)` with `j > nrows` is always out-of-bounds when `ncols > nrows`.
- **Fix 3:** The early return with sentinel values (`huge(1.0_dp)`) makes dimension mismatches immediately visible rather than producing silent wrong answers deep in the comparison loop.
- **Fix 4:** The array size is the single source of truth. Hardcoded literals work only as long as the array dimension stays fixed.
- **Fix 5:** A test that cannot fail provides zero coverage. Removing the tautology eliminates false confidence.
- **Fix 6:** Disjoint sparsity is a boundary case for `csr_add`'s index-merge path. Without it, only overlapping patterns were tested.

## Prevention

- **Never use `stop` in test helpers.** All assertion helpers should return a logical result that pFUnit reports through `@assertTrue`/`@assertFalse`. The `stop 1` pattern is appropriate for production error handling but defeats test frameworks.
- **Use `size(array)` instead of hardcoded integers** for loop bounds over allocatable arrays. Magic numbers become stale the moment the array dimension changes.
- **Guard array dimensions before indexing in generic helpers.** Any routine that accepts arrays and indexes them should validate dimensions early, especially when handling multiple matrix shapes (square, rectangular, tall, wide).
- **Avoid tautological assertions.** When writing a regression test, verify that the expected and actual values derive from independent sources. If both come from the same constants, the test is vacuous.
- **Test boundary cases for sparse operations.** Include cases with disjoint sparsity, empty rows, and single-entry matrices alongside typical overlapping patterns.

## Related Issues

- [topological-magnetic-index-logic-errors](../logic-errors/topological-magnetic-index-logic-errors-2026-05-08.md) -- the originating bugs that motivated the CSR structural testing infrastructure
- Plan: `docs/plans/archive/2026-05-09-001-feat-csr-structure-testing-plan.md`
- Commit: `7254c29` on `feature/bdg-topological-superconductivity`
