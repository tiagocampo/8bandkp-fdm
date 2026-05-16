---
date: 2026-05-09
topic: csr-structure-testing
---

# CSR Structure Testing

## Summary

Two-layer CSR structure testing: reusable property-based invariant fixtures that assert structural validity of any `csr_matrix` after build/transform operations, plus Krylov-action snapshot tests that verify full Hamiltonian assembly correctness across all CSR-using code paths (wire, topology, optics, g-factor wire, SC loop) by applying fixed SpMV chains and comparing against committed reference vectors.

---

## Problem Frame

CSR index-arithmetic errors in Hamiltonian assembly can corrupt the sparse structure while preserving eigenvalues — the Hamiltonian is a different matrix but with the same spectrum. Of the five topological/magnetic bugs documented in `docs/solutions/logic-errors/`, one (BHZ same-site hopping) was a genuine CSR structural bug invisible to eigenvalue tests; the others were Hermiticity breaks, OOB crashes, and loop logic errors caught by different test mechanisms. The proposed testing complements rather than replaces Hermiticity checks and crash-detecting integration tests. Current testing validates CSR operations via eigenvalues and SpMV correctness, but never independently asserts that the assembled Hamiltonian has the correct sparse structure. This gap means future CSR index-arithmetic regressions can pass all existing tests.

---

## Requirements

**Property-based structural invariants**

- R1. Reusable assertion fixture that validates CSR structural invariants on any `csr_matrix`: `rowptr` monotonic and 1-based, `colind` sorted within each row, `nnz = rowptr(nrows+1) - 1`, all `colind` entries in `[1, ncols]`. Diagonal entry in every row is required only for square matrices — first-derivative and cross-derivative operators legitimately lack diagonal entries.
- R2. Dedicated unit tests for currently untested CSR operations that include structural invariant assertions after the operation: `csr_conjugate_transpose`, `csr_conjugate_transpose_to_preallocated`, `csr_find_diag_positions`, `csr_clone_structure`, `csr_add` with overlapping but non-identical sparsity patterns (existing `test_csr_add_basic` covers same-pattern diagonal only), `kron_dense_eye`, `kron_dense_dense_1d`, and `csr_scale`.
- R3. Existing CSR operations that already have tests (`csr_build_from_coo`, `csr_build_from_coo_cached`, `kron_dense_dense`, `kron_eye_dense`, `csr_spmv`, `csr_apply_variable_coeff`) get structural invariant assertions added to their existing test wrappers or via new wrapper tests.

**Krylov-action snapshots**

- R4. For each CSR-using code path — wire, topology/QHE, topology/BHZ, topology/BdG, optics wire mode, g-factor wire mode, SC loop — a canonical small-grid config that builds the full Hamiltonian as CSR and applies a fixed-seed SpMV chain (k iterations, deterministic). For the SC loop path, snapshot the initial Hamiltonian before any SC iteration to test CSR assembly independently of SC convergence behavior.
- R5. Reference vector snapshots committed alongside the tests, with a documented one-command regeneration procedure for when Hamiltonian construction intentionally changes.
- R6. Comparison against reference within a specified tolerance, with failure output identifying which iteration diverged and the magnitude of divergence.

---

## Success Criteria

- Any structural CSR bug in Hamiltonian construction (wrong `colind`, missing diagonal, corrupted `rowptr`, off-by-one in indexing) is caught by at least one test layer — either property-based invariant on the CSR primitive or Krylov snapshot on the full assembly.
- The regeneration workflow for Krylov reference data is documented and one-command.
- All existing tests continue to pass without modification.

---

## Scope Boundaries

- Golden CSR structure files (row pointers, column indices, element values as reference data) — rejected in favor of Krylov snapshots which are more sensitive and less brittle under intentional sparsity pattern changes.
- Performance benchmarking of CSR operations.
- CSR testing for dense-solver (non-CSR) code paths.
- Changes to production code — this is test-only.
- Richardson extrapolation convergence testing (separate ideation item).
- Refactoring existing CSR tests to use the new invariant fixtures (new tests only).

---

## Key Decisions

- **Krylov snapshots over golden CSR files:** An SpMV chain at fixed seed propagates structural errors across the matrix's connected components, making it highly sensitive to index-arithmetic bugs that redirect nonzero entries. Detection depends on the seed vector having nonzero support in the affected row/column subspace — a random seed covers this with high probability. Reference vectors require regeneration whenever Hamiltonian values or sparsity pattern change (both affect SpMV output). The advantage over golden CSR files is that Krylov snapshots test the full assembly pipeline end-to-end, and the reference data is a compact vector per iteration rather than the full CSR structure.
- **Property-based invariants as reusable fixtures:** Any new `csr_matrix` operation or test gets structural validation with zero additional effort. The fixture catches the class of bugs (corrupted indices, missing diagonals) that are invisible to value-only tests.
- **Both layers rather than one:** CSR primitives and Hamiltonian assembly are different bug surfaces. The topological bugs were in assembly logic, not CSR infrastructure, but infrastructure bugs would be equally invisible to eigenvalue tests.

---

## Dependencies / Assumptions

- pFUnit supports the assertion patterns needed for structural invariant checks (array comparison, index validation).
- Canonical configs for Krylov snapshots are small enough (grid size, number of materials) to run within normal test execution time (under 30 seconds each).
- The existing `csr_to_dense` helper in `test_hamiltonian_2d.pf` can serve as a reference pattern for the property-based fixture, but the fixture itself should be in a shared location accessible to all CSR test files.
- A shared test-support library (e.g., `8bandkp_test_support`) or a Fortran include file (e.g., `tests/unit/csr_invariant.inc`) is needed for the reusable fixture. pFUnit test modules cannot import from other test modules — only from production code libraries. The existing `LINK_LIBRARIES 8bandkp_common` pattern in `tests/CMakeLists.txt` supports either approach with minimal change.

---

## Outstanding Questions

### Deferred to Planning

- [Affects R4, R5][Technical] What SpMV chain length (k) provides sufficient sensitivity without excessive test runtime? Planning should experiment with k=3..10 on representative configs.
- [Affects R4][Technical] What canonical configs best exercise each CSR path while staying small enough for fast test runs? Planning should survey existing regression configs and select minimal representatives.
- [Affects R1][Technical] Where should the reusable structural invariant fixture live? Options: (a) a `tests/unit/test_csr_helpers.f90` compiled as a module and linked via `LINK_LIBRARIES` in CMake, or (b) a `tests/unit/csr_invariant.inc` include file. Existing codebase uses option (a) for production modules; no existing test uses shared includes.
