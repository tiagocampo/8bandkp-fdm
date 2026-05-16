---
title: "feat: CSR structure testing — invariant fixtures and Krylov snapshots"
type: feat
status: active
date: 2026-05-09
origin: docs/brainstorms/2026-05-09-csr-structure-testing-requirements.md
---

# feat: CSR structure testing — invariant fixtures and Krylov snapshots

## Summary

Two-layer CSR testing infrastructure: a reusable structural-invariant fixture in a new `8bandkp_test_support` library consumed by all pFUnit test executables, plus Krylov-action snapshot tests that apply fixed-seed SpMV chains to Hamiltonians assembled via each CSR code path and compare against committed reference vectors. Six implementation units cover the fixture library, tests for 8 untested CSR operations, refactoring existing tests to use the shared fixture, invariant assertions on existing tests, Krylov infrastructure, and per-path snapshot tests. No production code is modified.

---

## Problem Frame

CSR index-arithmetic errors in Hamiltonian assembly can corrupt the sparse structure while preserving eigenvalues — the matrix is different but has the same spectrum. Of the five topological/magnetic bugs documented in `docs/solutions/logic-errors/`, the BHZ same-site hopping bug was a genuine CSR structural error invisible to eigenvalue tests. Current testing validates CSR operations via eigenvalues and SpMV correctness, but never independently asserts that the assembled Hamiltonian has the correct sparse structure. Eight CSR operations have zero dedicated tests, and structural invariants (rowptr monotonicity, colind sorting, nnz consistency) are never checked. This plan closes that gap with two complementary testing layers.

---

## Requirements

- R1. Reusable assertion fixture validating CSR structural invariants on any `csr_matrix`: rowptr monotonic and 1-based, colind sorted within each row, nnz = rowptr(nrows+1) - 1, all colind entries in [1, ncols]. Diagonal entry required only for square matrices.
- R2. Dedicated unit tests for 8 untested CSR operations that include structural invariant assertions: `csr_conjugate_transpose`, `csr_conjugate_transpose_to_preallocated`, `csr_find_diag_positions`, `csr_clone_structure`, `csr_add` with overlapping patterns, `kron_dense_eye`, `kron_dense_dense_1d`, `csr_scale`.
- R3. Existing CSR operations that already have tests (`csr_build_from_coo`, `csr_build_from_coo_cached`, `kron_dense_dense`, `kron_eye_dense`, `csr_spmv`, `csr_apply_variable_coeff`) get structural invariant assertions added to their existing tests.
- R4. For each CSR-using code path — wire, topology/QHE, topology/BHZ, topology/BdG, optics wire mode, g-factor wire mode, SC loop — a canonical small-grid config that builds the full Hamiltonian as CSR and applies a fixed-seed SpMV chain.
- R5. Reference vector snapshots committed alongside the tests, with a documented one-command regeneration procedure for when Hamiltonian construction intentionally changes.
- R6. Comparison against reference within a specified tolerance, with failure output identifying which iteration diverged and the magnitude of divergence.

---

## Scope Boundaries

- Golden CSR structure files (rejected in favor of Krylov snapshots)
- Performance benchmarking of CSR operations
- CSR testing for dense-solver code paths
- Changes to production code — this is test-only
- Richardson extrapolation convergence testing

### Deferred to Follow-Up Work

- Remove `kron_dense_dense_1d` dead code (tested here for coverage; removal is separate)
- Extend Krylov snapshots to test velocity matrix assembly independently (covered by structural invariant assertions on CSR operations used to build them)

---

## Context & Research

### Relevant Code and Patterns

- CSR type and all operations: `src/math/sparse_matrices.f90` — type at lines 37-48, operations at lines 61-1150+
- Wire Hamiltonian assembly (primary CSR consumer): `src/physics/hamiltonian_wire.f90`
- BdG Hamiltonian assembly: `src/physics/bdg_hamiltonian.f90`
- Topological analysis assembly (QWZ, BHZ): `src/physics/topological_analysis.f90`
- Confinement init (Kronecker products + variable coeff): `src/physics/confinement_init.f90`
- Existing CSR unit tests: `tests/unit/test_csr_spmv.pf` (6 tests), `tests/unit/test_hamiltonian_2d.pf` (21 tests with local helpers: `csr_to_dense`, `assert_csr_interior_symmetric`, `assert_csr_interior_hermitian`, `assert_csr_hermitian`)
- Test build configuration: `tests/CMakeLists.txt` — `add_pfunit_ctest(... LINK_LIBRARIES 8bandkp_common)` pattern
- Library build: `src/CMakeLists.txt` — `8bandkp_common` static library

### Institutional Learnings

- Five index-arithmetic bugs in topological/magnetic modules; BHZ same-site hopping was a CSR structural bug invisible to eigenvalue tests. Prevention: structural CSR entry tests verifying specific entries exist or are absent.
- 8-band verification ladder (Rung 4: Wire + sparse) tests CSR assembly via dense-sparse consistency but not structural invariants independently.
- Build with `-fcheck=all` during development catches OOB at point of occurrence.

### External References

None — strong local testing patterns, well-understood domain.

---

## Key Technical Decisions

- **Separate `8bandkp_test_support` static library** over polluting `8bandkp_common`: keeps production code clean. CMake adds the library and links it to all test executables via the existing `LINK_LIBRARIES` pattern. pFUnit test modules cannot import from each other — only from linked libraries — so this is the only clean way to share code across test executables.
- **SpMV chain length k=5**: for 3x3 grids (graph diameter ~4), five iterations propagate errors through all connected components. Each iteration is O(nnz), negligible for small test grids. Chain length is tuneable in the Krylov helper.
- **Programmatic config construction** over `input.cfg` files: tests call `confinementInitialization_2d`, `ZB8bandGeneralized`, `build_bdg_hamiltonian`, etc. directly, following the pattern in `test_hamiltonian_2d.pf`. No file I/O, fully deterministic.
- **Fortran module constants** for reference vectors over external files: eliminates I/O complexity, file-path issues, and format parsing. Reference data lives in `tests/support/krylov_reference_data.f90` as parameter arrays. Regeneration overwrites this file.
- **Test `kron_dense_dense_1d` despite dead code**: coverage value justifies the small effort. Removal deferred to follow-up.

---

## Open Questions

### Resolved During Planning

- Where should the reusable fixture live? Separate `8bandkp_test_support` library — see Key Technical Decisions.
- What SpMV chain length? k=5 default — see Key Technical Decisions.
- How to build configs for Krylov snapshots? Programmatically in test code — see Key Technical Decisions.
- How to store reference vectors? Fortran module constants — see Key Technical Decisions.

### Deferred to Implementation

- Exact small-grid parameters for each code path's Krylov snapshot (grid size, materials, physics parameters): implementation discovers the smallest config that exercises each path within 30 seconds per test.
- Exact tolerance for Krylov reference comparison: depends on accumulated floating-point error across 5 SpMV iterations. Start with 1e-12 and relax if needed.
- Which local helpers in other test files (beyond `test_hamiltonian_2d.pf`) are candidates for migration to shared library: implementation surveys all `.pf` files during U3.

---

## Implementation Units

- U1. **Test-support library and structural invariant fixture**

**Goal:** Create the shared test infrastructure — a CMake library target and a Fortran module providing reusable CSR structural invariant assertions and common test helpers.

**Requirements:** R1

**Dependencies:** None

**Files:**
- Create: `tests/support/csr_test_helpers.f90`
- Create: `tests/support/CMakeLists.txt`
- Modify: `tests/CMakeLists.txt`
- Test: `tests/unit/test_csr_helpers.pf`

**Approach:**
- Create `tests/support/` directory with `CMakeLists.txt` that builds `8bandkp_test_support` from `csr_test_helpers.f90`
- The `csr_test_helpers` module provides `assert_csr_structural_invariants(csr, require_diagonal)` as the primary fixture, plus individual boolean check functions for each invariant (rowptr monotonic, colind sorted, colind bounds, nnz consistency, diagonal presence)
- Move commonly-needed helpers here: `csr_to_dense` (from `test_hamiltonian_2d.pf`), `assert_csr_hermitian` (from `test_hamiltonian_2d.pf`)
- Update `tests/CMakeLists.txt` to add subdirectory and link `8bandkp_test_support` to all existing test executables
- Write `test_csr_helpers.pf` validating the fixture against CSR matrices with known structural violations

**Patterns to follow:**
- Module structure: `private` default with explicit `public` exports (per CLAUDE.md conventions)
- CMake: `add_pfunit_ctest(... LINK_LIBRARIES 8bandkp_common 8bandkp_test_support)`
- Existing local helpers: `csr_to_dense` and `assert_csr_hermitian` in `test_hamiltonian_2d.pf`

**Test scenarios:**
- Happy path: valid 3x3 CSR built via `csr_build_from_coo` passes all invariants
- Happy path: valid CSR with diagonal present passes with `require_diagonal=.true.`
- Edge case: empty CSR (0 nonzeros, nrows=3, ncols=3) passes all invariants
- Edge case: rectangular CSR (nrows /= ncols) passes without diagonal requirement
- Edge case: CSR with empty rows passes invariants (rowptr(i) == rowptr(i+1))
- Error path: rowptr not monotonic — fixture detects and reports
- Error path: colind unsorted within a row — fixture detects
- Error path: colind entry outside [1, ncols] — fixture detects
- Error path: nnz inconsistent with rowptr — fixture detects
- Error path: missing diagonal when `require_diagonal=.true.` — fixture detects

**Verification:**
- `test_csr_helpers.pf` passes with all 10 test cases
- All existing tests continue to pass (no regression from CMake changes)

---

- U2. **Unit tests for untested CSR operations**

**Goal:** Dedicated unit tests for the 8 CSR operations with zero test coverage, each including structural invariant assertions after the operation.

**Requirements:** R1, R2

**Dependencies:** U1

**Files:**
- Create: `tests/unit/test_csr_structural.pf`
- Modify: `tests/CMakeLists.txt`

**Approach:**
- Single test file with one `@test` subroutine per untested operation
- Each test: build known input CSR(s), perform the operation, assert output correctness (values), then call `assert_csr_structural_invariants` on the result
- `csr_find_diag_positions` gets an additional error-path test for missing diagonal
- `kron_dense_dense_1d` flagged as dead code in comment; tested for coverage

**Patterns to follow:**
- Test structure: `test_csr_spmv.pf` (one operation per `@test`, cleanup with `call csr%free()`)
- CSR construction: `csr_build_from_coo` for test matrices
- Invariant assertion: `call assert_csr_structural_invariants(result, .true.)` from U1

**Test scenarios:**
- `csr_conjugate_transpose`: known 3x3 complex matrix, A^H matches hand-computed result, structural invariants hold
- `csr_conjugate_transpose_to_preallocated`: fast-path result matches slow path, structural invariants hold
- `csr_find_diag_positions`: known 4x4 matrix, returned positions match expected indices
- `csr_find_diag_positions` error path: non-square or missing diagonal triggers error
- `csr_clone_structure`: same sparsity pattern, all values zero, structural invariants hold
- `csr_add` overlapping: union sparsity, correct values at overlap, structural invariants hold
- `kron_dense_eye`: known 2x3 dense x I_2, correct structure and values, structural invariants hold
- `kron_dense_dense_1d`: known 1D vectors, outer product matches hand-computed, structural invariants hold
- `csr_scale`: values scaled correctly, structure unchanged, structural invariants hold

**Verification:**
- All 9 test cases in `test_csr_structural.pf` pass (8 operations; `csr_find_diag_positions` has happy-path and error-path tests)
- Structural invariant assertions hold on every operation output

---

- U3. **Refactor existing tests to use shared fixture**

**Goal:** Replace local helper subroutines in existing test files with calls to the shared `csr_test_helpers` module, eliminating code duplication across test executables.

**Requirements:** R1 (extended scope)

**Dependencies:** U1

**Files:**
- Modify: `tests/unit/test_hamiltonian_2d.pf`
- Modify: `tests/unit/test_bdg_hamiltonian.pf` (has `csr_to_dense_local` duplicate at line 347)
- Modify: other `.pf` files with local CSR helpers (implementation surveys remaining test files)
- Modify: `tests/CMakeLists.txt` (ensure `8bandkp_test_support` linked to all targets)

**Approach:**
- Local helpers in `test_hamiltonian_2d.pf` (`csr_to_dense`, `assert_csr_interior_symmetric`, `assert_csr_interior_hermitian`, `assert_csr_hermitian`) already moved to `csr_test_helpers` in U1
- Replace local definitions and calls with `use csr_test_helpers, only: ...` imports
- Survey all other `.pf` files for similar local CSR helpers and consolidate into the shared module
- Ensure all test executables link `8bandkp_test_support` in `tests/CMakeLists.txt`

**Patterns to follow:**
- Existing local helper signatures and behavior preserved exactly
- Module structure: `private` default with explicit `public` exports

**Test scenarios:**
- All existing tests pass identically after refactoring (zero behavioral change)
- `test_hamiltonian_2d.pf` no longer contains local CSR helper definitions
- All test executables link `8bandkp_test_support` in CMake

**Verification:**
- `ctest --test-dir build -L unit` green — no test behavior changes
- Grep confirms no local CSR helper duplicates remain in `.pf` files

---

- U4. **Invariant assertions on existing CSR tests**

**Goal:** Add structural invariant assertions to existing tests for the 6 CSR operations that already have test coverage but lack invariant checks.

**Requirements:** R1, R3

**Dependencies:** U1 (hard). U3 is recommended as workflow ordering (avoid modifying `test_hamiltonian_2d.pf` twice in succession), not a technical dependency — U4's `assert_csr_structural_invariants` calls come from U1, not U3.

**Files:**
- Modify: `tests/unit/test_csr_spmv.pf`
- Modify: `tests/unit/test_hamiltonian_2d.pf`

**Approach:**
- After each CSR build/transform operation in existing tests, add `call assert_csr_structural_invariants(csr)` (no diagonal requirement unless the test context is a square matrix)
- Operations targeted per R3: `csr_build_from_coo` (7 uses in `test_csr_spmv.pf`), `csr_build_from_coo_cached` (1 use), `kron_dense_dense` (3 uses in `test_hamiltonian_2d.pf`), `kron_eye_dense` (1 use), `csr_apply_variable_coeff` (2 uses)
- `csr_spmv` does not modify CSR structure — for the 6 SpMV tests in `test_csr_spmv.pf`, add invariant assertions on the CSR matrix before SpMV (verifying input structure is valid), not after

**Patterns to follow:**
- Insert `call assert_csr_structural_invariants(csr)` immediately after matrix construction, before value assertions

**Test scenarios:**
- All existing tests pass unchanged with the added invariant assertions
- Structural invariants confirmed for every CSR built in existing tests

**Verification:**
- `ctest --test-dir build -L unit` green — no existing test behavior changes
- Purely additive assertions, no test logic modification

---

- U5. **Krylov snapshot infrastructure**

**Goal:** Build the reusable SpMV chain helper, reference data module, comparison routine, and one-command regeneration procedure.

**Requirements:** R5, R6

**Dependencies:** None (independent of U1-U4)

**Files:**
- Create: `tests/support/krylov_helpers.f90`
- Create: `tests/support/krylov_reference_data.f90`
- Create: `tests/support/regenerate_krylov_references.f90`
- Modify: `tests/support/CMakeLists.txt`

**Approach:**
- `krylov_helpers` module provides `krylov_chain(csr, seed, k, vectors)` — applies k SpMV iterations from seed vector, stores all intermediate vectors; and `krylov_compare(vectors, reference, tolerance, failing_iteration, max_divergence)` — compares against reference, reports failure details
- `krylov_reference_data` module stores committed reference vectors as `complex(dp), parameter` arrays — one named constant per code path and iteration (e.g., `wire_krylov_k1`, `bdg_krylov_k5`)
- Regeneration program: standalone Fortran executable in `tests/support/` that builds all Hamiltonians, runs SpMV chains, writes `krylov_reference_data.f90` with updated parameter arrays. Added as a CMake custom target in `tests/support/CMakeLists.txt`. Invoked as `cmake --build build --target regenerate_krylov_references`
- **Naming note:** "Krylov chain" is shorthand for "fixed-seed SpMV iteration chain" throughout — the implementation is repeated SpMV (power iteration), not a Krylov subspace method. The module/function names retain "krylov" for brevity.

**Patterns to follow:**
- SpMV call: `call csr_spmv(A, x, y, ONE, ZERO)` (from production code in `gfactor_functions.f90`, `optical_spectra.f90`)
- Parameter array: `complex(dp), parameter :: name(n) = [(...)]`

**Test scenarios:**
- Happy path: Krylov chain on known 3x3 identity produces expected vectors at each iteration
- Happy path: comparison at tolerance 1e-12 passes for exact-match reference
- Error path: comparison detects divergence at a specific iteration and reports magnitude
- Error path: divergence in final iteration only — earlier iterations pass
- Integration: regeneration program produces identical reference data on two consecutive runs (determinism)

**Verification:**
- Krylov chain helper produces deterministic output for known CSR
- Reference comparison correctly passes/fails with informative output
- Regeneration program runs end-to-end and produces valid Fortran source

---

- U6. **Krylov snapshot tests per code path**

**Goal:** Krylov snapshot tests for all 7 CSR-using code paths, each building the Hamiltonian programmatically with a minimal config and comparing SpMV chain output against committed reference.

**Requirements:** R4, R5, R6

**Dependencies:** U5

**Files:**
- Create: `tests/unit/test_krylov_snapshots.pf`
- Modify: `tests/CMakeLists.txt`
- Modify: `tests/support/krylov_reference_data.f90`

**Approach:**
- One `@test` subroutine per code path:
  1. Programmatically construct minimal config (grid, materials, physics parameters)
  2. Call the appropriate Hamiltonian assembly routine
  3. Run `krylov_chain` with deterministic seed (e.g., all `cmplx(1.0, 0.0)`) and k=5
  4. Compare against committed reference via `krylov_compare`
- Code paths and assembly routines:
  - **Wire**: `ZB8bandGeneralized` — minimal GaAs wire (3x3 grid)
  - **Wire with magnetic field (Peierls phase)**: `ZB8bandGeneralized` with `ExternalField` magnetic field enabled — exercises the Peierls phase branch that applies position-dependent phases to CSR hopping terms (documented OOB bug source)
  - **BHZ**: `build_bhz_wire_hamiltonian` — BHZ wire in topological phase
  - **BdG**: `build_bdg_hamiltonian` — s-wave pairing + Zeeman splitting
  - **SC loop**: `ZB8bandGeneralized` with SC-enabled config — snapshot initial Hamiltonian before SC iteration
  - **Optics wire**: `ZB8bandGeneralized` — wire Hamiltonian with optics config
  - **G-factor wire**: `ZB8bandGeneralized` — wire Hamiltonian with g-factor config
- Reference vectors generated by running regeneration program once and committing the output

**Execution note:** Generate reference data first (run regeneration program), then write comparison tests. Reference data must be committed before tests can pass.

**Patterns to follow:**
- Programmatic config: `test_hamiltonian_2d.pf` pattern (direct calls to initialization routines)
- Krylov chain: `krylov_helpers` module from U5

**Test scenarios:**
- Happy path: wire Hamiltonian Krylov snapshot matches reference
- Happy path: wire-with-magnetic-field Hamiltonian Krylov snapshot matches reference (Peierls phase)
- Happy path: BHZ wire Hamiltonian Krylov snapshot matches reference
- Happy path: BdG Hamiltonian Krylov snapshot matches reference
- Happy path: SC loop initial Hamiltonian Krylov snapshot matches reference
- Happy path: optics wire mode Hamiltonian Krylov snapshot matches reference
- Happy path: g-factor wire mode Hamiltonian Krylov snapshot matches reference
- Edge case: tolerance relaxation from 1e-12 to 1e-8 still passes (reference is not at noise floor)
- Integration: deliberately corrupted colind entries detected by Krylov divergence

**Verification:**
- All 7 code path tests pass against committed reference
- `ctest --test-dir build -L unit` green
- Deliberate Hamiltonian corruption caught by at least one Krylov snapshot

---

## System-Wide Impact

- **Interaction graph:** No production code changes. CMake test infrastructure gains one new library target and one regeneration target. Existing test executables gain an additional link dependency.
- **Error propagation:** Invariant fixture failures surface as pFUnit assertion failures with descriptive messages (which invariant, which row). Krylov comparison failures report iteration number and divergence magnitude.
- **State lifecycle risks:** None — test-only, no persistent state.
- **API surface parity:** All CSR operations get consistent test coverage. No production API changes.
- **Integration coverage:** Krylov snapshots test the full config -> initialization -> Hamiltonian assembly -> SpMV pipeline end-to-end, catching integration bugs that unit-level invariant tests cannot.
- **Unchanged invariants:** Production CSR operations, Hamiltonian assembly routines, and all existing test behaviors remain identical. Only additive test code and CMake changes.

---

## Risks & Dependencies

| Risk | Mitigation |
|------|------------|
| Reference vectors become stale when Hamiltonian construction intentionally changes | Documented one-command regeneration (`regenerate_krylov_references` target); reference data module has header comment explaining regeneration |
| pFUnit assertion patterns insufficient for structural invariant checks | Research confirmed pFUnit supports `@assertTrue` with message — sufficient for all invariant checks |
| Small-grid Krylov snapshots insufficiently sensitive to structural bugs | k=5 chain propagates errors across connected components; edge-case test validates detection of deliberate corruption |
| Refactoring existing tests introduces regressions | U3 is purely structural (move helpers); all tests must pass identically before and after |
| `kron_dense_dense_1d` dead code tested but never used | Flagged in code comment; removal deferred to follow-up |

---

## Sources & References

- **Origin document:** [docs/brainstorms/2026-05-09-csr-structure-testing-requirements.md](docs/brainstorms/2026-05-09-csr-structure-testing-requirements.md)
- Related code: `src/math/sparse_matrices.f90` (CSR type and operations)
- Related code: `src/physics/hamiltonian_wire.f90` (wire Hamiltonian assembly)
- Related code: `src/physics/bdg_hamiltonian.f90` (BdG assembly)
- Related code: `src/physics/topological_analysis.f90` (QWZ, BHZ assembly)
- Related tests: `tests/unit/test_csr_spmv.pf`, `tests/unit/test_hamiltonian_2d.pf`
- Documented bugs: `docs/solutions/logic-errors/topological-magnetic-index-logic-errors-2026-05-08.md`
- Verification ladder: `docs/solutions/best-practices/8band-verification-ladder-2026-05-09.md`
