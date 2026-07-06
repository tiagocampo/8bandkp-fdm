**Status**: COMPLETE (2026-07-05)

Status: ready-for-agent

# PRD: Unify Hamiltonian Block Structure (C4) + Wire Convergence Fix

## Problem Statement

The 8-band zinc-blende k.p Hamiltonian block structure — which bands couple, with what prefactors, and how k.p terms map to matrix elements — is encoded independently in two modules: `hamiltonianConstructor.f90` (dense, 52 HT() assignments for k.p + 38 inline strain + 8 Zeeman) and `hamiltonian_wire.f90` (sparse/COO, 52 insert_csr_block* calls for k.p + 32-entry table-driven strain + 8 Zeeman). The strain table already exists in `strain_solver.f90` but only the wire path uses it; the dense path has two inline copies of the same Bir-Pikus entries. A bug fix in one representation must be manually replicated in the other, and the two copies have already diverged in structure. Additionally, the wire convergence test (U6) has a path-handling bug that prevents standalone execution.

## Solution

Three phases of unification, plus a wire convergence test fix:

**Wire convergence fix (U6):** Make the wire convergence test script robust to both absolute and relative paths for the build directory, matching the pattern used by other convergence tests. Verify wire convergence results are produced.

**Phase A (strain table for dense path):** Replace the two inline strain implementations in `hamiltonianConstructor.f90` with a new `apply_strain_table_dense` subroutine that reads `get_strain_table()` from `strain_solver` and applies the 32 entries to the dense matrix. This eliminates the most dangerous duplication: Bir-Pikus sign conventions are easy to get wrong, and CLAUDE.md's boundaries section warns against changing them.

**Phase B (Zeeman table):** Extract the 8-entry Zeeman diagonal into a data-driven table in `strain_solver` (alongside the strain table), and make both the dense and COO paths use it. Follows the identical pattern as Phase A.

**Phase C (k.p block table):** Create a new `hamiltonian_blocks.f90` module that defines the full 8-band k.p block structure as data — a table of 52 entries mapping (band_pair, kp_term, prefactor). Both the dense builder (`hamiltonianConstructor`) and the COO builder (`hamiltonian_wire`) read this table and emit into their representation. Derivative modes (dH/dk for g-factor velocity operators) are handled by the builder, not the table.

## User Stories

### Wire Convergence (U6)

1. As a test developer, I want the wire convergence test to work when invoked directly with a relative build directory, so that I can run it standalone for debugging without going through ctest.
2. As a test developer, I want the wire convergence test to produce non-empty JSON results for both subband energy and g-factor, so that Richardson extrapolation validates wire physics.

### Phase A: Strain Table Consolidation

3. As a physics researcher, I want the dense Hamiltonian builder to use the same strain table as the wire builder, so that Bir-Pikus sign conventions are verified once against Chuang/Winkler, not twice.
4. As a developer, I want the inline strain code in `ZB8bandQW` (the copy made to avoid compiler temporaries) to be replaced by the table-driven approach, so that there is exactly one representation of strain in the codebase.
5. As a developer, I want the standalone `add_bp_strain_dense` subroutine to be replaced by the table-driven approach, so that the two strain codepaths (inline + standalone) collapse into one.
6. As a developer, I want the table-driven dense strain to avoid array temporaries, so that the original -O3 + OpenMP stack corruption issue does not recur.
7. As a developer, I want the scalar bulk strain (`apply_bp_strain_inline`) to also adopt the table, so that all three strain variants (QW inline, QW standalone, bulk scalar) are unified.
8. As a codebase maintainer, I want the strain table in `strain_solver` to be the single source of truth for Bir-Pikus block entries, so that any future change (e.g., adding crystal-field splitting) touches one place.

### Phase B: Zeeman Table

9. As a physics researcher, I want the Zeeman splitting to be defined in a data-driven table, so that the 8 band-dependent g-factors (HH=-1.5, LH=+0.5, SO=-0.5, CB=+1.0) are tabulated once.
10. As a developer, I want the dense Zeeman insertion (8 inline HT() assignments in `ZB8bandQW`) to use the same table as the COO Zeeman insertion (`insert_zeeman_coo`), so that Zeeman entries are not duplicated.
11. As a developer, I want the bulk Zeeman (in `ZB8bandBulk`) to also use the table, so that all three geometries (bulk, QW, wire) share the Zeeman definition.

### Phase C: k.p Block Table

12. As a physics researcher, I want the 8-band zinc-blende k.p block structure defined as data in a dedicated module, so that the coupling pattern (52 band-pair entries with k.p terms Q, R, S, T, P) is verified once against Winkler/Chuang.
13. As a developer, I want a `hamiltonian_blocks.f90` module that exports the k.p block table, so that both the dense and COO builders import from a single source.
14. As a developer, I want each k.p entry to identify the k.p term by a named constant (KP_Q, KP_R, KP_S, KP_T, KP_P, KP_DIAGONAL), so that the table says "which term, which band pair, with what prefactor" and the builder translates that to the appropriate matrix element.
15. As a developer, I want the dense builder to read the k.p block table and emit HT() assignments, so that the 52 hardcoded assignments in `ZB8bandQW` are replaced by a data-driven loop.
16. As a developer, I want the COO builder to read the k.p block table and emit CSR block inserts, so that the 52 `insert_csr_block*` calls in `insert_main_blocks` are replaced by a data-driven loop.
17. As a developer, I want the builder (not the table) to handle derivative modes (dH/dk for g-factor velocity operators), so that the table encodes only the standard physics.
18. As a developer, I want adding a new term (e.g., Dresselhaus spin-orbit coupling) to require adding one table entry, not modifying two modules, so that the codebase is maintainable as the physics expands.
19. As a developer, I want band offsets (the 8 diagonal entries from material parameters) to be included in the k.p block table, so that the diagonal structure is also data-driven.
20. As a codebase maintainer, I want `hamiltonianConstructor.f90` to shrink by replacing 52 + 38 + 8 = 98 hardcoded assignments with table-driven loops, so that the module stays under the 300-line guideline.
21. As a codebase maintainer, I want `hamiltonian_wire.f90` to import the k.p block table for its block insertion, so that the 52 `insert_csr_block*` calls in `insert_main_blocks` are replaced.
22. As a developer, I want the existing convergence tests (U4, U5, U6) to validate that the block table refactoring does not change physics results, so that regression is caught automatically.

## Implementation Decisions

### Module structure: new `hamiltonian_blocks.f90`

A new module in `src/physics/` that defines:
- Named constants for k.p terms: `KP_Q`, `KP_R`, `KP_S`, `KP_T`, `KP_P`, `KP_PZ`, `KP_PM`, `KP_PP`, `KP_A`, `KP_DIAGONAL`
- A `kp_entry` derived type with fields: `row_band`, `col_band` (0-indexed), `kp_term` (named constant), `prefactor` (real or complex), `use_conjg` (logical)
- A `get_kp_block_table()` function that returns the 52-entry table (cached, same pattern as `get_strain_table`)
- Band offset entries included in the table as `KP_DIAGONAL` entries

The module depends only on `defs.f90` (for `dp` kind). No physics module dependencies.

### Phase A: Dense strain uses the strain table

The `hamiltonianConstructor` module currently imports `compute_bp_scalar` and `bir_pikus_blocks_free` from `strain_solver`. After Phase A, it additionally imports `strain_entry`, `get_strain_table` (same imports as `hamiltonian_wire`).

A new subroutine `apply_strain_table_dense(HT, ii, N, bp)` reads `get_strain_table()` and applies all 32 entries to the dense matrix at spatial index `ii`. This replaces:
- The inline strain block in `ZB8bandQW` (lines 224-284, 61 lines of inline code)
- The standalone `add_bp_strain_dense` (lines 864-919, 56 lines)
- The scalar `apply_bp_strain_inline` (lines 922-970, 49 lines)

The original inline was done to avoid compiler temporaries with `-O3 + OpenMP`. The table-driven approach uses scalar field lookups from the `bir_pikus_blocks` type (no array temporaries), so the original reason for inlining should not apply. This must be verified: run `OMP_NUM_THREADS=4 ctest -L unit` and a regression test to confirm no stack corruption.

### Phase B: Zeeman table in `strain_solver`

Add a `zeeman_entry` type and `get_zeeman_table()` to `strain_solver`, following the same pattern as the strain table. The table has 8 entries (one per band) with fields: `band_index`, `g_multiplier` (HH=-1.5, LH=+0.5, SO=-0.5, CB=+1.0). Both dense and COO paths call `compute_zeeman_vz` internally, which already encapsulates the band-dependent g-factors — the table formalizes this data and makes the 8 values visible as data.

### Phase C: k.p block table in `hamiltonian_blocks.f90`

The 52-entry k.p block table captures the full 8x8 block structure of the zinc-blende Hamiltonian. Each entry maps a band pair to a k.p term with a prefactor. The builders (dense adapter, COO adapter) read the table and emit into their representation.

**Dense adapter:** Iterates over entries, computes the k.p term value from precomputed `kpterms` arrays (Q, R, S, T, P from `confinementInitialization`), applies prefactor and conjugation, and writes to `HT(row_band*N + ii, col_band*N + ii)`.

**COO adapter:** Iterates over entries, gets the CSR block for each k.p term from the precomputed blocks, and inserts via `insert_csr_block*`.

**Derivative modes:** The builder has a mode flag. For derivative mode (dH/dk), the builder applies the power rule to the k.p term (e.g., k^2 becomes 2k) before emitting. Only the subset of entries with non-zero k-dependence contributes. The existing `insert_g3_blocks` (20 entries for dH/dkz) is handled by the derivative-mode builder.

**Band offsets:** Included in the table as `KP_DIAGONAL` entries. The dense builder applies them to the diagonal blocks; the COO builder inserts them via the existing profile diagonal mechanism.

### Phasing and commit strategy

Three commits (after the wire convergence fix):
1. **Phase A+B combined**: Strain table + Zeeman table for dense path. One commit because they follow the same pattern and should be validated together.
2. **Phase C**: k.p block table + `hamiltonian_blocks.f90` module. One commit covering the new module + both adapters.

### Validation strategy

The existing Richardson convergence tests (U4: QW grid, U5: QW order, U6: wire) provide physics regression validation. After each phase, run `ctest -L convergence` and verify Richardson-extrapolated values are unchanged within machine precision. The 33 unit tests (including `test_hamiltonian` and `test_hamiltonian_2d`) must continue to pass.

## Testing Decisions

### What makes a good test

Tests verify external behavior (physics output) not implementation details (which subroutines are called). A good test for the block table refactoring produces identical eigenvalues, g-factors, and absorption spectra before and after the change.

### Modules to test

- **`hamiltonian_blocks.f90`** (Phase C): Unit tests for `get_kp_block_table()` — verify entry count (52), verify all 52 band pairs are covered, verify no duplicate (row, col) pairs, verify Hermitian symmetry (entry for (i,j) has conjugate entry for (j,i)). New pFUnit test file `tests/unit/test_hamiltonian_blocks.pf`.
- **`apply_strain_table_dense`** (Phase A): Tested implicitly via existing `test_hamiltonian` pFUnit test. Additionally, a targeted unit test comparing the output of `apply_strain_table_dense` against the old `add_bp_strain_dense` for identical input, verifying element-wise agreement to machine precision.
- **Phase A regression**: Run full convergence test suite before and after. Richardson-extrapolated values must agree within 1e-12 eV.
- **Wire convergence test (U6)**: Fix the path handling so `python3 test_wire_convergence.py build-nolto .` works. Verify JSON output contains CB1_energy and gz results.

### Prior art

- `tests/unit/test_hamiltonian.pf` — tests dense Hamiltonian construction including strain
- `tests/unit/test_hamiltonian_2d.pf` — tests wire Hamiltonian construction
- `tests/unit/test_fd_convergence.pf` — operator-level convergence test
- `tests/integration/convergence_helpers.py` — Richardson extrapolation infrastructure
- `.scratch/optics-engine-encapsulation/PRD.md` — PRD format precedent

## Out of Scope

- C5 (split `simulation_config` into domain configs) — speculative, high blast radius, deferred
- Magnetic field / Landau level terms in the block table — wire does not implement Landau levels; defer until needed
- Foreman renormalization terms in the block table — disabled by default, not a duplication source
- S5 (InAs/GaSb broken-gap) negative convergence rate investigation — flagged as "needs calibration"
- SC CB1_shift 40% GCI failure investigation — flagged as "needs calibration"
- Dresselhaus spin-orbit coupling — mentioned as future use case for the table, not implemented now
- Renaming or restructuring `hamiltonianConstructor.f90` beyond what Phase A/B/C requires
- Any changes to the eigensolver, Poisson solver, or SC loop modules

## Further Notes

### Known convergence test issues

- S5 (InAs/GaSb broken-gap QW) shows negative convergence rates for CB1_energy (-4.81) and subband_spacing (-1.96). The Richardson plan explicitly deferred S5 calibration. The convergence data is informational; rate assertions will be tuned after empirical calibration.
- SC CB1_shift has 40.75% GCI and fails the convergence assertion. The grid spacing range may be too coarse for SC convergence. Flagged for future investigation.
- Wire convergence (U6) produces valid results when called with absolute paths (CB1 Richardson -0.8445 eV, GCI 0.64%, rate 4.57; g-factor Richardson 9.64). The test has a relative-path bug preventing standalone execution.

### Entry counts summary

| Category | Dense (inline) | COO (table-driven) | Unified (target) |
|----------|---------------|-------------------|-----------------|
| k.p blocks | 52 HT() assignments | 52 insert_csr_block* calls | 52-entry kp_entry table |
| Strain | 38 HT() x 2 copies + scalar | 32-entry table | 32-entry strain_entry table (shared) |
| Zeeman | 8 HT() x 3 copies | 8 COO entries | 8-entry zeeman_entry table (shared) |

### Branch

All work on `feat/richardson-observables-convergence`. The Richardson convergence tests and C4 refactoring are complementary: convergence tests validate that C4 doesn't change physics.
