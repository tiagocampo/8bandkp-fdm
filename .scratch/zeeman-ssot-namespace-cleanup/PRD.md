# PRD: Zeeman Single-Source-of-Truth Fix + Namespace Pollution Cleanup (C3 + C4)

## Problem Statement

The Zeeman diagonal coefficients are defined in **two independent places**: the canonical `get_zeeman_table()` in `strain_solver.f90` (declared as single source of truth in CLAUDE.md) and a hard-coded copy in `compute_zeeman_vz()` in `magnetic_field.f90`. The BdG module uses the `magnetic_field` copy, while the dense and COO Hamiltonian builders use the table. These must be kept in sync manually — a real SSOT violation. Additionally, the Bir-Pikus `field_id` dispatch (7-way `select case`) is character-for-character duplicated between `hamiltonianConstructor.f90` and `hamiltonian_wire.f90`, creating a maintenance hazard.

Separately, 22 modules import `definitions` without `only:` clauses — any new public symbol forces recompilation of all 22. `finitedifferences.f90` exports 16 symbols but 10 are never called externally. `get_unit` and `ensure_output_dir` live in `outputFunctions.f90`, causing 4 physics modules and the input parser to depend on the I/O subsystem for basic file utilities.

## Solution

Two coordinated changes in a single PR:

1. **C3 — Zeeman SSOT fix:** Move the Zeeman table (`zeeman_entry` type, `build_zeeman_table`, `get_zeeman_table`, cache, `init_zeeman_cache`) from `strain_solver.f90` to `magnetic_field.f90` — its domain-correct home. Rewrite `compute_zeeman_vz` to read from the table instead of hard-coding. Extract a shared `lookup_bp_field` pure function into `strain_solver.f90` to eliminate the duplicated `field_id` dispatch.

2. **C4 — Namespace pollution cleanup:** Add `only:` clauses to all 22 broad `use definitions` imports. Make 10 unused FD symbols PRIVATE. Move `get_unit` and `ensure_output_dir` from `outputFunctions.f90` to `utils.f90`.

All changes are structural — zero behavioral change.

## User Stories

### Zeeman SSOT (C3)

1. As a developer modifying the Zeeman diagonal coefficients (e.g. adding spin-orbit corrections), I want them defined in exactly one place, so that I don't have to remember to update two independent copies
2. As a developer reading `magnetic_field.f90`, I want `compute_zeeman_vz` to use the same table as the Hamiltonian builders, so that I can trust the BdG Zeeman splitting matches the dense/COO paths
3. As a developer reading the BdG Hamiltonian code, I want the Zeeman coefficients to come from the same source as the rest of the codebase, so that inconsistencies are impossible
4. As a developer adding a new Bir-Pikus strain field (e.g. field_id 8), I want the `field_id` → `bir_pikus_blocks` mapping defined in one function, so that I update one place instead of two identical `select case` blocks
5. As a developer reading the CLAUDE.md invariants, I want the Zeeman SSOT location to accurately reflect where the table lives in code, so that the documentation is trustworthy
6. As a developer running BdG simulations, I want identical physics output before and after this refactor, so that I can trust the change is structural only

### Namespace Pollution (C4)

7. As a developer adding a new public symbol to `definitions`, I want only the modules that actually use it to recompile, so that my build stays fast
8. As a developer reading a module's imports, I want to see exactly which symbols from `definitions` it uses (via `only:`), so that I can understand the module's real dependencies at a glance
9. As a developer reading `finitedifferences.f90`, I want the public interface to contain only the 6 symbols that external callers use, so that the module's API surface reflects its real contract
10. As a developer working on `input_parser.f90`, I want it to import `get_unit` and `ensure_output_dir` from a utility module rather than the I/O subsystem, so that the parser doesn't depend on output formatting code
11. As a developer working on a physics module that needs `get_unit`, I want to import it from `utils` rather than `outputFunctions`, so that my module's dependency graph reflects its actual needs (file utility, not formatted output)
12. As a developer adding a new executable, I want the file utilities (`get_unit`, `ensure_output_dir`) available without pulling in the output formatting module, so that my new executable has minimal dependencies

## Implementation Decisions

### C3: Zeeman Table Relocation

- **Move Zeeman table to `magnetic_field.f90`:** The `zeeman_entry` type, `build_zeeman_table` function, `get_zeeman_table` function, `zeeman_table_cache`/`zeeman_table_cached` SAVE variables, and `init_zeeman_cache` subroutine all move from `strain_solver.f90` to `magnetic_field.f90`. This is the domain-correct home — Zeeman is a magnetic effect, not a strain effect. Pre-empts the Zeeman portion of C8 (strain split).
- **Rewrite `compute_zeeman_vz` to use the table:** Instead of hard-coding 8 coefficients, loop over `get_zeeman_table()` results: `Vz(b) = table(b)%g_multiplier * E0`. Loses `pure` attribute (table uses SAVE cache), but the only caller (`bdg_hamiltonian.f90`) is non-pure — no impact.
- **Update import sites:** `hamiltonianConstructor.f90` and `hamiltonian_wire.f90` move `zeeman_entry, get_zeeman_table` from their `use strain_solver` to a new `use magnetic_field` import. `main.f90` moves `init_zeeman_cache` from `use strain_solver` to `use magnetic_field`. No circular dependency: `strain_solver` does not import `magnetic_field`.
- **No change to `init_strain_cache`:** The strain table cache initializer stays in `strain_solver.f90` — it's strain-related.
- **Extract `lookup_bp_field` into `strain_solver.f90`:** A new `pure function lookup_bp_field(bp, field_id, ii) result(field_val)` encapsulates the 7-way `select case (field_id)` mapping from `bir_pikus_blocks` components. Both `apply_strain_table_dense` (hamiltonianConstructor) and `insert_strain_coo` (hamiltonian_wire) call this instead of inlining the switch. Natural home: the `field_id` values are defined by the strain table in `strain_solver.f90`.
- **`strain_solver.f90` retains:** Bir-Pikus formulas (`compute_bp_scalar`, `compute_bir_pikus_blocks`), strain table (`strain_entry`, `get_strain_table`, `build_strain_table`), Navier-Cauchy solver (`compute_strain_wire`), and the new `lookup_bp_field` helper. Loss: ~40 lines of Zeeman code.

### C4: Namespace Cleanup

- **Add `only:` to all 22 `use definitions` imports:** The heaviest user (`gfactor_functions.f90`) references only 16 of 88 public symbols — no pragmatic cutoff needed. The median is ~6 symbols. All 22 modules get explicit `only:` lists.
- **Make 10 FD symbols PRIVATE:** `FDmatrixDense`, `FDstencil`, `toeplitz`, `FDcentralCoeffs2nd`, `FDcentralCoeffs1st`, `FDforwardCoeffs2nd`, `FDbackwardCoeffs2nd`, `FDforwardCoeffs1st`, `FDbackwardCoeffs1st`, `vandermonde_2nd_deriv`, `vandermonde_1st_deriv`, `vandermonde_interp`. The 6 that stay public: `Identity`, `buildFD2ndDerivMatrix`, `buildFD1stDerivMatrix`, `buildStaggeredD1Inner`, `buildStaggeredD1Outer`, `interpolateToHalfPoints`.
- **Move `get_unit` and `ensure_output_dir` to `utils.f90`:** Both are file-level utilities with no dependency on output formatting logic. `outputFunctions.f90` imports them from `utils` after the move. Four physics modules (`optical_spectra`, `scattering`, `exciton`, `main_optics`) already use `only:` — they change their import target. Four broad importers (`main.f90`, `main_gfactor.f90`, `main_topology.f90`, `input_parser.f90`) also update.
- **Ordering:** Do the Zeeman move and FD visibility changes first, then the `only:` sweep last — because the `only:` sweep touches all 22 files and should include the updated imports from the Zeeman relocation.

### Documentation Updates

- **Update CLAUDE.md Zeeman SSOT invariant:** Change "Single source of truth: [...] Zeeman table (`strain_solver.f90`)" to reference `magnetic_field.f90`. Update the "NEVER change the strain or Zeeman block tables" boundary to reference `magnetic_field.f90` for Zeeman.
- **Update CLAUDE.md dependency graph:** The `magnetic_field` module entry gains Zeeman table exports. The `strain_solver` entry loses Zeeman.
- **Update `src/physics/AGENTS.md`:** Adjust module descriptions for `magnetic_field.f90` (now exports Zeeman table) and `strain_solver.f90` (no longer exports Zeeman). Add `lookup_bp_field` to strain_solver's public API.

### Commit Structure

Three commits on a single branch:

1. `refactor: move Zeeman table to magnetic_field and extract lookup_bp_field` (C3)
2. `refactor: namespace cleanup — only imports, private FD exports, file utils to utils` (C4)
3. `docs: update CLAUDE.md and AGENTS.md for Zeeman SSOT and namespace changes`

### Branch Strategy

Branch `refactor/zeeman-ssot-namespace-cleanup` forked from `main` (after PR #36 merge).

## Testing Decisions

### What makes a good test

- **C3:** No new tests needed. The Zeeman table values are tested transitively through the Hamiltonian unit tests (`test_hamiltonian.pf`, `test_hamiltonian_blocks.pf`) and the BdG unit test (`test_bdg.pf`). The `lookup_bp_field` helper is a pure function with identical behavior to the inlined code it replaces. Behavioral equivalence is the primary test — all 113 existing tests must pass.
- **C4:** No new tests needed. All changes are mechanical (import adjustments, visibility changes, subroutine relocation). The existing test suite validates compilation correctness and behavioral equivalence.

### Modules tested

- Both C3 and C4 rely on the existing 113-test suite (35 unit + 44 regression + 18 verification + 9 standard-star + 6 convergence + 2 strain + 1 coverage + misc).

### Prior art

- PR #36 (C6+C1) established the pattern of structural refactoring verified by the full test suite
- `test_hamiltonian.pf` tests Zeeman diagonal placement via the Hamiltonian builders
- `test_magnetic_field.pf` was deleted in C6 (tested the now-removed `add_zeeman_coo`)

## Out of Scope

- **Wire topology subroutines routing through simulation_setup** — deferred to C5 (wire_setup type extraction)
- **green_functions re-initialization** — deferred to C5
- **Eigensolve/output absorption** (C2) — deferred until C5 lands
- **Strain_solver decomposition** (C8 — splitting into bir_pikus.f90, navier_cauchy.f90) — this PR pre-empts the Zeeman portion of C8 but does not split the remaining 3 concerns
- **Pipeline stage assertions** (C9) — independent PR
- **BdG COO workspace reuse** (C10) — independent PR
- **Confinement init unit tests** (C7) — independent PR
- **Any behavioral changes to physics output** — explicitly zero

## Further Notes

- This is the second PR in the architecture cleanup campaign. PR #36 (C6+C1) addressed dead code and the Landau simulation_setup gap.
- The Zeeman table move to `magnetic_field.f90` pre-empts the Zeeman portion of C8, meaning when C8 is eventually done, the Zeeman work is already complete. C8 would then only need to split `strain_solver.f90` into `bir_pikus.f90` and `navier_cauchy.f90`.
- The `compute_zeeman_vz` function loses its `pure` attribute when reading from the SAVE-cached table. This is acceptable because: (a) the only caller (`bdg_hamiltonian.f90`) is non-pure, (b) the cache is pre-warmed by `init_zeeman_cache()` before any OpenMP fork, so the SAVE variable is never actually written in a parallel context.
