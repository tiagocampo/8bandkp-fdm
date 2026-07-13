**Status**: COMPLETE (2026-07-05)

# PRD: Architecture Cleanup — Dead Code Removal + Landau Gap in Simulation Setup (C6 + C1)

## Problem Statement

The `simulation_setup` module handles 3 of 4 confinement modes (bulk, QW, wire) but bypasses Landau, forcing `main.f90` to duplicate ~43 lines of initialization and workspace query inline. Meanwhile, the codebase accumulates dead code (`add_zeeman_coo` with zero callers, `setup_alloc_sweep` with zero callers, stale `use confinement_init` import), misleading error messages (`'FEAST convergence failed at k-point 1'` in `zheevx` LAPACK paths), and deprecated `stop` usage that should be `error stop` per project convention.

These issues actively harm debugging — a developer seeing "FEAST convergence failed" in a LAPACK code path wastes time investigating the wrong eigensolver. Dead public exports pollute module interfaces and force unnecessary recompilation.

## Solution

Two coordinated changes in a single PR with distinct commits:

1. **C6 — Dead code removal and error message cleanup:** Mechanical deletion of unused code, correction of misleading error strings, and replacement of deprecated `stop` with `error stop`. Zero behavioral change. Seven items across five files.

2. **C1 — Close Landau gap in simulation_setup:** Add `case('landau')` to `simulation_setup_init` and `setup_build_H`, absorbing the Landau init block and workspace query from `main.f90`. Restructure `main.f90` so Landau becomes a self-contained block with early `stop` (mirroring the existing wire pattern). The k-sweep (uses `zheevx` + OpenMP + thread-local arrays) and B-sweep (Landau-specific fan diagram) remain inline in `main.f90`.

The B-sweep exclusion is documented in ADR 0003.

## User Stories

### Dead Code and Error Messages (C6)

1. As a developer debugging a diagonalization failure in the QW k-sweep, I want the error message to say "zheevx diagonalization failed" instead of "FEAST convergence failed at k-point 1", so that I investigate the correct eigensolver path
2. As a developer debugging a diagonalization failure in the Landau k-sweep, I want the error message to say "zheevx diagonalization failed" instead of "FEAST convergence failed at k-point 1", so that I investigate the correct eigensolver path
3. As a developer debugging a diagonalization failure in the Landau B-sweep, I want the error message to say "zheevx diagonalization failed in B-sweep" instead of "FEAST convergence failed at k-point 1", so that I know which sweep failed
4. As a developer debugging a wire FEAST failure at k-point 50, I want the error stop to say "FEAST convergence failed" without hardcoding "k-point 1", so that the printed detail above (which correctly shows k) is the authoritative source
5. As a developer reading `magnetic_field.f90`, I want dead code (`add_zeeman_coo`, zero callers, `stop 1` instead of `error stop`) removed, so that the module contains only live code
6. As a developer reading `simulation_setup.f90`, I want dead code (`setup_alloc_sweep`, zero callers) removed, so that the public interface reflects what is actually used
7. As a developer reading `hamiltonianConstructor.f90`, I want the stale `use confinement_init` import removed, so that module dependencies are accurate
8. As a developer hitting a fatal error in `main.f90` bulk diagonalization, I want `stop "error diag"` to be `error stop "error diag"`, so that the exit follows project convention
9. As a developer hitting a fatal error in `main_gfactor.f90`, I want `stop 'evnum...'` to be `error stop 'evnum...'`, so that the exit follows project convention

### Landau Gap in Simulation Setup (C1)

10. As a developer adding a new confinement-mode feature, I want `simulation_setup` to handle all 4 modes (bulk, QW, wire, Landau), so that I don't need to add parallel init paths in `main.f90`
11. As a developer modifying the Landau initialization logic, I want the init code in `simulation_setup_init` rather than inline in `main.f90`, so that the change is localized to one place
12. As a developer reading `main.f90`, I want Landau structured as a self-contained block with early stop (like wire), so that the flow is clear: wire → stop, Landau → stop, bulk/QW → shared path
13. As a developer working on the bulk/QW shared setup path, I want Landau branches (`else if (conf_direction == 'x')`) removed from that path, so that the shared code only handles the modes that actually share it
14. As a developer building a Landau Hamiltonian, I want `setup_build_H` to dispatch to `ZB8bandLandau` via the same interface as bulk/QW, so that Hamiltonian construction is uniform across modes
15. As a future developer wondering why the Landau B-sweep isn't in `simulation_setup`, I want ADR 0003 to explain that the B-sweep is specialized physics output (fan diagram, B-field parameterization) that belongs in the application layer
16. As a developer running `bandStructure` with `confinement='landau'`, I want the program to produce identical output before and after this refactor, so that I can trust the change is structural only

## Implementation Decisions

### C6: Dead Code and Error Messages

- **Delete `add_zeeman_coo` from `magnetic_field.f90`:** The subroutine has zero callers. Its `stop 1` (not `error stop`) is moot since the whole subroutine is deleted. Remove from public exports.
- **Delete `setup_alloc_sweep` from `simulation_setup.f90`:** The subroutine has zero callers. Remove from public exports.
- **Delete `use confinement_init` from `hamiltonianConstructor.f90`:** The import is only referenced in a comment on line 295. Remove the import; the comment remains.
- **Fix error messages in `main.f90` zheevx paths:** Lines 645, 709, 786 all say `'FEAST convergence failed at k-point 1'` but are in `zheevx` (LAPACK) code paths. Replace with mode-specific messages: `'zheevx diagonalization failed'` for QW and Landau k-sweeps, `'zheevx diagonalization failed in B-sweep'` for the Landau B-sweep. The detailed context (k value, info code) is already printed above each `error stop`.
- **Fix hardcoded "k-point 1" on line 327:** Inside the wire k-sweep loop, `error stop 'FEAST convergence failed at k-point 1'` fires for all k values. Change to `'FEAST convergence failed'` — the k value is already printed on line 325. Fortran `error stop` accepts only character literals, not interpolated strings, so the detail must live in the `print` above.
- **`stop` → `error stop` in error exits:** Line 586 (`stop "error diag"`) and `main_gfactor.f90` line 272 (`stop 'evnum...'`) are fatal error exits that should use `error stop` per project convention.
- **Preserve `stop` on line 380:** This is a normal successful exit ("wire mode complete"), not an error. The early-return pattern for `program` units requires `stop`. Leave unchanged.

### C1: Landau Gap in Simulation Setup

- **Add `case('landau')` to `simulation_setup_init`:** Calls `confinementInitialization_landau` with single-material params, allocates and stores `profile` and `kpterms`, sets `N = landau%nx * 8`, computes `il`/`iuu` from band selection, performs `zheevx` workspace query. No strain, no electric field, no self-consistent loop — Landau is single-material only. Reuses existing `setup%profile`, `setup%kpterms`, `setup%HT`, `setup%work` fields (no new type components).
- **Add `case('landau')` to `setup_build_H`:** Calls `ZB8bandLandau(HT_out, kvec, setup%profile, setup%kpterms, cfg%grid%x, cfg=cfg)` via the existing `HT_out` dense-matrix path, same as QW.
- **No `case('landau')` in `setup_solve_kpoint_serial`:** The Landau k-sweep uses `zheevx` (range-based) with OpenMP parallel thread-local workspace — a different eigensolver pattern than the full-spectrum `zheev`/`zheevd` used by `setup_solve_kpoint_serial` for QW. The k-sweep stays inline in `main.f90`.
- **Restructure `main.f90` into self-contained mode blocks:** After wire block and its `stop`, add a Landau block that calls `simulation_setup_init`, then performs the k-sweep and B-sweep inline, then `stop`. The shared bulk/QW setup path (matrix dims, workspace, SC, profile output) loses its `else if (conf_direction == 'x')` branches and becomes bulk/QW-only.
- **B-sweep stays in `main.f90` (ADR 0003):** The B-sweep is 81 lines of specialized Landau physics output (magnetic field parameterization, fan diagram). Moving it into `simulation_setup` would relocate complexity without reducing it. The method would be Landau-specific, violating the generic intent of the module.
- **Wire topology subroutines deferred to C5:** The three topology subroutines in `main_topology.f90` that bypass `simulation_setup` by calling `confinementInitialization_2d` directly are a separate concern. They need a lightweight wire-only init path, which is C5's `wire_setup` type extraction.
- **green_functions re-init deferred to C5:** The spectral function module's direct call to `confinementInitialization_2d` is the same wire duplication pattern. Also C5.

### Documentation Updates

- **Update CLAUDE.md line 176:** `simulation_setup.f90` description to note it handles all 4 confinement modes.
- **ADR 0003:** Documents the B-sweep exclusion decision. Written at `docs/adr/0003-landau-b-sweep-stays-in-main.md`.
- **No changes needed to `src/physics/AGENTS.md`** or other AGENTS files — the pattern description for wiring through `simulation_setup` is already correct.

### Commit Structure

Three commits on a single branch:

1. `refactor: dead code removal and error message cleanup` (C6)
2. `refactor: close Landau gap in simulation_setup` (C1)
3. `docs: update CLAUDE.md and add ADR 0003`

### Branch Strategy

Branch `refactor/architecture-cleanup` forked from `main` after `feat/publishable-benchmarks-phase21` merges.

## Testing Decisions

### What makes a good test

- **C6:** No new tests needed. Changes are deletions and string replacements with zero behavioral impact. The existing test suite (113 tests) validates that nothing breaks. The error message changes only fire on LAPACK failures, which don't occur in normal test runs.
- **C1:** Behavioral equivalence is the primary test. Running `bandStructure` with a Landau config should produce identical output files before and after the refactor. The existing Landau regression tests (if any) and verification ladder tests cover this.

### Modules tested

- **C1:** The `simulation_setup_init` and `setup_build_H` Landau branches should be exercised by existing Landau integration tests. If no dedicated Landau regression test exists, one should be added to prevent future regressions of the setup pipeline.

### Prior art

- `tests/regression/` golden-output tests validate end-to-end output stability
- `tests/integration/` verification ladder tests validate physics correctness
- ADR 0001 (`simulation_setup` fat type) established the pattern of routing modes through `simulation_setup_init` — C1 extends this to Landau

## Out of Scope

- **Wire topology subroutines routing through simulation_setup** — deferred to C5 (wire_setup type extraction)
- **green_functions re-initialization** — deferred to C5
- **Zeeman single-source-of-truth violation** (C3) — independent PR
- **Namespace cleanup / broad imports** (C4) — independent PR
- **Strain_solver decomposition** (C8) — independent PR
- **Pipeline stage assertions** (C9) — independent PR
- **BdG COO workspace reuse** (C10) — independent PR
- **Confinement init unit tests** (C7) — independent PR
- **Eigensolve/output absorption** (C2) — deferred until C1 lands
- **Adding `setup_solve_kpoint_serial` for Landau** — not needed (k-sweep uses different eigensolver pattern)
- **Any behavioral changes to physics output** — explicitly zero

## Further Notes

- This PR is the first of a 5-PR architecture cleanup campaign. The recommended campaign order is: C6+C1 (this PR) → C3+C8 (Zeeman SSOT + strain split) → C4 (namespace) → C5 (wire setup extraction, depends on C1) → C2+C7+C9+C10 (eigensolve absorption, confinement tests, pipeline assertions, BdG workspace).
- C5 depends on C1 being merged first, because C5's `wire_setup` type can live inside `simulation_setup` once Landau is also a case in the setup module.
- The architecture review (2026-06-04 deep pass) produced 47 raw findings → 36 surviving (15 Strong, 11 Worth exploring, 10 Speculative), grouped into 10 candidates. This PR addresses the top 2 Strong candidates.
