**Status**: COMPLETE (2026-07-05)

Status: ready-for-agent

# PRD: Unified Simulation Setup Orchestration

## Problem Statement

The three application programs (`bandStructure`, `gfactorCalculation`, `opticalProperties`) each independently implement the same physics initialization pipeline: confinement initialization, strain computation, eigensolver configuration, and self-consistent Schrodinger-Poisson loop integration. This duplication has two consequences. First, certain physics combinations are architecturally impossible -- for example, computing optical properties of a strained and doped quantum well requires both the SC loop and strain setup, but the optics app has neither. Second, bugs fixed in one app's setup must be manually replicated in the others, and the copies have already begun to diverge (wire SC uses different DIIS warm-up, different Fermi search).

## Solution

Extract a `simulation_setup` module that encapsulates the full physics initialization pipeline behind a single entry point. Each application calls `simulation_setup_init(cfg, setup)` and receives a fully prepared simulation state -- confinement operators, strain, SC-modified profile, eigensolver, and solver workspace. The setup type dispatches on geometry (bulk/QW/wire) internally, so applications do not need to know which confinement operators are dense or sparse. Applications only provide their domain-specific solve logic (k-sweep, g-factor extraction, optical accumulation).

New physics combinations (optics + strain + SC, g-factor + strain for QW, optics of doped wire) become available to all apps without code changes beyond calling the setup.

## User Stories

1. As a physics researcher, I want to compute optical absorption of a strained, doped quantum well, so that I can predict laser gain spectra for device design.
2. As a physics researcher, I want to compute the Landau g-factor of a strained quantum well, so that I can compare against experimental magneto-optical measurements.
3. As a physics researcher, I want to compute optical spectra of a doped wire, so that I can study intersubband transitions in 1D-confined systems with self-consistent charge.
4. As a codebase maintainer, I want the setup pipeline (confinement init, strain, SC, eigensolver config) to exist in one place, so that a bug fix applies everywhere.
5. As a codebase maintainer, I want each app to only contain its domain-specific solve logic, so that the apps are short and readable.
6. As a codebase maintainer, I want new physics features (e.g., Dresselhaus spin-orbit coupling) to be added to one setup module, so that all apps gain the feature automatically.
7. As a developer, I want to write a new application (e.g., transport calculation) without reimplementing confinement, strain, and SC setup, so that I can focus on the transport physics.
8. As a developer, I want the setup module to print progress messages, so that I can monitor long-running SC iterations.
9. As a developer, I want the setup type to own and clean up all allocated resources, so that I do not need to track deallocation manually in each app.
10. As a developer, I want per-thread LAPACK workspace for OpenMP-parallel QW k-sweeps, so that diagonalization scales across cores.
11. As a developer, I want a unified velocity matrix builder, so that g-factor and optics apps share the same construction logic.
12. As a developer, I want `read_and_setup` to be a pure parser, so that I can test physics setup without needing an input file on disk.
13. As a developer, I want `solve_kpoint` to work for any geometry, so that k-sweep loops are geometry-agnostic.
14. As a developer, I want the setup type to provide dimension helpers for k-sweep storage, so that I do not repeat the `NUM_VB_STATES*fdStep - numvb + 1` calculation.
15. As a developer, I want Landau mode to remain outside the setup type, so that its special B-sweep fan diagram logic does not complicate the common path.
16. As a developer, I want the adapter pattern for Hamiltonian construction to remain a future possibility, so that Candidate 4 (unified Hamiltonian block structure) can be implemented later without changing the setup API.

## Implementation Decisions

### Module structure

A new module `simulation_setup_mod` in `src/core/simulation_setup.f90`. It depends on: `definitions`, `parameters`, `strain_solver`, `confinement_init`, `hamiltonianConstructor`, `hamiltonian_wire`, `sc_loop`, `eigensolver`, `sparse_matrices`, `linalg`. It is consumed by the three app programs and any future apps.

### Type design: fat type with allocatable components

The `simulation_setup` type uses a single fat derived type (not polymorphic subtypes) with allocatable components for both the dense (bulk/QW) and sparse (wire) paths. Only one path's components are allocated based on `cfg%confinement`. This avoids abstract types and deferred procedures while keeping the allocated footprint lean.

Key fields:
- **Common**: `confinement` (0/1/2), `grid`, `strain_blocks`, `has_strain`, `sc_was_run`, `fermi_level`
- **Dense path** (bulk/QW): `profile(:,:)`, `kpterms(:,:,:)`, `N`, `il`, `iuu`, `HT(:,:)`, LAPACK workspace arrays (`work`, `rwork`, `iwork`, `ifail`, `lwork`)
- **Sparse path** (wire): `profile_2d(:,:)`, `kpterms_2d(:)` CSR, `HT_csr`, `Ngrid`, `Ntot`, `nev_wire`, `eigen_cfg`, `eigen_solver`, `coo_cache`, `wire_ws`
- **Velocity matrices** (opt-in): `vel(3)` CSR, `vel_built` flag

### Init routine: one-shot, driven by cfg flags

`simulation_setup_init(cfg, setup)` performs all steps in sequence:

1. Dispatch on `cfg%confinement`:
   - **Bulk (0)**: Set `N=8`, no confinement init needed
   - **QW (1)**: Call `confinementInitialization(cfg, profile, kpterms)`, compute `N`, `il`, `iuu`, allocate LAPACK workspace, run workspace query
   - **Wire (2)**: Call `confinementInitialization_2d(cfg%grid, cfg%params, cfg%regions, profile_2d, kpterms_2d, cfg%FDorder)`, compute `Ngrid`, `Ntot`, `nev_wire`, configure `eigen_cfg`, create eigensolver
2. If `cfg%strain%enabled`: call `compute_strain` + `compute_bir_pikus_blocks`, store in `setup%strain_blocks`
3. If `cfg%sc%enabled`: call `self_consistent_loop` (QW) or `self_consistent_loop_wire` (wire), modifying profile in-place
4. For wire after SC: free COO cache and workspace (profile changed)

The init routine prints progress messages at each step (consistent with current app output).

### `read_and_setup` becomes a pure parser

The existing `read_and_setup` subroutine in `input_parser.f90` is split:
- The parsing + validation logic (first ~1365 lines) becomes `read_config(cfg)` -- reads `input.cfg`, populates `cfg`, returns.
- The confinement initialization and electric field setup (last ~20 lines) move into `simulation_setup_init`.
- The old `read_and_setup` signature remains as a thin wrapper for backward compatibility during transition.

### `build_H`: dispatch wrapper (not adapter pattern)

`setup_build_H(setup, kvec, HT_out)` dispatches on `setup%confinement`:
- Bulk: `ZB8bandBulk`
- QW: `ZB8bandQW`
- Wire: `ZB8bandGeneralized`

This is a select-case wrapper, not a polymorphic adapter. The adapter pattern (Candidate 4 from the architecture review) can be introduced later by replacing the dispatch body without changing the API.

### `solve_kpoint`: geometry-agnostic solve

`setup_solve_kpoint_serial(setup, kvec, evals, evecs, ws)` builds H and diagonalizes:
- QW/bulk: uses per-thread `thread_workspace` for LAPACK scratch
- Wire: uses `setup%eigen_solver%solve`

The `thread_workspace` type holds per-thread LAPACK arrays (`HT_loc`, `work_loc`, `rwork_loc`, `iwork_loc`, `ifail_loc`). Apps allocate one per OpenMP thread.

### Velocity matrices: opt-in

`setup_build_velocity_matrices(setup, cfg)` builds `setup%vel(3)` for any geometry:
- Bulk: three calls to `ZB8bandBulk(g='g')`, convert to CSR
- QW: commutator-based z-velocity via `build_velocity_matrices`, `ZB8bandQW(g='g')` for x/y
- Wire: `build_velocity_matrices` for x/y, `ZB8bandGeneralized(g='g3')` for z

Only called by apps that need velocity matrices (gfactor, optics). The unified `vel(3)` field means consumers don't branch on geometry.

### K-sweep storage helper

`setup_alloc_sweep(setup, npts, eig, eigv)` allocates correctly-sized eigenvalue and eigenvector arrays for the full k-sweep, using the setup type's pre-computed dimensions (`N`, `il`, `iuu`, `nev_wire`). The arrays are owned by the caller.

### Cleanup

A public `simulation_setup_free(setup)` deallocates all components. A finalizer delegates to it. This follows the project convention of explicit `*_free` routines + finalizers.

### Landau mode exclusion

Confinement mode 3 (Landau) is not handled by the setup type. `main.f90` continues to special-case Landau with its own init path (`confinementInitialization_landau`, `ZB8bandLandau`, B-sweep fan diagram). If Landau gains SC or strain support in the future, it can be added to the setup type.

### Wire branch tracking exclusion

The `reorder_wire_branches` subroutine stays as a `contains` procedure in `main.f90`. Only `bandStructure` uses it for continuous band tracing across k-points. It is not part of the setup pipeline.

### File structure after implementation

- `src/core/simulation_setup.f90` -- new module
- `src/io/input_parser.f90` -- `read_config` extracted, `read_and_setup` becomes wrapper
- `src/apps/main.f90` -- refactored to use setup type (minus Landau path)
- `src/apps/main_gfactor.f90` -- refactored to use setup type
- `src/apps/main_optics.f90` -- refactored to use setup type (gains SC + strain)
- `CMakeLists.txt` -- add new source file

## Testing Decisions

### What makes a good test

Tests should verify external behavior of the setup module: given a config, does the setup produce correct initialized state? Tests should not assert on internal dispatch details. The existing pFUnit test framework and regression test infrastructure are used.

### Modules to test

1. **`simulation_setup_mod`** (new, highest priority):
   - Unit test: `simulation_setup_init` for each geometry (bulk, QW, wire) produces correct dimensions and allocations
   - Unit test: strain setup populates `strain_blocks` when enabled, skips when disabled
   - Unit test: SC setup modifies profile when enabled, skips when disabled
   - Integration test: full pipeline (read config -> setup -> solve_kpoint) produces correct eigenvalues for a known system (use existing regression test configs)

2. **`read_config`** (modified):
   - Existing regression tests cover parsing. Verify they still pass after splitting `read_and_setup`.

3. **App regression tests** (existing):
   - All ~91 existing tests must pass after refactoring. No new test logic needed -- the existing suite covers correctness.

4. **Velocity matrix builder**:
   - Verify `setup_build_velocity_matrices` produces matrices with the correct sparsity pattern and nonzero values for each geometry.

### Prior art

- `tests/unit/sc_loop.pf` -- tests SC loop in isolation, good model for testing setup module
- `tests/unit/strain_solver.pf` -- tests strain computation in isolation
- `tests/regression/configs/` -- existing configs for integration testing
- `validation/` -- cross-code validation against kdotpy for correctness

## Out of Scope

- **Candidate 4 (unified Hamiltonian block structure)**: The adapter pattern for dense/COO builders is deferred. The setup type's `build_H` dispatch wrapper provides the seam for future Candidate 4 integration without changing the API.
- **Candidate 3 (module-level state encapsulation)**: Converting `optical_spectra`'s 14 `save` arrays into a derived type is orthogonal and can be done independently.
- **Candidate 5 (split simulation_config)**: Splitting the god struct into domain configs is speculative and has high blast radius.
- **Landau mode**: Remains a special case in `main.f90`.
- **New executable**: No new app program is created. The refactoring enables future apps but does not add one.
- **OpenMP restructuring**: The parallel k-sweep pattern in the apps remains as-is, just using `thread_workspace` from the setup module instead of manual allocation.

## Further Notes

This PRD addresses Candidate 2 ("Deepen the app orchestration") from the 2026-05-25 architecture review. It is the highest-leverage deepening after Candidate 1 (QW/Wire duplication collapse, partially completed in the same session).

The design decisions were validated through a 19-question grilling session with the codebase owner. All decisions have user sign-off.

The implementation should be done incrementally: first create the setup module with the type and init routine, then refactor one app at a time (starting with `main_gfactor.f90` as the simplest, then `main_optics.f90` as the one that gains the most, then `main.f90` as the most complex). Each app refactoring should be tested against the full regression suite before moving to the next.
