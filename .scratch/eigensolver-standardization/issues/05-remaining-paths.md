# Remaining paths — Landau, optics, SC, topology, spectral through unified interface

**Type:** AFK

## Parent

PRD: `.scratch/eigensolver-standardization/PRD.md`

## What to build

Route all remaining eigensolve call sites through the polymorphic dispatch established in #02 and #03. No new architecture — this is pure pattern-following. Each path is migrated individually from direct LAPACK calls to `solver%solve_dense` or `solver%solve_sparse`.

**Paths to migrate:**

- **Landau k-sweep** (`main.f90`): Replace inline `zheevx` in OMP parallel k-sweep with per-thread solver instances. Smart default: DENSE + INDEX. Landau B-sweep: replace `zheevx('N',...)` with solver dispatch.
- **Optics bulk k-sweep** (`main_optics.f90`): Replace inline `zheevx` in OMP parallel sweep. Smart default: DENSE + INDEX.
- **Optics QW k-sweep** (`main_optics.f90`): Same pattern. Smart default: DENSE + INDEX.
- **QW SC loop** (`sc_loop.f90`): Replace inline `zheevx` in nested SC×k_par loop. Workspace previously pre-allocated by SC loop — now managed by solver object.
- **Wire SC loop** (`sc_loop.f90`): Already uses polymorphic FEAST — update from `cfg%feast` to `cfg%solver` if not done in #02.
- **QW BdG** (`main_topology.f90`): Replace direct `zheev` with `solver%solve_dense`.
- **QW Z2 TRIM loop** (`topological_analysis.f90`): Replace direct `zheev` at 4 TRIM points with solver dispatch.
- **Bulk spectral function** (`green_functions.f90`): Replace `zheev('N',...)` with `solver%solve_dense`.
- **QW spectral function** (`green_functions.f90`): Same pattern.
- **Topology sweep** (`main_topology.f90`): Wire BdG sweep already uses FEAST. QW Fu-Kane sweep: replace repeated `zheev` calls with solver dispatch.

**Simulation setup**: Ensure `eigen_solver` is allocated for Landau mode in `simulation_setup_init`. Landau smart defaults: DENSE + INDEX.

## Acceptance criteria

- [ ] No direct `zheev`, `zheevd`, or `zheevx` calls remain in `main.f90`, `main_optics.f90`, `main_topology.f90`, `sc_loop.f90`, `green_functions.f90`, `topological_analysis.f90`, or `simulation_setup.f90`
- [ ] All eigensolve calls go through `solver%solve_dense` or `solver%solve_sparse`
- [ ] `eigen_solver` is allocated for all 4 confinement modes in `simulation_setup_init`
- [ ] OpenMP parallel k-sweeps work correctly (Landau, optics bulk, optics QW)
- [ ] SC loop works correctly (QW dense, wire FEAST)
- [ ] BdG and Z2 topology paths work correctly
- [ ] Spectral function paths work correctly
- [ ] `cmake --build build` succeeds
- [ ] `ctest --test-dir build -j4 --output-on-failure` passes all tests

## Blocked by

- #03 (dense solver + bulk/QW paths established)
- #04 (QW CSR builder for QW sparse path)
