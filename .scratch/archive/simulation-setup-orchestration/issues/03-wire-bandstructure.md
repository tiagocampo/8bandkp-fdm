Status: ready-for-agent

## Parent

PRD: `.scratch/simulation-setup-orchestration/PRD.md`

## What to build

Add wire (confinement=2) support to the setup module and refactor the bandStructure app to use it. This is the most complex slice because bandStructure handles all three geometries plus the Landau special case.

Extend `simulation_setup_init` with the wire path: call `confinementInitialization_2d`, compute `Ngrid`/`Ntot`/`nev_wire`, configure `eigen_cfg`, create eigensolver via `make_eigensolver(config)`. Wire path allocates the sparse components: `profile_2d`, `kpterms_2d` CSR, `HT_csr`, `coo_cache`, `wire_ws`. If SC is enabled for wire, call `self_consistent_loop_wire`, then free COO cache and workspace (profile changed).

Extend `setup_build_H` with wire dispatch: `ZB8bandGeneralized`.

Extend `setup_solve_kpoint_serial` with wire path: uses `setup%eigen_solver%solve`.

Add `thread_workspace` type holding per-thread LAPACK arrays (`HT_loc`, `work_loc`, `rwork_loc`, `iwork_loc`, `ifail_loc`). Apps allocate one per OpenMP thread. Add `setup_alloc_sweep(setup, npts, eig, eigv)` helper that allocates correctly-sized eigenvalue/eigenvector arrays using the setup type's pre-computed dimensions.

Refactor `main.f90` to call `read_config` → `simulation_setup_init` → use setup type for k-sweep loop → `simulation_setup_free`. The Landau path (confinement=3) remains special-cased in main.f90, untouched. The wire branch tracking (`reorder_wire_branches`) stays as a `contains` procedure in main.f90.

All existing bandStructure regression tests must pass, including wire configs with OpenMP.

## Acceptance criteria

- [ ] `simulation_setup_init` handles wire path (confinement=2): confinement init, sparse allocations, eigensolver creation
- [ ] `setup_build_H` dispatches wire via `ZB8bandGeneralized`
- [ ] `setup_solve_kpoint_serial` handles wire via `eigen_solver%solve`
- [ ] `thread_workspace` type with per-thread LAPACK arrays exists
- [ ] `setup_alloc_sweep` allocates correctly-sized sweep arrays
- [ ] `main.f90` refactored to use setup module for bulk, QW, and wire paths
- [ ] Landau path (confinement=3) remains untouched in main.f90
- [ ] `reorder_wire_branches` stays as `contains` procedure in main.f90
- [ ] All existing bandStructure regression tests pass (including wire with OpenMP)

## Blocked by

- Issue 01 (setup module skeleton + gfactor migration)

## User stories covered

2, 10, 14, 15
