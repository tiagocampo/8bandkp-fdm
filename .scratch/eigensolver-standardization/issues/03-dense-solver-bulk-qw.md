# Dense solver + bulk/QW through unified interface

**Type:** AFK

## Parent

PRD: `.scratch/eigensolver-standardization/PRD.md`

## What to build

Complete the dense solver end-to-end path: implement `dense_lapack_solver_t` with mode dispatch (FULL→`zheev`, INDEX→`zheevx range='I'`, ENERGY→`zheevx range='V'`), route bulk and QW paths in `main.f90` through `solver%solve_dense`, set smart defaults for bulk (DENSE+FULL) and QW (DENSE+INDEX), and allocate `eigen_solver` for bulk/QW in `simulation_setup_init`.

This slice extends the architecture established in #02 to cover the dense LAPACK solver and the two most common dense-geometry paths (bulk and QW k-sweeps). The OpenMP parallel k-sweep pattern is refactored: instead of thread-private LAPACK workspace arrays, each OMP thread gets its own polymorphic solver instance via `make_eigensolver`.

**`dense_lapack_solver_t%solve_dense`** dispatches on `config%mode`:
- `EIGEN_MODE_FULL` → `zheev` (all eigenvalues/vectors), workspace cached in solver object
- `EIGEN_MODE_INDEX` → `zheevx range='I'` with `il=config%il, iu=config%iu`
- `EIGEN_MODE_ENERGY` → `zheevx range='V'` with `vl=config%emin, vu=config%emax`

**`simulation_setup.f90`**: Allocate `eigen_solver` for bulk and QW modes (currently only wire allocates it). Set up `eigensolver_config` from `cfg%solver` with smart defaults: bulk→DENSE+FULL, QW→DENSE+INDEX. Compute `il`/`iu` from `cfg%bands%num_cb/num_vb` for INDEX mode (existing formulas from main.f90).

**`main.f90` bulk path**: Replace inline `zheevx` workspace query + solve with `solver%solve_dense(H_dense, config, result)`. Bulk is serial (8×8 trivial).

**`main.f90` QW path**: Replace the OMP parallel k-sweep's inline `zheevx` with per-thread solver instances. Each thread creates its own `make_eigensolver(setup%eigen_cfg)`, calls `solver%solve_dense(H_loc, config, result)` per k-point, frees result after processing.

**`main_gfactor.f90` bulk/QW**: `setup_solve_kpoint_serial` already handles bulk/QW via simulation_setup — ensure it uses the new solver dispatch.

## Acceptance criteria

- [ ] `dense_lapack_solver_t%solve_dense` dispatches on EIGEN_MODE_FULL/INDEX/ENERGY with correct LAPACK routine
- [ ] `dense_lapack_solver_t%solve_sparse` converts CSR→dense then calls solve_dense
- [ ] `eigen_solver` is allocated for bulk and QW modes in `simulation_setup_init`
- [ ] Bulk path in `main.f90` uses `solver%solve_dense` (no direct `zheevx` call)
- [ ] QW k-sweep in `main.f90` uses per-thread solver instances in OMP parallel region (no direct `zheevx` call)
- [ ] Smart defaults work: bulk→DENSE+FULL, QW→DENSE+INDEX
- [ ] `il`/`iu` computed from `cfg%bands` for INDEX mode (same formulas as before)
- [ ] No direct `zheev`/`zheevd`/`zheevx` calls remain in `main.f90` bulk/QW paths
- [ ] Unit tests: `test_dense_full_mode`, `test_dense_index_mode`, `test_dense_energy_mode` pass
- [ ] All existing bulk and QW regression tests pass
- [ ] `cmake --build build` succeeds
- [ ] `ctest --test-dir build -j4 --output-on-failure` passes all tests

## Blocked by

- #02 (unified interface + solver_config must exist)
