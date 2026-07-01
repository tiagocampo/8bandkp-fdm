Status: ready-for-agent

## Parent

PRD: `.scratch/simulation-setup-orchestration/PRD.md`

## What to build

Create the `simulation_setup` module (`src/core/simulation_setup.f90`) and refactor the gfactor app to use it. This is the thinnest vertical slice: a working setup type that serves one real app end-to-end.

The module provides a fat derived type `simulation_setup` with allocatable components for both dense (bulk/QW) and sparse paths (only dense is populated in this slice). The public API is:

- `simulation_setup_init(cfg, setup)` — dispatches on `cfg%confinement`, runs confinement init, strain (if enabled), SC (if enabled). Bulk sets N=8, no confinement init. QW calls `confinementInitialization`, computes N/il/iuu, allocates LAPACK workspace.
- `setup_build_H(setup, kvec, HT_out)` — select-case dispatch: bulk → `ZB8bandBulk`, QW → `ZB8bandQW`.
- `setup_solve_kpoint_serial(setup, kvec, evals, evecs, ws)` — builds H and diagonalizes via LAPACK (bulk/QW).
- `setup_build_velocity_matrices(setup, cfg)` — builds `setup%vel(3)` CSR velocity matrices. Bulk: three `ZB8bandBulk(g='g')` calls, convert to CSR. QW: commutator z-velocity + `ZB8bandQW(g='g')` for x/y.
- `simulation_setup_free(setup)` — deallocates all components. Finalizer delegates to it.

Simultaneously, split `read_and_setup` in `input_parser.f90`: extract a pure `read_config(cfg)` that reads `input.cfg` and populates the config type (the first ~1365 lines of parsing + validation). The old `read_and_setup` signature becomes a thin wrapper during transition. The confinement initialization and electric field setup (last ~20 lines) move into `simulation_setup_init`.

Then refactor `main_gfactor.f90` to call `read_config` → `simulation_setup_init` → `setup_build_velocity_matrices` → use the setup type for its g-factor extraction loop → `simulation_setup_free`. Remove the app's direct calls to `confinementInitialization`, strain, SC, and velocity matrix construction.

Unit tests in `tests/unit/test_simulation_setup.pf`: verify `simulation_setup_init` produces correct dimensions (N=8 for bulk, correct N/il/iuu for QW with a test config), verify `setup_build_H` returns a matrix of the expected size for each geometry. Use a simple GaAs bulk config.

All existing gfactor regression tests must pass after refactoring. Add `simulation_setup.f90` to `CMakeLists.txt`.

## Acceptance criteria

- [ ] `src/core/simulation_setup.f90` exists with `simulation_setup` type, init, build_H, solve_kpoint, velocity matrices, free+finalizer
- [ ] `read_config` extracted from `read_and_setup` in `input_parser.f90`; `read_and_setup` remains as backward-compat wrapper
- [ ] `main_gfactor.f90` refactored to use setup module — no direct confinement/strain/SC/velocity calls remain
- [ ] `test_simulation_setup.pf` covers: bulk init dimensions, QW init dimensions, build_H matrix size for each geometry
- [ ] All existing gfactor regression tests pass
- [ ] `CMakeLists.txt` updated with new source file

## Blocked by

None - can start immediately.
