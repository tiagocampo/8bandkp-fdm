# Issue 2: Eliminate `ngrid` — replace with `grid%npoints()`

**Type**: AFK
**Blocked by**: Issue 1 (both modify `simulation_config` and `input_parser.f90`)
**User stories**: 1, 3, 7, 8, 9

## What to build

Remove the `ngrid` field from `simulation_config`. Replace all ~20 consumers with `cfg%grid%npoints()`, which already holds the correct spatial DOF count for all modes (bulk=1, QW=fd_step, wire=nx*ny, landau=landau.nx). Remove the four mode-dependent `cfg%ngrid = ...` overwrite assignments in `input_parser.f90`.

The user-facing TOML key `fd_step` is unchanged — it remains parsed and stored as `cfg%fd_step`, consumed by `init_grid_from_config` for QW grid sizing.

Wire code paths never used `cfg%ngrid` (they already use `grid_ngrid()` directly), so this is purely a cleanup of bulk/QW/Landau paths. It also fixes the dormant wire bug where `cfg%ngrid` was `ny` instead of `nx*ny`.

Add a unit test for `conf_direction()` (from Issue 1) alongside the grid-related test updates.

## Acceptance criteria

- [ ] `ngrid` field removed from `simulation_config` type in `defs.f90`
- [ ] All `cfg%ngrid` references replaced with `cfg%grid%npoints()` in `main.f90`, `main_gfactor.f90`, `main_optics.f90`, `simulation_setup.f90`, `confinement_init.f90`, and any other consumers
- [ ] Parser overwrite assignments (`cfg%ngrid = 1`, `cfg%ngrid = cfg%fd_step`, `cfg%ngrid = cfg%wire%ny`, `cfg%ngrid = cfg%landau%nx`) removed from `input_parser.f90`
- [ ] `validate()` and `validate_semantic()` updated to use `cfg%grid%npoints()`
- [ ] Unit test added for `conf_direction()` covering all four modes
- [ ] All 34+ unit tests pass
- [ ] All regression tests pass (`ctest --test-dir build`)
- [ ] No remaining references to `cfg%ngrid` in `src/` (grep confirms zero hits, excluding comments)
- [ ] `fd_step` TOML key remains parseable — no config file changes needed

## Blocked by

- Issue 1 (both modify `simulation_config` and `input_parser.f90`)
