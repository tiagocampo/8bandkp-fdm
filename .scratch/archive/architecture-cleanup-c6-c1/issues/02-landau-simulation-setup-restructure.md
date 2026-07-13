# Close Landau Gap in Simulation Setup + Restructure main.f90 (C1)

**Type:** AFK
**Blocked by:** #01 (C6 must land first — shares `simulation_setup.f90` and `main.f90`)

## What to Build

Add Landau as a fourth confinement mode to the `simulation_setup` module and restructure `main.f90` so Landau becomes a self-contained block. Commit as `refactor: close Landau gap in simulation_setup`.

### Part A: Add `case('landau')` to `simulation_setup_init`

Landau is single-material only — no strain, no electric field, no self-consistent loop. The branch is simple (~20 lines):

- Set `setup%N = cfg%landau%nx * 8`
- Compute `setup%il` and `setup%iuu` from band selection (same formula as QW: `NUM_VB_STATES * nx - num_vb + 1` to `NUM_VB_STATES * nx + num_cb`)
- Allocate and zero `setup%kpterms(cfg%landau%nx, cfg%landau%nx, 10)`
- Allocate `setup%profile` (3 columns, like QW)
- Call `confinementInitialization_landau(cfg%grid%x, cfg%material_names(1:1), cfg%params(1:1), setup%profile, setup%kpterms, cfg%FDorder)`
- Allocate `setup%HT(setup%N, setup%N)`, zero it
- Perform zheevx workspace query: build H at k=0 via `ZB8bandLandau`, call `zheevx` with `lwork=-1`, read optimal lwork, reallocate `setup%work`
- Allocate `setup%rwork(7*N)`, `setup%iwork(5*N)`

No new type fields needed — reuses `setup%profile`, `setup%kpterms`, `setup%HT`, `setup%work`, `setup%rwork`, `setup%iwork`.

The module `simulation_setup.f90` will need `use confinement_init, only: ..., confinementInitialization_landau` added to its imports. Check whether `ZB8bandLandau` is already accessible via the existing `use hamiltonianConstructor` import.

### Part B: Add `case('landau')` to `setup_build_H`

Follow the QW pattern exactly:

```
case('landau')
  if (.not. present(HT_out)) error stop 'setup_build_H: landau requires HT_out'
  call ZB8bandLandau(HT_out, kvec, setup%profile, setup%kpterms, cfg%grid%x, cfg=cfg)
```

Note: `ZB8bandLandau` takes an extra `x_grid` argument that `ZB8bandQW` does not. Use `cfg%grid%x`.

### Part C: Restructure `main.f90`

After the wire block (which ends with `stop` at line 380), add a self-contained Landau block:

1. `if (trim(cfg%confinement) == 'landau') then` — open block
2. Call `simulation_setup_init(cfg, setup)` to get profile, kpterms, workspace
3. Extract `N`, `il`, `iuu`, `work`, `rwork`, `iwork`, `HT`, `profile`, `kpterms` from `setup` into local variables (or use setup components directly)
4. Move the Landau profile output (current lines 468–477) into this block
5. Move the Landau k-sweep (current lines 665–728) into this block
6. Move the Landau B-sweep (current lines 730–810) into this block
7. Call `simulation_setup_free(setup)`
8. `stop  ! Landau mode complete`
9. `end if` — close block

Then remove Landau branches from the shared bulk/QW path:
- Remove `else if (conf_direction(cfg%confinement) == 'x')` blocks from the workspace query (lines 555–570)
- Remove `else if (conf_direction(cfg%confinement) == 'x')` from profile output (lines 468–477)
- Remove `else if (conf_direction(cfg%confinement) == 'x')` from the k-sweep dispatch (lines 665+)
- Remove `if (trim(cfg%confinement) == 'landau')` init block (lines 387–403) — now handled by simulation_setup
- Remove `use confinement_init, only: confinementInitialization_landau` from main.f90 imports — no longer needed

The shared code (lines 408–570 after deletions) becomes purely bulk/QW.

### Behavioral equivalence

The Landau output files (`potential_profile.dat`, eigenvalue files, `landau_fan.dat`) must be bit-identical before and after this refactor. The only change is structural — where the code lives, not what it computes.

## Acceptance Criteria

- [ ] `case('landau')` branch added to `simulation_setup_init` in `simulation_setup.f90`
- [ ] `case('landau')` branch added to `setup_build_H` in `simulation_setup.f90`
- [ ] `confinementInitialization_landau` imported in `simulation_setup.f90`
- [ ] `confinementInitialization_landau` removed from `main.f90` imports
- [ ] `main.f90` restructured: Landau is a self-contained block between wire and bulk/QW
- [ ] All `else if (conf_direction == 'x')` branches removed from shared bulk/QW path
- [ ] `cmake --build build` succeeds
- [ ] `ctest --test-dir build -j4 --output-on-failure` passes all tests
- [ ] Landau output files are bit-identical to pre-refactor output (if Landau regression test exists, verify it passes)
- [ ] `setup_solve_kpoint_serial` has NO `case('landau')` — k-sweep stays inline
- [ ] Single commit with message `refactor: close Landau gap in simulation_setup`

## Blocked by

- #01 (C6 dead code cleanup must land first — shared files)

## User Stories Covered

- #10: simulation_setup handles all 4 modes
- #11: Landau init localized to simulation_setup
- #12: main.f90 structured with self-contained mode blocks
- #13: Shared bulk/QW path free of Landau branches
- #14: setup_build_H dispatches Landau uniformly
- #16: Identical output before and after refactor
