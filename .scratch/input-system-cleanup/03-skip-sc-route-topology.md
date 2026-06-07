# Issue 3: Add skip_sc flag and route topology through simulation_setup

## Parent

PRD: Input System Architecture Cleanup (`.scratch/input-system-cleanup/PRD.md`)

## What to build

Two related changes to the initialization pipeline:

**Add `skip_sc` optional argument to `simulation_setup_init`.** Currently main.f90 disables SC by saving `cfg%sc%enabled`, setting it to 0, calling `simulation_setup_init`, then restoring the saved value. This is config mutation for control flow. Instead, add `logical, intent(in), optional :: skip_sc` to `simulation_setup_init`. When present and true, the routine skips all SC-related initialization. Update main.f90 to call `simulation_setup_init(cfg, setup, skip_sc=.true.)`.

**Route main_topology through `simulation_setup_init`.** Currently main_topology.f90 inlines QW initialization (confinementInitialization call, z=0 electric field check, external field setup) instead of using the shared pipeline. Replace this inline block with `call simulation_setup_init(cfg, setup, skip_sc=.true.)`. This gives topology strain integration for free and ensures all four executables share one initialization path.

The `simulation_setup` type is already used by main.f90, main_gfactor.f90, and main_optics.f90. main_topology.f90 needs to `use simulation_setup_mod` and declare a `type(simulation_setup)` variable.

## Acceptance criteria

- [ ] `simulation_setup_init` accepts optional `skip_sc` argument
- [ ] When `skip_sc=.true.`, `setup%sc_was_run` is false and no SC iteration occurs
- [ ] main.f90 uses `skip_sc=.true.` instead of save/restore pattern on `cfg%sc%enabled`
- [ ] main_topology.f90 calls `simulation_setup_init(cfg, setup, skip_sc=.true.)` for QW mode
- [ ] main_topology.f90 no longer inlines confinementInitialization, z=0 check, or external field setup
- [ ] main_topology QW initialization produces identical physics results (verified by topology regression tests)
- [ ] All unit tests pass
- [ ] Topology regression tests pass

## Blocked by

None — can start immediately. Note: may have minor merge overlap with Issue 1 if both modify simulation_setup.f90, but no hard dependency.
