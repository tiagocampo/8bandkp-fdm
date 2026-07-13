Status: ready-for-agent

## Parent

PRD: `.scratch/simulation-setup-orchestration/PRD.md`

## What to build

Refactor the optics app to use the setup module, unlocking SC and strain support for optical calculations for the first time.

The velocity matrix builder (`setup_build_velocity_matrices`) was created in issue 01 for bulk/QW paths. This issue consumes it in the optics app. Refactor `main_optics.f90` to call `read_config` → `simulation_setup_init` → `setup_build_velocity_matrices` → use the setup type's `vel(3)` for optical accumulation → `simulation_setup_free`. Remove the app's direct calls to confinement init, strain, and velocity matrix construction.

Optics now automatically gains SC and strain support because `simulation_setup_init` handles those pipelines internally. No optics-specific SC or strain code is needed — the setup module dispatches based on `cfg` flags.

Unit tests in `test_simulation_setup.pf`: verify `setup_build_velocity_matrices` produces CSR matrices with correct sparsity pattern (same nonzero count as manually-constructed velocity matrices) for bulk GaAs and a QW config.

All existing optics regression tests must pass. The new physics combinations (optics + SC, optics + strain) are implicitly available but not explicitly tested in this slice — they will be tested when users exercise those configs.

## Acceptance criteria

- [ ] `main_optics.f90` refactored to use setup module — no direct confinement/strain/velocity calls remain
- [ ] Optics app uses `setup%vel(3)` from setup module for all velocity matrix access
- [ ] Velocity matrix unit tests cover bulk and QW sparsity patterns
- [ ] All existing optics regression tests pass
- [ ] Optics app gains SC and strain support through setup module (no new optics-specific SC/strain code)

## Blocked by

- Issue 01 (setup module skeleton + gfactor migration)

## User stories covered

1, 3, 11
