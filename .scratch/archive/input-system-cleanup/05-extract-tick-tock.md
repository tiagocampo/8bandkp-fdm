# Issue 5: Extract tick/tock timing utilities to utils module

## Parent

PRD: Input System Architecture Cleanup (`.scratch/input-system-cleanup/PRD.md`)

## What to build

Move the `tick` and `tock` timing subroutines from defs.f90 to utils.f90. These utilities (`tick(t)` captures `cpu_time`, `tock(t, label)` prints elapsed time) have no dependency on `simulation_config` or any other type in defs.f90. They belong in utils.f90 alongside other standalone utility routines (dense-to-sparse conversion, Simpson integration).

Update `use` statements in consuming modules to import `tick`/`tock` from utils instead of definitions. This is a pure reorganization — no logic changes.

## Acceptance criteria

- [ ] `tick` and `tock` subroutines are defined in utils.f90, not defs.f90
- [ ] All `use definitions, only: ... tick, tock` imports updated to `use utils, only: tick, tock`
- [ ] defs.f90 no longer contains timing utilities
- [ ] All unit tests pass
- [ ] All ~29 regression tests pass (timing output is cosmetic, not compared in golden files)

## Blocked by

- Issue 1 (Eliminate legacy aliases) — both modify defs.f90; avoid merge conflicts.
