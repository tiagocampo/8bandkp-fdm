# Issue 1: Eliminate legacy alias fields from simulation_config

## Parent

PRD: Input System Architecture Cleanup (`.scratch/input-system-cleanup/PRD.md`)

## What to build

Remove the 6 legacy alias fields from `simulation_config` (`materialN`, `startPos`, `endPos`, `intStartPos`, `intEndPos`) and update all consumers to use the canonical names (`material_names`, `z_min`/`z_max`, `int_start_pos`/`int_end_pos`). Rename the `confinementInitialization_raw` interface parameters to use the new names. Remove the defensive fallback logic in `confinement_init.f90` that copies new→old when old is unallocated. Remove the parallel alias writes in `input_parser.f90` (the lines like `cfg%materialN(i) = trim(name_val)` alongside `cfg%material_names(i) = trim(name_val)`). Remove the alias field deallocations from the `simulation_config` finalizer. Update `validate()` to check `material_names` instead of `materialN`.

This is a mechanical rename across ~10 files. No logic changes — only field names change.

## Acceptance criteria

- [ ] `materialN`, `startPos`, `endPos`, `intStartPos`, `intEndPos` fields no longer exist in `simulation_config`
- [ ] `confinementInitialization_raw` dummy parameters use `z_min`/`z_max`/`material_names`/`int_start_pos`/`int_end_pos`
- [ ] No defensive fallback logic remains in `confinement_init.f90` (the `if (allocated(cfg%int_start_pos) .and. .not. allocated(cfg%intStartPos))` block is gone)
- [ ] `validate()` checks `material_names` not `materialN`
- [ ] Finalizer deallocates only canonical fields
- [ ] `input_parser.f90` writes only canonical fields (no parallel alias writes)
- [ ] All 34 unit tests pass
- [ ] All ~29 regression tests pass (some pre-existing failures are unrelated)
- [ ] `grep -rn 'materialN\|startPos\|endPos\|intStartPos\|intEndPos' src/` returns zero matches (excluding comments and string literals in output messages)

## Blocked by

None — can start immediately.
