# PRD: Input System Architecture Cleanup

## Problem Statement

The TOML parser refactor replaced the 1,369-line sequential parser with a clean 925-line TOML-based parser, addressing 4 of 5 architecture review candidates. However, the migration left behind legacy alias fields in `simulation_config`, scattered validation across 3 layers, implicit ordering dependencies between parse routines, an executable that bypasses the shared initialization pipeline, and silent fallthrough on type-mismatched TOML input. These remnants create maintenance debt: every field write is doubled, validation rules must be updated in multiple places, and the bridge code in confinement_init contains defensive fallback logic for aliases that should no longer exist.

## Solution

Consolidate the remaining architectural debt from the TOML parser migration into a single cleanup pass: eliminate legacy alias fields, consolidate validation into two clear layers, fix the b_field/BdG ordering dependency, route all executables through the shared initialization pipeline, and add diagnostics for type-mismatched TOML input.

## User Stories

### Legacy alias elimination (F1+F2)

1. As a developer, I want `simulation_config` to have one canonical name per field, so that I never have to remember whether to use `materialN` or `material_names`.
2. As a developer, I want `confinementInitialization_raw` to accept the new field names as parameters, so that the bridge code does not need fallback logic.
3. As a developer, I want the finalizer in defs.f90 to deallocate only the canonical fields, so that cleanup is simple and complete.
4. As a developer, I want `validate()` to check `material_names` not `materialN`, so that validation uses the canonical names.
5. As a developer, I want `sc_loop` to read `int_start_pos`/`int_end_pos` instead of `intStartPos`/`intEndPos`, so that the physics modules use the same naming convention as the parser.
6. As a developer, I want `main_gfactor` to read `z_min`/`z_max` instead of `startPos`/`endPos`, so that all executables use TOML-native field names.

### Validation consolidation (F3)

7. As a developer, I want range-bounds validation (fd_step >= 3, wire nx >= FDorder+1, etc.) to live in exactly one place, so that changing a constraint requires updating only one check.
8. As a developer, I want `validate()` in defs.f90 to be the single source of truth for structural and state invariants, so that I can trust a config that passes validation.
9. As a developer, I want app-specific constraints (g-factor requires k0 mode, optics requires [optics] section) to live in a `validate_semantic(cfg, app_name)` routine, so that parser-level validation does not need to know which executable is running.
10. As a developer, I want the parser to keep only structural checks (required sections, TOML parse failures, enum validation), so that parsing and validation have clear separation of concerns.
11. As a user running the solver, I want to see clear error messages when my config violates a semantic constraint (e.g., "g-factor calculation requires wave_vector mode = k0"), so that I can fix my input quickly.

### b_field/BdG ordering fix (F5)

12. As a developer, I want `parse_b_field` to not write to `cfg%bdg` fields, so that each parse routine only populates its own sub-config.
13. As a developer, I want `parse_bdg` to self-containedly set its defaults from b_field before applying TOML overrides, so that the data flow is explicit within one routine.
14. As a developer, I want to be able to reorder the `parse_*` calls in `read_config` without breaking BdG default values, so that the orchestrator has no implicit ordering constraints.

### simulation_setup skip_sc flag (F6)

15. As a developer, I want `simulation_setup_init` to accept an optional `skip_sc` argument, so that callers can explicitly opt out of SC rather than mutating `cfg%sc%enabled`.
16. As a developer, I want `main.f90` to call `simulation_setup_init(cfg, setup, skip_sc=.true.)` instead of saving/restoring `cfg%sc%enabled`, so that config mutation is not used for control flow.

### Route main_topology through simulation_setup (F7)

17. As a developer, I want `main_topology` to call `simulation_setup_init` instead of inlining QW initialization, so that topology benefits from strain integration and future setup improvements.
18. As a developer, I want the z=0 electric field check in main_topology to be handled by simulation_setup_init, so that validation is not duplicated across executables.
19. As a developer, I want all four executables to share the same initialization pipeline, so that changes to QW setup only need to be made in one place.

### Type-mismatch diagnostics (F8)

20. As a user writing TOML configs, I want to see a warning when my value has the wrong type (e.g., a string where a number is expected), so that I know my input is being ignored.
21. As a developer, I want a reusable `check_optional_stat` helper that warns on stat /= 0 for optional field lookups, so that I don't have to write the same boilerplate for every optional key.
22. As a developer, I want required-field lookups to keep using the existing `require_string`/`require_int`/`require_real` helpers (which error-stop), so that missing required fields are hard errors.
23. As a user, I want warnings to include the section name and key name, so that I can quickly locate the problem in my config file.

### Grid and timing extraction (F4, partial)

24. As a developer, I want `tick`/`tock` timing utilities to live in utils.f90, so that defs.f90 focuses on type definitions and validation.
25. As a developer, I want `init_grid_from_config` to remain in defs.f90 since it couples tightly to simulation_config internals, so that I avoid circular dependencies or parameter bags.

## Implementation Decisions

### Decision 1: Rename-in-place for legacy aliases

All legacy alias fields (`materialN`, `startPos`, `endPos`, `intStartPos`, `intEndPos`) are removed from `simulation_config`. Consumers are updated to use the canonical names (`material_names`, `z_min`/`z_max`, `int_start_pos`/`int_end_pos`). The `confinementInitialization_raw` interface is updated to accept the new parameter names. This is a mechanical rename across ~10 files, done in a single commit.

No two-phase migration. The test suite (34 unit tests + ~29 regression tests) provides immediate feedback on any missed rename.

### Decision 2: Two-layer validation

Validation is consolidated into two layers:

- **Syntactic validation** (input_parser.f90): TOML parse failures, required section presence, enum validation, shape validation. These fail fast with `print + stop 1`.
- **Semantic validation** (defs.f90): A `validate_semantic(cfg, app_name)` subroutine that checks structural invariants (field ranges, array allocations, cross-field consistency) AND app-specific constraints dispatched on `app_name` (e.g., g-factor needs k=0, optics needs [optics] section).

The overlapping range checks currently duplicated in both the parser and `validate()` are removed from the parser. `validate()` remains the single source of truth for "is this config valid?"

The `validate_semantic` interface:

```fortran
subroutine validate_semantic(cfg, app_name)
  type(simulation_config), intent(in) :: cfg
  character(len=*), intent(in) :: app_name  ! 'bandStructure', 'gfactor', 'opticalProperties', 'topologicalAnalysis'
```

App-specific checks currently in main*.f90 that move to `validate_semantic`:
- g-factor: wave_vector%nsteps == 0 (k0 mode required)
- optics: optics%enabled must be true
- topology: topo%enabled must be true, topo%mode must be non-empty

### Decision 3: b_field copy moved into parse_bdg

The unconditional copy from `parse_b_field` (`cfg%bdg%B_vec = cfg%b_field%components` and `cfg%bdg%g_factor = cfg%b_field%g_factor`) moves to the top of `parse_bdg` as default initialization, before the `[bdg]` TOML section is read. This makes parse_bdg self-contained: it sets defaults from b_field, then overrides with explicit `[bdg]` values if present.

### Decision 4: skip_sc optional argument for simulation_setup_init

`simulation_setup_init` gains `skip_sc` optional logical argument. When present and true, the routine skips all SC-related initialization. This replaces the save/restore pattern in main.f90 where `cfg%sc%enabled` was temporarily zeroed.

```fortran
subroutine simulation_setup_init(cfg, setup, strain_out, sc_phi_out, sc_ne_out, sc_nh_out, skip_sc)
  type(simulation_config), intent(inout) :: cfg
  type(simulation_setup), intent(inout) :: setup
  type(strain_result), intent(out), optional, allocatable :: strain_out
  real(kind=dp), allocatable, intent(out), optional :: sc_phi_out(:,:), sc_ne_out(:,:), sc_nh_out(:,:)
  logical, intent(in), optional :: skip_sc
```

### Decision 5: Route main_topology through simulation_setup_init

main_topology.f90's inline QW initialization block (confinementInitialization, electric field check, external field setup) is replaced with `call simulation_setup_init(cfg, setup, skip_sc=.true.)`. This makes all four executables share the same initialization pipeline, bringing strain integration to topology for free.

### Decision 6: Type-mismatch warning helper

A private helper `check_optional_stat(stat, key, section)` prints a warning when `stat /= 0` on optional field lookups. Applied incrementally to optional `get_value` calls (those with default values), not to required-field lookups (which use `require_*` helpers that already error-stop). The helper is defined once and used throughout all parse_* routines.

### Decision 7: tick/tock moved to utils.f90

The `tick` and `tock` timing subroutines are moved from defs.f90 to utils.f90, with `use utils, only: tick, tock` added to consuming modules. This reduces defs.f90's responsibility to type definitions, constants, grid initialization, and validation.

### Decision 8: init_grid_from_config stays in defs.f90

`init_grid_from_config` couples tightly to `simulation_config` internals (reads cfg%confinement, cfg%z, cfg%int_start_pos, cfg%wire%nx, etc.). Extracting it would require either importing the full config type into another module or creating a parameter bag. For 70 lines, the extraction cost exceeds the benefit. It stays.

## Testing Decisions

### What makes a good test

Tests verify external behavior, not implementation details. A config that passes validation should produce correct results when consumed by any executable. A config that violates a constraint should produce a clear error message. Field rename tests verify that the physics results are unchanged, not that a specific variable name was used.

### Modules to test

- **validate_semantic** — New unit tests in pFUnit for each app-specific constraint: g-factor with non-k0 mode (should error), optics without [optics] section (should error), topology without [topology] section (should error), valid configs for each app (should pass).
- **skip_sc flag** — Unit test calling `simulation_setup_init` with `skip_sc=.true.` and verifying SC was not run (sc_was_run is false).
- **Field rename regression** — All 34 existing unit tests and ~29 regression tests must pass without modification (except updating verifier field-name references if any check legacy names). This proves the rename is behavior-preserving.
- **b_field/BdG ordering** — Existing regression tests with BdG configs should pass. The TOML override of b_field defaults by [bdg] B_vec should still work.

### Prior art

- pFUnit unit tests in `tests/unit/` follow the pattern of calling a subroutine with test inputs and asserting on outputs or error conditions.
- Regression tests in `tests/regression/` compare numerical output against golden reference data.
- The TOML parser refactor used the same approach: mechanical rename verified by full test suite pass.

## Out of Scope

- **Extracting physics sub-types from defs.f90** — Moving `topology_config`, `optics_config`, `bdg_config`, etc. into their respective physics modules would create circular dependencies (those modules already `use definitions`). This is deferred until defs.f90 grows past 1000 lines or until a forward-type pattern is adopted.
- **Extracting init_grid_from_config from defs.f90** — Tightly coupled to simulation_config internals. Extraction cost exceeds benefit at 70 lines.
- **Fixing pre-existing OpenMP reentrancy bug** — The "Recursive call to nonrecursive procedure" in hamiltonianConstructor.f90 affects 4+ regression tests but is unrelated to the input system.
- **Fixing wire kz NaN with nsteps=1** — Division by zero in main.f90:90 is a separate bug.
- **Adding strict parsing mode** — A flag to make type mismatches into errors (not just warnings) is deferred. The warning helper (F8) is the incremental first step.

## Further Notes

### Execution order

The work items have the following dependency chain:

1. **F1+F2 (legacy aliases)** — First. Unblocks F3 and F4 since validate() and grid init reference the alias fields.
2. **F5 (b_field copy)** — Independent, quick win.
3. **F4 (tick/tock extraction)** — After F1 since both touch defs.f90.
4. **F3 (validation consolidation)** — After F1 since validate() references the renamed fields.
5. **F6+F7 (skip_sc + topology routing)** — Independent of F1-F5.
6. **F8 (warning helper)** — Lowest priority, purely additive.

Items 2, 5, and 6 can potentially run in parallel after item 1 lands.

### ADR alignment

This PRD aligns with ADR-0001 (simulation_setup fat type). The `skip_sc` optional argument extends the existing `simulation_setup_init` interface without changing the fat-type pattern. Routing main_topology through the shared init pipeline completes the vision of "one type, one init routine" across all four executables.

### Blast radius summary

| Work item | Files touched | Risk |
|---|---|---|
| F1+F2 (aliases) | ~10 files, mechanical rename | Medium — test suite catches missed renames |
| F3 (validation) | 6 files | Low — moves checks, doesn't change logic |
| F4 (timing) | 2 files | Low — pure reorganization |
| F5 (b_field copy) | 1 file, 3 lines moved | Low |
| F6+F7 (skip_sc + topology) | 3 files | Low-medium — topology gets strain for free but needs testing |
| F8 (warnings) | 1 file, additive | Minimal |
