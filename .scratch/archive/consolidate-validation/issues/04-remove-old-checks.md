# Issue 4: Remove duplicated validation from executables and simulation_setup

**Type**: AFK**Blocked by**: Issue 2, Issue 3
**User stories**: 20
**GitHub issue**: #23

## What to build

Remove old validation checks from the four executables and `simulation_setup.f90` that are now covered by the consolidated `validate()` and `validate_semantic()` in `defs.f90`. After issues 2 and 3 land, these checks fire at parse time — the executable-level checks are dead code.

Checks to remove:

**From bandStructure (`main.f90`):**
- V1 old: evnum > 8 bulk cap with warning (replace the silent cap block)
- V2 old: num_cb > max QW cap with warning
- V3 old: num_vb > max QW cap with warning
- V4 old: wave vector mode `stop "no such direction"`
- V5 old: B-sweep step <= 0 check

**From gfactorCalculation (`main_gfactor.f90`):**
- V4 old: wave vector mode `stop "no such direction"` (gfactor only allows k0 — now validated by `validate_semantic`)
- S1 old: bandIdx range check for wire

**From topologicalAnalysis (`main_topology.f90`):**
- S2 old: QSHE Z2 confinement check
- S3 old: BdG mode requires [bdg] check
- S4 old: BdG confinement check
- S5 old: spectral_eta > 0
- S6 old: spectral_nk >= 1
- S7 old: spectral_nE >= 1
- S8 old: sweep_model confinement match
- S9 old: conductance_method enum
- S10 old: topology mode enum

**From simulation_setup (`simulation_setup.f90`):**
- V7 old: SC + bulk warning (line ~90)
- V8 old: Electric field z(1) == 0 stop (line ~106)

After removal, each executable should call only `validate_semantic(cfg, app_name)` for config validation. The `simulation_setup_init` should have no config validation — only initialization and LAPACK/runtime error handling.

## Acceptance criteria

- [ ] All listed old checks removed from the 4 executables and `simulation_setup.f90`
- [ ] No change to error behavior — invalid configs still caught at parse time by consolidated validators
- [ ] All 34 unit tests pass
- [ ] No regression in executable behavior for valid configs

## Blocked by

- Issue 2 (structural validation in validate())
- Issue 3 (app-semantic validation in validate_semantic())
