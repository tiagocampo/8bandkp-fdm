# Issue 4: Consolidate validation into two clear layers

## Parent

PRD: Input System Architecture Cleanup (`.scratch/input-system-cleanup/PRD.md`)

## What to build

Establish two validation layers with clear ownership:

**Layer 1 — Syntactic validation (input_parser.f90):** Keep only structural/parsing checks: required section presence, TOML parse failures, enum validation (unknown confinement, unknown wire shape), required-field presence via `require_string`/`require_int`/`require_real` helpers. Remove range-bounds checks that are duplicated in `validate()` (fd_step >= 3, wire nx/ny >= FDorder+1, landau nx >= FDorder+1, landau width > 0).

**Layer 2 — Semantic validation (defs.f90):** The existing `validate()` routine remains the single source of truth for structural and state invariants. Additionally, add a new `validate_semantic(cfg, app_name)` subroutine that checks app-specific constraints dispatched on `app_name`:
- `'gfactor'`: wave_vector%nsteps must be 0 (k0 mode required)
- `'opticalProperties'`: optics%enabled must be true
- `'topologicalAnalysis'`: topo%enabled must be true, topo%mode must be non-empty
- `'bandStructure'`: no additional constraints beyond what validate() checks

Each main*.f90 executable calls `validate_semantic(cfg, '<app_name>')` after `read_config`. App-specific checks currently in main*.f90 are removed and move to `validate_semantic`.

Add pFUnit unit tests for `validate_semantic` covering: each app with valid config (passes), g-factor with non-k0 mode (errors), optics without enabled section (errors), topology without enabled section (errors), topology without mode (errors).

## Acceptance criteria

- [ ] `validate_semantic(cfg, app_name)` subroutine exists in defs.f90
- [ ] Range-bounds checks (fd_step >= 3, wire nx/ny >= FDorder+1, landau nx >= FDorder+1, landau width > 0) are removed from input_parser.f90
- [ ] Those checks remain in `validate()` in defs.f90
- [ ] App-specific checks (g-factor k0, optics enabled, topology enabled/mode) move from main*.f90 to `validate_semantic`
- [ ] Each main*.f90 calls `validate_semantic(cfg, '<app_name>')` after read_config
- [ ] pFUnit unit tests for validate_semantic cover: valid config per app, invalid g-factor k=0, invalid optics missing, invalid topology missing
- [ ] All 34 unit tests pass (existing + new)
- [ ] All ~29 regression tests pass

## Blocked by

- Issue 1 (Eliminate legacy aliases) — `validate()` references `materialN` which is renamed to `material_names` by Issue 1.
