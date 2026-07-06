**Status**: COMPLETE (2026-07-05)

# PRD: Consolidate Config Validation (C4)

## Problem Statement

Config validation is split across four layers — the TOML parser, two validators in `defs.f90`, the four executables, and `simulation_setup.f90`. Invalid configurations can pass the parser and structural validator, then crash deep in physics code with unhelpful errors. Worse, some invalid configs are silently corrected: a user who sets `num_cb = 20` for bulk confinement gets results for `num_cb = 2` with only a warning buried in stdout. There is no single place to find or modify validation rules.

## Solution

Consolidate all config-level validation — checks catchable from `simulation_config` alone at parse time — into the two existing validators in `defs.f90`. The parser retains its existing key-existence and type-correctness checks. Runtime checks (LAPACK results, FEAST convergence, file I/O) stay in the executables where they belong. After consolidation, the old duplicated checks in executables and `simulation_setup.f90` are removed.

Architecture decision record: `docs/adr/0002-consolidate-validation.md`.

## User Stories

### Structural validator (`validate()`) — config-agnostic checks

1. As a user running bulk bandStructure, I want the solver to reject `evnum > 8` immediately with a clear error, so that I don't get results for a different band count than I requested.
2. As a user running QW bandStructure, I want the solver to reject `num_cb` exceeding the maximum available conduction bands, so that I know my configuration is physically impossible.
3. As a user running QW bandStructure, I want the solver to reject `num_vb` exceeding the maximum available valence bands, so that I know my configuration is physically impossible.
4. As a user, I want the solver to reject unrecognized `wave_vector` modes immediately, so that I see a clear error instead of a cryptic crash later.
5. As a user running B-sweep calculations, I want the solver to reject zero or negative B-sweep step size, so that I don't enter an infinite loop or get nonsensical results.
6. As a user specifying material layers, I want the solver to reject `z_min >= z_max`, so that I catch inverted or zero-width layers before the solver silently produces wrong results.
7. As a user enabling self-consistent mode, I want the solver to reject SC with bulk confinement, so that I'm told the SC loop requires a confined geometry.
8. As a user enabling an external electric field, I want the solver to reject `z(1) == 0`, so that the field calculation doesn't silently produce zero field.

### App-semantic validator (`validate_semantic()`) — app-specific checks

9. As a user running gfactorCalculation with a wire, I want the solver to reject out-of-range `bandIdx`, so that I get a clear error instead of an array bounds crash.
10. As a user running topologicalAnalysis with QSHE Z2 mode, I want the solver to reject non-QW/non-wire confinement, so that I'm told this topology mode requires a confined geometry.
11. As a user running topologicalAnalysis with BdG mode, I want the solver to reject missing `[bdg]` section, so that I'm told to add the required configuration block.
12. As a user running topologicalAnalysis with BdG mode, I want the solver to reject non-QW/non-wire confinement, so that I'm told BdG requires a confined geometry.
13. As a user running topologicalAnalysis spectral mode, I want the solver to reject non-positive `spectral_eta`, so that I get a clear error about the broadening parameter.
14. As a user running topologicalAnalysis spectral mode, I want the solver to reject `spectral_nk < 1`, so that I get a clear error about the k-point count.
15. As a user running topologicalAnalysis spectral mode, I want the solver to reject `spectral_nE < 1`, so that I get a clear error about the energy count.
16. As a user running topologicalAnalysis with gap sweeps, I want the solver to reject mismatched sweep_model and confinement, so that I'm told which confinement the model requires.
17. As a user running topologicalAnalysis, I want the solver to reject unrecognized `conductance_method` values, so that I see valid options in the error message.
18. As a user running topologicalAnalysis, I want the solver to reject unrecognized `topology mode` values, so that I see valid options in the error message.

### Parser guard

19. As a user configuring a QW with `fd_step = 1`, I want the parser to reject this before the division-by-zero crash, so that I get a clear error about the minimum grid size.

### Developer stories

20. As a developer adding a new validation rule, I want to add it in exactly one place (`validate()` or `validate_semantic()`), so that I don't have to remember which executables need updating.
21. As a developer debugging a validation failure, I want error messages to include the actual value and what was expected, so that I can fix the config without reading source code.
22. As a developer testing new validation rules, I want to write a pFUnit test that constructs a valid config and confirms the validator accepts it, so that I can verify the acceptance criterion.

## Implementation Decisions

### Two validators, two responsibilities

The existing `validate()` (type-bound on `simulation_config`) and `validate_semantic(cfg, app_name)` (standalone subroutine) are expanded in place. No new subroutines, no new files, no new types.

- `validate()` absorbs checks that require no app context: evnum caps, wave vector enum, band count ranges, B-sweep step, material layer ordering, SC+bulk mismatch, electric field+z(1).
- `validate_semantic(cfg, app_name)` absorbs checks that require knowing which executable is running: gfactor bandIdx range, topology mode/confnement/parameter validation.

This preserves the structural/semantic boundary established by the existing code.

### Error mechanism

All new checks use `error stop` with contextual error messages. Each message includes the actual value and what was expected:

```
Error: evnum (=20) exceeds 8 for bulk confinement (8x8 Hamiltonian)
Error: wave_vector mode 'kz' invalid for gfactor (requires k0)
Error: SC loop requires confinement='qw', got 'bulk'
```

Existing `stop 1` calls in the parser are not refactored in this pass.

### Hard errors, no silent corrections

All moved checks are hard errors. The current behavior of silently capping `evnum` to 8 for bulk, or capping band counts for QW, is eliminated. Users get a clear error telling them exactly what's wrong, not results for a different configuration than they requested.

### Old checks removed

After consolidation, the original checks in the executables and `simulation_setup.f90` are removed. No duplication — one place to look, one place to maintain.

### Parser guard for fd_step division-by-zero

The division `cfg%delta = cfg%totalSize / real(cfg%fd_step - 1, dp)` in the QW material parser runs before `validate()` is called. A one-line guard is added in the parser immediately before this division to fail fast on `fd_step < 2`. The existing `validate()` check (`fd_step >= 3` for QW) remains as belt-and-suspenders.

### Moderate scope

Only crash-level gaps are addressed: the fd_step division-by-zero and `z_min >= z_max` for material layers. Softer validation gaps (negative temperature, negative grid spacings, unrecognized fermi_mode fallback, ldos_E_range ordering, wire region inner/outer radius) are deferred to a future pass.

### Implementation slices (all AFK)

| Issue | Slice | Blocked by |
|-------|-------|------------|
| #20 | Parser guard: reject fd_step < 2 (P1) | None |
| #21 | Structural validation in validate() (V1–V8) | None |
| #22 | App-semantic validation in validate_semantic() (S1–S10) | None |
| #23 | Remove old duplicated checks from executables | #21, #22 |
| #24 | pFUnit passing-case tests for new validators | #21, #22 |
| #25 | Mechanical documentation updates | #23, #24 |

### Modules modified

| Module | Change | Slice |
|--------|--------|-------|
| `input_parser.f90` — `parse_materials_qw` | Add P1 guard before fd_step division | #20 |
| `defs.f90` — `simulation_config_validate` | Add V1–V8 checks | #21 |
| `defs.f90` — `validate_semantic` | Add S1–S10 checks | #22 |
| `main.f90` | Remove V1–V5 old checks | #23 |
| `main_gfactor.f90` | Remove V4, S1 old checks | #23 |
| `main_topology.f90` | Remove S2–S10 old checks | #23 |
| `simulation_setup.f90` | Remove V7–V8 old checks | #23 |
| `test_defs.pf` | Add passing-case tests for V1–V8 | #24 |
| `test_parameters.pf` | Add passing-case tests for S1–S10 | #24 |
| `docs/adr/0002-consolidate-validation.md` | Status: Proposed → Accepted | #25 |
| `CLAUDE.md` | Add validation bullet, remove stale language | #25 |
| `docs/reference/input-reference.md` | Remove stale silent-correction language if present | #25 |

## Testing Decisions

### Test philosophy

Tests verify external behavior (validator accepts valid configs), not implementation details (which specific `if` branch runs). Following the existing pattern established in commit `1e8ce9f`.

### Testability constraint

pFUnit 4.x cannot catch `error stop`. Tests cover passing cases only: construct a valid config, call the validator, assert it does not crash. Failure cases are verified manually during development.

### Modules tested

- `validate()` new checks — tested in `test_defs.pf` (same file as existing validate tests)
- `validate_semantic()` new checks — tested in `test_parameters.pf` (same file as existing validate_semantic tests)

### Prior art

- `test_parameters.pf:228–269` — five existing `validate_semantic` passing-case tests
- `test_defs.pf:174–210` — existing `validate()` passing-case test for wire multi-region

Each new check gets one test: a minimal valid config that exercises the acceptance path for that specific rule.

## Out of Scope

- **Refactoring existing `stop 1` to `error stop`** in `input_parser.f90` or `simulation_setup.f90` — separate cleanup pass
- **Adding missing range checks** for negative temperature, negative grid spacings, wire region inner/outer, ldos_E_range ordering — deferred
- **Silent enum fallback for fermi_mode** — deferred
- **Runtime error handling** (LAPACK info checks, FEAST convergence, file I/O, Simpson parity) — stays in executables
- **Creating per-domain sub-validators** — over-engineering for ~20 checks
- **Merging `validate()` and `validate_semantic()`** into one subroutine — would violate the parser↔app boundary

## Further Notes

- This work resolves candidate C4 from the 2026-05-27 architecture review. Candidates C1–C3 and C5 are already resolved by the TOML migration and normalize-config work.
- The ADR is at `docs/adr/0002-consolidate-validation.md`.
- Implementation follows a 4-phase order: add checks → write tests → remove old checks → docs. Each phase should compile and pass all tests independently.
- The `electric field + z(1) == 0` check (V8) needs a guard for `allocated(cfg%z)` — bulk mode may not allocate the z array.
