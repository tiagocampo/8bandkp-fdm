# ADR 0002: Consolidate Config Validation

## Status

Accepted

## Context

Config validation is split across four layers:

1. **Parser** (`input_parser.f90`, ~50+ checks): Required keys, TOML type correctness, enum values. Uses `stop 1`.
2. **Structural** (`defs.f90` `validate()`, 19 checks): Grid size, band counts, material allocations, mode-specific grid constraints. Uses `error stop`.
3. **App-semantic** (`defs.f90` `validate_semantic()`, 5 checks): gfactor↔k0, optics↔enabled, topology↔enabled+mode. Uses `error stop`.
4. **Runtime** (4 executables + `simulation_setup.f90`, 48+ checks): evnum caps, wave vector mode, LAPACK results, FEAST convergence, file I/O. Mixed `stop`/`stop 1`/warnings.

An invalid config can pass the parser and structural validator, then crash deep in physics code with an unhelpful error. For example, `evnum = 20` with `confinement = 'bulk'` passes the parser, but `main.f90` silently caps it to 8 with a warning — the user gets results for a different configuration than they requested.

The architecture review (2026-05-27) identified this as candidate C4. Candidates C1–C3 and C5 are resolved by the TOML migration and normalize-config work.

## Decision

Consolidate config-level validation into the existing two validators. The dividing line is **catchable from config alone at parse time** — if a check only depends on `simulation_config` fields (including computed fields like `grid%npoints()`), it moves. If it depends on runtime state (LAPACK results, FEAST convergence, file I/O), it stays.

### Checks added to `validate()` (config-agnostic, no app knowledge)

| ID | Check | Source |
|----|-------|--------|
| V1 | `evnum > 8` for bulk | main.f90:415 |
| V2 | `num_cb > max` for QW | main.f90:425 |
| V3 | `num_vb > max` for QW | main.f90:429 |
| V4 | Wave vector mode is valid enum | main.f90:116 |
| V5 | B-sweep step > 0 when present | main.f90:759 |
| V6 | `z_min(i) < z_max(i)` for all materials | New |
| V7 | SC + bulk confinement mismatch | simulation_setup.f90:90 |
| V8 | Electric field + `z(1) == 0` | simulation_setup.f90:106 |

### Checks added to `validate_semantic()` (app-specific)

| ID | App | Check | Source |
|----|-----|-------|--------|
| S1 | gfactor | bandIdx range for wire | main_gfactor.f90:135 |
| S2 | topologicalAnalysis | confinement must be QW/wire for QSHE Z2 | main_topology.f90:166 |
| S3 | topologicalAnalysis | BdG mode requires [bdg] section | main_topology.f90:181 |
| S4 | topologicalAnalysis | BdG confinement must be QW/wire | main_topology.f90:191 |
| S5 | topologicalAnalysis | spectral_eta > 0 | main_topology.f90:925 |
| S6 | topologicalAnalysis | spectral_nk >= 1 | main_topology.f90:933 |
| S7 | topologicalAnalysis | spectral_nE >= 1 | main_topology.f90:937 |
| S8 | topologicalAnalysis | sweep_model matches confinement | main_topology.f90:1137 |
| S9 | topologicalAnalysis | conductance_method is valid enum | main_topology.f90:1099 |
| S10 | topologicalAnalysis | topology mode is valid enum | main_topology.f90:227 |

### Parser-level guard

| ID | Check | Location |
|----|-------|----------|
| P1 | `fd_step >= 2` before division | input_parser.f90, before line 233 |

### Key design choices

1. **All moved checks are hard errors** (`error stop`). No silent corrections. If a user requests `num_cb = 20` for bulk, they get a clear error, not silently capped results.

2. **Error messages include context**: the actual value and what was expected. E.g., `evnum (=20) exceeds 8 for bulk confinement (8x8 Hamiltonian)`.

3. **Old checks removed from executables** after consolidation. No duplication — one place to look.

4. **Moderate scope**: only crash-level gaps addressed (fd_step division-by-zero, z_min >= z_max). Softer gaps (negative temperature, unrecognized fermi_mode, negative grid spacings) deferred.

5. **Passing-case pFUnit tests only**, following the existing pattern. `error stop` cannot be caught by pFUnit 4.x, so failure cases are tested via shell script (`test_validate_rejects_bad_configs.sh`) that feeds bad TOML configs to the executable and checks for non-zero exit codes and expected error patterns.

6. **Implementation order**: add checks → write tests → remove old checks → docs. Each phase compiles and passes tests independently.

## Alternatives Considered

### Per-domain sub-validators (`validate_wave_vector`, `validate_bands`, etc.)

Create 5-6 new subroutines, each validating one config section.

**Rejected**: Over-engineering for ~20 checks. Adds call-site complexity with no benefit — the caller still needs to call all of them in sequence. Two validators (structural + app-semantic) are simpler and already established.

### Single unified `validate_all(cfg, app_name)`

Merge `validate()` and `validate_semantic()` into one subroutine.

**Rejected**: Loses the structural/semantic separation. `validate()` runs at parse time (no app context); `validate_semantic()` runs in each executable (has app context). Merging requires threading `app_name` through the parser, violating the parser↔app boundary established in ADR 0001.

### Soft corrections (warn and continue)

Keep the existing warn-and-cap behavior for evnum, band counts, etc.

**Rejected**: Silent corrections are a footgun. Users get results for a different configuration than requested, and may not notice the warning in stdout. Hard errors with clear messages are more honest.

### Aggressive scope (add all validation gaps)

Add negative temperature, negative grid spacing, unrecognized fermi_mode, ldos_E_range ordering, wire region inner/outer, etc.

**Deferred**: These are "wrong physics" risks, not crash risks. Some have legitimate edge cases (temperature = 0 is physically meaningful in some contexts; negative strain is real). Deserve their own pass with more domain-specific thought.

## Consequences

### Positive

- **One place to look**: all config validation in two subroutines in `defs.f90`.
- **Fail early**: invalid configs caught at parse time, not deep in physics code.
- **Clear error messages**: users see what's wrong and what's expected.
- **No silent corrections**: users never get results for a different config than requested.

### Negative

- **Harder to bypass**: code that constructs configs programmatically must call `cfg%validate()` explicitly. The old executable checks were a safety net that's now removed.
- **Moderate scope leaves gaps**: negative temperature, fermi_mode fallback, and other "wrong physics" issues are not addressed. These could produce wrong results silently.
- **Test coverage limited**: ~47 `error stop` branches cannot be unit-tested with pFUnit 4.x. Failure paths are covered by `test_validate_rejects_bad_configs.sh` (8 test cases covering the most impactful rejection checks). Remaining branches are verified manually.
