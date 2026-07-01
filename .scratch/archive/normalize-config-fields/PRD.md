# PRD: Normalize config type fields — eliminate mode-dependent overwrites

## Problem Statement

The `simulation_config` type has three fields (`ngrid`, `num_layers`, `conf_dir`) whose meaning silently changes depending on the `confinement` mode. These are "poly-fields" — one name, multiple semantics. The parser overwrites user input without warning (e.g., `ngrid` set to `wire%ny` for wire, `num_layers` set to region count for wire). Downstream code that branches on these fields cannot distinguish QW material layers from wire regions. The wire `ngrid` value is semantically wrong (`ny` instead of `nx*ny`), masked only because wire code paths avoid `cfg%ngrid` entirely and use the grid system directly.

## Solution

Eliminate all three poly-fields in one pass:

1. **Remove `cfg%ngrid`** — replace all consumers with `cfg%grid%npoints()`. The grid system already holds the correct spatial DOF count for all modes. The user-facing TOML key `fd_step` remains unchanged.
2. **Split `cfg%num_layers`** — it means material layers only (1 for bulk/landau, `[[material]]` count for QW, 1 for wire). Wire consumers switch to `cfg%wire%num_regions` with explicit confinement-mode branches.
3. **Remove `cfg%conf_dir`** — add a pure function `conf_direction(conf)` that computes the direction character from the confinement string. Replace all stored-field checks with function calls.

Result: each concept has one name, one source of truth, no silent overwrites.

## User Stories

1. As a developer reading `simulation_config`, I want each field to have one stable meaning regardless of confinement mode, so that I don't need to mentally track mode-dependent overrides.
2. As a developer writing a new confinement mode, I want the config type to have no mode-specific overwrite traps, so that adding a mode doesn't require auditing every existing field for silent reassignment.
3. As a developer debugging wire physics, I want `cfg%grid%npoints()` to always return the correct total grid size (`nx*ny`), so that I don't encounter a stale `ny`-only value through `cfg%ngrid`.
4. As a developer reading downstream code, I want to see `cfg%wire%num_regions` explicitly when iterating over wire regions, so that the code says what it means.
5. As a developer branching on confinement direction, I want a named function `conf_direction()` instead of a stored character, so that the mapping from mode to direction is centralized and the stored field cannot go stale.
6. As a developer writing validation logic, I want `num_layers` to always mean material layers, so that validation checks don't accidentally validate region arrays against a layer count.
7. As a developer running the test suite, I want all 34+ unit tests and regression tests to pass after the refactor, so that I have confidence the normalization preserved behavior.
8. As a user writing TOML configs, I want `fd_step` to remain the TOML key name for grid points, so that my existing configs don't break.
9. As a user, I want the README TOML examples to remain valid, so that I can copy-paste them without adjustment.
10. As a developer reading the lecture documentation, I want internal field names (`ngrid`, `confDir`, `fdStep`) updated to reflect the new architecture, so that docs match the code.
11. As a developer reading the input-reference docs, I want the notes section to explain the new relationship between `fd_step` (user input) and `cfg%grid%npoints()` (computed), so that I understand the system without reading source.
12. As a developer reading the output-reference docs, I want eigenvalue count formulas to reference `grid%npoints()` instead of `ngrid`, so that the docs match the implementation.
13. As a developer reading CLAUDE.md, I want the architecture section to reflect the normalized config type, so that onboarding documentation is accurate.
14. As a developer reading the refactor plan, I want the completed normalization recorded, so that future work knows this is done.

## Implementation Decisions

### Module 1: `defs.f90` — type changes and new accessor

- Remove `ngrid` field from `simulation_config`. Remove the comment "computed grid size (was fdStep)".
- Remove `conf_dir` field from `simulation_config`.
- `num_layers` stays but its meaning is narrowed: material layers only. Wire mode sets it to 1.
- Add `public :: conf_direction` pure function:
  ```
  pure function conf_direction(conf) result(d)
    character(len=*), intent(in) :: conf
    character(len=1) :: d
    select case(trim(conf))
    case('bulk');   d = 'n'
    case('qw','wire'); d = 'z'
    case('landau'); d = 'x'
    end select
  end function
  ```
- Update `validate()` and `validate_semantic()`: replace any `cfg%ngrid` references with `cfg%grid%npoints()`. Wire validation for array sizes uses `cfg%wire%num_regions` instead of `cfg%num_layers`.

### Module 2: `input_parser.f90` — remove overwrite assignments

- Remove all `cfg%ngrid = ...` assignments (bulk line ~91, QW line ~237, wire line ~346, landau line ~403). Grid initialization already happens in `init_grid_from_config`.
- Wire mode: set `cfg%num_layers = 1` instead of `cfg%num_layers = nreg`. The `nreg` value stays in `cfg%wire%num_regions`.
- Remove `cfg%conf_dir = ...` assignments (all four mode branches). Callers now use `conf_direction(cfg%confinement)`.
- `fd_step` remains parsed from TOML and stored in `cfg%fd_step`. It is consumed only by `init_grid_from_config` which sets up the grid for QW mode.

### Module 3: Four executables — replace field references

- `main.f90`: Replace all `cfg%ngrid` with `cfg%grid%npoints()`. Replace `cfg%conf_dir` with `conf_direction(cfg%confinement)`.
- `main_gfactor.f90`: Same replacements.
- `main_optics.f90`: Same replacements.
- `main_topology.f90`: Same replacements. Wire branches already use `grid_ngrid()` — no change needed in wire paths.

### Module 4: Physics modules — replace field references

- `simulation_setup.f90`: Replace `cfg%ngrid * 8` with `cfg%grid%npoints() * 8`. For eigensolver choice, branch on `confinement` instead of `num_layers == 1`.
- `confinement_init.f90`: Replace `cfg%ngrid` fallback with `cfg%grid%npoints()`. Replace `confDir` parameter with `conf_direction(cfg%confinement)` where needed.
- `sc_loop.f90`: Replace `cfg%ngrid` with `cfg%grid%npoints()`. Wire doping iteration uses `cfg%wire%num_regions` with explicit confinement guard.
- `bdg_hamiltonian.f90`: Already uses Hamiltonian-derived size — verify no `cfg%ngrid` remains, audit only.
- `gfactor_functions.f90`, `optical_spectra.f90`, `exciton.f90`, `scattering.f90`, `outputFunctions.f90`: These receive `fdstep`/`ngrid` as parameters from callers. After the executables switch to `grid%npoints()`, the parameter names in these routines can stay as-is (they're local parameter names, not config field accesses).

### Module 5: `topological_analysis.f90` — wire midpoint fix

- Replace `layer = max(1, (cfg%num_layers + 1) / 2)` with confinement-aware logic: wire uses `max(1, (cfg%wire%num_regions + 1) / 2)`, QW uses `max(1, (cfg%num_layers + 1) / 2)`.

### Documentation updates

- `docs/reference/input-reference.md`: Update the Notes section (lines ~916-918). Replace "computed `ngrid` field" with explanation that grid size is accessible via `cfg%grid%npoints()`. Remove mention of `fd_step` being ignored for wire/landau (it's still parsed but the grid system handles mode differences internally).
- `docs/reference/output-reference.md`: Replace `ngrid` in eigenvalue count formulas (lines 31, 55) with `grid%npoints()`.
- `docs/lecture/00-quickstart.md`: Update narrative about `fdStep` warning (line 122). Update table entry (line 99) to clarify `fd_step` is used by the grid system, not stored as `ngrid`.
- `docs/lecture/02-quantum-well.md`: Replace `confDir = 'z'` reference (line 1046) with `conf_direction()` explanation.
- `docs/lecture/05-gfactor.md`: Replace `fdStep`/`ngrid` prose references (lines 320, 392, 471) with `grid%npoints()` terminology.
- `docs/lecture/08-quantum-wire.md`: Update backward-compatibility note (line 325) to state that `fdStep`/`confDir` fields no longer exist. Update `confDir` reference (line 633).
- `docs/lecture/12-extending-the-code.md`: Replace `Ngrid`/`fdStep` references (line 208) with `grid%npoints()`.
- `CLAUDE.md`: Update architecture section — remove `ngrid` from dual-mode description (line 167), update `fd_step` description (line 182), update `grid_ngrid` exception note (line 197). Add `conf_direction` to module dependency concepts.
- `README.md`: TOML examples with `fd_step` are unchanged (still valid TOML key). Prose around examples may need `ngrid` → `grid%npoints()` if any exists.

### What does NOT change

- The TOML key `fd_step` remains. All 88 TOML configs in `tests/regression/configs/` are unchanged.
- All lecture scripts (`scripts/lecture_*.py`) that generate TOML config strings are unchanged — `fd_step = <value>` is still valid TOML.
- All integration test config strings are unchanged.
- `cfg%fd_step` remains in the type — it's the raw user input, used by grid initialization for QW mode.

## Testing Decisions

- **Good test**: verify external behavior (eigenvalues, grid sizes, validation errors) not internal field names. After the refactor, the same TOML input must produce the same physics output.
- **Unit tests** (pFUnit in `tests/unit/`):
  - Add unit test for `conf_direction()` — verify all four modes return the correct character.
  - Add unit test for `validate_semantic` with wire configs — verify it validates against `wire%num_regions`, not `num_layers`.
  - Update existing `test_parameters.pf` if it references `cfg%ngrid` or `cfg%conf_dir`.
- **Regression tests** (`tests/regression/`): All 88 configs must produce identical output after the refactor. No config changes needed — `fd_step` TOML key is unchanged.
- **Integration tests** (`tests/integration/`): Lecture scripts and convergence tests must pass without modification.
- **Prior art**: The `validate_semantic` unit tests added in commit `1e8ce9f` follow the pattern — create a minimal config, call the validator, check error conditions.

## Out of Scope

- Removing `cfg%fd_step` from the type — it remains as the user-facing input for QW grid sizing.
- Changing any TOML key names — `fd_step` stays as-is in all configs.
- Updating legacy `.cfg` format files — those are pre-TOML and not maintained.
- Updating archived plan documents in `docs/plans/archive/` — historical reference only.
- Changing physics module parameter names (e.g., `fdstep` argument in `gfactor_functions.f90`) — these are local parameter names passed by value, not config field accesses.
- Wire mode grid initialization refactoring — already correct via `init_grid_from_config`.
- Adding new confinement modes.
- Any changes to the k.p block table, Hamiltonian construction, or FD stencils.

## Further Notes

- The wire `cfg%ngrid = ny` bug (set to `wire%ny` instead of `nx*ny`) is effectively harmless because no wire code path uses `cfg%ngrid`. But eliminating the field entirely removes the trap for future code.
- The `conf_direction` function introduces a central mapping that must be updated if a new confinement mode is added. This is intentional — a single source of truth is better than scattered character assignments.
- `simulation_setup.f90`'s `num_layers == 1` heuristic for choosing between `zheev` and `zheevd` needs a new criterion after the split. The natural replacement is `cfg%confinement == 'bulk'` — bulk is the only mode where the 8x8 dense path uses `zheev`.
- The 88 TOML configs and all lecture scripts remain unchanged because `fd_step` as a TOML key is not being removed — only the internal `ngrid` computed field is being eliminated in favor of the grid system accessor.
