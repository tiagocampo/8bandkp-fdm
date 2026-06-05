# PRD: TOML Input Parser Refactor

**Labels:** enhancement, ready-for-agent
**Branch:** feat/input-parser-toml (from main)

## Problem Statement

The input parser (`input_parser.f90`) is a 1,369-line monolithic subroutine with 36 `backspace` calls for peek-ahead control flow, fragile positional ordering constraints, and mode-dependent field overwrites (e.g., `cfg%fdStep` gets silently overwritten for wire and Landau modes). The parser uses a custom `label: value` format that requires strict line ordering, making it error-prone to add new fields and nearly impossible to write comprehensive parser tests without fragile file-position assumptions.

Adding a new optional block (e.g., exciton, scattering, strain) currently requires inserting peek/backspace logic at exactly the right position in the sequential read stream. The topology block alone is 220 lines of deeply nested optional sub-fields, each using backspace for conditional parsing.

## Solution

Replace the entire custom parser with a TOML-based parser using the `toml-f` library (pure Fortran 2008, MIT/Apache dual license, used by fpm itself). All configuration files will be converted from the custom `label: value` format to standard TOML. The `simulation_config` derived type will be restructured into sub-types mirroring the TOML section structure. A computed `ngrid` field replaces the mode-dependent `fdStep` overwrites.

This is a **full migration** — no backward compatibility, no dual-format hacks, no compatibility shims. The old `label: value` format is dead. Every config file, every consumer, every test, every lecture document, and every reference doc will be updated to the new TOML format.

## User Stories

1. As a physicist running simulations, I want to write config files in a standard, well-documented format (TOML), so that I don't need to memorize line ordering rules
2. As a physicist running simulations, I want config sections to be order-independent, so that I can organize parameters in a way that makes sense to me
3. As a physicist running simulations, I want to see clear error messages when a required field is missing, so that I can fix config errors quickly
4. As a physicist running simulations, I want `confinement = "qw"` instead of `confinement = 1`, so that I don't need to memorize integer codes
5. As a physicist running simulations, I want material layers defined as `[[material]]` tables with named fields, so that I can easily see what each value means
6. As a physicist running simulations, I want optional physics blocks (topology, optics, SC, etc.) to be enabled by their presence in the config, so that I don't need to add `T/F` enable flags
7. As a physicist using wire geometry, I want wire parameters in a `[wire]` section with named sub-fields, so that the config is self-documenting
8. As a physicist using Landau levels, I want Landau parameters in a `[landau]` section, so that they are clearly separated from QW parameters
9. As a physicist using external fields, I want the electric field and magnetic field to be in separate `[external_field]` and `[b_field]` sections, so that they don't collide with each other
10. As a physicist using strain, I want strain parameters in a `[strain]` section, so that I can easily find and modify them
11. As a physicist using self-consistent Schrödinger-Poisson, I want SC parameters in an `[sc]` section with human-readable keys, so that the config is easier to understand
12. As a physicist using topology analysis, I want gap sweep, conductance, and spectral parameters in nested `[topology.*]` sections, so that complex options are organized hierarchically
13. As a physicist using doping profiles, I want doping defined as `[[doping]]` tables, so that per-layer uniform and delta doping are handled uniformly
14. As a physicist using g-factor calculations, I want `which_band` and `band_idx` as named keys, so that I know what they control
15. As a physicist using optics, I want optical parameters in an `[optics]` section, so that they are grouped together
16. As a physicist using FEAST eigensolver, I want FEAST tuning parameters in a `[feast]` section, so that they are clearly separated from simulation parameters
17. As a developer maintaining the parser, I want no `backspace` calls in the codebase, so that the parser is linear and easy to reason about
18. As a developer adding new config fields, I want to add a single `get_value` call for optional fields, so that new features don't require peek/backspace logic
19. As a developer writing parser tests, I want to write TOML content to a temp file and parse it, so that tests are order-independent and easy to write
20. As a developer reviewing the parser, I want each optional block parsed by its own sub-parser function, so that the code is modular and each block's logic is isolated
21. As a developer debugging config issues, I want the parser to report the TOML section and key name in error messages, so that I can find the problem immediately
22. As a developer working on simulation_config, I want sub-types matching TOML sections, so that the type structure mirrors the config structure
23. As a developer, I want a computed `ngrid` field instead of mode-dependent `fdStep` overwrites, so that the original user input is preserved alongside the computed grid size
24. As a developer, I want snake_case TOML keys, so that the config follows TOML community conventions
25. As a developer, I want no per-field `print *` echo of parsed values, so that stdout is clean and reserved for simulation output
26. As a developer running tests, I want config files named `input.toml`, so that the file extension signals the format
27. As a developer running CI, I want all 88 regression configs, 32 test scripts, 15 lecture scripts, and 12 validation configs updated atomically, so that the codebase is consistent
28. As a developer, I want toml-f integrated as a git submodule, so that builds work offline and are reproducible
29. As a reader of lecture documentation, I want all config examples shown in TOML format, so that I can copy-paste them directly into `input.toml`
30. As a reader of the input reference, I want a complete TOML schema reference with all sections, keys, types, defaults, and examples, so that I can write configs without reading source code
31. As a reader of lecture docs, I want inline config snippets updated from the old format to TOML, so that the docs match what I actually need to write
32. As a reader of the output reference, I want references to `input.cfg` updated to `input.toml`, so that I'm not confused about which file to edit
33. As a reader of the README, I want quickstart examples in TOML format, so that my first experience with the project uses the current format
34. As a reader of brainstorm/ideation docs, I want to know that old `.cfg` references are historical, so that I don't try to use outdated config formats

## Implementation Decisions

### TOML library: toml-f via git submodule

- `toml-f` v0.4.2 added as a git submodule in `subprojects/toml-f/`
- Pure Fortran 2008, compatible with `-std=f2018` enforcement
- CMake integration via `add_subdirectory(subprojects/toml-f)`, linking against `toml-f-lib`
- API: `toml_parse(table, io, error)` to load file, `get_value(table, key, val, default)` to extract values

### Config filename: `input.toml`

- Hard rename from `input.cfg` to `input.toml`
- All four executables hard-code the new filename
- No `--config` CLI argument (out of scope)

### No backward compatibility, no hacks

- The old `label: value` format is completely removed
- No dual-format support, no compatibility shims, no parser that handles both formats
- Every config file in the repo is converted to TOML
- The conversion is atomic — the feature branch will not mix old and new formats

### TOML key naming: snake_case

- All TOML keys use `snake_case` (e.g., `wave_vector_max`, `num_cb`, `fd_step`)
- Internal Fortran field names follow the same convention
- No camelCase or abbreviated keys in TOML files

### Confinement: string enum

- `confinement = "bulk"` instead of `confinement = 0`
- Four values: `"bulk"`, `"qw"`, `"wire"`, `"landau"`
- All consumer `select case` and `if` comparisons updated to use strings directly — no integer compatibility layer

### Section presence = enabled

- Optional physics blocks (topology, optics, SC, strain, etc.) are enabled by their presence in the TOML file
- No `enabled = true` keys inside sections
- If `[topology]` exists in the file, topology is enabled

### TOML schema structure

Every multi-field group is a TOML section. Only true singletons are top-level keys:

- Top-level: `confinement`, `FDorder`, `fd_step`
- `[wave_vector]`: `mode`, `max`, `step`
- `[bands]`: `num_cb`, `num_vb`
- `[[material]]`: `name`, `z_min`, `z_max` (QW/bulk)
- `[wire]`: `nx`, `ny`, `dx`, `dy` + `[wire.geometry]` + `[[region]]`
- `[landau]`: `nx`, `width`, `sweep`, `material`
- `[external_field]`: `type`, `value`
- `[b_field]`: `components` (3-element array), `g_factor`
- `[strain]`: `substrate`
- `[sc]`: all SC parameters
- `[[doping]]`: `layer`, `ND`, `NA` (uniform) or `layer`, `type`, `NS`, `fwhm`, `pos` (delta)
- `[topology]`: all topology parameters, with nested `[topology.gap_sweep]`, `[topology.conductance]`, `[topology.spectral]`
- `[bdg]`: all BdG parameters
- `[optics]`: all optical spectra parameters
- `[exciton]`: exciton solver parameters
- `[scattering]`: phonon scattering parameters
- `[feast]`: `emin`, `emax`, `m0`

### simulation_config type restructure

New sub-types mirroring TOML sections:

```
wave_vector_config: mode, max, step
bands_config: num_cb, num_vb
external_field_config: type, value
b_field_config: components(3), g_factor
wire_config: nx, ny, dx, dy, geom, regions
landau_config: nx, width, sweep, material
feast_config: emin, emax, m0
```

The `simulation_config` type gains a computed `ngrid` field that replaces mode-dependent `fdStep` overwrites:
- bulk: `ngrid = 1`
- QW: `ngrid = fd_step`
- wire: `ngrid = wire%ny`
- Landau: `ngrid = landau%nx`

The `confinement` field changes from `integer` to `character(len=8)` — all `select case` and `if` comparisons in `simulation_setup`, physics modules, and apps are updated to compare against string values directly.

### Parser architecture: sub-parser functions

The single `read_config` subroutine is replaced by a coordinator that calls section-specific sub-parsers:
- `parse_wave_vector(table, cfg)`
- `parse_bands(table, cfg)`
- `parse_materials(table, cfg)` — bulk/QW
- `parse_wire(table, cfg)` — wire geometry + regions
- `parse_landau(table, cfg)` — Landau grid
- `parse_external_field(table, cfg)`
- `parse_b_field(table, cfg)`
- `parse_sc(table, cfg)` + `parse_doping(table, cfg)`
- `parse_topology(table, cfg)` + nested gap_sweep/conductance/spectral
- `parse_bdg(table, cfg)`
- `parse_optics(table, cfg)`
- `parse_exciton(table, cfg)`
- `parse_scattering(table, cfg)`
- `parse_feast(table, cfg)`
- `parse_strain(table, cfg)`

### Error handling: fail-fast + defaults

- `require_value(table, key, val, section_name)`: prints `"Error: [section] requires key 'name'"` and stops
- `get_value(table, key, val, default)`: silently uses default if key is missing
- No `print *` echo of parsed values — stdout is reserved for simulation output

### Consumer field migration mapping

Every reference to old field names is updated. No compatibility layer.

Key renames (252 references across 12 source files):
- `cfg%fdStep` → `cfg%ngrid` (computed field)
- `cfg%confinement == 0/1/2/3` → `cfg%confinement == "bulk"/"qw"/"wire"/"landau"`
- `cfg%waveVector` → `cfg%wave_vector%mode`
- `cfg%waveVectorMax` → `cfg%wave_vector%max`
- `cfg%waveVectorStep` → `cfg%wave_vector%step`
- `cfg%numcb` → `cfg%bands%num_cb`
- `cfg%numvb` → `cfg%bands%num_vb`
- `cfg%materialN` → `cfg%material_names`
- `cfg%startPos` → `cfg%z_min`
- `cfg%endPos` → `cfg%z_max`
- `cfg%intStartPos` → `cfg%int_start_pos`
- `cfg%intEndPos` → `cfg%int_end_pos`
- `cfg%confDir` → `cfg%conf_dir`
- `cfg%ExternalField`/`cfg%EFtype` → `cfg%external_field%type`
- `cfg%Evalue` → `cfg%external_field%value`
- `cfg%wire_nx/ny/dx/dy` → `cfg%wire%nx/ny/dx/dy`
- `cfg%landau_nx/width/sweep` → `cfg%landau%nx/width/sweep`
- `cfg%feast_emin/emax/m0` → `cfg%feast%emin/emax/m0`
- `cfg%whichBand` → `cfg%which_band`
- `cfg%bandIdx` → `cfg%band_idx`

### Wire polygon vertices in TOML

Polygon vertex data uses TOML arrays of arrays:
```toml
[wire.geometry]
shape = "polygon"
vertices = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]
```

### Doping uniform representation

Both uniform and delta doping use `[[dopping]]` tables. The `type` key defaults to `"uniform"`:
```toml
[[doping]]
layer = 1
ND = 1e18
NA = 0.0

[[doping]]
layer = 2
type = "delta"
NS = 1e12
fwhm = 5.0
pos = 50.0
```

### Documentation overhaul

All documentation is updated from old `label: value` format to TOML. No doc references the old format as current.

**Scope (~90 references in current docs):**

- `docs/reference/input-reference.md` — complete rewrite as TOML schema reference with all sections, keys, types, defaults, and TOML examples
- `docs/reference/output-reference.md` — update `input.cfg` → `input.toml` references
- `docs/reference/benchmarks.md` — update config filename references
- `README.md` — rewrite quickstart section with TOML config examples
- `docs/lecture/*.md` (15 files) — update all inline config snippets from old format to TOML, update `input.cfg` → `input.toml`, update `cat config.cfg > input.cfg` → `cp config.toml input.toml`
- `CLAUDE.md` — update "Running" section, "Input file" section, and "Architecture" section with TOML format details

Note: old plans (`docs/plans/`), brainstorm archives (`docs/brainstorms/archive/`), and ideation docs (`docs/ideation/`) are NOT updated — they are historical records. Only current active documentation is in scope.

## Testing Decisions

### What makes a good test

Parser tests should test **external behavior** (TOML file in → `simulation_config` field values out), not implementation details (no assertions about internal toml-f API calls). Tests write TOML content to temp files, parse them, and assert field values.

### Module 1: TOML Config Parser (highest priority)

New comprehensive parser tests in pFUnit format:
- **Bulk config**: minimal valid TOML, verify defaults, verify required field errors
- **QW config**: material layers, fd_step, grid computation (z, delta, dz, int_start_pos)
- **Wire config**: wire grid, each shape variant (circle, rectangle, hexagon, polygon), regions
- **Landau config**: grid params, sweep modes, material
- **Each optional block**: topology (all sub-sections), SC + doping, optics, exciton, scattering, FEAST, strain, BdG, b_field, external_field
- **Validation errors**: missing required fields, invalid confinement, fd_step too small for FDorder, invalid FDorder
- **Defaults**: verify that omitting optional fields produces correct default values
- **Mode-specific validation**: bulk must have 1 material, QW fd_step >= 3, wire nx/ny >= 3, Landau nx >= 3

Prior art: `test_topology_parser.pf` and `test_bdg_config.pf` — these write config to disk and parse. The new tests follow the same pattern but with TOML content.

### Module 2: simulation_config type defaults

- Verify all new sub-types have correct default values
- Verify `validate()` catches: invalid confinement string, fd_step < 3 for QW, missing materials, conflicting options

### Module 3: Config File Converter (Python script)

- Convert known old-format configs and verify TOML output matches expected schema
- Round-trip test: convert → parse → verify same numerical values as old parser would produce
- Test all 4 geometry modes
- Test each optional block

### Consumer tests (no new tests needed)

Existing unit tests (`test_hamiltonian*.pf`, `test_sc_loop.pf`, `test_simulation_setup.pf`, etc.) construct `simulation_config` objects directly — field renames will be caught at compile time. No new tests needed for mechanical renaming.

## Out of Scope

- `--config` CLI argument for custom filenames
- Backward compatibility with old `input.cfg` format (no dual-format support)
- TOML generation/serialization (parser is read-only)
- Schema validation beyond what the parser already checks (no external JSON Schema or TOML schema validator)
- Changes to the physics modules' algorithms or interfaces beyond field renaming
- Changes to the output format or file structure
- Any modifications to material parameters in `parameters.f90`
- Migration of the `fpm.toml` experimental build system
- GUI or interactive config editors

## Further Notes

### Relationship to ADR 0001

ADR 0001 (`simulation-setup-fat-type`) established the `read_config` / `simulation_setup_init` split. This PRD deepens that separation by making `read_config` a clean TOML parser with no physics logic. The `simulation_setup` fat type and its `select case` dispatch remain — only the comparison operand changes from integer to string.

### Migration script

A Python conversion script (`scripts/convert_cfg_to_toml.py`) will automate the bulk of the 88 config file conversions. Manual spot-checking is still required for configs with unusual feature combinations.

### Build dependency chain

The implementation must follow a strict order:
1. Add toml-f submodule and verify build (Phase 1)
2. Restructure types in defs.f90 (Phase 2) — breaks all consumers until Phase 4
3. Rewrite parser (Phase 3)
4. Update consumers (Phase 4) — restores compilation
5. Convert configs (Phase 5)
6. Update test scripts (Phase 6)
7. Documentation overhaul (Phase 7) — 121 doc references across 20+ files
8. Full verification (Phase 8)

Phases 2-4 should be done as a single atomic commit (or feature branch) to avoid leaving the codebase in a broken state.

### Estimated scope

- ~175 files touched
- 252 consumer field references to rename (across 12 source files)
- 12 unit test files with direct `simulation_config` construction
- 2 parser test files to rewrite completely (`test_topology_parser.pf`, `test_bdg_config.pf`)
- 88 regression configs to convert
- 32 integration test shell scripts to update
- 15 lecture scripts to update
- 12 validation configs to convert
- ~90 documentation references in current docs (15 lecture docs, input reference, output reference, benchmarks, README, CLAUDE.md)
