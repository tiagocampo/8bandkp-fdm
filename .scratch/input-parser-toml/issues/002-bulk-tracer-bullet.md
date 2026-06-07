# Issue 002: Bulk mode tracer bullet — type restructure + TOML parser + consumer migration

**Type:** AFK
**Blocked by:** Issue 001 (toml-f build integration)
**User stories:** US 1-4, 6, 10, 16-18, 20-28

## Parent

PRD: `.scratch/input-parser-toml/PRD.md`

## What to build

This is the **keystone tracer bullet** that cuts through every integration layer. It restructures the `simulation_config` derived type, rewrites the parser for **bulk mode only** (simplest geometry — single material, no grid computation), renames ALL consumer fields across the entire codebase, and converts bulk configs to TOML. After this issue, the codebase compiles and all tests pass with the new TOML-based bulk mode.

This is the biggest single slice because phases 2-4 of the PRD must be atomic — the codebase won't compile between type restructure and consumer updates.

### Type restructure (defs.f90)

Add new sub-types mirroring TOML sections:
- `wave_vector_config`: mode, max, step
- `bands_config`: num_cb, num_vb
- `external_field_config`: type, value
- `b_field_config`: components(3), g_factor
- `feast_config`: emin, emax, m0

Update `simulation_config`:
- `confinement` changes from `integer` to `character(len=8)` with values `"bulk"`, `"qw"`, `"wire"`, `"landau"`
- `fdStep` renamed to `fd_step` (user input) + new computed `ngrid` field
- `waveVector/waveVectorMax/waveVectorStep` → `wave_vector` sub-type
- `numcb/numvb` → `bands` sub-type
- `ExternalField/EFtype/Evalue` → `external_field` sub-type
- `b_field` → `b_field` sub-type
- `feast_emin/emax/m0` → `feast` sub-type
- `materialN` → `material_names`
- `startPos/endPos` → `z_min/z_max`
- `intStartPos/intEndPos` → `int_start_pos/int_end_pos`
- `confDir` → `conf_dir`
- `whichBand` → `which_band`
- `bandIdx` → `band_idx`
- Wire fields (`wire_nx/ny/dx/dy`, `wire_geom`, `numRegions`, `regions`) → `wire_config` sub-type
- Landau fields (`landau_nx/width/sweep`) → `landau_config` sub-type

### Parser rewrite (input_parser.f90)

Replace the entire 1,369-line `read_config` with a toml-f based parser. For this slice, only **bulk mode** must be fully functional. Other modes (QW, wire, Landau) can be stubbed or throw "not yet implemented" errors — they'll be filled in by subsequent issues.

The parser reads `input.toml` (not `input.cfg`). Architecture:
- `read_config(filename)` returns `type(simulation_config)`
- Coordinator calls sub-parser functions per TOML section
- `require_value(table, key, val, section)` for required fields (fail-fast)
- `get_value(table, key, val, default)` for optional fields (silent defaults)
- No `print *` echo of parsed values
- Post-parse computation: `ngrid`, `evnum`, `conf_dir`

### Consumer updates (12 source files, 252 references)

Update every file that references old `simulation_config` field names. The full mapping is in the PRD under "Consumer field migration mapping". Key files:
- `src/apps/main.f90` (79 refs)
- `src/apps/main_gfactor.f90` (55 refs)
- `src/apps/main_optics.f90` (40 refs)
- `src/apps/main_topology.f90` (15 refs)
- `src/core/defs.f90` (29 refs — type definition itself)
- `src/core/simulation_setup.f90` (10 refs)
- `src/physics/sc_loop.f90` (8 refs)
- `src/physics/confinement_init.f90` (7 refs)
- `src/physics/green_functions.f90` (4 refs)
- `src/physics/hamiltonianConstructor.f90` (2 refs)
- `src/io/outputFunctions.f90` (2 refs)
- `src/physics/bdg_hamiltonian.f90` (1 ref)

### Unit test updates (12 test files + 2 parser test rewrites)

Update all pFUnit tests that construct `simulation_config` directly. Rewrite `test_topology_parser.pf` and `test_bdg_config.pf` to write TOML content instead of old format. Add new `test_input_parser.pf` with bulk mode parser tests.

### Config conversion (bulk configs only)

Convert all bulk-mode regression configs to TOML format. This slice only requires bulk configs to work end-to-end.

### Executable filename change

All four executables read `input.toml` instead of `input.cfg`.

## Acceptance criteria

- [ ] `simulation_config` uses new sub-types (`wave_vector_config`, `bands_config`, etc.)
- [ ] `confinement` is `character(len=8)` with string values
- [ ] `ngrid` computed field exists alongside `fd_step` user input
- [ ] No `backspace` calls remain in `input_parser.f90`
- [ ] `input_parser.f90` uses toml-f API exclusively
- [ ] All four executables read `input.toml`
- [ ] No `print *` echo of parsed values
- [ ] `cmake --build build` succeeds (zero compilation errors)
- [ ] Bulk mode configs parse correctly and produce identical physics results to before
- [ ] All existing tests pass (bulk mode end-to-end; other modes may be stubbed)
- [ ] `test_input_parser.pf` has tests for: bulk config parsing, defaults, required field errors, validation errors
- [ ] Old `test_topology_parser.pf` and `test_bdg_config.pf` rewritten for TOML (topology/BdG parsing tested with TOML configs)
- [ ] All consumer files use new field names (no old names remain)

## Blocked by

- Issue 001 (toml-f build integration)
