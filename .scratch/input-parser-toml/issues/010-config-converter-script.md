# Issue 010: Config converter script

**Type:** AFK
**Blocked by:** Issue 002 (bulk tracer bullet)
**User stories:** US 27

## Parent

PRD: `.scratch/input-parser-toml/PRD.md`

## What to build

Write a Python script (`scripts/convert_cfg_to_toml.py`) that automates converting old `label: value` config files to TOML format. This script handles all 4 geometry modes and all optional blocks.

The converter must:
- Parse the old positional `label: value` format
- Detect confinement mode from the `confinement` field (0/1/2/3)
- Generate correct TOML output using the schema from the PRD:
  - Top-level singletons: `confinement` (integer→string), `FDorder`, `fd_step`
  - `[wave_vector]` section
  - `[bands]` section
  - `[[material]]` or `[wire]`/`[landau]` based on mode
  - All optional sections as detected
- Handle edge cases: polygon vertices, delta doping, topology nested sub-sections
- Accept input filename as argument, output to stdout or specified file
- Support batch mode: convert all `.cfg` files in a directory

This can be developed in parallel with issues 3-9 but is needed for full verification.

## Acceptance criteria

- [ ] `python3 scripts/convert_cfg_to_toml.py tests/regression/configs/bulk_gaas_k0.cfg` produces valid TOML
- [ ] Converter handles all 4 confinement modes correctly
- [ ] Converter handles all optional blocks (SC, topology, optics, etc.)
- [ ] Converter handles polygon vertices, delta doping
- [ ] Batch mode converts all 88 regression configs
- [ ] Round-trip verification: converted configs produce identical simulation results

## Blocked by

- Issue 002 (bulk tracer bullet) — needs the TOML parser to exist for round-trip verification
