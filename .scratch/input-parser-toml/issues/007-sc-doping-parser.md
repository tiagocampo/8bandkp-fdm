# Issue 007: SC + doping parser

**Type:** AFK
**Blocked by:** Issue 002 (bulk tracer bullet)
**User stories:** US 11, 13

## Parent

PRD: `.scratch/input-parser-toml/PRD.md`

## What to build

Add `[sc]` section and `[[doping]]` array-of-tables parsing for the self-consistent Schrödinger-Poisson solver.

### `[sc]` section
- All SC parameters: `max_iterations`, `tolerance`, `mixing_alpha`, `diis_history`, `temperature`, `fermi_mode` (string: `"charge_neutrality"` or `"fixed"`), `fermi_level`, `num_kpar`, `kpar_max`, `bc_type`, `bc_left`, `bc_right`
- Section presence = enabled (no `enabled` key)
- `fermi_mode` changes from integer (0/1) to string (`"charge_neutrality"`/`"fixed"`)

### `[[doping]]` array-of-tables
- Uniform doping: `layer`, `ND`, `NA`
- Delta doping: `layer`, `type = "delta"`, `NS`, `fwhm`, `pos`
- `type` defaults to `"uniform"` when omitted
- Number of doping entries can differ from number of material layers

Convert all SC-mode regression configs to TOML.

## Acceptance criteria

- [ ] `[sc]` section parses correctly, all 12 parameters extracted
- [ ] `fermi_mode` is string enum `"charge_neutrality"` or `"fixed"`
- [ ] `[[doping]]` parses uniform and delta doping types
- [ ] Delta doping correctly sets `dtype`, `NS`, `delta_fwhm`, `delta_pos`
- [ ] All SC regression configs converted to TOML
- [ ] SC simulation outputs match pre-migration golden data
- [ ] Parser test has SC-specific test cases including both doping types

## Blocked by

- Issue 002 (bulk tracer bullet)
