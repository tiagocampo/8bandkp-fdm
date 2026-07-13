# Issue 005: Landau geometry parser — grid, sweep, material

**Type:** AFK
**Blocked by:** Issue 002 (bulk tracer bullet)
**User stories:** US 8

## Parent

PRD: `.scratch/input-parser-toml/PRD.md`

## What to build

Add `[landau]` section parsing for Landau level mode (1D x-discretization for orbital quantization).

The parser must:
- Parse `[landau]` section: `nx`, `width`, `sweep` (optional, default `"ky"`), `material`
- Compute `ngrid = landau%nx` for Landau mode
- Compute `dz = width / (nx - 1)`
- Compute `start_pos`/`end_pos` from width
- Set `conf_dir = 'x'` for Landau mode
- Validate: `nx >= 3`, `width > 0`, `nx >= FDorder + 1`

Convert all Landau-mode regression configs to TOML.

## Acceptance criteria

- [ ] Landau TOML configs with `[landau]` parse correctly
- [ ] `ngrid` computed as `landau%nx` for Landau mode
- [ ] Default sweep mode is `"ky"` when omitted
- [ ] Validation catches: nx < 3, width <= 0, nx < FDorder + 1
- [ ] All Landau regression configs converted to TOML
- [ ] Landau simulation outputs match pre-migration golden data
- [ ] Parser test has Landau-specific test cases

## Blocked by

- Issue 002 (bulk tracer bullet)
