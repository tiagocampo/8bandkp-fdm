# Issue 004: Wire geometry parser — grid, shapes, regions

**Type:** AFK
**Blocked by:** Issue 002 (bulk tracer bullet)
**User stories:** US 7

## Parent

PRD: `.scratch/input-parser-toml/PRD.md`

## What to build

Add `[wire]` section parsing with all shape variants and `[[region]]` array-of-tables. This enables 2D wire confinement mode.

The parser must:
- Parse `[wire]` section: `nx`, `ny`, `dx`, `dy`
- Parse `[wire.geometry]` sub-section with shape-dependent fields:
  - `circle`/`hexagon`: `radius`
  - `rectangle`: `width`, `height`
  - `polygon`: `vertices` (array of [x, y] pairs)
- Parse `[[region]]` entries with `material`, `inner`, `outer`
- Compute `ngrid = wire%ny` for wire mode
- Compute dummy `start_pos`/`end_pos` from wire geometry
- Validate: `nx >= 3`, `ny >= 3`, `dx > 0`, `dy > 0`, `nx >= FDorder + 1`, `ny >= FDorder + 1`

Convert all wire-mode regression configs to TOML.

## Acceptance criteria

- [ ] Wire TOML configs with `[wire]` + `[[region]]` parse correctly
- [ ] All 4 shape variants (circle, rectangle, hexagon, polygon) parse correctly
- [ ] Polygon vertex arrays converted from TOML to `wire_geometry%verts(2, nverts)`
- [ ] `ngrid` computed as `wire%ny` for wire mode
- [ ] Validation catches: nx/ny < 3, dx/dy <= 0, nx/ny < FDorder + 1
- [ ] All wire regression configs converted to TOML
- [ ] Wire simulation outputs match pre-migration golden data
- [ ] Parser test has wire-specific test cases for each shape variant

## Blocked by

- Issue 002 (bulk tracer bullet)
