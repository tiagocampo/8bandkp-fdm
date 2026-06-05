# Issue 003: QW geometry parser — materials, fd_step, grid computation

**Type:** AFK
**Blocked by:** Issue 002 (bulk tracer bullet)
**User stories:** US 5, 23

## Parent

PRD: `.scratch/input-parser-toml/PRD.md`

## What to build

Add `[[material]]` array-of-tables parsing and QW grid computation to the TOML parser. This enables the most common simulation mode — quantum well confinement.

The parser must:
- Parse `[[material]]` entries with `name`, `z_min`, `z_max` fields
- Compute `num_layers` from the count of `[[material]]` entries
- Compute `ngrid = fd_step` for QW mode
- Compute the z-grid: `z(:)`, `delta`, `dz`, `total_size`
- Compute `int_start_pos` / `int_end_pos` from layer positions
- Call `paramDatabase` to populate `params(:)` for each material
- Validate: `fd_step >= 3`, `fd_step >= FDorder + 1`, at least 1 material layer

Convert all QW-mode regression configs to TOML. Verify QW simulations produce identical results to before.

## Acceptance criteria

- [ ] QW TOML configs with `[[material]]` entries parse correctly
- [ ] Grid computation (z, delta, dz, int_start_pos, int_end_pos) matches old parser output
- [ ] `ngrid` computed as `fd_step` for QW mode
- [ ] Validation catches: fd_step < 3, fd_step < FDorder + 1
- [ ] All QW regression configs converted to TOML
- [ ] QW simulation outputs match pre-migration golden data
- [ ] Parser test `test_input_parser.pf` has QW-specific test cases

## Blocked by

- Issue 002 (bulk tracer bullet)
