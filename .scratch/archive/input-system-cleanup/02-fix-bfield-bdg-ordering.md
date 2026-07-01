# Issue 2: Fix b_field/BdG implicit ordering dependency

## Parent

PRD: Input System Architecture Cleanup (`.scratch/input-system-cleanup/PRD.md`)

## What to build

Move the unconditional b_field→BdG copy from `parse_b_field` into `parse_bdg` as default initialization. Currently `parse_b_field` writes `cfg%bdg%B_vec = cfg%b_field%components` and `cfg%bdg%g_factor = cfg%b_field%g_factor`. This creates an implicit ordering dependency: parse_b_field must be called before parse_bdg, and reordering the parse calls in `read_config` would silently break BdG defaults.

After this change, `parse_b_field` only writes to `cfg%b_field`. `parse_bdg` sets its own defaults from `cfg%b_field%components` and `cfg%b_field%g_factor` before reading the `[bdg]` TOML section, then overrides with explicit TOML values if present. The existing comment "Read B_vec from [bdg] if present (overrides b_field copy)" already describes this pattern — just move the defaults to the right place.

## Acceptance criteria

- [ ] `parse_b_field` does not write to `cfg%bdg` fields
- [ ] `parse_bdg` sets `cfg%bdg%B_vec` and `cfg%bdg%g_factor` from `cfg%b_field` as defaults before reading TOML overrides
- [ ] Swapping the order of `parse_b_field` and `parse_bdg` calls in `read_config` produces identical results
- [ ] BdG regression tests pass (configs with [bdg] section still work)
- [ ] BdG configs without explicit B_vec still inherit from [b_field]

## Blocked by

None — can start immediately.
