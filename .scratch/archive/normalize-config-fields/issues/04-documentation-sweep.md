# Issue 4: Documentation sweep

**Type**: AFK
**Blocked by**: Issues 1–3 (all code must land first)
**User stories**: 10, 11, 12, 13, 14

## What to build

Update all documentation to reflect the normalized config type. No code changes — purely prose and comment updates.

Files to update:
- `docs/reference/input-reference.md` — Notes section: replace "computed `ngrid` field" with `grid%npoints()` explanation
- `docs/reference/output-reference.md` — eigenvalue count formulas: replace `ngrid` with `grid%npoints()`
- `docs/lecture/00-quickstart.md` — update `fdStep` warning narrative and table entry
- `docs/lecture/02-quantum-well.md` — replace `confDir = 'z'` reference with `conf_direction()` explanation
- `docs/lecture/05-gfactor.md` — replace `fdStep`/`ngrid` prose with `grid%npoints()` terminology
- `docs/lecture/07-self-consistent-sp.md` — update any `ngrid` references
- `docs/lecture/08-quantum-wire.md` — update backward-compatibility note (fdStep/confDir no longer exist as stored fields), update confDir reference
- `docs/lecture/10-qcse.md` — update any `ngrid` references
- `docs/lecture/12-extending-the-code.md` — replace `Ngrid`/`fdStep` with `grid%npoints()`
- `CLAUDE.md` — update architecture section: remove `ngrid` from dual-mode description, update `fd_step` description, add `conf_direction` to key concepts

What does NOT change:
- `README.md` TOML examples — `fd_step` remains a valid TOML key, examples unchanged
- All 88 TOML configs — unchanged
- All lecture scripts — unchanged (`fd_step` in generated TOML strings still valid)
- Archived plan documents — historical reference only

## Acceptance criteria

- [ ] `docs/reference/input-reference.md` Notes section updated
- [ ] `docs/reference/output-reference.md` formulas updated
- [ ] All 6 lecture chapters updated (00, 02, 05, 07, 08, 10, 12)
- [ ] `CLAUDE.md` architecture section updated
- [ ] No stale references to `cfg%ngrid`, `cfg%conf_dir`, or mode-dependent `num_layers` in any documentation file (grep confirms)
- [ ] README.md examples remain valid and unchanged
- [ ] No source code changes in this issue — documentation only

## Blocked by

- Issue 1 (conf_direction)
- Issue 2 (ngrid elimination)
- Issue 3 (num_layers split)
