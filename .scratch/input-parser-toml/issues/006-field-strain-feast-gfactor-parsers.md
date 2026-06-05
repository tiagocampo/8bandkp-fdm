# Issue 006: External field, B-field, strain, FEAST, g-factor parsers

**Type:** AFK
**Blocked by:** Issue 002 (bulk tracer bullet)
**User stories:** US 9, 10, 14, 16

## Parent

PRD: `.scratch/input-parser-toml/PRD.md`

## What to build

Add parsing for the cross-cutting optional sections that apply to multiple geometry modes.

### `[external_field]` section
- `type` (currently only `"EF"` for electric field)
- `value` (electric field strength in kV/cm)

### `[b_field]` section
- `components` — 3-element array `[Bx, By, Bz]` in Tesla
- `g_factor` — optional, default 2.0
- Enables Zeeman splitting when present (sets `bdg%enabled = .true.` for non-Landau modes)

### `[strain]` section
- `substrate` — material name for strain reference (default 0.0 = no strain)

### `[feast]` section
- `emin`, `emax` — energy window (0 = auto)
- `m0` — subspace size (0 = auto: 2*nev)

### G-factor keys (top-level)
- `which_band` — 0 = CB, 1 = VB
- `band_idx` — subband index (default 1)

This eliminates the `b_field`/`EFParams` collision hack from the old parser (65 lines of peek/backspace logic → 2 clean TOML sections).

## Acceptance criteria

- [ ] `[external_field]` parses correctly, value passed to simulation
- [ ] `[b_field]` parses correctly, components set in `b_field_config` and propagated to `bdg%B_vec`
- [ ] `[strain]` substrate parses correctly, enables strain when present
- [ ] `[feast]` parameters parse correctly, 0 = auto behavior preserved
- [ ] `which_band` and `band_idx` parse correctly from top-level keys
- [ ] All configs using these sections converted to TOML
- [ ] Simulations using electric field, B-field, strain produce identical results
- [ ] No backspace calls remain for field/b_field parsing

## Blocked by

- Issue 002 (bulk tracer bullet)
