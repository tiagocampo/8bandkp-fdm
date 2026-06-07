# Issue 012: Documentation overhaul

**Type:** AFK
**Blocked by:** Issue 009 (optics + exciton + scattering parsers — last parser issue)
**User stories:** US 29-33

## Parent

PRD: `.scratch/input-parser-toml/PRD.md`

## What to build

Update all current documentation to reflect the TOML format. No references to the old `label: value` format should remain in active docs. Historical docs (plans, brainstorm archives, ideation) are NOT updated.

### `docs/reference/input-reference.md` — complete rewrite
- Full TOML schema reference with all sections, keys, types, defaults
- TOML examples for each geometry mode
- TOML examples for each optional block
- Validation rules documented

### `docs/reference/output-reference.md`
- Update `input.cfg` → `input.toml` references

### `docs/reference/benchmarks.md`
- Update config filename references

### `README.md`
- Rewrite quickstart section with TOML config example
- Update `input.cfg` → `input.toml`
- Update config directory references

### `docs/lecture/*.md` (15 files)
- Update all inline config snippets from old format to TOML
- Update `input.cfg` → `input.toml`
- Update `cat config.cfg > input.cfg` → `cp config.toml input.toml`
- Update any config format descriptions

### `CLAUDE.md`
- Update "Running" section with `input.toml`
- Update "Input file" section with TOML format description
- Update "Architecture" section: module dependency graph, key design concepts
- Update "Code Conventions" if needed

## Acceptance criteria

- [ ] `docs/reference/input-reference.md` is a complete TOML schema reference
- [ ] `README.md` quickstart uses TOML examples
- [ ] All 15 lecture docs show TOML config snippets
- [ ] No `input.cfg` references remain in active docs (only in historical docs/plans)
- [ ] `CLAUDE.md` accurately describes the new TOML-based input system
- [ ] `docs/reference/output-reference.md` references `input.toml`

## Blocked by

- Issue 009 (optics + exciton + scattering parsers — the last parser issue, needed so docs describe the final complete schema)
