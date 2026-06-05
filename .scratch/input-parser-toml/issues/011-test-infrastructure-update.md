# Issue 011: Test infrastructure update

**Type:** AFK
**Blocked by:** Issue 010 (config converter script)
**User stories:** US 26, 27

## Parent

PRD: `.scratch/input-parser-toml/PRD.md`

## What to build

Update all test shell scripts and lecture companion scripts to use `input.toml` and the new TOML config files.

### Integration test shell scripts (32 files in `tests/integration/`)
- Change `cp $config input.cfg` → `cp $config input.toml`
- Update any config filename references from `.cfg` to `.toml`

### Lecture companion scripts (15 files in `scripts/`)
- Update config file references from `.cfg` to `.toml`
- Update any inline config generation to write TOML format
- Update `run_exe` helper if it hard-codes `input.cfg`

### Validation pipeline (`validation/`)
- Convert all 12 validation configs to TOML
- Update Python runners to reference `.toml` files

### Regression test infrastructure
- Ensure `compare_output.py` still works (it compares numerical output, not config format)
- Update any config-related assertions

## Acceptance criteria

- [ ] All 32 integration test shell scripts use `input.toml`
- [ ] All 15 lecture scripts reference `.toml` configs
- [ ] All 12 validation configs converted to TOML and runners updated
- [ ] `ctest --test-dir build -j4` passes all tests
- [ ] `ctest --test-dir build -L regression` passes all regression tests
- [ ] `ctest --test-dir build -L verification` passes all verification tests

## Blocked by

- Issue 010 (config converter script) — configs must be in TOML before scripts can reference them
