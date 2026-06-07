# Issue 013: Full verification + cleanup

**Type:** AFK
**Blocked by:** Issue 011 (test infrastructure), Issue 012 (documentation)
**User stories:** US 27

## Parent

PRD: `.scratch/input-parser-toml/PRD.md`

## What to build

Run the complete verification suite to confirm the migration is correct, then clean up any remaining artifacts.

### Full test suite
- `ctest --test-dir build -j4 --output-on-failure` — all ~91 tests
- `ctest --test-dir build -L unit` — all pFUnit unit tests
- `ctest --test-dir build -L regression` — all regression tests
- `ctest --test-dir build -L verification` — all verification tests
- `ctest --test-dir build -L standard-star` — standard-star benchmarks
- `ctest --test-dir build -L strain-validation` — strain validation
- `ctest --test-dir build -L coverage` — validation coverage matrix

### Validation pipeline
- `python3 validation/run_all.py` — all 12 cross-code validation tests

### Lecture scripts
- Run all 15 lecture scripts to verify they work with TOML configs

### Cleanup
- Remove any remaining `.cfg` files in `tests/regression/configs/` (should all be `.toml` now)
- Remove any dead code from the old parser
- Verify no `backspace` calls remain anywhere in the codebase
- Verify no `input.cfg` references remain in source code or active docs
- Update `.gitignore` if needed

## Acceptance criteria

- [ ] All unit tests pass
- [ ] All regression tests pass with identical numerical results
- [ ] All verification rungs pass
- [ ] All standard-star benchmarks pass
- [ ] All strain validation tests pass
- [ ] Validation coverage matrix reports correctly
- [ ] All 12 cross-code validation tests pass
- [ ] All 15 lecture scripts run successfully
- [ ] No `.cfg` files remain in `tests/regression/configs/`
- [ ] No `backspace` calls in `src/io/input_parser.f90`
- [ ] No `input.cfg` references in source code or active docs
- [ ] `grep -r 'backspace' src/io/input_parser.f90` returns empty
- [ ] `grep -r 'input\.cfg' src/` returns empty
