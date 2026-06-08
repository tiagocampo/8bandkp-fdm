# Config migration + docs + comprehensive tests

**Type:** AFK

## Parent

PRD: `.scratch/eigensolver-standardization/PRD.md`

## What to build

Final vertical slice: migrate all remaining `[feast]` TOML sections to `[solver]`, remove backward compatibility, update all documentation, and add comprehensive test coverage for all solver × mode × geometry combinations.

**Config migration** — Replace all `[feast]` TOML sections with `[solver]` across 17 regression configs, 3 integration test Python files (inline TOML string literals), and 4 Python scripts. Migration pattern:
- `m0 = -1` (dense fallback) → `method = "DENSE"`, `mode = "ENERGY"`
- `m0 = 0` (auto FEAST) → `method = "FEAST"`, `mode = "ENERGY"`
- `m0 > 0` (explicit FEAST subspace) → `method = "FEAST"`, `mode = "ENERGY"`, `m0 = <value>`

**Remove backward compat** — Delete `parse_feast` subroutine from `input_parser.f90`. Remove `[feast]` fallback in `read_config`. Only `[solver]` accepted.

**Documentation** — Update all docs to reflect the new `[solver]` section:
- `docs/reference/input-reference.md`: rewrite FEAST section → `[solver]` section with method/mode/emin/emax/m0, valid combinations table, smart defaults table
- `docs/lecture/08-quantum-wire.md`: TOML examples, FEAST→solver references
- `docs/lecture/09-numerical-methods.md`: eigensolver section, remove ARPACK
- `docs/lecture/05-gfactor.md`, `06-optical-properties.md`, `11-convergence.md`, `12-extending-the-code.md`, `13-topological-superconductivity.md`: solver references
- `CLAUDE.md`, `src/physics/AGENTS.md`: feast_config→solver_config, `[feast]`→`[solver]`, remove ARPACK, add mode/method
- `README.md`: TOML example

**Comprehensive tests:**
- `test_solver_config_validation`: FEAST+INDEX rejection, emin≥emax, il>iu, valid combos
- `verify_solver_defaults.py`: runs each confinement mode without `[solver]`, verifies smart defaults produce correct results
- `verify_qw_sparse_solver.py`: QW with FEAST matches QW with DENSE
- Coverage annotations on new tests

## Acceptance criteria

- [ ] All 17 regression TOML configs use `[solver]` (zero `[feast]` sections remain)
- [ ] All 3 integration test Python files use `[solver]` string literals
- [ ] `convert_cfg_to_toml.py` emits `[solver]` section
- [ ] `parse_feast` deleted; no `[feast]` backward compat in parser
- [ ] `docs/reference/input-reference.md` has complete `[solver]` documentation
- [ ] All 6 lecture docs updated (08, 09, 05, 06, 11, 12, 13)
- [ ] CLAUDE.md and AGENTS.md updated
- [ ] README.md updated
- [ ] `test_solver_config_validation` passes
- [ ] `verify_solver_defaults.py` passes for all 4 confinement modes
- [ ] `verify_qw_sparse_solver.py` passes
- [ ] Coverage annotations added
- [ ] `cmake --build build` succeeds
- [ ] `ctest --test-dir build -j4 --output-on-failure` passes all tests (113+ existing + new)

## Blocked by

- #05 (all call sites migrated first — configs and docs must reflect the final state)
