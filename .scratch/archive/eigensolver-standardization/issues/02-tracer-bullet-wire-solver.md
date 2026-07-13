# Tracer bullet: unified interface + wire FEAST with [solver] TOML

**Type:** AFK

## Parent

PRD: `.scratch/eigensolver-standardization/PRD.md`

## What to build

The thin end-to-end vertical slice that proves the entire architecture works for the most important path: wire mode with FEAST solver configured via the new `[solver]` TOML section.

This slice establishes the full architecture in one shot:

**Eigensolver interface redesign** — Redesign `eigensolver_base` with two deferred methods `solve_dense(H, config, result)` and `solve_sparse(H_csr, config, result)`. Add `EIGEN_MODE_FULL`, `EIGEN_MODE_INDEX`, `EIGEN_MODE_ENERGY` named constants. Extend `eigensolver_config` with `mode`, `il`, `iu` fields. Implement `feast_solver_t` for both methods: `solve_sparse` dispatches on mode (ENERGY=native FEAST, FULL=Gershgorin-bounded window, INDEX=error stop). `solve_dense` converts dense→CSR then calls `solve_sparse`. Keep legacy `solve` binding as alias for backward compat. Add `eigensolver_config_validate` helper.

**Config type rename** — Replace `feast_config` with `solver_config` in defs.f90: fields `method` ('AUTO'/'DENSE'/'FEAST'), `mode` ('AUTO'/'FULL'/'INDEX'/'ENERGY'), `emin`, `emax`, `m0`. Update `simulation_config` field name. Update validation checks.

**TOML parser** — Add `parse_solver` reading `[solver]` with method/mode/emin/emax/m0. Keep `parse_feast` with deprecation warning for backward compat. `read_config` tries `[solver]` first, falls back to `[feast]`.

**Wire path dispatch** — Update `simulation_setup.f90` wire path: replace `cfg%feast` with `cfg%solver`, replace `m0 < 0` hack with string method dispatch, implement smart defaults (wire → FEAST + ENERGY), call `eigensolver_config_validate`. Wire call sites in main.f90, main_gfactor.f90, main_optics.f90, main_topology.f90, green_functions.f90, sc_loop.f90 updated from `cfg%feast` to `cfg%solver`.

**Tests** — Unit tests: `test_feast_energy_mode` (FEAST + ENERGY on CSR), `test_feast_full_mode` (FEAST + FULL via Gershgorin), `test_feast_rejects_index` (error stop verification). All existing wire tests pass via legacy `[feast]` backward compat.

**Docs** — Update CLAUDE.md and src/physics/AGENTS.md: replace `feast_config` with `solver_config`, `[feast]` with `[solver]`, remove ARPACK references, add mode/method description.

## Acceptance criteria

- [ ] `eigensolver_base` has deferred `solve_dense` and `solve_sparse` methods; legacy `solve` exists as alias
- [ ] `feast_solver_t` implements both methods with mode dispatch (ENERGY, FULL, INDEX→error stop)
- [ ] `EIGEN_MODE_FULL/INDEX/ENERGY` named constants are public
- [ ] `eigensolver_config` has `mode`, `il`, `iu` fields with correct defaults
- [ ] `eigensolver_config_validate` rejects FEAST+INDEX, validates emin<emax, validates il≤iu
- [ ] `solver_config` type exists in `defs.f90` with method, mode, emin, emax, m0 fields; `feast_config` no longer exists
- [ ] `parse_solver` reads `[solver]` section; `parse_feast` maps `[feast]` with deprecation warning
- [ ] Wire path in `simulation_setup.f90` uses `cfg%solver` with smart defaults (wire → FEAST + ENERGY)
- [ ] All wire/bdG call sites updated from `cfg%feast` to `cfg%solver`
- [ ] 3 new unit tests pass (FEAST energy, FEAST full, FEAST rejects index)
- [ ] CLAUDE.md and AGENTS.md updated for new types and sections
- [ ] All existing tests pass (backward compat via `[feast]` still works)
- [ ] `cmake --build build` succeeds
- [ ] `ctest --test-dir build -j4 --output-on-failure` passes all tests

## Blocked by

- #01 (ARPACK removed to avoid merge conflicts in eigensolver.f90)
