# Implementation Backlog — Ordered by Priority

Consolidated from REVIEW.md on 2026-05-06. Updated 2026-06-03.
Piezoelectric explicitly excluded (ZB [001] = zero by symmetry, wires = negligible).
Phases 1-20 complete. Phase 21 in progress (PR #35). Phases 18-19 remain active backlog (architectural cleanup + extended cross-code).

---

## Phase 1: COMPLETED (2026-05-03)

Closed groups #45, #5, #16, #10, #8. Realistic InAsW(8nm)/GaSbW(8nm) broken-gap config now shows proper 34 meV anticrossing.

---

## Phase 1: Testing Infrastructure

Solidify what exists before adding new physics.

| Source | What | Effort |
|--------|------|--------|
| #4 | 5 gfactor regression tests (bulk CB/VB GaAs, GaAsW, InAsW, QW VB) | Low |
| #23 | Core-shell benchmark: create `wire_inas_gaas_core_shell.cfg` | Low |
| #37 | `validate_simulation_config` + close `contiguous` gaps (3 known sites) | Medium |
| #8 | 3 integration tests (wire hexagon, wire strain, SC wire) + 2 regression datasets | Medium |

**Result:** All existing features have full test coverage.

---

## Phase 2: COMPLETED (2026-05-03)

Closed group #6. Bulk EF shift added to ZB8bandBulk, delta-doping (Gaussian profile) implemented, gfactor SC wired for wire and QW modes. Bulk SC uses QW path (confDir=z, single material). b_field parsing for bulk Landau levels.

**Result:** All 4 executables support full feature set (SC, EF, gfactor, optics, wire).

---

## Phase 3: COMPLETED (2026-05-03)

All 21 review findings resolved: 7 correctness bugs (C1-C7), 9 code quality issues (H1-H12), 5 test coverage gaps (T1-T8), plus 5 pre-existing test failures traced to root cause and fixed. 53/53 tests pass.

**Plan:** `docs/plans/archive/2026-05-03-pr13-review-fixes-plan.md` (archived)
**Commits:** 20 commits from `4342513` to `67a653f`, plus 5 root-cause fixes to `c1fcee2`.

---

## Phase 4: COMPLETED (2026-05-04)

All 74 existing figures validated, ISBT cross-verified (z-dipole vs commutator velocity, rel_diff=0.007%), 5 new physics figures added (bulk E(k), QW subbands, wavefunctions, wire geometry, Zeeman fan). Ch06 + Ch08 rebuilt with verified captions. All 14 lecture chapters have zero broken figure references. Groups #22, #26, #50 -> COMPLETE.

**Plan:** `docs/plans/archive/2026-05-04-phase4-optics-figures-plan.md` (archived)
**Commits:** 12 commits from `88840ce` to `79688ab`.

---

## Phase 5: COMPLETED (2026-05-05)

New confinement=3 (Landau) mode with x-discretized 8NxN Hamiltonian, Landau gauge A=(0,Bz*x,-By*x), Zeeman splitting, orbital quantization, B-sweep fan diagram, analytical validation against E_n = E_C + hw_c(n+1/2). 58/58 tests pass.

**Plan:** `docs/plans/archive/2026-05-04-landau-bulk-phase5.md` (archived)
**Commits:** 82 commits from `70c817b` to `8bd2387`, plus sign fix and CI wiring.

---

## Phase 6: COMPLETED (2026-05-06)

Full topological suite completion repair delivered and pushed.

**Plan:** `docs/plans/archive/2026-05-05-phase6-topological-suite.md` (archived)
**Commits:** `c56fbd4` through `20c3f19`. 66/66 tests passed.

---

## Phase 7: COMPLETED (2026-05-08)

Topological and magnetic bug fixes from Codex review findings, plus ce-code-review follow-up. 5 Codex findings (F1-F5) and 5 ce-code-review findings (U1-U5) fixed. 66/66 tests passed.

**Commits:** `46be100` through `d8fcbbd`.

---

## Phase 8: COMPLETED (2026-05-09)

CSR structure testing — reusable invariant fixture, unit tests for 8 untested CSR operations, Krylov snapshot infrastructure.

**Brainstorm:** `docs/brainstorms/archive/2026-05-09-csr-structure-testing-requirements.md` (archived)
**Plan:** `docs/plans/archive/2026-05-09-001-feat-csr-structure-testing-plan.md` (archived)

*Note: Both brainstorm and plan are in their respective `archive/` directories.*

**Delivered:**
- `tests/support/csr_test_helpers.f90` — structural invariant fixture (5 invariants + helpers)
- `tests/support/krylov_helpers.f90` — Krylov chain + comparison with failure diagnostics
- `tests/support/krylov_reference_data.f90` — reference vectors for 4 code paths
- `tests/support/generate_krylov_reference.f90` + `regenerate_references.sh` — one-command regeneration
- `tests/unit/test_csr_helpers.pf` — 10 fixture tests
- `tests/unit/test_csr_structural.pf` — 10 CSR operation tests
- `tests/unit/test_krylov_helpers.pf` — 4 infrastructure tests
- `tests/unit/test_krylov_snapshots.pf` — 4 snapshot tests (wire, wire+Peierls, BHZ, BdG)
- Refactored existing tests (`test_csr_spmv.pf`, `test_hamiltonian_2d.pf`, `test_bdg_hamiltonian.pf`) to use shared fixture

**Remaining gaps (U6 partial):** Resolved in Phase 14. All 7/7 Krylov snapshot code paths now covered (SC loop, optics wire, g-factor wire added).

---

## Phase 9: COMPLETED (2026-05-09)

8-band verification ladder — 4-rung hierarchy validating the full 8-band Hamiltonian from bulk k=0 through dispersion through confinement. All 19 requirements (R1-R19) implemented.

**Brainstorm:** `docs/brainstorms/archive/2026-05-08-8band-verification-ladder-requirements.md` (archived)
**Plan:** `docs/plans/archive/2026-05-08-feat-8band-verification-ladder-plan.md` (archived)

**Delivered:**
- `tests/integration/verify_8band_rung1_bulk_k0.py` — R1-R5 (eigenvalues, degeneracies, basis ordering, T_d symmetry, normalization) for 5 materials
- `tests/integration/verify_8band_rung2_dispersion.py` — R6-R9 (effective masses vs Kane/Vurgaftman) for 4 materials
- `tests/integration/verify_8band_rung3_qw.py` — R10, R12, R13 (QW subbands, broken-gap, degeneracies)
- `tests/integration/verify_8band_rung4_wire.py` — R14-R16 (wire convergence, dense/sparse agreement, eigenvalue count)
- 14 config files in `tests/regression/configs/`
- Section 7 of `docs/reference/benchmarks.md` with full results
- `docs/solutions/best-practices/8band-verification-ladder-2026-05-09.md`
- All registered under ctest label `verification`

**Result:** COMPLETE. No remaining gaps.

---

## Phase 10: COMPLETED (2026-05-09)

Standard-star benchmark systems — 7 canonical material systems (S1-S7) validated against published observables. 24 code-review findings addressed.

**Brainstorm:** `docs/brainstorms/2026-05-09-standard-star-benchmarks-requirements.md`
**Plans:** `docs/plans/2026-05-09-002-feat-standard-star-benchmarks-plan.md`, `docs/plans/2026-05-09-003-fix-standard-star-review-findings-plan.md`

**Delivered:**
- `tests/integration/star_helpers.py` — shared utility module (run_exe, parse helpers, roth_gfactor, physical constants)
- 7 verification scripts (S1 GaAs bulk through S7 InAs wire, ~3,423 lines total)
- `tests/integration/aggregate_star_benchmarks.py` — aggregation script
- 6 new config files
- All 7 registered under ctest labels `standard-star;verification`

**Remaining gaps (~10-15%):** Resolved in Phase 15. S4 absorption onset has regression reference. S7 wire g-factor tightened. Run wrappers centralized in S3/S5/S6. Roth formula corrected in benchmarks.md.

---

## Phase 11: COMPLETED (2026-05-10)

Executable lecture-test pairs — 14 standalone Python scripts (L00-L13) serving as pedagogical companions and integration tests. Each runs a Fortran executable, validates physics, generates overlay plots.

**Brainstorm:** `docs/brainstorms/archive/2026-05-10-executable-lecture-test-pairs-requirements.md` (archived)
**Plan:** `docs/plans/archive/2026-05-10-004-feat-executable-lecture-test-pairs-plan.md` (archived)

**Delivered:**
- 14 scripts in `scripts/` (lecture_00 through lecture_13)
- 18 overlay plots in `docs/lecture/figures/`
- Makefile targets (`make lecture-XX`, `make lectures`)
- All 14 lecture docs updated with Verification sections and image references
- 4 legacy verify scripts absorbed and deleted
- Shared infrastructure via `tests/integration/star_helpers.py`

**Result:** COMPLETE. No remaining gaps.

---

## Phase 12: COMPLETED (2026-05-10)

Validation coverage matrix — YAML universe file + COVERAGE annotation convention + Python generator + ctest integration.

**Brainstorm:** `docs/brainstorms/archive/2026-05-10-validation-coverage-matrix-requirements.md` (archived)
**Plan:** `docs/plans/archive/2026-05-10-005-feat-validation-coverage-matrix-plan.md` (archived)

**Delivered:**
- `tests/integration/validation_universe.yml` — 59 cells (48 required, 11 aspirational), 28 observables
- `tests/integration/coverage_matrix.py` — generator with heat map + traceability output
- CTest target `coverage_matrix` under label `coverage`
- 50/53 test files annotated with `# COVERAGE:` lines (92 total annotations)
- 3 shell wrappers correctly excluded (delegating scripts carry annotations)

**Remaining gaps:** Resolved in Phase 14. 3 orphan annotations fixed (gfactor, E_sub, state_character).

**Result:** COMPLETE. No remaining gaps.

---

## Phase 13: COMPLETED (2026-05-12)

Strain validation across geometries — bulk, QW, and wire strain verification against published/analytical results.

**Brainstorm:** `docs/brainstorms/2026-05-11-strain-validation-requirements.md`
**Plan:** `docs/plans/2026-05-11-strain-validation-plan.md`

**Delivered:**
- `tests/integration/verify_strain_rung5_bulk.py` — R1-R3 (Bir-Pikus eigenvalue shifts, HH-LH splitting, additive modification)
- `tests/integration/verify_strain_rung6_qw.py` — R4a-R6 (strained QW confinement, HH-LH ordering, grid convergence)
- `tests/integration/verify_strain_wire_profile.py` — R7-R9 (Navier-Cauchy profile, hydrostatic strain, band edge shifts)
- CTest registration: `verification_rung5_strain_bulk`, `verification_rung6_strain_qw`, `strain_validation_wire`
- 4 cells added to `validation_universe.yml` with COVERAGE annotations
- Lecture 04 audit (docstring updated with cross-reference)
- Shared `bir_pikus_biaxial_001()` extracted to `star_helpers.py`
- Code-review fix pass (commit `aff4ca7`)

**Remaining gaps:** All items have documented justification (see Known Limitations below).

---

## Completion Sprint (Phases 14-16): COMPLETED (2026-05-13)

All remaining backlog items resolved. 91/91 tests pass.

### Phase 14: Quick Wins

- P14.1: Added `contiguous` attribute to remaining hot-path arrays (utils.f90, spin_projection.f90, hamiltonianConstructor.f90)
- P14.2: Fixed 3 coverage matrix orphan annotations (gfactor, E_sub, state_character)
- P14.3: Committed strain validation planning docs
- P14.4: Generated gfactor golden data for 3 configs (bulk VB, bulk Winkler CB, QW VB) + registered in CTest

### Phase 15: Validation Tightening

- P15.1: Added 3 CSR Krylov snapshot tests (SC loop, optics wire, gfactor wire) — now 7/7 code paths covered
- P15.2: Tightened standard-star assertions: S4 absorption onset reference, S7 wire g-factor reference, centralized run wrappers in S3/S5/S6, fixed Roth formula in benchmarks.md
- P15.3: Re-scoped docs physics revamp (Backlog #26) — all 12 tasks either completed or superseded by lecture-test pairs

### Phase 16: Integration + Rashba

- P16.1: Added wire hexagon integration test + SC wire integration test (2 new regression tests)
- P16.2: Calibrated Rashba BdG physics — mu at CB subband (0.638 eV), fixed FEAST window, updated lecture_13 B sweep, fixed regression test eigenvalue count extraction

### Additional fixes during the sprint

- Fixed 5 "pre-existing" test failures: CB spacing (112 meV not 9.92), wire convergence (FEAST sparse), band overlap (diagnostic-only), gfactor segfault (contiguous attribute)

---

## Phase 17: COMPLETED (2026-05-14)

Post-completion bug fixes and polish from code review.

**Commits:** `c97c86d` (29 code review findings), `a2ea66a` (failing test), `30c6895` (fit-tail fix), `6223273` (NaN guard + test cleanup).

### 17a: Code Review Findings (29 fixes)

**Plan:** `docs/plans/archive/2026-05-12-001-fix-code-review-findings-plan.md` (archived)
**Brainstorm:** `docs/brainstorms/archive/2026-05-12-code-review-fixes-requirements.md` (archived)

**Delivered:**
- BdG correctness: `pairing_sign_xi` constant, `delta_0 <= 0` guard, Hermiticity tests with Peierls
- Topological: `is_z2_transition` simplified to Z2-change-only, Fu-Kane parity pairing, cumulative Z2 tracking, spectral detection replacing width threshold
- Numerical: `n_fit_actual` regression bias fix, `converged` flag for near-transition, shared `fit_tail_exponential` helper
- Dead code: removed `diag_2x2`, duplicate conductance wrappers, `bdg_config%self_consistent`
- `contiguous` attributes: 30+ hot-path arguments in bdg/topological/magnetic modules
- Python: all `mkdtemp` → `TemporaryDirectory`, `trapz_fn` compat shim, `__file__`-relative paths
- Tests: BdG Hermiticity with Peierls (3 tests), Z2 transition/value tests, `fit_exponential_decay` accuracy

### 17b: Right-Edge Majorana xi Regression Fix

**Plan:** `docs/superpowers/plans/archive/2026-05-13-fix-fit-tail-regression-direction.md` (archived)
**Brainstorm:** `docs/brainstorms/archive/2026-05-13-fit-tail-regression-direction-bug.md` (archived)

**Delivered:**
- Added `forward_tail` boolean to track search direction in `fit_tail_exponential`
- Split regression loop into forward/backward branches with correct distance calculation
- Added fallback for forward tails with zero-density data
- Fixed early-return guard and `domain_extent` for backward tails
- New test: `test_fit_exponential_decay_right_edge` in `test_edge_states.pf`

### 17c: NaN Guard and Test Polish

**Delivered:**
- BdG `delta_0` guard now also rejects NaN (`delta_0 /= delta_0`)
- `is_z2_transition` tests simplified to use 2-argument signature

---

## Remaining Backlog — Phases 1-20 CLOSED, Phases 18-19 Active

| Source | What | Status | Resolution |
|--------|------|--------|------------|
| #4 | 5 gfactor regression tests | COMPLETE | All 5 configs covered: 2 by standard-star (GaAs CB, InAsW CB), 3 by new golden-data regression (bulk VB, GaAsW CB, QW VB) |
| #37 | `validate_simulation_config` | COMPLETE | Implemented as `simulation_config_validate` in `defs.f90` (lines 518-586) |
| #37 | `contiguous` gaps | COMPLETE | All 3 known sites fixed: ZB8bandQW, utils.f90 dns(:,:), spin_projection.f90 psi(:), hamiltonianConstructor.f90 velocity arrays |
| #8 | Integration tests: wire hexagon | COMPLETE | `regression_wire_hexagon` registered in CTest |
| #8 | Integration tests: wire strain | COMPLETE | `verify_strain_wire_profile.py` registered as `strain_validation_wire` |
| #8 | Integration tests: SC wire | COMPLETE | `regression_sc_wire` registered in CTest |
| #26 | Docs physics revamp tasks 3-12 | COMPLETE | All 12 tasks complete or superseded by lecture-test pairs, verification ladder, and standard-star benchmarks |
| #55 | CSR Krylov snapshots (7/7 paths) | COMPLETE | SC loop, optics wire, gfactor wire snapshots added (Phase 15) |
| #56 | Standard-star tightening | COMPLETE | S4 onset reference, S7 g-factor reference, wrapper centralization, Roth formula fix (Phase 15) |
| #38 | Rashba BdG calibration | COMPLETE | mu at CB subband (0.638 eV), FEAST window fixed, B sweep demonstrates phase transition in lecture_13 (Phase 16) |
| #59 | Strain validation | COMPLETE | Bulk/QW/wire scripts, CTest registered, code-review pass applied (Phase 13) |
| #60 | Code review findings (29 fixes) | COMPLETE | BdG/topological/numerical/dead-code/contiguous/Python/test fixes (Phase 17a) |
| #61 | Right-edge Majorana xi regression | COMPLETE | `forward_tail` direction fix + fallback + test (Phase 17b) |
| #62 | NaN guard + test polish | COMPLETE | BdG delta_0 NaN rejection, is_z2_transition test cleanup (Phase 17c) |

---

## Phase 18: Architectural Cleanup + Code Modernization

Source: Deferred review items from PR #27 (`/tmp/8bandkp-fdm-pr27-review-handoff.md`).

### Architectural Debt

| ID | What | Why Deferred | Effort | Priority |
|----|------|-------------|--------|----------|
| I1 | `main.f90` bypasses `simulation_setup` for Landau + k-sweep — refactor to route through it | Works correctly; refactor is architectural debt, not a bug | High | Medium |
| I2 | `defs.f90` at ~1,089 lines — 3.6× the 300-line guideline — decompose toward guideline | Natural split point unclear; don't split without a clear decomposition | Medium | Medium |
| I3 | Duplicate confinement string checks in multiple modules — deduplicate | Low priority — works correctly, just repetitive | Low | Medium |
| I4 | `grid_ngrid`/`grid%npoints()` redundancy — migrate 6 remaining source files (~11 call sites) | Legacy wrapper still called; migration is gradual | Low | Medium |

### Code Cleanup

| ID | What | Why Deferred | Effort | Priority |
|----|------|-------------|--------|----------|
| I13-ext | Replace ~42 remaining `stop 1` calls with `error stop` + messages across 14 files (physics/math/io modules) | Non-executable files; bulk sed replacement. Done for `simulation_setup.f90` + 4 executables in PR #27. | Low | Medium |

### Testing Gaps (incremental)

| ID | What | Why Deferred | Effort | Priority |
|----|------|-------------|--------|----------|
| I17 | ~~No wire optics path tested in any unit test~~ | DONE — `test_optical.pf:test_optical_matrix_wire_basic()` added | — | — |
| I18 | ~~`hamiltonian_blocks` has no dedicated unit test~~ | DONE — `test_hamiltonian_blocks.pf` (177 lines, entry count, duplicate pairs) | — | — |
| I19 | Expand convergence suite to cover all physics modules | Expanding incrementally per feature | Medium | Low |
| I20 | Harden ISBT test assertions further (C8 improved; further hardening incremental) | C8 added proper physics tests; further hardening is incremental | Low | Low |

---

## Phase 19: Extended Cross-Code Validation

Source: `docs/brainstorms/archive/kdotpy-cross-validation-requirements.md` (deferred items)
Blocked on: kdotpy API capabilities (LL symbolic mode, hzy() integration)

| Item | What | Effort | Priority |
|------|------|--------|----------|
| R9 | Cross-code Zeeman comparison (requires kdotpy LL symbolic mode) | High | Low |
| R13-R14 | Cross-code wire comparison (requires kdotpy hzy() integration) | High | Low |
| R15 | QW Landau fan diagram E(B) at 10+ B-field values | Medium | Low |
| R17 | BHZ parameter extraction comparison | Medium | Low |
| R18 | QW (not bulk) g-factor cross-code comparison | Medium | Low |
| R21-R22 | Berry curvature Omega(k) and LL Chern number cross-code comparison | High | Low |
| R24-R26 | Formal discrepancy logging, resolution protocol, layered verification gates | Medium | Low |
| R5/R6 | HgCdTe x=0.3 and CdZnTe alloy materials in parameters.f90 with interpolation | Low | Low |

---

## Phase 20: COMPLETED (2026-06-01)

Deep physics validation — TRK sum rule, Hermiticity, verification ladder rungs 7-8, wire strain quantitative.

**Commits:** `b8f16a2` (delivery), `0582ee6` (review fixes), `a36ef34` (test count, dynamic edge, refactor), `3596737` (wire strain test fix).

### Delivered:

| Item | What | Status |
|------|------|--------|
| Idea 7 | TRK sum rule / f-sum rule invariant checks | **DONE** — `tests/support/trk_helpers.f90` (6 routines) + `tests/unit/test_trk_sum_rule.pf` (7 tests: bulk GaAs/InAs CB+HH, isotropy, QW CSR, wire CSR) |
| Idea 8 | Hermiticity smoke test for all geometries | **DONE** — 20 hermiticity tests across 6 test files (bulk, QW, wire, BdG, Landau, BHZ, blocks) |
| Idea 5 (partial) | Verification ladder rungs 7-8 | **DONE** — `verify_8band_rung7_gfactor.py` (R7.1 bulk Roth 4 materials, R7.2 QW consistency) + `verify_8band_rung8_optical.py` (R8.1 Kane self-consistency, R8.2 absorption edge 10meV, R8.3 TE/TM ordering) |
| Idea 6 (partial) | Wire strain quantitative | **DONE** — `verify_strain_wire_quantitative.py` (R10a gap direction, R10b HH>LH, R10c measurable VB shift, R11 quantitative gap shift 15% tolerance) |

### Code review fixes applied:
- Dynamic absorption edge search window (replaces hard-coded 1.5–2.0 eV with `e_transition ± 0.3 eV`)
- Extracted `diagonalize_bulk_at_k` helper in `trk_helpers.f90` (eliminates 3× boilerplate)
- Wire strain test: 25×25 grid, narrowed FEAST window [-0.5, 1.0] eV, R10 qualitative + R11 quantitative

---

## Phase 21: IN_PROGRESS (PR #35)

Publishable validation benchmarks — ISBT oscillator strength, SC Schrödinger-Poisson, exciton Bastard reference.

**PRD:** `.scratch/publishable-benchmarks/PRD.md`
**Issues:** `.scratch/publishable-benchmarks/issues/` (6 issues)
**Branch:** `feat/publishable-benchmarks-phase21`
**PR:** #35 (not yet merged)
**Commits:** `c36cfed` through `d28bd0f` (3 feature + 2 review-fix commits)

### Delivered:

| Item | What | Status |
|------|------|--------|
| Issue #01 | ISBT GaAs/AlGaAs tracer bullet (100A) | **DONE** — `isbt_gaas_algaas_w100.toml`, infinite-well analytical comparison |
| Issue #02 | SC charge neutrality hard check <1% | **DONE** — `verify_sc_benchmark.py` Check 2 |
| Issue #03 | Exciton Miller tolerance tightening (30%→15%) + Bastard 1982 reference | **DONE** — `test_exciton_convergence.py` updated |
| Issue #04 | ISBT width sweep (50, 100, 200 A GaAs/AlGaAs) | **DONE** — 3 TOML configs + 3-width loop in benchmark |
| Issue #05 | SC Fermi level, subband shift, potential profile checks | **DONE** — `verify_sc_benchmark.py` Checks 3-5 |
| Issue #06 | ISBT InGaAs/InAlAs material sweep (50, 100, 200 A) | **DONE** — 3 TOML configs + cross-material comparison |

### Verification scripts:

- `verify_isbt_benchmark.py` — 6-config sweep, 7 checks per config (B1-B7), cross-material comparison, summary table
- `verify_sc_benchmark.py` — 5 checks (SC convergence, charge neutrality, Fermi level, subband shift, potential profile)
- `test_exciton_convergence.py` — Updated: Miller 15%, Bastard 1982 (8.5 meV, 20%)

### Infrastructure updates:

- `validation_universe.yml`: +ISBT_dipole, +ISBT_oscillator_strength, +fermi_level QW cells, ISBT promoted to required
- `coverage_matrix.py`: Added `test_*.py` scan pattern
- `tests/CMakeLists.txt`: 2 new tests (TIMEOUT 1200/900)

### Code review fixes (2 rounds):

- Round 1 (`8f687fa`): Fixed f12_inf mass-dependence, added B6/B7 checks, fixed SC docstring, added universe cells
- Round 2 (`d28bd0f`): Fixed B6/B7 docstring sync, increased CTest TIMEOUTs, minor cleanup

### Deferred items (follow-up PRs):

| Item | Description | Recommended PR |
|------|-------------|----------------|
| US19/US25 | Benchmark figures (matplotlib) | `feat/publishable-figures` |
| US24 | ISBT golden data regression tests | `feat/isbt-golden-data` |
| US3/US4 | Well-width monotonic trend assertion | `feat/isbt-width-trends` |
| US18 | Multi-width exciton sweep | `feat/exciton-width-sweep` |
| US10/US11 | Tighter Fermi/subband tolerances | `feat/sc-tolerance-tightening` |

---

## Minor Carry-Over Items

| Source | What | Effort | Priority |
|--------|------|--------|----------|
| Richardson R9 gap | Investigate why absorption_edge observable absent from stored JSON results despite code being complete | Low | Medium |
| PR14 review T10 gap | Replace return None in test_bulk_zeeman.py:44 with RuntimeError (consistency) | Trivial | Low |
| Phase 21 PRD | US19/US25: Benchmark figures (matplotlib) matching lecture-script pattern | Medium | Medium |
| Phase 21 PRD | US24: ISBT golden data regression tests | Medium | Medium |
| Phase 21 PRD | US3/US4: Well-width monotonic trend + hard cross-material assertion | Low | Medium |
| Phase 21 PRD | US18: Multi-width exciton sweep (30, 50, 80, 100, 150, 200 A) | Medium | Low |
| Phase 21 PRD | US10/US11: Tighter Fermi/subband tolerances (need baseline 8-band values) | Low | Low |

---

## Known Limitations

Items that were explicitly deferred or relaxed with documented justification:

- **QW strain validation** uses qualitative checks instead of Bastard analytical formulas (justified: narrow 20A well, Bastard formula assumes infinite barriers)
- **Wire strain R7 tolerance** relaxed from 5-10% to 60% (justified: free-surface relaxation produces non-trivial displacement profiles)
- **`config_validation_result` helper type** not implemented (low priority — `simulation_config_validate` covers the functional need)

---

## Summary

| Phase | Groups resolved | Net new capability | Status |
|-------|----------------|-------------------|--------|
| 1. Testing | 4 | Full regression coverage | DONE |
| 2. Physics wiring | 1 | Bulk SC, bulk EF, gfactor SC | DONE |
| 3. PR review fixes | — | Correctness, compat, quality | DONE |
| 4. Figures | 3 | Complete documentation | DONE |
| 5. Peierls/Landau | 2 | Magnetic field for all modes | DONE |
| 6. Topological | 3 | Full topo suite | DONE |
| 7. Bug fixes + review | 2 | Codex + ce-code-review findings closed | DONE |
| 8. CSR testing | — | Structural invariant fixture + Krylov snapshots (7/7 paths) | DONE |
| 9. Verification ladder | — | 4-rung 8-band hierarchy (R1-R19) | DONE |
| 10. Standard-star benchmarks | — | 7 material systems S1-S7, assertions tightened | DONE |
| 11. Lecture-test pairs | — | 14 executable lectures | DONE |
| 12. Coverage matrix | — | YAML universe + annotations + generator | DONE |
| 13. Strain validation | — | Bulk/QW/wire strain verification | DONE |
| 14. Quick wins | — | Contiguous, coverage orphans, gfactor golden data | DONE |
| 15. Validation tightening | — | Krylov 7/7, standard-star tightening, docs revamp re-scope | DONE |
| 16. Integration + Rashba | — | Wire hexagon, SC wire, Rashba BdG calibrated | DONE |
| 17. Post-completion fixes | — | 29 code review findings, fit-tail regression fix, NaN guard | DONE |
| 18. Architectural cleanup | — | Landau refactor, defs decomposition, grid accessor migration, stop 1 replacement | PENDING |
| 19. Extended cross-code | — | Wire/Landau/g-factor/strain/Berry/Chern cross-code validation | PENDING (blocked on kdotpy API) |
| 20. Deep physics validation | — | TRK sum rule (7 tests), Hermiticity (20 tests), rungs 7-8 (g-factor + optical), wire strain quantitative | DONE |
| 21. Publishable benchmarks | — | ISBT 6-config sweep (42 checks), SC benchmark (5 checks), exciton tightening + Bastard 1982 | IN_PROGRESS (PR #35) |

**Phases 1-20 complete.** 113 tests pass. Phase 21 in progress (PR #35, branch `feat/publishable-benchmarks-phase21`). Phases 18-19 are active backlog.
Source of truth: `docs/plans/BACKLOG.md` (this file) and `docs/plans/REVIEW.md`.
