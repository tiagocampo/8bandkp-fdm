# Plans Review — Implementation Status

Review date: 2026-06-01 (updated from 2026-05-30)

| # | Group | Files | Status | Notes |
|---|-------|-------|--------|-------|
| 1 | 2026-03-28-codebase-improvement | design | COMPLETE | Archived |
| 2 | 2026-03-29-cmake-testing | design | COMPLETE | Archived |
| 3 | 2026-03-29-gfactor-improvements | design | COMPLETE | Archived |
| 4 | 2026-03-29-gfactor-validation | standalone | COMPLETE | All 5 configs covered: 2 by standard-star benchmarks, 3 by golden-data regression tests (Phase 14). See details below. |
| 5 | 2026-03-29-pr-review-fixes | design | COMPLETE | Archived |
| 6 | 2026-03-29-self-consistent-sp | design | COMPLETE | Archived |
| 7 | 2026-03-30-sc-validation-report | standalone | COMPLETE | Archived |
| 8 | 2026-04-02-quantum-wire | design | COMPLETE | Archived. Integration tests added in Phase 16 (wire hexagon, SC wire). |
| 9 | 2026-04-03-linalg-backend | design | COMPLETE | Archived |
| 10 | 2026-04-04-documentation-overhaul | design | COMPLETE | Archived |
| 11 | 2026-04-04-pr-review-fixes | design | COMPLETE | Archived |
| 12 | 2026-04-06-bulk-chapter-overhaul | design | COMPLETE | Archived |
| 13 | 2026-04-06-documentation-deep-rewrite | standalone | COMPLETE | Archived |
| 14 | 2026-04-06-quickstart-fixes | design | COMPLETE | Archived |
| 15 | 2026-04-07-strain-unification | design + plan | COMPLETE | Archived |
| 16 | 2026-04-12-qw-documentation-overhaul | design + 2 plans | COMPLETE | Archived |
| 17 | 2026-04-17-fd-variable-coeff-fix | design + plan | COMPLETE | Archived |
| 18 | 2026-04-18-ch01-03-review-fix | design + plan | COMPLETE | Archived |
| 19 | 2026-04-18-ch04-06-review-fix | plan | COMPLETE | Archived |
| 20 | 2026-04-18-ch04-strain-figures | design + plan | COMPLETE | Archived |
| 21 | 2026-04-19-ch03-ch04-figures | plan | COMPLETE | Archived |
| 22 | 2026-04-19-qw-tutorials | design | COMPLETE | All optics figure functions + 5 new Group #50 figures. ISBT cross-validated |
| 23 | 2026-04-20-remaining-work | design + plan | COMPLETE | Archived |
| 24 | 2026-04-20-strain-sign-and-buffer-fix | standalone | COMPLETE | Archived |
| 25 | 2026-04-21-benchmark-matrix | standalone | COMPLETE | Archived |
| 26 | 2026-04-21-docs-physics-revamp | design + plan | COMPLETE | All 12 tasks DONE or superseded by lecture-test pairs (Phase 15 re-scope confirmed). See details below. |
| 27 | 2026-04-21-failure-ledger | standalone | COMPLETE | Archived |
| 28 | 2026-04-21-figure-provenance | standalone | COMPLETE | Archived |
| 29 | 2026-04-21-nextnano-bulk-qw-comparison | standalone | COMPLETE | Archived |
| 30 | 2026-04-21-nextnano-qcse-optics-comparison | standalone | COMPLETE | Archived |
| 31 | 2026-04-21-nextnano-wire-comparison | standalone | COMPLETE | Archived |
| 32 | 2026-04-21-validation-summary | standalone | COMPLETE | Archived |
| 33 | 2026-04-21-wire-validation-summary | standalone | COMPLETE | Archived |
| 34 | 2026-04-25-optical-properties | redesign + implementation | COMPLETE | Archived |
| 35 | 2026-04-25-wire-gfactor-commutator | design | COMPLETE | Archived |
| 36 | 2026-04-26-hamiltonian-performance-refactor | design + plan | COMPLETE | Archived |
| 37 | 2026-04-26-modern-fortran-migration | design + plan | COMPLETE | `validate_simulation_config` DONE (as `simulation_config_validate`). All contiguous gaps closed in Phase 14 (utils.f90, spin_projection.f90, hamiltonianConstructor.f90). |
| 38 | 2026-04-27-bdg-topological-sc | design + impl plan | COMPLETE | Archived. Rashba BdG calibrated in Phase 16 (mu at CB subband, FEAST window, B sweep). |
| 39 | 2026-04-27-pr11-review-fixes | design + plan | COMPLETE | Archived |
| 40 | 2026-04-27-wire-fastpath-followup | doc + plan | COMPLETE | Archived |
| 41 | 2026-04-28-modern-fortran-phase2 | plan | COMPLETE | Archived |
| 42 | 2026-04-28-pr11-cache-safety-fixes | plan | COMPLETE | Archived |
| 43 | 2026-04-29-modern-fortran-remaining | plan | COMPLETE | Archived |
| 44 | 2026-04-29-phase4-exploratory | standalone | COMPLETE | Archived |
| 45 | 2026-04-29-pr12-fixes | plan | COMPLETE | Archived |
| 46 | 2026-05-01-bhz-remaining-work | plan | COMPLETE | Archived |
| 47 | 2026-05-02-magnetic-field-landau | design + impl plan | COMPLETE | Archived |
| 48 | 2026-05-02-physics-figures-implementation | plan | COMPLETE | Archived |
| 49 | 2026-05-02-topological-suite-verification | spec + plan | COMPLETE | Archived |
| 50 | 2026-05-03-physics-figures-extended | plan | COMPLETE | Archived |
| 51 | 2026-05-05-phase6-completion-repair | plan | COMPLETE | Archived |
| 52 | 2026-05-07-topological-magnetic-bugfixes | spec + plan | COMPLETE | Archived |
| 53 | 2026-05-08-review-findings-fixes | plan | COMPLETE | Archived |
| 54 | 2026-05-08-8band-verification-ladder | brainstorm + plan | COMPLETE | Archived. 4 rungs (R1-R19), 14 configs, benchmarks.md Section 7 |
| 55 | 2026-05-09-csr-structure-testing | brainstorm + plan | COMPLETE | Archived. U1-U5 complete. U6: 7/7 Krylov paths done (SC loop, optics wire, gfactor wire added in Phase 15). |
| 56 | 2026-05-09-standard-star-benchmarks | brainstorm + 2 plans | COMPLETE | 7 scripts S1-S7, 24 review findings addressed. S4 onset tightened, S7 g-factor tightened, wrappers centralized, Roth formula fixed (Phase 15). |
| 57 | 2026-05-10-executable-lecture-test-pairs | brainstorm + plan | COMPLETE | Archived. 14 scripts, 18 plots, Makefile targets |
| 58 | 2026-05-10-validation-coverage-matrix | brainstorm + plan | COMPLETE | Archived. 59 cells, 92 annotations, ctest `coverage` label |
| 59 | 2026-05-11-strain-validation | brainstorm + plan | COMPLETE | Bulk/QW/wire scripts done, CTest registered. Code-review pass applied. Docs committed. Tolerances have documented justification. |
| 60 | 2026-05-12-code-review-fixes | brainstorm + plan | COMPLETE | 29 findings resolved: BdG correctness (pairing_sign_xi, delta_0 guard, Hermiticity), topological (Z2-only transition, Fu-Kane parity, cumulative Z2), numerical (n_fit_actual, converged flag), dead code removal, contiguous attributes, Python cleanup. Commit `c97c86d`. |
| 61 | 2026-05-13-fit-tail-regression | brainstorm + superpowers plan | COMPLETE | Right-edge Majorana xi regression fix: `forward_tail` direction tracking, split regression loop, fallback for zero-density forward tails, backward domain_extent. Commits `a2ea66a` (test) + `30c6895` (fix). |
| 62 | 2026-05-14-nan-guard-polish | — | COMPLETE | BdG delta_0 NaN guard, is_z2_transition test cleanup for 2-arg signature. Commit `6223273`. |
| 63 | 2026-05-16-figure-validation-fixes | brainstorm + plan | COMPLETE | All 18 requirements (R1-R18) implemented. Fortran optics (bulk isotropy, BZ normalization, ISBT), Python scripts (exe routing, column mapping, width sweep), lecture text alignment. 74 PNGs regenerated. Commit `0f5adc7`. |
| 64 | 2026-05-17-kdotpy-cross-validation | brainstorm + plan | COMPLETE | 12/12 pipeline tests passing. 10 implementation units delivered (some restructured). Cross-code: bulk k=0, dispersion, QW subbands/dispersion/convergence, strain bandedge. Analytical: Zeeman, Landau, wire, g-factor, strain QW, SC QW. Commits `33859cc` + `0f5adc7`. Remaining items are extensions (Berry, Chern, BHZ, wire cross-code). |
| 65 | 2026-05-19-richardson-observables-convergence | brainstorm + plan | COMPLETE | All 18 requirements (R1-R18) delivered. convergence_helpers.py (416 lines), 6 ctest tests under `convergence` label, JSON diagnostics. Observables: CB1, subband spacing, effective mass, g-factor, absorption edge, Fermi level, subband shift, charge integral, exciton binding energy. Commit `0f5adc7`. |
| 66 | 2026-05-19-fix-pr14-review-issues | plan | COMPLETE | All 12 tasks (T1-T12) verified. Doc comments, empty-array comparison, skipped-line warnings, kdotpy wrapping, run_all.py truncation, parse_gfactor, DRY strain runners, tensile strain config, physics parameter verification, return None replacement, SC tolerance, CLAUDE.md update. Commit `33859cc`. |
| 67 | 2026-05-27-input-parser-toml-refactor | plan | COMPLETE | All 8 phases complete. toml-f submodule, 7 sub-types in defs.f90, 956-line parser (zero BACKSPACE), all consumers updated, 88 TOML configs, test infrastructure updated. Commit `0f5adc7`. |
| 68 | 2026-05-07-test-coverage-completion | spec | COMPLETE | All 3 tasks (T1-T3) done. gfactor regression (3 golden + 2 standard-star), integration tests (wire hexagon, wire strain, SC wire), contiguous attributes (all 3 sites). Phase 14. |
| 69 | 2026-05-08-core-kp-validation-ideation | ideation | PARTIAL | 6 of 8 ideas fully implemented (cross-code, Richardson, exciton, SC loop, TRK sum rule, Hermiticity). 2 partially done (verification ladder rungs 7-8, wire strain quantitative). Remaining items tracked in BACKLOG Phase 19. |
| 70 | 2026-05-30-backlog-deferred-review-items | backlog doc | ACTIVE | 2 of 9 deferred items completed (I17 wire optics unit test, I18 hamiltonian_blocks unit test). 2 partially done (I19 convergence expansion, I20 ISBT hardening). 5 not done (I1, I2, I3, I4, I13-ext). Tracked in BACKLOG Phase 18. |
| 71 | 2026-05-31-deep-physics-validation-phase20 | implementation | COMPLETE | All 4 items delivered: TRK sum rule (7 pFUnit tests), Hermiticity (20 tests across 6 files), verification ladder rungs 7-8 (g-factor + optical), wire strain quantitative (4 checks). Commits `b8f16a2` through `3596737`. |

---

## Detailed Findings

### 4. 2026-03-29-gfactor-validation

**Status: COMPLETE**

All code-level implementation is done. All 5 configs now have automated test coverage:

| Config | Golden data | CTest | Standard-star coverage |
|--------|-------------|-------|----------------------|
| `gfactor_bulk_gaas_cb.cfg` | No | No | Yes — `verify_star_gaas_bulk.py` (Roth analytical) |
| `gfactor_bulk_gaasw_cb.cfg` | Yes | Yes | Yes — golden-data regression (Phase 14) |
| `gfactor_bulk_inasw_cb.cfg` | No | No | Yes — `verify_star_inas_bulk.py` (Roth analytical) |
| `gfactor_bulk_gaas_vb.cfg` | Yes | Yes | Yes — golden-data regression (Phase 14) |
| `gfactor_qw_vb.cfg` | Yes | Yes | Yes — golden-data regression (Phase 14) |

Pre-existing `gfactor_qw_cb.cfg` retains full regression coverage.

### 26. 2026-04-21-docs-physics-revamp

**Status: COMPLETE**

Design doc archived (defines scope only). Plan had 12 tasks; all completed or superseded by lecture-test pairs, verification ladder, and standard-star benchmarks:

- **Tasks 1-2 (benchmark matrix, figure provenance):** DONE. Archived as #25, #27, #28.
- **Tasks 3-5 (figure fixes):** Lecture 05 validates g-factors vs Roth formula + Landau levels + Zeeman fan. Lecture 10 validates Stark shift with overlay plot. Lecture 03 validates wavefunction band character + normalization. All figure references corrected (`7b80237`).
- **Tasks 6-7 (solver repair):** ISBT dipole sign addressed: commutator-based velocity with cross-validation (`88840ce`, rel_diff=0.007%). Gain quasi-Fermi fully implemented with bisection solver, separate mu_e/mu_h, gain_reset(). Bulk/QW/strain benchmarks covered by verification ladder rungs 1-3, 5-6 and standard-stars S1-S3.
- **Tasks 8-9 (wire hardening):** FEAST window fix (`f9a3f6f`). Dense-sparse consistency verified in lecture 08 and verification rung 4 (R14-R16). CSR structural testing (Phase 8). Wire strain validated (Phase 13).
- **Tasks 10-12 (chapter rebuilds, verification):** All 14 chapters have Verification sections (`49ff0d3`). All 74+ figures validated (`79688ab`). Validation summaries created as #32, #33.

### 37. 2026-04-26-modern-fortran-migration

**Status: COMPLETE**

`validate_simulation_config` is DONE — implemented as `simulation_config_validate` in `defs.f90:518-586` (not `input_parser.f90` as originally planned, but functionally equivalent as a type-bound method). `config_validation_result` helper type still missing (low priority — functional need covered).

Contiguous attribute gaps — all closed in Phase 14:
- ~~ZB8bandQW profile/kpterms~~: CLOSED (already had `contiguous`)
- ~~`utils.f90 dns(:,:)`~~: CLOSED (Phase 14)
- ~~`spin_projection.f90 psi(:)`~~: CLOSED (Phase 14)
- ~~`hamiltonianConstructor.f90` velocity arrays~~: CLOSED (Phase 14)

### 54. 2026-05-08-8band-verification-ladder

**Status: COMPLETE**

4 verification scripts (rung1-rung4), 14 config files, benchmarks.md Section 7, solutions doc. CTest label `verification`. All 19 requirements implemented (some revised with documented justification: R7/R9 tolerance relaxation, R11 dropped, R10 uses self-consistency instead of Bastard).

### 55. 2026-05-09-csr-structure-testing

**Status: COMPLETE**

U1-U5 complete (invariant fixture, 10 CSR operation tests, 4 Krylov snapshot tests, refactored existing tests, regeneration infrastructure). U6 complete: all 7 of 7 Krylov paths done. Remaining 3 paths added in Phase 15:
- SC loop snapshot
- Optics wire mode snapshot
- G-factor wire mode snapshot

### 56. 2026-05-09-standard-star-benchmarks

**Status: COMPLETE**

All 7 scripts implemented and registered under ctest labels `standard-star;verification`. 24 code-review findings addressed (P0 ctest label fix, P1 S5 overlap check). Phase 15 tightening:
- S4 absorption onset: now uses regression reference value
- S7 wire g-factor: tightened with reference value
- Run wrappers centralized in S3/S5/S6 via `star_helpers.run_exe`
- `benchmarks.md` Roth formula documentation error corrected

### 59. 2026-05-11-strain-validation

**Status: COMPLETE**

All 3 verification scripts implemented and CTest registered. Shared infrastructure extracted. Code-review fix pass applied. Planning docs committed in Phase 14. Remaining items have documented justification:
- QW tests use qualitative checks instead of Bastard analytical (justified for narrow 20A well)
- Wire R7 tolerance relaxed from 5-10% to 60% (justified by free-surface relaxation)
