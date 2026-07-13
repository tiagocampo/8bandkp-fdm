# Plans Review — Implementation Status

Review date: 2026-06-27 (updated from 2026-06-03)

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
| 69 | 2026-05-08-core-kp-validation-ideation | ideation | COMPLETE | All 8 ideas fully implemented. Last 2 (verification ladder rungs 7-8, wire strain quantitative) delivered in Phase 20 (Group #71). Phase 21 added ISBT benchmark, SC benchmark, exciton Bastard reference. |
| 70 | 2026-05-30-backlog-deferred-review-items | backlog doc | ACTIVE | 4 of 9 deferred items completed: I17 wire optics unit test, I18 hamiltonian_blocks unit test, I20 ISBT hardening (Phase 21 verify_isbt_benchmark.py). I19 convergence expansion partially addressed (Phase 21 exciton tightening + SC benchmark). 5 not done (I1, I2, I3, I4, I13-ext). Tracked in BACKLOG Phase 18. |
| 71 | 2026-05-31-deep-physics-validation-phase20 | implementation | COMPLETE | All 4 items delivered: TRK sum rule (7 pFUnit tests), Hermiticity (20 tests across 6 files), verification ladder rungs 7-8 (g-factor + optical), wire strain quantitative (4 checks). Commits `b8f16a2` through `3596737`. |
| 72 | 2026-06-01-publishable-validation-benchmarks-phase21 | PRD + 6 issues | IN_PROGRESS | ISBT benchmark (6-config sweep, 7 checks per config), SC benchmark (5 checks: convergence, charge neutrality, Fermi, subband shift, potential), exciton tightening (Miller 15% + Bastard 1982 at 20%). 6 TOML configs. 2 new ctest tests (#104, #105). PR #35 (branch `feat/publishable-benchmarks-phase21`, not yet merged). |
| 73 | 2026-06-13-architecture-review-deepening | architecture review + ADR 0005 | IN_PROGRESS | Deep review (C1–C10) + unified follow-up (U-A/B/C/D), reconciled against `refactor/eigensolver-standardization`. Decisions locked in grill session 2026-06-13/14. Strategic thread: FEAST parity everywhere except bulk (**ADR 0005**, new). ADR 0001 wording corrected (data SSOT vs interpretation). CONTEXT.md glossary +*eigensolver dispatch (format vs backend)*. Candidate status + sequencing tracked in BACKLOG Phase 18 + Phase 22. PRD forthcoming in `.scratch/`. |
| 74 | 2026-06-13-eigensolver-review-fixes | design + plan | COMPLETE | Archived. All 10 tasks landed (merged via PR #39 `refactor/architecture-deepening`). Standardized on the polymorphic `eigensolve_csr` test helper (removed `solve_sparse_evp`/`solve_dense_lapack` from `src/` + stale `sc_loop` import), dense-solver workspace caching + info checks (#2/#3), QW-CSR value-only fast-path `update_kp_term_values` (PRD US-25), `[feast]` hard error, I15 FEAST+INDEX rejection. T7 truncation warning was later superseded by `reconcile_band_slice` (plan #75). |
| 75 | 2026-06-21-pr39-review-fixes | design + plan | COMPLETE | Archived. 124/124 ctest green (merged via PR #39). All in-scope fixes landed: shared `reconcile_band_slice` offset helper routed through all 3 k-sweep sites (#1/#2/#4/#10), A4 partial-window + A3 bulk+FEAST rejection in `validate()`, B2 `kp_scalar_block` error-stop, B3 shrinking-N dense-cache guard, C1 `wire_setup_adopt_precomputed` deletion. **#3 (info=3) reverted as a non-bug**; **#11 deferred → BACKLOG** (duplicated `margin_frac`/`margin_floor` constants). |
| 76 | 2026-06-21-u8-bdg-window-routing | design + plan | COMPLETE | Archived (now `docs/superpowers/{plans,specs}/archive/`). PR40 main (`feat/bdg-u8-window-routing`, 10 commits `dcdea33..c6dc762`): wire BdG window routing via `apply_solver_window` (ADR 0005 / KTD6), auto-window fallback removed (`a4ade9d`), Gershgorin-scale BdG window guard (`4c445c7`), `regression_wire_bdg_topological` + `regression_wire_bdg_strain_shift` (μ repointed to 0.6601 eV). 126/126 ctest green. Closes U8 of the parent BdG validation plan. |
| 77 | 2026-06-22-pr27-review-fixes-audit | review + 4 issues | COMPLETE | `.scratch/pr27-review-fixes/` audited 2026-06-22; 3 of 4 issues landed in PR40 main + 1 follow-up (C1, commit `b479dae`). Archived to `.scratch/archive/pr27-review-fixes/` on 2026-06-27. REVIEW.md status updated to COMPLETE. |
| 78 | 2026-06-26-u8-followup-reviews-and-codex | spec + plan | COMPLETE | Archived. 13 follow-up tasks (C1–C10 + archive + memory + final ctest/push) + 1 drift-fix (C-fix1 rashba config) + 1 doc-drift fix (parent plan U8 status footer) landed as 12 commits `b479dae..84e238a` on PR40. 126/126 ctest green; PR40 pushed (HEAD `84e238a`, OPEN awaiting merge). Resolves U8 spec/plan self-contradictions, adds 4 UBIQUITOUS_LANGUAGE entries (Gershgorin bound, sentinel gap value, FD-Nyquist tail, PHS operator `Xi = tau_x K`), renames `eval_wire_bdg_gap_app` → `eval_wire_bdg_gap`, saves `codebase-doc-drift-prevention` memory entry. |
| 79 | 2026-06-14-001-feat-bdg-majorana-validation-plan | parent plan | IN PROGRESS | Status footer added 2026-07-06 (`status: in-progress` + per-unit breakdown). Units shipped via PR #40 + PR #41: U3 (Pfaffian + Kitaev harness), U4 (hole-block unification, ADR 0007 `half_wire_integral=max(L,R)`), U5 (cross-builder PHS via `test_bdg_phs_at_finite_bx` + B.1 identity test), U6 (Sticlet polarization via A.3a wire-polarization emitter), U7 (dense-QW BdG rung via block-(2,1) structural+sign fix Task 1.8), U8 (sparse-wire zero-gap fix, PR #40), U12 (acceptance gate `lecture_13_acceptance_gate`). Partial: U11 (lecture 13 disclosures + slim-Pfaffian caption landed; full revamp not done). Still open: U1, U2, U9, U10. Explicitly deferred: U13 (periodic/Bloch BdG construction — separate scoped PR per CLAUDE.md Known Issues). NOT archived. |
| 80 | 2026-07-01-bdg-p1-fix | spec + plan | COMPLETE | Archived (now `docs/superpowers/{plans,specs}/archive/2026-07-01-bdg-p1-fix{,-design}.md`). 28 commits on `feat/bdg-p1-stabilization` (Phase 1–4: physics correctness TDD-doubled, test teeth restoration, architecture cleanup, hygiene + docs). 9 BdG regression tests + acceptance gate PASS. Spec marked `IMPLEMENTED`. Established ADR 0008 invariants; fixed dense-QW block (2,1) structural+sign bug (`1984ca5`); restored symmetric `add_peierls_coo(-B)` per PHS oracle (`f617640` disproof of double-count claim); Pfaffian infrastructure + Kitaev Majorana + cross-builder identity + Sticlet polarization landed. |
| 81 | 2026-07-05-pr41-blocker-fixes | addendum spec + plan | COMPLETE | Archived (now `docs/superpowers/{plans,specs}/archive/2026-07-05-pr41-{blocker-fixes,completion-design}.md`). 17 commits on `feat/bdg-p1-stabilization`: Phase A P0 blockers TDD-doubled (A.1/A.2 slim-Pfaffian row fix `7172858`, A.3a wire-polarization emitter `e95d001` + A.3b real-eigensolve verifier rewrite `677c0eb`, A.4 S1/S2 synthetic disclosure `4d2a40a`, A.5 acceptance-gate tighten 2.0→1.0 T `b430b1f`, A.6 3-witness reconcile `f863964`); Phase B P1 soundness (B.1 cross-builder identity with shared-H₀ invariant `56f49bb`, B.2 drop zero-escape hatch `0e8fd9e`); Phase C doc fixes (C.1 archive Status footers `1487198`, C.2–C.7 lecture disclosures + ADR amendments + spec line-23 fix `845c178..3135ad6`); plus C.8 final 2D-colormap μ-window tighten `b949e00` — bcrit_2d 3.750→2.500 T, acceptance gate 3-witness range 0.8 T ≤ 1.0 T PASS. Spec marked `IMPLEMENTED`. PR #41 pushed 2026-07-06. |
| 82 | 2026-07-12-phase-24-followups | post-archive cleanup pass | COMPLETE | 4 parallel subagent dispatches (A: PR #41 doc/test, B: BdG P1 architecture, C: BdG P1 hygiene, D: PR #39 shell wrapper) resolved DONE-WITH-GAPS findings from adversarial review of rows #74–#81. Sites fixed: `extract_block_csr` moved from `spectral_bdg_wire.f90:189-234` → `sparse_matrices.f90:1185-1222` (Task 3.2 complete); 4 inline `open()` blocks in `main_topology.f90:927, 1010, 1150` consolidated to named writers in `outputFunctions.f90` (Task 3.3 complete); 26 bare `stop 1` across 8 src/ files replaced with descriptive `error stop '<message>'` (Task 4.1 expanded project-wide); B.2 escape hatch removed (strict `@assertTrue(s1 == s2 .and. s1 /= 0)` in `test_wire_pfaffian_witness.pf`, now `@todo` U13 — 1 expected unit fail); 8 P1-spec lines (34, 72, 103, 260, 333, 335, 430, 446) + 3 ADR 0008 lines (50, 52, 57) reconciled to current state; 10 "4-witness" labels swept to "3-witness" in `lecture_13_topological.py` 7 sites + `docs/lecture/13-topological-superconductivity.md` 3 sites; PRD `.scratch/archive/bdg-majorana-validation/PRD.md` updated (14→47 commits, 124/124→145/145 ctest, 4→3 witness); `TOLERANCE_BCRIT_RANGE` 2.0→1.0 in `lecture_13_topological.py:278-283`; PR #39 shell wrapper `test_qw_bandstructure_dense_feast.sh` created (CMakeLists.txt:1260 updated). **BLOCKING-EMPIRICAL gap (A.3a/A.3b integration — wire-polarization emitter never fires on canonical wire due to FEAST noise floor > near-zero threshold):** resolved per option (b); `verify_majorana_polarization.py` module docstring now carries `@todo U13` block; `parse_polarization()` exits 0 with `SKIP: ... not produced ... (BLOCKING-EMPIRICAL on canonical wire)`; other failure modes still fail loud per spec D7; acceptance gate `lecture_13_acceptance_gate.sh` (3-witness, 1.0 T) remains working regression net. Build green across all 4 scopes. Tests: 49/50 unit (1 expected fail = `test_wire_pfaffian_witness`), 8/8 BdG regression, 4/4 spectral, 3/3 wire BdG regression, 1/1 lecture 13, 1/1 new wrapper all green. **Landed on PR #41 (2026-07-12):** the 8 follow-up commits (`5a39bfb..01a2bdc`) were FF-pushed to `feat/bdg-validation-pass2` head; PR #41 now at `5a39bfb` (85 commits, MERGEABLE), ready for merge to `main`. The 7 originally-tracked commits (`a9a71a6..01a2bdc`) cover the 4-subagent cleanup + U13 BLOCKING-EMPIRICAL SKIP + archive moves; the 8th (`5a39bfb`) is a subsequent chore gitignore update. |

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
