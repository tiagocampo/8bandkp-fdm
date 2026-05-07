# Plans Review — Implementation Status

Review date: 2026-05-06

| # | Group | Files | Status | Notes |
|---|-------|-------|--------|-------|
| 1 | 2026-03-28-codebase-improvement | design | COMPLETE | Archived |
| 2 | 2026-03-29-cmake-testing | design | COMPLETE | Archived |
| 3 | 2026-03-29-gfactor-improvements | design | COMPLETE | Archived |
| 4 | 2026-03-29-gfactor-validation | standalone | INCOMPLETE | Missing regression tests for 5/6 configs |
| 5 | 2026-03-29-pr-review-fixes | design | COMPLETE | Archived (Makefile target + z(1) guard added) |
| 6 | 2026-03-29-self-consistent-sp | design | COMPLETE | Archived (bulk EF shift, delta-doping, gfactor SC wired; bulk SC uses QW path with confDir=z) |
| 7 | 2026-03-30-sc-validation-report | standalone | COMPLETE | Archived |
| 8 | 2026-04-02-quantum-wire | design | COMPLETE | Archived (piezoelectric excluded: ZB [001]=zero by symmetry; FFTW Poisson not needed, PARDISO used) |
| 9 | 2026-04-03-linalg-backend | design | COMPLETE | Archived |
| 10 | 2026-04-04-documentation-overhaul | design | COMPLETE | Archived (wire_optical_spectrum figure added; 7 renamed figures have better names) |
| 11 | 2026-04-04-pr-review-fixes | design | COMPLETE | Archived |
| 12 | 2026-04-06-bulk-chapter-overhaul | design | COMPLETE | Archived |
| 13 | 2026-04-06-documentation-deep-rewrite | standalone | COMPLETE | Archived |
| 14 | 2026-04-06-quickstart-fixes | design | COMPLETE | Archived |
| 15 | 2026-04-07-strain-unification | design + plan | COMPLETE | Archived |
| 16 | 2026-04-12-qw-documentation-overhaul | design + 2 plans | COMPLETE | Archived (anticrossing prose updated with 34 meV gap from realistic InAsW(8nm)/GaSbW(8nm) config) |
| 17 | 2026-04-17-fd-variable-coeff-fix | design + plan | COMPLETE | Archived (staggered-grid deviation justified) |
| 18 | 2026-04-18-ch01-03-review-fix | design + plan | COMPLETE | Archived |
| 19 | 2026-04-18-ch04-06-review-fix | plan | COMPLETE | Archived |
| 20 | 2026-04-18-ch04-strain-figures | design + plan | COMPLETE | Archived |
| 21 | 2026-04-19-ch03-ch04-figures | plan | COMPLETE | Archived |
| 22 | 2026-04-19-qw-tutorials | design | COMPLETE | All optics figure functions + 5 new Group #50 figures. ISBT cross-validated |
| 23 | 2026-04-20-remaining-work | design + plan | INCOMPLETE | Task 8: core-shell benchmark config missing |
| 24 | 2026-04-20-strain-sign-and-buffer-fix | standalone | COMPLETE | Archived |
| 25 | 2026-04-21-benchmark-matrix | standalone | COMPLETE | Archived (reference doc) |
| 26 | 2026-04-21-docs-physics-revamp | design + plan | INCOMPLETE | 10/12 tasks missing, design archived |
| 27 | 2026-04-21-failure-ledger | standalone | COMPLETE | Archived (reference doc) |
| 28 | 2026-04-21-figure-provenance | standalone | COMPLETE | Archived (reference doc) |
| 29 | 2026-04-21-nextnano-bulk-qw-comparison | standalone | COMPLETE | Archived (reference doc) |
| 30 | 2026-04-21-nextnano-qcse-optics-comparison | standalone | COMPLETE | Archived (reference doc) |
| 31 | 2026-04-21-nextnano-wire-comparison | standalone | COMPLETE | Archived (reference doc, 5 follow-ups pending) |
| 32 | 2026-04-21-validation-summary | standalone | COMPLETE | Archived (reference doc) |
| 33 | 2026-04-21-wire-validation-summary | standalone | COMPLETE | Archived (reference doc) |
| 34 | 2026-04-25-optical-properties | redesign + implementation | COMPLETE | Archived |
| 35 | 2026-04-25-wire-gfactor-commutator | design | COMPLETE | Archived |
| 36 | 2026-04-26-hamiltonian-performance-refactor | design + plan | COMPLETE | Archived |
| 37 | 2026-04-26-modern-fortran-migration | design + plan | MOSTLY COMPLETE | validate_simulation_config missing |
| 38 | 2026-04-27-bdg-topological-sc | design + impl plan | COMPLETE | Archived; Phase 6 completion repair delivered QW BdG, Fu-Kane QW, conductance, spectral, LDOS, and sweep paths |
| 39 | 2026-04-27-pr11-review-fixes | design + plan | COMPLETE | Archived |
| 40 | 2026-04-27-wire-fastpath-followup | doc + plan | COMPLETE | Archived |
| 41 | 2026-04-28-modern-fortran-phase2 | plan | COMPLETE | Archived |
| 42 | 2026-04-28-pr11-cache-safety-fixes | plan | COMPLETE | Archived |
| 43 | 2026-04-29-modern-fortran-remaining | plan | COMPLETE | Archived |
| 44 | 2026-04-29-phase4-exploratory | standalone | COMPLETE | Archived |
| 45 | 2026-04-29-pr12-fixes | plan | COMPLETE | Archived (zdotc decision documented in CLAUDE.md) |
| 46 | 2026-05-01-bhz-remaining-work | plan | COMPLETE | Archived |
| 47 | 2026-05-02-magnetic-field-landau | design + impl plan | COMPLETE | Archived; confinement=3 Landau mode, gauge shifts, fan diagram, analytical validation |
| 48 | 2026-05-02-physics-figures-implementation | plan | COMPLETE | Landau verification, input-reference updated, CI wired |
| 49 | 2026-05-02-topological-suite-verification | spec + plan | COMPLETE | Archived; Phase 6 completion repair delivered Fu-Kane QW, gap sweep, conductance/spectral regressions; Berry/phase figures already present |
| 51 | 2026-05-05-phase6-completion-repair | plan | COMPLETE | Archived; 66/66 tests passed; pushed through commit 20c3f19 |
| 50 | 2026-05-03-physics-figures-extended | plan | COMPLETE | All 5 phases done: bulk E(k), QW subbands, wavefunctions, wire geometry, Zeeman fan |

---

## Detailed Findings

(Filled in as each group is reviewed)

### 4. 2026-03-29-gfactor-validation

**Status: INCOMPLETE**

All code-level implementation is done (g-factor functions, sigma matrices, self-term exclusion, near-zero guards). Missing:

- **5/6 test configs lack reference data and automated regression tests:**
  - `gfactor_bulk_gaas_cb.cfg` — no reference data, no CMake test
  - `gfactor_bulk_gaasw_cb.cfg` — no reference data, no CMake test
  - `gfactor_bulk_inasw_cb.cfg` — no reference data, no CMake test
  - `gfactor_bulk_gaas_vb.cfg` — no reference data, no CMake test
  - `gfactor_qw_vb.cfg` — no reference data, no CMake test
- Only `gfactor_qw_cb.cfg` has full regression coverage
- Document values are stale (QW CB g-factor values differ from current reference data)

### 6. 2026-03-29-self-consistent-sp

**Status: COMPLETE**

All items resolved in PR #13 (Phase 2 physics wiring):
- **Bulk EF shift**: Uniform diagonal shift added to `ZB8bandBulk` via `cfg%Evalue` and `cfg%sc_potential_shift` (TDD test: `test_bulk_ef_shift`)
- **Delta-doping**: Extended `doping_spec` with Gaussian profile (`delta<N>: NS FWHM POS` syntax), implemented in `build_doping_charge` (TDD test: `test_delta_doping_gaussian`)
- **gfactor SC**: `self_consistent_loop` (QW) and `self_consistent_loop_wire` (wire) wired into `main_gfactor.f90`
- **Bulk SC**: Uses QW path (`confDir='z'`, single material + doping) — no separate bulk SC routine needed. Warning printed if user tries `confDir='n'` with SC enabled.
- **b_field parsing**: Bulk Landau level magnetic field support with peek/backspace pattern

### 22. 2026-04-19-qw-tutorials

**Status: COMPLETE**

Fortran code fully done. All 15 optics figure functions implemented in generate_all_figures.py. All PNG files generated and verified in docs/figures/. Verified on 2026-05-04: absorption edge positions, ISBT peak energies, gain blue-shift behavior all physically correct.

ISBT cross-validation: z-dipole vs commutator velocity consistency confirmed (see test_isbt_dipole_velocity_consistency in test_optical_qw.pf, rel_diff = 0.007%).

### 26. 2026-04-21-docs-physics-revamp

**Status: INCOMPLETE**

Design doc archived (defines scope only). Plan has 12 tasks; Tasks 1-2 (control docs, provenance) DONE. Missing Tasks 3-12:

- **Tasks 3-5**: Figure fixes (broken refs, missing optics/gain/ISBT figures, PNG generation)
- **Tasks 6-7**: Solver repair (ISBT dipole sign, gain quasi-Fermi integration)
- **Tasks 8-9**: Wire hardening (FEAST convergence, branch tracking edge cases)
- **Tasks 10-11**: Chapter rebuilds (Ch06 optics, Ch08 wire with real data)
- **Task 12**: Final verification sweep

### 23. 2026-04-20-remaining-work

**Status: INCOMPLETE**

13/14 tasks done. Missing:

- **Task 8 (Stier & Bimberg core-shell wire benchmark)**: Integration test `test_wire_core_shell.sh` exists but regression config `wire_inas_gaas_core_shell.cfg` was never created — test cannot function without it

### 37. 2026-04-26-modern-fortran-migration

**Status: MOSTLY COMPLETE**

Phases 1-2 done (forall→do concurrent, goto→named blocks, dsqrt→sqrt, external→interfaces, elemental pure, private default, do concurrent in hot paths, finalizers, type-bound methods). Missing:

- **`validate_simulation_config` subroutine** (Plan Task 10): Post-parse validation routine not implemented in input_parser.f90
- **`config_validation_result` helper type**: Not present in defs.f90
- **Known `contiguous` gaps** documented in CLAUDE.md: ZB8bandQW profile/kpterms, utils.f90 dns(:,:), spin_projection.f90 psi(:)

### 38. 2026-04-27-bdg-topological-sc

**Status: COMPLETE**

**Archive:** `docs/plans/archive/2026-04-27-bdg-topological-superconductivity-design.md`, `docs/plans/archive/2026-04-29-bdg-topological-sc-implementation-plan.md`

Core physics and Phase 6 repair items are complete:

- Dense QW BdG assembly and QW app path implemented.
- QW Fu-Kane parity invariant implemented with band-major inversion operator and pair-subspace parity sign.
- Berry/Kubo conductance and Landauer helper implemented.
- QW, bulk, and wire spectral functions implemented; LDOS shifted-matrix bug fixed.
- Gap sweep now uses real evaluators: module-level BHZ analytic, app-level QW Fu-Kane, and app-level wire BdG.
- Unit and regression coverage added for BdG, Z2/Fu-Kane, Berry/Kubo, Green/LDOS/Landauer, parser, spectral, conductance, and sweep paths.

Verification on 2026-05-06: fresh configure/build passed; full suite `66/66` passed; manual topologicalAnalysis smoke configs passed.

### 47. 2026-05-02-magnetic-field-landau

**Status: COMPLETE**

**Archive:** `docs/plans/archive/2026-05-02-magnetic-field-landau-design.md`, `docs/plans/archive/2026-05-02-magnetic-field-landau-implementation-plan.md`, `docs/plans/archive/2026-05-04-landau-bulk-phase5.md`

All phases done. New confinement=3 (Landau) mode replaces the original "Peierls in ZB8bandBulk" approach with a proper x-discretized 8NxN Hamiltonian. Landau gauge A=(0,Bz*x,-By*x) via `compute_gauge_shifts`. B-sweep fan diagram, analytical validation (E_n = E_C + hw_c(n+1/2)), ky-degeneracy check. `b_field` config workaround used for ExternalField (float parsing deferred — no user impact). Fu-Kane Z2 deferred to Phase 6.

**Commits:** 82 commits on feature/bdg-topological-superconductivity.

### 48. 2026-05-02-physics-figures-implementation

**Status: COMPLETE**

All tasks done. Landau level verification via analytical formula and B-sweep fan diagram. input-reference.md updated for confinement=3. Landau configs registered as regression tests.

### 49. 2026-05-02-topological-suite-verification

**Status: COMPLETE**

**Archive:** `docs/plans/archive/2026-05-02-topological-suite-verification-spec.md`, `docs/plans/archive/2026-05-02-topological-suite-verification-plan.md`

Phase 6 completion repair closed the remaining verification gaps:

- Fu-Kane Z2 for QW implemented and regression-tested (`regression_topology_qw_fukane_z2`).
- Gap sweep implemented with BHZ analytic, QW Fu-Kane, and wire BdG executable paths; sweep regressions registered.
- Conductance and spectral modes registered as regressions.
- Berry curvature / phase-diagram / edge / spectral figure support was delivered earlier in the branch.

Verification on 2026-05-06: `ctest --test-dir build -j4 --output-on-failure` passed `66/66`.

### 51. 2026-05-05-phase6-completion-repair

**Status: COMPLETE**

**Archive:** `docs/plans/archive/2026-05-05-phase6-completion-repair.md`

Executed task-by-task with review gates. Final branch state pushed to `origin/feature/bdg-topological-superconductivity` through:

- `c56fbd4 docs(topo): design phase 6 completion repair`
- `ad63187 fix(topo): repair phase 6 foundation defects`
- `bd3b98f feat(topo): complete Phase 6 topological suite wiring`
- `b35a979 feat(topo): add Phase 6 parser fields and regression configs`
- `889109d test(topo): add Phase 6 integration test scripts`
- `20c3f19 fix(topo): resolve phase 6 verification issues`

Final verification:

- Fresh configure/build passed.
- Full suite passed: `66/66`.
- Manual smoke passed for conductance QWZ, QW Fu-Kane, QW BdG, spectral QW/wire/bulk, sweep BHZ, and sweep QW.

### 50. 2026-05-03-physics-figures-extended

**Status: COMPLETE**

All 5 phases delivered in Phase 4 (commit 79688ab):
- Phase 1: Bulk E(k) for GaAs, InAs, InSb (3-panel figure)
- Phase 2: QW wavefunctions with band edge profile
- Phase 3: Wire material map + radial potential cut
- Phase 4: Analytical Zeeman fan diagram (spin splitting vs B)
- Phase 5: Figures integrated into `generate_all_figures.py` (not separate script)
