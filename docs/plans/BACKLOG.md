# Implementation Backlog — Ordered by Priority

Consolidated from REVIEW.md on 2026-05-06.
Piezoelectric explicitly excluded (ZB [001] = zero by symmetry, wires = negligible).

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

**Additional fixes (root-cained from "pre-existing" failures):**
- `compute_chern_qwz`: non-Hermitian QWZ Hamiltonian, eigenvector cross-contamination, wrong FHS plaquette formula (fixed test_chern_number + regression_topology_qwz_chern)
- BHZ Z2 configs: missing `compute_hall`/`qwz_u` fields caused parser field-order misalignment (fixed regression_topology_bhz_z2)
- `mu_B = e*hbar/(2*m0)` gave 9.274e-4 instead of 5.788e-5 eV/T (16x too large); added CODATA constant (fixed regression_landau_inas)
- Rashba config: `bdg:` block inside topology block, missing `B_vec`, `ldos_E_range` consumed `ldos_num_E` (fixed regression_topology_rashba_phase)

---

## Phase 4: COMPLETED (2026-05-04)

All 74 existing figures validated, ISBT cross-verified (z-dipole vs commutator velocity, rel_diff=0.007%), 5 new physics figures added (bulk E(k), QW subbands, wavefunctions, wire geometry, Zeeman fan). Ch06 + Ch08 rebuilt with verified captions. All 14 lecture chapters have zero broken figure references. Groups #22, #26, #50 -> COMPLETE.

The historical `timing_dense_vs_sparse` timeout note is archived in `docs/plans/archive/phase4-discrepancy-log.md`; `fig_timing_dense_vs_sparse` now uses the 600s timeout identified in that log. The `qw_strained_bands` segfault was fixed (commit 9afc205: gfortran -O3 array temporaries in derived-type components).

**Plan:** `docs/plans/archive/2026-05-04-phase4-optics-figures-plan.md` (archived)
**Commits:** 12 commits from `88840ce` to `79688ab`.

---

## Phase 5: COMPLETED (2026-05-05)

New confinement=3 (Landau) mode with x-discretized 8NxN Hamiltonian, Landau gauge A=(0,Bz*x,-By*x), Zeeman splitting, orbital quantization, B-sweep fan diagram, analytical validation against E_n = E_C + hw_c(n+1/2). 58/58 tests pass.

**Plan:** `docs/plans/archive/2026-05-02-magnetic-field-landau-design.md` (archived)
**Implementation plan:** `docs/plans/archive/2026-05-04-landau-bulk-phase5.md` (archived)
**Commits:** 82 commits from `70c817b` to `8bd2387`, plus sign fix and CI wiring.

**Code review fixes applied:**
- `compute_gauge_shifts` added to `magnetic_field.f90` (was missing, blocked build)
- Pi_z sign error fixed (latent: By=0 in all configs)
- B_vec mutation in B-sweep, contiguous attribute, variable renames (kx2→piy2 etc.)
- landau_sweep validation, nB guard

**Tests added:**
- `test_landau.pf` — pFUnit unit tests (hermiticity, bulk recovery, gauge shifts with By≠0)
- `regression_landau_bulk_gaas/inas/inas_bsweep` — golden data comparison via `test_landau_bulk.sh`
- `regression_landau_analytical` — Python verification against analytical Landau levels
- SC golden data regenerated (stale from Zeeman/mu_B/sc_potential_shift changes)

**Result:** Magnetic field works for all geometries. Groups 47, 48 -> COMPLETE.

---

## Phase 6: COMPLETED (2026-05-06)

Full topological suite completion repair delivered and pushed.

**Plan:** `docs/plans/archive/2026-05-05-phase6-topological-suite.md` (archived)
**Repair plan:** `docs/plans/archive/2026-05-05-phase6-completion-repair.md` (archived)
**Commits:** `c56fbd4` through `20c3f19`.

**Delivered:**
- Dense QW BdG path and app dispatch.
- QW Fu-Kane parity invariant using band-major inversion and Kramers-pair parity sign.
- Berry/Kubo conductance and Landauer helper.
- Bulk/QW/wire spectral functions and LDOS shifted-matrix fix.
- Real gap sweep evaluators: BHZ analytic, QW Fu-Kane, wire BdG.
- Phase 6 parser fields, output fields, docs, and regression wiring.
- New/expanded regressions for QW Fu-Kane, conductance, spectral bulk/QW/wire, and sweep BHZ/QW.

**Verification:** fresh configure/build passed; full suite `66/66` passed; manual topologicalAnalysis smoke configs passed.

**Result:** Full topological analysis suite. Groups #38, #49, #51 -> COMPLETE.

---

## Remaining Backlog

Only non-Phase-6 items remain from the review:

| Source | What | Effort |
|--------|------|--------|
| #4 | 5 gfactor regression tests (bulk CB/VB GaAs, GaAsW, InAsW, QW VB) | Low |
| #37 | `validate_simulation_config` + close `contiguous` gaps (3 known sites) | Medium |
| #8 | 3 integration tests (wire hexagon, wire strain, SC wire) + 2 regression datasets | Medium |
| #26 | Docs physics revamp remaining tasks 3-12 | Medium-High |

---

## Summary

| Phase | Groups resolved | Net new capability | Effort |
|-------|----------------|-------------------|--------|
| 1. Testing | 4 | Full regression coverage | ~2-3 days |
| 2. Physics wiring | 1 | Bulk SC, bulk EF, gfactor SC | DONE |
| 3. PR review fixes | — | Correctness, compat, quality | DONE |
| 4. Figures | 3 | Complete documentation | DONE |
| 5. Peierls/Landau | 2 | Magnetic field for all modes | DONE |
| 6. Topological | 3 | Full topo suite | DONE |
| **Remaining** | **5 groups** | Non-Phase-6 backlog above | TBD |
