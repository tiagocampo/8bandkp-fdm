# Implementation Backlog — Ordered by Priority

Consolidated from REVIEW.md on 2026-05-03.
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

## Phase 4: Optics & Documentation Figures

Generate all missing Python figures. Independent of Phases 5-6.

| Source | What | Effort |
|--------|------|--------|
| #22 | 11 optics figure functions in generate_all_figures.py (absorption, gain, ISBT, exciton, spin-resolved) | Medium |
| #50 P1-3 | Bandstructure figures (bulk, QW, wire) + wavefunctions + wire geometry | Medium |
| #26 | Figure fixes, solver repairs (ISBT dipole sign), chapter rebuilds (Ch06, Ch08) | Medium |

**Result:** Complete documentation with all figures. Groups 22, 26, 50 -> COMPLETE.

---

## Phase 5: Peierls / Landau Levels (New Confinement Mode)

New physics capability. Peierls in bulk = new confinement mode (discretize x, 8Nx8N).

| Source | What | Effort |
|--------|------|--------|
| New | Design: `confinement=3` mode (Landau) -- x-discretized, gauge choice, 8Nx8N assembly reusing QW machinery | Design |
| #47 | Implement Peierls in ZB8bandBulk -> calls confinement=3 path | High |
| #47 | Fix ExternalField float parsing in input_parser.f90 | Low |
| #48 | Landau level verification: InAs, GaAs benchmarks vs analytical E_n = hw_c(n+1/2) | Medium |
| #50 P4 | Landau fan diagram (E_n vs B) + effective mass extraction | Medium |

**Result:** Magnetic field works for all geometries. Groups 47, 48 -> COMPLETE.

---

## Phase 6: Topological Suite Completion

Most advanced features. Fu-Kane is independent of Peierls but logically last.

| Source | What | Effort |
|--------|------|--------|
| #38 | `build_bdg_hamiltonian_qw` (dense variant for QW BdG) | Medium |
| #49 | `compute_z2_fukane` for QW (Fu-Kane parity method) | Medium |
| #49 | `compute_z2_gap_sweep` for wire (phase diagram M vs width) | Medium |
| #49 | Berry curvature heatmap figures + BHZ phase transition figures | Low-Medium |
| #38 | 4 unit test files (magnetic_field, z2, edge_states, green_functions) | Medium |
| #38 | `compute_longitudinal_conductance` + `compute_spectral_function` | Medium-High |

**Result:** Full topological analysis suite. Groups 38, 49 -> COMPLETE.

---

## Summary

| Phase | Groups resolved | Net new capability | Effort |
|-------|----------------|-------------------|--------|
| 1. Testing | 4 | Full regression coverage | ~2-3 days |
| 2. Physics wiring | 1 | Bulk SC, bulk EF, gfactor SC | DONE |
| 3. PR review fixes | — | Correctness, compat, quality | DONE |
| 4. Figures | 3 | Complete documentation | ~3-4 days |
| 5. Peierls/Landau | 2 | Magnetic field for all modes | ~5-7 days |
| 6. Topological | 2 | Full topo suite | ~5-7 days |
| **Total** | **10 remaining** | | **~12-18 days** |
