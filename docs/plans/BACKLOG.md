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

## Phase 3: Optics & Documentation Figures

Generate all missing Python figures. Independent of Phases 4-5.

| Source | What | Effort |
|--------|------|--------|
| #22 | 11 optics figure functions in generate_all_figures.py (absorption, gain, ISBT, exciton, spin-resolved) | Medium |
| #50 P1-3 | Bandstructure figures (bulk, QW, wire) + wavefunctions + wire geometry | Medium |
| #26 | Figure fixes, solver repairs (ISBT dipole sign), chapter rebuilds (Ch06, Ch08) | Medium |

**Result:** Complete documentation with all figures. Groups 22, 26, 50 -> COMPLETE.

---

## Phase 4: Peierls / Landau Levels (New Confinement Mode)

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

## Phase 5: Topological Suite Completion

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
| 3. Figures | 3 | Complete documentation | ~3-4 days |
| 4. Peierls/Landau | 2 | Magnetic field for all modes | ~5-7 days |
| 5. Topological | 2 | Full topo suite | ~5-7 days |
| **Total** | **10 remaining** | | **~15-21 days** |
