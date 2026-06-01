---
date: 2026-05-19
topic: richardson-observables
---

# Richardson Extrapolation on Physics Observables

## Summary

A convergence validation fixture that runs standard-star systems, SC loop benchmarks, and exciton benchmarks at multiple grid spacings and FD orders, extracts physics observables (subband energy, effective mass, g-factor, absorption edge, SC physics outputs, exciton binding energy), performs Richardson extrapolation with Grid Convergence Index (GCI) uncertainty quantification, and asserts convergence rates match theoretical predictions. Runs as an on-demand ctest label.

---

## Problem Frame

The codebase supports FD orders 2-10 and arbitrary grid refinement, but no test validates that physics observables converge at the theoretically expected rate. The existing Richardson extrapolation fixture covers only FD stencil coefficients at the operator level — a silently wrong matrix assembly (like the order-10 D1 coefficient bug) would produce plausible-looking eigenvalues that no test catches.

Standard-star benchmarks (S1-S7) report numbers at single resolutions with no discretization uncertainty. Without convergence certification, these numbers cannot distinguish "correct to N digits" from "looks right by coincidence." The SC loop test only checks that convergence happens, not that the converged physics is correct. The exciton module has zero test coverage.

Richardson extrapolation on physics observables closes these gaps: every observable gets a convergence certificate (rate + GCI), the extrapolated continuum value becomes a self-generated reference, and future regressions are caught by convergence rate deviation — a more sensitive detector than absolute-value comparison.

---

## Requirements

**Grid convergence (Richardson + GCI)**

- R1. For each supported system and observable, run the simulation at 4-5 grid spacings (fixed FD order), extract the observable at each resolution, fit E(h) = E_exact + c_2 h^2 + c_4 h^4 + ..., compute the Richardson-extrapolated continuum value and the Grid Convergence Index (Roache 1998).
- R2. Assert that the observed convergence rate (from log-log slope of consecutive grid pairs) matches the theoretical rate for the FD order used, within a tolerance calibrated per system. Material interfaces (abrupt barrier-well boundaries) degrade smooth h^p convergence; the tolerance must be empirically established for each QW system before asserting. The theoretical rate serves as a reference baseline, not a hard gate.
- R3. Report GCI with Roache safety factors: Fs=3 for two-grid Richardson, Fs=1.25 for three-grid. With 4-5 grid levels, compute GCI from the finest 3-grid subset using Fs=1.25. Each standard-star benchmark number gets a companion GCI uncertainty band.

**FD order convergence**

- R4. For each QW system (S4-S6) and observable, run the simulation at FD orders {2, 4, 6, 8, 10} (fixed grid spacing, fixing h and varying N so each order has sufficient boundary padding), extract the observable at each order, and verify monotonic convergence toward the Richardson limit. For wire (S7), restrict to orders {2, 4, 6} until grid convergence is proven stable at those orders. Bulk systems (S1-S3) are excluded from order convergence — the 8x8 Hamiltonian has no spatial discretization.
- R5. Assert that the observed convergence rate as a function of FD order is consistent with the expected acceleration (higher order => faster convergence when grid is resolved). Orders 8 and 10 are informational at QW material interfaces — monotonic convergence may not hold due to stencil boundary effects, and assertion tolerance is relaxed accordingly.

**Observable coverage**

- R6. Extract and converge subband energy (CB1, CB2, subband spacing) from `bandStructure` eigenvalue output.
- R7. Extract and converge effective mass from `bandStructure` k-sweep output via parabolic fit near k=0.
- R8. Extract and converge g-factor (gz component) from `gfactorCalculation` output.
- R9. Extract and converge absorption edge position from `opticalProperties` spectrum output. The output energy grid (num_energy_points, E_min, E_max) must be held constant across all grid spacings in the convergence sweep to isolate FD convergence from post-processing discretization.
- R10. Extract and converge SC physics outputs: Fermi level, subband energy shift relative to flat-band, and charge density integral — from `bandStructure` with SC enabled. SC convergence is empirically determined (composed of FD discretization + Poisson solver error + DIIS mixing), so R2 rate assertion does not apply to SC observables. Instead, assert monotonic convergence + Richardson extrapolation + GCI band. SC_tolerance must be set sufficiently tight (at least 10x finer than expected FD error at the finest grid) to ensure Richardson measures FD convergence, not SC iteration residual.
- R11. Extract and converge exciton binding energy from `bandStructure` or `opticalProperties` with exciton enabled. Validated against published experimental values (e.g., Miller et al., Phys. Rev. B 1985 for GaAs/AlGaAs) as the primary reference, with the 2D hydrogen model E_b ~ Rydberg/4 as a sanity check during development.
- R11a. Prerequisite smoke test for the exciton module: a single-resolution test that verifies the binding energy is physically plausible (e.g., 1-20 meV for a GaAs/AlGaAs QW) and the module runs without errors. This must pass before R11 convergence work begins.

**System coverage**

- R12. Grid and order convergence tests for QW and wire systems: S4 (GaAs/AlGaAs QW), S5 (InAs/GaSb QW), S6 (InAs/GaAs strained QW), S7 (InAs wire). Bulk systems (S1-S3) are excluded from convergence testing — the 8x8 dense Hamiltonian has no spatial grid, so Richardson extrapolation on h is undefined.
- R13. SC loop convergence benchmark using the existing GaAs/AlGaAs SC configuration, at multiple grid resolutions.
- R14. Exciton convergence benchmark using a GaAs/AlGaAs QW configuration with exciton enabled, at multiple grid resolutions and FD orders.

**Test infrastructure and reporting**

- R15. All convergence tests run under a dedicated ctest label (`convergence`), on-demand only, not gating CI. Convergence tests carry ONLY the `convergence` label — not `standard-star` or `verification` — to avoid bloating existing label runs.
- R16. Convergence results (Richardson-extrapolated values, GCI, observed convergence rates, pass/fail status) written to JSON per system per observable.
- R17. Shared convergence helper module provides reusable Richardson fitting (variable-order, not hardcoded p=2), GCI computation, convergence rate extraction, and observable extraction. The module should be extracted from the first few convergence test scripts once the actual code duplication pattern is clear, not designed upfront. Existing Richardson implementations in `validation/qw/test_qw_convergence.py` and `scripts/lecture_11_convergence.py` provide starting points to extend with GCI support.
- R18. If a convergence test fails (wrong rate, non-monotonic convergence, GCI exceeds a threshold), diagnostics identify which grid level or FD order caused the failure and the specific observable value at that point.

---

## Acceptance Examples

- AE1. **Covers R1, R2, R6, R12.** Given the S4 (GaAs/AlGaAs QW) config at FDorder=2 with 4 grid spacings, when the convergence fixture runs, it extracts CB1 energy at each resolution, computes Richardson extrapolation, reports GCI, and asserts the convergence rate matches h^2 scaling within tolerance.
- AE2. **Covers R4, R5, R6, R12.** Given the S4 config at a fixed FDstep=201 with FD orders {2, 4, 6, 8, 10}, when the order sweep runs, CB1 energy converges monotonically toward the Richardson limit and the convergence rate accelerates with higher FD order.
- AE3. **Covers R10, R13.** Given the SC benchmark config at 4 grid resolutions, when the SC convergence fixture runs, it extracts Fermi level, subband shift, and charge integral at each resolution, performs Richardson extrapolation, and asserts monotonic convergence with a GCI band. SC rate is empirically measured, not compared to a theoretical FD order.
- AE4. **Covers R11, R14.** Given a GaAs/AlGaAs QW config with exciton enabled at 4 grid resolutions, when the exciton convergence fixture runs, it extracts exciton binding energy, performs Richardson extrapolation, and the extrapolated value is checked against published experimental values first, falls back to 2D hydrogen model comparison if no published value matches the config, and always verifies internal consistency of the Richardson-extrapolated reference.
- AE5. **Covers R15, R18.** When `ctest -L convergence` is run, all convergence tests execute and produce JSON results. If the S7 wire convergence test fails due to non-monotonic CB1 energy, the diagnostics report which grid level diverged and the CB1 values at each level.

---

## Success Criteria

- Every QW and wire standard-star observable (subband energy, effective mass, g-factor, absorption edge) has a convergence certificate: observed rate, Richardson-extrapolated value, and GCI uncertainty band.
- SC loop physics outputs (Fermi level, subband shift, charge integral) converge monotonically across grid resolutions with GCI bands — the first quantitative physics validation of the self-consistent solver.
- Exciton binding energy converges and the extrapolated value matches published experimental values — the first automated test of the exciton module (after smoke test R11a passes).
- A downstream implementer can run `ctest -L convergence` and get actionable pass/fail results with JSON diagnostics for every system and observable.

---

## Scope Boundaries

- CI gating — convergence tests run on demand only, never blocking merges
- BdG / topological invariant convergence — separate physics modules with different convergence semantics
- Cross-code convergence — already covered by `validation/qw/test_qw_convergence.py`
- Automated reference update from Richardson-extrapolated values — the extrapolated values inform future reference setting but no automatic golden file update
- FD operator-level convergence — already covered by `tests/unit/test_fd_convergence.pf`
- Bulk system convergence — bulk (S1-S3) uses 8x8 dense Hamiltonian with no spatial discretization; grid convergence is undefined
- New physics observables beyond the six listed (subband energy, effective mass, g-factor, absorption edge, SC outputs, exciton binding energy)

---

## Key Decisions

- **Both grid convergence and order convergence axes:** Grid convergence quantifies discretization error (GCI); order convergence verifies FD stencil correctness at the physics level. Both are needed for full convergence certification.
- **QW and wire only for grid convergence:** Bulk systems (S1-S3) have no spatial FD grid (8x8 dense Hamiltonian), so Richardson extrapolation on h is undefined. Convergence testing is restricted to QW (S4-S6) and wire (S7) systems.
- **FD orders 2-10 for QW, 2-6 for wire:** Higher FD orders on wire CSR are untested territory. Wire order sweep is restricted until grid convergence proves stable. QW systems use all five orders.
- **On-demand ctest label:** Multi-resolution sweeps are 10-25x slower than single-resolution runs. An on-demand label makes convergence data available without slowing CI.
- **Exciton smoke test before convergence:** The exciton module has zero test coverage. A prerequisite smoke test (R11a) verifies basic correctness before convergence certification (R11).
- **Published experimental values as primary exciton reference:** Miller et al. (Phys. Rev. B 1985) for GaAs/AlGaAs, not a three-layer chain. The 2D hydrogen model serves as a development sanity check, not a test target.
- **Shared helper extracted from first tests:** The convergence helper module (R17) emerges from the first few test scripts rather than being designed upfront. Existing Richardson implementations in `validation/qw/test_qw_convergence.py` and `scripts/lecture_11_convergence.py` provide starting points.

---

## Dependencies / Assumptions

- The existing `star_helpers.py` parsers (`parse_eigenvalues`, `parse_gfactor`, `parse_absorption`, `extract_effective_mass`) are reusable for observable extraction at each grid/FD-order level. Exciton parsers from `scripts/lecture_14_excitons_scattering.py` must be migrated or imported into the shared helper.
- The SC benchmark config used by `verify_sc_benchmarks.py` is suitable for multi-resolution convergence testing (scales to different FDstep values).
- SC Fermi level is output only to per-iteration stdout (`sc_loop.f90:270`), not to a file. The convergence test must either parse stdout (fragile) or a Fortran file output must be added (e.g., `output/sc_summary.dat`). This is a planning-level decision.
- SC charge density integral must be computed numerically from `sc_charge.dat` using trapz integration over the actual z-coordinates (which vary with grid resolution).
- An exciton-enabled config must be created — no existing test config has exciton enabled (only `wire_gaas_optical_window.cfg` has `exciton: F`). The `_EXCITON_CONFIG_TEMPLATE` in `scripts/lecture_14_excitons_scattering.py` provides a working starting point.
- Physics-appropriate convergence rate tolerance must be calibrated empirically per system — material interface discontinuities and non-parabolic bands degrade the clean h^p scaling seen in operator-level tests.
- Wire convergence tests may be significantly more expensive per grid level than QW tests due to CSR assembly and sparse eigensolve. Estimated total wall-clock for the full suite: 30+ minutes (existing wire tests timeout at 600s each; 4-5 grid levels per system).
- Effective mass extraction (R7) depends on k-range. The k-sweep parameters (waveVectorMax, waveVectorStep) must be held constant across grid resolutions.

---

## Outstanding Questions

### Deferred to Planning

- [Affects R1, R2][Needs research] What grid spacing range provides clean convergence for each system? Lecture_11 uses FDstep = [51, 101, 201, 401] for GaAs/AlGaAs QW; other systems may need different ranges.
- [Affects R2][Needs research] What physics-appropriate tolerance for convergence rate assertion? The operator-level tests use 0.05-0.5 (order-adaptive); physics observables will be looser due to material interfaces.
- [Affects R17][Technical] How to structure the shared convergence helper — as an extension of `star_helpers.py`, a new module alongside it, or a separate convergence framework module.
- [Affects R4, R5][Technical] Whether FD order 10 is reliable enough for physics-observable convergence given the past coefficient bug. The operator-level tests pass, but physics-level testing at order 10 may reveal different behavior.
