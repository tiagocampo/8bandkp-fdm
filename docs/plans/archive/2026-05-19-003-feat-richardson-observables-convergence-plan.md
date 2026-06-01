---
title: Richardson Extrapolation on Physics Observables
type: feat
status: active
date: 2026-05-19
origin: docs/brainstorms/richardson-observables-requirements.md
---

# Richardson Extrapolation on Physics Observables

## Summary

A convergence validation suite that runs QW and wire standard-star systems at multiple grid spacings and FD orders, extracts physics observables (subband energy, effective mass, g-factor, absorption edge, SC outputs, exciton binding energy), performs Richardson extrapolation with GCI uncertainty quantification, and asserts convergence rates. Eight implementation units: shared convergence helper module, exciton smoke test, SC Fortran output enhancement, then five convergence test groups. All tests run under an on-demand `convergence` ctest label with JSON diagnostics.

---

## Problem Frame

The codebase supports FD orders 2-10 and arbitrary grid refinement, but no test validates that physics observables converge at the theoretically expected rate. Standard-star benchmarks report numbers at single resolutions with no discretization uncertainty. The SC loop test checks that convergence happens, not that converged physics is correct. The exciton module has zero test coverage. Richardson extrapolation on physics observables closes these gaps with convergence certificates (rate + GCI) for every observable. (see origin: `docs/brainstorms/richardson-observables-requirements.md`)

---

## Requirements

- R1. Grid convergence: 4-5 grid spacings per system, Richardson extrapolation, GCI with Roache safety factors.
- R2. Convergence rate assertion: observed rate matches theoretical FD order rate within empirically calibrated tolerance per system. Material interfaces degrade clean h^p scaling; tolerance is empirically established.
- R3. GCI with Roache safety factors: Fs=1.25 for three-grid (from finest 3 of 4-5 levels).
- R4. FD order convergence: QW systems at orders {2, 4, 6, 8, 10}, wire at {2, 4, 6}. Verify monotonic convergence toward Richardson limit.
- R5. Higher FD order => faster convergence when grid is resolved. Orders 8/10 are informational at QW material interfaces.
- R6. Extract and converge subband energy (CB1, CB2, subband spacing).
- R7. Extract and converge effective mass via parabolic fit near k=0.
- R8. Extract and converge g-factor (gz component).
- R9. Extract and converge absorption edge position. Energy grid held constant across grid spacings.
- R10. Extract and converge SC physics: Fermi level, subband shift, charge density integral. SC rate is empirically measured (not compared to theoretical FD order). SC_tolerance set 10x finer than expected FD error.
- R11. Extract and converge exciton binding energy. Validated against published experimental values.
- R11a. Prerequisite exciton smoke test: single-resolution, verify binding energy 1-20 meV for GaAs/AlGaAs QW.
- R12. Systems: S4 (GaAs/AlGaAs QW), S5 (InAs/GaSb QW), S6 (InAs/GaAs strained QW), S7 (InAs wire). Bulk excluded.
- R13. SC convergence benchmark at multiple grid resolutions.
- R14. Exciton convergence benchmark at multiple grid resolutions and FD orders.
- R15. All convergence tests under `convergence` ctest label only, on-demand, not gating CI.
- R16. JSON results per system per observable: Richardson-extrapolated values, GCI, observed rates, pass/fail.
- R17. Shared convergence helper module with reusable Richardson fitting, GCI, rate extraction. Extracted iteratively from first test scripts.
- R18. Failure diagnostics identify which grid level or FD order caused failure and the specific observable value.

**Origin actors:** Test developer (writes and maintains convergence tests), CI system (runs tests on demand)
**Origin flows:** F1 (grid convergence sweep), F2 (FD order convergence sweep), F3 (SC convergence), F4 (exciton convergence)
**Origin acceptance examples:** AE1 (covers R1,R2,R6,R12 — S4 grid convergence), AE2 (covers R4,R5,R6,R12 — S4 order convergence), AE3 (covers R10,R13 — SC benchmark), AE4 (covers R11,R14 — exciton benchmark), AE5 (covers R15,R18 — ctest label + diagnostics)

---

## Scope Boundaries

- CI gating — convergence tests run on demand only, never blocking merges
- BdG / topological invariant convergence — separate physics modules with different convergence semantics
- Cross-code convergence — already covered by `validation/qw/test_qw_convergence.py`
- Automated golden file updates from Richardson-extrapolated values
- FD operator-level convergence — already covered by `tests/unit/test_fd_convergence.pf`
- Bulk system convergence (S1-S3) — 8x8 dense Hamiltonian has no spatial discretization
- New physics observables beyond the six listed (subband energy, effective mass, g-factor, absorption edge, SC outputs, exciton binding energy)
- Refactoring existing `lecture_11` or `validation/qw` Richardson code (they remain as-is for their respective purposes)
- Fixing the 4 SC configs with overlapping layer patterns (separate concern flagged in `docs/solutions/logic-errors/standard-star-benchmark-suite-logic-errors-2026-05-09.md`)

### Deferred to Follow-Up Work

- Convergence rate calibration for S5 (InAs/GaSb broken-gap) and S6 (InAs/GaAs strained): these systems may exhibit non-standard convergence behavior due to interband coupling and strain discontinuities. Initial implementation uses S4 as the calibration baseline; S5/S6 tolerances are tuned empirically during implementation.
- Wire FD order convergence beyond {2, 4, 6}: once grid convergence proves stable at orders 2-6, higher orders can be added.
- Absorption edge convergence for wire systems: wire optical output is less mature than QW; defer until wire optics are validated.

---

## Context & Research

### Relevant Code and Patterns

- `tests/integration/star_helpers.py` — shared parsers: `parse_eigenvalues`, `parse_gfactor`, `parse_absorption`, `extract_effective_mass`, `compare_value`, `run_exe`. Physical constants and tolerance tiers. All parsers are reusable for convergence tests.
- `validation/qw/test_qw_convergence.py` — existing Richardson extrapolation: `richardson_extrapolate(h_vals, E_vals, order=2)` using two-finest-grid formula. JSON output to `validation/qw/results/qw_convergence.json`. No GCI, no variable order, no pairwise rate extraction.
- `scripts/lecture_11_convergence.py` — alternative Richardson: `richardson_extrapolation(dz_arr, e_arr)` with different formula, `estimate_convergence_rate()` via log-log fit. FDstep sweep [51, 101, 201, 401] for GaAs/AlGaAs QW. No GCI.
- `scripts/lecture_14_excitons_scattering.py` — exciton parsers: `parse_exciton_stdout()`, `parse_exciton_file()`, `_EXCITON_CONFIG_TEMPLATE`. Tests well widths [30-200] A, validates Eb in [1, 20] meV.
- `tests/integration/verify_sc_benchmarks.py` — SC benchmark verifier. Parses SC convergence from stdout. No Fermi level file output (only per-iteration `mu:` lines in stdout). No charge density integration.
- `tests/CMakeLists.txt` — test registration: Python-based tests use `add_test(NAME ... COMMAND python3 ... ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR})` with `LABELS` and `TIMEOUT` properties.
- `validation/shared/comparison.py` — `write_json_report()` pattern for structured JSON output.
- `docs/solutions/best-practices/2026-05-11-richardson-convergence-fixture-fd-operators.md` — max-rate convergence strategy: take maximum observed rate across consecutive pairs (captures asymptotic sweet spot where truncation error dominates but round-off has not contaminated). Order-adaptive tolerance: 0.05 for orders 2-4, 0.2 for 6-8, 0.5 for 10.
- Wire configs at varying resolution: `wire_gaas_21x21.cfg`, `wire_gaas_26x26.cfg`, `wire_gaas_31x31.cfg`, `wire_gaas_41x41.cfg` — but these vary wire SIZE, not grid spacing. New configs needed for FD convergence.
- QW config pattern: `qw_gaas_algaas.cfg` (subband/k-sweep), `qw_gaas_algaas_optics.cfg` (absorption), `gfactor_qw_cb.cfg` (g-factor) — 2-layer pattern (barrier full domain, well center).

### Institutional Learnings

- **Max-rate convergence strategy** (`docs/solutions/best-practices/2026-05-11-richardson-convergence-fixture-fd-operators.md`): Use maximum observed rate across consecutive grid pairs. Round-off causes rates to drop, never increase, so the maximum is the best estimator. Adopt for physics observables.
- **FD order 10 coefficient bug** (`docs/solutions/logic-errors/2026-05-11-wrong-d1-order10-fd-coefficients.md`): Silently wrong coefficients at order 10 — no runtime error. Operator-level tests now catch this, but physics-level convergence provides an additional safety net.
- **8-band verification ladder** (`docs/solutions/best-practices/8band-verification-ladder-2026-05-09.md`): CB effective mass from 8-band is NOT the Vurgaftman value (29% deviation for GaAs). Compare against Kane model `m*_kane = Eg/(EP+Eg)`. Use numerical differentiation at first nonzero k-point, not polynomial fits.
- **SC convergence detection** (`docs/solutions/logic-errors/cross-validation-review-logic-errors-2026-05-19.md`): Use exact substring `"SC loop converged"` for convergence detection, never loose substring matching that could match "NOT converged".
- **Wire test performance** (`docs/solutions/test-failures/lecture-test-pair-review-findings-2026-05-10.md`): Grid size is a test parameter, not a physics requirement. Use smallest grid that validates the property. Wire tests timeout at 600s; convergence sweeps at 4-5 levels need generous timeouts.

### External References

- Roache, P.J. (1998). *Verification and Validation in Computational Science and Engineering*. Hermosa Publishers. — GCI methodology with safety factors.
- Lecture 11 convergence chapter (`docs/lecture/11-convergence.md`): Richardson extrapolation formula (section 11.3.3), GCI formula with Roache safety factors (section 11.3.4), observed convergence rates for QW systems.

---

## Key Technical Decisions

- **New `convergence_helpers.py` module alongside `star_helpers.py`**: Avoids bloating the existing import chain. The helper starts as a minimal skeleton with core Richardson/GCI functions (per R17: extract iteratively from first test scripts), then grows as U4-U8 reveal duplication patterns. Rationale: `star_helpers.py` is imported by all standard-star tests and lecture scripts; adding convergence-specific infrastructure there would couple unrelated concerns.

- **Add Fortran `sc_summary.dat` output**: The SC Fermi level is currently only available from per-iteration stdout (`sc_loop.f90:270` — `'  mu:', fermi_level`). Stdout parsing is fragile (format-dependent, breaks if SC loop output changes). A structured file output is cleaner and more maintainable. The Fortran change requires: (1) add optional `fermi_level_out` argument to both SC subroutine signatures, (2) update call sites in main.f90, (3) write `sc_summary.dat` from main.f90. Note: CLAUDE.md flags sc_loop.f90 changes as requiring approval.

- **Wire convergence uses new configs at fixed geometry**: Existing wire configs (`wire_gaas_21x21.cfg` etc.) vary wire SIZE while keeping spacing fixed — they test physics convergence (wider wire = less confinement), not FD convergence. New configs keep `wire_width`/`wire_height` constant and vary `wire_nx`/`wire_ny`/`wire_dx`/`wire_dy`. Rationale: Richardson extrapolation on h requires the same physical domain with different grid resolutions.

- **Max-rate convergence strategy for physics observables**: Adopted from the FD operator convergence fixture. Take the maximum observed convergence rate across all consecutive grid pairs. This captures the asymptotic "sweet spot" where truncation error dominates but round-off has not contaminated. Rationale: physics observables at material interfaces show non-monotonic convergence rates; the max-rate estimator is more robust than a single finest-pair rate.

- **Each convergence test script is self-contained**: Each script handles its own config generation (parameterized by grid spacing or FD order), Fortran execution, observable extraction, Richardson fitting, assertion, JSON output, and ctest registration. The helper module provides reusable primitives; the scripts compose them. Rationale: avoids a monolithic test framework; each system's convergence test can be run and debugged independently.

- **Exciton parsers remain in U2 scope**: `parse_exciton_stdout()` and `parse_exciton_file()` are first consumed by U2 (smoke test). Per R17, they stay in U2's scope until U8 reveals whether they should be extracted to the shared helper. U2 imports directly from `lecture_14` or copies the parsers inline.

---

## Open Questions

### Resolved During Planning

- **SC Fermi level extraction method**: Resolved — add Fortran `sc_summary.dat` output (see Key Technical Decisions). Stdout parsing rejected as fragile.
- **Shared helper module structure**: Resolved — new `convergence_helpers.py` alongside `star_helpers.py`, not extending it.
- **Wire convergence config strategy**: Resolved — new configs at fixed geometry with varying grid spacing.

### Deferred to Implementation

- **Exact grid spacing ranges per system**: Lecture_11 uses FDstep = [51, 101, 201, 401] for GaAs/AlGaAs QW. Other systems need calibration. The first implementation run will establish which ranges produce clean convergence. Start with [51, 101, 201, 401] for all QW systems and adjust empirically.
- **Physics-appropriate convergence rate tolerance per system**: Material interface discontinuities and non-parabolic bands degrade clean h^p scaling. The operator-level tests use 0.05-0.5 (order-adaptive); physics observables will be looser. Establish empirically: run the convergence sweep, compute observed rates, then set tolerance as max-rate * margin. Lecture 11 shows rates of 1.2-1.9 for FDorder=2 (theoretical 2.0) — so initial tolerance might be ~0.3-0.5 for order 2.
- **FD order 10 reliability at physics level**: The past coefficient bug was fixed and operator-level tests pass, but physics-level convergence at order 10 may reveal stencil boundary effects at material interfaces. R5 already marks orders 8/10 as informational — assertion tolerance is relaxed accordingly. Defer further investigation to implementation results.
- **Wire convergence test timeout**: 4-5 grid levels of wire CSR eigensolve at 600s per level gives ~30+ minutes. May need TIMEOUT 3600 or individual test timeouts in ctest. Calibrate from first wire convergence run.
- **Exciton config parameters for convergence**: The `_EXCITON_CONFIG_TEMPLATE` in lecture_14 uses FDstep=51 and well widths [30-200] A. The convergence sweep needs FDstep parameterization and a fixed well width. Exact parameters deferred.

### From Document Review (2026-05-20)

- **Calibration Phase A/B split**: Consider splitting the convergence suite into Phase A (informational: produce Richardson values and observed rates without asserting) and Phase B (assertive: rate assertions with calibrated tolerances). Add a minimum-observed-rate floor (e.g., rate > 1.0 for order 2) as a non-calibrated hard gate that catches gross failures.
- **Richardson order estimation at material interfaces**: Variable-order Richardson requires known p. At material interfaces, the effective order is degraded. Add an order-estimation mode using the Roache three-grid formula (estimate p from finest three levels). Use estimated p when observed rate deviates significantly from theoretical.
- **Effective mass extraction method**: R7 specifies parabolic fit but institutional learnings recommend numerical differentiation at first nonzero k-point. The adaptive fit in `extract_effective_mass` uses different k-ranges at different resolutions, introducing systematic bias. Consider fixing the number of k-points or using numerical differentiation.
- **Absorption edge extraction algorithm**: Specify the method: find the energy where absorption crosses a fixed fraction of the peak value using linear interpolation. No edge extraction function currently exists in star_helpers.py.
- **SC config verification (U7)**: Verify that `sc_gaas_alas_qw.cfg` well boundaries at z=+/-50 A align with grid points at ALL FDstep values in the sweep. If not, adjust domain width for exact alignment. Add an assertion in the convergence test checking boundary alignment.
- **SC_tolerance scaling**: SC_tolerance is on potential |dPhi|, not eigenvalue error. Verify empirically that SC eigenvalues don't change when tolerance tightens from 1e-8 to 1e-10. Consider grid-spacing-adaptive tolerance.
- **U5 merge into U4**: U5 depends heavily on U4 output and covers the same systems/observables. Consider merging as a second mode within U4's test script, using a `--mode` flag for separate ctest registration.
- **Exciton smoke test decoupling**: U2 could ship independently as a standalone PR (import parsers from lecture_14 directly) for immediate exciton module regression protection, without waiting for convergence_helpers.py.
- **Wire grid sizing**: The proposed [11x11, 16x16, 21x21, 26x26] range gives grid ratio ~1.36, below the typical Richardson recommendation of >= 1.5. Consider starting at 16x16 or adding a 31x31 level.
- **Lecture 11 consistency**: Lecture 11 section 11.3.4 already references the convergence fixture as if it exists. Verify consistency after implementation.

---

## Implementation Units

- U1. **Shared convergence helper module**

**Goal:** Create a minimal `convergence_helpers.py` skeleton with core Richardson extrapolation and GCI. Additional functions extracted iteratively as duplication patterns emerge from U4-U8 (per R17).

**Requirements:** R17

**Dependencies:** None

**Files:**
- Create: `tests/integration/convergence_helpers.py`
- Reference: `tests/integration/star_helpers.py` (parser imports)
- Reference: `scripts/lecture_11_convergence.py` (Richardson formula starting point)
- Reference: `validation/qw/test_qw_convergence.py` (Richardson + JSON pattern)

**Approach:**
- Start with only `richardson_extrapolate(h_vals, observable_vals, order)` with variable-order support and `compute_gci(h_vals, observable_vals, order, safety_factor=1.25)`. These are the proven duplication points from the two existing implementations.
- Add `extract_convergence_rates()`, `max_convergence_rate()`, observable extraction wrappers, JSON output, and failure diagnostics incrementally as U4-U8 reveal the actual duplication pattern (per R17: "extracted from the first few convergence test scripts once the actual code duplication pattern is clear, not designed upfront").
- Exciton parsers remain in U2's scope (U2 is the first and only consumer until U8). Do not migrate to convergence_helpers upfront.

**Execution note:** Write unit tests for Richardson/GCI functions using known analytical solutions. Do not pre-design the full API — let it emerge from U4-U8.

**Patterns to follow:**
- JSON output from `validation/shared/comparison.py`: `os.makedirs` + `json.dump(results, f, indent=2)`
- Max-rate strategy from `docs/solutions/best-practices/2026-05-11-richardson-convergence-fixture-fd-operators.md`
- Parser structure from `tests/integration/star_helpers.py`

**Test scenarios:**
- Happy path: Richardson extrapolation with exact h^2 data produces E_exact within machine epsilon
- Happy path: GCI computation with known three-grid data matches hand-calculated Roache GCI
- Happy path: Convergence rate extraction from h^2 data returns rate ~2.0
- Edge case: Richardson with non-monotonic data (observables not decreasing with h) returns gracefully with diagnostic
- Edge case: Two-grid Richardson uses Fs=3 safety factor; three-grid uses Fs=1.25
- Edge case: Variable-order Richardson with p=4 produces correct extrapolation
- Error path: Empty or single-point input raises informative error
- Integration: JSON output matches existing schema from validation pipeline

**Verification:**
- All helper function unit tests pass
- Richardson extrapolation on analytical test data (e.g., E(h) = E_exact + c*h^2) recovers E_exact to machine precision

---

- U2. **Exciton smoke test**

**Goal:** Single-resolution test verifying the exciton module runs correctly and produces a physically plausible binding energy (1-20 meV for GaAs/AlGaAs QW). Must pass before R11 convergence work begins.

**Requirements:** R11a

**Dependencies:** U1 (imports exciton parsers from convergence_helpers)

**Files:**
- Create: `tests/integration/test_exciton_smoke.py`
- Create: `tests/regression/configs/exciton_gaas_algaas.cfg` (exciton-enabled config based on `_EXCITON_CONFIG_TEMPLATE`)
- Modify: `tests/CMakeLists.txt` (add test registration)
- Reference: `scripts/lecture_14_excitons_scattering.py` (config template reference)

**Approach:**
- Create exciton-enabled config: GaAs/AlGaAs QW, well width ~80 A, `exciton: T`, `method: variational`, FDstep=51. Based on the `_EXCITON_CONFIG_TEMPLATE` from lecture_14.
- Run `bandStructure` (or `opticalProperties` — whichever produces exciton output) with the config.
- Parse output using migrated exciton parsers from convergence_helpers.
- Assert: binding energy in [1, 20] meV, no runtime errors, exciton output file exists.
- Register under `convergence` label (per R15 — even the smoke test carries this label). Note: consider using a separate label (e.g., `exciton`) since this is a smoke test, not a convergence test. If shipped independently before the convergence suite, use the existing label structure.

**Patterns to follow:**
- Test script pattern from standard-star verifiers: `sys.argv[1:]` for build_dir/source_dir, `run_exe()` from star_helpers, parsers for output, `sys.exit(0/1)`.
- Config template from `scripts/lecture_14_excitons_scattering.py:_EXCITON_CONFIG_TEMPLATE`.
- ctest registration pattern from `tests/CMakeLists.txt`.

**Test scenarios:**
- Happy path: Exciton binding energy falls in [1, 20] meV for GaAs/AlGaAs QW
- Happy path: Exciton output file (`exciton.dat`) exists and is parseable
- Happy path: Variational parameters (lambda_opt, mu, eps_r) are physically reasonable
- Error path: Fortran executable fails to run → test reports clear error, exits 1
- Edge case: Binding energy outside [1, 20] meV → test reports value and fails

**Verification:**
- `ctest -L convergence -R exciton_smoke` passes
- Exciton binding energy matches lecture_14 result for same config within tolerance

---

- U3. **SC Fortran output enhancement**

**Goal:** Add structured `sc_summary.dat` file output to the Fortran SC loop, containing Fermi level, converged flag, and iteration count. This replaces fragile stdout parsing for Fermi level extraction.

**Requirements:** R10 (indirect — enables reliable SC observable extraction)

**Dependencies:** None (Fortran change independent of Python test infrastructure)

**Files:**
- Modify: `src/physics/sc_loop.f90` (add optional `fermi_level_out` argument to `self_consistent_loop` and `self_consistent_loop_wire` signatures)
- Modify: `src/apps/main.f90` (update call sites at lines ~613 and ~235 to receive fermi_level; add sc_summary.dat writing after SC loop completes)
- Test: `tests/unit/test_sc_loop.pf` (update if it calls self_consistent_loop directly)
- Test: `tests/integration/test_sc_convergence.py` (U7 will use this output)

**Approach:**
- The Fortran change is larger than initially estimated: `fermi_level` is a local variable in both SC subroutine signatures with no output parameter. Requires: (1) add optional `fermi_level_out` argument to both `self_consistent_loop` (QW) and `self_consistent_loop_wire`, (2) assign at end of each subroutine, (3) update call sites in main.f90, (4) write `sc_summary.dat` from main.f90 using the returned value.
- Note: CLAUDE.md flags this as requiring approval ("changes to sc_loop.f90 convergence logic").
- Format: single-line header + single data line. Example: `# converged  iterations  |dPhi|  fermi_level(eV)` followed by `T  23  9.87e-07  0.567`.
- Keep the existing stdout output unchanged — the new file is supplementary.

**Execution note:** Build and run existing SC tests to verify no regression before and after the change.

**Patterns to follow:**
- Fortran output file pattern from `src/apps/main.f90` (existing `output/sc_potential_profile.dat`, `output/sc_charge.dat` writing).
- Single-header + single-data-line format for simple parsing.

**Test scenarios:**
- Happy path: `sc_summary.dat` exists after SC-enabled run, contains correct Fermi level
- Happy path: Fermi level in file matches `mu:` value in stdout
- Edge case: SC does not converge — file still written with `converged=F` and final Fermi level
- Edge case: Non-SC run — no `sc_summary.dat` created (no file pollution)

**Verification:**
- Existing SC tests (`ctest -R sc`) still pass after the change
- Manual run of `bandStructure` with SC config produces `sc_summary.dat` with plausible values

---

- U4. **QW grid convergence tests**

**Goal:** Grid convergence tests for S4 (GaAs/AlGaAs QW), S5 (InAs/GaSb QW), S6 (InAs/GaAs strained QW) at fixed FD order, extracting subband energy, effective mass, g-factor, and absorption edge. Covers AE1.

**Requirements:** R1, R2, R3, R6, R7, R8, R9, R12

**Dependencies:** U1 (convergence_helpers)

**Files:**
- Create: `tests/integration/test_qw_grid_convergence.py`
- Create: config templates for S4-S6 grid sweeps (parameterized by FDstep)
- Modify: `tests/CMakeLists.txt` (add test registration)
- Test: `tests/integration/test_qw_grid_convergence.py` (self-testing script)

**Approach:**
- For each QW system (S4, S5, S6), generate configs at 4-5 grid spacings (starting from [51, 101, 201, 401], adjusted per system).
- Run `bandStructure` for subband energy and effective mass extraction. Run `gfactorCalculation` for g-factor. Run `opticalProperties` for absorption edge.
- Energy grid for absorption held constant across resolutions (R9).
- k-sweep parameters (waveVectorMax, waveVectorStep) held constant across grid resolutions for effective mass extraction (R7), analogous to the constant energy grid for absorption.
- Extract observables via convergence_helpers wrappers around star_helpers parsers.
- Compute Richardson extrapolation (variable order matching FDorder), GCI, and convergence rates.
- Assert convergence rate matches theoretical FD order rate within empirically calibrated tolerance.
- Write JSON results per system per observable.

**Execution note:** Start with S4 (GaAs/AlGaAs) to calibrate grid spacings and tolerances. S5 and S6 follow the calibrated pattern.

**Patterns to follow:**
- Grid sweep from `scripts/lecture_11_convergence.py`: FDstep = [51, 101, 201, 401].
- Observable extraction from `star_helpers.py`: `parse_eigenvalues`, `extract_effective_mass`, `parse_gfactor`, `parse_absorption`.
- Config generation: parameterize existing QW configs by FDstep (write to temp file per resolution).

**Test scenarios:**
- Covers AE1. S4 GaAs/AlGaAs at FDorder=2, 4 grid spacings → Richardson extrapolation, GCI, rate assertion for CB1 energy
- Happy path: S4 CB1 convergence rate ~2.0 for FDorder=2 (within tolerance ~0.3)
- Happy path: S4 effective mass convergence with Richardson extrapolation
- Happy path: S4 g-factor convergence with Richardson extrapolation
- Happy path: S4 absorption edge convergence with Richardson extrapolation
- Happy path: GCI uncertainty band is positive and reasonable (fraction < 5% for well-resolved grids)
- Edge case: Non-monotonic convergence at coarsest grid → max-rate strategy still produces valid rate
- Integration: All S4/S5/S6 observables produce JSON output with Richardson, GCI, rate, pass/fail
- Integration: JSON schema matches convergence_helpers format

**Verification:**
- `ctest -L convergence -R qw_grid` passes for all three QW systems
- JSON results files contain Richardson-extrapolated values and GCI for all observables

---

- U5. **QW FD order convergence tests**

**Goal:** FD order convergence tests for S4-S6 at fixed grid spacing, running at orders {2, 4, 6, 8, 10}. Verify monotonic convergence toward Richardson limit and rate acceleration with higher order. Covers AE2.

**Requirements:** R4, R5, R6, R12

**Dependencies:** U1 (convergence_helpers), U4 (calibrated grid spacings and Richardson limit from grid convergence)

**Files:**
- Create: `tests/integration/test_qw_order_convergence.py`
- Modify: `tests/CMakeLists.txt` (add test registration)

**Approach:**
- For each QW system, fix FDstep at the finest grid from U4. Vary FDorder in {2, 4, 6, 8, 10}.
- Extract CB1 energy (primary observable) at each FD order.
- The convergence target is the finest-grid result at the highest reliable FD order — not U4's order-specific Richardson limit. Each FD order has its own discretization error model; using U4's p=2 Richardson limit as a target for higher orders is not theoretically sound.
- Verify monotonic convergence: observable values approach the common answer as FD order increases.
- Verify rate acceleration: convergence rate (log-log slope vs order) increases with FD order.
- Orders 8 and 10 are informational — relaxed assertion tolerance per R5.

**Patterns to follow:**
- FD order sweep from `scripts/lecture_11_convergence.py`: FDorder sweep at fixed FDstep.
- Richardson limit from U4 grid convergence as the convergence target.

**Test scenarios:**
- Covers AE2. S4 at FDstep=401 with orders {2, 4, 6, 8, 10} → monotonic convergence toward Richardson limit
- Happy path: CB1 energy monotonically approaches Richardson limit across orders
- Happy path: Convergence rate accelerates (order 4 converges faster than order 2, etc.)
- Edge case: Order 10 may not improve over order 8 at material interfaces → informational, not a failure
- Edge case: Order 2 vs 4 gap is larger than order 4 vs 6 gap (diminishing returns expected)
- Integration: JSON results include per-order values and monotonicity flag

**Verification:**
- `ctest -L convergence -R qw_order` passes
- Order convergence JSON shows monotonic approach to Richardson limit for S4

---

- U6. **Wire convergence tests**

**Goal:** Grid convergence tests for S7 (InAs wire) at fixed geometry, extracting subband energy and g-factor. FD order sweep restricted to {2, 4, 6}.

**Requirements:** R1, R2, R3, R6, R8, R12

**Dependencies:** U1 (convergence_helpers)

**Files:**
- Create: `tests/integration/test_wire_convergence.py`
- Create: wire convergence configs at fixed geometry, varying grid spacing (e.g., InAs 55x55 nm rectangle at 11x11, 16x16, 21x21, 26x26 grids)
- Modify: `tests/CMakeLists.txt` (add test registration with generous TIMEOUT)

**Approach:**
- New wire configs at fixed physical dimensions (same wire_width, wire_height) with varying grid resolution (wire_nx, wire_ny, wire_dx, wire_dy). Start with InAs rectangle ~55x55 nm at grids [11x11, 16x16, 21x21, 26x26].
- Run `bandStructure` for subband energy, `gfactorCalculation` for g-factor.
- FEAST eigensolver: use `feast_m0=-1` (dense LAPACK fallback) for smaller grids to avoid FEAST non-determinism. For larger grids, use FEAST with consistent `feast_emin`/`feast_emax`.
- Only use k=0 eigenvalues (FEAST non-deterministic at k>0).
- Compute Richardson, GCI, rates. Set generous ctest TIMEOUT (1800s per grid level).
- FD order sweep: orders {2, 4, 6} only, per R4. Skip higher orders until grid convergence proves stable.

**Patterns to follow:**
- Wire config pattern from existing `wire_inas_rectangle.cfg`: `confinement: 2`, `wire_nx`, `wire_ny`, etc.
- k=0-only extraction from `tests/regression/compare_output.py` wire test (first data line of eigenvalues.dat).

**Test scenarios:**
- Happy path: S7 subband energy converges across grid levels with positive GCI
- Happy path: S7 g-factor converges across grid levels
- Edge case: FEAST returns different eigenvalue count at different resolutions → handle gracefully (compare only available bands)
- Edge case: Coarsest wire grid (11x11) may have poor resolution → convergence rate may be lower than theoretical
- Integration: JSON results for wire grid and order convergence

**Verification:**
- `ctest -L convergence -R wire` passes (with generous timeout)
- Wire convergence JSON contains Richardson values for subband energy and g-factor

---

- U7. **SC convergence benchmark**

**Goal:** Multi-resolution SC loop runs extracting Fermi level, subband energy shift, and charge density integral. Monotonic convergence with GCI — no theoretical rate assertion (SC convergence is empirically determined). Covers AE3.

**Requirements:** R10, R13

**Dependencies:** U1 (convergence_helpers), U3 (sc_summary.dat for Fermi level)

**Files:**
- Create: `tests/integration/test_sc_convergence.py`
- Modify: `tests/CMakeLists.txt` (add test registration)

**Approach:**
- Use the existing SC benchmark config (`sc_gaas_alas_qw.cfg`) as template, parameterized by FDstep.
- Run at 4 grid resolutions (e.g., FDstep = [51, 101, 201, 401]).
- SC_tolerance set to 1e-8 (at least 10x finer than expected FD error at finest grid) to ensure Richardson measures FD convergence, not SC iteration residual.
- Extract: Fermi level from `sc_summary.dat` (U3), subband energy shift relative to flat-band from eigenvalues, charge density integral from `sc_charge.dat` using numpy trapezoid integration over z-coordinates.
- Assert monotonic convergence, compute Richardson extrapolation and GCI.
- SC rate is empirically measured (not compared to theoretical FD order) per R10.
- Write JSON results.

**Patterns to follow:**
- SC config from `tests/regression/configs/sc_gaas_alas_qw.cfg`.
- Charge density file format: `sc_charge.dat` with `# z(A) n_e(cm^-3) n_h(cm^-3)` per line.
- Convergence detection: exact substring `"SC loop converged"` (not loose matching).

**Test scenarios:**
- Covers AE3. SC benchmark at 4 resolutions → Fermi level, subband shift, charge integral with Richardson + GCI
- Happy path: Fermi level converges monotonically across resolutions
- Happy path: Subband energy shift converges monotonically
- Happy path: Charge density integral converges monotonically
- Happy path: Richardson extrapolation produces finite GCI band for all three observables
- Edge case: SC loop does not converge at coarsest grid → skip that resolution, report diagnostic
- Edge case: SC_tolerance too loose → Richardson measures SC residual, not FD error → diagnostic warning
- Integration: JSON results with SC-specific fields (sc_tolerance, iterations_per_resolution)

**Verification:**
- `ctest -L convergence -R sc_conv` passes
- SC convergence JSON contains Richardson-extrapolated Fermi level, subband shift, and charge integral with GCI bands

---

- U8. **Exciton convergence benchmark**

**Goal:** Multi-resolution exciton runs extracting binding energy, performing Richardson extrapolation, and validating against published experimental values (Miller et al., Phys. Rev. B 1985). Covers AE4.

**Requirements:** R11, R14

**Dependencies:** U1 (convergence_helpers), U2 (smoke test must pass first)

**Files:**
- Create: `tests/integration/test_exciton_convergence.py`
- Modify: `tests/CMakeLists.txt` (add test registration)

**Approach:**
- Exciton-enabled config based on U2's smoke test config, parameterized by FDstep.
- Run at 4 grid resolutions (e.g., FDstep = [51, 101, 201, 401]).
- Extract binding energy at each resolution using exciton parsers from convergence_helpers.
- Richardson extrapolation + GCI.
- Validate extrapolated value against published experimental values as primary reference.
- Fallback to 2D hydrogen model (Rydberg/4) as sanity check if no published value matches the config.
- Always verify internal consistency of Richardson-extrapolated reference.
- Also run FD order sweep at fixed FDstep (orders {2, 4, 6}) for order convergence (R14).

**Execution note:** U2 smoke test must pass before this work begins. If the smoke test reveals exciton module issues, fix those first.

**Patterns to follow:**
- Config template from U2 (exciton_gaas_algaas.cfg).
- Exciton output parsing from `parse_exciton_file()` in convergence_helpers.
- Validation against published values: Miller et al. binding energy for GaAs/AlGaAs QW at similar well width.

**Test scenarios:**
- Covers AE4. GaAs/AlGaAs exciton at 4 resolutions → binding energy Richardson extrapolation matches published values
- Happy path: Binding energy converges monotonically across resolutions
- Happy path: Richardson-extrapolated value within tolerance of Miller et al. experimental value
- Happy path: GCI band is reasonable (fraction < 10%)
- Edge case: Richardson-extrapolated value does not match published → fallback to 2D hydrogen model, report diagnostic
- Edge case: Binding energy non-monotonic at coarsest grid → diagnostic identifies which grid level diverged (R18)
- Integration: JSON results with exciton-specific fields (experimental_reference, hydrogen_model_comparison)

**Verification:**
- `ctest -L convergence -R exciton_conv` passes
- Exciton convergence JSON contains Richardson-extrapolated binding energy matching published values within tolerance

---

## System-Wide Impact

- **Interaction graph:** Convergence tests invoke the same Fortran executables (`bandStructure`, `gfactorCalculation`, `opticalProperties`) as standard-star tests. No new executables or Fortran modules except `sc_summary.dat` output in `main.f90`.
- **Error propagation:** Convergence test failures produce JSON diagnostics identifying the failing grid level, observable value, and expected rate. No error propagation to other test suites.
- **State lifecycle risks:** Config files generated in temp directories per run (same pattern as standard-star tests). No persistent state risk.
- **API surface parity:** The Fortran SC summary output is additive — no existing output format changes.
- **Integration coverage:** Each convergence test script is end-to-end (config generation → Fortran execution → parsing → Richardson → assertion → JSON). No external service dependencies.
- **Unchanged invariants:** Existing standard-star tests (`ctest -L standard-star`), verification tests (`ctest -L verification`), and unit tests (`ctest -L unit`) are unaffected. The `convergence` label is exclusive to new tests per R15.

---

## Risks & Dependencies

| Risk | Mitigation |
|------|------------|
| Wire convergence tests exceed 30 minutes total | Set generous TIMEOUT per test (1800s). Start with smallest grids that show convergence. Wire order sweep deferred to {2,4,6}. |
| Physics convergence rates at material interfaces don't match theoretical FD order | Use empirically calibrated tolerances per system (not hard theoretical gates). Max-rate strategy captures best-case rate. |
| FD order 10 produces non-monotonic physics convergence | Orders 8/10 are informational per R5. Relaxed assertion tolerance. Non-monotonic results reported but not a failure. |
| SC_tolerance too loose → Richardson measures SC residual instead of FD error | Set SC_tolerance to 1e-8 (10x finer than expected FD error). Diagnostic warning if SC iterations are excessive. |
| Exciton module bugs block convergence testing | R11a smoke test runs first. If it fails, exciton convergence is blocked until the module is fixed. |
| Exciton parsers fragile (regex-based stdout parsing) | Migrate parsers to convergence_helpers with error handling. Parse `exciton.dat` file (structured) preferentially over stdout. |
| Convergence rate calibration is circular (need convergence data to set tolerance) | Bootstrap: run convergence sweep first, compute rates, then set tolerance as observed_max_rate * margin. Initial runs are informational. |

---

## Documentation / Operational Notes

- Lecture 11 convergence chapter (`docs/lecture/11-convergence.md`) provides pedagogical context for Richardson extrapolation methodology. No update needed — convergence tests are orthogonal to the lecture content.
- `docs/reference/input-reference.md` may need update if exciton config parameters are not already documented.
- The `convergence` ctest label should be documented in CLAUDE.md's Testing section alongside existing labels.

---

## Sources & References

- **Origin document:** [docs/brainstorms/richardson-observables-requirements.md](docs/brainstorms/richardson-observables-requirements.md)
- Existing Richardson: [validation/qw/test_qw_convergence.py](validation/qw/test_qw_convergence.py), [scripts/lecture_11_convergence.py](scripts/lecture_11_convergence.py)
- Standard-star infrastructure: [tests/integration/star_helpers.py](tests/integration/star_helpers.py)
- Exciton infrastructure: [scripts/lecture_14_excitons_scattering.py](scripts/lecture_14_excitons_scattering.py)
- SC benchmark: [tests/integration/verify_sc_benchmarks.py](tests/integration/verify_sc_benchmarks.py)
- FD convergence methodology: [docs/solutions/best-practices/2026-05-11-richardson-convergence-fixture-fd-operators.md](docs/solutions/best-practices/2026-05-11-richardson-convergence-fixture-fd-operators.md)
- Convergence lecture: [docs/lecture/11-convergence.md](docs/lecture/11-convergence.md)
- Test registration: [tests/CMakeLists.txt](tests/CMakeLists.txt)
- SC loop output: [src/physics/sc_loop.f90](src/physics/sc_loop.f90) (Fermi level per-iteration stdout at line 270)
