---
date: 2026-05-12
topic: fix-code-review-findings
status: completed
origin: docs/brainstorms/2026-05-12-code-review-fixes-requirements.md
---

# Fix Code Review Findings — Topological/Superconducting Branch

## Summary

Fix 29 findings from the code review of `feature/bdg-topological-superconductivity`: BdG Hermiticity verification, topological invariant correctness (Z2, Fu-Kane, phase diagrams), numerical regression bias, dead code removal, missing Fortran attributes, and Python script robustness — with tests for each fix.

---

## Problem Frame

The branch adds topological invariant computation (Chern, Z2, Majorana via BdG) across 206 commits and ~11K lines of Fortran source/test changes. A multi-agent code review found 29 issues: 14 P1, 14 P2, and 1 P3. Several P1 findings affect physics correctness — the BdG Hamiltonian may lose Hermiticity when Peierls phases are applied, the Z2 invariant is tracked incorrectly in phase diagrams, and Majorana detection produces false positives when `delta_0=0`. The branch's purpose is new physics; these correctness issues must be resolved for the physics to be trustworthy.

---

## Key Technical Decisions

- **Dead code first:** Remove dead code (U9) before physics fixes to shrink the diff surface and avoid fixing code that will be deleted.
- **Physics fixes ordered by dependency:** Z2 transition logic (U2) before phase diagram cumulative tracking (U5) and QSHE wire Z2 (U6), since the latter two depend on correct transition detection.
- **R1 is a false positive:** Research confirmed the BdG hole block construction `-H0^T` is mathematically correct because H0 contains both (i,j) and (j,i) entries. No Fortran fix — only a verification test (R25).
- **R8 uses gap-closing + cumulative tracking:** QSHE wire Z2 will use gap-closing detection with cumulative Z2 tracking rather than the Pfaffian method, consistent with R7's approach and simpler to implement.
- **R9 spatial localization:** `compute_z2_gap` will check spatial localization (edge vs bulk weight) of near-zero eigenvalues, not just energy threshold.
- **R17 shared helper:** Extract tail-search and log-linear fit from `fit_exponential_decay` and `compute_majorana_profile` into a private helper to eliminate duplication.
- **File splits deferred:** `topological_analysis.f90` (1287 lines) and `main_topology.f90` (1234 lines) will be split in a follow-up task per scope boundaries.

---

## Implementation Units

### U1: delta_0 guard + pairing sign constant

**Requirements:** R3, R4
**Files:**
- `src/physics/bdg_hamiltonian.f90`
- `src/apps/main_topology.f90`
**Risk:** low — guard adds a check; constant extraction is mechanical

Define the pairing sign antidiagonal array `[+1,-1,+1,-1,-1,+1,+1,-1]` as a module-level constant in `bdg_hamiltonian.f90`. Replace the two duplicated definitions in `build_bdg_hamiltonian_1d` and `build_bdg_hamiltonian_qw`.

Add a guard at the entry of `build_bdg_hamiltonian_1d` (and `_qw`): if `delta_0 <= 0`, print an error message and stop. This prevents silent production of degenerate matrices that trigger false Majorana detection.

In `main_topology.f90` `run_bdg_wire`, the Majorana threshold `0.001*delta_0` becomes unreachable when `delta_0 <= 0` because the BdG build will have already errored.

**Test scenarios:**
- Build BdG with `delta_0 = 0.0`: confirm error stop with message, no Majorana detection runs
- Build BdG with `delta_0 < 0` (e.g., -0.001): same error stop
- Build BdG with valid `delta_0 > 0`: pairing signs match the module constant

---

### U2: is_z2_transition fix

**Requirements:** R5, R28
**Files:**
- `src/physics/topological_analysis.f90`
**Risk:** medium — changes transition detection semantics

Fix `is_z2_transition` (line ~1265) to check that Z2 actually changes between adjacent points. Current logic marks transition when gap is small even if Z2 is identical. The fix: return `.true.` only when `z2_a /= z2_b` (Z2 changed), regardless of gap size.

This is a prerequisite for U5 (phase diagram cumulative tracking) and U6 (QSHE wire Z2), which rely on correct transition detection.

**Test scenarios:**
- Two adjacent points with Z2=0, Z2=0 but small gap: `is_z2_transition` returns `.false.`
- Two adjacent points with Z2=0, Z2=1 with any gap: `is_z2_transition` returns `.true.`
- Two adjacent points with Z2=1, Z2=0: `is_z2_transition` returns `.true.`

---

### U3: Linear regression bias + fitting extraction

**Requirements:** R10, R11, R17, R29
**Files:**
- `src/physics/topological_analysis.f90`
**Risk:** medium — changes regression behavior, may affect existing golden outputs

**R17 first:** Extract the shared tail-search and log-linear fit logic from `fit_exponential_decay` (~line 508) and `compute_majorana_profile` into a private helper. Both routines: (1) find the tail start where density drops below threshold, (2) take log of nonzero values, (3) perform linear regression on (position, log(density)).

**R10:** In the shared helper, track the actual count of points contributing to the regression sums. Skip zero-density points during accumulation and decrement `n_fit` for each skipped point. This eliminates the bias where `n_fit` includes zero-density positions that contribute nothing to the sums.

**R11:** When the localization length xi is comparable to or larger than the wire length (near-transition regime), return a result flagged as "unconverged" with the raw xi value, rather than a sentinel that is indistinguishable from failure. The caller can check convergence status.

**Test scenarios:**
- Synthetic exponential with xi=5.0 and zero-valued tail entries: recovered xi within 10% of true value (R29)
- Near-transition case (xi ≈ wire_length): returns result with unconverged flag, not sentinel
- Existing Majorana localization tests still pass

---

### U4: Fu-Kane Z2 robustness

**Requirements:** R6, R27
**Files:**
- `src/physics/topological_analysis.f90`
**Risk:** high — changes Z2 computation for 2D BHZ model, may affect existing results

Fix `compute_z2_fukane_qw_result` (~lines 992-1098) to pair Kramers partners by parity/momentum structure rather than by consecutive zheev index. Current approach pairs `(istate, istate+1)` which fails when zheev reorders within a near-degenerate subspace at TRIM points.

Approach: compute parity eigenvalue for each band individually, then pair bands with opposite parity at the same energy. If the parity pattern is ambiguous (more than two bands near the same energy), fall back to the minimum-energy-gap pairing criterion but with explicit tie-breaking.

**Test scenarios:**
- Z2 value (not just allocation) verified at phase transition boundary — must be 0 on trivial side, 1 on topological side
- Near-degenerate bands at TRIM point: pairs correctly despite zheev reordering
- Known BHZ model: Z2 matches published results at representative parameter points

---

### U5: Phase diagram cumulative Z2 tracking

**Requirements:** R7
**Files:**
- `src/physics/topological_analysis.f90`
- `src/apps/main_topology.f90`
**Risk:** high — changes phase diagram output semantics

Fix `compute_phase_diagram` (~line 935) to track Z2 cumulatively along the B-sweep. Start at Z2=0 (trivial at B=0) and flip each time the gap closes. The current per-point gap-closing flag approach incorrectly marks individual points as topological without cumulative context.

Depends on U2 for correct transition detection.

**Test scenarios:**
- BHZ model swept B=0 to B=5 with one gap closing at B=2: Z2=0 for B<2, Z2=1 for B>2
- Multiple gap closings: Z2 flips each time (0→1→0→1...)
- No gap closing: Z2=0 throughout

---

### U6: QSHE wire Z2 from spectrum

**Requirements:** R8, R9
**Files:**
- `src/physics/topological_analysis.f90`
- `src/apps/main_topology.f90`
**Risk:** medium — replaces heuristic with physics-based detection

**R8:** Replace the hardcoded wire-width threshold (width >= 70Å) in `run_qshe_wire` with gap-closing detection from the actual eigenspectrum. Use cumulative Z2 tracking consistent with U5.

**R9:** Fix `compute_z2_gap` to distinguish topological edge states from numerical noise. Add a spatial localization check: compute the weight of each near-zero eigenstate at the wire edges vs. the bulk. Only count eigenvalues that are both near-zero energy AND spatially localized at edges as evidence of nontrivial topology. Remove the unused `N` parameter.

Depends on U2 and U5 for cumulative tracking infrastructure.

**Test scenarios:**
- QSHE wire below critical width: Z2=0 (trivial), no edge states detected
- QSHE wire above critical width: Z2=1 (topological), edge states with spatial localization at wire ends
- Near-zero eigenvalue from numerical noise (not localized at edges): not counted as topological

---

### U7: BdG Hermiticity verification test

**Requirements:** R1, R25
**Files:**
- `tests/unit/test_bdg.pf`
**Risk:** low — test only, no Fortran source changes

R1 was confirmed as a false positive through research. The BdG hole block construction `-H0^T` is mathematically correct because H0 stores the full matrix with both (i,j) and (j,i) entries. When H0 is Hermitian, `-H0^T` is also Hermitian (since `-H0^T = -(H0^T) = -(H0*) = -conjg(H0)`, and for Hermitian H0, `conjg(H0) = H0^T`, so the hole block is `-H0^T` which equals `-conjg(H0)` which is Hermitian). Peierls phase is zero for single-column wires (Bx has no orbital effect when nx=1), so no Hermiticity issue arises in practice.

Write a test that builds a multi-column wire BdG Hamiltonian with nonzero Bx, computes `H - H†`, and checks the Frobenius norm is below 1e-12. This confirms Hermiticity for the Peierls-phased case and serves as a regression guard.

**Test scenarios:**
- Build wire BdG with B_vec=[0.5, 0, 0], check `||H - H†|| < 1e-12`
- Build wire BdG with B_vec=[0, 0, 0], check Hermiticity (baseline)
- Single-column wire with nonzero Bx: confirm Peierls phase is zero (no orbital effect)

---

### U8: QW BdG hole block k-dependence

**Requirements:** R2
**Files:**
- `src/physics/bdg_hamiltonian.f90`
**Risk:** medium — changes hole block construction for QW BdG

Fix the QW BdG hole block construction in `build_bdg_hamiltonian_qw` to use the mathematically correct operation for all k-values. The current implementation uses a shortcut that is correct only at k=0. The general operation requires proper handling of the k-dependent terms in the hole block.

Research the exact form: the hole block should be `-H(-k)^T + mu*I` where H(-k) means evaluating the normal-state Hamiltonian at the time-reversed momentum. For the QW case with Peierls phases, this requires conjugation of the hopping terms.

**Test scenarios:**
- QW BdG at k=0: unchanged from current (already correct at k=0)
- QW BdG at nonzero k: hole block is Hermitian
- Full BdG matrix Hermiticity at nonzero k with Peierls phases

---

### U9: Dead code cleanup

**Requirements:** R12, R13, R14, R15, R16
**Files:**
- `src/physics/topological_analysis.f90`
- `src/physics/green_functions.f90`
- `src/core/defs.f90`
- `src/io/input_parser.f90`
**Risk:** low — removing unused code

Remove the following dead code:

- `diag_2x2` private function in `topological_analysis.f90` (zero callers)
- `compute_berry_curvature` pass-through wrapper (zero external callers — the actual computation happens inline)
- `extract_edge_states` wrapper (zero callers, discards 2/3 of return values)
- Two of the three identical conductance wrappers (`compute_hall_conductance_from_chern`, `compute_conductance_kubo_chern`), keeping only `compute_hall_conductance`. Update any call sites and tests.
- `compute_conductance_landauer_from_transmission` identity function in `green_functions.f90` (G=T adds no value). Have callers assign directly.
- Unused `N` parameter from `compute_z2_gap` in `topological_analysis.f90`. Update call sites.
- Dead `bdg_config%self_consistent` field in `defs.f90` and its parser entry in `input_parser.f90`.

Grep for each removed symbol before deletion to confirm zero callers.

**Test scenarios:**
- All existing tests pass after removal
- Build succeeds without undefined reference errors
- `grep -r` for each removed symbol returns no hits outside the removal itself

---

### U10: Missing contiguous attributes

**Requirements:** R18
**Files:**
- `src/physics/bdg_hamiltonian.f90`
- `src/physics/topological_analysis.f90`
- `src/physics/magnetic_field.f90`
**Risk:** low — attribute-only change, no logic change

Add `contiguous` attribute to all assumed-shape hot-path array arguments:

- `bdg_hamiltonian.f90`: `profile_2d(:,:)` parameter
- `topological_analysis.f90`: ~19 assumed-shape arguments in analysis routines (profile arrays, eigenvalue arrays, state arrays)
- `magnetic_field.f90`: assumed-shape arguments where applicable

Do not add `contiguous` to optional or allocatable arguments.

**Test scenarios:**
- Build succeeds with `-std=f2018` (no new warnings)
- Existing tests pass (no behavior change expected)

---

### U11: Python script fixes

**Requirements:** R19, R20, R21, R22, R23, R24
**Files:**
- `scripts/lecture_03_wavefunctions.py`
- `scripts/lecture_04_optical.py`
- `scripts/lecture_07_scsp.py`
- `scripts/lecture_08_spin.py`
- `scripts/lecture_09_exciton.py`
- `scripts/lecture_10_strain.py`
- `scripts/lecture_13_topological.py`
- `scripts/generate_all_figures.py`
**Risk:** low — script-level fixes, no Fortran changes

- **R19:** Replace 16 manual `tempfile.mkdtemp`/`shutil.rmtree` pairs in 7 lecture scripts with `tempfile.TemporaryDirectory` context managers (already the pattern in lectures 00-02, 05-06).
- **R20:** Replace the local numpy trapezoid compat shim in `lecture_07_scsp.py` with the existing `trapz_fn` from `star_helpers`.
- **R21:** Move the inline `import glob` in `lecture_03_wavefunctions.py` (line ~373) to top-level imports.
- **R22:** Make `run_bhz_z2` and `run_bdg_sweep` in `lecture_13_topological.py` print stderr on subprocess failure instead of silently returning None.
- **R23:** Fix `generate_all_figures.py` to use `__file__`-relative paths instead of CWD-relative paths.
- **R24:** Add trailing newline to `generate_all_figures.py`.

**Test scenarios:**
- Each modified script runs without import errors
- `generate_all_figures.py` works from any CWD
- `lecture_13_topological.py` subprocess failures print stderr to console

---

## Execution Sequencing

```
U9 (dead code) ──► U1 (delta_0 guard + pairing constant)
                 ──► U2 (is_z2_transition)
                       ├──► U4 (Fu-Kane Z2)
                       ├──► U5 (phase diagram cumulative Z2)
                       │     └──► U6 (QSHE wire Z2 from spectrum)
                       └──► U3 (regression bias + fitting)
                 ──► U8 (QW BdG hole block)
                 ──► U7 (Hermiticity verification test)
                 ──► U10 (contiguous attributes)
                 ──► U11 (Python scripts)
```

U9 first (shrinks diff surface). U1, U2, U8, U7, U10 can proceed in parallel after U9. U4, U5, U3 depend on U2. U6 depends on U5. U11 is independent throughout.

---

## Test Plan

### Existing regression guard

All 82+ existing tests must pass after each implementation unit. Run `ctest --test-dir build -j4` after each unit.

### New tests (R25-R29)

| Test | Requirement | File | Description |
|------|------------|------|-------------|
| R25 | R1 | `tests/unit/test_bdg.pf` | BdG Hermiticity with nonzero Peierls phase |
| R26 | R3 | `tests/unit/test_bdg.pf` | delta_0=0 rejection with error message |
| R27 | R6 | `tests/unit/test_z2.pf` | Z2 value at phase transition boundary |
| R28 | R5 | `tests/unit/test_z2.pf` | is_z2_transition false when Z2 unchanged |
| R29 | R10 | `tests/unit/test_topology.pf` | fit_exponential_decay accuracy with known xi=5.0 |

---

## Scope Boundaries

- Python type annotations for `star_helpers.py` and lecture scripts
- REPO/BUILD_DIR/CONFIGS_DIR/FIGURES_DIR boilerplate consolidation into shared helper
- BHZ hardcoded model parameters (intentional HgTe/CdTe parameterization)
- Krylov snapshot brittleness mitigation (inherent to the approach)
- LDOS PARDISO symbolic factorization reuse (optimization, not correctness)
- File length refactoring (splitting `topological_analysis.f90` into 4 modules and `main_topology.f90` refactoring) — deferred to a follow-up
- New physics features beyond the 29 findings (no scope expansion)

---

## Dependencies / Assumptions

- The 5 index-arithmetic bug fixes from `docs/solutions/logic-errors/topological-magnetic-index-logic-errors-2026-05-08.md` must not be regressed
- The FEAST sentinel value fix (min_gap = -1.0 on solver failure) must survive changes to `main_topology.f90`
- Peierls phase is zero for single-column wires (Bx has no orbital effect when nx=1) — confirmed by research
- The existing test infrastructure (pFUnit, ctest labels) supports all new tests
- `csr_build_from_coo` sums duplicate (row,col) entries — confirmed in `sparse_matrices.f90`

---

## Risks

| Risk | Impact | Mitigation |
|------|--------|------------|
| U4 (Fu-Kane) changes Z2 for existing BHZ results | Regression in golden outputs | Run verification tests before/after; compare Z2 at known parameter points |
| U5 (cumulative Z2) changes phase diagram output | Existing phase diagram tests may need update | Compare old vs new output for a representative sweep |
| U8 (QW BdG hole block) may not be k-broken | Research suggests the shortcut may work for QW (no Peierls) | Test at multiple k-points before changing code |
| U3 (regression bias) changes Majorana xi values | Golden outputs for Majorana localization may shift | Tolerance-based comparison, not exact match |

---

## Review Checkpoints

- After U9: verify all tests pass with dead code removed
- After U1-U2: verify delta_0 guard and Z2 transition logic with new tests
- After U3-U6: full test suite + verification tests for topological invariants
- After U7-U8: BdG Hermiticity confirmed, QW hole block verified
- After U10-U11: build clean, all scripts run, full test suite green
