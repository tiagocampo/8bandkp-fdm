---
date: 2026-05-12
topic: code-review-fixes
---

# Code Review Fixes for Topological/Superconducting Branch

## Summary

Fix all 29 findings from the code review of `feature/bdg-topological-superconductivity`: BdG Hamiltonian Hermiticity, topological invariant correctness (Z2, Fu-Kane, phase diagrams), numerical regression bias, dead code removal, missing Fortran attributes, and Python script robustness — with tests for each fix.

---

## Problem Frame

The branch adds topological invariant computation (Chern, Z2, Majorana via BdG) across 206 commits and ~11K lines of Fortran source/test changes. A multi-agent code review found 29 issues: 14 P1, 14 P2, and 1 P3. Several P1 findings affect physics correctness — the BdG Hamiltonian may lose Hermiticity when Peierls phases are applied, the Z2 invariant is tracked incorrectly in phase diagrams, and Majorana detection produces false positives when `delta_0=0`. The branch's purpose is new physics; these correctness issues must be resolved for the physics to be trustworthy.

---

## Requirements

**BdG Hamiltonian correctness**

- R1. The BdG Hamiltonian must remain Hermitian when Peierls phases are applied to multi-column wires (Bx != 0). The hole block (2,2) construction must correctly handle complex Peierls-phased entries from the pre-merge COO.
- R2. The QW BdG hole block construction must use the mathematically correct operation for all k-values, not just k=0 where the current `-conjg(H(-k))` shortcut happens to be correct.
- R3. Building a BdG Hamiltonian with `delta_0 <= 0` must be rejected with a clear error message, not silently produce a degenerate matrix that triggers false Majorana detection.
- R4. The pairing sign array (antidiagonal pattern) must be defined once as a module-level constant, not duplicated across `build_bdg_hamiltonian_1d` and `_qw`.

**Topological invariant correctness**

- R5. `is_z2_transition` must report transitions only when Z2 actually changes between adjacent points, not when the gap is small but Z2 is unchanged.
- R6. The Fu-Kane Z2 computation must not assume Kramers partners are ordered consecutively by zheev. Band pairing must be robust against near-degenerate bands at TRIM points.
- R7. Phase diagram Z2 must be tracked cumulatively along the B-sweep: start at Z2=0 (trivial at B=0) and flip each time the gap closes. Per-point gap-closing flags must not replace cumulative tracking.
- R8. The QSHE wire Z2 invariant must be determined from the actual eigenspectrum (gap closing/reopening or parity analysis), not from a hardcoded wire-width threshold.
- R9. `compute_z2_gap` must distinguish between topological edge states and numerical noise near zero energy, not count all eigenvalues within a threshold as evidence of nontrivial topology.

**Numerical correctness**

- R10. The exponential decay linear regression in `fit_exponential_decay` and `compute_majorana_profile` must use the actual count of data points contributing to the sums, not the full tail range count. Zero-density points skipped during accumulation must not inflate `n_fit`.
- R11. The Majorana localization fit must handle the near-transition regime (xi comparable to wire length) gracefully, not silently return a sentinel value that is indistinguishable from failure.

**Dead code and redundant wrappers**

- R12. Remove `diag_2x2` (private, zero callers), `compute_berry_curvature` pass-through wrapper (zero external callers), `extract_edge_states` wrapper (zero callers, discards 2/3 of return values).
- R13. Remove two of the three identical conductance wrappers (`compute_hall_conductance_from_chern`, `compute_conductance_kubo_chern`), keeping only `compute_hall_conductance`. Update call sites and tests.
- R14. Remove the identity function `compute_conductance_landauer_from_transmission` (G=T adds no value). Have callers assign directly.
- R15. Remove the unused `N` parameter from `compute_z2_gap`. Update call sites.
- R16. Remove the dead `bdg_config%self_consistent` field and its parser entry.

**Duplicated logic extraction**

- R17. Extract the tail-search and log-linear fit logic shared between `fit_exponential_decay` and `compute_majorana_profile` into a shared private helper.

**Missing Fortran attributes**

- R18. Add `contiguous` attribute to all assumed-shape hot-path array arguments in `bdg_hamiltonian.f90` (`profile_2d`), `topological_analysis.f90` (19 arguments in analysis routines), and `magnetic_field.f90` where applicable.

**Python script fixes**

- R19. Replace the 16 manual `tempfile.mkdtemp`/`shutil.rmtree` pairs in 7 lecture scripts with `tempfile.TemporaryDirectory` context managers (already the pattern in lectures 00-02, 05-06).
- R20. Replace the local numpy trapezoid compat shim in `lecture_07_scsp.py` with the existing `trapz_fn` from `star_helpers`.
- R21. Move the inline `import glob` in `lecture_03_wavefunctions.py` to top-level imports.
- R22. Make `run_bhz_z2` and `run_bdg_sweep` in `lecture_13_topological.py` print stderr on subprocess failure instead of silently returning None.
- R23. Fix `generate_all_figures.py` to use `__file__`-relative paths instead of CWD-relative paths.
- R24. Add trailing newline to `generate_all_figures.py`.

**Testing**

- R25. Add a test verifying BdG Hamiltonian Hermiticity when Peierls phases are nonzero (build wire BdG with Bx != 0, check H - H† has zero norm).
- R26. Add a test for `delta_0 = 0` rejection (confirm error message, no false Majoranas).
- R27. Add a test verifying Z2 value (not just allocation) at the phase transition boundary in `test_z2_invariant.pf`.
- R28. Add a test for `is_z2_transition` with matching Z2 but small gap (the false-positive case must now return false).
- R29. Add a test for `fit_exponential_decay` accuracy: construct a known exponential with xi=5.0 and verify the recovered xi is within tolerance.

---

## Acceptance Examples

- AE1. **Covers R1, R25.** Given a multi-column wire with B_vec=[0.5, 0, 0], when the BdG Hamiltonian is built, the matrix satisfies H - H† = 0 (norm < 1e-12).
- AE2. **Covers R3, R26.** Given delta_0=0 in the BdG config, when run_bdg_wire is called, it prints an error and stops without attempting Majorana detection.
- AE3. **Covers R5, R28.** Given two adjacent phase diagram points with Z2=0, Z2=0 but gap_a < threshold, `is_z2_transition` returns false.
- AE4. **Covers R7.** Given a BHZ model swept from B=0 to B=5 with one gap closing at B=2, the phase diagram shows Z2=0 for B<2 and Z2=1 for B>2 (cumulative flip, not per-point).
- AE5. **Covers R10, R29.** Given a synthetic exponential density with xi=5.0 and some zero-valued tail entries, `fit_exponential_decay` recovers xi within 10% of the true value.

---

## Success Criteria

- All 82+ existing tests continue to pass (no regressions).
- New tests for R25-R29 pass.
- The code review findings P1-#1 through P1-#14 and P2-#15 through P2-#28 are resolved (fix verified or documented with rationale).
- `compute_chern_qwz` returns correct integer Chern numbers for QWZ model at u=0, ±1, ±2.
- BdG wire produces Majorana zero modes at phase boundary, not in the trivial phase.

---

## Scope Boundaries

- Python type annotations for `star_helpers.py` and lecture scripts
- REPO/BUILD_DIR/CONFIGS_DIR/FIGURES_DIR boilerplate consolidation into shared helper
- BHZ hardcoded model parameters (intentional HgTe/CdTe parameterization)
- Krylov snapshot brittleness mitigation (inherent to the approach)
- LDOS PARDISO symbolic factorization reuse (optimization, not correctness)
- File length refactoring (splitting `topological_analysis.f90` into 4 modules and `main_topology.f90` refactoring) — deferred to a follow-up to avoid disrupting active physics work

---

## Key Decisions

- **Physics fixes over documentation:** Advisory findings about Peierls limitations, QW BdG fragility, and hardcoded Z2 are treated as correctness fixes (the branch's purpose IS this physics), not documentation items.
- **File splits deferred:** The 1287-line and 1234-line file violations are real but the refactoring risk is high mid-feature. Splitting will be a separate task after the physics fixes stabilize.

---

## Dependencies / Assumptions

- The Peierls Hermiticity issue (R1) may be a false positive for single-column wires (the current primary use case) but needs a test to confirm.
- The existing 5 index-arithmetic bug fixes from `docs/solutions/logic-errors/topological-magnetic-index-logic-errors-2026-05-08.md` must not be regressed by the new fixes.
- The FEAST sentinel value fix (min_gap = -1.0 on solver failure) must survive the changes to `main_topology.f90`.

---

## Outstanding Questions

### Deferred to Planning

- [Affects R1][Needs research] Whether the Peierls phase requires conjugation in the hole block depends on whether the COO entries represent pre-merge or post-merge values. The planner should verify by reading the CSR merge path in `sparse_matrices.f90`.
- [Affects R6][Technical] The correct Fu-Kane band pairing strategy for near-degenerate bands needs research: is parity-based individual computation sufficient, or does the Kramers-pair structure need explicit reconstruction?
- [Affects R8][Technical] Whether the BHZ wire Z2 can use the Pfaffian method or should use gap-closing counting with cumulative tracking.
