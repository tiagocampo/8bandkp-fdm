---
date: 2026-05-16
topic: figure-validation-fixes
---

# Figure Validation: Remaining Fixes

## Summary

Fix 6 HIGH/STALE and 10 MEDIUM-severity issues found during the figure validation audit of 104 PNGs against 15 lecture markdowns. Fixes span Fortran optics code (2 bugs), Python figure-generation scripts (4 bugs + 1 gap), and lecture markdown text (10 mismatches).

---

## Problem Frame

Phases 1-2 of the figure validation audit classified all 104 figures by severity. Phases 3-4 partially completed: 2 strain figures regenerated, 6 text edits in lecture 04, 5 orphaned figures removed. Investigation revealed clear root causes for all remaining issues. The codebase currently ships figures with physically incorrect content (bulk TE≠TM, ISBT reading stale data, scattering lifetimes off by orders of magnitude) and lecture text that contradicts its own figures.

---

## Requirements

**Fortran code fixes**

- R1. Bulk absorption TE=TM: for confinement=0 (bulk), `optics_accumulate` must treat all three velocity components identically, producing equal TE and TM spectra. Currently it uses `px+py` vs `pz` with eigenstates from a z-only k-sweep, which violates cubic isotropy.
- R2. Bulk absorption magnitude: the 3D k-space weight in `main_optics.f90` must include proper `(2*pi)^3` BZ normalization. Currently missing a `(2*pi)^2` factor, causing ~40x underestimation.
- R3. Wire `compute_intersubband_transitions` into `main_optics.f90` so that `isbt_transitions.dat` is produced when ISBT mode is active, unblocking `fig_isbt_dipole_moments`.

**Python script fixes**

- R4. `fig_isbt_absorption` and `fig_gain_strained_comparison` both call `EXE_BAND` (bandStructure) instead of `EXE_OPTICS` (opticalProperties). ISBT and gain computation only exist in the optics executable. Both must use `EXE_OPTICS`.
- R5. `fig_scattering_lifetime_vs_width` has wrong column mapping: reads col 1 (subband index j) as energy, col 2 (transition energy in meV) as rate, then inverts it for lifetime. Must read cols 3-4 (rates in 1/s) and use cols 5-6 (lifetimes in ps) directly, summing emission + absorption contributions.
- R6. `fig_scattering_lifetime_vs_field` has correct column mapping but relies on stale cached data from before the Fortran scattering-unit fix (commit c160eff). The sweep must be re-run to produce physically correct magnitudes and trends.
- R7. `fig_qw_absorption_vs_width` must implement a width-sweep loop following the `fig_exciton_binding_vs_width` pattern: loop over well widths, write config per width, run bandStructure, read absorption output, cache results, plot multi-curve comparison.
- R8. `fig_qw_absorption_strained` config (`qw_ingaas_algaas_strained_absorption.cfg`) must use a material system where strain effects are significant. The current InGaAs/AlInAs on InP is lattice-matched (near-zero strain). Switch to GaAs/InGaAs where the InGaAs well is compressively strained.

**Lecture text fixes (MEDIUM severity)**

- R9. `docs/lecture/03-wavefunctions.md`: `vb_hh_lh_mixing` text claims majority-LH crossover not visible in figure — update text to describe what the figure actually shows.
- R10. `docs/lecture/03-wavefunctions.md`: `cb_parts_evolution` shows CB purity ~95% but text says 67% — align text with figure.
- R11. `docs/lecture/09-numerical-methods.md` / `docs/lecture/11-convergence.md`: `convergence_fd_order` order-2 error value in text disagrees with figure — fix text to match computed value.
- R12. `docs/lecture/09-numerical-methods.md`: `timing_dense_vs_sparse` wire timing is 2.3x larger than text claims — update text.
- R13. `docs/lecture/06-optical-properties.md`: `absorption_excitonic_TE` magnitudes too low relative to text expectations — investigate and align.
- R14. `docs/lecture/06-optical-properties.md`: `isbt_dipole_moments` shows oscillator strength values that violate Thomas-Reiche-Kuhn sum rule — investigate and fix figure or text.
- R15. `docs/lecture/06-optical-properties.md`: `exciton_bohr_vs_width` shows Bohr radius exceeding bulk for wide wells — verify physics and update text if figure is correct.
- R16. `docs/lecture/07-self-consistent-sp.md`: `sc_delta_doped_potential` shows 9 meV notch but text says 100 meV — update text to match figure.
- R17. `docs/lecture/05-gfactor.md`: `qw_gfactor_vs_width` data range is 20-80 A but text says 10-120 A — update text.
- R18. `docs/lecture/08-quantum-wire.md`: `wire_inas_gaas_subbands` limited display (2 kz points) — update text to acknowledge limited k-resolution.

---

## Success Criteria

- All HIGH/STALE figures regenerated with correct physics (TE=TM for bulk, ISBT in correct energy window, scattering lifetimes in 1-10 ps range, strained absorption showing visible HH/LH splitting, width sweep showing multiple curves)
- All 15 lecture markdowns have text consistent with their referenced figures
- `OMP_NUM_THREADS=12 python3 scripts/plotting/generate_all_figures.py` completes without errors
- Regression test suite (`ctest --test-dir build`) still passes with no new failures

---

## Scope Boundaries

- No changes to material parameters in `parameters.f90`
- No changes to basis ordering (bands 1-4 valence, 5-6 split-off, 7-8 conduction)
- No refactoring bulk optics to full 3D spherical k-integration (isotropic averaging is sufficient)
- No new test cases for fixed figures
- No changes to the input parser's sequential reading behavior (document the limitation)

---

## Key Decisions

- Bulk TE=TM fix uses isotropic averaging `(px+py+pz)` for both TE and TM when confinement=0, rather than implementing directional k-sweep integration. Simpler, correct, and sufficient.
- Strained absorption switches from InGaAs/AlInAs-on-InP (lattice-matched) to GaAs/InGaAs system where strain effects are physically significant.
- Width sweep follows the existing `fig_exciton_binding_vs_width` caching pattern (write once, reuse on subsequent runs).
- MEDIUM-severity text fixes update lecture markdown to match figures rather than regenerating figures to match text (figures are computed from the code; text is the thing that can be wrong).

---

## Dependencies / Assumptions

- The Fortran scattering-unit fix (commit c160eff) is correct; re-running the field sweep will produce physically valid lifetimes.
- The `compute_intersubband_transitions` subroutine in `optical_spectra.f90` is functionally complete and only needs to be called from the right place.
- The strained absorption config change (GaAs/InGaAs) will require verifying that GaAs and InGaAs material parameters exist in `parameters.f90` for the chosen composition.
