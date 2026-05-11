---
date: 2026-05-08
topic: 8band-verification-ladder
---

# 8-Band Verification Ladder

## Summary

A 4-rung validation hierarchy for the 8-band zinc-blende k.p Hamiltonian, testing 4 canonical materials (GaAs, InAs, InSb, GaSb/AlSb) at increasing complexity — bulk k=0 structural checks, bulk dispersion effective masses, QW subband energies, and wire internal consistency — using a hybrid of regression tests and a Python analytical validation script with tiered tolerances.

---

## Problem Frame

The codebase has 24 regression tests that compare output against golden files, but the coverage is uneven and undocumented. When a test fails, there is no diagnostic hierarchy — the developer must manually trace whether the error is in parameters, Hamiltonian construction, FD stencils, confinement, or sparse assembly. The benchmarks.md reference doc lists 6 validation regimes but lacks a systematic progression from simplest to most complex physics.

A physics code without layered validation cannot distinguish "the solver is correct but the model has limitations" from "something is broken at an unknown layer." The 20-30% g-factor shortfall vs experiment is attributed to the 8-band model — but without independent validation at each layer of complexity, that attribution is unverified. The ladder creates a "ladder of trust" where each rung is independently confirmed, making residual errors unambiguously attributable to the model rather than the implementation.

---

## Requirements

**Rung 1: Bulk k=0 structural validation**

- R1. For each canonical material, compute the 8x8 bulk Hamiltonian at k=0 and verify eigenvalues match the band edges (EV, EC, Delta_SO) from `parameters.f90` within machine precision.
- R2. Verify eigenvalue degeneracies are correct: HH+LH (4x degenerate at k=0), SO (2x), CB (2x) for the spin-degenerate basis.
- R3. Verify basis ordering: bands 1-4 are valence (HH, LH, LH, HH), bands 5-6 are split-off, bands 7-8 are conduction — confirmed by eigenfunction inspection.
- R4. Verify T_d symmetry at k=0: all off-diagonal k-dependent blocks (Q, R, S, T) vanish.
- R5. Verify eigenfunction normalization (unit sum of squares).

**Rung 2: Bulk dispersion — effective masses**

- R6. For each canonical material, compute the bulk band structure along at least one k-direction with sufficiently small k-step, then extract effective masses (m*_e, m*_hh, m*_lh, m*_so) by parabolic fitting near Gamma.
- R7. Compute the expected effective masses analytically from the code's own Kane parameters (P, Eg, Delta_SO, gamma1/2/3) using the standard 8-band formulas, and verify the numerically extracted masses agree within 1e-8.
- R8. Compare the analytically predicted masses against the Vurgaftman 2001 tabulated values for the same material, and report the relative deviation.
- R9. Flag any material where the self-consistency check (R7) passes but the literature comparison (R8) deviates by more than 0.1% — this indicates a parameter-set inconsistency, not a solver bug.

**Rung 3: QW subband energies**

- R10. For a GaAs/AlGaAs single QW with sufficiently high barriers, compute subband energies and compare against the Bastard infinite-barrier analytical formula E_n = n^2 pi^2 hbar^2 / (2m* L^2) within the expected accuracy for finite barriers.
- R11. For the same GaAs/AlGaAs QW with realistic barrier height, compare subband energies against published nextnano (or equivalent) reference data, agreeing within 0.5 meV.
- R12. For an InAs/GaSb broken-gap QW, compute the band overlap and subband structure, comparing against the existing benchmarks.md validation (142 meV overlap) and published results.
- R13. Verify that QW subband degeneracies and ordering are physically correct (CB subbands, HH/LH subbands, SO subbands).

**Rung 4: Wire internal consistency**

- R14. For a GaAs rectangular wire with large transverse dimensions, verify that wire eigenvalues converge toward the QW result for the same material as the wire width/height increase.
- R15. For the same wire geometry, run with both dense (LAPACK) and sparse (CSR/FEAST) solver paths and verify eigenvalues agree within solver tolerance.
- R16. Verify that wire eigenvalue count matches the expected number given the grid size and numcb/numvb settings.

**Cross-cutting**

- R17. All 4 canonical materials (GaAs, InAs, InSb, GaSb/AlSb) are tested at rungs 1 and 2 (bulk). Material-geometry mapping for rungs 3-4: GaAs/AlGaAs for R10-R11 (standard QW), InAs/GaSbW for R12 (broken-gap QW, using W-variant GaSbW since GaSb lacks EV/EC), GaAs for R14-R16 (wire). Bulk rungs use both W and non-W variants where applicable.
- R18. Tolerances are specified per-requirement: machine precision for rung 1 (exact by construction), 1e-8 for rung 2 (analytical fits), and per-requirement tolerances for rungs 3-4: R10 (TBD based on validation approach), R11 (0.5 meV absolute), R12 (1% relative from benchmarks.md), R14 (1% relative), R15 (1e-10 solver agreement).
- R19. Validation results are documented in `docs/reference/benchmarks.md` alongside the existing 6 regimes, with config file names, expected values, and tolerances.

---

## Acceptance Examples

- AE1. **Covers R1, R2, R3, R4, R5.** Given GaAs parameters from `parameters.f90`, when the bulk Hamiltonian is computed at k=0, the 8 eigenvalues are [EV, EV, EV, EV, EV-Delta_SO, EV-Delta_SO, EC, EC] with degeneracies HH+LH(4x degenerate at k=0), SO(2x), CB(2x) in the spin-resolved basis, off-diagonal blocks are zero, and eigenfunctions have unit norm.
- AE2. **Covers R6, R7, R8.** Given GaAs dispersion computed at k = 0.001 nm^-1 steps near Gamma, the parabolic fit yields m*_e = 0.067 m0, which matches the analytical prediction from (P, Eg, Delta_SO) within 1e-8, and matches Vurgaftman Table III within 0.1%.
- AE3. **Covers R10.** Given a 10 nm GaAs QW with Al0.3Ga0.7As barriers at FD order 2 with 200 grid points, the E1-E2 CB subband spacing agrees with the Bastard analytical result within the expected finite-barrier correction (~5%).
- AE4. **Covers R14, R15.** Given a 100x100 nm GaAs rectangular wire computed with 21x21 grid, the ground-state eigenvalue is within 1% of the QW ground state for the same material; dense and sparse solver paths agree within 1e-10.

---

## Success Criteria

- All 4 rungs pass for all targeted materials with no tolerance violations.
- When a regression is introduced (e.g., wrong sign on a Bir-Pikus term), the ladder identifies the failing rung, narrowing the diagnosis to a specific physics layer.
- The validation results are documented in benchmarks.md with sufficient detail that an independent researcher can reproduce them.
- A downstream planner can implement this without inventing physics decisions — all reference values, tolerances, and material choices are specified.

---

## Scope Boundaries

- Temperature-dependent validation (all tests at default temperature)
- Optical properties validation (absorption, gain, ISBT — separate ideation item)
- g-factor validation via Lowdin partitioning (covered by Standard-Star Benchmarks)
- Topological invariants (Chern, Z2 — separate module, not on this ladder)
- Strain validation (covered by ideation item #7 — Strain Validation Across Geometries)
- Richardson convergence testing (covered by ideation item #6)
- CSR structure testing (covered by ideation item #2)
- Expanding to all 25+ materials (start with 4 canonical, expand later)

---

## Key Decisions

- **4 canonical materials, not all 25+:** GaAs, InAs, InSb, GaSb/AlSb cover the range from large gap to narrow gap to broken-gap. Deep multi-property validation per material is more valuable than shallow validation across all materials.
- **No intermediate band-count models:** The ladder validates the 8-band model only — no 2-band Kane or 4-band Luttinger rungs. The 8-band is the core of this codebase.
- **Wire rung uses internal consistency only:** No external published wire references — validate via wire→QW limit and dense↔sparse cross-check.
- **Hybrid test format:** Regression tests for golden-output pass/fail (rungs 1, 3, 4) plus Python script for analytical comparisons that compute references on-the-fly (rung 2).
- **QW rung covers three sub-cases:** Infinite-barrier Bastard, finite-barrier vs nextnano, and InAs/GaSb broken-gap.

---

## Dependencies / Assumptions

- Vurgaftman 2001 tabulated effective masses and band parameters are the authoritative literature reference.
- The Bastard 1981 infinite-barrier formula provides a sufficiently accurate analytical reference for the QW infinite-barrier sub-case.
- nextnano reference data for GaAs/AlGaAs QW is available (either from benchmarks.md cross-validation or from published nextnano output).
- The Python validation script can parse the code's output format (space-delimited ASCII, g14.6 formatting).
- Parabolic fitting near Gamma provides reliable effective mass extraction for the 8-band dispersion.

---

## Outstanding Questions

### Resolve Before Planning

None — all scope decisions were resolved in dialogue.

### Deferred to Planning

- [Affects R6][Technical] What k-range and k-step size gives reliable parabolic fitting for each material's effective mass extraction?
- [Affects R11][Needs research] Where is the nextnano reference data for GaAs/AlGaAs QW — in the existing benchmarks.md cross-validation, or does it need to be generated?
- [Affects R14][Technical] What wire dimensions constitute a "large" wire that approaches the QW limit? How many grid points are needed for convergence?
- [Affects R19][Technical] Should the Python validation script produce a summary report, or should it integrate with the existing `compare_output.py` framework?

## Deferred / Open Questions

### From 2026-05-08 review

- **Bastard formula comparison invalid for 8-band QW** — Rung 3 / R10, AE3 (P0, all 5 reviewers, confidence 100)

  R10 and AE3 compare 8-band QW subband energies against the Bastard infinite-barrier formula, claiming ~5% agreement. The Bastard formula is single-band; benchmarks.md documents a 5.4x discrepancy (302 meV vs 56.2 meV). This comparison cannot work as stated. Needs a replacement validation approach: finite-barrier 1-band model, published 8-band QW results, or nextnano cross-validation.

  <!-- dedup-key: section="rung 3 r10 ae3" title="bastard formula comparison invalid for 8band qw" evidence="r10 and ae3 compare 8band qw subband energies against the bastard infinitebarrier formula claiming 5 agreement" -->

- **Analytical effective mass formula undefined for 8-band** — Rung 2 / R7 (P0, feasibility, confidence 100)

  R7 claims analytical effective masses from Kane parameters at 1e-8 tolerance, but no simple closed-form 8-band formula achieves this. The code's CB diagonal uses A=1/meff PLUS explicit P coupling, so total CB curvature is non-trivial. Numerical mass for GaAs ~0.046 m0 vs Vurgaftman 0.067 vs Kane formula 0.053. Needs reworked validation: numerical extraction from dispersion, self-consistent analytical prediction, or direct Vurgaftman comparison with relaxed tolerance.

  <!-- dedup-key: section="rung 2 r7" title="analytical effective mass formula undefined for 8band" evidence="r7 claims analytical effective masses from kane parameters at 1e8 tolerance but no simple closedform 8band formula achieves this" -->

- **R9 deviation threshold too tight at 0.1%** — Rung 2 / R9 (P2, scope-guardian, adversarial, confidence 100)

  R9 flags deviations >0.1% between analytical and Vurgaftman tabulated values. This is unrealistically tight — 8-band model parameters don't exactly reproduce experimental tabulated values. Needs relaxed threshold (1-5%) and clarified purpose as a parameter-consistency check, not a solver-correctness test.

  <!-- dedup-key: section="rung 2 r9" title="r9 deviation threshold too tight at 01" evidence="r9 flags deviations 01 between analytical and vurgaftman tabulated values this is unrealistically tight" -->

- **Problem frame g-factor motivation disconnect** — Problem Frame (P2, product-lens, adversarial, confidence 100)

  The problem frame motivates the ladder by referencing the 20-30% g-factor shortfall, but the ladder (rungs 1-4) doesn't validate g-factors at all. This attribution is unverified without the ladder covering g-factors. Needs either: remove g-factor from motivation and replace with general validation-gap argument, or note that g-factor validation is out of scope and referenced as context only.

  <!-- dedup-key: section="problem frame" title="problem frame gfactor motivation disconnect" evidence="the problem frame motivates the ladder by referencing the 2030 gfactor shortfall but the ladder rungs 14 doesnt validate gfactors at all" -->
