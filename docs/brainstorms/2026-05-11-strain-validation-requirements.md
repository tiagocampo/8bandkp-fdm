---
date: 2026-05-11
topic: strain-validation
---

# Strain Validation Across Geometries

## Summary

End-to-end strain validation suite for InAs/GaAs across bulk, QW, and wire geometries. Validates the full chain (strain computation → Bir-Pikus → Hamiltonian → eigenvalues) against published analytical formulas and reference data. Bulk and QW tests extend the verification ladder; wire tests get a dedicated `strain-validation` ctest label.

---

## Problem Frame

The Bir-Pikus strain implementation is one of the most fragile code paths in the codebase. The sign convention (`P_eps = -av * Tr(eps)`) is a known hazard — CLAUDE.md flags it as a "NEVER change" boundary. A gfortran -O3 segfault previously hit the strained QW path, confirming the code is undertested under optimization. The existing 10 unit tests validate Bir-Pikus formulas and COO assembly mechanics in isolation, but no test verifies that the end-to-end physics (strained eigenvalues) matches published results.

The gap is integration-level: if a sign error, unit mismatch, or Hamiltonian insertion bug exists in the strain path, the unit tests pass but every strained calculation is silently wrong. This affects band structure, g-factors, and optical properties — all downstream of the strained Hamiltonian.

Published reference data is available for all three geometries: Vurgaftman 2001 deformation potentials for analytical Bir-Pikus shifts (used for bulk validation in R1), nextnano published QW transition energies, and analytical plane-strain solutions for wire strain profiles.

---

## Requirements

**Bulk strain validation**

- R1. Strained bulk eigenvalues at k=0 match analytical Bir-Pikus band edge shifts (delta_Ec, delta_EHH, delta_ELH, delta_ESO) computed from Vurgaftman deformation potentials for InAs biaxially strained on GaAs substrate, within 1% tolerance. The existing `lecture_04_strain.py` Section 2 performs a similar validation for strained bulk GaAs; this requirement extends it to InAs-on-GaAs and promotes the result to the verification ladder with ctest integration. A new `tests/regression/configs/bulk_inas_gaas_strained_k0.cfg` is needed (material1: InAs, confinement=0, strain: F, strainSubstrate: 5.65325 for GaAs lattice constant, no external field). Uses the `strainSubstrate` mechanism rather than the `strain:` block (bulk bypasses `compute_strain` entirely).
- R2. HH-LH splitting in strained bulk matches the analytical shear deformation formula `b * (eps_zz - eps_xx)` (equivalently `2 * Q_eps` where `Q_eps = b/2 * (eps_zz - eps_xx)` for biaxial strain) using Vurgaftman `b` parameter for InAs, within 1% tolerance.
- R3. The strained-minus-unstrained eigenvalue difference at k=0 equals the analytical Bir-Pikus shift exactly (to numerical precision), confirming strain is a pure additive modification. Requires two separate bulk InAs runs (unstrained: strainSubstrate=0, strained: strainSubstrate=GaAs a0), following the pattern in `lecture_04_strain.py` Section 2.

**QW strain validation**

- R4a. Strained InAs/GaAs QW subband energies at k=0 (ground-state electron, ground-state HH, ground-state LH) match Bastard/Winkler analytical approximations, within 3-5% tolerance accounting for analytical model limitations (infinite barrier, parabolic bands).
- R4b. Alternatively, if nextnano numerical reference data is available for a comparable InAs/GaAs QW config, strained subband energies match within ~1% tolerance (covering only numerical discretization differences between the two FD implementations). The planner should evaluate both options and commit to one before implementation.
- R5. HH-LH splitting behavior in the strained QW is validated against analytical expectations. The HH-LH ordering under compressive strain must be verified against `compute_bp_scalar` during planning (see Outstanding Questions). The splitting magnitude matches the shear deformation prediction within 3% tolerance.
- R6. Strained QW results are insensitive to a 2x increase in grid density (convergence check), confirming the strain+FD combination is numerically stable.

**Wire strain validation**

- R7. Navier-Cauchy PDE strain profile (eps_xx, eps_yy, eps_zz components) for an InAs square core in GaAs shell matches the analytical plane-strain solution at representative cross-section points, within 5% tolerance for strain magnitude and 10% for profile shape (decay rate from core to boundary). Uses square inclusion geometry to avoid staircase-interface artifacts inherent in Cartesian grid discretization of circular boundaries. The planner should derive or cite the closed-form analytical expressions for the square-inclusion strain field.
- R8. Hydrostatic strain (Tr(eps)) is negative (compressive) in the InAs core and decays toward zero in the GaAs shell, matching the physical expectation for lattice-mismatch-driven strain.
- R9. Band edge positions in the strained wire cross-section show quantitatively correct behavior: CB shift in the InAs core is positive (matching the Bir-Pikus hydrostatic prediction within 20%), and HH-LH splitting in the core exceeds a minimum threshold derived from the R2 shear deformation formula scaled to the wire geometry.

**Lecture and documentation**

- R10. Audit `scripts/lecture_04_strain.py` and its companion markdown doc against the validated R1-R9 numerical values. Update any hardcoded constants, expected values, or assertions that differ from the validated results. If no discrepancies are found, record that fact.

**Test infrastructure**

- R11. Bulk (R1-R3) and QW (R4a/R4b-R6) tests are integrated into the verification ladder with appropriate rung numbering and ctest labels. These complement (not duplicate) the existing S6 standard-star strained QW benchmark, which validates HH-LH splitting and g-factor against published references using `qw_inas_gaas_strained.cfg`. R4a uses analytical references; S6 uses regression golden output — they validate different aspects.
- R12. Wire tests (R7-R9) use a dedicated `strain-validation` ctest label, separate from the verification ladder.
- R13. All new tests include `# COVERAGE:` annotations following the established convention, tagging observable, geometry, material, and reference source.
- R14. If validation exposes bugs in the strain solver (wrong signs, unit mismatches, incorrect Hamiltonian insertion), the underlying code is fixed — not just the tests.

---

## Success Criteria

- Running `ctest -L verification` includes strained bulk and QW checks that pass against published reference data.
- Running `ctest -L strain-validation` includes wire strain checks that pass against analytical strain profiles.
- The Bir-Pikus sign convention (`P_eps = -av * Tr(eps)`) is validated end-to-end for the first time, confirming the "NEVER change" boundary is correct in practice.
- The strain lecture script produces output consistent with the validated numbers.
- Any strain solver bugs discovered during validation are fixed and covered by new or updated unit tests.

---

## Scope Boundaries

- Additional material systems beyond InAs/GaAs (GaInAs/InP, GaSb/AlSb, InSb) — follow-up after the framework is proven
- Multi-property strain cross-checks (g-factor shift under strain, optical absorption edge shift) — follow-up once strained eigenvalues are validated
- Piezoelectric effects (orthogonal to Bir-Pikus, already has a config flag)
- Modifications to the wire PDE solver algorithm itself (validation may surface bugs that get fixed, but no solver redesign)
- New lecture scripts beyond lecture_04_strain.py

---

## Key Decisions

- **End-to-end eigenvalue chain rather than formula isolation.** The existing 10 unit tests already cover Bir-Pikus formulas and COO assembly in isolation. The new suite validates the integration — strain through to eigenvalues — which is the gap that matters. (Rationale: formula tests passing while end-to-end physics is wrong is exactly the failure mode to prevent.)
- **InAs/GaAs single system, full geometry coverage.** Proves the framework works across all three geometry tiers before adding material diversity. InAs/GaAs has a 6.7% lattice mismatch and is chosen for its rich published reference data (Vurgaftman deformation potentials, Bastard/Winkler analytical formulas, nextnano cross-validation data). (Rationale: best availability of independent reference results for validation, plus significant strain effects.)
- **Hybrid test framework.** Bulk and QW extend the verification ladder because they validate the same Hamiltonian at increasing complexity. Wire gets a dedicated label because it exercises the Navier-Cauchy PDE solver, a distinct subsystem. (Rationale: keeps the verification ladder conceptually coherent while giving wire strain appropriate visibility.)
- **Architectural split in bulk vs QW/wire strain paths.** Bulk strain uses `strainSubstrate` applied inline in `ZB8bandBulk` via `apply_bp_strain_inline`, while QW and wire use the `compute_strain` + `compute_bir_pikus_blocks` pipeline. Both paths call `compute_bp_scalar` for the Bir-Pikus formulas (shared single source of truth), but the integration differs. R1-R3 validate the bulk-specific path; R4-R9 validate the QW/wire pipeline.

---

## Dependencies / Assumptions

- Vurgaftman 2001 deformation potentials for InAs (a_c, a_v, b, d) are correctly implemented in `parameters.f90` — this is a CLAUDE.md boundary ("NEVER modify without verifying against published references") and is assumed correct.
- The analytical plane-strain solution for a square inclusion on a Cartesian grid is a valid reference for the wire strain profile. Assumes elastic constants (c11, c12) in `parameters.f90` match those used in the analytical calculation. Note: the `strain_ref` config field is decorative only — the reference lattice constant is always `params(1)%a0`. All QW and wire strain configs must list the substrate (GaAs) as material1/region1 to produce correct InAs-on-GaAs strain.
- Bastard/Winkler analytical QW subband energies are approximate (infinite barrier, parabolic bands) — the 3% tolerance accounts for the analytical model's limitations, not just numerical error.
- The existing `tests/integration/star_helpers.py` infrastructure (run_exe, parse helpers) can be reused for the new verification scripts.
- Tolerance values (1% bulk, 3-5% QW analytical, ~1% QW numerical, 5% wire) are starting points — actual achievable tolerances will be determined during implementation and may need adjustment based on the analytical model fidelity vs. the 8-band numerical result.
- The existing S6 standard-star benchmark (`verify_star_inas_gaas_qw.py`) validates strained InAs/GaAs QW HH-LH splitting and g-factor against published references using `qw_inas_gaas_strained.cfg`. R4a/R4b tests use analytical or numerical references (not regression golden output), providing orthogonal validation to S6.

---

## Outstanding Questions

### Deferred to Planning

- [Affects R4a/R4b, R6][Needs research] Which specific Bastard or Winkler analytical formulas and parameter sets give the best match for InAs/GaAs QW subband energies? The planner should verify whether Bastard infinite-barrier or Winkler 8-band analytical corrections are more appropriate, then commit to R4a or R4b.
- [Affects R7][Technical] How to compute the analytical plane-strain solution for a square inclusion on a Cartesian grid? The planner should derive or find the closed-form expressions for eps_xx, eps_yy, eps_zz as functions of position, core size, and elastic constants. If a square-inclusion solution is not readily available, consider using a cylindrical Eshelby solution evaluated along the x-axis (where eps_xx = eps_rr) as an approximate reference with appropriately relaxed tolerance. Risk: square-inclusion strain fields are typically non-uniform and may require Fourier series evaluation rather than simple closed-form expressions. Tolerance may need relaxation to 10-15% if the reference solution has numerical integration error. Interior core points should be used for comparison (avoiding boundary cells where first-order FD extraction reduces accuracy).
- [Affects R10][Needs research] What is the current state of `scripts/lecture_04_strain.py` and its companion doc? The planner should audit the existing content for outdated or incorrect values before determining the scope of updates needed.

### From 2026-05-11 Review

- [Affects R5][Physics verification] Verify the HH-LH ordering for compressively strained InAs on GaAs. The Bir-Pikus formulas with InAs parameters (b = -1.8 eV, eps_xx = -0.067, eps_zz = +0.073) give Q_eps = -0.126 eV, which implies HH shifts below LH for in-plane compression (delta_EHH = -P_eps + Q_eps < delta_ELH = -P_eps - Q_eps when Q_eps < 0). If this is correct, R5 needs revision; if the code uses a different sign convention for the VB diagonal terms, document which convention applies.
