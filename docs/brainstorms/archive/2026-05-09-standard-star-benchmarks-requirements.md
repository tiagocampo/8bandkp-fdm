---
date: 2026-05-09
topic: standard-star-benchmarks
---

# Standard-Star Benchmark Systems

## Summary

A suite of 7 per-material Python verification scripts that validate multiple published observables per canonical material system — g-factors, optical absorption edges, strain shifts, and topological invariants that the verification ladder doesn't cover. Each script runs as a ctest target under a `standard-star` label and generates a publication-ready benchmark table with literature references, tolerances, and computed-vs-expected deltas.

---

## Problem Frame

The verification ladder (rungs 1-4) validates solver machinery: bulk eigenvalues, dispersion effective masses, QW subband energies, and wire internal consistency. It covers the "does the solver work?" question but leaves the "does the physics match published results?" question partially answered. G-factors, optical absorption spectra, strain shifts, and topological invariants are currently validated by running executables and manually comparing output against literature — no automated regression, no traceability, no repeatable benchmark table.

This matters because a methods paper needs systematic, reproducible benchmark results. A reviewer or reader should be able to run `ctest --test-dir build -L standard-star` and see every standard material system validated against Vurgaftman 2001, Winkler 2003, Bastard 1981, and other authoritative references. Today, that table doesn't exist in automated form.

---

## Requirements

**Standard-star system definitions**

- R1. Seven standard-star material systems, each representing a distinct physics regime:

  | # | System | Regime | Geometry |
  |---|--------|--------|----------|
  | S1 | GaAs | Standard III-V | Bulk |
  | S2 | InAs | Narrow-gap | Bulk |
  | S3 | InSb | Extreme spin-orbit | Bulk |
  | S4 | GaAs/AlGaAs | Type-I QW | Quantum well |
  | S5 | InAs/GaSb | Type-III broken-gap / topological QW | Quantum well |
  | S6 | InAs/GaAs | Strained QW | Quantum well |
  | S7 | InAs nanowire | Wire geometry | Nanowire |

- R2. Each system validates 2-5 independent published observables (minimum 2 for systems with thin literature coverage). The observable set per system is:

  | System | Eg | ΔSO | m* | g* | Subbands | Overlap | Stark | Optics | Strain | Topo | Wire g* |
  |--------|:--:|:---:|:--:|:--:|:--------:|:-------:|:-----:|:------:|:------:|:----:|:-------:|
  | GaAs bulk | x | x | x | x | | | | x | | | |
  | InAs bulk | x | | x | x | | | | | | | |
  | InSb bulk | x | x | x | x | | | | | | | |
  | GaAs/AlGaAs QW | | | | | x | | x | x | | | |
  | InAs/GaSb QW | | | | x | | x | | | | x | |
  | InAs/GaAs QW | | | | x | x | | | | x | | |
  | InAs wire | | | | x | | | | | | | x |

- R3. Every validated observable has a published literature reference (author, year, source) and a named tolerance regime.

**Observable validation**

- R4. Analytical observables (band gap, spin-orbit splitting, Roth g-factor) validated at exact or near-exact tolerance — determined by construction from the parameter set.

- R5. Kane-model observables (effective masses) validated via parabolic fitting near Gamma. Tolerance reflects the 8-band model's known non-parabolicity deviation from Vurgaftman tabulated masses.

- R6. Confinement observables (subband spacing, band overlap, Stark shift, HH-LH splitting) validated against published analytical or computational references with a stated tolerance per observable.

- R7. Topological observables (band inversion, Z2 invariant) validated against the InAs/GaSb broken-gap QW's topological phase transition and inverted regime.

**Test infrastructure**

- R8. One Python verification script per standard-star system, following the naming convention `verify_star_<material>.py` and the existing `verify_8band_rung*.py` pattern.

- R9. All scripts registered as ctest targets under a `standard-star` label, runnable via `ctest --test-dir build -L standard-star`.

- R10. Scripts share a common utility layer for running executables, parsing output, and comparing against literature values with tolerance. This avoids duplicating boilerplate across 7 scripts.

- R11. Standard-stars complement the verification ladder — organized by material rather than by complexity level, and covering observables the ladder doesn't. Overlapping observables (bulk Eg/ΔSO for GaAs/InAs/InSb from rung 1; effective masses from rung 2; GaAs/AlGaAs subbands from rung 3) are cross-referenced, not duplicated. Each standard-star script should note which observables are already covered by the ladder and link to the relevant rung.

**Publication output**

- R12. Each script outputs a benchmark table row in markdown: `| Material | Observable | Computed | Expected | Reference | Tolerance | Delta | Status |`.

- R13. Running the full standard-star suite produces a complete benchmark table aggregating all systems. An aggregation step (script or make target) concatenates per-system output into a publication-ready table.

- R14. Literature references (author, year, journal/book) are embedded in each script and reproduced verbatim in the benchmark table output.

---

## Acceptance Examples

- AE1. **Covers R4.** For GaAs bulk at k=0, the script verifies Eg = 1.519 eV against Vurgaftman 2001 at machine precision and reports PASS with delta < 1e-14 eV.

- AE2. **Covers R6.** For GaAs/AlGaAs QW, the script verifies CB subband spacing against Bastard 1981 within a stated absolute tolerance, reporting the computed value, expected value, reference, and PASS/FAIL.

- AE3. **Covers R9, R12, R13, R14.** Each standard-star script outputs a benchmark table row. Running the full suite and concatenating output produces a complete markdown benchmark table with material name, observable, computed value, expected literature value, literature citation, tolerance, delta, and PASS/FAIL.

- AE4. **Covers R10.** Adding a new standard-star system requires only a new script that imports the shared utilities — no modification to existing scripts.

---

## Success Criteria

- All 7 standard-star systems pass in CI with automated comparisons against published literature values.
- The benchmark table generated by the standard-star suite (with aggregation) is directly usable in a methods paper.
- G-factors, optical absorption edges, strain shifts, and topological invariants — observables currently lacking automated regression — are covered.
- A reader unfamiliar with the codebase can understand what each standard-star validates and why, from the test output alone.

---

## Scope Boundaries

- Not replacing or refactoring the existing verification ladder
- Not adding Richardson extrapolation convergence testing (ideation idea #6)
- Not adding executable lecture-test pairs (ideation idea #4)
- Not adding a validation coverage matrix (ideation idea #5)
- Not adding cross-code validation pipeline with nextnano (rejected ideation idea #13)
- Not auditing Vurgaftman parameter provenance (rejected ideation idea #4)
- Standard-stars validate against published references; they do not compute new physics or novel predictions
- Performance benchmarking (time-to-solution) is out of scope — standard-stars validate physics accuracy only

---

## Key Decisions

- **Per-material organization (not per-observable):** A reviewer reads results by material system ("does GaAs work?"), not by observable type ("do all g-factors work?"). Publication tables are organized the same way.
- **Dual output (test + publication table):** Making the test output the publication artifact eliminates the divergence risk where the test passes but the published table is stale.
- **Complement, not extend, the ladder:** The ladder's "one rung = one complexity level" principle is orthogonal to standard-stars' "one system = many observables." Merging them would compromise both organizational principles.
- **7 systems for full regime coverage:** GaAs/InAs/InSb span bulk parameter space; three QW types (Type-I, Type-III, strained) span confinement physics; InAs/GaSb QW also covers topological regime (inverted bands, phase transition); InAs wire covers wire geometry. HgTe/CdTe excluded because parameters.f90 lacks II-VI materials.
- **Phased implementation order:** Phase 1 — S1-S3 bulk (validates parameter database); Phase 2 — S4 GaAs/AlGaAs QW (validates confinement); Phase 3 — S5/S6 advanced QW (broken-gap, strain); Phase 4 — S7 wire (most complex, requires wire-specific configs).
- **Three tolerance tiers:** Exact (machine precision, for analytical quantities), analytical (sub-percent, for Kane/Roth formulas), numerical (percent-level, for published computational/experimental comparisons).

---

## Dependencies / Assumptions

- The `gfactorCalculation` and `opticalProperties` executables must support the material systems listed. Where a config file doesn't exist in `tests/regression/configs/`, one must be created. Known missing configs: `gfactor_bulk_inas_cb.cfg`, `gfactor_bulk_insb_cb.cfg`, `qw_inas_gaas_strained.cfg`, `wire_inas_gfactor.cfg`.
- S6 was changed from "InGaAs/GaAs" to "InAs/GaAs" because `parameters.f90` has no InGaAs ternary alloy. InAs/GaAs strained QW uses materials that exist in the database and provides natural lattice-mismatch strain.
- S5 topological observables (band inversion, Z2 invariant) require the `topologicalAnalysis` executable to support real 8-band material QW configs in topology mode — a path currently untested (all existing topology configs use GaAs-based model Hamiltonians). This capability must be verified during planning; if unavailable, the "topo" column for S5 is deferred.
- Published reference values for each observable must be available in the cited literature. Systems where references are thin may have fewer observables.
- The shared utility layer must be created from scratch — no shared Python utilities exist in `tests/support/` (Fortran-only). The existing verification ladder Python scripts each contain inline parsing/comparison helpers (~15-25 lines each) that serve as the extraction source.
- S5 (InAs/GaSb) covers both broken-gap physics and topological observables since HgTe/CdTe parameters are unavailable in `parameters.f90`.

---

## Outstanding Questions

### Deferred to Planning

- [Affects R2] [Needs research] What specific published values are available for each observable per system? The exact literature citations and expected values need research during planning.
- [Affects R10] [Technical] The shared utility layer must be new (tests/support/ is Fortran-only). Planning should extract common patterns from existing rung scripts (~15-25 lines each of parsing/comparison logic) into a new Python module under tests/integration/.
- [Affects R6] [Technical] What tolerance is appropriate for strained QW observables (InAs/GaAs) — depends on what published reference is available.
- [Affects R13] [Technical] How should the benchmark table handle systems where some observables are not yet testable (e.g., wire optics lacking published reference)?
- [Affects R4-R7] [Technical] Planning must assign concrete tolerance values per tier and map each observable to its tier. Each tolerance should document what source of error it accounts for: parameter reproduction, model limitation, or numerical discretization.
