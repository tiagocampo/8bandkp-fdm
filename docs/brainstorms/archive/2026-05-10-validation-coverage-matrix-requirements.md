---
date: 2026-05-10
topic: validation-coverage-matrix
---

# Validation Coverage Matrix

## Summary

A dual-layer coverage tracking system for the k.p test suite. A YAML universe file declares which (observable, geometry, material) cells should be covered; tests self-report what they cover via structured annotations; a Python generator cross-checks both directions and produces a physics heat map plus an infrastructure traceability table. CI gates on required-tier cells; aspirational gaps are reported as info. The generated output doubles as a publication-ready coverage summary.

---

## Problem Frame

The codebase has ~80 regression configs, 4 verification ladder rungs, 7 standard-star scripts, and ~15 integration tests — but no systematic view of what physics is validated where. Test creation is ad hoc: when a gap is discovered, it's discovered incidentally, not by consulting a coverage map. The archived benchmark matrix tracks lecture-chapter status, not physics-observable coverage across geometries and materials.

A methods paper needs a reproducible statement of what's validated and how. A reviewer running the test suite should see not just "tests pass" but "GaAs bulk band gap validated against Vurgaftman 2001; InAs wire g-factor has regression but no published benchmark." Today that statement doesn't exist in automated form.

---

## Requirements

**Universe definition**

- R1. A YAML file declares the full set of (observable, geometry, material) cells that should be covered, with each cell annotated as `required` or `aspirational`.
- R2. Each cell in the universe file carries: observable name, geometry (bulk / QW / wire), material name (individual material, e.g. GaAs, InAs, InSb), tier, and optional published reference.
- R3. The universe file is the authoritative definition of "what should be tested." Adding or removing a cell is an explicit scope decision.

**Test self-reporting**

- R4. Tests self-report coverage via structured comment lines in their source files, specifying observable, geometry, and material per cell.
- R5. A single test file can report multiple cells.
- R6. The self-report format is machine-parseable by the generator (exact syntax determined during planning).

**Generator**

- R7. A Python generator reads the universe file and scans test source files for self-report annotations, then cross-checks both directions.
- R8. The generator detects three structural conditions: (a) required cells with no test coverage (CI failure), (b) aspirational cells with no test coverage (info), (c) test annotations claiming cells not in the universe (orphan warning).
- R9. When a test annotation claims a cell that exists in the universe, the generator resolves the cell's status to green or yellow based on whether a published reference is cited. When no test covers a cell, the status is red.

**CI integration**

- R10. The generator runs as a ctest target. CI exits non-zero when any required-tier cell has red status.
- R11. Aspirational-tier red cells, orphan annotations, and the full coverage summary are printed as info and do not cause CI failure.
- R12. When a required cell transitions from green/yellow to red (e.g., a covering test was removed), CI fails — this is a coverage regression.

**Output**

- R13. The generator produces a physics-layer heat map: a table with axes (observable x geometry x material), cells colored green/yellow/red.
- R14. The generator produces an infrastructure-layer traceability table: each covered cell maps to test label, executable, config file, and published reference.
- R15. The combined output is formatted as markdown suitable for inclusion in a methods paper or documentation.

**Retrofitting existing tests**

- R16. All existing verification ladder rungs, standard-star scripts, and integration tests that validate physics observables are annotated with self-report comments during implementation.

---

## Acceptance Examples

- AE1. **Covers R7, R8, R10.** The generator runs on the current test suite, finds that Eg/bulk/GaAs is covered by `verify_8band_rung1_bulk_k0.py` (required, green), that Chern/wire/GaAs has no test (aspirational, red), and that a test claims coverage for "strain/wire/GaAs" which is not in the universe file (orphan warning). CI exits zero because all required cells are green.

- AE2. **Covers R10, R12.** A test covering g-factor/bulk/GaAs (required, green) is removed from the suite. The next CI run detects the cell is now red and fails, printing "COVERAGE REGRESSION: required cell g-factor/bulk/GaAs lost coverage."

- AE3. **Covers R13, R14, R15.** Running the generator standalone produces a markdown table showing all cells with their status, plus a traceability section listing test names, configs, and references per cell. The output is directly pasteable into a methods paper.

---

## Success Criteria

- Running the generator on the current test suite produces an accurate heat map with no false positives (covered cells marked uncovered) or false negatives (uncovered cells marked covered).
- CI fails if any required-tier cell loses coverage.
- A reviewer can read the generated output and understand exactly what's validated, against what reference, and where the gaps are — without reading the test source code.
- Adding a new test requires only: (a) writing the test with a self-report annotation, (b) optionally adding the cell to the universe file if it's new.

---

## Scope Boundaries

- Not replacing the archived lecture-chapter benchmark matrix (orthogonal tracking dimension)
- Not auto-discovering test coverage from config filenames or output parsing (fragile)
- Not generating coverage-over-time drift reports or historical tracking
- Not including Richardson extrapolation convergence testing (ideation idea #6)
- Not including performance benchmarking or timing measurements
- Not auditing Vurgaftman parameter provenance (rejected ideation idea #4)
- Not generating new tests or configs (the matrix surfaces gaps; filling them is separate work)

---

## Key Decisions

- **Dual-layer output (physics + infrastructure):** The physics heat map serves gap analysis and publication; the infrastructure traceability table serves CI and developer debugging. A single-layer approach would force the reader to infer one from the other.
- **Hybrid data model (universe file + self-report):** Pure universe-file maintenance creates update friction; pure auto-discovery is unreliable. The hybrid gives explicit scope control (universe file) with low per-test maintenance (self-report annotations).
- **Required/aspirational tiers:** Starting with all cells as required would make CI impossible to pass (many cells are uncovered). Tiered gating lets the matrix ship immediately with current coverage as required and future targets as aspirational.
- **Individual materials, not families:** Standard-stars and verification ladder already use individual materials. Family-level grouping (III-As, III-Sb) would lose granularity needed for gap analysis.
- **Green = published benchmark, yellow = regression-only, red = untested:** Three-status model matches the ideation doc's traffic-light scheme. The CI gate (R10) tests existence only — a yellow cell (regression-only, no published reference) still counts as "covered" for CI purposes. The green/yellow distinction serves the publication artifact, not the CI gate.
- **Unit tests excluded from coverage universe:** pFUnit unit tests (.pf files in tests/unit/) validate code-level contracts (matrix symmetry, solver convergence) rather than physics observables against reference values. They are out of scope for the coverage matrix.

---

## Dependencies / Assumptions

- The self-report annotation format must be chosen during planning and then applied retroactively to all existing verification ladder rungs, standard-star scripts, and relevant integration tests.
- The universe file's initial contents require an upfront audit of all existing test configs and their coverage — a one-time cost.
- The generator depends on being able to statically parse test source files for annotations (no runtime introspection needed).
- The standard-star aggregation script (`tests/integration/aggregate_star_benchmarks.py`) may serve as a starting point for the publication-output formatting.
- PyYAML is a new Python dependency required by the generator for parsing the universe file. It is not currently used elsewhere in the project (existing scripts depend only on numpy and the standard library).
- Existing test labels (unit, regression, verification, standard-star) are used as-is in the infrastructure layer — no label restructuring required.

---

## Outstanding Questions

### From 2026-05-10 review

- [Affects R8, R9] R8 defines three structural conditions (red, info, orphan) but R9 adds green/yellow status resolution without a unified status model. Consider consolidating into a single requirement or adding a status enumeration that explicitly defines all four statuses (green, yellow, red, orphan) and their CI semantics. (coherence, P1)
- [Affects R10, R12] R12's "coverage regression" framing implies temporal transition detection (comparing against a previous run's status), but no artifact stores previous status. Either clarify R12 as a restatement of R10's state check, or specify a baseline storage mechanism. (coherence, P1)
- [Affects Problem Frame] Whether a simpler static-table approach (extending the existing benchmark-matrix.md or aggregate_star_benchmarks.py) suffices for the methods paper before investing in the full annotation + generator infrastructure. The three stated purposes (gap analysis, CI gate, publication) may not all require the same solution. (product-lens, P1)
- [Affects R4-R5] No mechanism verifies that annotations match actual test behavior. The system is entirely trust-based — consider spot-checking annotation accuracy (e.g., verifying config file existence, cross-referencing material names). (product-lens, P1)
- [Affects Problem Frame] Whether filling the top 5-10 highest-impact coverage gaps first delivers more value for the methods paper than building the tracking infrastructure. (product-lens, P1)
- [Affects R4-R7, R10-R12] Test annotation drift is undetectable: CI gate checks annotation presence, not test pass/fail. Consider consuming ctest JUnit/XML output to correlate annotations with pass/fail status. (adversarial, P1)
- [Affects R1-R3] Universe cells are coarse-grained (observable x geometry x material) but tests validate at specific k-points. "Green" means validated at one scenario, not comprehensively. Consider adding a regime/conditions dimension. (adversarial, P1)
- [Affects R1-R3, R10-R12] Universe file has no staleness detection relative to codebase capabilities. Consider cross-referencing against parameters.f90 materials list or executable capabilities. (adversarial, P1)
- [Affects R4-R5, R12] No failure mode for semantic annotation drift (test modified, annotation unchanged). Consider self-validating annotations or placing annotations in assertion code. (adversarial, P1)
- [Affects Requirements] Maintenance surface is substantial for a solo-developer codebase with no staleness detection. Consider whether a simpler first version (ctest labels + filenames, no annotations) provides 80% of the value. (product-lens, P2)

### Deferred to Planning

- [Affects R4, R6] [Technical] What exact syntax should the self-report annotation use? Inline comments (`# COVERAGE: observable=Eg geometry=bulk material=GaAs`) vs. a structured header block vs. a sidecar file. Planning should evaluate parseability vs. maintenance burden. Should shell scripts (.sh) and Python scripts (.py) use identical annotation syntax, or language-appropriate variants? Both use `#` for comments, but bash has no multi-line string literals. The parser must handle both file types.
- [Affects R1, R2] [Technical] What is the definitive list of observables for the initial universe file? A starting set would be: Eg, Delta_SO, m*_e, m*_hh, m*_lh, g*_cb, g*_vb, subband spacing, Stark shift, optical absorption edge, ISBT, strain shift, HH-LH splitting, Chern number, Z2 invariant, LDOS. Planning should enumerate from parameters.f90 and the four executables' capabilities.
- [Affects R7] [Technical] How should the generator discover test source files? Directory scan of `tests/integration/` with filename pattern? Explicit file list? CMake-level registration?
- [Affects R13, R14] [Technical] For the publication output, should the generator produce a single combined document or separate files (heat map + traceability)?
- [Affects R16] [Effort estimate] How many existing test files need retroactive annotation? Actual count: ~43 files (20 Python verify_*.py scripts + 23 inline shell scripts with physics PASS/FAIL assertions in tests/integration/). Note: 8 shell scripts delegate to Python verify scripts — decide whether to annotate the shell script, the Python script, or both.
