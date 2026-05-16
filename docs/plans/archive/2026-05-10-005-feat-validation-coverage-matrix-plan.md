---
title: "feat: Validation Coverage Matrix"
type: feat
status: active
date: 2026-05-10
origin: docs/brainstorms/2026-05-10-validation-coverage-matrix-requirements.md
---

# feat: Validation Coverage Matrix

## Summary

A Python generator (`tests/integration/coverage_matrix.py`) reads a YAML universe file declaring (observable, geometry, material) cells, scans ~51 test source files (20 verify_*.py + 31 test_*.sh in `tests/integration/`) for `# COVERAGE:` annotations, and cross-checks both directions to produce a dual-layer markdown report (physics heat map + infrastructure traceability). A ctest target gates on required-tier red cells. Existing verification ladder rungs, standard-star scripts, and integration tests are retroactively annotated.

*(see origin: docs/brainstorms/2026-05-10-validation-coverage-matrix-requirements.md)*

---

## Problem Frame

The codebase has ~80 tests across 4 ctest labels but no systematic view of what physics is validated where. Test creation is ad hoc; gaps are discovered incidentally. A methods paper needs a reproducible statement of what's validated and how — today that statement doesn't exist in automated form.

*(see origin: Problem Frame)*

---

## Requirements

- R1. YAML universe file declares the full set of (observable, geometry, material) cells, each annotated as `required` or `aspirational`
- R2. Each cell carries: observable name, geometry (bulk/QW/wire), material name, tier, optional published reference
- R3. Universe file is the authoritative definition of what should be tested
- R4. Tests self-report coverage via structured `# COVERAGE:` comment lines specifying observable, geometry, and material
- R5. A single test file can report multiple cells
- R6. Self-report format is machine-parseable
- R7. Generator reads universe file and scans test source files for annotations, cross-checks both directions
- R8. Generator detects: (a) required cells with no coverage (red, CI failure), (b) aspirational cells with no coverage (red, info), (c) annotations claiming cells not in universe (orphan warning)
- R9. Status resolution: green = covered + published reference, yellow = covered + regression-only, red = not covered, orphan = annotation references unknown cell
- R10. Generator runs as a ctest target; CI exits non-zero when any required-tier cell is red
- R11. Aspirational red, orphan annotations, and full summary are printed as info
- R12. Removing a covering test causes CI failure on next run (state check, not temporal comparison)
- R13. Physics-layer heat map: table with axes (observable x geometry x material), cells colored green/yellow/red
- R14. Infrastructure-layer traceability table: each covered cell maps to test label, executable, config file, published reference
- R15. Combined output formatted as markdown suitable for methods paper
- R16. All existing verification ladder rungs, standard-star scripts, and integration tests annotated during implementation

**Origin acceptance examples:** AE1 (covers R7, R8, R10), AE2 (covers R10, R12), AE3 (covers R13, R14, R15)

---

## Scope Boundaries

- Not replacing the archived lecture-chapter benchmark matrix
- Not auto-discovering coverage from config filenames or output parsing
- Not generating coverage-over-time drift reports or historical tracking
- Not including Richardson extrapolation convergence testing
- Not including performance benchmarking or timing measurements
- Not auditing Vurgaftman parameter provenance
- Not generating new tests or configs (matrix surfaces gaps; filling them is separate work)
- Not creating a CI/CD pipeline (no `.github/workflows/` — ctest gating is the local equivalent)
- Not annotating lecture scripts (`scripts/lecture_*.py`) — they are pedagogical companions, not CI-gated tests
- Not including pFUnit unit tests in the coverage universe — they validate code-level contracts, not physics observables

### Deferred to Follow-Up Work

- Semantic annotation drift detection (verifying annotations match actual test behavior): defer to a future validation pass
- Historical coverage tracking / trend reporting: separate effort after the matrix stabilizes
- CI/CD pipeline creation: independent infrastructure work

---

## Context & Research

### Relevant Code and Patterns

- `tests/integration/aggregate_star_benchmarks.py` (120 lines) — existing aggregation pattern: hardcoded script list, subprocess invocation, markdown table assembly. Not a ctest target itself; model ctest registration after the verify_8band_rung* pattern instead.
- `tests/integration/star_helpers.py` — shared helpers including `compare_value()`, `format_benchmark_row()`, `run_exe()`. Used by verification and star scripts.
- `tests/CMakeLists.txt` (812 lines) — central test registration with 4-label taxonomy (unit, regression, verification, standard-star). All `add_test()` calls and `set_tests_properties()` labels live here. Standard-star tests carry compound label `"standard-star;verification"`.
- `tests/integration/verify_8band_rung1_bulk_k0.py` — canonical verify script: `MATERIALS` dict, `check_*()` functions, `[R<n>]` PASS/FAIL reporting. Template for annotation placement.
- `tests/integration/test_bulk_bandstructure.sh` — canonical shell test: `set -euo pipefail`, temp directory, config copy, executable invocation, output comparison. Template for `.sh` annotation placement.
- `docs/solutions/best-practices/8band-verification-ladder-2026-05-09.md` — 4-rung verification ladder design, maps naturally to universe cell hierarchy.

### Institutional Learnings

- **Subprocess timeout pattern** (from `docs/solutions/test-failures/lecture-test-pair-review-findings-2026-05-10.md`): all `subprocess.run` with `timeout` must catch `TimeoutExpired`; established pattern in `star_helpers.py`.
- **Aggregate reporting** (from `docs/solutions/logic-errors/standard-star-benchmark-suite-logic-errors-2026-05-09.md`): never nest pipe-delimited markdown tables; use `##` headings with raw rows.
- **Avoid tautological assertions** (from `docs/solutions/test-failures/csr-test-infra-stop1-oob-tautology-fixes-2026-05-09.md`): verify expected and actual derive from independent sources.

### External References

- PyYAML 6.0.3 already installed system-wide — no new dependency installation needed.

---

## Key Technical Decisions

- **Inline comment annotations** (`# COVERAGE: observable=Eg geometry=bulk material=GaAs`): Both `.sh` and `.py` files use `#` comments, so one syntax works across all test types. Avoids sidecar files and keeps annotations co-located with the test code they describe. *(resolves deferred Q1)*
- **Unified 4-status model** (green/yellow/red/orphan): Consolidates R8's three structural conditions with R9's green/yellow distinction into a single enumeration. CI semantics are unambiguous: red+required → fail; all other combinations → info. *(resolves P1 coherence question on R8/R9)*
- **R12 as state check, not temporal comparison**: Every CI run checks the current state of (universe file + annotations). No baseline storage or comparison against previous runs. If a required cell is red, CI fails — regardless of whether it was green yesterday. *(resolves P1 coherence question on R12)*
- **Generator at `tests/integration/coverage_matrix.py`**: Lives alongside `aggregate_star_benchmarks.py` which it extends. Keeps coverage infrastructure with the test infrastructure.
- **Universe file at `tests/integration/validation_universe.yml`**: Co-located with the test scripts it governs. YAML format for human-editability and PyYAML parsing.
- **ctest label `"coverage"` as meta-check**: Separate from the existing 4 physics-test labels. The coverage matrix is a structural check on the test suite, not a physics test itself.
- **Lecture scripts excluded from coverage universe**: They run outside ctest, serve a pedagogical purpose, and would conflate two different quality dimensions if included.
- **Observable enum derived from executable capabilities**: Initial set of ~20 observables enumerated from what the four executables actually compute, not from a theoretical taxonomy. Observable names match the language already used in verification scripts and the verification ladder doc.

---

## Open Questions

### Resolved During Planning

- **Annotation syntax**: Inline `# COVERAGE:` comments in both .sh and .py files. Multiple lines for multiple cells. *(resolves deferred Q1)*
- **Observable list**: Eg, Delta_SO, m*_e, m*_hh, m*_lh, g*_cb, g*_vb, subband_spacing, stark_shift, HH_LH_splitting, absorption_edge, gain_peak, ISBT, strain_shift, chern_number, z2_invariant, ldos, majorana_modes, landau_levels, fermi_level. Finalized during universe file population. *(resolves deferred Q2)*
- **File discovery**: glob scan of `tests/integration/verify_*.py` and `tests/integration/test_*.sh`. No verify scripts exist in `scripts/`. *(resolves deferred Q3)*
- **Output format**: Single markdown document with two sections (physics heat map + infrastructure traceability). *(resolves deferred Q4)*
- **Retrofitting scope**: ~51 files: 20 verify_*.py (including 4 verification ladder rungs and 7 standard-star scripts) + 31 test_*.sh. Of the 31 shell scripts, ~7 delegate to verify scripts (annotate the verify script, not the wrapper); ~24 have inline physics assertions (annotate the .sh file directly). Shell scripts delegating to verify scripts get annotations in the verify script, not the shell wrapper. *(resolves deferred Q5)*
- **Status model**: green/yellow/red/orphan with CI semantics defined per tier. *(resolves P1 coherence on R8/R9)*
- **R12 semantics**: State check (no temporal baseline). *(resolves P1 coherence on R12)*

### Deferred to Implementation

- **Exact universe cell population**: The initial audit of which tests cover which observables requires reading each test file and identifying its physics assertions. The universe file schema is defined here; the cell content is an execution-time enumeration.
- **Spot-check annotation validity**: Whether to validate that referenced config files exist or material names match `parameters.f90`. Deferred to implementation as a robustness enhancement.
- **Annotation drift detection**: No mechanism to verify annotations match actual test behavior. The system is trust-based by design (acknowledged P1 concern). Spot-checks are a best-effort mitigation.

---

## Implementation Units

- U1. **Universe File Schema and Initial Population**

**Goal:** Create `validation_universe.yml` with a well-defined schema and populate it with the initial set of (observable, geometry, material, tier) cells based on an audit of existing test coverage.

**Requirements:** R1, R2, R3

**Dependencies:** None

**Files:**
- Create: `tests/integration/validation_universe.yml`
- Reference: `tests/integration/verify_8band_rung1_bulk_k0.py` (audit for initial cells)
- Reference: `tests/integration/verify_star_gaas_bulk.py` (audit for initial cells)
- Reference: `docs/solutions/best-practices/8band-verification-ladder-2026-05-09.md` (observable hierarchy)
- Reference: `src/physics/parameters.f90` (material name enumeration)

**Approach:**
- Define YAML schema: list of cells, each with `observable`, `geometry`, `material`, `tier` (required/aspirational), optional `reference` (published citation string)
- Add a `metadata` block at the top with `observables` enum, `geometries` enum, `materials` enum for validation
- Audit existing verification ladder rungs (4), standard-star scripts (7), and integration verify scripts (~20) to determine which cells are already covered
- Mark covered cells as `required`; mark known-gap cells as `aspirational`
- Observable names drawn from the executable capabilities (band gap, effective mass, g-factor, etc.)

**Patterns to follow:**
- `tests/integration/aggregate_star_benchmarks.py` for the idea of a hardcoded universe (but YAML makes it externalized)
- Material names must match `parameters.f90` exactly (GaAs, InAs, InSb, AlAs, GaSb, AlSb, etc.)

**Test scenarios:**
- Happy path: YAML file parses without error and contains >= 30 cells
- Edge case: All cells have required fields (observable, geometry, material, tier); generator schema validation catches missing fields
- Edge case: Observable values all appear in the `observables` metadata enum
- Edge case: Material values match names in `parameters.f90`
- Edge case: Tier values are `required` or `aspirational` only

**Verification:**
- `python3 -c "import yaml; yaml.safe_load(open('tests/integration/validation_universe.yml'))"` succeeds
- At least 30 cells defined; all required cells correspond to existing test coverage

---

- U2. **Annotation Convention and Parser**

**Goal:** Define the `# COVERAGE:` annotation syntax and build a reusable parser that extracts annotations from `.sh` and `.py` test source files.

**Requirements:** R4, R5, R6

**Dependencies:** None

**Files:**
- Create: `tests/integration/coverage_matrix.py` (parser module within the generator)
- Reference: `tests/integration/verify_8band_rung1_bulk_k0.py` (annotation target example)
- Reference: `tests/integration/test_bulk_bandstructure.sh` (annotation target example)

**Approach:**
- Define annotation format: `# COVERAGE: observable=<name> geometry=<type> material=<mat> [ref=<citation>]`
- Multiple cells per file: multiple `# COVERAGE:` lines
- Optional `ref=` field for published reference (drives green vs yellow status)
- Parser: regex-based line scanner, reads each file, extracts all COVERAGE lines into structured dicts
- Returns list of `(file_path, observable, geometry, material, ref)` tuples
- Validates annotation fields against the observable/geometry/material enums from the universe file
- Reports malformed annotations (missing fields, unknown values) as warnings

**Technical design:**

The annotation line format is positional for required fields, with `ref=` as an optional trailing key-value pair:
```
# COVERAGE: observable=Eg geometry=bulk material=GaAs ref=Vurgaftman2001
# COVERAGE: observable=m*_e geometry=bulk material=InAs
```

Parser flow: for each file, read all lines, match against `^#\s*COVERAGE:\s*(.+)$`, parse key-value pairs from the capture group, validate fields, return structured records. Invalid annotations are collected as warnings, not errors (the generator reports them as orphan/malformed).

**Patterns to follow:**
- `aggregate_star_benchmarks.py` uses regex parsing of stdout for PASS/FAIL rows — similar line-scanning approach
- `star_helpers.py`'s `compare_value()` return pattern (structured dict) for test result data

**Test scenarios:**
- Happy path: parse single COVERAGE line from a .py file → correct (observable, geometry, material) tuple
- Happy path: parse single COVERAGE line from a .sh file → same result
- Happy path: parse multiple COVERAGE lines from one file → list of tuples
- Edge case: file with zero COVERAGE lines → empty list, no error
- Edge case: non-COVERAGE comment lines → ignored
- Error path: malformed annotation (missing `geometry=`) → collected as warning, not crash
- Error path: unknown observable name → collected as warning
- Integration: parse all existing verify_*.py files → list of annotations (initially empty since U5 hasn't run yet)

**Verification:**
- Parser module extracts annotations from a synthetic test file with known COVERAGE lines
- Parser handles both .sh and .py files identically
- Malformed annotations produce warnings, not exceptions

---

- U3. **Coverage Matrix Generator (Core Logic)**

**Goal:** Build the main generator that reads the universe file, scans test files for annotations, cross-references both, and produces the dual-layer markdown output.

**Requirements:** R7, R8, R9, R13, R14, R15

**Dependencies:** U1, U2

**Files:**
- Modify: `tests/integration/coverage_matrix.py` (add main logic and output formatting)
- Read: `tests/integration/validation_universe.yml`
- Read: `tests/integration/verify_*.py`, `tests/integration/test_*.sh`

**Approach:**
- Main entry: `coverage_matrix.py` takes `--source-dir` (repo root) and `--output` (markdown file path) arguments
- Reads universe file (U1) via PyYAML
- Scans test files for annotations using parser (U2)
- Cross-reference logic:
  1. For each universe cell, find matching annotations → green (if ref) or yellow (if no ref)
  2. For each universe cell with no matching annotation → red
  3. For each annotation with no matching universe cell → orphan
- Output section 1 — Physics heat map: markdown table with rows = observables, columns = geometry x material combinations. Cell contents: green/yellow/red indicator + reference citation for green cells
- Output section 2 — Infrastructure traceability: table mapping each covered cell to test file path, ctest label, config file (parsed from annotation or inferred from test), and reference
- CI semantics: exit code 1 if any required cell is red; exit code 0 otherwise
- Summary to stderr: counts per status, list of red required cells

**Technical design:**

Cross-reference algorithm:
1. Load universe cells into a dict keyed by `(observable, geometry, material)`
2. Load annotations into a list of `(file_path, observable, geometry, material, ref)` tuples
3. For each universe cell: search annotations for matching `(observable, geometry, material)`. If found → green (annotation has ref) or yellow (no ref). If not found → red.
4. For each annotation: if `(observable, geometry, material)` not in universe keys → orphan.
5. Report: heat map + traceability + summary.

The heat map table uses observable as the primary axis and (geometry, material) pairs as columns. For compact display, group by geometry first (bulk section, QW section, wire section) with material sub-columns.

**Patterns to follow:**
- `aggregate_star_benchmarks.py` for subprocess invocation pattern and markdown assembly
- `docs/solutions/logic-errors/standard-star-benchmark-suite-logic-errors-2026-05-09.md`: use `##` headings with raw rows, never nested tables
- ctest label extraction: read `tests/CMakeLists.txt` to map test names to labels, or accept label as an optional annotation field

**Test scenarios:**
- Happy path: universe has 3 cells, annotations cover all 3 → all green/yellow, exit 0
- Happy path: required cell covered with ref → green in heat map, appears in traceability
- Happy path: required cell covered without ref → yellow in heat map, appears in traceability
- Edge case: required cell with no annotation → red, exit 1, "UNCOVERED REQUIRED CELL" message to stderr
- Edge case: aspirational cell with no annotation → red, info only, exit 0
- Edge case: annotation references cell not in universe → orphan warning to stderr
- Integration: run against actual test suite (after U5 annotations) → produces accurate heat map with no false positives or false negatives
- Integration: verify AE1 scenario — Eg/bulk/GaAs covered, Chern/wire/GaAs aspirational red, orphan annotation for strain/wire/GaAs

**Verification:**
- `python3 tests/integration/coverage_matrix.py --source-dir .` runs without error
- Output markdown contains both heat map and traceability sections
- Exit code is 0 when all required cells are covered, 1 when any required cell is red

---

- U4. **ctest Registration**

**Goal:** Register the coverage matrix generator as a ctest target with a `coverage` label, so `ctest -L coverage` runs the coverage check.

**Requirements:** R10, R11, R12

**Dependencies:** U3

**Files:**
- Modify: `tests/CMakeLists.txt`

**Approach:**
- Add `add_test(NAME coverage_matrix COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/coverage_matrix.py --source-dir ${CMAKE_SOURCE_DIR} --output ${CMAKE_BINARY_DIR}/coverage_report.md WORKING_DIRECTORY ${CMAKE_BINARY_DIR})`
- Set label: `set_tests_properties(coverage_matrix PROPERTIES LABELS "coverage")`
- The test fails (nonzero exit) when required cells are red — ctest interprets nonzero as FAIL
- Add `TIMEOUT 120` property to prevent hanging on large test suites

**Patterns to follow:**
- Existing `add_test()` registration pattern for pure Python tests in `tests/CMakeLists.txt` (e.g., `verification_8band_rung1`)
- Label assignment via `set_tests_properties(... PROPERTIES LABELS "...")`

**Test scenarios:**
- Happy path: `ctest -R coverage_matrix` runs the generator and reports PASS when all required cells are covered
- Happy path: `ctest -L coverage` includes the coverage_matrix test
- Edge case: Generator finds required red cells → ctest reports FAIL
- Integration: `ctest --test-dir build` includes coverage_matrix alongside existing tests

**Verification:**
- `ctest --test-dir build -R coverage_matrix` succeeds (or fails with clear red-cell report)
- Coverage report appears at `build/coverage_report.md`

---

- U5. **Retrofit Existing Test Files with Annotations**

**Goal:** Add `# COVERAGE:` annotations to all existing verification ladder rungs, standard-star scripts, and integration verify scripts.

**Requirements:** R16

**Dependencies:** U1 (universe file for cell enumeration), U2 (annotation syntax definition)

**Files:**
- Modify: `tests/integration/verify_8band_rung1_bulk_k0.py`
- Modify: `tests/integration/verify_8band_rung2_dispersion.py`
- Modify: `tests/integration/verify_8band_rung3_qw.py`
- Modify: `tests/integration/verify_8band_rung4_wire.py`
- Modify: `tests/integration/verify_star_gaas_bulk.py`
- Modify: `tests/integration/verify_star_inas_bulk.py`
- Modify: `tests/integration/verify_star_insb_bulk.py`
- Modify: `tests/integration/verify_star_gaas_algaas_qw.py`
- Modify: `tests/integration/verify_star_inas_gasb_qw.py`
- Modify: `tests/integration/verify_star_inas_wire.py`
- Modify: `tests/integration/verify_star_inas_gaas_qw.py`
- Modify: `tests/integration/verify_bulk_benchmarks.py`
- Modify: `tests/integration/verify_qw_benchmarks.py`
- Modify: `tests/integration/verify_sc_benchmarks.py`
- Modify: `tests/integration/verify_qw_absorption_polarization.py`
- Modify: `tests/integration/verify_qw_state_character.py`
- Modify: `tests/integration/verify_landau_analytical.py`
- Modify: `tests/integration/verify_stark_shift.py`
- Modify: `tests/integration/verify_scattering_transition_sanity.py`
- Modify: `tests/integration/verify_wire_optical_selection.py`
- Modify: Additional verify/test scripts as discovered during audit (~8 more)
- Modify: Shell scripts with inline physics assertions and no paired verify script (~24 files)

**Approach:**
- For each verify script: read the script, identify what observables/geometries/materials it validates (from its MATERIALS dict and check functions), add `# COVERAGE:` lines after the docstring
- For shell scripts without a paired verify script: add `# COVERAGE:` lines in the header comment block
- Each annotation includes `ref=` when the test validates against a published value (Vurgaftman 2001, Winkler 2003, etc.)
- Annotations placed early in the file (after imports, before functions) for discoverability

**Execution note:** Each annotation requires reading the test's assertions to determine the correct observable/geometry/material triple — this is semantically non-trivial for tests that validate multiple observables across multiple materials. Budget review time proportional to the semantic difficulty, not the line count.

**Patterns to follow:**
- Annotation format: `# COVERAGE: observable=Eg geometry=bulk material=GaAs ref=Vurgaftman2001`
- Place after docstring, before imports or after imports, before function definitions
- Use the MATERIALS dict in verify scripts as the source of truth for which materials are tested

**Test scenarios:**
- Happy path: each annotated file has at least one valid COVERAGE line
- Edge case: no file has annotations referencing cells not in the universe (no orphans after retrofit)
- Integration: after retrofit, running the generator produces a heat map where all existing test coverage is correctly represented
- Integration: verify AE1 scenario passes end-to-end with the retrofitted annotations

**Verification:**
- `python3 tests/integration/coverage_matrix.py --source-dir .` reports no orphan warnings for the retrofitted files
- All previously-covered physics is represented as green or yellow in the heat map
- No previously-covered physics shows as red (no false negatives from missing annotations)

---

## System-Wide Impact

- **Interaction graph:** The coverage matrix generator is a read-only tool — it scans source files and the universe YAML, produces output, and exits. No callbacks, no middleware, no side effects on test execution.
- **Error propagation:** Generator exit code propagates to ctest. Nonzero → test FAIL. Stderr messages carry human-readable diagnostics (red required cells, orphan warnings).
- **State lifecycle risks:** None. The generator is stateless — no persistent state between runs. R12 semantics are state-check, not temporal comparison.
- **API surface parity:** Not applicable — this is a new standalone tool with no shared API surface.
- **Integration coverage:** The ctest registration (U4) is the sole integration point. It must not interfere with existing test execution timing or ordering.
- **Unchanged invariants:** Existing test files continue to function identically — annotations are comments that do not affect execution. The 4 existing ctest labels (unit, regression, verification, standard-star) are unchanged. Test discovery and execution are unaffected.

---

## Risks & Dependencies

| Risk | Mitigation |
|------|------------|
| Annotation drift — tests modified but annotations unchanged | Accept as known limitation; document as trust-based system. Future: spot-check config file existence in generator. |
| Universe file staleness — new code capabilities not reflected | Accept; universe file is an explicit scope decision (R3). Adding cells is a deliberate act. |
| Retrofitting effort underestimated — ~51 files, each requiring semantic analysis of test assertions | Prioritize by tier: required-tier coverage first, aspirational second. Run generator incrementally during retrofit to catch mismatches early. |
| PyYAML not available in all environments | Already confirmed available (6.0.3 system-wide). Add a `try: import yaml` guard with clear error message. |
| No CI/CD pipeline — coverage gating is local only | Accept as current-state limitation. ctest gating is the local equivalent; CI pipeline is deferred work. |
| Observable enumeration incomplete on first pass | Universe file is version-controlled and iteratively refined. Missing observables are caught by orphan warnings when tests add them. |

---

## Sources & References

- **Origin document:** [docs/brainstorms/2026-05-10-validation-coverage-matrix-requirements.md](docs/brainstorms/2026-05-10-validation-coverage-matrix-requirements.md)
- Existing aggregation: `tests/integration/aggregate_star_benchmarks.py`
- Test registration: `tests/CMakeLists.txt`
- Verification ladder design: `docs/solutions/best-practices/8band-verification-ladder-2026-05-09.md`
- Standard-star bug patterns: `docs/solutions/logic-errors/standard-star-benchmark-suite-logic-errors-2026-05-09.md`
