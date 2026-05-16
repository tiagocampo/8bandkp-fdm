---
title: Fix Standard-Star Benchmark Code Review Findings
type: fix
status: active
date: 2026-05-09
origin: Code review of docs/plans/2026-05-09-002-feat-standard-star-benchmarks-plan.md implementation
---

# Fix Standard-Star Benchmark Code Review Findings

## Summary

Fix 24 findings from an 8-reviewer code review of the standard-star benchmark suite (S1-S7). One P0 (wrong ctest label), one P1 (S5 overlap check uses hardcoded Python constants instead of Fortran output), and 22 P2-P3 issues spanning dead code, duplicated helpers, inconsistent patterns, loose assertions, and missing deliverables. Organized into 4 dependency-ordered units: refactor shared helpers first, normalize all 7 scripts to use them, fix physics assertions, then fix infrastructure.

---

## Problem Frame

The standard-star benchmark suite (commit `1eed461`) passed basic functionality but a thorough 8-reviewer code review found 24 issues. The most critical: `ctest -L standard-star` finds zero tests because the label is `"verification"` instead. S5's band overlap test validates hardcoded Python constants against a literature value (two independent published sources, so not tautological — but the Fortran code output is never checked). S6's Bir-Pikus diagnostic omits the Poisson correction (~4% error in the printed diagnostic value, not ~2x as initially reported). The remaining findings are code quality, consistency, and assertion-tightening issues.

---

## Requirements

- R1. CTest label must be `"standard-star"` so `ctest -L standard-star` discovers all 7 tests (plan KD8, R9)
- R2. S5 band overlap must validate Fortran code output (eigenvalue-based), not just hardcoded Python constants vs literature (plan KD4)
- R3. S6 Bir-Pikus HH-LH diagnostic must include Poisson correction (~4% correction for InAs, not ~2x as initially claimed)
- R4. Remove dead code from `star_helpers.py` (`run_star_test`, `cleanup_tmpdir`)
- R5. Centralize duplicated helpers into `star_helpers.py` (run wrappers, Roth formula, HBAR2_OVER_2M0)
- R6. All 7 scripts must use consistent patterns (tempdir management, row format, import strategy)
- R7. S1 must report FAIL when a subprocess crashes (not silently return empty rows)
- R8. S4 absorption onset, S6 CB spacing, and S7 wire g-factor assertions must be tightened with regression references
- R9. Add ctest TIMEOUT properties for expensive simulations
- R10. Create aggregation script placeholder (plan U9, R13)
- R11. Add numpy 2.0 compatibility shim for `np.trapezoid`

---

## Scope Boundaries

- This plan fixes review findings only — no new benchmark systems or physics capabilities
- S5 Z2 topological invariant test remains SKIP (plan KD5 capability gate) — deferred to follow-up
- S5 band overlap literature value sourcing (Liu et al. PRL 2008) is already established — this plan adds eigenvalue-based validation alongside it
- Effective mass extraction method unification (two-point vs adaptive fit) is advisory — noted but not changed
- Adding type hints to all public functions is advisory — noted but not changed in this fix pass
- No Fortran code changes — all fixes are in Python test scripts and CMake configuration

### Deferred to Follow-Up Work

- Z2 topology benchmark for InAsW/GaSbW broken-gap QW (requires new topology config)
- Machine-readable `--json` output format for benchmark scripts
- Unit tests for `star_helpers.py` utility functions
- Provenance tracking for regression reference values (separate data files with commit hashes)

---

## Context & Research

### Relevant Code and Patterns

- `tests/integration/star_helpers.py` — shared utility module (currently ~330 lines with dead code)
- `tests/integration/verify_star_*.py` — 7 benchmark scripts, each 100-200 lines
- `tests/CMakeLists.txt:741-811` — ctest registration for standard-star benchmarks
- `tests/regression/configs/` — canonical location for test configs (CLAUDE.md)
- Existing verification scripts in `tests/integration/` (verify_8band_rung*.py) — pattern reference for tempdir and error handling
- `docs/solutions/test-failures/csr-test-infra-stop1-oob-tautology-fixes-2026-05-09.md` — tautological assertion antipattern

### Institutional Learnings

- CSR test infra learning: tautological assertions (comparing code output against values derived from same parameters) can never fail. S5 band overlap compares Winkler 2003 parameters against Liu et al. 2008 (two independent published sources) — not tautological, but still doesn't validate Fortran code output.
- Verification ladder learning: CB effective mass must be compared against Kane 2-band formula, not Vurgaftman. The Roth g-factor comparison is analogous — compare against analytical formula, not experiment.
- Topological index learning: parameter sweeps must compute inside the loop, not before it. Bir-Pikus sign convention is a known source of subtle errors (CLAUDE.md boundary).

---

## Key Technical Decisions

- **KD1. Shared run_exe wrappers in star_helpers.** Four scripts define near-identical `run_bandstructure`/`run_gfactor` wrappers. Consolidate into `run_exe(build_dir, name, config, work_dir, timeout=300)` in `star_helpers.py` with executable existence check centralized. Avoids ~60 lines of duplication.

- **KD2. Roth g-factor as shared function.** Six scripts compute `g = 2 - 2*EP*DeltaSO/(3*Eg*(Eg+DeltaSO))` inline. Add `roth_gfactor(ep, eg, delta_so)` to `star_helpers.py` with Winkler 2003 Eq. 6.42 reference in docstring. Single source of truth for the formula.

- **KD3. S5 band overlap: add eigenvalue-based check alongside existing parameter check.** The current check compares Winkler 2003 parameters (`INASW_EC - GASBW_EV = -142 meV`) against Liu et al. PRL 2008 (`-150 meV`) — two independent published sources, not tautological. However, it validates parameter transcription, not Fortran code output. Add a second check that parses CB ground state and VB top from `bandStructure` eigenvalue output and compares the code-computed overlap against the literature value. Keep the existing parameter-vs-literature check as a separate self-consistency row. The bandStructure run (currently diagnostic-only in S5) is repurposed to provide eigenvalue data for both checks.

- **KD4. S6 Bir-Pikus: include Poisson correction in diagnostic.** Replace simplified `2*b*eps_biaxial` (0.241 eV) with full formula `2*b*(eps_zz - eps_xx)/2` where `eps_zz = -2*(C12/C11)*eps_xx` (gives 0.251 eV, a ~4% correction). The actual test (regression reference) is unaffected, but the diagnostic now prints the correct value for anyone debugging Bir-Pikus sign issues.

- **KD5. Assertion tightening uses regression references.** S4 absorption onset, S6 CB spacing, and S7 wire g-factor currently use wide range checks. Replace with `compare_value` against known 8-band regression references with appropriate tolerances (2-5%). This requires running the executables once to establish references, then hardcoding them (same pattern as S5/S6 g-factor).

- **KD6. Skip row format normalization.** S2/S7 store formatted strings while others store dicts. Normalizing would be cosmetic churn across all 7 scripts with no behavioral change. Leave as-is; future aggregation can handle both formats.

- **KD7. numpy compat shim for np.trapezoid.** Use `trapz_fn = getattr(np, 'trapezoid', np.trapz)` in `star_helpers.py` since `np.trapezoid` was introduced in numpy 2.0 and the project doesn't specify a minimum version.

---

## Open Questions

### Resolved During Planning

- How to fix S5 overlap: keep parameter-vs-literature check (two independent published sources), add eigenvalue-based check from Fortran output (KD3)
- Bir-Pikus correction scope: diagnostic output only (KD4) — test PASS/FAIL unaffected. Corrected value is 0.251 eV (~4% correction), not 0.502 eV as initially claimed
- S4 config path: copy to `tests/regression/configs/` rather than referencing `docs/benchmarks/`

### Deferred to Implementation

- Exact regression reference values for S4 absorption onset, S6 CB spacing, and S7 wire g-factor — must be established by running executables during implementation
- S5 eigenvalue parsing: the config has `numvb=8`, so CB starts at index 8 (not 4 as in the current diagnostic). The parser must use the correct band count to separate VB and CB eigenvalues

---

## Implementation Units

- U1. **Refactor star_helpers.py shared module**

**Goal:** Remove dead code, add shared wrappers and functions that eliminate duplication across the 7 scripts, copy stray config file.

**Requirements:** R4, R5, R11

**Dependencies:** None

**Files:**
- Modify: `tests/integration/star_helpers.py`
- Create: `tests/regression/configs/qw_gaas_algaas.cfg` (copy from `docs/benchmarks/`)

**Approach:**
- Remove `run_star_test` (lines 299-324) and `cleanup_tmpdir` (lines 327-331) — unused dead code
- Add `run_exe(build_dir, name, config_path, work_dir, timeout=300)` — resolves exe path from `build_dir/src/<name>`, checks `os.path.isfile`, calls `run_executable`. Replaces 4 copies of `run_bandstructure`/`run_gfactor`/`run_optical` across S3-S6
- Add `roth_gfactor(ep, eg, delta_so)` returning `2 - 2*ep*delta_so/(3*eg*(eg+delta_so))` with Winkler 2003 Eq. 6.42 reference. Replaces 6 inline copies
- Add `trapz_fn = getattr(np, 'trapezoid', np.trapz)` compatibility shim; update S4 to use it
- Fix numpy import strategy: replace `np = None` soft-fail with raising `ImportError` immediately with a clear message: `"numpy is required for standard-star benchmarks. Install with: pip install numpy"`
- Ensure `HBAR2_OVER_2M0` is documented as the single source of truth (already defined at line 155)
- Document the `compare_value` absolute/relative threshold behavior in docstring

**Config copy:**
- Copy `docs/benchmarks/qw_gaas_algaas.cfg` to `tests/regression/configs/qw_gaas_algaas.cfg` (S4 references this; U2 will update the path)
- Also update `tests/integration/verify_8band_rung3_qw.py` if it references the same `docs/benchmarks/` path

**Patterns to follow:**
- Existing `run_executable` function signature and error handling in `star_helpers.py`
- Existing `TOL_EXACT`, `TOL_ANALYTICAL`, `TOL_NUMERICAL` constants as public exports

**Test scenarios:**
- Happy path: `run_exe` correctly resolves and calls executables
- Error path: `run_exe` returns meaningful error when executable not found
- Happy path: `roth_gfactor` returns correct values for known materials (GaAs: -0.317, InAs: -14.6)
- Edge case: `roth_gfactor` with zero Eg raises or returns inf (division by zero guard)

**Verification:**
- All 7 scripts can import from updated `star_helpers.py` without errors
- Dead functions removed, new functions accessible

---

- U2. **Normalize all 7 benchmark scripts**

**Goal:** Make all scripts use consistent patterns — shared helpers, tempdir management, row format, error handling, and remove dead code paths.

**Requirements:** R6, R7

**Dependencies:** U1

**Files:**
- Modify: `tests/integration/verify_star_gaas_bulk.py` (S1)
- Modify: `tests/integration/verify_star_inas_bulk.py` (S2)
- Modify: `tests/integration/verify_star_insb_bulk.py` (S3)
- Modify: `tests/integration/verify_star_gaas_algaas_qw.py` (S4)
- Modify: `tests/integration/verify_star_inas_gasb_qw.py` (S5)
- Modify: `tests/integration/verify_star_inas_gaas_qw.py` (S6)
- Modify: `tests/integration/verify_star_inas_wire.py` (S7)

**Approach:**
Per-script changes (apply systematically):

**All scripts:**
- Replace local `HBAR2_OVER_2M0` with import from `star_helpers`
- Replace inline Roth formula with `star_helpers.roth_gfactor(ep, eg, delta_so)` call
- Replace local `run_bandstructure`/`run_gfactor`/`run_optical` wrappers with `star_helpers.run_exe` calls
- Remove redundant `import numpy` guards (centralized in star_helpers now)
- Normalize tempdir to `with tempfile.TemporaryDirectory(...) as tmpdir:` pattern (S1/S4 already use this)

**S1 (`verify_star_gaas_bulk.py`):**
- Fix silent failure swallowing: when `run_exe` returns non-zero or parsing fails, append a FAIL row (status `'FAIL'`) instead of returning empty list. Current behavior: empty rows → zero failures counted → PASS reported
- Add note in Eg/DeltaSO benchmark row output: `'self-consistency: eigenvalue == parameter at k=0'`

**S2 (`verify_star_inas_bulk.py`):**
- No changes beyond shared helper adoption (row format left as-is per KD6)

**S3 (`verify_star_insb_bulk.py`):**
- Remove local `TOL_EXACT = 1e-12` on line 65, import from `star_helpers` instead
- Remove local `extract_effective_mass` wrapper, call `star_helpers.extract_effective_mass` directly

**S4 (`verify_star_gaas_algaas_qw.py`):**
- Change `CONFIG_SUBBAND` to reference `tests/regression/configs/qw_gaas_algaas.cfg` (config copied in U1)
- Use `star_helpers` numpy compat shim for `np.trapezoid`
- Fix TE>TM polarization row: use actual ratio as delta, not binary 0.0/1.0

**S5 (`verify_star_inas_gasb_qw.py`):**
- Replace `check_z2` 65-line dead stub with a single TODO comment and one-line SKIP row emission
- Keep the `bandStructure` run (lines 275-296) — U3 will repurpose its output for eigenvalue-based overlap validation
- Fix CB eigenvalue index: config has `numvb=8`, so CB starts at index 8, not 4 (currently hardcoded as `cb_start=4` at line 288)
- Fix Z2 SKIP row: pass `'N/A'` string to format_benchmark_row instead of 0.0 floats, or skip the row entirely
- Update script header and inline comments to accurately describe 3-layer AlSbW/GaSbW/InAsW structure

**S6 (`verify_star_inas_gaas_qw.py`):**
- Import `HBAR2_OVER_2M0` from `star_helpers`
- Define `CB_INDEX = 6` as named constant (currently hardcoded magic number)
- Remove duplicate `bandStructure` run — parse both VB and CB eigenvalues from single run

**S7 (`verify_star_inas_wire.py`):**
- Define `CB_INDEX` constant for wire mode (row format left as-is per KD6)

**Test scenarios:**
- Happy path: each script runs end-to-end with the refactored helpers
- Error path: S1 reports FAIL when subprocess crashes (not silent PASS)
- Cleanup: no tempdir leaks on exception

**Verification:**
- All 7 scripts import successfully from updated `star_helpers.py`
- No local redefinitions of `HBAR2_OVER_2M0`, `TOL_EXACT`, Roth formula, or run helpers
- S5 `check_z2` is ≤ 5 lines (TODO + SKIP row)
- S5 bandStructure run is preserved (not removed) for U3 eigenvalue extraction

---

- U3. **Fix physics assertions and benchmark validity**

**Goal:** Fix tautological checks, wrong formulas, and overly loose assertions so benchmarks actually validate physics.

**Requirements:** R2, R3, R8

**Dependencies:** U2

**Files:**
- Modify: `tests/integration/verify_star_inas_gasb_qw.py` (S5 — overlap)
- Modify: `tests/integration/verify_star_inas_gaas_qw.py` (S6 — Bir-Pikus, CB spacing)
- Modify: `tests/integration/verify_star_gaas_algaas_qw.py` (S4 — absorption onset)
- Modify: `tests/integration/verify_star_inas_wire.py` (S7 — wire g-factor)

**Approach:**

**S5 band overlap (`verify_star_inas_gasb_qw.py`):**
- Keep the existing parameter-vs-literature check (`INASW_EC - GASBW_EV` vs Liu et al. 2008) as a "parameter self-consistency" benchmark row
- Add a new "band overlap (eigenvalue)" benchmark row that parses CB ground state and VB top from the `bandStructure` eigenvalue output (already running in S5)
- Eigenvalue extraction: with `numvb=8`, CB starts at index 8. Parse `eigenvalues[numvb]` (CB ground) and `eigenvalues[numvb-1]` (VB top), compute overlap as `(CB_ground - VB_top) * 1000` in meV
- Compare eigenvalue-based overlap against `OVERLAP_LITERATURE = -150.0` (Liu et al. PRL 2008) with `TOL_OVERLAP = 0.10`
- Note: QW eigenvalues include confinement shifts, so the eigenvalue-based overlap will differ from the bulk parameter overlap. The literature value is for bulk band alignment, so use a generous tolerance (10%) and note in the script that this validates the sign and order of magnitude

**S6 Bir-Pikus diagnostic (`verify_star_inas_gaas_qw.py`):**
- Replace `HH_LH_SPLITTING_BULK = 2.0 * INAS_B * EPS_BIAXIAL` (gives 0.241 eV) with full biaxial formula:
  ```
  C11 = 832.9  # GPa
  C12 = 452.6  # GPa
  eps_zz = -2.0 * (C12/C11) * eps_biaxial
  HH_LH_SPLITTING_BULK = abs(INAS_B * (eps_zz - eps_biaxial))
  ```
  This gives ~0.251 eV (a ~4% correction from 0.241 eV, not ~2x as initially claimed)
- The regression test (HHLH_REF = 0.082710) remains unchanged — only the printed diagnostic value is corrected
- Add comment: "Full biaxial Bir-Pikus including Poisson correction (Chuang 2003, Ch. 4)"

**S4 absorption onset (`verify_star_gaas_algaas_qw.py`):**
- Replace 458 meV range check `[onset_min, onset_max]` with `compare_value` against a regression reference onset energy
- Regression reference TBD: run `opticalProperties` with the absorption config and extract onset from output, then hardcode as `ONSET_REF` with `TOL_ONSET = 0.02`
- Keep the range check as a secondary sanity bound

**S6 CB spacing (`verify_star_inas_gaas_qw.py`):**
- Replace `[0.001, 0.600]` bounds with `compare_value` against regression reference
- Known 8-band result: ~409 meV for the 20 nm InAs/GaAs strained well
- Hardcode as `CB_SPACING_REF` with `TOL_CB_SPACING = 0.05`
- Remove the fake `compare_value` row with midpoint expected value and binary delta

**S7 wire g-factor (`verify_star_inas_wire.py`):**
- Replace `tolerance_frac = 0.50` (factor-of-4 range) with regression reference
- Run `gfactorCalculation` with wire config, extract gz, hardcode as `G_WIRE_REF` with `TOL_G_WIRE = 0.10`
- Keep bulk Roth value as a sanity bound: `abs(g_wire)` should be within factor 2 of `abs(g_roth)`

**Execution note:** Regression reference values must be established by running executables first. Run each relevant config, capture output, extract values, then hardcode as references. This is execution-time discovery.

**Test scenarios:**
- Happy path: S5 overlap comparison uses parsed eigenvalues, not hardcoded constants
- Happy path: S6 Bir-Pikus diagnostic prints ~0.251 eV (not 0.241 eV)
- Edge case: S5 eigenvalue parsing handles missing or empty output file
- Regression: S4 onset within 2% of reference
- Regression: S6 CB spacing within 5% of reference
- Regression: S7 wire g-factor within 10% of reference

**Verification:**
- S5 has two overlap checks: parameter-vs-literature (existing) and eigenvalue-vs-literature (new)
- S6 diagnostic prints correct Bir-Pikus value with Poisson correction
- S4, S6, S7 use `compare_value` with regression references instead of wide range checks

---

- U4. **Infrastructure fixes and missing deliverables**

**Goal:** Fix ctest configuration, add timeouts, create aggregation script.

**Requirements:** R1, R9, R10

**Dependencies:** U2

**Files:**
- Modify: `tests/CMakeLists.txt`
- Create: `tests/integration/aggregate_star_benchmarks.py`

**Approach:**

**CMakeLists.txt fixes:**
- Change line 810 from `PROPERTIES LABELS "verification"` to `PROPERTIES LABELS "standard-star;verification"` (dual labels: `ctest -L standard-star` runs only benchmarks, `ctest -L verification` runs all verification including benchmarks)
- Add `TIMEOUT 600` to `set_tests_properties` for all 7 standard-star tests (wire FEAST and QW optics can take several minutes on slow CI)

**Aggregation script:**
- Create `tests/integration/aggregate_star_benchmarks.py` per plan U9/R13
- Script runs each of the 7 standard-star scripts via subprocess and parses their stdout benchmark table output
- Concatenates parsed rows into a single markdown table with header
- Writes to stdout or a configurable output path
- No new CLI flags needed on the benchmark scripts — parse their existing output format

**Test scenarios:**
- Happy path: `ctest --test-dir build -L standard-star` discovers and runs all 7 tests
- Happy path: `python3 tests/integration/aggregate_star_benchmarks.py` produces a complete markdown table
- Error path: aggregation script handles individual test failures gracefully (includes FAIL rows in output)

**Verification:**
- `ctest -N -L standard-star` lists exactly 7 tests
- `ctest -L standard-star` runs only the standard-star tests
- `ctest -L verification` runs both standard-star and verification ladder tests
- Aggregation script exists and produces valid markdown

---

## Risks & Dependencies

| Risk | Mitigation |
|------|------------|
| Regression reference values (S4 onset, S6 spacing, S7 wire g) may change with any Fortran code update | Document provenance: record git commit hash and parameter set alongside each reference. Accept that regression references are snapshots, not ground truth |
| S5 eigenvalue-based overlap differs from bulk parameter overlap due to QW confinement shifts | Use generous 10% tolerance; note in script that this validates sign and order of magnitude. The existing parameter-vs-literature check remains as a separate row |
| S6 6.7% lattice mismatch exceeds Bir-Pikus validity — regression references may shift with solver changes | Acknowledged in script docstring as stress-test. Tolerance (5%) allows small solver variations |
| numpy compat shim may not cover all numpy 1.x differences | The only numpy 2.0-only function used is `np.trapezoid`; the shim is targeted and sufficient |

---

## Sources & References

- **Origin plan:** `docs/plans/2026-05-09-002-feat-standard-star-benchmarks-plan.md`
- **Review artifacts:** `/tmp/compound-engineering/ce-code-review/20260509-164959-e0255952/`
- **Learning:** `docs/solutions/test-failures/csr-test-infra-stop1-oob-tautology-fixes-2026-05-09.md` (tautological assertion pattern)
- **Reference:** Winkler 2003, Eq. 6.42 (Roth g-factor)
- **Reference:** Chuang 2003, Ch. 4 (Bir-Pikus strain with Poisson correction)
