---
status: active
origin: docs/brainstorms/2026-05-09-standard-star-benchmarks-requirements.md
created: 2026-05-09
---

# Standard-Star Benchmark Systems

## Summary

Seven per-material Python verification scripts validating published observables (band gap, spin-orbit splitting, effective mass, g-factor, subband spacing, optical absorption, strain HH-LH splitting, topological Z2, wire g-factor) against literature references. Scripts produce ctest PASS/FAIL plus publication-ready markdown benchmark tables. Organized in 4 phases: bulk (S1-S3), standard QW (S4), advanced QW (S5-S6), wire (S7). Includes a shared `star_helpers.py` utility module, ~5 new executable configs, and ctest registration under a `standard-star` label. Failing benchmarks trigger investigation and Fortran code fixes where physics is not reproducible.

## Problem Frame

The verification ladder (rungs 1-4) validates solver machinery but leaves the "does the physics match published results?" question partially answered. G-factors, optical absorption spectra, strain shifts, and topological invariants lack automated regression with traceability to literature references. A methods paper needs systematic, reproducible benchmark results. The standard-star suite creates per-material scripts that each validate multiple published observables, generating a publication-ready benchmark table runnable via `ctest --test-dir build -L standard-star`.

## Requirements

Traced from origin `docs/brainstorms/2026-05-09-standard-star-benchmarks-requirements.md`.

### System definitions (R1-R3)

- **R1.** Seven standard-star material systems spanning distinct physics regimes (see table in origin).
- **R2.** Each system validates 2-5 independent published observables with named literature references and tolerance tiers.
- **R3.** Every validated observable has a published literature reference (author, year, source) and a named tolerance tier.

### Observable validation (R4-R7)

- **R4.** Analytical observables (band gap, spin-orbit splitting, Roth g-factor) validated at exact or near-exact tolerance.
- **R5.** Kane-model observables (effective masses) validated via parabolic fitting near Gamma, comparing against 2-band Kane formula `m*_kane = Eg / (EP + Eg)` at ~10% tolerance. NOT compared against Vurgaftman tabulated masses (8-band model is inherently non-parabolic; GaAs deviates 29%).
- **R6.** Confinement observables (subband spacing, band overlap, Stark shift, HH-LH splitting) validated against published or analytical references with stated tolerance per observable.
- **R7.** Topological observables (band inversion, Z2 invariant) validated against InAs/GaSb broken-gap QW (capability-gated: requires topologicalAnalysis to support real 8-band QW materials in Fu-Kane mode).

### Test infrastructure (R8-R11)

- **R8.** One Python verification script per system, named `verify_star_<material>.py`.
- **R9.** All scripts registered as ctest targets under a `standard-star` label.
- **R10.** Shared utility module (`tests/integration/star_helpers.py`) for running executables, parsing output, comparing values, formatting benchmark table rows.
- **R11.** Standard-stars complement the verification ladder. Overlapping observables cross-referenced, not duplicated.

### Publication output (R12-R14)

- **R12.** Each script outputs benchmark table rows: `| Material | Observable | Computed | Expected | Reference | Tolerance | Delta | Status |`.
- **R13.** Aggregation script concatenates per-system output into a publication-ready table.
- **R14.** Literature references embedded in each script and reproduced verbatim in benchmark table output.

## Key Technical Decisions

### KD1. Shared utility module is new infrastructure

`tests/support/` is Fortran-only. No shared Python utilities exist. Each existing verification script duplicates ~15-25 lines of parsing/comparison helpers. Create `tests/integration/star_helpers.py` with common functions: `parse_eigenvalues`, `parse_gfactor`, `parse_absorption`, `parse_topology_result`, `run_executable`, `compare_value`, `format_benchmark_row`.

**Why:** R10 explicitly requires this. Avoids duplicating boilerplate across 7 new scripts. Existing scripts remain unchanged (no retroactive refactoring).

### KD2. g-factor comparison against Roth analytical formula, not experiment

The 8-band Lowdin partitioning systematically underestimates g-factors vs experiment by 20-30% (GaAs: 8-band gives -0.315, experiment is -0.44). The Roth formula computed from the same parameters gives -0.315 — matching the 8-band result at near-exact tolerance. Comparing 8-band output against Roth analytical prediction validates the solver; comparing against experiment would always fail due to model limitations.

**Why:** Distinguishes solver correctness from model limitations. A failure against Roth means a code bug; a shortfall vs experiment is the expected 8-band behavior.

**Note on benchmarks.md Roth formula:** The Roth formula printed in `docs/reference/benchmarks.md` (line 206) has an error — `g0 - (2Ep/3)(1/Eg + 1/(Eg+Delta))` computes to ~-20.96 for GaAs, not the claimed -0.44. Standard-star scripts must use the correct Winkler 2003 Eq. (6.42): `g = 2 - 2*EP*Delta/(3*Eg*(Eg+Delta))` which gives -0.317 for GaAs (matching 8-band output of -0.315). The benchmarks.md error should be corrected in a separate commit.

### KD3. Effective mass comparison against 2-band Kane formula, not Vurgaftman

The 8-band CB effective mass for GaAs is ~0.0476 m0 (numerical), Vurgaftman tabulated is 0.067, and the 2-band Kane formula gives ~0.053. The 29% Vurgaftman deviation is the 8-band model's known non-parabolicity. The correct self-consistency check is against the Kane formula at ~10% tolerance. Observed deviations: GaAs 4.9%, InAs 1.2%, InSb 7.6%, GaSb 1.2%.

**Why:** Parabolic fitting near Gamma using numerical differentiation `d2E/dk2 = 2*(E(k1) - E(0))/k1^2` gives the 8-band parabolic limit. This must match the analytical Kane prediction from the same parameters.

### KD4. Literature values from independent published references, never from parameters.f90

The tautological R12 overlap check in the verification ladder (comparing code output against a value derived from the same parameter set) could never fail. For standard-star benchmarks, expected values must come from published literature (Vurgaftman 2001, Winkler 2003, Chuang 2003, Bastard 1981) and be cited in both the script and the benchmark table.

**Why:** A benchmark that reproduces the code's own input is a circular check. Meaningful validation requires independently sourced expected values.

**Exception:** R4 observables (Eg, DeltaSO at k=0) are exact by construction — the code must reproduce its own parameter values. These are self-consistency checks, not literature comparisons, and are already covered by verification ladder rung 1.

### KD5. Topological Z2 benchmark is capability-gated

The `topologicalAnalysis` executable's QSHE Fu-Kane mode constructs the real 8-band QW Hamiltonian for Z2 parity calculation. This path is currently untested with real 8-band materials (all existing topology configs use model Hamiltonians). Before writing the S5 Z2 test, verify this capability works. If it doesn't, defer the "topo" column for S5 to a follow-up.

**Why:** The origin document flags this as an open dependency. Building a test against unverified infrastructure wastes effort.

### KD6. Three tolerance tiers with explicit error-source documentation

| Tier | Tolerance | Error source | Observables |
|------|-----------|--------------|-------------|
| Exact | ~1e-12 (machine precision) | Parameter reproduction at k=0 | Eg, DeltaSO |
| Analytical | 0.5-10% | 8-band model non-parabolicity, Roth vs full Lowdin difference | Kane effective masses (~10%), Roth g-factor vs 8-band (~1%) |
| Numerical | 1-5% | Published computational/experimental reference uncertainty | Subband spacing, Stark shift, absorption edge, HH-LH splitting |

**Why:** Each tier documents what source of error it accounts for, making tolerance violations immediately attributable: parameter error, model limitation, or solver bug.

### KD7. Failing benchmarks trigger investigation and code fixes

When a benchmark fails, the phased implementation order provides diagnostic context:
- Phase 1 failures (bulk S1-S3) → parameter database or bulk Hamiltonian construction
- Phase 2 failures (QW S4) → confinement, FD stencils, QW Hamiltonian assembly
- Phase 3 failures (advanced QW S5-S6) → strain (Bir-Pikus), broken-gap alignment, topology
- Phase 4 failures (wire S7) → wire Hamiltonian, CSR assembly, wire geometry

Each failure must be investigated and resolved before proceeding to the next phase. Resolution may require Fortran code fixes in the solver, Hamiltonian constructor, or physics modules.

**Why:** Benchmarks that silently pass with inflated tolerances are worse than no benchmarks. Tight tolerances with investigation on failure create a trustworthy validation suite.

### KD8. CTest label is `standard-star`, not `verification`

The verification ladder uses `verification`. Standard-stars are a distinct category: per-material literature benchmarks vs per-complexity-level solver checks. A new label keeps `ctest -L standard-star` separate from `ctest -L verification` while allowing `ctest -L verification -L standard-star` for combined runs.

**Why:** Origin R9 specifies the `standard-star` label. Separate labels serve different audiences (solver developers vs methods-paper readers).

## Implementation Units

### U1: Shared utility module

**Files:**
- Create: `tests/integration/star_helpers.py`
- Reference: `tests/integration/verify_8band_rung1_bulk_k0.py` (pattern source)

**Description:**
Shared Python utilities for standard-star scripts. Extracts common patterns from existing verification scripts into a reusable module.

**Functions to provide:**
- `run_executable(exe_path, config_path, work_dir)` — copies config as `input.cfg`, runs executable, returns output directory path
- `parse_eigenvalues(filepath)` — returns `list[tuple[float, list[float]]]` of (k, eigenvalues)
- `parse_gfactor(filepath)` — returns `(gx, gy, gz)` tuple from `output/gfactor.dat`
- `parse_absorption(filepath)` — returns `list[tuple[float, float]]` of (energy, alpha) via numpy loadtxt
- `parse_topology_result(filepath)` — returns `dict[str, str]` of key-value pairs from comment-prefixed lines
- `compare_value(actual, expected, tolerance, name, unit="")` — returns `(passed, delta, row_dict)` with formatted benchmark row data
- `format_benchmark_row(material, observable, computed, expected, reference, tolerance, delta, status)` — returns markdown table row string
- `print_benchmark_header()` — prints `| Material | Observable | ... |` header

**Test scenarios:**
- `run_executable` correctly copies config, runs in tempdir, returns output path
- `parse_eigenvalues` handles multi-k files and single-k files
- `parse_gfactor` parses 3 space-separated floats
- `compare_value` correctly computes absolute and relative tolerance
- `format_benchmark_row` produces valid markdown table row

**Blocks:** U2-U8 (all subsequent scripts import from this module)

---

### U2: S1 — GaAs bulk standard-star

**Files:**
- Create: `tests/integration/verify_star_gaas_bulk.py`
- Reuse configs: `tests/regression/configs/bulk_gaas_k0.cfg`, `tests/regression/configs/gfactor_bulk_gaas_cb.cfg`, `tests/regression/configs/bulk_gaas_optics.cfg`
- Reference: `docs/reference/benchmarks.md` Section 1

**Description:**
Validates GaAs bulk observables against Vurgaftman 2001 and Winkler 2003.

**Observables and literature references:**

| Observable | Expected | Reference | Tolerance tier |
|-----------|----------|-----------|----------------|
| Eg | 1.519 eV | Vurgaftman 2001, Table I | Exact (1e-12) |
| DeltaSO | 0.341 eV | Vurgaftman 2001, Table I | Exact (1e-12) |
| m*_e (Kane) | Eg/(EP+Eg) = 1.519/(28.8+1.519) = 0.0501 m0 | 2-band Kane formula from code's own P,Eg | Analytical (~10%) |
| g* (Roth vs 8-band) | Roth formula from (Eg, DeltaSO, EP) | Winkler 2003, Eq. (6.42) | Analytical (~1%) |
| Absorption edge | ~1.519 eV (at Eg) | Vurgaftman 2001 + Elliott formula | Numerical (2%) |

**Cross-reference with verification ladder:** Eg and DeltaSO are covered by rung 1 (R1). Effective mass is covered by rung 2 (R6-R7). Standard-star adds: g-factor validation, optical absorption edge, and the publication-ready benchmark table format. Note the overlap in the script header.

**Test scenarios:**
- Eg at k=0 matches parameters.f90 within machine precision
- DeltaSO at k=0 matches parameters.f90 within machine precision
- Parabolic fit near Gamma gives m* within 10% of 2-band Kane prediction
- 8-band g-factor matches Roth analytical formula within ~1%
- Optical absorption onset occurs at Eg within 2%

---

### U3: S2 — InAs bulk standard-star

**Files:**
- Create: `tests/integration/verify_star_inas_bulk.py`
- Reuse configs: `tests/regression/configs/bulk_inas_k0.cfg`, `tests/regression/configs/gfactor_bulk_inasw_cb.cfg`
- Reference: `docs/reference/benchmarks.md` Section 1

**Description:**
Validates InAs bulk observables. Uses InAsW variant for g-factor (InAsW has Winkler parameters, used in existing g-factor configs).

**Observables:**

| Observable | Expected | Reference | Tolerance tier |
|-----------|----------|-----------|----------------|
| Eg | 0.417 eV | Vurgaftman 2001, Table I | Exact (1e-12) |
| m*_e (Kane) | 0.417/(21.5+0.417) = 0.0190 m0 | 2-band Kane formula | Analytical (~10%) |
| g* (Roth vs 8-band) | Roth formula from (Eg, DeltaSO, EP) | Winkler 2003, Eq. (6.42) | Analytical (~1%) |

**Cross-reference with verification ladder:** Eg and m* covered by rung 1 and rung 2. Standard-star adds: g-factor.

**Test scenarios:**
- Eg at k=0 matches parameters.f90 within machine precision
- Parabolic m* within 10% of Kane prediction (observed deviation: 1.2%)
- 8-band g-factor matches Roth formula within ~1%

---

### U4: S3 — InSb bulk standard-star

**Files:**
- Create: `tests/integration/verify_star_insb_bulk.py`
- Reuse configs: `tests/regression/configs/bulk_insb_k0.cfg`, `tests/regression/configs/bulk_insb_kx_dispersion.cfg`
- Create config: `tests/regression/configs/gfactor_bulk_insb_cb.cfg`
- Reference: `docs/reference/benchmarks.md` Section 1

**Description:**
Validates InSb bulk observables. InSb has extreme spin-orbit coupling (DeltaSO = 0.81 eV), making it the most stringent g-factor test. New g-factor config needed (modeled on `gfactor_bulk_gaas_cb.cfg`).

**Observables:**

| Observable | Expected | Reference | Tolerance tier |
|-----------|----------|-----------|----------------|
| Eg | 0.235 eV | Vurgaftman 2001, Table I | Exact (1e-12) |
| DeltaSO | 0.81 eV | Vurgaftman 2001, Table I | Exact (1e-12) |
| m*_e (Kane) | 0.235/(23.3+0.235) = 0.00997 m0 | 2-band Kane formula | Analytical (~10%) |
| g* (Roth vs 8-band) | Roth formula from (Eg, DeltaSO, EP) | Winkler 2003, Eq. (6.42) | Analytical (~1%) |

**Config `gfactor_bulk_insb_cb.cfg`:**
```
waveVector: k0
waveVectorStep: 0
confinement: 0
FDstep: 1
FDorder: 2
numLayers: 1
material1: InSb
numcb: 2
numvb: 6
ExternalField: 0
EFParams: 0.0005
whichBand: 0
bandIdx: 1
```

**Test scenarios:**
- Eg and DeltaSO at k=0 match parameters.f90 within machine precision
- Parabolic m* within 10% of Kane prediction (observed deviation: 7.6%)
- 8-band g-factor matches Roth formula within ~1%
- New config runs without error

---

### U5: S4 — GaAs/AlGaAs QW standard-star

**Files:**
- Create: `tests/integration/verify_star_gaas_algaas_qw.py`
- Reuse configs: `docs/benchmarks/qw_gaas_algaas.cfg`, `tests/regression/configs/qw_gaas_algaas_absorption.cfg`, `tests/regression/configs/sc_qcse_gaas_algaas_ef.cfg`
- Reference: `docs/reference/benchmarks.md` Sections 2, 6

**Description:**
Validates GaAs/AlGaAs QW observables: subband spacing, Stark shift, and optical absorption. Subband spacing and absorption use the established Al30Ga70As(200A)/GaAs(100A)/Al30Ga70As(200A) benchmark structure. Stark shift uses the existing QCSE config (`sc_qcse_gaas_algaas_ef.cfg`) which has a different structure: Al20Ga80As(460A)/GaAs(60A)/Al20Ga80As(460A) — the expected Stark shift value must correspond to this structure, not the 100A/30% Al structure.

**Observables:**

| Observable | Expected | Reference | Tolerance tier |
|-----------|----------|-----------|----------------|
| CB subband spacing (E1-E2) | 9.92 meV | benchmarks.md Section 2; nextnano cross-validation | Numerical (1%) |
| Stark shift at -70 kV/cm | TBD (must match 60A/20%Al QCSE config) | benchmarks.md Section 6 (value for corresponding structure) | Numerical (5%) |
| Absorption TE onset | At CB1 edge (~1.021 eV) | benchmarks.md Section 2 | Numerical (2%) |

**Cross-reference with verification ladder:** Subband spacing covered by rung 3 (R10). Standard-star adds: Stark shift, optical absorption, publication table.

**Test scenarios:**
- CB spacing (E1-E2) = 9.92 meV within 1% (absolute tolerance ~0.1 meV)
- QCSE Stark shift = -1.79 meV at -70 kV/cm within 5%
- TE absorption onset energy matches CB1 subband within 2%
- Absorption TE integral > TM integral (polarization sanity check from existing verify_qw_absorption_polarization.py)

---

### U6: S5 — InAs/GaSb QW standard-star

**Files:**
- Create: `tests/integration/verify_star_inas_gasb_qw.py`
- Reuse configs: `tests/regression/configs/qw_inasw_gasbw_broken_gap.cfg`, `tests/regression/configs/gfactor_qw_cb.cfg`
- Reference: `docs/reference/benchmarks.md` Section 3

**Description:**
Validates InAs/GaSb broken-gap QW. Covers band overlap, g-factor, and (capability-gated) topological Z2 invariant. Uses InAsW/GaSbW variants (Winkler parameters with defined EV/EC).

**Observables:**

| Observable | Expected | Reference | Tolerance tier |
|-----------|----------|-----------|----------------|
| g* CB (Roth vs 8-band) | Roth formula from (Eg_InAsW, DeltaSO, EP) | Winkler 2003 | Analytical (~1%) |
| Band overlap (EC-EV) | From published literature (not code parameters) | Search for published InAs/GaSb overlap value | Numerical (5%) |
| Z2 invariant | 1 (inverted regime) | Fu-Kane method | Numerical (exact value) |

**Capability gate for Z2:** Before writing the Z2 test, run `topologicalAnalysis` with the InAsW/GaSbW QW config in QSHE Fu-Kane mode. If it produces a valid result, include the Z2 test. If it fails or produces garbage, defer the topo column and document the gap.

**Literature sourcing for band overlap:** The 142 meV value from benchmarks.md is derived from the code's own EC(InAsW) and EV(GaSbW) — this is a tautological check (origin doc's learning). Instead, source the expected overlap from a published paper (e.g., Liu et al. PRL 2008, or a review article on InAs/GaSb broken-gap systems). The script should note that the code's parameter-derived overlap is a separate self-consistency check (covered by ladder rung 3), while the literature comparison is the standard-star benchmark.

**Test scenarios:**
- 8-band QW g-factor matches Roth analytical prediction within ~1%
- Band overlap agrees with published literature within 5%
- Z2 = 1 in inverted regime (capability-gated)
- Script notes ladder rung 3 overlap cross-reference

---

### U7: S6 — InAs/GaAs strained QW standard-star

**Files:**
- Create: `tests/integration/verify_star_inas_gaas_qw.py`
- Create configs: `tests/regression/configs/qw_inas_gaas_strained.cfg`, `tests/regression/configs/gfactor_qw_inas_gaas_strained.cfg`
- Reference: Vurgaftman 2001 Tables XIV-XV (deformation potentials); Chuang 2003 Chapter 4 (strained band edges)

**Description:**
Validates InAs/GaAs strained QW. Natural lattice mismatch (a_InAs=6.058 A vs a_GaAs=5.653 A, ~6.7% mismatch) provides strain. Uses Bir-Pikus analytical formulas as reference.

**Observables:**

| Observable | Expected | Reference | Tolerance tier |
|-----------|----------|-----------|----------------|
| g* CB (Roth vs 8-band) | Roth formula from strained parameters | Winkler 2003 | Analytical (~2%) |
| CB subband spacing | Finite-barrier analytical (InAs well in GaAs barrier) | Bastard finite-barrier formula (for strained well width) | Numerical (5%) |
| HH-LH splitting (strain) | 2*b*eps_biaxial (Bir-Pikus) | Vurgaftman 2001 deformation potentials + Chuang Ch. 4 | Analytical (2%) |

**Config `qw_inas_gaas_strained.cfg`:**
```
waveVector: k0
waveVectorStep: 0
confinement: 1
FDstep: 101
FDorder: 4
numLayers: 3
material1: GaAs 0 50 0
material2: InAs 50 70 0
material3: GaAs 70 120 0
numcb: 4
numvb: 12
ExternalField: 0
EFParams: 0.0
strain: T
strain_ref: GaAs
strain_solver: pardiso
```

**Config `gfactor_qw_inas_gaas_strained.cfg`:** Same structure with `waveVector: k0`, `waveVectorStep: 0`, `whichBand: 0`, `bandIdx: 1`.

**HH-LH splitting calculation:**
- Biaxial strain: `eps_xx = (a_barrier - a_well) / a_well`
- Bir-Pikus: `delta_EHH - delta_ELH = 2*b*eps_biaxial` where `b = -1.8` for InAs
- With `a_GaAs=5.65325`, `a_InAs=6.0583`: `eps = (5.65325-6.0583)/6.0583 = -0.0669`
- `HH-LH splitting = 2*(-1.8)*(-0.0669) = 0.241 eV` (large, easily detectable)

**Test scenarios:**
- 8-band g-factor matches Roth prediction within ~2%
- CB subband spacing reasonable for 20 nm InAs well with GaAs barriers
- HH-LH splitting matches Bir-Pikus analytical within 2%
- Strain config runs without error (stress test for Bir-Pikus sign convention)

**Failure diagnostic:** If HH-LH splitting is wrong, check `compute_bp_scalar` in `src/physics/strain_solver.f90`. The `P_eps = -av * Tr(eps)` sign flip is a known source of subtle errors (CLAUDE.md boundary).

---

### U8: S7 — InAs nanowire standard-star

**Files:**
- Create: `tests/integration/verify_star_inas_wire.py`
- Create configs: `tests/regression/configs/wire_inas_rectangle.cfg`, `tests/regression/configs/wire_inas_gfactor.cfg`
- Reference: Winkler 2003 (wire g-factor theory)

**Description:**
Validates InAs nanowire g-factor. Most complex system: wire geometry requires CSR sparse solver + commutator-based velocity matrices. Wire g-factor depends on transverse confinement via the commutator `[r_alpha, H]`.

**Observables:**

| Observable | Expected | Reference | Tolerance tier |
|-----------|----------|-----------|----------------|
| g* CB wire | Roth formula (bulk limit for large wire) | Winkler 2003 | Analytical (5%) |
| Wire eigenvalues → QW limit | Converge toward bulk InAs as wire size increases | Internal consistency | Numerical (5%) |

**Config `wire_inas_rectangle.cfg`:**
```
waveVector: k0
waveVectorMax: 0.01
waveVectorStep: 2
confinement: 2
FDstep: 1
FDorder: 2
numLayers: 1
wire_nx: 11
wire_ny: 11
wire_dx: 5.0
wire_dy: 5.0
wire_shape: rectangle
wire_width: 55.0
wire_height: 55.0
numRegions: 1
region: InAs  0.0  100.0
numcb: 4
numvb: 8
ExternalField: 0  EF
EFParams: 0.0
SC: 0
feast_emin: -1.5
feast_emax: 2.0
feast_m0: -1
```

**Config `wire_inas_gfactor.cfg`:** Same structure with `whichBand: 0`, `bandIdx: 1`, `EFParams: 0.0005` (magnetic perturbation).

**Wire → bulk convergence check:** For a large wire (e.g., 50x50 nm), the ground-state eigenvalue should approach the bulk InAs band edge. This reuses the wire-convergence pattern from verification ladder rung 4 (R14).

**Test scenarios:**
- Wire g-factor is computed without error (CSR + commutator velocity path)
- g* is physically reasonable (same sign as bulk, magnitude within factor of 2 of Roth bulk value)
- Large-wire eigenvalue converges toward bulk InAs band edge

---

### U9: CTest registration and aggregation

**Files:**
- Modify: `tests/CMakeLists.txt`
- Create: `tests/integration/aggregate_star_benchmarks.py`

**Description:**
Register all 7 standard-star scripts as ctest targets under `standard-star` label. Create aggregation script for publication-ready benchmark table.

**CTest registration pattern** (per script):
```cmake
add_test(
    NAME standard_star_<material>
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_star_<material>.py
        ${CMAKE_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(standard_star_<material>
    PROPERTIES LABELS "standard-star"
)
```

**Aggregation script:**
- Runs each of the 7 standard-star scripts with `--benchmark-only` flag (no ctest exit codes)
- Captures benchmark table rows from stdout
- Concatenates into a single markdown table with header
- Writes to `docs/reference/standard_star_benchmarks.md` or stdout

**Test scenarios:**
- `ctest --test-dir build -L standard-star` discovers and runs all 7 tests
- Each test passes with all observables within tolerance
- `python3 tests/integration/aggregate_star_benchmarks.py` produces a complete markdown table
- Aggregated table has correct column headers per R12

**Blocks:** Nothing (last unit, registers all previous work)

---

## Failure Diagnostics

When a benchmark fails, use the phased structure to narrow the root cause:

| Phase | System | Failing observable | Likely root cause | Investigation target |
|-------|--------|-------------------|-------------------|---------------------|
| 1 | GaAs/InAs/InSb bulk | Eg, DeltaSO | Parameter database | `src/core/parameters.f90` |
| 1 | GaAs/InAs/InSb bulk | m* (Kane) | Bulk Hamiltonian or dispersion | `src/physics/hamiltonianConstructor.f90` |
| 1 | GaAs/InAs/InSb bulk | g* (Roth) | Lowdin partitioning | `src/physics/gfactor_functions.f90` |
| 2 | GaAs/AlGaAs QW | Subband spacing | Confinement, FD stencils | `src/physics/hamiltonianConstructor.f90` QW path |
| 2 | GaAs/AlGaAs QW | Stark shift | Electric field implementation | `src/apps/main.f90` EF setup |
| 2 | GaAs/AlGaAs QW | Absorption edge | Optical spectra | `src/physics/optical_spectra.f90` |
| 3 | InAs/GaSb QW | Band overlap | Band alignment, W-variant EV/EC | `src/core/parameters.f90` GaSbW/InAsW |
| 3 | InAs/GaSb QW | Z2 | Topological analysis, Fu-Kane | `src/physics/topological_analysis.f90` |
| 3 | InAs/GaAs QW | HH-LH splitting | Bir-Pikus sign convention | `src/physics/strain_solver.f90` `compute_bp_scalar` |
| 4 | InAs wire | g* wire | Wire velocity matrices, CSR | `src/physics/hamiltonianConstructor.f90` `build_velocity_matrices` |
| 4 | InAs wire | Eigenvalue convergence | Wire Hamiltonian, geometry | `src/physics/hamiltonian_wire.f90` |

**Investigation workflow:**
1. Identify failing phase and observable
2. Check the diagnostic table for likely root cause
3. Add targeted debugging to the failing Fortran module
4. Fix the code error
5. Re-run the failing benchmark to confirm the fix
6. Re-run all previous phases to confirm no regression
7. Document the fix in `docs/solutions/`

## Scope Boundaries

### In scope
- 7 Python verification scripts with shared utility module
- ~5 new executable config files
- CTest registration under `standard-star` label
- Publication-ready benchmark table aggregation
- Fortran code fixes when benchmarks reveal solver bugs

### Out of scope
- Refactoring existing verification ladder scripts to use star_helpers.py
- Performance benchmarking (time-to-solution)
- Adding new Fortran physics capabilities (only fixing existing ones)
- HgTe/CdTe or II-VI materials (not in parameters.f90)
- Richardson convergence testing (separate ideation item)
- Validation coverage matrix (separate ideation item)
- Executable lecture-test pairs (separate ideation item)
- Cross-code validation with nextnano (too expensive)

### Deferred to follow-up work
- S5 topological Z2 test if capability gate fails
- Additional observables for systems with thin literature coverage
- Retroactive refactoring of existing verification scripts to use star_helpers.py
- Adding S8-S10 for additional materials (InP, GaP, etc.)

## Dependencies / Assumptions

- All four executables (`bandStructure`, `gfactorCalculation`, `opticalProperties`, `topologicalAnalysis`) must build successfully from the current branch
- Existing configs in `tests/regression/configs/` are correct and produce valid output
- Python 3 with numpy is available (already required by existing verification scripts)
- Published reference values exist for all targeted observables (Vurgaftman 2001, Winkler 2003, Chuang 2003, benchmarks.md)
- The Roth g-factor formula provides a reliable analytical reference for 8-band g-factor self-consistency
- The 2-band Kane formula `m* = Eg/(EP+Eg)` provides a reliable analytical reference for 8-band effective mass
- The Bir-Pikus sign convention in `compute_bp_scalar` is correct (CLAUDE.md: NEVER change)
- Wire velocity matrices via `build_velocity_matrices` are correct (commutator-based)
- InAs/GaAs strained QW works with existing strain implementation (no new Fortran strain features needed)
- S5 topological Z2 requires Fu-Kane mode to work with real 8-band materials (unverified, gated)
- **Material variant consistency:** S2 (InAs bulk) uses non-W InAs (EP=21.5, Vurgaftman parameters). S5 (InAs/GaSb QW) and S6 (InAs/GaAs g-factor) use InAsW (EP=22.2, Winkler parameters). S7 (InAs wire) uses non-W InAs. The EP difference (3.2%) affects Roth g-factor predictions. Scripts must use the variant matching each config's existing material specification. Bulk dispersion configs use InAs; QW/gfactor configs use InAsW where existing configs already specify it.

## Deferred / Open Questions

### From 2026-05-09 review

- **S5 band overlap has no concrete published value** — The plan states the expected band overlap should come from "published literature (not code parameters)" but does not specify a concrete value or citation. Implementation requires literature research to find a published InAs/GaSb overlap value (e.g., Liu et al. PRL 2008, or a review article). Without this, the S5 overlap benchmark cannot be implemented. (feasibility, adversarial, confidence 100)

- **S6 InAs/GaAs 6.7% mismatch may exceed Bir-Pikus validity** — The lattice mismatch between InAs (6.058 A) and GaAs (5.653 A) is ~6.7%, which exceeds the ~1-2% range where Bir-Pikus linear elasticity is accurate. The resulting HH-LH splitting of 0.241 eV is very large and may be outside the perturbative regime. Tradeoff: this is the only strained QW system available with existing materials in `parameters.f90` (no InGaAs ternary). Accepting the large strain tests the code at its limits but may produce physically inaccurate results that aren't cleanly attributable to bugs vs model limitations. (feasibility, adversarial, confidence 75)
