---
title: 8-Band Verification Ladder
type: feat
status: active
date: 2026-05-08
deepened: 2026-05-08
origin: docs/brainstorms/2026-05-08-8band-verification-ladder-requirements.md
---

# 8-Band Verification Ladder

## Summary

Implements a 4-rung validation hierarchy for the 8-band k.p Hamiltonian using Python verification scripts that run the existing Fortran executables and check physical correctness. Each rung covers one complexity level: bulk k=0 structural checks (R1-R5), bulk dispersion effective masses (R6-R9), QW subband energies (R10, R12-R13), and wire internal consistency (R14-R16). Resolves review findings: R10 uses benchmarks.md self-consistency with matching config, R7 compares against Kane model prediction (not Vurgaftman) at 10%, R11 dropped (no independent nextnano data), R15 tolerance aligned with existing test at 1e-8, R14 uses constant grid spacing to avoid confounding FD error. No Fortran source changes.

---

## Problem Frame

The codebase has 24 regression tests and 9 integration verification scripts, but the coverage is uneven and undocumented as a diagnostic hierarchy. When a test fails, the developer must manually trace whether the error is in parameters, Hamiltonian construction, FD stencils, confinement, or sparse assembly. The verification ladder creates a structured progression where each rung validates one layer, making residual errors unambiguously attributable.

---

## Requirements

- R1. For each canonical material, compute 8x8 bulk Hamiltonian at k=0 and verify eigenvalues match band edges (EV, EC, Delta_SO) from parameters.f90 within machine precision.
- R2. Verify eigenvalue degeneracies: HH+LH (4x degenerate at k=0), SO (2x), CB (2x).
- R3. Verify basis ordering: bands 1-4 valence (HH, LH, LH, HH), 5-6 split-off, 7-8 conduction — via eigenfunction weight inspection.
- R4. Verify T_d symmetry at k=0: eigenvalues are exactly the band edges (indirect check — if any k-dependent coupling existed, eigenvalues would deviate).
- R5. Verify eigenfunction normalization (unit sum of squares).
- R6. For each canonical material, compute bulk dispersion along one k-direction, extract effective masses by parabolic fitting near Gamma.
- R7. (Revised v2) Extract effective masses numerically from dispersion curves and compare against the Kane model prediction m* = Eg/(EP + Eg) within 10% (accounts for higher-order band-mixing corrections beyond the simple 2-band Kane formula). The 8-band code produces non-parabolic dispersion — GaAs gives ~0.046 m0 vs Vurgaftman 0.067 — so direct comparison against Vurgaftman is informational only (see R8/R9).
- R8. Report the relative deviation between numerically extracted masses and both the Kane model prediction and Vurgaftman 2001 tabulated values. Both are informational.
- R9. (Revised v2) Report Vurgaftman deviation as informational — the 8-band model intentionally includes band-mixing effects that make CB masses lighter than single-band experimental values. Do not flag as a parameter-set inconsistency; this is a known model feature documented in benchmarks.md and the existing verify_bulk_benchmarks.py (which prints "INFO: Not a failure -- non-parabolic dispersion is a feature of the 8-band model").
- R10. (Revised v2) For a GaAs/AlGaAs QW, verify subband energies match benchmarks.md values (CB spacing = 9.92 meV for 100A GaAs well) using the exact config that produced those values (`docs/benchmarks/qw_gaas_algaas.cfg`, FDstep=101, FDorder=2).
- R11. (Revised) Merge with R10 — benchmarks.md does not contain independent nextnano data for the standard GaAs/AlGaAs QW structure. R10 already serves as the regression guard. R11 is dropped.
- R12. For InAs/GaSbW broken-gap QW, verify band overlap (~142 meV) against benchmarks.md.
- R13. Verify QW subband degeneracies and ordering are physically correct.
- R14. For a GaAs rectangular wire, verify eigenvalues converge toward QW result as transverse dimensions increase.
- R15. For the same wire, verify dense (LAPACK) and sparse (CSR/FEAST) solver paths agree within 1e-8 (matching the existing test_wire_dense_sparse_consistency.sh tolerance).
- R16. Verify wire eigenvalue count matches expected number given grid size and numcb/numvb.
- R17. 5 materials tested at rungs 1-2 (GaAs, InAs, InSb, GaSb, GaSbW — GaSbW included because its W-variant parameters are independently defined and QW rungs need W-variants). Bulk mode does not use EV/EC (only Eg, DeltaSO on the diagonal), so GaSb non-W works correctly for rungs 1-2. Material-geometry mapping for rungs 3-4: GaAs/AlGaAs for R10, InAs/GaSbW for R12, GaAs for R14-R16.
- R18. Per-requirement tolerances: machine precision (rung 1), 10% Kane model comparison (rung 2), R10 (benchmarks.md values with matching config), R12 (benchmarks.md, 1%), R14 (1%), R15 (1e-8).
- R19. Validation results documented in benchmarks.md.

**Origin acceptance examples:** AE1 (R1-R5), AE2 (R6-R8, revised v2 — origin specified 1e-8 analytical and 0.1% Vurgaftman thresholds, replaced by 10% Kane model comparison per P0 review finding that 8-band GaAs mass is ~0.046 vs Vurgaftman 0.067), AE3 (R10, revised v2 — uses benchmarks.md config, R11 dropped), AE4 (R14-R15).

---

## Scope Boundaries

- No changes to Fortran source code — all validation is external Python scripts
- No CSR structural testing (ideation item #2)
- No Richardson convergence testing (ideation item #6)
- No strain validation (ideation item #7)
- No pFUnit unit tests for Hamiltonian internals
- Temperature-dependent validation (all tests at default temperature)
- g-factor validation (covered by Standard-Star Benchmarks ideation item)
- Topological invariants (separate module)
- Expanding beyond 4 canonical materials

### Deferred to Follow-Up Work

- Additional materials beyond GaAs, InAs, InSb, GaSb/AlSb (expand later)
- Automated regression-cause attribution (the ladder identifies the rung, not the root cause)

---

## Context & Research

### Relevant Code and Patterns

- `tests/integration/verify_bulk_benchmarks.py` — existing bulk validation (Eg, Delta_SO, effective mass for GaAs/InAs at 1% tolerance)
- `tests/integration/verify_qw_benchmarks.py` — QW eigenvalue sanity checks
- `tests/integration/verify_qw_state_character.py` — broken-gap QW state ordering
- `tests/regression/compare_output.py` — golden-file numeric comparison (1e-8 relative, 1e-14 absolute near-zero)
- `tests/integration/test_*.sh` — shell scripts that run executables, check exit codes, compare output
- `docs/reference/benchmarks.md` — 6 validated benchmark systems with quantitative values
- `docs/reference/input-reference.md` — complete config file format reference
- `docs/reference/output-reference.md` — output file format reference
- `src/physics/hamiltonianConstructor.f90` — `ZB8bandBulk` (line 593) builds 8x8 bulk; `ZB8bandQW` (line 42) builds 8N×8N QW
- `src/core/parameters.f90` — material parameters; GaSb non-W has EV=0/EC=0 (unset)

### Institutional Learnings

- Eigenvalue comparison alone is insufficient for catching Hamiltonian construction bugs — CSR structural testing is needed (docs/solutions/topological-magnetic-index-logic-errors). The verification ladder complements but does not replace structural tests.
- Wire regression tests only compare k=1 eigenvalues due to FEAST non-determinism at higher k-points.

### External References

- Vurgaftman 2001, J. Appl. Phys. 89, 5815 — tabulated III-V material parameters (band gaps, effective masses, Luttinger parameters)
- Winkler 2003, Spin-Orbit Coupling Effects in Two-Dimensional Electron and Hole Systems — W-variant parameter sets

---

## Key Technical Decisions

- **R10 replacement: benchmarks.md self-consistency check.** The Bastard infinite-barrier formula is single-band; benchmarks.md documents a 5.4x discrepancy vs 8-band. Replacing with the code's own documented values (CB spacing 9.92 meV for 100A GaAs/AlGaAs QW) validates that the implementation hasn't regressed from its validated state. Uses the exact config that produced the benchmark (`docs/benchmarks/qw_gaas_algaas.cfg`, FDstep=101, FDorder=2). This is a regression guard — the original validation against nextnano is not independently available for this specific QW structure, so R11 (nextnano comparison) is dropped.
- **R7 rework: Kane model self-consistency, not Vurgaftman comparison.** The 8-band code produces non-parabolic CB dispersion — GaAs gives ~0.046 m0 vs Vurgaftman 0.067 (31% deviation). This is a known feature of the 8-band model (documented in benchmarks.md and verify_bulk_benchmarks.py: "INFO: Not a failure -- non-parabolic dispersion is a feature of the 8-band model"). Comparing against Vurgaftman at any tight tolerance would always fail. Instead, compare against the Kane model prediction m* = Eg/(EP + Eg) within 10% tolerance, which accounts for higher-order band-mixing corrections. Report Vurgaftman deviation as informational only.
- **R11 dropped.** benchmarks.md section 2 contains only the code's own computed values (CB1=1.02133, CB2=1.03125, spacing=9.92 meV). There is no independent nextnano data for this specific GaAs/AlGaAs QW structure. R10 already serves as the regression guard. An independent physics validation for this QW is a gap to address in a future iteration.
- **R15 tolerance aligned with existing test.** Relaxed from 1e-10 to 1e-8 to match the existing `test_wire_dense_sparse_consistency.sh` which uses 1e-8 and passes.
- **Python-only, no Fortran changes.** The verification scripts run the existing executables and parse their output. This avoids touching production code and matches the existing verify_*.py pattern.
- **Single material database per script.** Each script contains a hardcoded dict of parameters from parameters.f90 and Vurgaftman 2001 rather than parsing the Fortran source at runtime. Simpler, more portable, and the values change rarely.

---

## Open Questions

### Resolved During Planning

- **k-range for parabolic fitting:** Use k = 0.001 to 0.05 nm⁻¹ for all materials, with adaptive fitting that narrows the range until R² > 0.9999. InSb and InAs (narrow gap) will naturally use a smaller range due to stronger non-parabolicity.
- **R11 nextnano reference:** No independent nextnano data exists for the standard GaAs/AlGaAs QW in benchmarks.md. R11 dropped; R10 serves as the regression guard. Future iteration can add independent QW validation.
- **Wire dimensions for QW convergence:** Use 63×63 Å (21×21, dx=dy=3 Å), 93×93 Å (31×31, dx=dy=3 Å), and 123×123 Å (41×41, dx=dy=3 Å) — constant grid spacing, varying grid points and physical size. Matches the pattern in wire_gaas_rectangle.cfg vs wire_gaas_rectangle_large.cfg.
- **Python script format:** Standalone scripts following the existing `verify_*.py` pattern — no integration with compare_output.py needed.
- **GaSb non-W bulk:** Bulk mode (ZB8bandBulk) does not use EV/EC — only Eg and DeltaSO go on the diagonal. GaSb non-W works correctly for rung 1 bulk tests. QW rungs use GaSbW which has EV/EC defined.

### Deferred to Implementation

- Exact parabolic fitting algorithm (polynomial degree, fitting window selection) — implementation-time tuning
- Whether to include VB effective masses (m*_hh, m*_lh) in rung 2 or CB only — the VB masses are anisotropic and more complex to validate
- Whether the dense-sparse comparison (R15) needs a separate wire config or can reuse the existing wire_gaas_rectangle.cfg
- Exact Kane model prediction values for each material (need EP from parameters.f90)

---

## Implementation Units

- U1. **Rung 1 — Bulk k=0 Structural Validation**

**Goal:** Create configs and a Python verification script that validates bulk eigenvalues, degeneracies, basis ordering, T_d symmetry, and eigenfunction normalization at k=0 for all 4 canonical materials.

**Requirements:** R1, R2, R3, R4, R5. Covers AE1.

**Dependencies:** None.

**Files:**
- Create: `tests/integration/verify_8band_rung1_bulk_k0.py`
- Create: `tests/regression/configs/bulk_inas_k0.cfg`
- Create: `tests/regression/configs/bulk_insb_k0.cfg`
- Create: `tests/regression/configs/bulk_gasb_k0.cfg`
- Create: `tests/regression/configs/bulk_gasbw_k0.cfg`
- Reference: `tests/regression/configs/bulk_gaas_kx.cfg` (existing pattern)

**Approach:**
Create minimal bulk k=0 configs for each material (waveVector: k0, waveVectorMax: 0, confinement: 0). The Python script runs `bandStructure` for each config, parses `output/eigenvalues.dat` and `output/parts.dat`, then checks:

1. Eigenvalues match expected band edges in ascending order: [-DeltaSO, -DeltaSO, 0, 0, 0, 0, Eg, Eg] at machine precision (1e-12 relative)
2. Degeneracy: pairs 1-2 are degenerate (SO), elements 3-6 are degenerate (HH/LH, 4-fold), pairs 7-8 are degenerate (CB)
3. Basis ordering: eigenfunction weights from parts.dat show correct band character — bands 1-4 have HH+LH weight > 0.9, bands 5-6 have SO weight > 0.9, bands 7-8 have CB weight > 0.9
4. T_d symmetry (indirect): eigenvalues are exactly the band edges — any non-zero off-diagonal coupling would cause deviations
5. Eigenfunction normalization: sum of squares of each eigenfunction ≈ 1.0 within 1e-10

Material database in the script includes: GaAs (Eg=1.519, DeltaSO=0.341), InAs (Eg=0.417, DeltaSO=0.390), InSb (Eg=0.235, DeltaSO=0.810), GaSb (Eg=0.812, DeltaSO=0.760), GaSbW (Eg=0.812, DeltaSO=0.760). Note: bulk mode does not use EV/EC (only Eg, DeltaSO on the diagonal), so GaSb non-W works correctly despite having unset EV/EC.

**Patterns to follow:**
- `tests/integration/verify_bulk_benchmarks.py` — script structure, material dict, eigenvalue parsing
- `tests/integration/verify_qw_state_character.py` — eigenfunction weight parsing

**Test scenarios:**
- Happy path: GaAs bulk k=0 → eigenvalues [-0.341,-0.341,0,0,0,0,1.519,1.519] at machine precision (ascending order)
- Happy path: InSb bulk k=0 → verify large DeltaSO (0.810) eigenvalue splitting
- Edge case: GaSb (non-W, EV=0) → bulk still works correctly
- Degeneracy check: all materials → HH+LH(4x), SO(2x), CB(2x) degeneracy pattern
- Basis ordering: eigenfunction weights identify HH, LH, SO, CB character correctly
- Normalization: all eigenfunctions have unit norm

**Verification:**
- Script exits with code 0 when all checks pass, non-zero on any failure
- Output includes per-material pass/fail with specific values and tolerances

---

- U2. **Rung 2 — Bulk Dispersion Effective Mass Validation**

**Goal:** Create configs and a Python verification script that computes bulk dispersion, extracts effective masses by parabolic fitting, and compares against Vurgaftman 2001 tabulated values.

**Requirements:** R6, R7, R8, R9 (revised v2). Covers AE2.

**Dependencies:** None.

**Files:**
- Create: `tests/integration/verify_8band_rung2_dispersion.py`
- Create: `tests/regression/configs/bulk_gaas_kx_dispersion.cfg`
- Create: `tests/regression/configs/bulk_inas_kx_dispersion.cfg`
- Create: `tests/regression/configs/bulk_insb_kx_dispersion.cfg`
- Create: `tests/regression/configs/bulk_gasb_kx_dispersion.cfg`
- Reference: `tests/integration/verify_bulk_benchmarks.py` (existing mass extraction pattern)

**Approach:**
Create bulk k-sweep configs (waveVector: kx, waveVectorMax: 0.05, waveVectorStep: 51) for each material. The Python script runs `bandStructure` for each config, parses `output/eigenvalues.dat`, then:

1. Identifies the CB band (highest eigenvalue pair) and extracts E(k) near k=0
2. Fits a parabola E(k) = E₀ + ℏ²k²/(2m*) to the first few k-points (adaptive range until R² > 0.9999)
3. Extracts m* from the fit coefficient
4. Computes the Kane model prediction m*_kane = Eg/(EP + Eg) using EP from parameters.f90 for each material
5. Compares numerically extracted m* against m*_kane within 10% tolerance — accounts for higher-order corrections beyond the 2-band Kane formula
6. Reports Vurgaftman deviation as informational (not a test failure) — the 8-band model intentionally includes band-mixing that makes CB masses lighter than experimental values

Also extracts HH, LH, SO band curvatures where possible (noting that VB masses are anisotropic — only the kx direction is measured).

The effective mass is extracted in units of m₀ using: m* = ℏ²/(2 × d²E/dk²) where ℏ² = 7.61996 eV·Å² and k is in Å⁻¹.

Kane model predictions using EP from parameters.f90:
- GaAs: m*_kane = 1.519/(28.8 + 1.519) = 0.0501 (8-band gives ~0.046, ~8% deviation from Kane)
- InAs: m*_kane = 0.417/(21.5 + 0.417) = 0.0190
- InSb: m*_kane = 0.235/(23.3 + 0.235) = 0.0100
- GaSb: m*_kane = 0.812/(27.0 + 0.812) = 0.0292

Vurgaftman informational targets: GaAs 0.067, InAs 0.026, InSb 0.0135, GaSb 0.041.

**Patterns to follow:**
- `tests/integration/verify_bulk_benchmarks.py` — effective mass extraction (lines with numpy polyfit)
- `scripts/verify_landau_levels.py` — material database pattern

**Test scenarios:**
- Happy path: GaAs CB dispersion → m*_e extracted and compared against Kane model prediction Eg/(EP+Eg) within 10%
- Happy path: InAs CB dispersion → same Kane model comparison
- Edge case: InSb (narrow gap, strong non-parabolicity) → adaptive fitting narrows range
- Edge case: GaSb → same Kane model comparison
- Informational: report Vurgaftman deviation for each material (expected ~30% for GaAs)

**Verification:**
- All 4 materials pass with < 10% deviation from Kane model prediction
- Adaptive fitting produces R² > 0.9999 for all materials
- Report includes k-range used, fit quality, extracted mass, Kane prediction, and Vurgaftman deviation for each material

---

- U3. **Rung 3 — QW Subband Energy Validation**

**Goal:** Create QW configs and a Python verification script that validates subband energies against benchmarks.md and published results, plus subband ordering and degeneracies.

**Requirements:** R10 (revised v2), R12, R13. Covers AE3 (revised).

**Dependencies:** None.

**Files:**
- Create: `tests/integration/verify_8band_rung3_qw.py`
- Reference: `docs/benchmarks/qw_gaas_algaas.cfg` (exact config that produced 9.92 meV benchmark — FDstep=101, FDorder=2, 3-layer Al30Ga70As/GaAs/Al30Ga70As)
- Reference: `tests/regression/configs/qw_inasw_gasbw_broken_gap.cfg` (existing broken-gap config — 3 layers: AlSbW/InAsW/GaSbW)
- Reference: `tests/integration/verify_qw_benchmarks.py` (existing QW checks)
- Reference: `tests/integration/verify_qw_state_character.py` (existing state character checks)

**Approach:**
Use the benchmark config for R10 (not the regression config — different parameters), reuse existing broken-gap config for R12. The script:

1. **R10 (revised v2):** Copies `docs/benchmarks/qw_gaas_algaas.cfg` to a temp dir, runs `bandStructure`, extracts CB subband energies, verifies E2-E1 spacing matches benchmarks.md value (9.92 meV) within a tight tolerance (0.1 meV). Uses the exact config that produced the benchmark to avoid FD order/grid resolution discrepancies.

2. **R11 is dropped.** No independent nextnano reference exists for this QW structure in benchmarks.md.

3. **R12:** Reuses existing `tests/regression/configs/qw_inasw_gasbw_broken_gap.cfg` (3-layer AlSbW/InAsW/GaSbW structure, matching the benchmarks.md setup). Runs it, extracts the CB-VB overlap from eigenvalues. Verifies overlap ≈ 142 meV as documented in benchmarks.md.

4. **R13:** For both QW configs at k=0, verifies:
   - CB subbands are above VB subbands in energy (standard QW only — broken-gap has VB above CB)
   - Spin degeneracy: each subband has a Kramers partner within 1e-6 eV
   - Eigenfunction weights from parts.dat show correct band character

**Patterns to follow:**
- `tests/integration/verify_qw_benchmarks.py` — QW eigenvalue sanity checks
- `tests/integration/verify_qw_state_character.py` — broken-gap state ordering
- `docs/benchmarks/qw_gaas_algaas.cfg` — exact benchmark config for R10

**Test scenarios:**
- Happy path: GaAs/AlGaAs QW CB spacing matches benchmarks.md (9.92 meV) using benchmark config
- Happy path: InAs/GaSbW broken-gap overlap ≈ 142 meV using existing 3-layer config
- Degeneracy: QW eigenvalues at k=0 are spin-degenerate (pairs within 1e-6 eV)
- Ordering: CB eigenvalues above VB (standard QW), correct band character from parts.dat
- Edge case: broken-gap regime — some VB states above CB states (verified by state character)

**Verification:**
- GaAs/AlGaAs CB spacing within 0.1 meV of benchmarks.md value (9.92 meV)
- InAs/GaSbW overlap within 1% of 142 meV
- All subbands correctly ordered and degenerate

---

- U4. **Rung 4 — Wire Internal Consistency Validation**

**Goal:** Create wire configs and a Python verification script that validates wire-to-QW convergence, dense vs sparse solver agreement, and eigenvalue count.

**Requirements:** R14, R15, R16. Covers AE4.

**Dependencies:** None.

**Files:**
- Create: `tests/integration/verify_8band_rung4_wire.py`
- Create: `tests/regression/configs/wire_gaas_31x31.cfg` (93×93 Å, dx=dy=3 Å)
- Create: `tests/regression/configs/wire_gaas_41x41.cfg` (123×123 Å, dx=dy=3 Å)
- Reference: `tests/regression/configs/wire_gaas_rectangle.cfg` (existing 21×21, 63×63 Å, dx=dy=3 Å)
- Reference: `tests/regression/configs/wire_gaas_rectangle_large.cfg` (existing 31×31 pattern)

**Approach:**
Create GaAs rectangular wire configs at multiple sizes to test convergence. The script:

1. **R14:** Runs wire configs with increasing transverse dimensions while keeping grid spacing constant (dx=dy=3 Å). Configs: 63×63 Å (21×21), 93×93 Å (31×31), 123×123 Å (41×41). Extracts ground-state eigenvalue at k=0, verifies it converges toward the bulk band edge (less confinement energy as wire gets wider). The convergence trend should be monotonic: larger wire → lower ground state energy. Keeping dx/dy constant avoids confounding FD discretization error with physical size changes.

2. **R15:** For the 21×21 grid wire config (small enough for dense LAPACK solve), runs with both dense solver path and sparse solver path (CSR/FEAST), verifies eigenvalues agree within 1e-8 (matching existing test_wire_dense_sparse_consistency.sh tolerance).

3. **R16:** Verifies that the number of eigenvalues returned equals numcb + numvb and that the wire grid produces the expected matrix size (8 × nx × ny).

Wire configs use: confinement=2, wire_shape=rectangle, constant dx=dy=3 Å, varying nx=ny (21, 31, 41) to change physical dimensions while keeping grid spacing constant.

**Patterns to follow:**
- `tests/integration/test_wire_bandstructure.sh` — wire test execution (only k=1 comparison due to FEAST non-determinism)
- `tests/regression/configs/wire_gaas_rectangle.cfg` — wire config format

**Test scenarios:**
- Happy path: wire eigenvalues at k=0 for 63×63 Å grid are physically reasonable
- Convergence: ground-state energy decreases monotonically as wire dimensions increase (constant grid spacing)
- Dense vs sparse: eigenvalues from both solver paths agree within 1e-8 (for 21×21 grid)
- Eigenvalue count: returned eigenvalues = numcb + numvb for each config

**Verification:**
- Monotonic convergence of ground-state energy with increasing wire size
- Dense-sparse agreement within 1e-8 (where both paths are available)
- Correct eigenvalue count for all wire configs

---

- U5. **Shell Scripts and CMake Registration**

**Goal:** Create integration shell scripts for all new configs and register them in CMakeLists.txt.

**Requirements:** R17, R19 (partial). R11 is dropped — R10 handles QW regression guard.

**Files:**
- Create: `tests/integration/test_8band_rung1_bulk_k0.sh`
- Create: `tests/integration/test_8band_rung2_dispersion.sh`
- Create: `tests/integration/test_8band_rung3_qw.sh`
- Create: `tests/integration/test_8band_rung4_wire.sh`
- Modify: `tests/CMakeLists.txt` — register new tests with "regression" label
- Reference: `tests/integration/test_bulk_bandstructure.sh` (existing pattern)

**Approach:**
Each shell script follows the existing pattern: create temp dir, copy config, run executable, check exit code, compare output, run verification script. One script per rung that loops over all materials/configs for that rung.

Register all new tests in CMakeLists.txt with the "regression" label so they run with `ctest -L regression`.

**Patterns to follow:**
- `tests/integration/test_bulk_bandstructure.sh` — template for shell scripts
- `tests/integration/test_wire_bandstructure.sh` — wire-specific handling

**Test scenarios:**
- Happy path: each shell script runs cleanly and exits 0
- Each test appears in `ctest -N` output with correct label
- `ctest -L regression` runs all new tests alongside existing ones

**Verification:**
- All new tests pass with `ctest --test-dir build -L regression`
- Test output shows per-material pass/fail from verification scripts

---

- U6. **Documentation Update**

**Goal:** Document validation results in benchmarks.md and update the requirements doc with resolved open questions.

**Requirements:** R19.

**Dependencies:** U1, U2, U3, U4, U5.

**Files:**
- Modify: `docs/reference/benchmarks.md` — add verification ladder results section
- Modify: `docs/brainstorms/2026-05-08-8band-verification-ladder-requirements.md` — update deferred open questions

**Approach:**
After all verification scripts pass, add a "Verification Ladder" section to benchmarks.md documenting:
- Rung 1 results: eigenvalue tables for each material, degeneracy confirmation
- Rung 2 results: effective mass values, parabolic fit quality, Vurgaftman deviations
- Rung 3 results: QW subband energies, broken-gap overlap
- Rung 4 results: wire convergence data, dense-sparse agreement

Update the requirements doc's Open Questions section to reflect resolved items.

**Test scenarios:**
- Happy path: benchmarks.md contains quantitative validation results with config file names, expected values, and tolerances
- Documentation is sufficient for independent reproduction

**Verification:**
- benchmarks.md updated with all 4 rungs
- Config file names and expected values are documented for each test

---

## System-Wide Impact

- **Interaction graph:** No changes to Fortran executables or shared modules. New Python scripts are standalone consumers of existing output.
- **Error propagation:** Verification scripts exit non-zero on failure, which the shell integration scripts and ctest report correctly.
- **Unchanged invariants:** All existing regression tests continue to pass. The new tests are additive.

---

## Risks & Dependencies

| Risk | Mitigation |
|------|------------|
| Effective mass extraction fails for narrow-gap materials (InSb) due to strong non-parabolicity | Adaptive fitting narrows k-range until R² > 0.9999; if fitting fails, report informational warning rather than test failure |
| 8-band CB effective masses deviate significantly from Vurgaftman (GaAs ~0.046 vs 0.067) | Compare against Kane model prediction m*=Eg/(EP+Eg) at 10% tolerance instead; Vurgaftman comparison is informational only |
| GaSb non-W lacks EV/EC for QW mode | Bulk mode does not use EV/EC (only Eg, DeltaSO on diagonal) — GaSb non-W works correctly for rungs 1-2. Use GaSbW/AlSbW W-variant for all QW configs (rungs 3-4) where EV/EC are needed for band offsets |
| Wire dense-sparse comparison may not apply for large grids (dense solve too slow) | Use small grid (21×21) for dense-sparse comparison at 1e-8 tolerance; convergence test uses larger grids with sparse only |
| FEAST non-determinism at k>1 in wire mode | Only compare k=1 eigenvalues for wire tests (existing pattern) |

---

## Sources & References

- **Origin document:** [docs/brainstorms/2026-05-08-8band-verification-ladder-requirements.md](docs/brainstorms/2026-05-08-8band-verification-ladder-requirements.md)
- **Material parameters:** src/core/parameters.f90
- **Existing benchmarks:** docs/reference/benchmarks.md
- **Input format:** docs/reference/input-reference.md
- **Vurgaftman 2001:** J. Appl. Phys. 89, 5815 (2001) — III-V parameter tables
- **Winkler 2003:** Spin-Orbit Coupling Effects in 2D Electron and Hole Systems
