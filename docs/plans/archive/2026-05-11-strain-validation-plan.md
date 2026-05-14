---
date: 2026-05-11
status: active
origin: docs/brainstorms/2026-05-11-strain-validation-requirements.md
depth: standard
---

# Strain Validation Across Geometries — Implementation Plan

## Overview

End-to-end strain validation suite for InAs/GaAs across bulk, QW, and wire geometries. Validates the full chain (strain computation -> Bir-Pikus -> Hamiltonian -> eigenvalues) against published analytical formulas. Three implementation units produce testable verification scripts; a fourth unit audits the lecture companion; a fifth updates coverage metadata.

The suite addresses the critical gap between existing unit tests (which validate Bir-Pikus formulas in isolation) and end-to-end physics (which no test currently verifies). If a sign error, unit mismatch, or Hamiltonian insertion bug exists in the strain path, every strained calculation is silently wrong.

## Requirements Trace

| Req | Description | Unit |
|-----|-------------|------|
| R1 | Strained bulk InAs/GaAs eigenvalues match analytical Bir-Pikus shifts at k=0, within 1% | U1 |
| R2 | HH-LH splitting matches shear deformation formula `b*(eps_zz - eps_xx)` within 1% | U1 |
| R3 | Strained-minus-unstrained eigenvalue difference equals Bir-Pikus shift exactly | U1 |
| R4a | Strained QW subband energies match Bastard analytical approximations within 3-5% | U2 |
| R5 | HH-LH splitting in strained QW matches shear deformation prediction within 3% | U2 |
| R6 | Strained QW results insensitive to 2x grid density increase | U2 |
| R7 | Wire strain profile matches algebraic biaxial formula at interior core points within 5-10% | U3 |
| R8 | Hydrostatic strain is compressive in InAs core, decays toward zero in GaAs shell | U3 |
| R9 | Band edge shifts in wire core match Bir-Pikus predictions within 20% | U3 |
| R10 | Audit and update lecture_04_strain.py against validated R1-R9 values | U4 |
| R11 | Bulk and QW tests integrated into verification ladder with ctest labels | U1, U2 |
| R12 | Wire tests use dedicated `strain-validation` ctest label | U3 |
| R13 | All tests include COVERAGE annotations | U1, U2, U3 |
| R14 | Strain solver bugs discovered during validation are fixed, not just tested around | U1-U3 |

## Key Technical Decisions

### D1. R4a (analytical) chosen over R4b (nextnano numerical)

No nextnano reference data is readily available for a comparable InAs/GaAs QW config. The Bastard infinite-barrier analytical formulas provide an independent physics check at 3-5% tolerance, which accounts for model limitations (infinite barrier, parabolic bands). This is the more valuable validation because it tests physics correctness, not just numerical agreement between two FD implementations. S6 standard-star benchmark already provides regression-reference validation for strained QW.

### D2. Algebraic biaxial formula as wire strain reference (R7)

The Eshelby square-inclusion analytical solution requires Fourier series evaluation and is complex to implement independently. The interior strain of a finite square inclusion on a Cartesian grid approaches the algebraic biaxial formula at points far from boundaries. Validation compares PDE solver output at interior core points against this formula with 5-10% tolerance. This is simpler, more robust, and directly tests the physics that matters: does the Navier-Cauchy PDE produce the correct lattice-mismatch-driven strain in the core?

### D3. HH-LH ordering under compressive strain (R5 resolution)

For InAs on GaAs (compressive in-plane strain), the Bir-Pikus formulas give:
- eps_xx = eps_yy = -0.0669 (compressive)
- eps_zz = +0.0727 (tensile, Poisson relaxation)
- Q_eps = b_dp/2 * (eps_zz - eps_xx) = -1.8/2 * 0.1396 = -0.1257 eV (negative)
- delta_EHH = -P_eps + Q_eps = -0.0611 + (-0.1257) = -0.1868 eV (HH shifts down)
- delta_ELH = -P_eps - Q_eps = -0.0611 - (-0.1257) = +0.0646 eV (LH shifts up)

HH is **below** LH under compressive strain. The HH-LH splitting is -0.251 eV (HH below LH by 0.251 eV). The test must verify this ordering explicitly.

### D4. LH-SO mixing required for bulk eigenvalue comparison (R1)

The naive diagonal Bir-Pikus shifts (delta_EHH, delta_ELH, delta_ESO) are NOT the actual eigenvalues because QT2 coupling mixes LH and SO bands into a 2x2 sub-problem. The reference eigenvalues must be computed by solving this 2x2 system, following the pattern already established in `scripts/lecture_04_strain.py` `bir_pikus_bulk()`. The coupled eigenvalues differ significantly from the diagonal-only shifts (e.g., LHSO_high = +0.120 eV vs diagonal ELH = +0.065 eV).

### D5. Verification ladder extension: rungs 5 and 6

Bulk strain (R1-R3) becomes verification rung 5, QW strain (R4a, R5, R6) becomes rung 6. This follows the established "simplest to most complex" hierarchy. Wire strain gets a separate `strain-validation` label because it exercises the Navier-Cauchy PDE solver — a distinct subsystem from the Hamiltonian validation that the ladder tracks.

### D6. Minimal inline configs (not committed config files)

Following the lecture_04 pattern, all test configs are written as Python strings and generated at test time. This avoids the parser confusion caused by extra fields in committed configs (e.g., `bulk_gaas_strained.cfg` has gWhichBand, SC, feast fields that are not needed). Each test writes only the fields it needs.

## InAs/GaAs Bir-Pikus Reference Values

Computed from `parameters.f90` InAs entries (Vurgaftman 2001) and GaAs substrate a0 = 5.65325 Angstrom:

```
InAs parameters: a0=6.0583, C11=832.9, C12=452.6, ac=-5.08, av=1.00, b_dp=-1.8, d_dp=-3.6
GaAs substrate:  a_sub=5.65325

eps_xx = eps_yy = (5.65325 - 6.0583) / 6.0583 = -0.06691
eps_zz = -2 * 452.6/832.9 * (-0.06691) = +0.07273
Tr_eps = -0.06109
P_eps = -av * Tr_eps = +0.06109
Q_eps = b_dp/2 * (eps_zz - eps_xx) = -0.9 * 0.13964 = -0.12568

Diagonal Bir-Pikus shifts:
  delta_Ec  = ac * Tr_eps = -5.08 * (-0.06109) = +0.31034
  delta_EHH = -P_eps + Q_eps = -0.18677
  delta_ELH = -P_eps - Q_eps = +0.06459
  delta_ESO = -P_eps = -0.06109

HH-LH splitting = delta_EHH - delta_ELH = -0.25135 eV (= b_dp * (eps_zz - eps_xx))

LH-SO 2x2 coupling: QT2 = 2*Q_eps = -0.25136, coupling = QT2/sqrt(2) = -0.17772
  a = delta_ELH = 0.06459, b = -0.390 + delta_ESO = -0.45109
  root = sqrt((0.25784)^2 + (0.17772)^2) = 0.31316
  E_LHSO_low  = -0.50641 (mixed SO/LH, deeper in VB)
  E_LHSO_high = +0.11991 (mixed LH/SO, above HH)

Strained eigenvalues at k=0 (ascending):
  -0.50641, -0.50641, -0.18677, -0.18677, +0.11991, +0.11991, +0.72734, +0.72734
```

## Implementation Units

### U1: Bulk Strain Verification (R1, R2, R3, R11, R13, R14)

**Script:** `tests/integration/verify_strain_rung5_bulk.py`

**What:** Validates that the bulk strain path (strainSubstrate -> apply_bp_strain_inline) produces correct eigenvalues, HH-LH splitting, and that strain is a pure additive modification.

**Files created:**
- `tests/integration/verify_strain_rung5_bulk.py` — verification script

**Files modified:**
- `tests/CMakeLists.txt` — register as `verification_rung5_strain_bulk` with label `"verification"`

**Config pattern:** Three inline configs generated in Python:
1. Unstrained InAs bulk at k=0: `confinement=0`, `material1: InAs`, `strainSubstrate: 0`, `waveVectorStep: 1`, `numcb: 2`, `numvb: 6`
2. Strained InAs bulk at k=0: same but `strainSubstrate: 5.65325` (GaAs lattice constant)
3. (R2 uses the same strained run, extracting the HH-LH splitting from eigenvalues)

**Test scenarios:**

1. **R1 — Strained eigenvalue match.** Run strained InAs bulk at k=0. Parse eigenvalues from `output/eigenvalues.dat`. Compare all 8 eigenvalues against the LH-SO mixed reference values computed by `bir_pikus_bulk()` (imported or reimplemented). Tolerance: 1% (`TOL_ANALYTICAL`). The script must compute reference values from InAs parameters (read from a shared constants dict matching `parameters.f90`).

2. **R2 — HH-LH splitting.** From the same strained run, extract the HH eigenvalue pair (bands 3-4) and the LHSO_high pair (bands 5-6). Compute HH-LH splitting = E_HH - E_LHSO_high. Compare against `b_dp * (eps_zz - eps_xx) = -0.25135 eV`. Tolerance: 1%.

3. **R3 — Strain as additive modification.** Run both unstrained and strained InAs bulk at k=0. Compute difference eigenvalue-by-eigenvalue (sorted ascending). Compare each difference against the Bir-Pikus diagonal shifts (delta_Ec, delta_EHH, delta_ELH, delta_ESO). Since LH-SO mixing changes the eigenvalue ordering, the comparison must pair eigenvalues by band character (not by sorted position). Simpler approach: compare CB shift (= delta_Ec, no mixing), HH shift (= delta_EHH, no mixing), and the LH+SO block shift as a sum. Tolerance: numerical precision (1e-10) for CB and HH; the LH+SO block requires solving the 2x2 system.

4. **R14 — Bug discovery protocol.** If any test fails, investigate the root cause in the strain solver code path before adjusting tolerances. Document any bugs found.

**COVERAGE annotations:**
```
# COVERAGE: observable=strain_shift geometry=bulk material=InAs ref=Vurgaftman2001
# COVERAGE: observable=HH_LH_splitting geometry=bulk material=InAs/GaAs ref=Vurgaftman2001
```

**Dependencies:** None (first unit).

**Estimated time:** 30-40 seconds per ctest run (two bulk runs, ~1s each).

---

### U2: QW Strain Verification (R4a, R5, R6, R11, R13, R14)

**Script:** `tests/integration/verify_strain_rung6_qw.py`

**What:** Validates that the QW strain path (strain block -> compute_strain_qw -> compute_bir_pikus_blocks -> add_bp_strain_dense) produces correct subband energies and HH-LH splitting.

**Files created:**
- `tests/integration/verify_strain_rung6_qw.py` — verification script

**Files modified:**
- `tests/CMakeLists.txt` — register as `verification_rung6_strain_qw` with label `"verification"`, `TIMEOUT 600`

**Config pattern:** Inline configs following the 2-layer pattern:
1. Strained InAs/GaAs QW: `confinement=1`, `material1: GaAs -60 60 0`, `material2: InAs -10 10 1`, `strain: T`, `strain_ref: GaAs`, `strain_solver: pardiso`, `piezo: F`, `waveVector: k0`, `waveVectorStep: 1`, `FDstep: 201`, `FDorder: 2`, `numcb: 4`, `numvb: 8`
2. (R6) Same config with `FDstep: 401` (2x grid density)

**Key constraint:** Use the 2-layer config pattern (barrier full domain, well overwrites center). The existing `qw_inas_gaas_strained.cfg` already uses this pattern correctly.

**Test scenarios:**

1. **R4a — Strained QW subband energies.** Run strained InAs/GaAs QW at k=0. Parse eigenvalues. Identify ground-state electron (lowest CB eigenvalue), ground-state HH (highest VB eigenvalue), ground-state LH (next-highest VB eigenvalue after HH). Compare against Bastard infinite-barrier analytical formulas:
   - Electron: `E_e1 = Eg + hbar^2 * pi^2 / (2 * m*_e * L^2)` where L = 20 Angstrom (well width), m*_e from InAs parameters
   - HH: `E_hh1 = hbar^2 * pi^2 / (2 * m*_hh * L^2)` where m*_hh from gamma1, gamma2
   - LH: `E_lh1 = hbar^2 * pi^2 / (2 * m*_lh * L^2)` where m*_lh from gamma1, gamma2
   Tolerance: 3-5% (accounts for infinite-barrier approximation vs finite barriers in the 8-band model).

2. **R5 — HH-LH splitting in strained QW.** Extract HH and LH ground-state energies. Compute splitting. Compare against the bulk shear deformation prediction (0.251 eV from D3) modified by confinement effects. The QW confinement plus strain HH-LH splitting should exceed the strain-only splitting because confinement typically increases HH-LH separation. At minimum, verify the splitting is > 0.251 eV and < 0.5 eV. Verify HH is below LH (Q_eps < 0 ordering confirmed in D3).

3. **R6 — Grid convergence.** Run with FDstep=201 and FDstep=401. Compare ground-state eigenvalues between the two runs. All subband energies must agree within 1% (demonstrating numerical stability of strain+FD combination).

4. **R14 — Bug discovery protocol.** Same as U1.

**COVERAGE annotations:**
```
# COVERAGE: observable=strain_shift geometry=QW material=InAs/GaAs ref=Bastard1981
# COVERAGE: observable=HH_LH_splitting geometry=QW material=InAs/GaAs ref=Vurgaftman2001
```

**Dependencies:** U1 (bulk reference values inform QW expectations).

**Estimated time:** ~2 minutes per ctest run (QW solve is slower; R6 doubles it).

**Bastard analytical formulas (reference):**

For an infinite-barrier QW of width L, the nth subband energy for carrier type with effective mass m* is:
```
E_n = n^2 * pi^2 * hbar^2 / (2 * m* * L^2)
```

InAs QW parameters (from parameters.f90):
- L = 20 Angstrom (well width from config: -10 to +10)
- m*_e = 0.026 m0 (from parameters.f90; but 8-band model gives different m*)
- m*_hh (along [001]) = m0 / (gamma1 - 2*gamma2) = 1/(20.0 - 2*8.5) = 1/3.0 = 0.333 m0
- m*_lh (along [001]) = m0 / (gamma1 + 2*gamma2) = 1/(20.0 + 2*8.5) = 1/37.0 = 0.027 m0

Note: The 8-band model m*_e differs from the parameters.f90 tabulated value (29% deviation for GaAs per verification ladder learnings). For the Bastard comparison, use the 2-band Kane m*_e = Eg/(EP + Eg) for better agreement with the 8-band model:
- m*_kane = 0.417 / (21.5 + 0.417) = 0.417/21.917 = 0.01903 m0

This gives a more meaningful comparison at 3-5% tolerance.

---

### U3: Wire Strain Validation (R7, R8, R9, R12, R13, R14)

**Script:** `tests/integration/verify_strain_wire_profile.py`

**What:** Validates that the wire Navier-Cauchy PDE solver produces correct strain profiles and that the resulting Bir-Pikus band edge shifts are physically correct.

**Files created:**
- `tests/integration/verify_strain_wire_profile.py` — validation script

**Files modified:**
- `tests/CMakeLists.txt` — register as `strain_validation_wire` with label `"strain-validation"`, `TIMEOUT 600`

**Config pattern:** Inline config based on existing `wire_inas_gaas_strain.cfg`:
- `confinement=2`, `wire_nx: 30`, `wire_ny: 30`, `wire_dx: 5.0`, `wire_dy: 5.0`
- `wire_shape: rectangle`, `wire_width: 150.0`, `wire_height: 150.0`
- `numRegions: 2`, `region: GaAs 40.0 150.0`, `region: InAs 0.0 40.0`
- `strain: T`, `strain_ref: GaAs`, `strain_solver: pardiso`, `piezo: F`
- `waveVector: kz`, `waveVectorMax: 0.01`, `waveVectorStep: 2`
- `numcb: 4`, `numvb: 8`

**Wire strain output parsing:** The wire solver writes strain data to `output/strain.dat` with columns for position and strain tensor components. The script must parse this file to extract eps_xx, eps_yy, eps_zz at each grid point. If the strain output format is not documented, inspect the actual output or read `outputFunctions.f90` to determine the format.

**Test scenarios:**

1. **R7 — Interior strain profile match.** Parse the wire strain output. At the center of the InAs core (grid point ~15,15 for a 30x30 grid with core 0-40 out of 150 Angstrom), compare eps_xx, eps_yy, eps_zz against the algebraic biaxial formula:
   - eps_xx_ref = eps_yy_ref = (a_GaAs - a_InAs) / a_InAs = -0.0669
   - eps_zz_ref = -2 * C12/C11 * eps_xx_ref = +0.0727
   Tolerance: 10% at center (accounts for finite core size, boundary effects, first-order FD at boundaries). Also check that the strain profile is approximately uniform across the interior core (variance < 5% across the central half of the core region).

2. **R8 — Hydrostatic strain sign and decay.** Compute Tr(eps) = eps_xx + eps_yy + eps_zz at all grid points. Assert:
   - Tr(eps) < 0 in the InAs core (compressive hydrostatic strain)
   - |Tr(eps)| in the GaAs shell at the outer boundary is < |Tr(eps)| at the interface (decay toward zero)
   - The interface region shows a sign change or sharp gradient in strain components

3. **R9 — Band edge shift consistency.** At the InAs core center, compute the expected Bir-Pikus shifts from the strain tensor (using `compute_bp_scalar` formulas). Compare against the eigenvalue shifts relative to unstrained InAs bulk. CB shift should be positive (delta_Ec = +0.31 eV) within 20%. HH-LH splitting magnitude should exceed the minimum threshold of 0.20 eV (from the bulk prediction of 0.25 eV, relaxed to account for finite-size effects).

4. **R14 — Bug discovery protocol.** Same as U1.

**COVERAGE annotations:**
```
# COVERAGE: observable=strain_shift geometry=wire material=InAs/GaAs-core-shell ref=biaxial_analytical
```

**Dependencies:** U1 (bulk reference values establish the analytical baseline).

**Estimated time:** ~2 minutes per ctest run (PARDISO solve for Navier-Cauchy).

**Risk: strain output format.** The wire solver may not output the full strain tensor to a parseable file. If not, the script may need to:
- Parse eigenvalue shifts instead (less direct but still validates the physics chain)
- Or add a strain output path to the Fortran code (minor, controlled change)

This risk should be assessed at implementation time by running the wire config and inspecting `output/` contents.

---

### U4: Lecture Audit (R10)

**Script modification:** `scripts/lecture_04_strain.py`

**What:** Audit the lecture script against the validated R1-R9 numerical values. Update any hardcoded constants, expected values, or assertions that differ.

**Files modified:**
- `scripts/lecture_04_strain.py` — update values and assertions
- `docs/lecture/lecture_04_strain.md` — update companion markdown if it exists

**Approach:**

1. Run `scripts/lecture_04_strain.py` and capture all PASS/FAIL results.
2. Compare Section 2 (strained bulk GaAs) values against the GaAs Bir-Pikus reference (already validated in the existing test). Verify GaAs analytical values match the code's `bir_pikus_bulk()` function.
3. Compare Section 3 (strained InAs/GaAs QW) values against the U2 validated results. If the lecture uses hardcoded expected values, update them to match.
4. If no discrepancies are found, record that fact as a comment in the script.

**Dependencies:** U1, U2 (validated reference values needed for comparison).

---

### U5: Validation Universe and CTest Registration (R11, R12, R13)

**What:** Update the validation universe with new strain cells and register all tests in CMake.

**Files modified:**
- `tests/integration/validation_universe.yml` — add cells for strain observables
- `tests/CMakeLists.txt` — register all test scripts with appropriate labels

**New validation universe cells:**
```yaml
  # --- Strain validation ---
  - observable: strain_shift
    geometry: bulk
    material: InAs/GaAs
    tier: required
    reference: Vurgaftman2001

  - observable: HH_LH_splitting
    geometry: bulk
    material: InAs/GaAs
    tier: required
    reference: Vurgaftman2001

  - observable: strain_shift
    geometry: QW
    material: InAs/GaAs
    tier: required
    reference: Bastard1981

  - observable: strain_shift
    geometry: wire
    material: InAs/GaAs-core-shell
    tier: required
    reference: biaxial_analytical
```

**CTest registration pattern:**
```cmake
# --- Strain Verification (Rungs 5-6) ---
add_test(
    NAME verification_rung5_strain_bulk
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_strain_rung5_bulk.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(
    NAME verification_rung6_strain_qw
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_strain_rung6_qw.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    verification_rung5_strain_bulk verification_rung6_strain_qw
    PROPERTIES LABELS "verification" TIMEOUT 600
)

# --- Wire Strain Validation ---
add_test(
    NAME strain_validation_wire
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/integration/verify_strain_wire_profile.py
        ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set_tests_properties(
    strain_validation_wire
    PROPERTIES LABELS "strain-validation" TIMEOUT 600
)
```

**Dependencies:** U1, U2, U3 (scripts must exist before registration).

---

## Sequencing

```
U1 (bulk) ────────> U2 (QW) ────────> U4 (lecture audit)
      \                                  ^
       ──────> U3 (wire) ───────────────┘
                     \
                      ──> U5 (ctest + universe) ── done
```

U1 must come first because it establishes the analytical reference baseline. U2 and U3 can proceed in parallel after U1. U4 depends on both U1 and U2 (needs validated bulk and QW values). U5 is the final integration step.

## Risk Register

| Risk | Impact | Mitigation |
|------|--------|------------|
| Wire strain output format unknown | R7-R9 blocked | Assess at U3 start by running wire config; eigenvalue-based fallback if needed |
| Bastard analytical tolerance too tight | R4a fails at 3% | Relax to 5% (documented in requirements as acceptable) |
| gfortran -O3 segfault in strained QW path | U2 blocked | Build with `-O2` for development; ensure existing inline-scalar fix is in place |
| LH-SO mixing changes eigenvalue ordering | R3 pairing fails | Use band-character pairing (CB, HH, LH+SO block) not sorted-position pairing |
| Existing lecture_04 values are stale | R10 scope expands | Audit first, then update only where values differ from validated references |

## Shared Infrastructure

All scripts import from `tests/integration/star_helpers.py`:
- `run_exe(build_dir, config_content, work_dir)` — write config, run bandStructure, return output dir
- `parse_eigenvalues(filepath)` — parse eigenvalues.dat
- `compare_value(actual, expected, tolerance, label)` — assert with formatted message
- `TOL_EXACT`, `TOL_ANALYTICAL`, `TOL_NUMERICAL` — tolerance tiers

Each script follows the verification script pattern:
- Signature: `script.py <build_dir> <source_dir>`
- Uses `tempfile.mkdtemp()` for working directory
- Returns 0 on all-pass, 1 on any failure
- Prints per-section PASS/FAIL with actual vs expected values
