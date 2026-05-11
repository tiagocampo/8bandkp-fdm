---
title: "Lecture-test pair robustness: physics correctness, performance, false-positive fixes, and FEAST fallback"
date: "2026-05-10"
last_updated: "2026-05-11"
category: test-failures
module: lecture-test-pairs
problem_type: test_failure
component: testing_framework
severity: high
symptoms:
  - "L04 summary hardcodes PASS regardless of actual validation results"
  - "L05 Landau test uses InAs (55% tolerance) instead of GaAs, inflating tolerance to mask non-parabolicity"
  - "L06 ISBT test asserts alpha > 0 on positive-definite quantity (always true)"
  - "L08 wire test uses 31x31 grid taking >10 min for a subband-count validation"
  - "L13 BdG FEAST returns 0 eigenvalues with no fallback; QW config has unrealistic pairing parameters"
root_cause: missing_validation
resolution_type: test_fix
tags:
  - lecture-test-pairs
  - false-positive
  - assertion-quality
  - subprocess-timeout
  - constant-deduplication
  - feast-eigensolver
  - landau-levels
  - bdg-hamiltonian
  - grid-scaling
  - physics-correctness
related_components:
  - documentation
  - tooling
---

# Lecture-test pair robustness: physics correctness, performance, false-positive fixes, and FEAST fallback

## Problem

14 Python lecture-test pair scripts (`scripts/lecture_00_quickstart.py` through `scripts/lecture_13_topological.py`) serve dual roles as pedagogical companions and integration tests. Five commits (`5a0984f..6923863`) addressed multiple classes of issues: incorrect material choices for physics validation (L05), oversized grids causing excessive runtime (L08), FEAST eigensolver failures at finite magnetic field (L13, Fortran fix), false-positive risks from hardcoded PASS and vacuous assertions (L04, L06), duplicate constants and missing timeouts, plus documentation mismatches across all 14 lecture chapters. Since these scripts validate the Fortran solver, false positives or fragile tests silently hide real regressions.

## Symptoms

- L04 `main()` always prints `[PASS]` for every section -- functions return data but not pass/fail status
- L05 Landau test uses InAs with 55% tolerance; strong non-parabolicity masks failure to validate against Kane effective mass
- L06 ISBT peak test asserts `alpha > 0` which is always true for a computed absorption spectrum
- L08 wire test uses 31x31 grid (7688-dim Hamiltonian, >10 min) for a subband-count check
- L13 BdG Majorana section uses QW config with unrealistic parameters; FEAST returns 0 eigenvalues with no fallback
- L13 subprocess timeout crashes script with uncaught `TimeoutExpired`
- L13 Z2 invariant tested at exact critical width where noise flips result
- Rashba shell test runs `topologicalAnalysis` with no timeout guard
- Physical constants duplicated across 4+ scripts with inconsistent precision
- All 14 lecture docs reference wrong script names and figure filenames

## What Didn't Work

- **Wrong material for Landau validation (L05)**: InAs has EP=21.5 eV, Eg=0.417 eV -- strong 8-band mixing at B=5T causes Landau spacings to deviate up to 55% from simple cyclotron formula. The 55% tolerance was embarrassingly loose and masked the inability to validate against a clean analytical formula. (session history)
- **Oversized grid for simple validation (L08)**: 31x31 grid (7688-dim) took >10 min but the test only validates eigenvalue existence and subband count, not spatial resolution. The strict assertion `n_evals == expected_nev` also broke because FEAST with `feast_m0=-1` auto-detects and may find more eigenvalues than requested. (session history)
- **FEAST eigenvalue window failure in BdG mode (L13)**: The Rashba wire BdG config produced zero eigenvalues at high B-fields. The initial approach switched L13 to use the QW BdG config entirely, but the proper fix was adding a Gershgorin-based auto-window fallback in the Fortran code. (session history)
- **Hardcoded PASS (L04)**: Section functions returned data arrays but no pass/fail boolean. `main()` printed PASS without checking results.
- **Vacuous assertion (L06)**: `alpha > 0` is always true for matrix-element-squared quantities. The assertion never catches a wrong peak position or zero signal.

## Solution

Thirteen fixes applied across 5 commits covering 25 files:

### Core limitation fixes (commit `5a0984f`)

#### 1. Match material to the validation formula (L05 Landau)

Switched from InAs to GaAs (weak non-parabolicity). Used Kane effective mass `m* = Eg/(EP+Eg)` instead of bare `m*_e`. Added spin-splitting filter: keep only spacings > 50% of `hbar*omega_c` to separate inter-Landau-level spacing (~12 meV) from spin-orbit sub-splitting (~1 meV). Tolerance tightened from 55% to 15%.

```python
# Before: InAs with bare mass and loose tolerance
m_star = MATERIALS['InAs']['meff']  # 0.026 m0
mean_spacing = np.mean(spacings[:5])
TOL_LANDAU = 0.55  # 55%

# After: GaAs with Kane mass and spin filtering
p = MATERIALS['GaAs']
m_star_kane = p['Eg'] / (p['EP'] + p['Eg'])  # 0.0501 m0
threshold = 0.5 * hbar_omega_c_meV
landau_spacings = spacings_meV[spacings_meV > threshold]
TOL_LANDAU = 0.15  # 15%
```

#### 2. Scale grid to the physics being tested (L08 wire)

Reduced from 31x31 (7688-dim, >10 min) to 21x21 (3528-dim, ~3 min). Fixed FEAST assertion: `n_evals >= min_expected_nev` instead of strict equality.

```python
# Before
cfg = CONFIGS_DIR / "wire_gaas_31x31.cfg"
assert n_evals_first == expected_nev  # strict equality

# After
cfg = CONFIGS_DIR / "wire_gaas_rectangle.cfg"
assert n_evals_first >= min_expected_nev  # FEAST may find more
```

#### 3. Add FEAST auto-window fallback (L13, Fortran)

Added Gershgorin-based auto-window in `main_topology.f90` when initial FEAST call finds 0 eigenvalues. Switched L13 to InAs Rashba wire (mu=0.1 meV, Delta=0.1 meV, B_crit ~ 1.22 T). Reduced grid to 11x11 for B-sweep (~30 s per point). Reduced sweep to 0-3 T.

```fortran
! main_topology.f90: auto-window retry on FEAST miss
if (eigen_res_local%nev_found == 0) then
  print *, 'Warning: FEAST found no eigenvalues in the search window'
  print *, '  Retrying with auto-computed energy window...'
  call auto_compute_energy_window(H_bdg_csr, eigen_cfg_local%emin, eigen_cfg_local%emax)
  call eigen_solver_local%solve(H_bdg_csr, eigen_cfg_local, eigen_res_local)
end if
```

### Code review fixes (commit `c48e533`)

#### 4. Return pass/fail tuples from section functions (L04)

Each section function now returns `(all_pass, data)`. `main()` collects booleans and prints actual PASS/FAIL, exiting with `sys.exit(1)` on failure.

```python
# Before: returns data only
def test_unstrained_reference():
    ...
    return evals

# After: returns status + data
def test_unstrained_reference():
    ...
    return all_pass, evals
```

#### 5. Replace vacuous assertions with meaningful bounds (L06)

```python
isbt_e_min, isbt_e_max = 0.03, 0.30  # eV, generous ISBT range
min_alpha = 1.0  # cm^-1, above numerical noise floor
in_range = isbt_e_min <= peak_energy <= isbt_e_max
above_noise = peak_alpha > min_alpha
passed = in_range and above_noise
```

#### 6. Catch subprocess.TimeoutExpired (L13)

```python
try:
    result = subprocess.run([str(exe)], cwd=work_dir,
                            capture_output=True, text=True, timeout=120)
except subprocess.TimeoutExpired:
    print(f"    TIMEOUT after 120s for width={width_angstrom}A")
    return None
```

#### 7. Add timeout to shell integration tests

```bash
timeout 300 "$EXE" > test_output.log 2>&1 || RC=$?
if [ $RC -eq 124 ]; then
    echo "FAIL: topologicalAnalysis timed out after 300s"
    exit 1
fi
```

#### 8. Centralize physical constants (star_helpers)

```python
# star_helpers.py (single source of truth)
HBAR_J_S  = 1.054571817e-34   # hbar in J*s
E_CHARGE  = 1.602176634e-19   # elementary charge in C
M0_KG     = 9.1093837015e-31  # free electron mass in kg

# lecture_05_gfactor.py, lecture_12_extending.py
from star_helpers import HBAR_J_S, E_CHARGE, M0_KG
```

#### 9. Skip phase-boundary critical point (L13 Z2)

```python
# Before: tests at exactly 70A (phase boundary)
expected_z2 = 0 if w < 70 else 1

# After: skip transition region
if w <= 60: expected_z2 = 0   # firmly trivial
elif w >= 80: expected_z2 = 1 # firmly topological
else: continue                 # skip transition region near ~70A
```

#### 10. Document physics rationale for numerical thresholds (L05)

```python
# The 0.5 threshold cleanly separates the two regimes for GaAs at B=5T
# (spin-split ~1 meV << threshold ~6 meV << Landau spacing ~12 meV).
threshold = 0.5 * hbar_omega_c_meV
```

### Documentation fixes (commits `7b80237`, `6a16c88`)

#### 11. Correct script names and figure references across all 14 lectures

All 14 `docs/lecture/*.md` files referenced wrong script names (e.g., `lecture_01_bulk_band_structure.py` instead of `lecture_01_bulk.py`) and wrong figure filenames. Batch find-and-replace across all docs. Three docs (05, 08, 13) received additional updates for physics changes.

### Minor fixes (commit `6923863`)

#### 12. Remove duplicate `if __name__` guard (L04)

`lecture_04_strain.py` had the guard duplicated after the pass/fail propagation fix.

## Why This Works

The unifying root causes are **missing validation** and **physics-test mismatch**:

- **Material must match the validation formula**: Testing a parabolic-mass formula with a highly non-parabolic material (InAs) inflates tolerance to the point of meaninglessness. GaAs with Kane effective mass gives tight, meaningful validation.
- **Grid size is a test parameter, not a physics requirement**: Use the smallest grid that validates the property under test. Subband count doesn't need 31x31 resolution.
- **FEAST eigenvalue searches need graceful fallback**: Search-window misses are expected in parameter sweeps; Gershgorin estimation provides reliable automatic recovery.
- **Test functions must return pass/fail, not just data**: A script that always prints PASS regardless of results is worse than no test at all.
- **Assertions on physical quantities must check meaningful ranges**: Positive-definite quantities are always > 0. Check energy range AND noise floor.
- **Physical constants belong in one module**: When `HBAR_J_S` appears in three files, one will eventually drift.

(session history: the initial implementation session created all 14 scripts via parallel subagent dispatch, with 3 of 14 requiring retries. The Landau tolerance was first widened as a workaround, then properly fixed by switching material. The BdG FEAST window issue was a pre-existing limitation in the Fortran code that was first worked around by switching configs, then properly fixed with auto-window fallback.)

## Prevention

- **Section functions must return `(bool, data)` tuples.** The `main()` summary must use the returned booleans, not hardcode PASS.
- **Every `subprocess.run` with `timeout` must catch `TimeoutExpired`.** The timeout parameter raises an exception, it does not return an error code.
- **Assertions on physical quantities must check meaningful ranges, not sign alone.** Positive-definite quantities are always > 0. Check energy range AND noise floor.
- **Never test at exactly the critical point of a phase transition.** Test firmly on each side, skipping the transition region.
- **Physical constants belong in one module.** Import from `star_helpers` -- never define `HBAR`, `M0`, `E_CHARGE` locally.
- **Match the material to the validation formula.** If tolerance exceeds 20%, suspect material mismatch or wrong effective mass.
- **Scale grids to the physics being tested.** Use the smallest grid that validates the property. Subband count doesn't need high resolution.
- **FEAST sweeps need auto-window fallback.** Add Gershgorin-based window recomputation when initial search finds zero eigenvalues.
- **Keep script names and figure references synchronized.** When renaming scripts or generating new figures, batch-update all lecture docs.

## Related Issues

- `docs/solutions/test-failures/csr-test-infra-stop1-oob-tautology-fixes-2026-05-09.md` -- HIGH overlap. Same class of test quality problems (vacuous assertions, hardcoded values, tautological checks) in Fortran pFUnit test infrastructure.
- `docs/solutions/logic-errors/standard-star-benchmark-suite-logic-errors-2026-05-09.md` -- MODERATE overlap. Established the project pattern for subprocess timeout handling in Python test scripts.
- `docs/solutions/logic-errors/topological-magnetic-index-logic-errors-2026-05-08.md` -- LOW-MODERATE. BdG/Majorana path in `main_topology.f90` is the same subsystem where the auto-window fallback was added.
- GitHub Issue #13 (OPEN): Parent feature issue for BdG/topological work stream.
- Plan: `docs/plans/2026-05-10-004-feat-executable-lecture-test-pairs-plan.md` -- implementation plan for all 14 scripts.
