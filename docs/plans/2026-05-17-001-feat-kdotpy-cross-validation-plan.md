---
title: kdotpy Cross-Code Validation Pipeline
type: feat
status: active
date: 2026-05-17
origin: docs/brainstorms/kdotpy-cross-validation-requirements.md
---

# kdotpy Cross-Code Validation Pipeline

## Summary

Build an automated cross-code validation pipeline comparing our Fortran 8-band k.p code against kdotpy (Wuerzburg, Python, plane-wave discretization). The pipeline uses kdotpy as a Python library (not CLI) for clean programmatic comparison. Parameter mapping resolves the convention differences: F=(1/meff-1)/2 for the CB diagonal, kappa=0 to eliminate interface-only discrepancies, and explicit strain deformation potential conversions. Verification proceeds in layers from bulk (exact agreement) through QW (convergence studies) to extended physics (wire, Landau, g-factor, strain, Berry, SC), with all discrepancies investigated strictly.

---

## Problem Frame

Our 91 tests validate against analytical formulas and literature, but no independent k.p code provides ground truth. kdotpy v1.3.1 is now installed locally with a programmable Python API, making systematic cross-code comparison feasible for the first time. The critical blocker — the F/kappa parameter convention gap — has been resolved through source-code analysis: our `A = 1/meff` maps to kdotpy's `2F+1`, and kappa is irrelevant for bulk and can be zeroed for heterostructures. (See origin: `docs/brainstorms/kdotpy-cross-validation-requirements.md`)

---

## Requirements

- R1. kdotpy venv installed (DONE). R2. validation/ directory structure. R3. Parameter mapping module with F/kappa/ge derivation. R4. Comparison runner.
- R5. III-V materials ported to kdotpy INI. R6. II-VI materials ported to parameters.f90.
- R7-R9. Bulk verification (k=0, dispersion, Zeeman).
- R10-R12. QW verification (subbands, dispersion, convergence).
- R13-R14. Wire verification. R15-R16. Landau levels.
- R17-R18. BHZ/g-factor. R19-R20. Strain.
- R21-R22. Berry curvature + LL Chern. R23. Self-consistent.
- R24-R26. Discrepancy logging, resolution protocol, layered verification order.

**Origin actors:** None (infrastructure task)
**Origin flows:** None (layered verification order defined in R26)
**Origin acceptance examples:** None

---

## Scope Boundaries

- Physics only our code has: BdG/Majorana, Z2 (Fu-Kane), optical spectra, exciton, scattering, spin-resolved spectra, LDOS via Green's functions, topological phase diagrams, spectral functions
- Physics only kdotpy has: 6-band mode, BIA/Dresselhaus, exchange (Mn), cubic Zeeman, arbitrary crystal orientation, lattice regularization, symbolic Hamiltonian, GPU solvers
- Code quality, style, or performance comparison
- Publishing a joint benchmark paper (future work)
- Parameter optimization or fitting to experimental data

### Deferred to Follow-Up Work

- Symbolic Hamiltonian matrix-element comparison (kdotpy-only capability, could provide deeper validation)
- BIA/Dresselhaus terms (kdotpy-only, no comparator in our code)
- Optical transition comparison in LL mode (different from our optical module)

---

## Context & Research

### Relevant Code and Patterns

- `tests/integration/star_helpers.py` — reusable runners (`run_executable`, `run_exe`), parsers (`parse_eigenvalues`, `parse_gfactor`), physics utilities (`extract_effective_mass`, `roth_gfactor`, `bir_pikus_biaxial_001`), comparison (`compare_value`, `format_benchmark_row`), tolerance tiers (`TOL_EXACT=1e-12`, `TOL_ANALYTICAL=0.01`, `TOL_NUMERICAL=0.05`)
- `tests/integration/verify_star_*.py` — standard-star benchmark pattern: run exe in tempdir, parse output, compare against reference, print markdown table
- `tests/regression/configs/` — existing config format (`.cfg`, label-value with material names)
- `src/io/outputFunctions.f90:147-189` — eigenvalue output format: `|k| eval_1 eval_2 ... eval_N` (g14.6, ascending order)
- `src/core/parameters.f90:741` — `A = 1/meff` definition
- `src/physics/hamiltonianConstructor.f90:729` — bulk CB diagonal: `HT(7,7) = A * k2`

### kdotpy API (installed at `validation/kdotpy_env/`)

- **Library usage**: `from kdotpy.diagonalization.diagonalization import hbulk, hz` — returns `DiagDataPoint` with `.eival` (numpy array, meV)
- **Materials**: `from kdotpy.materials import Material, allMaterials` — custom materials via `Material("name", param={...})` or INI files in `~/.kdotpy/materials/`
- **Parameters**: `from kdotpy.physparams import PhysParams` — `PhysParams(kdim=3, m_layers=[mat], l_layers=[1.0], norbitals=8, zres=0.25)`
- **k-points**: `from kdotpy.vector import Vector` — `Vector(kx, ky, kz)` in 1/nm
- **Bulk**: `hbulk(k, params)` → eigenvalues in meV
- **QW**: `hz(k, b, params)` via sparse solver from `hamiltonian/full.py`

### Parameter Mapping (Resolved)

| Our Parameter | kdotpy Parameter | Mapping | Exact? |
|---|---|---|---|
| `A = 1/meff` | `F` | `F = (1/meff - 1) / 2` | Yes (CB diagonal equality) |
| (absent) | `kappa` | Set `kappa = 0` in all INI files | Eliminates interface-only discrepancy |
| `EP` (eV) | `P` (meV·nm) | `P = sqrt(EP * 1000 * hbarm0_kdotpy)` where `hbarm0_kdotpy ≈ 38.1` meV·nm² | Yes (unit conversion) |
| `gamma1/2/3` | `gamma1/2/3` | Direct copy (dimensionless) | Yes |
| `Eg` (eV) | `Ec - Ev` (meV) | `Ec = Eg * 1000`, `Ev = 0` (our EV convention) | Needs careful alignment with our EC/EV |
| `delta_SO` (eV) | `delta_so` (meV) | `× 1000` | Yes |
| `ac` (eV) | `strain_C1` (meV) | `ac * 1000` | Yes |
| `av` (eV) | `strain_Dd` (meV) | `av * 1000` | Yes |
| `b_dp` (eV) | `strain_Du` (meV) | `-1.5 * b_dp * 1000` | Yes (convention) |
| `d_dp` (eV) | `strain_Duprime` (meV) | `-0.5 * sqrt(3) * d_dp * 1000` | Yes (convention) |
| `ge` | `ge` | `2.0` (default, adjust if Zeeman fails) | Approximate |
| (absent) | `q` | `0.0` | N/A (BIA parameter, not used) |

### Institutional Learnings

- CB effective mass from 8-band model ≠ Vurgaftman (29% deviation is expected model physics, not a bug) — `docs/ideation/2026-05-08-core-kp-validation-ideation.md`
- Last-layer-wins overlapping material layer pitfall — `CLAUDE.md`
- Richardson convergence: max-rate strategy, order-adaptive tolerance — `docs/solutions/`
- Our basis ordering: bands 1-4 = VB (HH,LH,LH,HH), 5-6 = SO, 7-8 = CB. kdotpy: 0-1 = CB, 2-5 = VB, 6-7 = SO — eigenvalue-only comparison is permutation-invariant

---

## Key Technical Decisions

- **Use kdotpy as Python library, not CLI**: Cleaner programmatic access via `hbulk(k, params)` / `hz(k, b, params)` returning numpy arrays directly. Avoids CLI parsing, temp files, and subprocess overhead. (Rationale: CLI would require parsing CSV/XML output; library gives structured data.)
- **Set kappa=0 in all kdotpy INI files**: Our code has no kappa parameter. kdotpy's kappa only contributes via `dzkappa` (spatial derivative at interfaces). For bulk, kappa is irrelevant. For QW, setting kappa=0 makes the VB off-diagonal blocks identical between codes. (Rationale: introduces zero error for bulk, eliminates interface discrepancy for QW.)
- **F = (1/meff - 1)/2 for all materials**: Exact mapping that guarantees identical CB diagonal terms. For GaAs: F=6.963. (Rationale: `A * k² = (1/meff) * k²` in our code vs `hbarm0 * (2F+1) * k²` in kdotpy, so `2F+1 = 1/meff`.)
- **Ev = 0, Ec = Eg*1000 for kdotpy**: Our code uses EC = EV + Eg with EV as the energy reference. kdotpy uses Ec and Ev independently. Setting Ev=0 and Ec=Eg*1000 meV aligns the band gap. (Rationale: energy reference is arbitrary; what matters is the gap.)
- **Reuse star_helpers.py pattern**: Tolerance tiers, benchmark table format, executable runners, effective mass extraction. (Rationale: consistency with existing test infrastructure; aggregate_star_benchmarks.py can include cross-code results.)
- **Our code runs as subprocess, kdotpy as library**: Our Fortran code requires compilation and has no Python bindings. kdotpy is Python and can be called directly. (Rationale: avoids coupling build systems.)

---

## Open Questions

### Resolved During Planning

- **F mapping**: `F = (1/meff - 1)/2` is exact for CB diagonal equality. Verified against source: our `HT(7,7) = A*k²` vs kdotpy `t0 = Ec + hbarm0*(2F+1)*k²`.
- **kappa**: Our code has no kappa. Set kappa=0 in kdotpy. Zero effect on bulk, eliminates interface-only discrepancy for QW.
- **kdotpy invocation**: Use as library via `hbulk(k, params)` and `hz(k, b, params)`, not CLI.
- **Strain_C1 convention**: kdotpy's `strain_C1` = our `ac` (same sign, unit conversion eV→meV). Verified against kdotpy blocks.py `hstrain`: `tt = cs * tr_e` matches our `delta_Ec = ac * Tr(eps)`.

### Deferred to Implementation

- **EC/EV alignment**: Our code uses absolute EC/EV values per material (from Vurgaftman). kdotpy's default materials use Ec/Ev with Ev as reference. The param_mapper must handle this correctly — verify by comparing bulk k=0 eigenvalues.
- **P conversion precision**: `P = sqrt(EP * 1000 * hbarm0_kdotpy)` — must use kdotpy's exact hbarm0 value (not ours) to avoid floating-point mismatch at the meV level.
- **QW grid resolution matching**: Our FDstep vs kdotpy's zres are fundamentally different convergence parameters. Richardson extrapolation (R12) resolves this but needs concrete grid spacings determined during implementation.
- **BHZ vs commutator velocity**: kdotpy's BHZ uses numerical dH/dk; our g-factor uses -i[r,H]. May differ by O((P/Eg)^4) for narrow QWs. Tolerance adjusted to 10% for narrow wells.

---

## Output Structure

```
validation/
├── kdotpy_env/                          # Python venv (DONE)
├── shared/
│   ├── param_mapper.py                  # Parameter mapping module
│   ├── kdotpy_runner.py                 # kdotpy library wrapper
│   ├── fortran_runner.py                # Fortran exe subprocess runner
│   ├── comparison.py                    # Comparison + reporting
│   └── materials/                       # kdotpy INI material files
│       ├── gaas.ini
│       ├── inas.ini
│       ├── insb.ini
│       ├── alas.ini
│       ├── gasb.ini
│       ├── alsb.ini
│       ├── inp.ini
│       ├── al20ga80as.ini
│       ├── al30ga70as.ini
│       └── inas_sb_alloys.ini
├── bulk/
│   ├── test_bulk_k0.py
│   ├── test_bulk_dispersion.py
│   ├── test_bulk_zeeman.py
│   └── results/                         # JSON + markdown reports
├── qw/
│   ├── test_qw_subbands.py
│   ├── test_qw_dispersion.py
│   ├── test_qw_convergence.py
│   └── results/
├── wire/
│   ├── test_wire_subbands.py
│   └── results/
├── landau/
│   ├── test_landau_qw.py
│   ├── test_landau_bulk.py
│   └── results/
├── gfactor/
│   ├── test_gfactor_qw.py
│   ├── test_bhz_extraction.py
│   └── results/
├── strain/
│   ├── test_strain_bandedge.py
│   ├── test_strain_qw.py
│   └── results/
├── berry/
│   ├── test_berry_curvature.py
│   ├── test_chern_ll.py
│   └── results/
├── selfconsistent/
│   ├── test_sc_qw.py
│   └── results/
└── run_all.py                           # Full pipeline runner
```

---

## Implementation Units

- U1. **Infrastructure Foundation**

**Goal:** Set up validation directory structure, kdotpy smoke test, and shared runner modules.

**Requirements:** R1, R2, R4

**Dependencies:** None

**Files:**
- Create: `validation/shared/__init__.py`
- Create: `validation/shared/kdotpy_runner.py`
- Create: `validation/shared/fortran_runner.py`
- Create: `validation/shared/comparison.py`
- Create: `validation/shared/materials/` directory
- Create: `validation/bulk/`, `validation/qw/`, `validation/wire/`, `validation/landau/`, `validation/gfactor/`, `validation/strain/`, `validation/berry/`, `validation/selfconsistent/` directories with `results/` subdirs

**Approach:**
- `kdotpy_runner.py`: Wraps kdotpy library calls. `init_kdotpy()` calls `initialize_config()` + `initialize_materials()`. `run_bulk(mat, k_points)` creates PhysParams(kdim=3), calls `hbulk()` for each k-point, returns eigenvalues in meV. `run_qw(mats, layers, k_points, zres)` creates PhysParams(kdim=2), calls sparse diagonalization, returns eigenvalues.
- `fortran_runner.py`: Wraps our Fortran executables via subprocess. Reuses `star_helpers.run_exe()` pattern. Parses `output/eigenvalues.dat` with `star_helpers.parse_eigenvalues()`. Returns eigenvalues in eV (our native unit).
- `comparison.py`: Core comparison logic. `compare_eigenvalues(ours_eV, kdotpy_meV, tolerance_meV)` converts units, aligns band indices (permutation from different basis ordering), computes per-band deltas. `format_report(results)` produces JSON + markdown output following `star_helpers.format_benchmark_row()` pattern.

**Patterns to follow:**
- `tests/integration/star_helpers.py` for executable running, output parsing, and benchmark table formatting
- Tolerance tiers: `TOL_EXACT = 0.01 meV`, `TOL_ANALYTICAL = 0.1 meV`, `TOL_NUMERICAL = 1.0 meV`

**Test scenarios:**
- Happy path: kdotpy init succeeds, HgTe bulk at k=0 returns 8 eigenvalues, Fortran bulk GaAs at k=0 returns 8 eigenvalues
- Edge case: kdotpy materials directory doesn't exist yet (create it)
- Integration: full init → material creation → bulk calculation → eigenvalue comparison pipeline runs without error

**Verification:**
- `python -c "from validation.shared.kdotpy_runner import init_kdotpy; init_kdotpy()"` succeeds
- HgTe bulk k=0 eigenvalues match kdotpy's own output within 0.001 meV

---

- U2. **Parameter Mapping Module**

**Goal:** Build and validate the parameter mapping between our format and kdotpy's format.

**Requirements:** R3

**Dependencies:** U1

**Files:**
- Create: `validation/shared/param_mapper.py`
- Test: `validation/shared/test_param_mapper.py`

**Approach:**
- `FortranToKdotpyMapper` class with methods:
  - `map_material(mat_name)` — reads our `parameters.f90` data (via a hardcoded dictionary extracted from the source), returns kdotpy `Material` object with all parameters converted
  - `write_ini(mat_name, filepath)` — writes kdotpy INI file
  - `get_physparams(mat_name, kdim, **kwargs)` — creates `PhysParams` from mapped material
- Parameter conversion formulas (from resolved research):
  - `F = (1/meff - 1) / 2`
  - `kappa = 0.0` (eliminates interface discrepancy)
  - `P = sqrt(EP_eV * 1000 * HBARM0_KDOTPY)` where `HBARM0_KDOTPY` is read from kdotpy's `physconst.py`
  - `Ec = (EC_paramStruct - EV_paramStruct + Eg) * 1000` in meV (aligns band gap)
  - `Ev = 0.0` (energy reference)
  - `delta_so = delta_SO * 1000` (eV→meV)
  - `strain_C1 = ac * 1000`, `strain_Dd = av * 1000`
  - `strain_Du = -1.5 * b_dp * 1000`, `strain_Duprime = -0.5 * sqrt(3) * d_dp * 1000`
  - `a = a0 / 10` (Å→nm)
  - `ge = 2.0` (default, adjustable)
  - `q = 0.0` (no BIA)
- Hardcode our material parameters from `parameters.f90` into a Python dict (single source of truth, avoids Fortran parsing)
- Validation: create material, build bulk Hamiltonian at k=0, compare 8 eigenvalues against our code's output

**Execution note:** Test-first — write param_mapper tests before implementation. Start with GaAs (best-known material), validate mapping at k=0, then expand.

**Test scenarios:**
- Happy path: GaAs material maps correctly (F=6.963, kappa=0, P≈9.4 meV·nm, Ec≈1519 meV)
- Happy path: Bulk k=0 eigenvalues from mapped GaAs match our code's k=0 eigenvalues within 0.01 meV
- Edge case: InSb (small gap, large EP) — verify no numerical overflow in P conversion
- Error path: unknown material name raises clear error
- Integration: mapped material creates valid PhysParams, hbulk runs without error

**Verification:**
- All 8 GaAs bulk eigenvalues at k=0 agree between mapped kdotpy material and our Fortran code to < 0.01 meV
- Same for InAs, InSb (at least 3 materials validated)

---

- U3. **III-V Material Files for kdotpy**

**Goal:** Create kdotpy INI material files for all our standard-star materials and validate each one.

**Requirements:** R5

**Dependencies:** U2

**Files:**
- Create: `validation/shared/materials/gaas.ini`
- Create: `validation/shared/materials/inas.ini`
- Create: `validation/shared/materials/insb.ini`
- Create: `validation/shared/materials/alas.ini`
- Create: `validation/shared/materials/gasb.ini`
- Create: `validation/shared/materials/alsb.ini`
- Create: `validation/shared/materials/inp.ini`
- Create: `validation/shared/materials/al20ga80as.ini`
- Create: `validation/shared/materials/al30ga70as.ini`
- Create: `validation/shared/materials/inas_sb_alloys.ini`
- Create: `validation/shared/load_materials.py`

**Approach:**
- Use param_mapper to generate INI files programmatically from our parameter data
- `load_materials.py`: copies INI files to `~/.kdotpy/materials/` and calls `allMaterials.load_from_file()`
- For alloy materials (AlGaAs, InAsSb): use kdotpy's `linearmix` feature in INI format or `Material` with composition parameter
- Each INI file includes all required kdotpy fields per R5
- Validation script: for each material, run bulk k=0 in both codes, compare band gap

**Test scenarios:**
- Happy path: each INI file loads without error in kdotpy
- Happy path: bulk k=0 band gap matches our code to < 0.01 meV for each material
- Edge case: alloy materials (AlGaAs) interpolate correctly
- Integration: all 10+ INI files loaded simultaneously, bulk calculation succeeds for each

**Verification:**
- `python -c "from validation.shared.load_materials import load_all; load_all(); from kdotpy.materials import allMaterials; assert 'GaAs' in allMaterials"` succeeds
- Bulk band gaps agree for all materials

---

- U4. **II-VI Material Entries for Our Code**

**Goal:** Add HgTe, CdTe, HgCdTe, CdZnTe to our `parameters.f90`.

**Requirements:** R6

**Dependencies:** U1

**Files:**
- Modify: `src/core/parameters.f90`
- Test: `tests/unit/parameters.pf` (existing unit test)

**Approach:**
- Add material entries following existing pattern (case in `setParams` subroutine)
- Parameters sourced from kdotpy's default material database, verified against primary references:
  - HgTe: Pfeuffer-Jeschke PhD thesis (kdotpy references this)
  - CdTe: same source
  - HgCdTe: Becker et al. PRB 62, 10353 (kdotpy's Eg reference)
  - CdZnTe: lattice constant only well-characterized
- Unit conversion: kdotpy meV → our eV (/1000), nm → Å (*10), meV·nm → eV·Å (/10000)
- For HgCdTe: implement as Vegard interpolation similar to AlGaAs pattern
- Requires CLAUDE.md boundary approval — must cite published references

**Execution note:** Requires approval per CLAUDE.md boundary on parameter modifications. Present proposed parameters with citations before implementation.

**Test scenarios:**
- Happy path: new materials recognized in input.cfg, bulk calculation runs
- Happy path: HgTe bulk k=0 band gap is negative (semimetal) as expected
- Happy path: CdTe bulk k=0 band gap matches literature (~1.6 eV)
- Integration: HgCdTe x=0.3 band gap interpolates correctly between HgTe and CdTe

**Verification:**
- `bandStructure` runs with `material1: HgTe` and produces eigenvalues
- Band gap for CdTe within 5% of literature value (1.606 eV from kdotpy's reference)
- Existing tests still pass (no regression)

---

- U5. **Bulk Verification Suite**

**Goal:** Implement and run bulk comparisons (k=0, dispersion, Zeeman) — the gate for all downstream work.

**Requirements:** R7, R8, R9

**Dependencies:** U2, U3

**Files:**
- Create: `validation/bulk/test_bulk_k0.py`
- Create: `validation/bulk/test_bulk_dispersion.py`
- Create: `validation/bulk/test_bulk_zeeman.py`

**Approach:**
- `test_bulk_k0.py`: For each III-V material, compare all 8 eigenvalues at Gamma. Reports per-band delta. Tolerance: 0.01 meV. This validates the entire parameter mapping pipeline.
- `test_bulk_dispersion.py`: For GaAs, InAs, InSb — sweep k along [100], [110], [111] (20+ points, |k| < 0.1 1/nm). Compare eigenvalues at each k-point (tolerance 0.1 meV). Extract effective masses via parabolic fit, compare to < 1%. This validates F/kappa mapping.
- `test_bulk_zeeman.py`: For GaAs with B along x, y, z (0.1 T, 1 T, 5 T). Compare Zeeman-split eigenvalues (tolerance 0.01 meV). Validates ge/kappa Zeeman mapping.
- All tests output markdown benchmark tables compatible with `aggregate_star_benchmarks.py`

**Execution note:** test_bulk_k0.py is the gate — it must pass before any other comparison proceeds. If it fails, the parameter mapping is wrong.

**Test scenarios:**
- Happy path: GaAs k=0 eigenvalues agree to < 0.01 meV (CB at ~1.519 eV, VB at 0, SO at ~0.34 eV)
- Happy path: InAs k=0 eigenvalues agree (smaller gap, larger SO splitting)
- Happy path: GaAs [100] dispersion shows matching CB effective mass to < 1%
- Edge case: InSb (very small gap) — dispersion comparison at small |k| only (large non-parabolicity)
- Error path: if k=0 fails for any material, report which eigenvalue(s) differ and by how much

**Verification:**
- All bulk k=0 tests pass for all ported materials
- Bulk dispersion effective masses agree to < 1%
- Zeeman splitting agrees to < 0.01 meV (or test skipped with documented reason if ge mapping unresolved)

---

- U6. **QW Verification Suite**

**Goal:** Implement QW subband, dispersion, and convergence comparisons.

**Requirements:** R10, R11, R12

**Dependencies:** U5

**Files:**
- Create: `validation/qw/test_qw_subbands.py`
- Create: `validation/qw/test_qw_dispersion.py`
- Create: `validation/qw/test_qw_convergence.py`

**Approach:**
- `test_qw_subbands.py`: For GaAs/AlGaAs (5, 7, 10, 15 nm wells), InAs/GaSb, InAs/GaAs strained QW. Run both codes, compare subband energies at k_par=0. Tolerance: 1 meV at converged resolution.
- `test_qw_dispersion.py`: In-plane k_par sweep for GaAs/AlGaAs QW. Compare subband curvature (effective masses) to < 2%.
- `test_qw_convergence.py`: Richardson extrapolation study. Our code at FD orders 2, 4, 6, 8 and FDstep sweep. kdotpy at zres = 0.5, 0.25, 0.125, 0.0625 nm. Extrapolate independently, compare continuum-limit values. Minimum 4 grid spacings per code.
- QW configuration in kdotpy: `PhysParams(kdim=2, m_layers=[barrier, well, barrier], l_layers=[lbarr, lwell, lbarr], zres=0.25, norbitals=8)`

**Test scenarios:**
- Happy path: GaAs/AlGaAs 10nm QW CB1 subband energy agrees to < 1 meV
- Happy path: HH1-LH1 splitting matches between codes
- Edge case: narrow QW (5 nm) — larger discrepancy expected due to discretization differences
- Integration: Richardson extrapolated values from both codes agree to < 0.5 meV

**Verification:**
- QW subband energies agree to < 1 meV at converged resolution
- Richardson-extrapolated continuum values agree between codes
- Effective masses from dispersion agree to < 2%

---

- U7. **Wire + Landau Level Verification**

**Goal:** Compare wire subbands and Landau level fan diagrams.

**Requirements:** R13, R14, R15, R16

**Dependencies:** U6

**Files:**
- Create: `validation/wire/test_wire_subbands.py`
- Create: `validation/landau/test_landau_qw.py`
- Create: `validation/landau/test_landau_bulk.py`

**Approach:**
- Wire: kdotpy 1D mode (`PhysParams(kdim=1, width=..., yres=...)`) vs our confinement=2. Rectangular wire, GaAs material. Run convergence pilot first to establish achievable agreement. Target: < 2 meV at convergence.
- Landau QW: kdotpy LL mode vs our confinement=3. GaAs QW with B perpendicular, 0.1-10 T, 10+ B values. Compare LL energies at matched grid resolution.
- Landau bulk: kdotpy bulk-ll vs our bulk with magnetic field. GaAs at representative B values.

**Test scenarios:**
- Happy path: wire subband energies agree to < 2 meV at convergence
- Happy path: LL fan diagram E(B) shows matching level crossings
- Edge case: high B (>5 T) — LL mixing may differ between discretizations
- Integration: LL degeneracy factor consistent between codes

**Verification:**
- Wire energies agree within established tolerance
- LL energies at 5+ B values agree to < 1 meV

---

- U8. **g-Factor, Strain, and Berry Curvature Verification**

**Goal:** Compare BHZ/g-factor, strained subbands, Berry curvature, and LL Chern numbers.

**Requirements:** R17, R18, R19, R20, R21, R22

**Dependencies:** U6

**Files:**
- Create: `validation/gfactor/test_gfactor_qw.py`
- Create: `validation/gfactor/test_bhz_extraction.py`
- Create: `validation/strain/test_strain_bandedge.py`
- Create: `validation/strain/test_strain_qw.py`
- Create: `validation/berry/test_berry_curvature.py`
- Create: `validation/berry/test_chern_ll.py`

**Approach:**
- g-factor: Our Lowdin g_x, g_y, g_z vs kdotpy's Zeeman splitting derivative. Standard-star QWs. Tolerance: < 5%.
- BHZ: Run kdotpy's `do_bhz()` on GaAs/AlGaAs QW, compare A, B, C, D, M parameters against our perturbation results. Tolerance: < 10% (accounts for commutator vs dH/dk methodology difference).
- Strain bandedge: Compare Bir-Pikus shifts against analytical formulas in both codes. Then compare strained subband energies between codes to < 1 meV.
- Berry: Compare Omega(k) at QW k-points. Tolerance: < 10% (numerical differentiation sensitivity). LL Chern numbers: same integer at matching B.
- Strain requires substrate material in kdotpy: `PhysParams(substrate_material=mat_gaas, rel_strain=...)`

**Test scenarios:**
- Happy path: GaAs/AlGaAs QW g_z agrees to < 5%
- Happy path: BHZ A parameter agrees to < 10% for 10 nm well
- Happy path: strained InAs/GaAs band edge shifts match analytical Bir-Pikus
- Edge case: narrow QW g-factor — larger tolerance acceptable
- Integration: strained QW subbands agree to < 1 meV

**Verification:**
- g-factors agree to < 5%
- BHZ parameters agree to < 10%
- Strained subband energies agree to < 1 meV
- Berry curvature pattern matches, Chern numbers agree

---

- U9. **Self-Consistent Verification**

**Goal:** Compare converged self-consistent solutions.

**Requirements:** R23

**Dependencies:** U6, U8 (strain)

**Files:**
- Create: `validation/selfconsistent/test_sc_qw.py`

**Approach:**
- Single test case: doped GaAs/AlGaAs QW. Match temperature, k_parallel grid, convergence tolerance, boundary conditions.
- Our code: SC loop with DIIS mixing. kdotpy: SC loop with dynamic time stepping.
- Compare converged V(z) profiles, subband energies, charge density n(z).
- Tolerances relaxed from initial requirements: < 5 meV potential, < 5 meV subband, < 10% charge density (different mixing strategies converge to slightly different states).
- Run both SC solvers on the same QW first to verify profile shape agreement before comparing absolute values.

**Test scenarios:**
- Happy path: converged potential profile shapes match (qualitative)
- Happy path: CB1 subband energy agrees to < 5 meV
- Edge case: convergence failure in either code — report as discrepancy
- Integration: full SC cycle completes in both codes, results comparable

**Verification:**
- V(z) profile shapes match qualitatively
- Subband energies agree to < 5 meV
- Charge density integrates to same total within 10%

---

- U10. **Pipeline Integration and Discrepancy Resolution**

**Goal:** Build the full pipeline runner that orchestrates all tests and produces a unified report.

**Requirements:** R24, R25, R26

**Dependencies:** U5 through U9

**Files:**
- Create: `validation/run_all.py`

**Approach:**
- `run_all.py` orchestrates the layered verification order from R26:
  1. Bulk (U5) — gate
  2. QW (U6) — after bulk passes
  3. Wire + Landau (U7), g-factor/strain/Berry (U8) — after QW, in parallel
  4. SC (U9) — after QW + strain
- Enforces strict discrepancy protocol: stops at each layer boundary, reports all discrepancies, requires resolution before proceeding
- Discrepancy report format: JSON with observable, our_value, kdotpy_value, delta, tolerance, status (PASS/FAIL/INVESTIGATE), resolution (if any)
- Aggregates all markdown benchmark tables into a unified report
- Exit code: 0 if all pass, 1 if any unresolved discrepancy

**Test scenarios:**
- Happy path: full pipeline runs, produces unified report, all comparisons pass
- Error path: bulk k=0 fails — pipeline stops, discrepancy logged, resolution required
- Integration: `run_all.py` can be called from CI (future)

**Verification:**
- Single command `python validation/run_all.py --build-dir build/` runs the full suite
- Unified markdown report produced at `validation/results/cross_code_report.md`
- Exit code reflects overall pass/fail status

---

## System-Wide Impact

- **Interaction graph:** New validation scripts read our Fortran output format and call kdotpy library. No changes to Fortran production code except U4 (new materials in parameters.f90).
- **Error propagation:** Parameter mapping errors in U2 propagate through all downstream units. The bulk k=0 test (U5) is the gate that catches mapping bugs.
- **State lifecycle risks:** kdotpy's `initialize_materials()` modifies global state (`allMaterials`). Tests must call init once at module level, not per-test.
- **API surface parity:** Not applicable — this is a validation pipeline, not an API.
- **Integration coverage:** Each unit tests the full pipeline from material creation through kdotpy invocation through eigenvalue comparison.
- **Unchanged invariants:** Our Fortran code's Hamiltonian construction, eigensolvers, and output format are unchanged. kdotpy's computation is unchanged. Only new materials (U4) and validation scripts are added.

---

## Risks & Dependencies

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| F/kappa mapping wrong despite analysis | Low | Critical — blocks all downstream | U5 bulk k=0 gate catches it immediately; per-band eigenvalue comparison isolates which parameter is wrong |
| kdotpy API changes across versions | Low | Medium | Pin to v1.3.1 in venv |
| P conversion precision mismatch | Medium | Medium — could cause ~0.1 meV bulk discrepancy | Use kdotpy's exact hbarm0 constant, not our computed value |
| QW convergence too slow for practical testing | Medium | Low | Start with coarse grids, refine only for discrepancy investigation |
| kappa=0 introduces measurable QW error | Low | Low | kappa=0 only affects interface terms; quantify impact by running with/without kappa for HgTe/CdTe (native kdotpy materials where kappa is nonzero) |
| Our code has bugs that kdotpy reveals | Possible | High (but desirable) | Strict discrepancy protocol: investigate, fix, document |
| EC/EV alignment differs from Eg | Medium | Medium — wrong band gap | Verify at k=0 first; Ec-Ev must equal Eg*1000 meV |

---

## Phased Delivery

### Phase 1: Foundation + Bulk (U1-U5)
Highest leverage — validates parameter mapping, establishes infrastructure, provides first physics results. After Phase 1, we know whether the codes agree on the 8x8 Hamiltonian.

### Phase 2: QW + Convergence (U6)
QW is the most important confined system. Convergence study (R12) is critical for establishing whether grid-dependent comparisons are meaningful.

### Phase 3: Extended Physics (U7-U9)
Wire, Landau, g-factor, strain, Berry, SC. These can proceed in parallel after Phase 2.

### Phase 4: Pipeline Integration (U10)
Automates the full suite for CI integration and reproducibility.

---

## Sources & References

- **Origin document:** `docs/brainstorms/kdotpy-cross-validation-requirements.md`
- Our Hamiltonian: `src/physics/hamiltonianConstructor.f90`
- Our parameters: `src/core/parameters.f90`
- kdotpy blocks: `validation/kdotpy_env/lib/python3.14/site-packages/kdotpy/hamiltonian/blocks.py`
- kdotpy materials: `validation/kdotpy_env/lib/python3.14/site-packages/kdotpy/materials/materials.py`
- Test infrastructure: `tests/integration/star_helpers.py`
- kdotpy paper: Beugeling et al., SciPost Phys. Codebases 47 (2025), arXiv:2407.12651
