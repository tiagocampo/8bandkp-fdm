# Phase 4: Optics & Documentation Figures — Design Spec

Date: 2026-05-04
Status: DRAFT
Groups: #22 (QW tutorials), #26 (docs physics revamp), #50 (physics figures extended)
Execution order: Validation-first (ISBT fix → validate existing → new figures → chapter rebuilds)

---

## Context

Phase 4 closes the documentation gap: all optics/physics capabilities need verified figures
and up-to-date lecture chapters. Current state:

- `generate_all_figures.py`: 74 figure functions, fully implemented, no stubs
- `docs/figures/`: 76 PNG files, all 70 lecture references resolve correctly
- `docs/lecture/figures/`: only 8 files (topology verification), Group #50 not started
- ISBT has two implementations with no cross-validation (z-dipole vs commutator velocity)
- Group #22 REVIEW.md entry is stale (says 0/11 figures but all exist)
- Groups #26 Tasks 6,10-11 and Group #50 entirely unbegun

---

## Work Package 1: ISBT Cross-Validation & Sign Fix

**Group:** #26 Task 6
**Effort:** Medium
**Blocks:** WP2, WP4

### Problem

`optical_spectra.f90` has two ISBT code paths:

1. `compute_intersubband_transitions()` (~line 663): uses `z_dipole()` helper
   - Computes `z_ij = dz * sum_n conjg(psi_i_n) * z_n * psi_j_n`
   - Outputs transition table to `isbt_transitions.dat`

2. `compute_isbt_absorption()` (~line 726): uses commutator velocity
   - Computes `v_z = -i[z, H]` via CSR SpMV on velocity matrix
   - Outputs absorption spectrum to `absorption_ISBT.dat`

Theoretical relationship: `v_z_ij = -i * E_ij * z_ij`, therefore
`|z_ij|^2 = |v_z_ij|^2 / E_ij^2`.

No cross-validation exists between these two methods. The sign or normalization
may be inconsistent.

### Validation Plan

1. Create GaAs/Al₀.₃Ga₀.₇As QW config (10nm well, L/2 barriers):
   - Well-known system, n=1→n=2 CB ISBT at ~150 meV
   - Use `opticalProperties` executable with `ISBT: T` flag

2. Extract from both methods:
   - From `isbt_transitions.dat`: z_ij dipole matrix elements and E_ij energies
   - From `absorption_ISBT.dat`: absorption peak position and intensity

3. Cross-check: `|z_ij|^2 * E_ij^2` vs `|v_z_ij|^2` for each transition
   - Tolerance: 1% (accounts for FD discretization)
   - If discrepancy > 1%: trace root cause and fix

4. Fix categories:
   - Sign error → flip convention in one method
   - Missing E_ij factor → add energy factor to commutator-based formula
   - Wrong normalization → fix dz sum in z_dipole

### Testing

- New pFUnit test: `test_isbt_dipole_velocity_consistency` in `test_optical_qw.pf`
  - Build QW Hamiltonian, compute both z-dipole and velocity for CB1→CB2
  - Assert `|z_ij|^2 * E_ij^2 ≈ |v_z_ij|^2` within 1%
- New regression config: `tests/regression/configs/optics_qw_isbt_crossval.cfg`
- Update `test_optical.pf` if any Fortran fixes are needed

### Deliverable

- ISBT absorption from both methods agrees within tolerance
- Unit test passes
- No ISBT figure changes needed (yet — WP4 will regenerate)

---

## Work Package 2: Validate Group #22 Figures

**Group:** #22 (QW tutorials optics figures)
**Effort:** Low
**Depends on:** WP1

### Scope

All 11 optics figure functions already exist in `generate_all_figures.py` and PNGs
exist in `docs/figures/`. This is a verification pass to catch regressions from WP1.

### Validation Steps

1. **Full run**: `python scripts/plotting/generate_all_figures.py` — verify no crashes
2. **Spot-check 11 optics figures**:
   - Absorption TE/TM (bulk GaAs, QW GaAs/AlGaAs)
   - ISBT absorption (QW)
   - Gain spectra (QW, TE+TM, with carrier density sweep)
   - Spontaneous emission
   - Spin-resolved absorption (up/down)
   - Exciton absorption (if exciton path works)
3. **Physics sanity checks**:
   - Absorption edge at E_g (GaAs ≈ 1.42 eV)
   - ISBT peak at expected energy (~150 meV for 10nm QW n=1→n=2)
   - Gain peak blue-shifts with increasing carrier density
4. **Update REVIEW.md**: Mark Group #22 as COMPLETE with verification notes

### Deliverable

- REVIEW.md Group #22 updated to COMPLETE
- All 11 optics figures verified (or issues documented)

---

## Work Package 3: Group #50 Physics Figures

**Group:** #50 (physics figures extended)
**Effort:** Medium
**Depends on:** WP1 (ISBT fix, for consistent data)
**Independent of:** WP2, WP4

### Approach

Extend existing `generate_all_figures.py` with 5 new figure functions.
Each follows the established pattern: build config → run executable → parse output → plot → save.

Output directory: `docs/lecture/figures/` (as per Group #50 plan).

### Figure 1: `bandstructure_bulk_ek`

**Purpose:** Bulk 8-band E(k) along high-symmetry directions for 3 materials.
**Executable:** `bandStructure` with `confinement=0`
**Layout:** 1×3 panel (GaAs, InAs, InSb)
**Content:**
- k along [100] (Γ→X), [110] (Γ→K), [111] (Γ→L)
- All 8 bands: HH (bands 1-2), LH (3-4), SO (5-6), CB (7-8)
- Spin-degenerate at k=0, splitting visible at finite k
- Axes: k in Å⁻¹, E in eV (relative to VBM)
**Config:** 3 bulk configs with wavevector sweep along each direction

### Figure 2: `bandstructure_qw_subbands`

**Purpose:** QW subband dispersion E(k_∥) showing non-parabolicity.
**Executable:** `bandStructure` with `confinement=1`
**Layout:** 2-panel: CB subbands (top), VB subbands (bottom)
**Content:**
- GaAs/Al₀.₃Ga₀.₇As QW, 10nm well
- k_∥ from 0 to 0.1 Å⁻¹
- First 3 CB subbands + first 6 VB subbands
- Show HH-LH mixing at finite k_∥ (anti-crossing)
**Config:** QW config with k_∥ sweep

### Figure 3: `wavefunctions_qw`

**Purpose:** Envelope functions |ψ(z)|² overlaid on band profile.
**Executable:** `bandStructure` with `confinement=1`, output eigenfunctions
**Layout:** 2-panel: CB states (top), VB states (bottom)
**Content:**
- First 3 CB states + first 3 VB states
- Band edge profile (CB/VB offset) as background shading
- Show confinement and barrier penetration
- Label each state with subband index and energy
**Config:** QW config with eigenfunction output enabled

### Figure 4: `wire_geometry_potential`

**Purpose:** Wire cross-section geometry and potential profile.
**Executable:** `bandStructure` in wire mode
**Layout:** 2-panel: material map (left), radial profile (right)
**Content:**
- InAs/GaAs hexagonal core-shell wire
- Material regions color-coded (core=InAs, shell=GaAs)
- Hexagonal grid overlay showing FD discretization
- Right panel: potential/strain along radial cut through center
**Config:** Wire config from existing regression test

### Figure 5: `landau_fan_diagram`

**Purpose:** Zeeman-split CB levels vs B showing spin splitting and effective mass.
**Executable:** `bandStructure` with `confinement=0` + B-field (Zeeman only)
**Note:** Orbital Landau levels (Peierls/confinement=3) are Phase 5 scope. This figure
shows Zeeman spin splitting of the CB edge as a function of B, which is currently implemented.
**Layout:** Single plot
**Content:**
- InAs CB Zeeman splitting: E_CB± = E_g ± g·μ_B·B/2
- B from 0 to 10 T
- Overlay analytical Zeeman formula with g* = 2.0 (default) and g* = -14.9 (InAs)
- Extract effective g-factor from slope, compare against known InAs g* ≈ -14.9
- Second panel (optional): QW Zeeman splitting at k=0 for GaAs/AlGaAs QW
**Config:** Bulk InAs with B-field sweep (loop over B values), parsed from eigenvalues.dat

### Registration

All 5 functions registered in `ALL_FIGURES` dictionary at end of
`generate_all_figures.py`, callable via `--only` flag.

### Deliverable

- 5 new PNG files in `docs/lecture/figures/`
- 5 new figure functions in `generate_all_figures.py`
- All callable via `python scripts/plotting/generate_all_figures.py --only FIG`

---

## Work Package 4: Chapter Rebuilds (Ch06 + Ch08)

**Group:** #26 Tasks 10-11
**Effort:** Medium
**Depends on:** WP1 (ISBT fix), WP3 (new figures)

### Ch06: Optical Properties

1. **Identify all figure references** in `docs/lecture/06-optical-properties.md`
2. **Regenerate each referenced figure** using corrected code (post-WP1)
3. **Verify captions match**:
   - Peak positions (absorption edge, ISBT energy, gain peak)
   - Units (cm⁻¹ vs eV — must be consistent)
   - TE/TM polarization labels match actual polarization
   - Line styles/colors match legend descriptions
4. **Fix text discrepancies**: update any caption or paragraph that describes
   spectral features that shifted after ISBT fix

### Ch08: Exciton & Scattering

1. **Identify all figure references** in `docs/lecture/08-exciton-scattering.md`
2. **Regenerate figures**: exciton absorption, LO-phonon scattering rates
3. **Verify physics**:
   - Exciton peak at expected binding energy (~4-10 meV for GaAs QW)
   - LO-phonon scattering rate peak at ℏω_LO ≈ 36 meV (GaAs)
   - Scattering rate magnitude order-of-magnitude correct
4. **Fix text discrepancies**: update descriptions to match actual figures

### Process

For each chapter:
1. `grep '!\[.*\](.*\.png)' chapter.md` to extract all figure references
2. For each reference: regenerate figure, compare old vs new PNG
3. Read surrounding text, verify description matches new figure
4. Edit markdown if needed

### Deliverable

- Ch06 + Ch08 markdown updated (captions, descriptions)
- All referenced figures regenerated with corrected code
- No broken figure references

---

## Summary

| WP | Group | Effort | Depends on | Deliverable |
|----|-------|--------|------------|-------------|
| 1. ISBT fix | #26 T6 | Medium | — | Cross-validated ISBT, unit test |
| 2. Validate #22 | #22 | Low | WP1 | REVIEW.md updated |
| 3. New figures | #50 | Medium | WP1 | 5 new physics figures |
| 4. Chapter rebuild | #26 T10-11 | Medium | WP1, WP3 | Ch06+Ch08 verified |

**Estimated total effort:** 3-4 days
