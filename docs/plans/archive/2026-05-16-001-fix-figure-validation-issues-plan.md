---
title: "fix: Figure Validation Remaining Fixes"
type: fix
status: active
date: 2026-05-16
origin: docs/brainstorms/figure-validation-fixes-requirements.md
---

# Fix: Figure Validation Remaining Fixes

## Summary

Fix 6 HIGH/STALE and 10 MEDIUM-severity figure issues identified in the validation audit. Eight implementation units span Fortran optics code (bulk TE=TM, magnitude normalization, ISBT transitions), Python figure-generation scripts (executable routing, column mapping, width sweep, strained config), and lecture markdown text alignment. All fixes have clear root causes from prior investigation.

---

## Problem Frame

The figure validation audit (Phases 1-2) classified all 104 figures by severity. Phase 3 partially completed: strain figures regenerated, lecture 04 HH/LH text fixed, orphaned figures removed. The remaining 16 issues ship physically incorrect figures (bulk TE≠TM, ISBT reading stale data, scattering lifetimes off by orders of magnitude) and lecture text contradicting its own figures.

---

## Requirements

**Origin requirements traced:**

- R1. Bulk absorption TE=TM for cubic materials
- R2. Bulk absorption magnitude (missing BZ normalization)
- R3. Wire `compute_intersubband_transitions` into `main_optics.f90`
- R4. ISBT/gain figures call wrong executable
- R5. Scattering lifetime vs width: wrong column mapping
- R6. Scattering lifetime vs field: stale cached data
- R7. Absorption vs width: implement sweep loop
- R8. Strained absorption config: switch to strained material system
- R9-R18. Ten MEDIUM-severity text-figure mismatches in lectures 03, 05, 06, 07, 08, 09, 11

---

## Scope Boundaries

- No changes to material parameters in `parameters.f90` without published reference
- No changes to basis ordering (bands 1-4 valence, 5-6 split-off, 7-8 conduction)
- No 3D spherical integration refactor for bulk optics (isotropic averaging sufficient)
- No new test cases for fixed figures
- No changes to input parser sequential reading behavior

---

## Context & Research

### Relevant Code and Patterns

- `src/physics/optical_spectra.f90` — `optics_accumulate` (line 105), `optics_finalize` (line 357), `compute_isbt_absorption` (line 726), `compute_intersubband_transitions` (line 663)
- `src/apps/main_optics.f90` — bulk k-sweep (line 116), 3D weight (line 219), ISBT call site (line 526), QW k-loop (line 508)
- `src/core/defs.f90` — `optics_config` type (line 240), `simulation_config` type
- `scripts/plotting/generate_all_figures.py` — `fig_isbt_absorption` (line 4243), `fig_gain_strained_comparison` (line 4289), `fig_scattering_lifetime_vs_width` (line 4754), `fig_scattering_lifetime_vs_field` (line 4824), `fig_qw_absorption_vs_width` (line 4085)
- `tests/regression/configs/` — existing config files for strained/unstrained comparisons
- Width-sweep pattern: `fig_exciton_binding_vs_width` (line 4374) with helper `_run_exciton_width` (line 4324)

### Institutional Learnings

- Scattering code had unit conversion bugs fixed in commit c160eff (J/eV coupling, k_lo units, form factor normalization). Cached sweep data predates this fix.
- The `optics_config` type has no confinement field — the subroutine cannot distinguish bulk from QW without modification.

---

## Key Technical Decisions

- **Isotropic averaging for bulk TE=TM**: Add `confinement` field to `optics_config` type. In `optics_accumulate`, when confinement=0, use `(px+py+pz)` for both TE and TM instead of `px+py` vs `pz`. This is physically correct for cubic zincblende and simpler than directional k-sweep integration. (see origin: key decisions)
- **Strained absorption config uses existing Ga47In53AsW on GaAs substrate**: The existing In0.53Ga0.47As material has a large lattice mismatch with GaAs (~3.8%), producing significant strain effects. This avoids adding new material parameters to `parameters.f90`. Use GaAs (or Al30Ga70As) as barrier with `strain: T` and `strain_ref: GaAs`.
- **ISBT transitions called outside k-loop**: `compute_intersubband_transitions` needs z_grid/dz/fdstep (not available in k-loop) and computes a table at k=0. Call it once after the k-loop, before `optics_finalize`.
- **Width sweep follows existing exciton pattern**: Loop over widths, write config per width, run executable, cache results, plot multi-curve comparison.

---

## Open Questions

### Resolved During Planning

- How to pass confinement to `optics_accumulate`: add field to `optics_config` rather than new argument (cleaner, no signature change needed)
- Material system for strained absorption: Ga47In53AsW on GaAs substrate (existing material, large strain, no new parameters)
- Where to call `compute_intersubband_transitions`: outside k-loop at k=1, before `optics_finalize`

### Deferred to Implementation

- Exact InGaAs well width for strained absorption config — must be within critical thickness for In0.53Ga0.47As on GaAs (very thin, ~5-10 nm max given 3.8% mismatch)
- Whether Ga47In53AsW at 3.8% mismatch produces numerically stable strain fields — may need to check strain solver convergence
- Exact text wording for the 10 MEDIUM-severity lecture edits — the approach is defined, wording is an execution detail

---

## Implementation Units

- U1. **Bulk optics isotropic averaging (R1, R2)**

**Goal:** Fix bulk TE=TM asymmetry and magnitude underestimation in the optics code.

**Requirements:** R1, R2

**Dependencies:** None

**Files:**
- Modify: `src/core/defs.f90` (add `confinement` field to `optics_config`)
- Modify: `src/physics/optical_spectra.f90` (isotropic averaging in `optics_accumulate`)
- Modify: `src/apps/main_optics.f90` (set confinement in optcfg, fix 3D weight)

**Approach:**
1. Add `integer :: confinement = 0` field to `optics_config` type in `defs.f90`
2. In `main_optics.f90`, after constructing `optcfg`, set `optcfg%confinement = cfg%confinement`
3. In `optics_accumulate` (`optical_spectra.f90` line 174): when `optcfg%confinement == 0`, use `v_total = px + py + pz` for BOTH TE and TM accumulation instead of `px+py` for TE and `pz` for TM
4. Fix 3D weight at `main_optics.f90` line 221: change `/(2.0_dp * pi_dp)` to `/(2.0_dp * pi_dp)**3` for proper `(2*pi)^3` BZ normalization

**Patterns to follow:**
- Existing `optics_config` field additions pattern in `defs.f90`
- QW weight computation at line 454 (`2*pi*k*dk` without BZ divisor) as reference for correct normalization

**Test scenarios:**
- Happy path: bulk GaAs absorption TE spectrum equals TM spectrum (within numerical precision)
- Happy path: bulk GaAs absorption magnitude is in physically reasonable range (~10^3-10^4 cm^-1 above bandgap)
- Edge case: QW optics unchanged (confinement=1 path still uses `px+py` vs `pz`)

**Verification:**
- Build succeeds
- Bulk absorption figure shows TE=TM and reasonable magnitude
- QW absorption figures unchanged from current behavior

---

- U2. **Wire ISBT transitions into optics executable (R3)**

**Goal:** Call `compute_intersubband_transitions` from `main_optics.f90` so `isbt_transitions.dat` is produced.

**Requirements:** R3

**Dependencies:** U1 (same file `main_optics.f90`)

**Files:**
- Modify: `src/apps/main_optics.f90` (add call after k-loop)

**Approach:**
1. After the QW k-loop (after line 532), before `optics_finalize` (line 537), add a conditional call:
   ```
   if (cfg%optics%isbt_enabled) then
     call compute_intersubband_transitions(eig(:, 1), eigv(:, :, 1), &
       cfg%z, cfg%dz, cfg%numcb, cfg%numvb, cfg%fdstep, &
       "output/isbt_transitions.dat")
   end if
   ```
2. This uses k=1 eigenvalues/eigenvectors (zone-center), which is where ISBT transitions are most meaningful

**Patterns to follow:**
- Existing ISBT call at line 526 for `compute_isbt_absorption` (same conditional, same scope)

**Test scenarios:**
- Happy path: `output/isbt_transitions.dat` is created when ISBT mode is active
- Edge case: no crash when `isbt_enabled = .false.`
- Integration: file contains CB-CB transition energies in the expected range (80-120 meV for 100A GaAs/AlGaAs QW)

**Verification:**
- Build succeeds
- `isbt_transitions.dat` file produced with reasonable values

---

- U3. **Fix ISBT and gain figure executables (R4)**

**Goal:** Fix `fig_isbt_absorption` and `fig_gain_strained_comparison` to call the correct executable.

**Requirements:** R4

**Dependencies:** None

**Files:**
- Modify: `scripts/plotting/generate_all_figures.py` (lines ~4253, ~4289)

**Approach:**
1. In `fig_isbt_absorption` (~line 4253): change `run_executable(EXE_BAND, ...)` to `run_executable(EXE_OPTICS, ...)`
2. In `fig_gain_strained_comparison` (~line 4289): same change `EXE_BAND` to `EXE_OPTICS`
3. Verify the config files used by these functions have correct `optics:` blocks

**Test scenarios:**
- Happy path: ISBT figure shows absorption in 20-300 meV range (not 1.3-2.0 eV)
- Happy path: gain figure shows gain spectrum, not empty/stale data

**Verification:**
- ISBT figure energy axis shows meV range
- Gain figure shows physically meaningful gain spectrum

---

- U4. **Fix scattering lifetime figures (R5, R6)**

**Goal:** Fix column mapping in vs_width and clear stale cache for vs_field.

**Requirements:** R5, R6

**Dependencies:** None

**Files:**
- Modify: `scripts/plotting/generate_all_figures.py` (lines ~4782-4793)
- Delete: `output/scattering_field_sweep.dat` (stale cache)

**Approach:**
1. In `fig_scattering_lifetime_vs_width`: fix column mapping
   - Column 2 = E_ij in meV (transition energy, for x-axis)
   - Columns 3-4 = rate_emission, rate_absorption in 1/s (sum for total rate)
   - Columns 5-6 = tau_emission, tau_absorption in ps (use directly for y-axis)
   - Replace the nonsensical `lifetime = 1.0 / rate` with direct use of columns 5-6
2. Delete stale `output/scattering_field_sweep.dat` so the field sweep re-runs with corrected Fortran scattering code

**Patterns to follow:**
- `_run_scattering_field_sweep` (lines 4874-4878) already has correct column mapping — follow its pattern for the vs_width function

**Test scenarios:**
- Happy path: scattering lifetimes plot shows 1-1000 ps range (not 10^8 ps)
- Happy path: x-axis shows transition energies in meV (not subband indices)
- Integration: field sweep shows physically reasonable trend (lifetime changes with field)

**Verification:**
- Scattering lifetime figures show physically reasonable magnitudes

---

- U5. **Implement absorption vs width sweep (R7)**

**Goal:** Replace passive file-reader with active width-sweep loop that generates multi-curve comparison.

**Requirements:** R7

**Dependencies:** U3 (needs working optics executable)

**Files:**
- Modify: `scripts/plotting/generate_all_figures.py` (lines 4085-4167)

**Approach:**
1. Add helper `_run_absorption_width_sweep` following the `_run_exciton_width` pattern:
   - Define widths: `widths_aa = [30, 50, 80, 100, 150, 200]`
   - For each width: write config with `Optics: T`, appropriate material layers, compute half/barrier/fdstep
   - Run `bandStructure`, read `output/absorption_TE.dat` and `output/absorption_TM.dat`
   - Cache each width's spectra as `absorption_TE_W{width}.dat` / `absorption_TM_W{width}.dat`
2. In `fig_qw_absorption_vs_width`: check for cached sweep data first, then run sweep if needed, then plot all curves

**Patterns to follow:**
- `fig_exciton_binding_vs_width` (line 4374) for cache-then-run pattern
- `_run_exciton_width` (line 4324) for per-width config generation

**Test scenarios:**
- Happy path: figure shows 4+ absorption curves at different well widths
- Edge case: cached data is reused on second run (no redundant computation)
- Integration: absorption peaks shift to lower energy with increasing width (quantum confinement trend)

**Verification:**
- Figure shows multiple color-coded curves with width legend

---

- U6. **Strained absorption config (R8)**

**Goal:** Create strained absorption config using Ga47In53AsW on GaAs substrate where strain effects are significant.

**Requirements:** R8

**Dependencies:** None

**Files:**
- Modify: `tests/regression/configs/qw_ingaas_algaas_strained_absorption.cfg`
- Modify: `scripts/plotting/generate_all_figures.py` (update `fig_qw_absorption_strained` if config path changes)

**Approach:**
1. Replace the current InGaAs/AlInAs-on-InP config with:
   - Barrier: `Al30Ga70As -200 200 0` (full domain, no strain)
   - Well: `Ga47In53AsW -25 25 1` (narrow InGaAs well, strain flag)
   - `strain: T`, `strain_ref: GaAs`, `strain_solver: pardiso`
   - This gives ~3.8% lattice mismatch → significant compressive strain in the InGaAs well
2. Keep the `Optics: T` block with appropriate energy range for InGaAs bandgap (~0.8-1.5 eV)
3. Update the unstrained comparison config to use the same materials but without `strain: T`

**Note:** The 3.8% mismatch is large. Well width must be kept thin (~5 nm) to stay within approximate critical thickness. If the strain solver produces issues at this mismatch, reduce well width or accept that this is a computational demonstration rather than an experimentally realizable structure.

**Test scenarios:**
- Happy path: strained absorption shows visible difference from unstrained (HH/LH splitting shifts absorption edge)
- Edge case: strain solver converges for the chosen well width

**Verification:**
- Strained vs unstrained figure shows clearly different absorption spectra

---

- U7. **Lecture text fixes (R9-R18)**

**Goal:** Align 10 lecture markdown text passages with their referenced figures.

**Requirements:** R9, R10, R11, R12, R13, R14, R15, R16, R17, R18

**Dependencies:** U1 (bulk optics fixes), U3 (ISBT fixes), U4 (scattering fixes) — to confirm what the corrected figures actually show before finalizing text

**Files:**
- Modify: `docs/lecture/03-wavefunctions.md` (R9: vb_hh_lh_mixing text ~line 199; R10: cb_parts_evolution ~line 535)
- Modify: `docs/lecture/05-gfactor.md` (R17: qw_gfactor_vs_width ~line 490)
- Modify: `docs/lecture/06-optical-properties.md` (R13: absorption_excitonic_TE ~line 615; R14: isbt_dipole_moments ~line 813; R15: exciton_bohr_vs_width ~line 979)
- Modify: `docs/lecture/07-self-consistent-sp.md` (R16: sc_delta_doped_potential ~line 670)
- Modify: `docs/lecture/08-quantum-wire.md` (R18: wire_inas_gaas_subbands ~line 692)
- Modify: `docs/lecture/09-numerical-methods.md` (R11: convergence_fd_order ~line 225; R12: timing_dense_vs_sparse ~line 568)
- Modify: `docs/lecture/11-convergence.md` (R11: convergence_fd_order ~line 124)

**Approach:**
For each MEDIUM issue, update text to match what the figure actually shows:

- **R9** (lec 03, ~line 199): The text says "majority-LH crossover" — the figure shows this crossover at ~67% LH. Text is approximately correct; add clarification that the crossover is at 67% LH (not 100%).
- **R10** (lec 03, ~line 535): Text says "CB purity 67% at k=0" but the figure shows ~81%. Update "67%" to match the figure value.
- **R11** (lec 09 ~line 225, lec 11 ~line 124): Fix Delta z = 3 A to Delta z = 5 A for FDstep=101 on 500 A domain. In lec 11, acknowledge order-2 is the only order with visible error.
- **R12** (lec 09, ~line 568): Fix "FEAST is approximately 15x faster" to match actual benchmark data (FEAST is ~4.4x slower for this problem size, not 15x faster).
- **R13** (lec 06, ~line 615): The caption already disclaims the excitonic marker is post-processed, not computed. Add a clearer note about the magnitudes being free-carrier (not excitonic).
- **R14** (lec 06, ~line 813): Clarify that the plotted values are velocity-gauge (dH/dk units), not length-gauge dipole moments.
- **R15** (lec 06, ~line 971/979): The 100-200 A Bohr radius for a 10 nm QW exceeds the bulk limit of 113 A. Verify figure and update the text range to match what the variational calculation actually produces.
- **R16** (lec 07, ~line 670): Replace "100 meV" with the actual computed notch depth from the figure (~9 meV).
- **R17** (lec 05, ~line 490): Update the sweep range text to match the actual data range in the figure.
- **R18** (lec 08, ~line 692): Add explicit note that the dispersion uses only 2 kz points — not suitable for branch tracking.

**Test scenarios:**
- Each text passage describes the corresponding figure accurately without contradiction
- No physics claims are introduced that aren't supported by the code output

**Verification:**
- Read each modified lecture section and confirm text matches figure content

---

- U8. **Regeneration and final verification**

**Goal:** Regenerate all affected figures and verify the full pipeline.

**Requirements:** All (R1-R18)

**Dependencies:** U1, U2, U3, U4, U5, U6, U7

**Files:**
- Generated: `docs/figures/*.png` (regenerated figures)
- Verify: all lecture markdowns

**Approach:**
1. Rebuild Fortran: `cmake --build build`
2. Run figure generation: `OMP_NUM_THREADS=12 python3 scripts/plotting/generate_all_figures.py`
3. Run regression tests: `OMP_NUM_THREADS=12 ctest --test-dir build`
4. Verify figure-to-text consistency for all corrected figures
5. Produce final consolidation report

**Test scenarios:**
- Integration: all figures generate without errors
- Integration: regression test suite passes (no new failures)
- Integration: no orphaned or missing figure references

**Verification:**
- `generate_all_figures.py` exits cleanly
- `ctest` reports all tests passed (or same pre-existing failures only)
- Visual spot-check of corrected figures

---

## Risks & Dependencies

| Risk | Mitigation |
|------|------------|
| 3.8% lattice mismatch exceeds strain solver stability | Reduce InGaAs well width to 25 A; if still unstable, use a thinner well or different material system |
| Bulk optics fix affects QW/wire paths | Confinement flag ensures only bulk path uses isotropic averaging; QW/wire unchanged |
| Width sweep is computationally expensive (~6 runs x ~2 min each) | Cache results for reuse; skip if cached data exists |
| `compute_intersubband_transitions` may crash on edge cases | Guard with `isbt_enabled` check; file I/O errors don't affect main optics output |
| MEDIUM text fixes may need re-adjustment after figure regeneration | Do text fixes after figures are finalized (U7 depends on U1-U6) |

---

## Sources & References

- **Origin document:** [docs/brainstorms/figure-validation-fixes-requirements.md](docs/brainstorms/figure-validation-fixes-requirements.md)
- Related code: `src/physics/optical_spectra.f90`, `src/apps/main_optics.f90`, `scripts/plotting/generate_all_figures.py`
- Related commit: c160eff (scattering unit fix)
