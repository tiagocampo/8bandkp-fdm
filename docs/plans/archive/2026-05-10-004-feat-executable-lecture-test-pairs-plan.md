---
name: executable-lecture-test-pairs
date: 2026-05-10
origin: docs/brainstorms/2026-05-10-executable-lecture-test-pairs-requirements.md
status: completed
---

# Executable Lecture-Test Pairs — Implementation Plan

## Problem Frame

The codebase has 14 lecture docs with detailed derivations, 20+ regression configs, 4 standalone verification scripts in `scripts/`, and 20+ integration verification scripts in `tests/integration/` (verification ladder rungs 1-4, standard-star benchmarks S1-S7). These are disconnected: reading a lecture doc doesn't verify the code, running a regression test doesn't reference the physics derivation, and the verification scripts aren't tied to any lecture. This plan creates 14 executable lecture-test pairs that connect each lecture's derivations to code output through standalone Python scripts.

*(see origin: docs/brainstorms/2026-05-10-executable-lecture-test-pairs-requirements.md)*

## Scope

**In scope:** 14 lecture-pair scripts, Makefile targets, overlay plots, lecture doc updates (embed plots, add code-output anchors), refactoring 4 existing `scripts/verify_*.py` scripts.

**Out of scope:** Jupyter notebooks, CI integration of lecture pairs, new shared framework module, new physics modules, changes to regression test infrastructure, changes to `generate_all_figures.py`.

*(see origin: Scope Boundaries)*

## Requirements Traceability

| Req | Description | Implementation Units |
|-----|-------------|---------------------|
| R1 | Standalone script per lecture (00-13) with config, run, parse, analytical reference, overlay plot, assertion | U1-U14 |
| R2 | Reuse existing regression configs; add new configs only for gaps | U1-U14 (per-lecture config mapping) |
| R3 | Fix-as-you-go: physics matches before advancing | Protocol in §Fix-as-You-Go |
| R4 | Makefile targets `make lecture-XX` and `make lectures` | U15 |
| R5 | Embed overlay plots in lecture docs | U16 |
| R6 | Update derivations with code-output anchors | U16 |
| R7 | Refactor 4 existing `scripts/verify_*.py` into lecture-pair format | U13 (lecture 13 absorbs all 4) |
| R8 | Methodological lectures (00, 12) get worked-example configs | U1, U14 |
| R9 | Lecture-specific tolerances, documented limitations | Per-unit tolerance specs |
| R10 | Rich diagnostic output: overlay plots + comparison tables | U1-U14 (each script) |
| R11 | Overlay plots committed under lecture figures directory | U1-U14 (output to `docs/lecture/figures/`) |

## Resolved Planning Questions

### Q1: Lecture-to-Infrastructure Mapping

| Lecture | Executable(s) | Existing Configs | Existing Verify Scripts | Analytical Reference Available? |
|---------|---------------|-----------------|------------------------|-------------------------------|
| 00 quickstart | `bandStructure` | `bulk_gaas_k0.cfg` | None | No (worked example per R8) |
| 01 bulk | `bandStructure` | `bulk_gaas_k0.cfg`, `bulk_gaas_kx.cfg`, `bulk_gaas_kx_dispersion.cfg`, `bulk_gaas_kxky.cfg`, `bulk_inas_kx_dispersion.cfg` | `verify_8band_rung1_bulk_k0.py`, `verify_8band_rung2_dispersion.py`, `verify_star_gaas_bulk.py`, `verify_star_inas_bulk.py`, `verify_star_insb_bulk.py` | Yes: k=0 eigenvalue self-check, Kane effective mass formula, Vurgaftman band gaps |
| 02 QW | `bandStructure` | `qw_gaas_algaas.cfg`, `qw_alsbw_gasbw_inasw.cfg`, `qw_gaas_algaas_double_qw.cfg` | `verify_8band_rung3_qw.py`, `verify_star_gaas_algaas_qw.py`, `verify_star_inas_gasb_qw.py` | Partial: QW subbands lack clean analytical formula (Bastard fails for 8-band); use inter-QW consistency checks and regression comparison |
| 03 wavefunctions | `bandStructure` | `qw_gaas_algaas.cfg` | `verify_qw_state_character.py` | Yes: band character fractions from eigenvector block decomposition |
| 04 strain | `bandStructure` | `bulk_gaas_strained.cfg`, `qw_inas_gaas_strained.cfg` | `verify_star_inas_gaas_qw.py` (strained mode) | Yes: Bir-Pikus band-edge shifts, Vurgaftman deformation potentials |
| 05 g-factor | `gfactorCalculation` | `gfactor_bulk_gaas_cb.cfg`, `gfactor_bulk_gaasw_cb.cfg`, `gfactor_bulk_inasw_cb.cfg`, `gfactor_bulk_insb_cb.cfg`, `gfactor_qw_cb.cfg` | `verify_landau_analytical.py`, `verify_star_gaas_bulk.py` (g-factor), `verify_star_insb_bulk.py` (g-factor) | Yes: Roth formula (GaAs CB g = -0.315), Landau level E_n = E_C + hbar*omega_c/2 |
| 06 optical | `opticalProperties` | `bulk_gaas_optics.cfg`, `qw_gaas_algaas_optics.cfg`, `qw_gaas_algaas_absorption.cfg`, `qw_gaas_algaas_isbt.cfg`, `qw_gaas_algaas_spin_resolved.cfg` | `verify_qw_absorption_polarization.py`, `verify_wire_optical_selection.py` | Partial: absorption edge at band gap, TE/TM polarization ratio, but no closed-form spectrum |
| 07 SC-SP | `bandStructure` | `sc_bulk_gaas_doped.cfg`, `sc_gaas_alas_qw.cfg`, `sc_delta_doped_gaas.cfg`, `sc_mod_doped_gaas_algaas.cfg` | `verify_sc_benchmarks.py`, `verify_star_gaas_algaas_qw.py` (SC mode) | Partial: convergence criteria, charge neutrality, but no closed-form for self-consistent potential |
| 08 wire | `bandStructure` | `wire_gaas_31x31.cfg`, `wire_inas_rectangle.cfg`, `wire_gaas_rectangle.cfg` | `verify_8band_rung4_wire.py`, `verify_star_inas_wire.py` | Partial: wire subband count vs confinement dimensions, internal consistency checks |
| 09 numerical methods | `bandStructure` | `bulk_gaas_k0.cfg` (vary FDorder) | `verify_8band_rung1_bulk_k0.py` (FD stencils) | Yes: Richardson extrapolation, convergence order vs theoretical FD order |
| 10 QCSE | `bandStructure` | `sc_qcse_gaas_algaas.cfg`, `sc_qcse_gaas_algaas_ef.cfg`, `qw_gaas_algaas_qcse_scattering.cfg` | `verify_stark_shift.py` | Yes: Stark shift = E(field) - E(0), quadratic vs linear field dependence |
| 11 convergence | `bandStructure` | `bulk_gaas_kx_dispersion.cfg`, `qw_gaas_algaas.cfg` (vary grid) | `verify_8band_rung2_dispersion.py` (mass extraction) | Yes: convergence rate = theoretical FD order (2, 4, 6, 8) |
| 12 extending | varies | None (demonstrative) | None | No (worked example per R8) |
| 13 topological | `topologicalAnalysis` | `topology_qhe_qwz.cfg`, `topology_qwz_chern_u-0.8.cfg`, `topology_bhz_z2_*.cfg`, `topology_qw_bdg.cfg`, `topology_rashba_phase.cfg`, `topology_spectral_*.cfg` | `scripts/verify_qwz_chern.py`, `scripts/verify_bhz_z2.py`, `scripts/verify_landau_levels.py`, `scripts/sweep_rashba_bdg.py` | Yes: QWZ Chern = integer, BHZ Z2 = 0 or 1, BdG gap = 2*Delta, Rashba phase boundary |

### Q2: Fix-as-You-Go Triage Protocol

When code output diverges from analytical derivation:

1. **Classify the divergence:**
   - **Code bug** — fix the Fortran code. Indicators: wrong sign, wrong basis ordering, indexing error, incorrect formula implementation.
   - **Model limitation** — document and adjust tolerance. Indicators: 8-band g-factor 20-30% shortfall vs experiment, spurious solutions at small grid spacing, non-parabolicity beyond Kane model.
   - **Analytical reference error** — update the lecture derivation. Indicators: derivation uses different sign convention, different parameter set, or simplified model that doesn't match the full 8-band treatment.

2. **Time-box:** If investigation exceeds 30 minutes without clear resolution, document the divergence as a known limitation with a TODO marker in the script and proceed to the next lecture. Return to unresolved items after all 14 pairs have initial scripts.

3. **Conflict resolution:** If a fix breaks an existing passing test, the existing test takes precedence. Investigate whether the lecture script's analytical reference is correct independently of the code. Do not modify existing passing tests to accommodate a new lecture pair.

### Q3: Shared Utilities — Import from star_helpers.py

The 10k-line `tests/integration/star_helpers.py` provides `run_executable`, `parse_eigenvalues`, `parse_gfactor`, `parse_absorption`, `parse_topology_result`, `compare_value`, tolerance tiers (`TOL_EXACT`, `TOL_ANALYTICAL`, `TOL_NUMERICAL`), and benchmark table formatting. This already solves the duplication problem that 14 standalone scripts would reproduce.

**Decision:** Lecture scripts import from `star_helpers.py` via `sys.path` prepend. No new shared module. Each lecture script remains standalone in the sense that it can be run independently via `make lecture-XX`, but it reuses existing parsers and runners rather than reimplementing them.

```python
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'tests', 'integration'))
from star_helpers import (
    run_exe, parse_eigenvalues, parse_gfactor, parse_absorption,
    parse_topology_result, compare_value, roth_gfactor,
    extract_effective_mass, format_benchmark_row, print_benchmark_header,
    TOL_EXACT, TOL_ANALYTICAL, TOL_NUMERICAL,
)
```

Key reusable functions: `run_exe` (dispatches by executable name), `roth_gfactor` (analytical Roth formula), `extract_effective_mass` (adaptive parabolic fit near k=0), plus all four parsers and `compare_value`.

This resolves the YAGNI tension: no new shared framework, but leverage the existing one that was already extracted for exactly this purpose.

### Q4: generate_all_figures.py Coexistence

`scripts/plotting/generate_all_figures.py` (6780 lines, 78 figure functions) produces publication-quality figures. Four functions target `docs/lecture/figures/`: `bandstructure_bulk_ek`, `bandstructure_qw_subbands`, `wavefunctions_qw`, `zeeman_fan_diagram`. The lecture-pair scripts produce **validation overlay plots** (code output vs analytical reference) — a different artifact.

**Decision:** Lecture-pair scripts coexist with `generate_all_figures.py`. They write overlay plots to `docs/lecture/figures/` with a `lecture_XX_` prefix to avoid name collisions. The existing figure functions in `generate_all_figures.py` remain untouched. No changes to `generate_all_figures.py`.

## Key Technical Decisions

1. **Script location:** `scripts/lecture_XX_<name>.py` — follows the existing `scripts/verify_*.py` pattern and is accessible from the Makefile.

2. **Config reuse:** Each lecture script references configs from `tests/regression/configs/` by absolute path (constructed from repo root). New configs are only added for lectures where no existing config exercises the right physics (L09: FD order sweep, L11: grid convergence).

3. **Output format per executable:** Each executable writes a different output file. Parsers in `star_helpers.py` already handle all four:
   - `bandStructure` → `output/eigenvalues.dat`
   - `gfactorCalculation` → `output/gfactor.dat`
   - `opticalProperties` → `output/absorption_TE.dat`, `output/absorption_TM.dat`
   - `topologicalAnalysis` → `output/topology_result.dat`

4. **Overlay plot format:** Each lecture script generates a matplotlib figure with two curves (code output, analytical reference) and saves to `docs/lecture/figures/lecture_XX_<name>.png`. Plots include axis labels, legends, tolerance bands where applicable, and a pass/fail indicator in the title.

5. **Tolerance tiers** (from `star_helpers.py` and verification ladder learnings):
   - `TOL_EXACT` (1e-12): k=0 eigenvalue self-check (Hamiltonian diagonal construction)
   - `TOL_ANALYTICAL` (1%): Roth g-factor, Kane effective mass, Bir-Pikus band shifts
   - `TOL_NUMERICAL` (5%): QW subbands, optical spectra, wire eigenvalues
   - `TOL_MODEL_LIMIT` (20-30%): 8-band g-factor vs experimental values (documented limitation, not a failure)

6. **Prerequisite fix:** `scripts/verify_landau_levels.py` reads from `band_results.dat` (doesn't exist) instead of `eigenvalues.dat`. Fixed in U13 as part of the refactoring.

7. **2-layer config pattern:** All QW configs in lecture scripts use the 2-layer pattern (barrier covers full domain, well overwrites center) to avoid the last-layer-wins pitfall documented in `docs/solutions/logic-errors/standard-star-benchmark-suite-logic-errors-2026-05-09.md`.

## Implementation Units

### U1: Lecture 00 — Quickstart (Worked Example)

**Physics:** Bulk GaAs band structure basics. No analytical assertions (R8).

**Executable:** `bandStructure`

**Config:** Reuse `tests/regression/configs/bulk_gaas_k0.cfg`

**Script:** `scripts/lecture_00_quickstart.py`

**Behavior:**
- Runs bulk GaAs k=0 config
- Parses eigenvalues, prints 8-band spectrum
- Generates band structure plot along kx direction (reuse `bulk_gaas_kx.cfg`)
- Labels bands by basis assignment (HH, LH, SO, CB)
- No numerical assertions — just demonstration
- Overlay plot shows code output with annotated band labels

**Test scenarios:**
- Script runs without error and produces a PNG figure
- Figure contains 8 bands with correct ascending ordering
- No assertion failures (methodological lecture)

**Files:** `scripts/lecture_00_quickstart.py`, `docs/lecture/figures/lecture_00_quickstart.png`

---

### U2: Lecture 01 — Bulk Band Structure

**Physics:** 8-band bulk Hamiltonian at k=0 and k≠0, effective mass extraction.

**Executable:** `bandStructure`

**Configs:** `bulk_gaas_k0.cfg`, `bulk_gaas_kx_dispersion.cfg`, `bulk_inas_kx_dispersion.cfg`

**Script:** `scripts/lecture_01_bulk.py`

**Behavior:**
- Runs GaAs k=0, asserts eigenvalues match `[-DeltaSO, -DeltaSO, 0, 0, 0, 0, Eg, Eg]` within `TOL_EXACT`
- Runs GaAs kx dispersion, extracts CB effective mass via single-point derivative `d2E_dk2 = 2*(E(k1)-E(0))/k1^2`
- Asserts GaAs CB mass matches Kane formula `m* = Eg/(EP+Eg)` within 10% (8-band higher-order mixing)
- Generates overlay: code E(k) vs parabolic Kane dispersion
- Runs InAs dispersion, asserts band gap matches Vurgaftman 0.417 eV within `TOL_EXACT`

**Analytical references:** GaAs Eg=1.519 eV, DeltaSO=0.341 eV (from `parameters.f90`); Kane mass from EP=28.8 eV; InAs Eg=0.417 eV.

**Tolerances:** k=0 eigenvalues: `TOL_EXACT`. Effective mass: 10%. Band gaps: `TOL_EXACT`.

**Test scenarios:**
- k=0 eigenvalues self-consistent with parameter values
- CB effective mass within 10% of Kane formula
- InAs band gap matches Vurgaftman
- Overlay plot shows parabolic vs code dispersion with visible non-parabolicity

**Files:** `scripts/lecture_01_bulk.py`, `docs/lecture/figures/lecture_01_bulk_dispersion.png`

---

### U3: Lecture 02 — Quantum Well

**Physics:** 8-band QW Hamiltonian, subband energies, confinement effect.

**Executable:** `bandStructure`

**Configs:** `qw_gaas_algaas.cfg`, `qw_gaas_algaas_double_qw.cfg`

**Script:** `scripts/lecture_02_qw.py`

**Behavior:**
- Runs GaAs/AlGaAs QW config, parses eigenvalues at k=0
- Verifies CB subband count matches expected (2 for this well width)
- Generates subband structure plot E(k_parallel)
- Runs double QW, verifies anticrossing splitting present (2 closely-spaced CB subbands)
- No clean analytical reference for 8-band QW subbands (Bastard formula has 5x discrepancy) — use regression comparison against known good output

**Tolerances:** Subband count: exact match. Subband energies: `TOL_NUMERICAL` (5%) vs regression reference.

**Test scenarios:**
- CB subband count is correct for well width
- Double QW shows anticrossing (2 closely-spaced levels)
- k_parallel dispersion plot generated successfully

**Files:** `scripts/lecture_02_qw.py`, `docs/lecture/figures/lecture_02_qw_subbands.png`

---

### U4: Lecture 03 — Wavefunctions

**Physics:** Eigenvector block decomposition, band character, probability density.

**Executable:** `bandStructure`

**Config:** `qw_gaas_algaas.cfg`

**Script:** `scripts/lecture_03_wavefunctions.py`

**Behavior:**
- Runs QW config, parses eigenvalues and eigenfunctions
- Decomposes CB ground state into 8-band basis: verifies CB character > 90%
- Verifies HH ground state has dominant HH character > 90%
- Generates wavefunction plot with per-band decomposition
- Verifies normalization: integral of |psi|^2 = 1 within `TOL_NORM` (1e-10)

**Tolerances:** Band character: >90% dominant. Normalization: `TOL_NORM` (1e-10).

**Test scenarios:**
- CB state has >90% CB character
- HH state has >90% HH character
- All states normalized
- Wavefunction plot shows confinement in well region

**Files:** `scripts/lecture_03_wavefunctions.py`, `docs/lecture/figures/lecture_03_wavefunctions.png`

---

### U5: Lecture 04 — Strain

**Physics:** Bir-Pikus strain Hamiltonian, band-edge shifts, HH-LH splitting.

**Executable:** `bandStructure`

**Configs:** `bulk_gaas_strained.cfg`, `qw_inas_gaas_strained.cfg`

**Script:** `scripts/lecture_04_strain.py`

**Behavior:**
- Runs strained bulk GaAs, compares eigenvalues at k=0 against unstrained
- Computes expected HH-LH splitting from Vurgaftman deformation potentials: `delta_E_HH-LH = 2*b*eps_xx` where b is the shear deformation potential
- Asserts computed splitting matches within `TOL_ANALYTICAL`
- Generates overlay: strained vs unstrained band structure with Bir-Pikus shifts annotated
- Runs strained InAs/GaAs QW, verifies strain-induced HH-LH crossover is present

**Analytical references:** GaAs b=-2.0 eV, a_v=1.16 eV, a_c=-7.17 eV (Vurgaftman 2001 Table XIV). Biaxial strain: eps_xx = (a_sub - a_layer) / a_layer.

**Tolerances:** Band-edge shifts: `TOL_ANALYTICAL` (1%). HH-LH splitting: `TOL_ANALYTICAL`.

**Test scenarios:**
- Strained eigenvalues shift correctly vs unstrained
- HH-LH splitting matches Bir-Pikus formula
- Strained QW shows correct band ordering

**Files:** `scripts/lecture_04_strain.py`, `docs/lecture/figures/lecture_04_strain_shifts.png`

---

### U6: Lecture 05 — g-Factor

**Physics:** Landau g-factor via Lowdin partitioning, Roth formula, Landau levels.

**Executable:** `gfactorCalculation`, `bandStructure` (Landau mode)

**Configs:** `gfactor_bulk_gaas_cb.cfg`, `gfactor_bulk_gaasw_cb.cfg`, `gfactor_bulk_insb_cb.cfg`, `landau_bulk_InAs.cfg`

**Script:** `scripts/lecture_05_gfactor.py`

**Behavior:**
- Runs bulk GaAs g-factor, asserts CB g matches Roth formula value (-0.315) within `TOL_ANALYTICAL`
- Runs bulk InSb g-factor, asserts extreme g-value regime (|g| > 40)
- Runs GaAsW g-factor, asserts Winkler parameter set result
- Generates overlay: code g-factors vs Roth formula for multiple materials
- Runs Landau mode InAs, verifies Landau level spacing matches `hbar*omega_c/2` within `TOL_NUMERICAL`
- Documents 20-30% g-factor shortfall as known model limitation (R9)

**Analytical references:** Roth formula GaAs CB g = -0.315, InAsW g = -14.858 (from past verification). Landau level: `E_n = E_C + (n+1/2)*hbar*omega_c`.

**Tolerances:** Roth g-factor: `TOL_ANALYTICAL` (1%). Landau levels: `TOL_NUMERICAL` (5%). Experimental comparison: `TOL_MODEL_LIMIT` (30%, documented limitation).

**Test scenarios:**
- GaAs CB g within 1% of Roth value
- InSb g in extreme regime
- Landau level spacing matches cyclotron frequency
- Known limitation documented in script output

**Files:** `scripts/lecture_05_gfactor.py`, `docs/lecture/figures/lecture_05_gfactor_comparison.png`

---

### U7: Lecture 06 — Optical Properties

**Physics:** Absorption, gain, spontaneous emission, ISBT, TE/TM polarization.

**Executable:** `opticalProperties`

**Configs:** `bulk_gaas_optics.cfg`, `qw_gaas_algaas_optics.cfg`, `qw_gaas_algaas_isbt.cfg`

**Script:** `scripts/lecture_06_optical.py`

**Behavior:**
- Runs bulk GaAs optics, verifies absorption onset at band gap energy
- Runs QW absorption, verifies TE/TM polarization ratio
- Runs ISBT config, verifies ISBT peak at expected transition energy
- Generates overlay: code spectra with absorption edge annotated
- No closed-form analytical spectrum — validate onset energies and peak positions

**Tolerances:** Absorption edge: `TOL_NUMERICAL` (5%). Polarization ratio: `TOL_NUMERICAL`. ISBT peak: `TOL_NUMERICAL`.

**Test scenarios:**
- Bulk absorption onset at Eg
- QW TE polarization dominant
- ISBT peak at expected energy
- Overlay plot shows spectra with annotated features

**Files:** `scripts/lecture_06_optical.py`, `docs/lecture/figures/lecture_06_absorption.png`

---

### U8: Lecture 07 — Self-Consistent SP

**Physics:** Schrödinger-Poisson iteration, charge neutrality, convergence, DIIS.

**Executable:** `bandStructure`

**Configs:** `sc_bulk_gaas_doped.cfg`, `sc_gaas_alas_qw.cfg`

**Script:** `scripts/lecture_07_scsp.py`

**Behavior:**
- Runs SC bulk GaAs doped, verifies Fermi level convergence
- Runs SC QW, verifies charge neutrality within tolerance
- Generates convergence plot: energy vs iteration number
- Verifies DIIS convergence faster than linear mixing (fewer iterations)
- No closed-form for self-consistent potential — validate convergence behavior and self-consistency

**Tolerances:** Charge neutrality: `TOL_NUMERICAL` (5%). Convergence: iteration count decrease vs linear mixing.

**Test scenarios:**
- SC loop converges (energy change < tolerance)
- Charge neutrality achieved
- DIIS faster than linear mixing
- Convergence plot generated

**Files:** `scripts/lecture_07_scsp.py`, `docs/lecture/figures/lecture_07_convergence.png`

---

### U9: Lecture 08 — Quantum Wire

**Physics:** 2D confinement, CSR sparse assembly, wire subbands.

**Executable:** `bandStructure`

**Configs:** `wire_gaas_31x31.cfg`, `wire_inas_rectangle.cfg`

**Script:** `scripts/lecture_08_wire.py`

**Behavior:**
- Runs GaAs wire, verifies subband count matches expectation for grid size
- Runs InAs wire, verifies wire subband spacing decreases with increasing cross-section
- Generates wire subband structure plot (E vs kz)
- Verifies dense-sparse consistency (same eigenvalues from both solvers)

**Tolerances:** Subband count: exact. Dense-sparse eigenvalue agreement: `TOL_EXACT` (1e-10). Subband spacing trend: qualitative.

**Test scenarios:**
- Wire subbands computed successfully
- Dense-sparse consistency within machine precision
- Subband count consistent with confinement geometry

**Files:** `scripts/lecture_08_wire.py`, `docs/lecture/figures/lecture_08_wire_subbands.png`

---

### U10: Lecture 09 — Numerical Methods

**Physics:** FD stencils, convergence order, Richardson extrapolation.

**Executable:** `bandStructure`

**Config:** `bulk_gaas_k0.cfg` (with varying FDorder 2, 4, 6, 8)

**Script:** `scripts/lecture_09_numerical.py`

**Behavior:**
- Runs bulk GaAs k=0 with FD orders 2, 4, 6, 8
- Verifies k=0 eigenvalues are order-independent (analytical by construction)
- Runs kx dispersion at each FD order, extracts convergence rate
- Asserts observed convergence rate matches theoretical order within tolerance
- Generates overlay: eigenvalue convergence vs FD order
- **New config needed:** `bulk_gaas_kx_fd_order_sweep.cfg` (or parameterized within script)

**Tolerances:** k=0 eigenvalues: `TOL_EXACT` (order-independent). Convergence rate: within 0.5 of theoretical order.

**Test scenarios:**
- k=0 eigenvalues identical across FD orders
- Convergence rate increases with FD order
- Overlay plot shows convergence behavior

**Files:** `scripts/lecture_09_numerical.py`, `docs/lecture/figures/lecture_09_convergence.png`, `tests/regression/configs/bulk_gaas_kx_fd2.cfg` (new, if needed)

---

### U11: Lecture 10 — QCSE

**Physics:** Quantum-confined Stark effect, field-dependent subband shift, scattering.

**Executable:** `bandStructure`

**Configs:** `sc_qcse_gaas_algaas.cfg`, `sc_qcse_gaas_algaas_ef.cfg`

**Script:** `scripts/lecture_10_qcse.py`

**Behavior:**
- Runs QW without field (reference), then with field (-70 kV/cm)
- Computes Stark shift = E_CB(field) - E_CB(0)
- Verifies Stark shift is negative (red shift) and within expected magnitude
- Generates overlay: band structure with and without field
- Verifies scattering lifetimes are physically reasonable (ps scale)

**Tolerances:** Stark shift direction: exact (must be negative). Stark shift magnitude: `TOL_NUMERICAL` (5%) vs reference data.

**Test scenarios:**
- Field causes red shift of CB subband
- Stark shift magnitude is physically reasonable
- Overlay plot shows field-induced band bending

**Files:** `scripts/lecture_10_qcse.py`, `docs/lecture/figures/lecture_10_stark_shift.png`

---

### U12: Lecture 11 — Convergence

**Physics:** Grid spacing convergence, FD order convergence, Richardson extrapolation.

**Executable:** `bandStructure`

**Configs:** `bulk_gaas_kx_dispersion.cfg` (varying grid), `qw_gaas_algaas.cfg` (varying grid)

**Script:** `scripts/lecture_11_convergence.py`

**Behavior:**
- Runs bulk dispersion at multiple grid spacings, plots eigenvalue convergence
- Runs QW at multiple grid spacings, plots subband convergence
- Verifies convergence rate matches FD order (Richardson extrapolation)
- Generates overlay: eigenvalue vs grid spacing with convergence order annotation
- **New config needed:** may need configs with different grid sizes (or parameterize within script)

**Tolerances:** Convergence rate: matches theoretical FD order within 0.5.

**Test scenarios:**
- Eigenvalues converge with decreasing grid spacing
- Convergence order matches FD stencil order
- Overlay plot shows convergence curves

**Files:** `scripts/lecture_11_convergence.py`, `docs/lecture/figures/lecture_11_grid_convergence.png`

---

### U13: Lecture 13 — Topological Superconductivity

**Physics:** Chern number (QWZ model), Z2 invariant (BHZ), BdG Majorana modes, Rashba phase boundary.

**Executable:** `topologicalAnalysis`

**Configs:** `topology_qwz_chern_u-0.8.cfg`, `topology_qwz_chern_u0.5.cfg`, `topology_qwz_chern_u2.5.cfg`, `topology_bhz_z2_trivial.cfg`, `topology_bhz_z2_topological.cfg`, `topology_qw_bdg.cfg`, `topology_rashba_phase.cfg`

**Script:** `scripts/lecture_13_topological.py`

**Behavior (absorbs R7 — all 4 existing verify scripts):**
- **QWZ Chern (from `verify_qwz_chern.py`):** Runs 3 QWZ configs (u=-0.8, 0.5, 2.5), asserts Chern numbers match theoretical predictions (C=+1, 0, -1). Generates Berry curvature heatmap.
- **BHZ Z2 (from `verify_bhz_z2.py`):** Runs trivial and topological configs, asserts Z2=0 and Z2=1. Generates phase diagram.
- **BdG Majorana (from `sweep_rashba_bdg.py`):** Runs Rashba phase sweep, verifies topological phase boundary. Generates phase diagram.
- **Landau levels (from `verify_landau_levels.py`, after fixing `band_results.dat` → `eigenvalues.dat` bug):** Runs Landau config, verifies level spacing. Plots Zeeman fan diagram.
- Overlay plots for each sub-section of lecture 13

**After refactoring, delete:** `scripts/verify_qwz_chern.py`, `scripts/verify_bhz_z2.py`, `scripts/verify_landau_levels.py`, `scripts/sweep_rashba_bdg.py`. Update `scripts/generate_all_figures.py` to call `lecture_13_topological.py` instead of the deleted scripts.

**Tolerances:** Chern number: exact integer match. Z2: exact integer match. BdG gap: within 1% of 2*Delta. Phase boundary: within 5% of analytical prediction.

**Test scenarios:**
- QWZ Chern numbers match theory for all 3 u values
- BHZ Z2 invariant is 0 (trivial) and 1 (topological)
- BdG gap matches Nambu convention (2*Delta)
- Rashba phase boundary at expected B-field
- All 4 legacy scripts fully absorbed
- Berry curvature heatmap, Z2 phase diagram, Majorana phase diagram generated

**Files:** `scripts/lecture_13_topological.py`, `docs/lecture/figures/lecture_13_*.png` (multiple), delete `scripts/verify_qwz_chern.py`, `scripts/verify_bhz_z2.py`, `scripts/verify_landau_levels.py`, `scripts/sweep_rashba_bdg.py`

---

### U14: Lecture 12 — Extending the Code (Worked Example)

**Physics:** Code architecture, adding new materials, extending the code. No analytical assertions (R8).

**Executable:** `bandStructure`

**Config:** `bulk_gaas_k0.cfg` (demonstrative)

**Script:** `scripts/lecture_12_extending.py`

**Behavior:**
- Runs bulk GaAs config as a worked example
- Demonstrates how to read and interpret output files
- Shows the relationship between input parameters and output eigenvalues
- Generates a simple demonstration plot
- No numerical assertions — pedagogical

**Test scenarios:**
- Script runs without error
- Output plot generated
- No assertion failures

**Files:** `scripts/lecture_12_extending.py`, `docs/lecture/figures/lecture_12_example.png`

---

### U15: Makefile Targets

**Files:** `Makefile`

**Behavior:**
- Add `lecture-00` through `lecture-13` targets, each running `python3 scripts/lecture_XX_<name>.py`
- Add `lectures` target that runs all 14 sequentially
- Each target builds the required executable first (dependency on `all`)
- Add `lectures` to `.PHONY`

**Pattern:**
```makefile
lecture-%: all
	python3 scripts/lecture_$*_*.py

lectures: all
	@for i in 00 01 02 03 04 05 06 07 08 09 10 11 12 13; do \
		echo "=== Lecture $$i ==="; \
		python3 scripts/lecture_$${i}_*.py || exit 1; \
	done

.PHONY: lecture-% lectures
```

Note: The glob pattern `lecture_$*_*.py` may need explicit mapping. Prefer explicit targets:

```makefile
lecture-00: all; python3 scripts/lecture_00_quickstart.py
lecture-01: all; python3 scripts/lecture_01_bulk.py
# ... etc for all 14
lectures: lecture-00 lecture-01 ... lecture-13
```

---

### U16: Lecture Doc Updates

**Files:** All 14 lecture docs in `docs/lecture/`

**Behavior:**
- For each lecture doc (U1-U14), add a section or update existing sections to embed the generated overlay plot
- Add code-output anchors: statements like "Running `bulk_gaas_k0.cfg` produces CB eigenvalue 1.519 eV at k=0" connecting derivations to code results
- Use relative image references: `![Overlay plot](figures/lecture_XX_name.png)`
- Add a "Verification" subsection to each lecture doc referencing the corresponding lecture-pair script

**Test scenarios:**
- Each lecture doc contains at least one image reference to a generated figure
- Each lecture doc contains at least one code-output anchor
- Image references use relative paths that resolve correctly

## Dependencies and Sequencing

```
Phase A — Infrastructure + Template
  U15 (Makefile) ─────────────────────────────┐
  U2  (L01 bulk) ← template for analytical    │
  U1  (L00 quickstart) ← template for worked  │
                                                │
Phase B — Core Physics Lectures                │
  U6  (L05 g-factor)                          │
  U3  (L02 QW)                                │
  U4  (L03 wavefunctions)                     │
  U5  (L04 strain)                            │
  U7  (L06 optical)                           │
  U8  (L07 SC-SP)                             │
                                                │
Phase C — Topological + Refactoring            │
  U13 (L13 topological) ← absorbs 4 scripts   │
                                                │
Phase D — Wire, Numerical, Specialized         │
  U9  (L08 wire)                              │
  U10 (L09 numerical)                          │
  U11 (L10 QCSE)                              │
  U12 (L11 convergence)                       │
  U14 (L12 extending)                          │
                                                │
Phase E — Documentation                        │
  U16 (lecture doc updates) ← after all U1-U14│
```

Within each phase, units can be implemented in any order. The phases enforce that:
- Phase A establishes the Makefile and script templates before Phase B
- Phase C (topological refactoring) is independent of Phase B
- Phase D can run in parallel with Phase C
- Phase E runs after all scripts are complete

## Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|-----------|
| 8-band g-factor 20-30% shortfall triggers false failures | Medium | Use R9 documented limitation flag; compare against Roth formula (code self-consistent) not experimental values |
| QW subband energies lack clean analytical reference | Medium | Use regression comparison and subband count checks; document that Bastard formula is inadequate for 8-band |
| Optical spectra lack closed-form validation | Low | Validate onset energies and peak positions; no spectrum shape comparison |
| `verify_landau_levels.py` bug (wrong filename) blocks L13 refactoring | Low | Fix as prerequisite in U13 |
| New configs needed for L09, L11 (FD order/grid sweep) | Low | Script can parameterize existing configs by modifying FDorder/grid at runtime |
| `star_helpers.py` import path breaks from `scripts/` directory | Low | Use `sys.path.insert(0, ...)` with repo-relative path; test in CI |
| 14 scripts take >30 min to run sequentially | Low | `make lectures` is manual-only per scope; individual `make lecture-XX` is fast |

## Success Criteria

- Running `make lectures` produces 14 passing pairs with overlay plots in `docs/lecture/figures/`
- Each lecture doc has at least one code-output anchor connecting derivation to computed value
- The 4 `scripts/verify_*.py` files are deleted after absorption into `scripts/lecture_13_topological.py`
- `scripts/generate_all_figures.py` continues to work (updated to reference lecture_13 script)
- Each lecture-pair script is independently runnable via `python3 scripts/lecture_XX_*.py`
