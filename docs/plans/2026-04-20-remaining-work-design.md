# Remaining Work Completion Design

**Goal:** Close all open gaps identified in the plan audit across 4 phases: cleanup, missing figures, wire gaps, and SC benchmarks.

**Architecture:** 4 phases ordered by dependency. Phase 1 is pure cleanup. Phase 2 is Python-only. Phase 3 requires Fortran development. Phase 4 is config + validation work.

**Tech Stack:** Fortran 90 (Phase 3), Python + matplotlib (Phase 2), shell configs (Phase 4).

---

## Phase 1: Cleanup Sweep

Delete stale files and fix cosmetic code items. No new functionality.

1. **Delete stale root files:** `rm modernization_plan.md plan.md checklist.md` — superseded by `docs/plans/` and lecture notes.
2. **Delete sub-READMEs:** `rm src/math/README.md src/io/README.md src/physics/README.md src/apps/README.md` — superseded by lecture docs and CLAUDE.md.
3. **Fix `2.00231` comment** in `src/apps/main_gfactor.f90` lines ~399/408/417: change `!+ 2.00231` to `!+ free-electron g-factor (included via sigma tensor)`.
4. **Add `dz=0` guard** in `pMatrixEleCalc` (`src/physics/gfactor_functions.f90`): add early return or error, matching existing guard in `sigmaElem`.

Commit: `chore: cleanup stale files and minor code fixes`

---

## Phase 2: Missing QW Tutorial Figures

4 figures from the QW tutorials design that were not yet generated. Python-driven parametric sweeps calling the Fortran executable.

### Figure 1: `exciton_bohr_vs_width.png`

**What:** Exciton Bohr radius vs well width (reproduction of Harrison Fig 6.5).
**Data source:** Run `bandStructure` for GaAs/AlGaAs QW with widths 40-200 AA in ~10 steps. Feed eigenvalues to `exciton.f90` (variational binding energy -> Bohr radius via `a_B = sqrt(hbar^2 / (2 * mu * E_b))`).
**Approach:** Python loop generates configs with varying `numLayers` (well width), calls executable, parses eigenvalues, computes exciton binding energy and Bohr radius using the same formulas as `exciton.f90` (or directly from exciton output if available).

### Figure 2: `scattering_lifetime_vs_field.png`

**What:** LO-phonon scattering lifetime vs applied electric field (Ferreira & Bastard PRB 1989, Fig 5).
**Data source:** Run `bandStructure` with QCSE configs at multiple field values. Feed wavefunctions to `scattering.f90` Frohlich rate computation.
**Approach:** Python loop over `ExternalField` values (0-50 kV/cm), parse scattering rate output, plot 1/rate vs field.

### Figure 3: `double_qw_anticrossing.png`

**What:** Subband anticrossing in an asymmetric double quantum well (Ferreira & Bastard PRB 1989, Fig 7).
**Data source:** Run `bandStructure` with a coupled double-QW config. Vary barrier width or asymmetry to show anticrossing.
**Approach:** Single config with double QW geometry, plot E vs k showing the avoided crossing between the two lowest CB subbands.

### Figure 4: `absorption_excitonic_TE.png`

**What:** Full excitonic absorption spectrum with TE polarization, comparing with and without excitonic correction.
**Data source:** Run optical spectrum computation (absorption) + exciton correction for a GaAs/AlGaAs QW.
**Approach:** Parse `absorption_spectrum.dat` output and overlay Sommerfeld-enhanced absorption from `exciton.f90` output.

### Implementation notes

- All figures go into `scripts/plotting/generate_all_figures.py` as new functions registered in `ALL_FIGURES`.
- Any new regression configs go into `tests/regression/configs/`.
- May need to check if existing Fortran output formats support the needed data; if not, minimal Fortran output additions are acceptable.

---

## Phase 3: Wire Gaps

Three items requiring Fortran development for the quantum wire path.

### 3a: Wire `parts.dat` (8-band decomposition)

**What:** Decompose each wire eigenstate into contributions from the 8 basis bands (HH, LH, SO, CB).
**Current state:** QW version exists (`writeParts` in `outputFunctions.f90`) which projects eigenvectors onto band subspaces and writes a multi-block file.
**Implementation:** Create `writeParts2d` that:
- Takes 2D eigenvector array (shape: `[8*nx*ny]` or `[8, nx, ny]`)
- Projects onto each of the 8 band subspaces
- Writes per-band densities as 2D arrays to `parts.dat`
- Follows the same multi-block format as the 1D version

### 3b: Wire SC diagnostics

**What:** Write `sc_phi.dat` and `sc_charge.dat` per SC iteration for wire mode.
**Current state:** 1D QW SC already writes these files in `sc_loop.f90`.
**Implementation:** Add analogous file output in the wire SC loop path (`self_consistent_loop_wire` in `sc_loop.f90`):
- `sc_phi.dat`: electrostatic potential profile at each iteration
- `sc_charge.dat`: charge density at each iteration
- Follow the same format as the 1D versions

### 3c: Tier 3 wire benchmarks

**What:** Two regression tests against published data.

**Benchmark 1: Stier & Bimberg 1997 core-shell wire**
- Material: InAs/GaAs core-shell nanowire
- Validates: strain field + band structure
- Config: `wire_inas_gaas_core_shell.cfg`
- Golden data: eigenvalues.dat from reference computation

**Benchmark 2: InSb wire g-factor**
- Material: InSb nanowire (extreme g-factor ~ -50)
- Validates: wire g-factor calculation at k=0
- Config: `wire_insb_gfactor.cfg`
- Golden data: gfactor.dat from reference computation

Both require:
1. New regression config files
2. Running the code to generate golden reference data
3. Shell test scripts in `tests/integration/`
4. Registration in `tests/regression/CMakeLists.txt` (or equivalent)

---

## Phase 4: SC Benchmarks + Validation

### 4a: Bulk doped GaAs regression test

**What:** Self-consistent bulk n-GaAs with doping, validating Fermi level against Sze (2007) carrier statistics.
**Config:** `sc_bulk_gaas_doped.cfg` — bulk mode with `SC=1`, n-type doping, temperature sweep.
**Validation:** Fermi level position from SC loop matches analytical Thomas-Fermi approximation within tolerance.

### 4b: InAs/AlSb QW regression test

**What:** Narrow-gap InAs/AlSb QW with self-consistent calculation, validating against Pfeffer (1999).
**Config:** `sc_qw_inas_alsb.cfg` — QW mode with InAs/AlSb materials, doping.
**Validation:** Subband energies and Fermi level match published values within tolerance.

### 4c: Validation documentation

**What:** Add comparison tables in relevant lecture chapters documenting agreement with published references.
- Ch06 (Optical): Comparison with Dumitras PRB 2002 absorption data, Harrison Figs 6.4-6.5 exciton binding energies, Ferreira & Bastard PRB 1989 scattering rates
- Ch05 (g-factor): Wire g-factor vs Stier & Bimberg 1997
- Ch07 (SC): Bulk carrier statistics vs Sze 2007, QW subbands vs Pfeffer 1999
- Format: Tables with "Published", "Computed", "Relative Error" columns

---

## Scope Exclusions

- **Drift-diffusion solver:** Explicitly deferred. Document as known limitation in Ch12.
- **GPU eigensolver (linalg Phase 2):** Deferred pending hardware.
- **I9 simulation_config type bloat:** Deferred from PR review, low priority.
