**Status**: COMPLETE (2026-07-05)

# PRD: Publishable Validation Benchmarks (Phase 21)

## Problem Statement

The solver has 111 tests across 20 phases of validation — verification ladder (8 rungs), standard-star benchmarks (S1-S7), cross-code kdotpy comparison (12 tests), Richardson convergence, TRK sum rule, and Hermiticity checks. Yet three physics modules that actively compute observables have no quantitative benchmark against an independent reference:

1. **ISBT (intersubband transitions):** `optical_spectra.f90` computes ISBT z-dipole matrix elements and oscillator strengths for QW systems. Two independent implementations (direct z-dipole and commutator velocity) agree to 0.007% — but nobody has checked whether either gives the *right answer* against an analytical formula. The only test is a smoke check that the absorption peak falls somewhere in a 10× energy window (30–300 meV).

2. **SC Schrödinger-Poisson:** `sc_loop.f90` iterates the self-consistent loop with DIIS acceleration. Convergence is validated (Fermi level Richardson convergence, |ΔΦ| decreases). But charge neutrality — an exact conservation law where ∫n_e dz must equal ∫N_D dz — is only a soft warning in a lecture script (5% WARN, 10% FAIL), not a hard CI-gated check. No test compares the self-consistent Fermi level, subband shifts, or potential profile against an analytical or published reference.

3. **Exciton binding energy:** `exciton.f90` implements the Bastard variational method. A convergence test compares against Miller et al. 1985 at 30% tolerance. The test is registered under the convergence label but not in the validation coverage matrix (`validation_universe.yml`). No test compares against the Bastard 1982 original quantitative values.

These gaps must be closed before a J. Comput. Phys. methods paper can claim comprehensive validation.

## Solution

Implement three publishable-quality benchmark suites, each comparing solver output against an independent analytical or published reference, with figures and tables suitable for inclusion in a methods paper:

1. **ISBT oscillator strength benchmark** — compare z-dipole and oscillator strength against the infinite square well analytical result for two type-I material systems (GaAs/AlGaAs and InGaAs/InAlAs) at three well widths each. The deviation from the infinite-well formula quantifies band-mixing and non-parabolicity corrections.

2. **SC Schrödinger-Poisson benchmark suite** — hard charge neutrality check (conservation law), Fermi level vs analytical parabolic approximation, CB1 subband shift vs flat-band with Bastard perturbative estimate, and qualitative potential profile validation.

3. **Exciton binding energy benchmark** — tighten Miller 1985 tolerance from 30% to ≤15%, add Bastard 1982 quantitative comparison at specific well widths, register in validation_universe.yml as required cell.

All three are promoted from aspirational/unregistered to required in the validation coverage matrix, with COVERAGE annotations and ctest registration.

## User Stories

### ISBT Oscillator Strength Benchmark

1. As a physicist running opticalProperties on a GaAs/AlGaAs QW, I want the computed e1→e2 ISBT z-dipole compared against the infinite-well analytical value `|⟨1|z|2⟩| = 16L/(9π²)`, so that I know the intersubband physics is quantitatively correct
2. As a physicist running opticalProperties on a GaAs/AlGaAs QW, I want the computed oscillator strength compared against the analytical `f₁₂ = (2m₀E₁₂/ℏ²)|⟨1|z|2⟩|²`, so that I know the dipole-velocity conversion is self-consistent
3. As a physicist studying ISBT across well widths, I want the benchmark to cover 50, 100, and 200 Å wells, so that I can see how confinement affects the 8-band vs infinite-well deviation
4. As a physicist studying material dependence of ISBT, I want the benchmark to cover GaAs/AlGaAs (wide gap) and InGaAs/InAlAs (narrow gap, QCL-relevant), so that I can see how non-parabolicity increases with band mixing
5. As a physicist reviewing a methods paper, I want a benchmark table showing |⟨1|z|2⟩|_8band / |⟨1|z|2⟩|_inf vs well width for two materials, so that I can evaluate whether the deviation is physically sensible (larger for narrower wells, larger for narrower-gap materials)
6. As a physicist running the InGaAs/InAlAs benchmark, I want the config to use the existing Ga₀.₄₇In₀.₅₃AsW / Al₀.₄₇In₀.₅₃AsW parameter sets (InP lattice-matched), so that no new material parameter entries are needed
7. As a developer modifying the optical accumulation pipeline, I want the ISBT benchmark to catch regressions in z-dipole or oscillator strength computation, so that I don't silently break ISBT physics
8. As a developer, I want the ISBT cell promoted from aspirational to required in validation_universe.yml, so that CI gates on ISBT correctness

### SC Schrödinger-Poisson Benchmark Suite

9. As a physicist running the SC loop on a doped GaAs/AlAs QW, I want charge neutrality enforced as a hard pass/fail check (`|∫n_e dz - ∫N_D dz| / ∫N_D < 1%`), so that I know the Schrödinger-Poisson iteration conserves charge
10. As a physicist running the SC loop on bulk n-GaAs, I want the self-consistent Fermi level compared against the analytical parabolic result `E_F = E_C + (ℏ²/2m*)(3π²n)^(2/3)`, so that I know the Fermi bisection converges to the correct chemical potential
11. As a physicist running the SC loop on a doped QW, I want the CB1 subband shift (SC eigenvalue vs flat-band) compared against the Bastard perturbative estimate, so that I know the self-consistent potential produces the expected confinement modification
12. As a physicist running the SC loop on a doped QW, I want the Hartree potential profile verified to have a V-shaped dip in the well and flat barriers, so that I know the Poisson solver produces a physically reasonable electrostatic potential
13. As a physicist running the SC loop on a doped QW, I want the potential energy scale verified (Hartree energy should be on the order of e²n^(1/3)/(4πε₀ε_r)), so that the magnitude is physically correct
14. As a developer modifying the Poisson solver or DIIS mixer, I want the SC benchmark to catch regressions in charge conservation or potential shape, so that I don't silently degrade SC accuracy
15. As a developer, I want `charge_neutrality`, `subband_shift`, and `potential_profile` registered as observable cells in validation_universe.yml, so that the coverage matrix tracks SC physics

### Exciton Binding Energy Benchmark

16. As a physicist running the exciton solver on a GaAs/AlGaAs QW, I want the binding energy compared against Miller et al. 1985 at ≤15% tolerance (not 30%), so that the benchmark provides meaningful quantitative validation
17. As a physicist running the exciton solver on a GaAs/AlGaAs QW, I want the binding energy compared against Bastard 1982 specific quantitative values at multiple well widths, so that the variational method is validated against its original reference
18. As a physicist studying exciton trends, I want the benchmark to cover the same well widths as lecture_14 (30, 50, 80, 100, 150, 200 Å), so that the width-dependent trend can be verified
19. As a physicist reviewing a methods paper, I want a figure showing Eb vs well width with both 8-band numerical results and the Bastard analytical trend overlaid, so that I can see the agreement visually
20. As a developer modifying the variational solver or Coulomb integral, I want the tightened tolerance to catch regressions, so that I don't silently worsen the exciton energy
21. As a developer, I want `exciton_Eb` registered as a required observable cell in validation_universe.yml (currently absent), so that the coverage matrix tracks exciton physics

### Developer Stories

22. As a developer, I want all new tests registered under appropriate ctest labels (verification for ISBT and SC, convergence for exciton update), so that they integrate with the existing 111-test CI structure
23. As a developer, I want all new verification scripts to include COVERAGE annotations, so that the coverage matrix (`ctest -L coverage`) tracks the new benchmark cells
24. As a developer, I want the ISBT golden data captured and registered in the regression test suite, so that future code changes that silently alter ISBT results are detected
25. As a developer, I want the benchmark figures generated automatically by the verification scripts (matching the lecture-script pattern), so that the paper figures are reproducible from CI

## Implementation Decisions

### ISBT Benchmark: infinite-well analytical comparison

**Analytical reference:** For an infinite square well of width L, the e1→e2 transition has:

- z-dipole: `|⟨1|z|2⟩| = 16L / (9π²)` (from the standard particle-in-a-box overlap integral)
- oscillator strength: `f₁₂ = (2m₀E₁₂ / ℏ²) |⟨1|z|2⟩|²` where `E₁₂ = 3π²ℏ² / (2m*L²)`
- combined: `f₁₂ = 0.9557` (universal constant for all infinite wells, independent of mass and width)

The 8-band result will deviate from these due to: (a) finite barrier height (wavefunction penetration into barriers), (b) band mixing (the 8 bands couple the CB states to VB states), (c) non-parabolicity (especially in narrow-gap InGaAs/InAlAs). The deviation ratio `|z₁₂|_8band / |z₁₂|_inf` is the key benchmark quantity — it should be < 1 (finite barriers reduce the dipole) and should decrease further for narrower wells and narrower-gap materials.

**Material systems (no new parameter entries needed):**

| System | Well | Barrier | Config basis |
|--------|------|---------|--------------|
| GaAs/AlGaAs | GaAs | Al30Ga70As | `qw_gaas_algaas_isbt.toml` (exists) |
| InGaAs/InAlAs | Ga47In53AsW | Al47In53AsW | `qw_ingaas_algaas_strained_optics.toml` (flip ISBT flag) |

**Config generation:** 6 TOML configs (3 widths × 2 materials). The well width is parameterized by the `[[material]]` layer extent. Reuse existing configs as templates.

**Verification script:** `verify_isbt_benchmark.py` following the pattern of `verify_8band_rung7_gfactor.py` and `verify_8band_rung8_optical.py` — run the opticalProperties executable, parse `output/isbt_transitions.dat` (z-dipole table) and `output/absorption_ISBT.dat` (spectrum), compare e1→e2 z-dipole and oscillator strength against infinite-well analytical.

**Key output:** Benchmark table with 6 rows (material × width), columns: well width, |z₁₂|_inf, |z₁₂|_8band, ratio, f₁₂_8band, E₁₂. Figure: ratio vs well width for both materials on the same axes.

### SC Benchmark Suite: four independent checks

**System:** GaAs/AlAs QW, uniformly n-doped at 5×10¹⁸ cm⁻³ in the well, matching `test_sc_convergence.py` config (51 grid points, FDorder=2, 41 k_∥ points).

**Charge neutrality check:** After SC convergence, integrate `sc_charge.dat` n_e(z) column via trapezoidal rule. Compare against `∫N_D dz` (doping × well thickness). Hard tolerance: 1% relative. This is an exact conservation law — any solver that doesn't conserve charge to < 1% has a bug.

**Fermi level check:** For bulk n-GaAs at doping n = 5×10¹⁸ cm⁻³, the analytical parabolic Fermi level is `E_F = E_C + (ℏ²/2m*)(3π²n)^(2/3)` where m* is the GaAs CB effective mass (0.067m₀). For the QW case, compare against the bulk value as a reference point — the QW Fermi level will differ due to subband quantization, but should be within a factor of 2. Use the bulk value as a sanity bound.

**Subband shift check:** Run bandStructure twice — once with SC enabled (self-consistent potential), once flat-band (no SC). Extract CB1 eigenvalue at k_∥=0 from both runs. The shift `ΔE₁ = E₁(SC) - E₁(flat)` should match the Bastard first-order perturbative estimate `ΔE ≈ -e²n_d/(2ε₀ε_r) × <ψ₁|z²|ψ₁>` to within ~20% (the Bastard estimate is approximate for finite barriers).

**Potential profile check:** Parse `output/sc_potential.dat` or equivalent. Verify: (a) potential is parabolic-like in the well (not flat, not linear), (b) potential is flat in the barriers (screened), (c) total potential swing `ΔV = V(center) - V(barrier)` is on the order of the expected Hartree energy `e²N_D × L_well² / (2ε₀ε_r)`.

**Verification script:** `verify_sc_benchmark.py` following the pattern of `test_sc_convergence.py` but with hard quantitative assertions instead of Richardson convergence only.

### Exciton Benchmark: tighten and extend

**Tolerance tightening:** The existing `test_exciton_convergence.py` compares Richardson-extrapolated Eb against Miller 1985 at 30%. The convergence data suggests the actual agreement is much better. Tighten to ≤15% (the Richardson extrapolant at 401 points should easily achieve this for a 100 Å GaAs/AlGaAs QW).

**Bastard 1982 comparison:** Add a second quantitative check comparing Eb against Bastard's original variational result. Bastard (PRB 1982) tabulates Eb for GaAs QWs at several well widths with Al₀.₃Ga₀.₇As barriers. Use these values as a second reference column. Expected tolerance: ≤20% (the Bastard model uses a simpler effective-mass approach, so some deviation from the 8-band result is expected).

**validation_universe.yml registration:** Add `exciton_Eb` as a new observable in the metadata.observables list. Add required cells:

- `exciton_Eb` / QW / GaAs/AlGaAs / required (reference: Miller 1985, Bastard 1982)

**Update `test_exciton_convergence.py`:** Add Miller tolerance constant (currently hardcoded), add Bastard reference values, add assertion. Keep Richardson convergence as the primary check.

### Config and infrastructure

**New TOML configs (6 total):**

| Config | Material | Width | Based on |
|--------|----------|-------|----------|
| `isbt_gaas_algaas_w50.toml` | GaAs/Al30Ga70As | 50 Å | `qw_gaas_algaas_isbt.toml` |
| `isbt_gaas_algaas_w100.toml` | GaAs/Al30Ga70As | 100 Å | (same, already exists as ISBT config for ~100 Å) |
| `isbt_gaas_algaas_w200.toml` | GaAs/Al30Ga70As | 200 Å | `qw_gaas_algaas_isbt.toml` |
| `isbt_ingaas_inalas_w50.toml` | Ga47In53AsW/Al47In53AsW | 50 Å | `qw_ingaas_algaas_strained_optics.toml` |
| `isbt_ingaas_inalas_w100.toml` | Ga47In53AsW/Al47In53AsW | 100 Å | (same) |
| `isbt_ingaas_inalas_w200.toml` | Ga47In53AsW/Al47In53AsW | 200 Å | (same) |

The InGaAs/InAlAs configs need `[strain]` section with reference to InP substrate (lattice-matched by construction, so strain is zero — but the strain section is needed to activate strain calculation for the optics pipeline).

**SC configs:** Reuse existing `sc_qw_gaas_alas.toml` from `test_sc_convergence.py`. May need a companion flat-band config (same geometry without `[sc]` section) for the subband shift comparison.

**ctest registration:**

- `verification_isbt_benchmark` under labels `verification;standard-star`
- `verification_sc_benchmark` under labels `verification;standard-star`
- Update `test_exciton_convergence` registration to include `verification` label

**COVERAGE annotations:** Each new verification script includes `# COVERAGE:` lines for its benchmark cells.

## Testing Decisions

### What makes a good test for publishable benchmarks

A publishable benchmark test must satisfy three criteria beyond a standard regression test:

1. **Independent reference:** The test compares against an analytical formula, published numerical data, or conservation law — never against the code's own output at a different resolution.
2. **Quantitative tolerance with physical justification:** The tolerance is not arbitrary — it reflects the expected accuracy of the reference (e.g., infinite-well vs finite-barrier, parabolic vs non-parabolic, Bastard variational vs 8-band).
3. **Reproducible figure or table:** The test produces a figure or data table that can be included in a publication without post-processing.

### Modules tested

| Module | Test type | Prior art |
|--------|-----------|-----------|
| `optical_spectra.f90` — `compute_intersubband_transitions` | ISBT z-dipole vs infinite-well analytical | `test_optical_qw.pf:test_isbt_dipole_velocity_consistency()` (internal consistency) |
| `optical_spectra.f90` — `compute_isbt_absorption` | ISBT spectrum peak position | `lecture_06_optical.py:test_isbt_peak()` (smoke test only) |
| `sc_loop.f90` — `self_consistent_loop` | Charge neutrality, Fermi level, subband shift | `test_sc_convergence.py` (Richardson convergence only) |
| `charge_density.f90` — `compute_charge_density_qw` | Integrated charge vs doping | `lecture_07_scsp.py` (soft warning only) |
| `poisson.f90` — box-integration Poisson | Potential profile shape | No existing shape test |
| `exciton.f90` — `compute_exciton_binding` | Eb vs Miller 1985, Bastard 1982 | `test_exciton_convergence.py` (30% tolerance) |

### Testing pattern

Follow the established verification script pattern (`verify_8band_rung*.py`, `verify_star_*.py`):

- Python script using `star_helpers.run_exe()` to invoke the Fortran executable
- Parse Fortran output files (`.dat` files in `output/`)
- Compare against analytical reference computed inline in the script
- Assert with quantitative tolerance
- Print PASS/FAIL per check
- Generate overlay plot (numerical vs analytical)
- Register as ctest target via shell wrapper in `tests/integration/`

### pFUnit considerations

No new pFUnit unit tests are planned for this phase — the benchmarks are integration-level tests that exercise the full Fortran executable against analytical references. The existing unit test coverage for optical_spectra (20 tests), sc_loop (17 tests), and exciton (indirectly via simulation_setup) is sufficient at the unit level.

## Out of Scope

- **VB g-factor validation:** Deferred to a second validation round. The Roth formula has no VB analog, and the comparison requires Winkler tabulated values or kdotpy cross-code — both higher effort.
- **Phase 18 architectural cleanup:** Separate workstream. The ISBT/SC/exciton benchmarks are additive and don't depend on architectural changes.
- **Phase 19 extended cross-code validation:** Blocked on kdotpy API capabilities. Not needed for these benchmarks (analytical comparisons suffice).
- **New material parameter entries:** All material systems use existing entries in `parameters.f90` (GaAs, Al30Ga70As, Ga47In53AsW, Al47In53AsW, AlAs).
- **Wire geometry benchmarks:** ISBT and SC benchmarks are QW-only. Wire ISBT and wire SC would require CSR-based eigensolver paths and are deferred.
- **Bulk ISBT:** ISBT is a confinement-induced phenomenon. No bulk ISBT benchmark is needed.
- **InAs/AlSb ISBT:** Dropped due to type-II band alignment making the infinite-well comparison misleading.
- **kdotpy cross-code for ISBT/SC/exciton:** The existing kdotpy pipeline (12 tests) doesn't cover these observables. Adding them would require kdotpy API extensions.
- **Sommerfeld enhancement validation:** The `sommerfeld_2d` function in `exciton.f90` is untested but is a supplementary correction, not the core variational method.
- **Excitonic corrections to absorption:** The `apply_excitonic_corrections` function modifies the absorption spectrum but is not the primary exciton observable.
- **Carrier-overlap corrections to ISBT beyond 1-band:**

## Further Notes

### Paper integration

These three benchmarks are designed to fill specific gaps in a J. Comput. Phys. methods paper's §4 results table. The expected paper structure is:

- §2 Method (FD discretization, Hamiltonian assembly, velocity matrices, SC loop, exciton variational)
- §3 Validation framework (verification ladder, standard stars, cross-code, Richardson, TRK sum rule, Hermiticity)
- §4 Results — the three new benchmarks fill gaps in this section:
  - §4.X ISBT: deviation from infinite-well formula as a function of well width and material system (new table + figure)
  - §4.Y SC: charge neutrality, Fermi level, subband shift, potential profile (new table)
  - §4.Z Exciton: tightened Miller comparison + Bastard reference (updated table)

### Effort estimate

| Item | Effort | Lines (est.) |
|------|--------|-------------|
| ISBT: 6 configs + verification script + figure | Medium | ~400 lines Python, 6 TOML configs |
| SC: verification script + quantitative checks | Medium | ~350 lines Python |
| Exciton: update convergence test + register in universe | Low | ~50 lines Python changes, universe.yml update |
| ctest registration + COVERAGE annotations | Low | ~30 lines CMakeLists.txt + shell wrappers |
| **Total** | ~2-3 sessions | ~830 lines Python, 6-7 TOML configs |

### Dependencies on existing infrastructure

- `tests/integration/star_helpers.py` — shared utility (run_exe, parse helpers, physical constants)
- `tests/integration/validation_universe.yml` — coverage matrix
- `tests/regression/configs/` — existing config templates
- `tests/regression/compare_output.py` — output comparison utilities
- `build/src/opticalProperties` — ISBT executable
- `build/src/bandStructure` — SC and exciton executable
