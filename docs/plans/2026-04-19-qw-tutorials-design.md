# QW Tutorial Coverage: Full Design (Tiers 1-3)

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add optical spectrum, exciton, and scattering capabilities to the 8-band k.p solver, enabling reproduction of all 13 nextnano QW tutorials.

**Architecture:** Three tiers of new physics modules, each building on the previous. Tier 1 (optical spectra) reuses existing momentum matrix element code. Tier 2 (excitons) adds a variational solver. Tier 3 (scattering) adds Fröhlich phonon rates. Each tier produces new Fortran subroutines, Python plotting functions, and lecture chapter sections.

**Tech Stack:** Fortran 90 (new physics), Python/matplotlib (figures), existing infrastructure (pMatrixEleCalc, fermi_dirac, Simpson integration, k_par grid).

---

## Tutorial Coverage Map

### Group A: Band Structure & Dispersion (COVERED)

| Tutorial | Topic | Our Coverage |
|----------|-------|-------------|
| 5.9.1 Broken-gap QW | InAs/GaSb/AlSb type-II, E(k_par), anticrossing | Done (Ch02, Ch03) |
| 5.9.5 k.p dispersion | VB subbands, strain, 6-band vs 8-band | Done (Ch02, Ch04) |

### Group B: Interband Optics (Tier 1a)

| Tutorial | Topic | Gap |
|----------|-------|-----|
| 5.9.7 Optics tutorial | Interband+intersubband absorption, Fermi's golden rule, polarization, occupation | k_par-integrated absorption spectrum alpha(E) |
| 5.9.8 Overlap integrals | Envelope function overlaps, delta-n selection rule | 1-band overlap code (small) |
| 5.9.9 InGaAs absorption | Strained QW absorption, TE/TM selection rules | Same as 5.9.7 |

### Group C: Intersubband Transitions (Tier 1b)

| Tutorial | Topic | Gap |
|----------|-------|-----|
| 5.9.10 ISBT infinite QW | ISBT absorption, dipole moments, doping effects | z-dipole extraction from existing code |
| 5.9.11 ISBT MQW | Single-band vs 8-band, coupled QWs, nonparabolicity | Same as 5.9.10 + multi-QW configs |

### Group D: PL & Gain (Tier 1c)

| Tutorial | Topic | Gap |
|----------|-------|-----|
| 5.9.12 Resonant PL | Spontaneous emission, generation rate, Beer's law | Spontaneous emission spectrum |
| 5.9.13 Non-resonant PL | Barrier absorption, carrier trapping, drift-diffusion | DEFERRED (needs drift-diffusion solver) |
| 5.9.14 Strained gain | TE/TM gain, quasi-Fermi levels, strain comparison | Gain spectrum from absorption + inversion |

### Group E: Excitons (Tier 2)

| Tutorial | Topic | Gap |
|----------|-------|-----|
| 5.9.2 Exciton binding energy | Variational solver, binding energy vs well width | Variational exciton solver |
| 5.9.15 Exciton absorption | Exciton peaks, Sommerfeld enhancement, 8-band k.p | Same + Sommerfeld factor |

### Group F: Scattering (Tier 3)

| Tutorial | Topic | Gap |
|----------|-------|-----|
| 5.9.3 Electron scattering | LO-phonon Fröhlich scattering, lifetimes, double QW | Fröhlich scattering rate code |

---

## Existing Infrastructure

### What Already Exists (Reusable)

| Component | Location | Notes |
|-----------|----------|-------|
| Momentum matrix elements (QW) | `gfactor_functions.f90`: `pMatrixEleCalc` | dH/dk in x,y,z directions |
| Optical transition type | `defs.f90`: `type optical_transition` | cb_idx, vb_idx, energy, px, py, pz, f_osc |
| Optical matrix at k=0 (QW) | `gfactor_functions.f90`: `compute_optical_matrix_qw` | CB-VB pairs at single k-point |
| Fermi-Dirac distribution | `charge_density.f90`: `fermi_dirac` | With overflow protection |
| Fermi level finder | `sc_loop.f90`: `find_fermi_level` | Bisection, charge neutrality |
| k_par grid construction | `sc_loop.f90` (lines 124-128) | Linear 0..kpar_max, odd count |
| Simpson integration | `utils.f90`: `simpson`, `simpson_real` | Complex + real, odd-point requirement |
| k_par integration pattern | `charge_density.f90`: `accumulate_band_density` | Pattern for accumulation loop |
| Spin matrices | `gfactor_functions.f90`: `SIGMA_X/Y/Z` | 8x8 in zincblende basis |
| Physical constants | `defs.f90` | hbar, m0, e, kB, pi, hbar2O2m0 |
| Strain (Bir-Pikus) | `hamiltonianConstructor.f90` | Full strain Hamiltonian |
| External electric field | `main.f90` | EFParams, QCSE support |
| Self-consistent SP | `sc_loop.f90` | Schrödinger-Poisson with DIIS |

### What Needs Writing

1. k_par-dependent optical matrix element loop (embed pMatrixEleCalc in k_sweep)
2. Absorption coefficient: alpha(E) = (e^2 pi)/(n_r c eps0 m0^2 E) * sum_CV int dk_par |M_CV|^2 * (f_V - f_C) * L(E - E_CV)
3. Lineshape broadening (Lorentzian + Gaussian -> Voigt)
4. Gain = negative absorption with separate quasi-Fermi levels
5. Intersubband z-dipole extraction (<psi_CB1 | z | psi_CB2>)
6. Variational exciton solver (Bastard 1982)
7. Sommerfeld enhancement factor (2D)
8. Fröhlich LO-phonon scattering rate
9. Output routines for all new spectra
10. New input.cfg parameters (linewidth, refractive index, carrier density, etc.)
11. Python plotting functions for all new figures
12. Lecture chapter sections

---

## Tier 1: Optical Spectra (Absorption, ISBT, Gain)

### Tier 1a: k_par-Integrated Absorption Spectrum

**New Fortran code:**

1. **`src/physics/optical_spectra.f90`** (new module)
   - `compute_absorption_qw(config, eigvals_k, eigvecs_k, k_grid, E_grid, alpha_TE, alpha_TM)`:
     For each k_par point:
     - Solve eigenproblem (already done in main k_sweep)
     - Compute momentum matrix elements between CB-VB pairs
     - Accumulate: alpha(E) += |M_CV|^2 * (f_V - f_C) * L(E - dE_CV) * k_par_weight
     - TE: M_TE^2 = |M_x|^2 + |M_y|^2, TM: M_TM^2 = |M_z|^2
   - `lorentzian(E, E0, gamma)` and `gaussian(E, E0, sigma)` lineshape functions
   - `voigt(E, E0, gamma_l, gamma_g)` approximation (numerical or Humlicek)

2. **`src/core/defs.f90`** additions:
   - New `optics_config` type: `linewidth_lorentzian`, `linewidth_gaussian`, `refractive_index`, `temperature_optics`, `num_energy_points`, `E_min`, `E_max`
   - Extend `simulation_config` with `optics_config`

3. **`src/io/input_parser.f90`** additions:
   - Parse new keywords: `Optics: T`, `OpticsLinewidthL: 0.030`, `OpticsLinewidthG: 0.005`, `OpticsNr: 3.3`, `OpticsEmin: 0.8`, `OpticsEmax: 2.0`, `OpticsNpts: 200`

4. **`src/apps/main.f90`** additions:
   - After k_sweep loop completes, if optics enabled:
     - Call `compute_absorption_qw` with all k_par eigenvectors
     - Write `output/absorption_TE.dat` and `output/absorption_TM.dat`

5. **Python plotting** (`generate_all_figures.py`):
   - `fig_qw_absorption_spectrum()`: TE + TM absorption vs photon energy
   - `fig_qw_absorption_strained()`: Comparison unstrained vs strained InGaAs
   - `fig_qw_absorption_vs_width()`: Well-width parameter sweep

6. **Lecture docs** (`docs/lecture/06-optical-properties.md` additions):
   - Section on absorption coefficient formula derivation
   - k_par integration methodology
   - TE/TM polarization selection rules with numerical examples
   - Comparison with nextnano InGaAs tutorial (Dumitras PRB 2002)

**Effort:** Medium (~500 lines Fortran, ~200 lines Python, ~1000 words docs)

### Tier 1b: Intersubband Transitions

**New Fortran code:**

1. Extend `optical_spectra.f90`:
   - `compute_intersubband_qw(config, eigvals_k, eigvecs_k, z_grid, z_dipole)`:
     Compute z-dipole matrix elements between CB subbands: z_ij = <psi_i | z | psi_j>
     This requires the position operator z, not dH/dk. For the envelope function:
     z_ij = sum_n psi_i(z_n) * z_n * psi_j(z_n) * dz
   - `compute_isbt_absorption_qw()`: ISBT absorption using z-dipole + Fermi occupation
   - Oscillator strength: f_ij = (2 m0 E_ij / hbar^2) |z_ij|^2

2. **Output:** `output/isbt_transitions.dat` (i, j, E_ij, z_ij, f_ij)

3. **Python plotting:**
   - `fig_isbt_dipole_moments()`: z_ij bar chart
   - `fig_isbt_absorption()`: ISBT absorption spectrum

4. **Lecture docs:** Section on intersubband selection rules, z-dipole, TM polarization

**Effort:** Small (~200 lines Fortran, ~100 lines Python)

### Tier 1c: Gain Spectrum

**New Fortran code:**

1. Extend `optical_spectra.f90`:
   - `compute_gain_qw()`: Same integrand as absorption but with separate quasi-Fermi levels
   - Requires: carrier density input (e.g., `GainCarrierDensity: 3e12` in cm^-2)
   - Find separate f_e and f_h for given 2D carrier density
   - Gain(E) = -alpha(E) when f_C > f_V (population inversion)
   - TE and TM gain separately

2. **Input:** `Gain: T`, `GainCarrierDensity: 3.0e12`, `GainStrain: compressive/tensile/unstrained`

3. **Output:** `output/gain_TE.dat`, `output/gain_TM.dat`

4. **Python plotting:**
   - `fig_gain_strained_comparison()`: Side-by-side gain for 3 strain conditions (reproduces nextnano 5.9.14)

5. **Lecture docs:** Section on gain, quasi-Fermi levels, strain-dependent TE/TM gain

**Effort:** Small (~150 lines Fortran, ~100 lines Python) -- builds on Tier 1a

---

## Tier 2: Excitons

### Tier 2a: Variational Exciton Solver

**New Fortran code:**

1. **`src/physics/exciton.f90`** (new module):
   - `compute_exciton_binding(config, eigvals, eigvecs, z_grid, E_binding, lambda_exc)`:
     Variational calculation following Bastard (PRB 1982):
     - Trial function: psi(r, ze, zh) = phi_e(ze) * phi_h(zh) * (2/pi) * (1/lambda) * exp(-r/lambda)
     - phi_e, phi_h = ground-state CB1 and HH1 envelope functions
     - Minimize <psi | H_ex | psi> w.r.t. lambda
     - H_ex = -hbar^2/(2*mu) * nabla^2_r - e^2/(4*pi*eps*eps0) * 1/sqrt(r^2 + (ze-zh)^2)
     - Return binding energy and exciton Bohr radius
   - Use 1D Simpson for z-integration, golden-section search for lambda minimization
   - Material parameters: dielectric constant from paramStruct, reduced mass from band character

2. **Output:** `output/exciton.dat` (E_binding in meV, lambda in nm)

3. **Python plotting:**
   - `fig_exciton_binding_vs_width()`: Binding energy vs well width (reproduces Harrison Fig 6.4)
   - `fig_exciton_bohr_vs_width()`: Bohr radius vs well width (reproduces Harrison Fig 6.5)

4. **Lecture docs:** New section in Ch06 on excitons in QWs, Bastard variational method

**Effort:** Medium (~400 lines Fortran, ~150 lines Python)

### Tier 2b: Sommerfeld Enhancement Factor

**New Fortran code:**

1. Extend `exciton.f90`:
   - `sommerfeld_2d(E_excess, E_binding)`: Returns S_2D = exp(pi/sqrt(D)) / cosh(pi/sqrt(D))
     where D = E_excess / (E_binding/4)
   - Apply to absorption continuum above band edge
   - Add exciton peak: delta-function peak at E = E_gap - E_binding, broadened by Voigt

2. Extend `compute_absorption_qw` to optionally include excitonic corrections

3. **Input:** `Exciton: T`, `ExcitonMethod: variational`

4. **Python plotting:**
   - `fig_absorption_with_exciton()`: Three curves: no exciton, Coulomb enhancement only, full excitonic

5. **Lecture docs:** Sommerfeld factor derivation, comparison of 2D vs 3D exciton limits

**Effort:** Small (~150 lines Fortran, ~100 lines Python) -- builds on Tier 2a

---

## Tier 3: Scattering

### Tier 3a: LO-Phonon Fröhlich Scattering

**New Fortran code:**

1. **`src/physics/scattering.f90`** (new module):
   - `compute_phonon_scattering(config, eigvals, eigvecs, z_grid, lifetimes)`:
     Fröhlich Hamiltonian for LO-phonon scattering (Ferreira & Bastard PRB 1989):
     - Scattering rate: 1/tau = (2*pi/hbar) * sum_q |M(q)|^2 * (N_q + 1/2 +/- 1/2) * delta(E_f - E_i +/- hbar*omega_LO)
     - M(q) = e^2 * hbar*omega_LO * (1/eps_inf - 1/eps_0) * (1/(V*q^2)) * <psi_f | e^{iq_z z} | psi_i>
     - For intersubband: integrate over in-plane q, sum over q_z modes
     - Parameters: LO phonon energy (material-dependent), eps_inf, eps_0
   - Output lifetimes in picoseconds
   - Works for single QW and double QW (coupled states)

2. **Input:** `Scattering: T`, `PhononEnergy: 0.036`, `EpsInf: 10.9`, `Eps0: 12.9`

3. **Output:** `output/scattering_rates.dat` (i, j, E_ij, rate, lifetime_ps)

4. **Python plotting:**
   - `fig_scattering_lifetime_vs_width()`: Lifetime vs QW width (reproduces Ferreira Fig 3)
   - `fig_scattering_lifetime_vs_field()`: Lifetime vs electric field (reproduces Ferreira Fig 5)
   - `fig_double_qw_anticrossing()`: Energy levels of asymmetric double QW

5. **Lecture docs:** New section on carrier scattering, Fröhlich interaction, LO-phonon emission/absorption

**Effort:** Medium (~400 lines Fortran, ~200 lines Python)

---

## Implementation Order

```
Phase 2 (Tier 1): Optical Spectra
  Task 1:  defs.f90 + input_parser.f90 -- optics_config type and input parsing
  Task 2:  optical_spectra.f90 -- absorption coefficient core
  Task 3:  main.f90 -- integrate optics into k_sweep loop
  Task 4:  Python figures -- absorption spectrum plots
  Task 5:  Lecture docs -- Ch06 absorption sections
  Task 6:  optical_spectra.f90 -- intersubband z-dipole (Tier 1b)
  Task 7:  Python figures -- ISBT plots
  Task 8:  optical_spectra.f90 -- gain spectrum (Tier 1c)
  Task 9:  Python figures -- gain plots
  Task 10: Lecture docs -- Ch06 gain + ISBT sections

Phase 3 (Tier 2): Excitons
  Task 11: exciton.f90 -- variational solver
  Task 12: Python figures -- binding energy plots
  Task 13: exciton.f90 -- Sommerfeld enhancement
  Task 14: Python figures -- excitonic absorption
  Task 15: Lecture docs -- Ch06 exciton sections

Phase 4 (Tier 3): Scattering
  Task 16: scattering.f90 -- Fröhlich LO-phonon rates
  Task 17: Python figures -- lifetime plots
  Task 18: Lecture docs -- scattering section
```

## New Input Parameters Summary

```
! --- Optical Spectra (Tier 1) ---
Optics: T                        ! Enable absorption/gain calculation
OpticsLinewidthL: 0.030          ! Lorentzian FWHM (eV)
OpticsLinewidthG: 0.005          ! Gaussian FWHM (eV)
OpticsNr: 3.3                    ! Refractive index
OpticsEmin: 0.5                  ! Spectrum energy minimum (eV)
OpticsEmax: 2.0                  ! Spectrum energy maximum (eV)
OpticsNpts: 200                  ! Number of energy grid points
OpticsTemperature: 300           ! Temperature for occupation (K)
OpticsCarrierDensity: 0.0        ! 2D carrier density (cm^-2), 0 = equilibrium

! --- Gain (Tier 1c) ---
Gain: F                          ! Enable gain calculation
GainCarrierDensity: 3.0e12       ! 2D carrier density for gain (cm^-2)

! --- Intersubband (Tier 1b) ---
ISBT: F                          ! Enable intersubband transitions

! --- Exciton (Tier 2) ---
Exciton: F                       ! Enable exciton calculation
ExcitonMethod: variational       ! Solver method

! --- Scattering (Tier 3) ---
Scattering: F                    ! Enable scattering rate calculation
PhononEnergy: 0.036              ! LO phonon energy (eV), material-dependent
EpsInf: 10.9                     ! High-frequency dielectric constant
Eps0: 12.9                       ! Static dielectric constant
```

## New Output Files Summary

```
output/absorption_TE.dat         ! Absorption coefficient vs E (TE polarization)
output/absorption_TM.dat         ! Absorption coefficient vs E (TM polarization)
output/isbt_transitions.dat      ! Intersubband transitions: i, j, E, z_ij, f_ij
output/gain_TE.dat               ! Optical gain vs E (TE polarization)
output/gain_TM.dat               ! Optical gain vs E (TM polarization)
output/exciton.dat               ! Exciton binding energy and Bohr radius
output/absorption_excitonic_TE.dat  ! Absorption with excitonic corrections
output/scattering_rates.dat      ! Phonon scattering rates and lifetimes
```

## New Python Figures Summary

| Figure | Tier | Reproduces Nextnano |
|--------|------|-------------------|
| `qw_absorption_spectrum.png` | 1a | 5.9.7 (QWIP), 5.9.9 (InGaAs) |
| `qw_absorption_strained.png` | 1a | 5.9.9 (InGaAs strained) |
| `isbt_dipole_moments.png` | 1b | 5.9.10 (infinite QW) |
| `isbt_absorption.png` | 1b | 5.9.10, 5.9.11 (MQW) |
| `gain_strained_comparison.png` | 1c | 5.9.14 (3 strain cases) |
| `exciton_binding_vs_width.png` | 2a | 5.9.2 (Harrison Fig 6.4) |
| `exciton_bohr_vs_width.png` | 2a | 5.9.2 (Harrison Fig 6.5) |
| `absorption_with_exciton.png` | 2b | 5.9.15 (3 curves) |
| `scattering_lifetime_vs_width.png` | 3a | 5.9.3 (Ferreira Fig 3) |
| `scattering_lifetime_vs_field.png` | 3a | 5.9.3 (Ferreira Fig 5) |
| `double_qw_anticrossing.png` | 3a | 5.9.3 (Ferreira Fig 7) |

## Key Design Decisions

1. **Optics module is separate from SC loop.** The absorption calculation reuses the k_par sweep pattern but does not modify the self-consistent loop. It runs as a post-processing step after the band structure sweep.

2. **Eigenvector storage.** Currently `main.f90` only stores eigenvectors at k_par=0 for output. For absorption, we need eigenvectors at ALL k_par points. This is the biggest memory concern. Solution: accumulate absorption integrand on-the-fly during the k_sweep loop (no need to store all eigenvectors simultaneously).

3. **Momentum matrix elements at each k_par.** The existing `pMatrixEleCalc` builds dH/dk from the Hamiltonian. For each k_par point, we already have eigenvectors from the k_sweep. We just need to call `pMatrixEleCalc` for each CB-VB pair at each k_par point. This is O(n_CB * n_VB * n_kpar) calls, each requiring a Hamiltonian build + dot product. Performance concern for large systems.

4. **Intersubband z-dipole.** Unlike interband transitions (which use dH/dk), intersubband transitions use the position operator z. The envelope function approximation gives z_ij = sum_n psi_i(z_n) * z_n * psi_j(z_n) * dz. This is a simple quadrature -- no Hamiltonian build needed.

5. **Exciton solver uses 1-band effective mass.** The variational solver takes the CB1 and HH1 envelope functions from the 8-band k.p calculation but treats the relative electron-hole motion with a 1s trial function. This is the Bastard (1982) approach and matches what nextnano does.

6. **Scattering uses Fröhlich model.** The LO-phonon scattering rate is computed from the Fröhlich coupling constant, integrating over in-plane phonon wavevectors. This follows Ferreira & Bastard (1989) exactly as in the nextnano tutorial.

7. **Non-resonant PL (Tutorial 5.9.13) is DEFERRED.** It requires a full drift-diffusion solver for carrier transport, which is a very large effort beyond the scope of this plan.

## Validation Strategy

Each tier validated against published references:
- **Tier 1a:** Dumitras et al., PRB 2002 (InGaAs absorption) -- nextnano tutorial 5.9.9
- **Tier 1b:** Chuang, Physics of Optoelectronic Devices (1995) -- nextnano tutorial 5.9.10
- **Tier 1c:** Chuang, Physics of Optoelectronic Devices (1995), Fig. 10.30 -- nextnano tutorial 5.9.14
- **Tier 2a:** Harrison, Quantum Wells, Wires and Dots (2005), Figs. 6.4-6.5 -- nextnano tutorial 5.9.2
- **Tier 2b:** 2D Sommerfeld factor analytical limits -- nextnano tutorial 5.9.15
- **Tier 3a:** Ferreira & Bastard, PRB 1989 -- nextnano tutorial 5.9.3
