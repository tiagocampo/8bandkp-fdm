# QW Documentation Overhaul + Physics Roadmap

## Context

The nextnano++ documentation has 13 quantum well tutorials covering band structure, optical absorption, intersubband transitions, photoluminescence, gain, excitons, and scattering. We want to cover equivalent physics in our documentation, identify what new code is needed, and plan a phased implementation.

**Goal:** A multi-phase roadmap that (1) overhauls QW documentation with what we can already compute, and (2) adds new physics modules to cover the remaining nextnano tutorials.

---

## Tutorial Grouping (13 tutorials → 6 clusters)

| Cluster | Tutorials | Nextnano Ref |
|---------|-----------|-------------|
| **A. Band Structure & Dispersion** | k.p dispersion QW, broken-gap QW | 5.9.4, 5.9.5 |
| **B. Interband Optics** | optics tutorial, overlap integrals, InGaAs absorption | 5.9.7, 5.9.8, 5.9.9 |
| **C. Intersubband Optics** | ISBT infinite QW, ISBT MQW | 5.9.10, 5.9.11 |
| **D. Luminescence & Gain** | resonant PL, non-resonant PL, strained gain | 5.9.12, 5.9.13, 5.9.14 |
| **E. Excitons** | binding energy, exciton absorption | 5.9.2, 5.9.15 |
| **F. Scattering** | electron scattering QW | 5.9.3 |

## Capability Gap Analysis

### What we already have (no new code needed)
- 8-band k.p Hamiltonian (zincblende, bulk + QW + wire)
- FD eigenvalue solver (dense LAPACK + sparse MKL)
- Bir-Pikus strain
- E(k_parallel) dispersion sweeps with eigenvector storage at every k
- Self-consistent Schrodinger-Poisson (charge neutrality + fixed Fermi)
- Broken-gap support (InAsW/GaSbW configs)
- **Optical momentum matrix elements** (`pMatrixEleCalc` in `gfactor_functions.f90`) — computes `<psi_i | dH/dk_dir | psi_j>` for QW (dense) and wire (CSR), axis-aligned x/y/z
- **k_parallel integration** infrastructure in `charge_density.f90` (Simpson rule, cylindrical weight k/(2pi), Fermi occupation)
- Material database (25+ semiconductors)
- External electric field (QCSE)

### What we need to add

| Tier | Capability | Enables | Effort |
|------|-----------|---------|--------|
| **1a** | k_parallel-integrated absorption spectrum alpha(E) | Clusters B, D (partial) | Medium |
| **1b** | Intersubband z_ij dipole moments | Cluster C | Small |
| **1c** | Gain spectrum with quasi-Fermi levels | Cluster D (gain) | Small (builds on 1a) |
| **2a** | Variational exciton solver | Cluster E | Medium |
| **2b** | Sommerfeld enhancement factor | Cluster E (absorption) | Small |
| **3a** | LO-phonon scattering rate (Frohlich) | Cluster F | Medium |
| **3b** | Drift-diffusion solver | Non-resonant PL only | Very Large (deferred) |

---

## Phase 1: Documentation Overhaul (no new Fortran)

**Scope:** Rewrite QW docs to cover Clusters A with existing code. Note gaps for B-F.

### Files to modify
- `docs/lecture/02-quantum-well.md` — add broken-gap E(k_parallel) example
- `docs/lecture/06-optical-properties.md` — add QW optical matrix element example (using existing `pMatrixEleCalc`)
- New tutorial chapters for dispersion and broken-gap walkthroughs

### New configs needed (in `tests/regression/configs/`)
- `qw_gaas_algaas_kpar.cfg` — GaAs/AlGaAs QW fine k_parallel sweep
- `qw_inas_gasb_broken_gap_kpar.cfg` — InAs/GaSb broken-gap E(k_parallel) showing anticrossing

### New figures needed (in `scripts/plotting/generate_all_figures.py`)
- `fig_qw_dispersion_gaas_algaas` — E(k_parallel) subbands
- `fig_qw_dispersion_broken_gap` — e-h anticrossing
- `fig_qw_optical_matrix_elements` — oscillator strength bar chart

### Enables fully: Cluster A (2 tutorials). Partially: Cluster B (1 tutorial).

---

## Phase 2: Absorption/Gain Spectrum Module (Tier 1)

**Scope:** Implement `alpha(E)` for interband + intersubband. Enables Clusters B, C, D.

### New module: `src/physics/optical_spectra.f90`

**Core algorithm** (mirrors `compute_charge_density_qw` pattern from `charge_density.f90`):
1. Loop over k_parallel points (Simpson integration, cylindrical weight)
2. At each k_par: compute matrix elements via existing `pMatrixEleCalc` (dense QW mode in `gfactor_functions.f90`)
3. Weight by Fermi occupation `(f_VB - f_CB)`
4. Broaden delta with configurable Lorentzian/Gaussian FWHM
5. Output TE = `alpha_x + alpha_y`, TM = `alpha_z`

**Intersubband z_ij:** `<psi_i | z | psi_j>` via Simpson integration of CB envelope overlap.

**Gain:** Same formula with quasi-Fermi levels for non-equilibrium injection.

### New input.cfg fields
- `Optics 1`, `Optics_energy_min/max/step`, `Optics_broadening`, `Optics_num_kpar`, `Optics_kpar_max`, `Optics_type`, `Optics_fermi_mode`, `Optics_injection_n2d`

### New output files
- `output/absorption_te.dat`, `output/absorption_tm.dat`, `output/intersubband.dat`

### Configs: `optics_gaas_algaas_interband.cfg`, `optics_ingaas_strained.cfg`, `optics_gaas_isbt.cfg`, `gain_ingaas_strained.cfg`

### Enables fully: Clusters B (3 tutorials), C (2 tutorials), D-partial (gain). Total after Ph1+2: 7 full, 4 partial.

---

## Phase 3: Exciton Solver (Tier 2)

**Scope:** Variational exciton + Sommerfeld enhancement. Enables Cluster E.

### New module: `src/physics/exciton.f90`

**Algorithm:** Bastard separable ansatz, minimize E_b(lambda) over Bohr radius lambda. Uses density-averaged reduced mass from envelope wavefunctions.

**Sommerfeld:** `S_2D = exp(pi/Delta) / cosh(pi/Delta)` as continuum enhancement.

### Configs: `exciton_cdte_infinite_qw.cfg`, `exciton_gaas_infinite_qw.cfg`

### Enables fully: Cluster E (2 tutorials). Total after Ph1+2+3: 10 full, 1 partial.

---

## Phase 4: Scattering (Tier 3)

**Phase 4a:** Frohlich LO-phonon scattering (`src/physics/scattering.f90`). Uses z_ij from Phase 2. Benchmark vs Ferreira & Bastard PRB 1989.

**Phase 4b (deferred):** Drift-diffusion for non-resonant PL. Very large effort, recommend documenting limitation.

### Enables fully: Cluster F (1 tutorial). Total: 11 full, 1 deferred.

---

## Tutorial Coverage Summary

| Tutorial | Phase 1 | Phase 2 | Phase 3 | Phase 4 |
|----------|---------|---------|---------|---------|
| 5.9.4 k.p Dispersion | **Full** | — | — | — |
| 5.9.5 Broken Gap | **Full** | — | — | — |
| 5.9.7 Optics (interband) | Partial | **Full** | **Full** | — |
| 5.9.8 Overlap Integrals | Partial | **Full** | — | — |
| 5.9.9 InGaAs Absorption | Partial | **Full** | **Full** | — |
| 5.9.10 ISBT Infinite QW | — | **Full** | — | — |
| 5.9.11 ISBT MQW | — | **Full** | — | — |
| 5.9.12 Resonant PL | — | Partial | **Full** | — |
| 5.9.13 Non-resonant PL | — | Partial | Partial | Deferred |
| 5.9.14 Strained Gain | — | **Full** | **Full** | — |
| 5.9.2 Exciton Binding | — | — | **Full** | — |
| 5.9.15 Exciton Absorption | — | Partial | **Full** | — |
| 5.9.3 Scattering | — | — | — | **Full** |

## Verification Benchmarks

- **Phase 2:** GaAs infinite QW E21=166.5 meV, z21=-1.82 nm (Chuang); InGaAs strained TE peak at ~1.303 eV (Dumitras PRB 2002)
- **Phase 3:** CdTe bulk E_ex=-10.0 meV, lambda=6.8 nm (Harrison); 2D limit = 4x bulk
- **Phase 4a:** GaAs scattering forbidden below 5.4 nm and above 18 nm well width (Ferreira-Bastard)
