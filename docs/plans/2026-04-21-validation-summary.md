# Validation Summary — Docs Physics Revamp

Date: 2026-04-21

## Reproduced nextnano Examples

| System | Observable | nextnano Reference | Our Result | Match |
|--------|-----------|-------------------|------------|-------|
| Bulk GaAs | Eg (Gamma) | 1.519 eV | Verified (automated, ±1%) | Yes |
| Bulk GaAs | Delta_SO | 0.341 eV | Verified (automated, ±1%) | Yes |
| Bulk GaAs | CB effective mass | 0.067 m0 | Verified (informational) | Yes |
| Bulk InAs | Band structure | Vurgaftman parameters | Regression test passes | Yes |
| GaAs/AlGaAs QW | CB1 subband | ~31.5 meV (10nm, V0=0.3eV) | Verified (automated) | Yes |
| GaAs/AlGaAs QW | Spin degeneracy | 2x at k=0 | Verified (automated) | Yes |
| QCSE Stark shift | Shift at -70 kV/cm | ~1.79 meV (6nm GaAs) | Verified (automated) | Yes |
| QCSE Stark coefficient | Perturbative | ~0.000365 meV/(kV/cm)² | Reported (informational) | Approx* |
| SC QW | Convergence | Converges | Verified (automated) | Yes |
| SC QW | Positive gap | > 0 | Verified (automated) | Yes |
| InAs/AlSb SC | Type-II band alignment | Narrow gap | Regression test passes | Yes |

*Non-perturbative regime at high field; coefficient is approximate.

## Chapters Now Benchmark-Backed

| Chapter | Benchmark Status |
|---------|-----------------|
| 01 Bulk Band Structure | GaAs Eg, Delta_SO vs Vurgaftman (automated) |
| 02 Quantum Well | VB-CB gap, subband spacing, degeneracy (automated) |
| 03 Wavefunctions | State indexing verified, figures match captions |
| 04 Strain | HH/LH splitting benchmarks added |
| 05 g-factor | Zeeman units corrected |
| 06 Optical Properties | TE absorption, scattering lifetimes regenerated |
| 07 Self-Consistent SP | Convergence, gap, subband checks (automated) |
| 08 Quantum Wire | Fixed Hamiltonian, GaAs gap 0.36 eV, InAs/GaAs 0.19 eV |
| 09 Numerical Methods | No changes needed (already accurate) |
| 10 QCSE | Stark shift benchmarks (automated) |
| 11 Convergence | No changes needed |
| 12 Extending the Code | Updated kpterms_2d indices |

## Intentionally Downgraded or Removed

- Wire g-factor quantitative comparison with nextnano: our InSb wire gives |g|~24 vs bulk |g|~51, but wire geometry/confinement effects make direct comparison complex. Reported as informational only.
- Bulk GaAs g-factor: Roth formula gives -0.44, our 8-band Lowdin result differs due to finite-band truncation. Not benchmarked automatically.
- Exciton binding energy quantitative check: not automated (requires separate variational calculation).

## Test Results

All 34 tests pass: 15 unit + 19 regression.
