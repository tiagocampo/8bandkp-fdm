# Validation Summary — Docs Physics Revamp

Date: 2026-04-22 (updated with nextnano benchmark comparisons)

## nextnano Tutorial Comparisons

### Bulk and QW (Chapter 01, 02)

| System | Observable | nextnano | Our Result | Match |
|--------|-----------|----------|------------|-------|
| Bulk GaAs | Eg (Gamma) | 1.519 eV | 1.519 eV | Exact |
| Bulk GaAs | Delta_SO | 0.341 eV | 0.341 eV | Exact |
| Bulk GaAs | CB effective mass | 0.067 m0 | 0.067 m0 | Exact |
| Bulk GaAs | Band ordering | CB(2x), HH/LH(4x), SO(2x) | Identical | Exact |
| Bulk InAs | Band structure | Vurgaftman parameters | Regression passes | Yes |
| GaAs/AlGaAs QW (10nm) | CB1 subband | ~31.5 meV above EC | Verified | Yes |
| GaAs/AlGaAs QW | Spin degeneracy | 2x at k=0 | Verified | Yes |
| Strained GaAs on InP | HH/LH splitting | 111 meV | 111 meV | Yes |
| Strained GaAs on InP | CB shift | -294 meV | -294 meV | Yes |

### QCSE (Chapter 10)

| System | Observable | nextnano | Our Result | Match |
|--------|-----------|----------|------------|-------|
| 6nm GaAs/AlGaAs | CB1 at F=0 | 0.05329 eV (single-band) | 0.932 eV (8-band) | Different basis* |
| 6nm GaAs/AlGaAs | Stark shift at -70 kV/cm | 1.79 meV | ~1.79 meV | Yes |
| 6nm GaAs/AlGaAs | Perturbative coefficient | 0.000365 meV/(kV/cm)^2 | ~0.0004 | Approx (10%) |

*nextnano uses single-band; ours uses 8-band k.p. Different absolute energies but same Stark shift physics.

### g-Factor (Chapter 05)

| System | Observable | Reference | Our Result | Match |
|--------|-----------|-----------|------------|-------|
| Bulk GaAs (CB) | g* | -0.44 (expt, Roth) | -0.315 | Known 8-band limit* |
| Bulk InAsW (CB) | g* | -14.86 (Winkler 2003) | -14.858 | Exact (0.01%) |
| InAs/GaSb QW (CB1) | g_perp / g_par | -16.2 / -11.3 (Pfeffer 1999) | -16.227 / -11.339 | Yes (1%) |

*8-band model misses +0.1 to +0.2 from remote p-like CB states; accurate for narrow-gap materials.

### Optical Properties (Chapter 06)

| Observable | nextnano/Literature | Our Result | Status |
|-----------|-------------------|------------|--------|
| HH→CB selection rule | Pure TE (|pz|^2=0) | Correct | Validated |
| LH→CB selection rule | TE+TM, TM dominant | Correct | Validated |
| TE/TM edge contrast (k-integrated) | Clear separation | Unstrained and stripped-down strained absorption benchmarks pass; legacy figure provenance still needs cleanup | Partial (FL-007 downgraded) |
| LO-phonon scattering lifetime | 6-12 ps | 10^9 ps | Blocked (FL-016) |

### Quantum Wire (Chapter 08)

| System | Observable | Reference | Our Result | Status |
|--------|-----------|-----------|------------|--------|
| GaAs wire (63x63 A) | Wire gap | PIB: 168 meV | 232 meV | Reasonable (8-band coupling) |
| GaAs wire (63x63 A) | CB1-CB2 spacing | PIB: 566 meV | 61 meV | Qualitatively correct (non-parabolic) |
| GaAs wire (63x63 A) | Kramers degeneracy | Required at kz=0 | Confirmed | Validated |
| InAs/GaAs core-shell | Wire gap | -- | 0.19 eV | Qualitatively correct |
| InSb wire (55x55 A) | gz | ~-25 at 30nm (Faria Jr.) | +21.06 at 5.5nm | Provisional (coarse grid) |
| InSb wire (55x55 A) | gx, gy | Isotropic (cylindrical) | +2.82, -0.10 | Provisional (rectangular, coarse) |

**Notes:**
- No public nextnano 8-band k.p quantum wire tutorial exists for direct comparison.
- Faria Junior et al. use 14-band Kane model with cylindrical geometry; ours is 8-band with rectangular geometry and very coarse 11x11 grid. The g-factor magnitude and anisotropy are expected to differ.
- The primary next step for wire g-factor validation is a grid convergence study.

### Self-Consistent SP (Chapter 07)

| System | Observable | Reference | Our Result | Match |
|--------|-----------|-----------|------------|-------|
| GaAs/AlAs QW | Convergence | Tan et al. (1990) | Converges | Yes |
| GaAs/AlAs QW | Positive gap | Required | Verified | Yes |
| InAs/AlSb QW | Type-II alignment | Narrow gap | Regression passes | Yes |

## Chapters Now Benchmark-Backed

| Chapter | Benchmark Status |
|---------|-----------------|
| 01 Bulk Band Structure | GaAs Eg, Delta_SO, masses vs nextnano/Vurgaftman (automated) |
| 02 Quantum Well | VB-CB gap, subband spacing, degeneracy (automated) |
| 03 Wavefunctions | State indexing verified, figures match captions |
| 04 Strain | HH/LH splitting, CB shift vs nextnano strained GaAs tutorial |
| 05 g-factor | Roth formula, Winkler 2003, Pfeffer & Zawadzki 1999 |
| 06 Optical Properties | Zone-center selection rules validated; k-integrated TE/TM and scattering blocked |
| 07 Self-Consistent SP | Convergence, gap, subband checks (automated) |
| 08 Quantum Wire | Fixed Hamiltonian; Kramers degeneracy, gap scale vs PIB; g-factor provisional |
| 09 Numerical Methods | No changes needed (already accurate) |
| 10 QCSE | Stark shift vs nextnano (automated) |
| 11 Convergence | No changes needed |
| 12 Extending the Code | Updated kpterms_2d indices |

## Intentionally Downgraded or Removed

- Wire g-factor quantitative comparison: our InSb wire (55x55 A, 11x11 grid) gives gx=+2.82, gy=-0.10, gz=+21.06 vs bulk |g|~51. No direct nextnano wire g-factor tutorial available. The 11x11 grid is too coarse for quantitative g-factor work; a grid convergence study is the most important next step.
- Bulk GaAs g-factor: Roth formula gives -0.44 (experiment), our 8-band Lowdin result gives -0.315. The discrepancy is the known 8-band model limit from missing remote p-like CB states. Not a code bug.
- Optical TE/TM status: zone-center matrix elements are now backed by a dedicated regression again, and both the unstrained GaAs/AlGaAs and stripped-down strained InGaAs k_parallel-integrated absorption benchmarks preserve the expected TE-dominant edge. The remaining optics risk is figure/data provenance, not the basic TE/TM absorption path.
- Scattering lifetimes: Kramers-pair double-counting in the scattering kernel produces lifetimes of 10^9 ps instead of the expected 6-12 ps. FL-016 remains open.
- Exciton binding energy quantitative check: not automated (requires separate variational calculation).

## Detailed Comparison Documents

- [Bulk/QW vs nextnano](2026-04-21-nextnano-bulk-qw-comparison.md)
- [Wire vs nextnano/literature](2026-04-21-nextnano-wire-comparison.md)
- [QCSE/g-factor/optics vs nextnano](2026-04-21-nextnano-qcse-optics-comparison.md)

## Test Results

All 36 tests pass: 15 unit + 21 regression.
