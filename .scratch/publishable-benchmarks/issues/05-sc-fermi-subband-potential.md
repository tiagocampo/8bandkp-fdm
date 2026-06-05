# SC Fermi Level + Subband Shift + Potential Profile — Quantitative SC Benchmarks

**Type:** AFK
**Blocked by:** #02 (SC charge neutrality hard check)
**GitHub:** #33

## What to Build

Extend the SC charge neutrality tracer bullet with three additional quantitative checks:

1. **Fermi level vs analytical parabolic:** For bulk n-GaAs at doping n = 5×10¹⁸ cm⁻³, the analytical parabolic Fermi level is `E_F = E_C + (ℏ²/2m*)(3π²n)^(2/3)` where m* = 0.067m₀. For the QW case, compare the SC Fermi level against this bulk value as a sanity bound (QW Fermi level will differ due to subband quantization, but should be within a factor of ~2). Assert the SC Fermi level is within a physically justified range.

2. **CB1 subband shift vs flat-band:** Run `bandStructure` twice — once with SC enabled, once flat-band (no `[sc]` section). Extract CB1 eigenvalue at k_∥=0 from both. The shift `ΔE₁ = E₁(SC) - E₁(flat)` should match the Bastard first-order perturbative estimate to within ~20%. The Bastard estimate is: `ΔE ≈ -e²N_D × L_well² / (2ε₀ε_r)` (uniform doping, infinite barrier limit).

3. **Potential profile shape check:** Parse the Hartree potential output. Verify: (a) potential has a V-shaped dip centered in the well (not flat, not linear), (b) potential is approximately flat in the barrier regions, (c) total potential swing `ΔV = V(center) - V(barrier)` is on the order of the expected Hartree energy `e²N_D × L_well² / (2ε₀ε_r)`. Generate a figure showing the potential profile with the expected Hartree scale overlaid.

Register `subband_shift` and `potential_profile` as observable cells in `validation_universe.yml`. Update COVERAGE annotations.

## Acceptance Criteria

- [ ] Fermi level assertion added: SC Fermi level within factor of 2 of analytical parabolic value
- [ ] Subband shift assertion added: ΔE₁ matches Bastard estimate within ~20%
- [ ] Potential profile shape check: V-dip in well, flat barriers, correct energy scale
- [ ] Figure generated: Hartree potential profile with expected scale annotation
- [ ] `subband_shift` and `potential_profile` observable cells added to `validation_universe.yml` as required
- [ ] COVERAGE annotations updated
- [ ] All existing tests still pass

## Blocked by

- #02 (SC charge neutrality tracer bullet must exist first)

## User Stories Covered

- #10: Fermi level vs analytical parabolic approximation
- #11: CB1 subband shift vs Bastard perturbative estimate
- #12: Potential profile shape verification (V-dip in well, flat barriers)
- #13: Potential energy scale verification
- #14: SC benchmark catches regressions
- #15: subband_shift + potential_profile registered in coverage matrix
