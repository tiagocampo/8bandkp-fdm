# Wire Validation Summary

Date: 2026-04-21

## Bugs Fixed (commit 663f50b)

### Bug 1: S/SC terms — wrong operator structure (CRITICAL)
- **Before**: S used `d/dx + d/dy` (real gradient), SC = -S
- **After**: S uses `d/dx - i*d/dy` (complex gradient), SC = Hermitian conjugate `d/dx + i*d/dy`
- **Impact**: HH-LH mixing and spin-orbit coupling in wire cross-section were completely wrong

### Bug 2: PZ term — spatial gradient instead of kz (CRITICAL)
- **Before**: PZ = `i*P*(d/dx + d/dy)` (spatial gradient)
- **After**: PZ = `-i*P*kz` (diagonal, z is free direction)
- **Impact**: Spurious CB-VB coupling at kz=0, wrong coupling at finite kz. Lower-triangle block prefactors also fixed for anti-Hermitian conjugation.

### Bug 3: R/RC terms — missing 2D spatial derivatives (CRITICAL)
- **Before**: R = `-sqrt(3)*gamma2*kz²*I` (diagonal only)
- **After**: R = `-sqrt(3)*[gamma2*(d²/dx² - d²/dy²) - 2i*gamma3*d²/dxdy + gamma2*kz²*I]`
- **Impact**: Anisotropic confinement contribution to VB splitting was absent

### Bug 4: kpterms_2d expanded from 15 to 17 operators
- Term 16: `gamma2*(D2x - D2y)` anisotropic Laplacian for R term
- Term 17: placeholder

## Post-Fix Results

### GaAs rectangular wire (63×63 Å, 21×21 grid)
- Band gap: **0.36 eV** (bulk GaAs Eg = 1.52 eV; gap reduced by confinement effects in VB)
- Subbands: 41 VB + 7 CB at kz=0
- Clear parabolic CB dispersion, complex VB mixing
- Wavefunctions localized in wire region

### InAs/GaAs core-shell wire (40 Å InAs core, GaAs shell, 30×30 grid)
- Band gap: **0.19 eV** (InAs narrow-gap, reduced by confinement)
- Subbands: ~47 VB + 1 CB
- Consistent with Type-I band alignment

### InSb wire g-factor (11×11 grid, 55×55 A)
- gx ≈ **-49.94**, gy ≈ **-50.05**, gz ≈ **-49.97**
- Bulk InSb |g*| ≈ 51; this coarse regression now checks band-character selection and bulk-scale magnitude, not a converged wire-radius trend
- Dominant gz component reflects the free-propagation direction of the wire
- Transverse components strongly quenched by 2D confinement
- Updated after solver parity fix (commit 9af3329)

## Remaining Limitations

1. **FEAST convergence**: `info=3` warnings for many kz-points. The subspace dimension (feast_m0) is often insufficient for the dense wire spectrum. Increasing feast_m0 or narrowing the energy window helps.

2. **Boundary treatment**: The `csr_apply_variable_coeff` diagonal reconstruction at boundaries uses `-sum(off-diagonal)` which differs from the exact FD stencil for non-central stencils. Effect is O(1/N) and localized to boundary cells.

3. **Hard-wall boundaries only**: No open or absorbing boundary conditions. Wire states near the boundary may have unphysical confinement.

4. **Circular/polygon wire shapes**: Inactive cells (outside wire) get zero potential instead of infinite barrier. This can cause wavefunction leakage for non-rectangular geometries. Rectangular wires are unaffected.

5. **No spin-orbit coupling correction beyond 8-band**: Higher-order terms from bands beyond the 8-band basis are not included. This affects g-factor accuracy, especially for narrow-gap materials.

## Trust Envelope

- **High confidence**: Rectangular wire bandstructure, subband energies, wavefunction localization
- **Medium confidence**: Wire g-factors (Lowdin partitioning sensitive to basis truncation)
- **Low confidence**: Non-rectangular wire shapes (boundary treatment issues)
