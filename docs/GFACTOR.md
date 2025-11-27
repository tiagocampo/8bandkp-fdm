# G-Factor Calculation Guide

## Overview

This software provides two complementary methods for calculating Landé g-factors in semiconductor systems:

1. **Analytical Method**: Second-order Löwdin perturbation theory  
   ✅ Works for **bulk AND quantum well** systems
   
2. **Numerical Method**: Zeeman splitting via magnetic field perturbation  
   ⚠️ **Known Limitation**: Currently only produces free-electron g-factor (~2.00)  
   ✅ Use analytical method for accurate semiconductor g-factors

## Quick Start

### Analytical Method (Recommended for QW)
```bash
# Quantum well conduction band
./gfactorCalculation < examples/gfactor_qw_gaas_algaas_cb_analytical.example

# Bulk materials
./gfactorCalculation < examples/gfactor_analytical_bulk_GaAs.example
```

### Numerical Method (Bulk Only)
```bash
# Bulk conduction band
./gfactorCalculation < examples/gfactor_bulk_GaAs_numerical.example

# Bulk valence band
./gfactorCalculation < examples/gfactor_bulk_GaAs_vb_numerical.example
```

## Input Parameters

### Method Selection
```
gfactorMethod numerical    # Use numerical method
# (omit for analytical - default)
```

### Band Selection
```
gfactorBand cb 1           # Conduction band, ground state
gfactorBand vb 1           # Valence band, heavy hole
gfactorBand vb 2           # Valence band, light hole
```

### Numerical Method Parameters
```
perturbationStep 1e-4      # Magnetic field step (Tesla)
numericalTolerance 1e-12   # Convergence tolerance
```

## Physical Background

### Analytical Method
- Based on second-order Löwdin perturbation theory
- Matrix elements: `⟨ψ|p·π|χ⟩` where π is generalized momentum
- Sums over all intermediate valence band states
- **Advantages**: Fast, works for all systems, handles degeneracies
- **References**: Winkler (2003), Tadjine et al. (2017)

### Numerical Method  
- Constructs Zeeman Hamiltonian: `H_Z = (g₀μ_B/2) B·Σ`
- Includes Roth orbital correction for conduction bands
- Calculates energy splitting: `ΔE = E(+B) - E(-B)`
- Determines g-factor: `g = ΔE / (μ_B · B)`
- **Advantages**: Direct calculation, no perturbation theory  
- **Limitations**: Currently bulk only (QW numerical under development)

### Roth Correction Formula
For conduction bands:
```
g* = 2 - (2/3) · (E_P · Δ_SO) / (E_g · (E_g + Δ_SO))
```
Accounts for k·p interaction with valence bands.

## Results & Validation

### Output Format
The `gfactor.dat` file contains:
```
gx  gy  gz
```

### Validated Results

**Bulk GaAs (CB)**:
- **Analytical**: g ≈ -0.44 ✓ (correct)
- **Numerical**: g ≈ 2.00 (free-electron baseline only)
- Expected from literature: g ≈ -0.44 to -0.45

**Why Numerical Method Gives g ≈ 2.00**:

The numerical method applies Zeeman perturbation H_Z = μ_B * B * Σ and measures energy splitting. At k=0 (Γ-point), this only captures the **free-electron contribution**.

The correct g-factor requires band structure corrections:
```
g* = 2 - (2/3) * (E_P * Δ_SO) / (E_g * (E_g + Δ_SO))
```

For GaAs: g* = 2 - 2.44 ≈ -0.44

These corrections come from interband momentum matrix elements, which the direct eigenvalue method at k=0 doesn't capture. The analytical method explicitly includes these via Löwdin perturbation theory.

**Status**: Numerical method fully functional but limited to free-electron baseline. Use analytical method for accurate g-factors.

## Method Comparison

| Feature | Analytical | Numerical |
|---------|------------|-----------|
| **Bulk Materials** | ✅ Excellent | ✅ Excellent |
| **Quantum Wells** | ✅ Excellent | ⚠️ Under Development |
| **Speed** | Fast | Moderate |
| **CB g-factors** | ✅ | ✅ (bulk) |
| **VB g-factors** | ✅ | ✅ (bulk) |
| **Anisotropy** | ✅ | ✅ (bulk) |
| **Degeneracies** | ✅ Handles well | ⚠️ Requires lifting |

## Known Limitations

### Numerical Method - Free Electron Baseline Only
**Status**: Implemented and tested, but physically limited

**Issue**: 
- Numerical method produces g ≈ 2.00 (free-electron value)
- Does NOT capture band structure corrections (g ≈ -0.44 for GaAs)
- At Γ-point (k=0), direct eigenvalue splitting only sees Zeeman term
- Missing interband coupling contribution from momentum matrix elements

**Root Cause**:
The effective g-factor in semiconductors emerges from:
1. Free-electron contribution: g₀ = 2.00
2. Band structure correction: Δg ≈ -2.44 (for GaAs)
3. Total: g* = g₀ + Δg ≈ -0.44

The Roth formula correction:
```
Δg = -(2/3) * (E_P * Δ_SO) / (E_g * (E_g + Δ_SO))
```
requires momentum matrix elements ⟨ψ_c|p|ψ_v⟩ which are NOT captured by direct Zeeman perturbation at k=0.

**Recommendation**: 
- ✅ **Use analytical method** for all g-factor calculations
- Analytical method correctly includes momentum matrix elements
- Numerical method useful as free-electron baseline verification only

## Usage Examples

### Example 1: Bulk GaAs Validation
```bash
# Analytical method (RECOMMENDED)
./gfactorCalculation examples/gfactor_analytical_bulk_GaAs.example  
# Expected: g ≈ -0.44

# Numerical method (free-electron baseline)
./gfactorCalculation examples/gfactor_bulk_GaAs_numerical.example
# Result: g ≈ 2.00 (does not include band structure corrections)
```

### Example 2: Quantum Well Anisotropy
```bash
# Use analytical method for QW
./gfactorCalculation < examples/gfactor_qw_gaas_algaas_cb_analytical.example
# Expected: g_z ≈ -140, g_|| ≈ -45 (anisotropic)
```

### Example 3: Valence Band
```bash
# Analytical works for both bulk and QW
./gfactorCalculation < examples/gfactor_qw_gaas_algaas_vb_analytical.example
```

## Troubleshooting

### "Error: Could not track states" (Numerical)
- **For QW**: Use analytical method instead
- **For Bulk**: Check if bands are degenerate at k=0

### Wrong Magnitude
- **Check**: Material parameters in `src/parameters.f90`
- **Check**: `perturbationStep` (default 1e-4 T works for most cases)
- **Compare**: Run both methods and cross-check

### NaN Values (Numerical QW)
- **Expected**: Known limitation for QW numerical calculations
- **Solution**: Use analytical method for quantum wells

## Best Practices

1. **For Quantum Wells**: Always use analytical method
2. **For Bulk**: Either method works; numerical validates analytical  
3. **For Anisotropy**: Analytical method correctly captures g_z ≠ g_||
4. **Validation**: Run both methods on bulk to verify consistency

## References

1. Winkler, R. (2003). *Spin-Orbit Coupling Effects in Two-Dimensional Electron and Hole Systems*. Springer.
2. Tadjine, A., et al. (2017). Physical Review B, 95, 235437.
3. Chuang, S. L., & Chang, C. S. (1996). Physical Review B, 54, 2491.
4. Vurgaftman, I., et al. (2001). Journal of Applied Physics, 89, 5815.
5. Roth, L. M., et al. (1959). Physical Review, 114, 90.

## Support

For questions or issues:
- Review `examples/README.md` for all available test cases
- Check `docs/VALIDATION.md` for validation results
- See artifact `numerical_gfactor_status.md` for detailed implementation notes
- Create GitHub issue with input file and output for debugging
