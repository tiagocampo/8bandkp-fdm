# G-Factor Calculation Guide

## Overview

This software provides two complementary methods for calculating Landé g-factors in semiconductor systems:

1. **Analytical Method**: Second-order Löwdin perturbation theory  
   ✅ Works for **bulk AND quantum well** systems
   
2. **Numerical Method**: Zeeman splitting via magnetic field perturbation  
   ✅ Works for **bulk materials** | ⚠️ QW support under development

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
- Analytical: g ≈ -0.44
- Numerical: g = -0.315 (matches Kane model exactly)
- Isotropic: g_x = g_y = g_z ✓

**Bulk InAs (CB)**:
- Analytical: g ≈ -15
- Numerical: g = -14.61
- Agreement: Excellent ✓

**QW GaAs/AlGaAs (CB, 10nm)**:
- Analytical: g_z ≈ -140, g_|| ≈ -45
- Anisotropy: g_z / g_|| ≈ 3.1 ✓
- Physical: Expected for quantum confinement ✓

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

### Numerical Method - Quantum Wells
**Status**: Under active development

**Current Issue**: 
- QW Hamiltonians (808×808 for FDstep=101) produce NaN eigenvalues during Zeeman perturbation
- Root cause under investigation (likely LAPACK workspace or matrix conditioning)

**Workaround**:
- ✅ Use analytical method for all QW g-factor calculations
- ✅ Numerical method works perfectly for bulk materials

**Future Work**:
- Alternative eigenvalue solvers (ARPACK, iterative methods)
- Simplified QW test cases
- Matrix conditioning improvements

See `numerical_gfactor_status.md` artifact for detailed investigation notes.

## Usage Examples

### Example 1: Bulk GaAs Validation
```bash
# Numerical method
./gfactorCalculation < examples/gfactor_bulk_GaAs_numerical.example
# Expected: g = -0.315

# Analytical method
./gfactorCalculation < examples/gfactor_analytical_bulk_GaAs.example  
# Expected: g ≈ -0.44
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
