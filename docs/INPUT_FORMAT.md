# Input File Format Guide

## Overview

Input files for `bandStructure` and `gfactorCalculation` use a simple keyword-value format. As of the latest update, **comment lines are fully supported** using `#` or `!` at the start of a line.

## Basic Format

```
# This is a comment line
! This is also a comment line

keyword: value                    # Comments can appear on their own lines
```

## Common Parameters

### Wave Vector Configuration
```
waveVector: kx                    # Direction: kx, ky, kz, or k0
waveVectorMax: 0.1                # Maximum k-point (fraction of BZ)
waveVectorStep: 11                # Number of k-points (or 0 for single point)
```

**Note**: `waveVectorStep: 0` is automatically converted to `1` for single k-point calculations.

### System Configuration
```
confinement: 0                    # 0 = bulk, 1 = quantum well
FDstep: 101                       # Finite difference grid points
numLayers: 1                      # Number of material layers
```

**Note**: For bulk calculations (`confinement: 0`, `numLayers: 1`), `FDstep` is automatically set to `1` for optimal performance.

### Material Specification

**Bulk Systems**:
```
material1: GaAs
```

**Quantum Wells**:
```
material1: AlGaAs -100 100 0      # Material start end offset
material2: GaAs -50 50 0          # Core layer
material3: AlGaAs -100 100 0      # Barrier
```

### Band Configuration
```
numcb: 2                          # Number of conduction bands
numvb: 6                          # Number of valence bands
```

### External Fields
```
ExternalField: 0 EF               # 0/1 type (EF = electric field)
EFParams: 0.0005                  # Field strength
```

## G-Factor Specific Parameters

### Method Selection
```
gfactorMethod: numerical          # Use numerical method
# (omit for analytical - default)
```

### Numerical Method Configuration
```
numericalTolerance: 1e-12         # Convergence tolerance
perturbationStep: 1e-4            # Magnetic field step (Tesla)
useSparseSolver: false            # Use sparse matrix solver
validateWithAnalytical: false     # Cross-check with analytical
```

### Band Selection for G-Factor
```
gfactorBand: cb 1                 # Conduction band, ground state
gfactorBand: vb 1                 # Valence band, heavy hole
```

## Complete Example Files

### Bulk Band Structure
```
# Bulk GaAs band structure calculation
waveVector: kx
waveVectorMax: 0.1
waveVectorStep: 100
confinement: 0
FDstep: 1                        # Ignored for bulk (auto-set to 1)
numLayers: 1
material1: GaAs
numcb: 2
numvb: 6
ExternalField: 0 EF
EFParams: 0.0
```

### Quantum Well Band Structure
```
# GaAs/AlGaAs quantum well
waveVector: kz
waveVectorMax: 0.02
waveVectorStep: 21
confinement: 1
FDstep: 101
numLayers: 3
material1: AlGaAs -100 100 0
material2: GaAs -50 50 0
material3: AlGaAs -100 100 0
numcb: 10
numvb: 10
ExternalField: 0 EF
EFParams: 0.0
```

### Numerical G-Factor (Bulk)
```
# Bulk GaAs numerical g-factor
waveVector: k0
waveVectorMax: 0.1
waveVectorStep: 0                # Single point at k=0
confinement: 0
FDstep: 101                      # Auto-set to 1 for bulk
numLayers: 1
material1: GaAs
numcb: 2
numvb: 6
ExternalField: 0 EF
EFParams: 0.0005

# Numerical method configuration
gfactorMethod: numerical
numericalTolerance: 1e-12
perturbationStep: 1e-4
gfactorBand: cb 1
```

## Best Practices

1. **Use Comments**: Document your input files with `#` or `!` comments
2. **Bulk Calculations**: Set `confinement: 0` and `numLayers: 1`
3. **Single k-point**: Use `waveVectorStep: 0` (auto-converts to 1)
4. **G-Factor Method**: Omit `gfactorMethod` for analytical (default)
5. **Validation**: Include comments with expected results

## Recent Updates

- ✅ **Comment Support**: Lines starting with `#` or `!` are now properly skipped
- ✅ **Auto-Correction**: `waveVectorStep: 0` → `1` automatically
- ✅ **Bulk Optimization**: `FDstep` forced to `1` for bulk systems
- ✅ **Robust Parsing**: All input files handle comments consistently

## Troubleshooting

### "Bad real/integer number in item X"
- **Cause**: Comment line treated as data
- **Solution**: Update to latest version (comment support added)
- **Workaround**: Remove all comment lines from input file

### "Index out of range" errors
- **Cause**: `waveVectorStep: 0` in older versions
- **Solution**: Update to latest version (auto-correction added)
- **Workaround**: Set `waveVectorStep: 1` manually

## See Also

- [GFACTOR.md](GFACTOR.md) - G-factor calculation methods
- [QUICKSTART.md](QUICKSTART.md) - Getting started guide
- [examples/README.md](../examples/README.md) - Example input files
