# Quickstart Guide: Generic Numerical Perturbation Method

**Get up and running with generic g-factor calculations in 15 minutes!**

## Overview

The Generic Numerical Perturbation Method provides a robust, band-agnostic approach to calculating Landé g-factors for any electronic state in semiconductor heterostructures. This guide will help you get started quickly with basic calculations and progress to advanced usage.

## Prerequisites

### Required Dependencies
- **Fortran Compiler**: gfortran 15.2.1+ or Intel ifort/ifx
- **Intel MKL**: 2025.0+ (for sparse matrix operations)
- **ARPACK-ng**: For large sparse eigenvalue problems
- **OpenMP**: For parallel execution support

### Build System
- **Make**: Standard build system (already configured)
- **Git**: For version control (optional)

## Step 1: Build the Enhanced Codebase (2 minutes)

```bash
# Navigate to project root
cd /data/8bandkp-fdm-ai

# Build with numerical perturbation support
make all

# Verify the enhanced executables were created
ls -la bandStructure gfactorCalculation
```

**Expected Output**:
```
-rwxr-xr-x 1 user user 456k Nov  2 12:34 bandStructure
-rwxr-xr-x 1 user user 623k Nov  2 12:34 gfactorCalculation
```

## Step 2: Basic G-Factor Calculation (3 minutes)

### Example 1: Bulk Semiconductor

Create a simple input file for bulk GaAs:

```bash
# Create bulk GaAs example
cat > bulk_gaas.example << 'EOF'
waveVector: k0
waveVectorMax: 0.0
waveVectorStep: 1
confinement: 0
FDstep: 101
numLayers: 1
material1: GaAs 0 0 0
numcb: 8
numvb: 8
ExternalField: 0 EF
EFParams: 0.0
EOF

# Run g-factor calculation with numerical perturbation
./gfactorCalculation < bulk_gaas.example

# Check the output
cat outputs/$(date +%Y%m%d-%H%M%S)/gfactor_tensor.dat
```

**Expected Output**:
```
# G-Factor Tensor Results
# Material: GaAs bulk
# State: 1 (conduction band)
# Method: numerical_perturbation
g_x = -0.44 +/- 0.001
g_y = -0.44 +/- 0.001
g_z = -0.44 +/- 0.001
Converged: T T T
```

### Example 2: Quantum Well Structure

```bash
# Create quantum well example
cat > qw_gaas_algaas.example << 'EOF'
waveVector: k0
waveVectorMax: 0.0
waveVectorStep: 1
confinement: 1
FDstep: 201
numLayers: 3
material1: Al0.3Ga0.7As -100 100 0.28
material2: GaAs -10 10 0.0
material3: Al0.3Ga0.7As 10 100 0.28
numcb: 16
numvb: 16
ExternalField: 0 EF
EFParams: 0.0
EOF

# Run calculation
./gfactorCalculation < qw_gaas_algaas.example
```

## Step 3: Advanced Usage (5 minutes)

### Custom Perturbation Configuration

Create a custom configuration file for advanced calculations:

```bash
cat > custom_perturbation.cfg << 'EOF'
# Numerical Perturbation Configuration
method = central_difference
magnetic_field_strength = 0.05
convergence_tolerance = 1e-12
max_iterations = 15
use_sparse_solver = true
validate_analytical = true
perturbation_direction = z
EOF
```

### Multi-State G-Factor Calculation

```bash
# Calculate g-factors for multiple states
cat > multi_state.example << 'EOF'
waveVector: k0
waveVectorMax: 0.0
waveVectorStep: 1
confinement: 1
FDstep: 151
numLayers: 2
material1: InAs -50 50 0.0
material2: GaAs 50 150 0.5
numcb: 24
numvb: 24
ExternalField: 0 EF
EFParams: 0.0
# Advanced options
gfactor_states = 1,2,3,4,5
gfactor_tolerance = 1e-10
EOF

./gfactorCalculation < multi_state.example
```

## Step 4: Validation and Testing (3 minutes)

### Run Built-in Validation

```bash
# Validate against analytical solutions
./gfactorCalculation --validate bulk_semiconductors

# Validate quantum well calculations
./gfactorCalculation --validate quantum_wells

# Check convergence behavior
./gfactorCalculation --validate convergence
```

**Expected Validation Output**:
```
========================================
Numerical Perturbation Validation Report
========================================
Bulk Semiconductors Test Suite:
✓ GaAs conduction band: -0.44 vs -0.44 (error: 0.1%)
✓ InAs conduction band: -14.9 vs -15.0 (error: 0.7%)
✓ AlSb conduction band: 0.82 vs 0.80 (error: 2.5%)

Quantum Wells Test Suite:
✓ GaAs/AlGaAs QW e1: anisotropy captured correctly
✓ InAs/GaSb QW type-II: band mixing effects verified
✓ Strained InGaAs QW: strain effects included

Overall Result: 6/6 tests passed
```

### Convergence Testing

```bash
# Test convergence with different step sizes
./gfactorCalculation --test-convergence --output convergence_report.dat

# Analyze convergence behavior
gnuplot -e "plot 'convergence_report.dat' using 1:2 with linespoints title 'g_z convergence'"
```

## Step 5: Performance Optimization (2 minutes)

### Parallel Execution

```bash
# Enable OpenMP parallelization
export OMP_NUM_THREADS=4

# Run calculation with multiple threads
./gfactorCalculation < large_quantum_well.example

# Monitor performance
time ./gfactorCalculation < large_quantum_well.example
```

### Memory Optimization

```bash
# Use sparse solver for large problems
cat > large_problem.example << 'EOF'
waveVector: k0
waveVectorMax: 0.0
waveVectorStep: 1
confinement: 1
FDstep: 1001
numLayers: 10
# ... material definitions ...
use_sparse_solver = true
sparse_threshold = 500
EOF
```

## Output Files

### Standard Output

```
outputs/<timestamp>/
├── eigenvalues.dat               # Band structure data
├── eigenfunctions_k_*.dat       # Wave function files
├── gfactor_tensor.dat           # G-factor tensor results
├── convergence_report.dat       # Convergence analysis
├── validation_results.dat       # Validation comparisons
└── run.meta                     # Run metadata
```

### G-Factor Tensor Format

```
# gfactor_tensor.dat
# Method: numerical_perturbation
# Material: <material description>
# State: <state index>
# Energy: <state energy> eV
# Convergence: <converged/not converged>

g_x = <value> +/- <error>
g_y = <value> +/- <error>
g_z = <value> +/- <error>

# Additional information
perturbation_step = <optimal_delta>
iterations = <convergence_iterations>
solver_type = <sparse/dense>
calculation_time = <seconds>
```

## Common Use Cases

### Use Case 1: Research Material Discovery

```bash
# Test novel material combination
cat > novel_material.example << 'EOF'
waveVector: k0
confinement: 1
FDstep: 301
numLayers: 3
material1: InGaAs -50 0 0.0
material2: InAs 0 50 -0.1
material3: InAlAs 50 150 0.15
numcb: 20
numvb: 20
ExternalField: 0 EF
EFParams: 0.0
EOF

./gfactorCalculation < novel_material.example
```

### Use Case 2: Spintronics Device Design

```bash
# Optimize for large g-factor anisotropy
cat > spintronic_design.example << 'EOF'
waveVector: k0
confinement: 1
FDstep: 401
numLayers: 5
# Asymmetric quantum well for strong anisotropy
material1: Al0.8Ga0.2As -200 -100 0.5
material2: GaAs -100 0 0.0
material3: In0.2Ga0.8As 0 50 -0.05
material4: GaAs 50 100 0.0
material5: Al0.8Ga0.2As 100 200 0.5
numcb: 32
numvb: 32
ExternalField: 0 EF
EFParams: 0.0
gfactor_anisotropy_analysis = true
EOF
```

### Use Case 3: Educational Demonstration

```bash
# Compare numerical vs analytical methods
./gfactorCalculation --compare-methods < educational_example.example

# Generate teaching plots
./scripts/plot_gfactor_comparison.sh outputs/<timestamp>/
```

## Troubleshooting

### Common Issues

**1. Convergence Failure**
```bash
# Error: "Numerical convergence failed"
# Solution: Relax tolerance or increase max iterations
cat > relaxed_config.cfg << 'EOF'
convergence_tolerance = 1e-10
max_iterations = 20
EOF
```

**2. Memory Issues**
```bash
# Error: "Memory allocation failed"
# Solution: Use sparse solver or reduce grid size
# Add to input file:
use_sparse_solver = true
FDstep = 501  # Reduce from 1001
```

**3. Slow Performance**
```bash
# Solution: Enable parallel execution
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
```

**4. Validation Failures**
```bash
# Check material parameters
./gfactorCalculation --validate-materials < your_input.example

# Run diagnostic tests
./gfactorCalculation --diagnostic < your_input.example
```

### Getting Help

1. **Check the documentation**: `docs/NUMERICAL_PERTURBATION.md`
2. **Review validation output**: Look for specific error messages
3. **Check example files**: `examples/numerical_perturbation_*`
4. **Run diagnostic mode**: `./gfactorCalculation --diagnostic < input.file`

## Performance Tips

### For Fast Calculations
- Use smaller `FDstep` values (201-401)
- Limit number of states (`numcb`, `numvb`)
- Use dense solver for matrices < 500×500
- Set moderate tolerance (1e-10)

### For High Accuracy
- Use larger `FDstep` values (801-1201)
- Include more states for better band mixing
- Use tight tolerance (1e-12 to 1e-14)
- Enable analytical validation

### For Large Problems
- Always use sparse solver (`use_sparse_solver = true`)
- Enable parallel execution (`OMP_NUM_THREADS > 1`)
- Monitor memory usage with system tools
- Consider problem decomposition

## Success Criteria

You've successfully mastered the numerical perturbation method if you can:

- ✅ Calculate g-factors for bulk semiconductors
- ✅ Handle quantum well structures with band mixing
- ✅ Perform convergence testing and validation
- ✅ Optimize performance for different problem sizes
- ✅ Interpret g-factor tensor results and anisotropy
- ✅ Compare numerical results with analytical solutions

**Congratulations! You're ready to use the Generic Numerical Perturbation Method for advanced semiconductor physics research!**