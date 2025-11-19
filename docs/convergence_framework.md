# Numerical Perturbation Convergence Framework

This document describes the convergence testing and step size optimization framework implemented for the numerical perturbation method in g-factor calculations.

## Overview

The convergence framework addresses numerical precision issues in g-factor calculations by implementing adaptive step sizing and convergence testing. This ensures that calculated g-factors are numerically stable and accurate within specified tolerances.

## Key Components

### 1. Configuration Structures

#### `ConvergenceConfig`
Controls the convergence testing behavior:

```fortran
type :: ConvergenceConfig
  real(kind=dp) :: energy_tolerance = 1e-6_dp        ! Energy difference < 1e-6 eV
  real(kind=dp) :: gfactor_tolerance = 1e-3_dp      ! g-factor change < 0.1%
  real(kind=dp) :: min_step_size = 1e-12_dp         ! Minimum perturbation step
  real(kind=dp) :: step_reduction_factor = 0.1_dp   ! Reduce step by factor of 10
  integer :: max_iterations = 10                    ! Maximum convergence iterations
  logical :: enable_adaptive_stepping = .true.      ! Enable adaptive step sizing
  logical :: verbose_convergence = .false.          ! Detailed convergence output
end type ConvergenceConfig
```

**Parameters:**
- `energy_tolerance`: Maximum allowed energy difference between iterations (default: 1e-6 eV)
- `gfactor_tolerance`: Maximum allowed relative g-factor change (default: 0.1% = 1e-3)
- `min_step_size`: Smallest allowed perturbation step size (default: 1e-12)
- `step_reduction_factor`: Factor by which step size is reduced (default: 0.1)
- `max_iterations`: Maximum number of convergence iterations (default: 10)
- `enable_adaptive_stepping`: Whether to automatically reduce step size (default: .true.)
- `verbose_convergence`: Whether to output detailed convergence information (default: .false.)

#### `ConvergenceResults`
Stores convergence metadata and results:

```fortran
type :: ConvergenceResults
  logical :: converged = .false.
  integer :: iterations_used = 0
  real(kind=dp) :: final_step_size = 0.0_dp
  real(kind=dp) :: final_energy_diff = 0.0_dp
  real(kind=dp) :: final_gfactor_diff = 0.0_dp
  real(kind=dp) :: initial_gfactor = 0.0_dp
  real(kind=dp) :: final_gfactor = 0.0_dp
  character(len=128) :: convergence_status = "Not started"
  logical :: numerical_stability_achieved = .false.
end type ConvergenceResults
```

### 2. Main Subroutine

#### `gfactorCalculationNumericalConverged`
Enhanced g-factor calculation with convergence testing:

```fortran
subroutine gfactorCalculationNumericalConverged(tensor, whichBand, bandIdx, numcb, numvb, &
  & cb_state, vb_state, cb_value, vb_value, nlayers, params, startz, endz, &
  & profile, kpterms, dz, num_config, conv_config, conv_results)
```

**Parameters:**
- `tensor`: Output g-factor tensor (2×2×3 complex array)
- `num_config`: Numerical configuration structure
- `conv_config`: Convergence configuration structure
- `conv_results`: Output convergence results structure

### 3. Helper Subroutines

#### `test_convergence_criteria`
Tests if convergence criteria are met:

```fortran
subroutine test_convergence_criteria(energy_diff, gfactor_diff, energy_tol, gfactor_tol, converged, stability_ok)
```

#### `output_convergence_metadata`
Outputs detailed convergence information:

```fortran
subroutine output_convergence_metadata(conv_results, field_dir, output_unit)
```

## Convergence Algorithm

### Step-by-Step Process

1. **Initialization**: Start with user-specified perturbation step size
2. **Hamiltonian Construction**: Build perturbed Hamiltonians with ±δB
3. **Eigenvalue Solution**: Solve eigenvalue problems for both Hamiltonians
4. **Energy Extraction**: Extract target state energies from both solutions
5. **G-factor Calculation**: Compute g-factor using central finite difference
6. **Convergence Testing**: Compare with previous iteration results
7. **Adaptive Stepping**: If not converged, reduce step size by reduction factor
8. **Termination**: Stop when converged or minimum step size reached

### Convergence Criteria

The framework tests two primary criteria:

1. **Energy Convergence**: `|ΔE_curr - ΔE_prev| < energy_tolerance`
2. **G-factor Convergence**: `|g_curr - g_prev|/|g_prev| < gfactor_tolerance`

Both criteria must be satisfied for convergence to be achieved.

### Step Size Adaptation

- **Automatic Reduction**: Step size reduced by `step_reduction_factor` (default 0.1)
- **Minimum Limit**: Stops when step size reaches `min_step_size` (default 1e-12)
- **User Control**: Can be disabled via `enable_adaptive_stepping = .false.`

## Usage Example

```fortran
! Configuration setup
type(NumericalConfig) :: num_config
type(ConvergenceConfig) :: conv_config
type(ConvergenceResults) :: conv_results

! Configure numerical parameters
num_config%perturbation_step = 1e-3_dp
num_config%verbose_output = .true.

! Configure convergence parameters
conv_config%energy_tolerance = 1e-6_dp
conv_config%gfactor_tolerance = 1e-3_dp
conv_config%max_iterations = 8
conv_config%verbose_convergence = .true.

! Run converged calculation
call gfactorCalculationNumericalConverged(tensor, whichBand, bandIdx, numcb, numvb, &
  & cb_state, vb_state, cb_value, vb_value, nlayers, params, startz, endz, &
  & dz, num_config, conv_config, conv_results)

! Check results
if (conv_results%converged) then
  print *, 'Converged after', conv_results%iterations_used, 'iterations'
  print *, 'Final g-factor:', conv_results%final_gfactor
else
  print *, 'Convergence failed:', trim(conv_results%convergence_status)
end if
```

## Output and Metadata

The framework provides comprehensive metadata:

- **Convergence Status**: Whether calculation converged and why
- **Iteration Count**: Number of iterations required
- **Final Step Size**: Optimal perturbation step size found
- **Energy Differences**: Final energy convergence metrics
- **G-factor Evolution**: Initial and final g-factor values
- **Stability Information**: Numerical stability verification

### Sample Output

```
=== Converged Numerical G-Factor Calculation ===
Target state index: 7
Initial perturbation step: 1.000E-03
Energy tolerance: 1.000E-06 eV
G-factor tolerance: 1.000E-03
Minimum step size: 1.000E-12
Maximum iterations: 10

Computing converged g-factor for direction: x
--- Iteration 1 for direction x ---
Current step size: 1.000E-03
E_plus: 0.500000 E_minus: 0.499000
Calculated g-factor: 5.123456

--- Iteration 2 for direction x ---
Current step size: 1.000E-04
Energy difference: 1.234E-07 eV (tolerance: 1.000E-06)
G-factor change: 2.345E-04 (tolerance: 1.000E-03)
Convergence achieved after 2 iterations

=== Converged Numerical G-Factor Calculation Complete ===
Final convergence status: Converged for direction x
Total iterations used: 2
Final step size: 1.000E-04
Final g-factor range: 5.123456 to 5.123456
```

## Error Handling and Robustness

### Numerical Stability

- **NaN Detection**: Checks for NaN and infinite values
- **Solver Failure Handling**: Graceful handling of eigenvalue solver failures
- **Step Size Limits**: Prevents numerical underflow/overflow

### Convergence Failure Modes

1. **Solver Failure**: Eigenvalue solver fails → step reduction
2. **Minimum Step Size**: Convergence not achieved before minimum step
3. **Maximum Iterations**: Convergence not achieved within iteration limit
4. **Numerical Instability**: Invalid numerical values detected

### Recovery Strategies

- **Automatic Step Reduction**: When solver fails or convergence issues detected
- **Fallback Results**: Uses best available result if convergence fails
- **Detailed Status Reporting**: Clear indication of failure reasons

## Integration with Existing Code

The convergence framework is designed to integrate seamlessly with existing code:

- **Backward Compatibility**: Original `gfactorCalculationNumerical` remains unchanged
- **Optional Usage**: Can be used selectively when convergence is needed
- **Same Interface**: Uses same input parameters as original numerical method
- **Metadata Output**: Provides additional information without changing core results

## Performance Considerations

- **Computational Cost**: Each iteration doubles the computational effort
- **Memory Usage**: Moderate additional memory for convergence testing
- **Typical Iterations**: 2-4 iterations usually sufficient for convergence
- **Speed vs Accuracy**: Trade-off controlled by tolerance settings

## Recommendations

### Tolerance Settings

- **High Accuracy**: `energy_tolerance = 1e-8`, `gfactor_tolerance = 1e-4`
- **Standard Use**: `energy_tolerance = 1e-6`, `gfactor_tolerance = 1e-3` (default)
- **Fast Computation**: `energy_tolerance = 1e-4`, `gfactor_tolerance = 1e-2`

### Step Size Strategy

- **Initial Step**: 1e-3 to 1e-4 typically works well
- **Minimum Step**: 1e-12 to 1e-14 for double precision
- **Reduction Factor**: 0.1 (factor of 10) recommended

### Usage Guidelines

1. **First Calculation**: Use default settings to establish baseline
2. **Convergence Testing**: Enable verbose output to understand behavior
3. **Performance Tuning**: Adjust tolerances based on accuracy requirements
4. **Validation**: Compare with analytical results when available

## Troubleshooting

### Common Issues

1. **No Convergence**: Check minimum step size and maximum iterations
2. **Large g-factors**: May indicate numerical precision issues
3. **Slow Convergence**: Consider adjusting initial step size
4. **Solver Failures**: Often resolved by step size reduction

### Debugging Options

- Enable `verbose_convergence` for detailed iteration information
- Check `convergence_status` for specific failure reasons
- Monitor `final_energy_diff` and `final_gfactor_diff` values
- Verify `numerical_stability_achieved` flag

This framework provides robust convergence testing for numerical g-factor calculations while maintaining compatibility with existing code and providing comprehensive diagnostics for troubleshooting and optimization.