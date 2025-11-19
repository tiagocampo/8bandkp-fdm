# API Contract: Generic Numerical Perturbation Method

**Version**: 1.0.0
**Date**: 2025-11-02
**Feature**: Generic Numerical Perturbation Method for G-Factor Calculations

## Overview

This API contract defines the interfaces for the generic numerical perturbation functionality that will be integrated into the existing k·p solver. The contract specifies function signatures, input/output formats, and error handling requirements for the new g-factor calculation capabilities.

## Module Interfaces

### 1. NumericalPerturbation Module

#### Public Interface

```fortran
module numerical_perturbation
  use iso_fortran_env, only: dp => real64
  use defs, only: paramStruct, wavevector
  implicit none

  ! Public types
  public :: PerturbationConfig
  public :: GFactorTensor
  public :: HamiltonianPerturbation
  public :: ConvergenceResult
  public :: ValidationTestCase

  ! Public procedures
  public :: initialize_perturbation_config
  public :: calculate_gfactor_tensor
  public :: validate_numerical_method
  public :: run_convergence_test
  public :: cleanup_perturbation_resources

contains

  ! Initialize perturbation configuration with defaults
  subroutine initialize_perturbation_config(config, material_params, k_vec, target_state)
    type(PerturbationConfig), intent(out) :: config
    type(paramStruct), intent(in) :: material_params
    type(wavevector), intent(in) :: k_vec
    integer, intent(in) :: target_state
  end subroutine

  ! Calculate complete g-tensor for specified state
  subroutine calculate_gfactor_tensor(config, g_tensor, convergence_info)
    type(PerturbationConfig), intent(in) :: config
    type(GFactorTensor), intent(out) :: g_tensor
    type(ConvergenceResult), intent(out) :: convergence_info
  end subroutine

  ! Validate numerical method against analytical solutions
  subroutine validate_numerical_method(test_cases, validation_results)
    type(ValidationTestCase), intent(in) :: test_cases(:)
    logical, intent(out) :: validation_results(size(test_cases))
  end subroutine

  ! Run step size convergence test
  subroutine run_convergence_test(config, direction, convergence_result)
    type(PerturbationConfig), intent(in) :: config
    character(len=1), intent(in) :: direction
    type(ConvergenceResult), intent(out) :: convergence_result
  end subroutine

  ! Clean up allocated resources
  subroutine cleanup_perturbation_resources()
  end subroutine

end module numerical_perturbation
```

#### Type Definitions

```fortran
! Configuration for numerical perturbation calculations
type :: PerturbationConfig
  ! Target specification
  integer :: target_state_index = 1
  real(kind=dp) :: k_point(3) = [0.0_dp, 0.0_dp, 0.0_dp]

  ! Magnetic field parameters
  real(kind=dp) :: magnetic_field_direction(3) = [0.0_dp, 0.0_dp, 1.0_dp]
  real(kind=dp) :: base_field_strength = 0.1_dp
  real(kind=dp) :: perturbation_step = 1e-4_dp

  ! Numerical method configuration
  character(len=32) :: numerical_method = 'central_difference'
  real(kind=dp) :: convergence_tolerance = 1e-12_dp
  integer :: max_iterations = 10
  logical :: use_sparse_solver = .true.
  logical :: auto_step_size = .true.

  ! Validation options
  logical :: enable_validation = .true.
  logical :: verbose_output = .false.

contains
  procedure :: validate => config_validate
  procedure :: to_string => config_to_string
end type

! Complete g-tensor result
type :: GFactorTensor
  real(kind=dp) :: components(3) = [0.0_dp, 0.0_dp, 0.0_dp]  ! gx, gy, gz
  real(kind=dp) :: uncertainties(3) = [0.0_dp, 0.0_dp, 0.0_dp]
  logical :: converged(3) = [.false., .false., .false.]
  real(kind=dp) :: target_energy = 0.0_dp
  character(len=32) :: calculation_method = 'unknown'
  real(kind=dp) :: calculation_time = 0.0_dp
  integer :: error_code = 0
  character(len=256) :: error_message = ''
contains
  procedure :: is_valid => tensor_is_valid
  procedure :: to_string => tensor_to_string
end type

! Hamiltonian with magnetic field perturbations
type :: HamiltonianPerturbation
  complex(kind=dp), allocatable :: H0(:,:)        ! Base Hamiltonian
  complex(kind=dp), allocatable :: H_plus(:,:)    ! +δB perturbation
  complex(kind=dp), allocatable :: H_minus(:,:)   ! -δB perturbation
  real(kind=dp) :: delta_B = 0.0_dp
  character(len=1) :: perturbation_direction = 'z'
  integer :: matrix_size = 0
  logical :: is_hermitian = .true.
contains
  procedure :: initialize => hamiltonian_init
  procedure :: cleanup => hamiltonian_cleanup
end type

! Convergence test results
type :: ConvergenceResult
  logical :: converged = .false.
  real(kind=dp) :: optimal_step_size = 0.0_dp
  real(kind=dp) :: optimal_g_factor = 0.0_dp
  real(kind=dp) :: error_estimate = 0.0_dp
  integer :: iterations_used = 0
  real(kind=dp) :: test_steps(4) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
  real(kind=dp) :: test_estimates(4) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
  character(len=128) :: status_message = ''
contains
  procedure :: is_successful => convergence_is_successful
end type
```

### 2. Input/Output Format Contract

#### Input File Extensions

Existing input configuration files will be extended with g-factor parameters:

```
# Existing parameters (unchanged)
waveVector: 0 0 1
waveVectorMax: 0.02
waveVectorStep: 10
confinement: 1
FDstep: 201
numLayers: 3
material1: AlSb -250 250 0.0
material2: InAs 0 50 0.0
material3: AlSb 50 300 0.0

# New g-factor calculation section
GFactorCalculation: 1                    ! Enable g-factor calculations (0/1)
GFactorTargetState: 1                     ! Target electronic state
GFactorMagneticField: 0 0 0.1            ! Magnetic field direction & magnitude [T]
GFactorNumericalMethod: central_difference ! Numerical method
GFactorStepSize: 1e-4                    ! Initial perturbation step [T]
GFactorTolerance: 1e-12                  ! Convergence tolerance
GFactorMaxIterations: 10                 ! Maximum refinement iterations
GFactorAutoStepSize: 1                   ! Automatic step size optimization
GFactorValidation: 1                     ! Enable validation
GFactorVerbose: 0                        ! Verbose output
```

#### Output File Extensions

##### 1. Extended eigenvalues.dat

```
# kx [1/Å]    ky [1/Å]    kz [1/Å]    E1 [eV]    E2 [eV]    ...    EN [eV]    g1x    g1y    g1z    g2x    g2y    g2z    ...
0.0000       0.0000      0.0000      0.7542     0.8765     ...    2.1345     -0.442  -0.441  -0.443  2.156   1.823   0.542   ...
0.0020       0.0000      0.0000      0.7589     0.8812     ...    2.1412     -0.438  -0.440  -0.442  2.148   1.819   0.539   ...
...
```

##### 2. New gfactor_results.dat

```
# G-Factor Calculation Results
# Generated: 2025-11-02 15:30:45
# Target State: 1
# Magnetic Field: (0.000, 0.000, 0.100) T
# Numerical Method: central_difference
# Auto Step Size: YES
# Convergence: YES

# Component    g-factor    Uncertainty    Converged    Step_Size [T]    Iterations    Time [s]
gx             -0.44213    1.2e-05        YES          1.0e-04           3             0.023
gy             -0.44189    1.3e-05        YES          1.0e-04           3             0.022
gz             -0.44227    1.1e-05        YES          1.0e-04           3             0.024

# Performance Summary
# Total Calculation Time: 0.082 seconds
# Memory Usage: 45.3 MB
# Solver Type: sparse (ARPACK)
# Convergence Test: PASSED

# Validation Results
# Analytical Solution (GaAs CB): -0.44
# Numerical Error: 0.05%
# Status: PASSED
```

##### 3. New convergence_test.dat

```
# Convergence Test Results
# Target State: 1
# Component: gx
# Test Date: 2025-11-02 15:30:42

# Step_Size [T]    g-factor    Uncertainty    Residual_Norm    Solver_Time [s]
1.0e-02            -0.42312    2.3e-03        1.2e-08          0.015
5.0e-03            -0.43758    1.1e-03        6.1e-09          0.014
1.0e-03            -0.44182    2.4e-04        2.8e-09          0.016
5.0e-04            -0.44201    1.2e-04        1.9e-09          0.015
1.0e-04            -0.44213    1.2e-05        1.1e-09          0.018
5.0e-05            -0.44214    6.1e-06        1.2e-09          0.019

# Convergence Analysis
# Optimal Step Size: 1.0e-04 T
# Converged: YES
# Convergence Error: 1.2e-05
# Richardson Extrapolation: -0.44214 ± 5.2e-06
# Extrapolation Available: YES
```

### 3. Error Handling Contract

#### Error Codes

```fortran
integer, parameter :: ERROR_NONE = 0
integer, parameter :: ERROR_INVALID_CONFIG = 1
integer, parameter :: ERROR_INVALID_TARGET_STATE = 2
integer, parameter :: ERROR_INVALID_FIELD = 3
integer, parameter :: ERROR_NUMERICAL_DIVERGENCE = 4
integer, parameter :: ERROR_EIGENVALUE_FAILURE = 5
integer, parameter :: ERROR_NON_CONVERGENCE = 6
integer, parameter :: ERROR_MEMORY_ALLOCATION = 7
integer, parameter :: ERROR_VALIDATION_FAILURE = 8
```

#### Error Messages

```fortran
character(len=256), parameter :: ERROR_MESSAGES(8) = [ &
  "Success", &
  "Invalid configuration parameters", &
  "Invalid target state index", &
  "Invalid magnetic field parameters", &
  "Numerical divergence detected", &
  "Eigenvalue solver failure", &
  "Non-convergence after maximum iterations", &
  "Memory allocation failure", &
  "Validation test failure" &
]
```

### 4. Integration Contract

#### Integration with Existing Modules

```fortran
! Integration with hamiltonianConstructor module
subroutine build_perturbed_hamiltonian(base_params, magnetic_field, direction, delta_B, H_perturbed)
  type(paramStruct), intent(in) :: base_params
  real(kind=dp), intent(in) :: magnetic_field(3)
  character(len=1), intent(in) :: direction
  real(kind=dp), intent(in) :: delta_B
  complex(kind=dp), intent(out) :: H_perturbed(:,:)
end subroutine

! Integration with outputFunctions module
subroutine write_gfactor_results(output_dir, g_tensor, config)
  character(len=*), intent(in) :: output_dir
  type(GFactorTensor), intent(in) :: g_tensor
  type(PerturbationConfig), intent(in) :: config
end subroutine

! Integration with utils module
subroutine solve_perturbed_eigenproblem(H_matrix, num_eigenvalues, eigenvalues, eigenvectors, use_sparse)
  complex(kind=dp), intent(in) :: H_matrix(:,:)
  integer, intent(in) :: num_eigenvalues
  complex(kind=dp), intent(out) :: eigenvalues(:)
  complex(kind=dp), intent(out) :: eigenvectors(:,:)
  logical, intent(in) :: use_sparse
end subroutine
```

### 5. Performance Contract

#### Performance Requirements

```fortran
! Maximum calculation time for standard quantum well
real(kind=dp), parameter :: MAX_CALCULATION_TIME = 30.0  ! seconds

! Maximum memory usage for typical problems
real(kind=dp), parameter :: MAX_MEMORY_USAGE = 1.0e9      ! bytes (1 GB)

! Minimum convergence tolerance
real(kind=dp), parameter :: MIN_TOLERANCE = 1e-15_dp

! Maximum problem size for dense solver
integer, parameter :: MAX_DENSE_SIZE = 1000              ! matrix dimension
```

#### Performance Metrics

```fortran
type :: PerformanceMetrics
  real(kind=dp) :: total_time = 0.0_dp
  real(kind=dp) :: solver_time = 0.0_dp
  real(kind=dp) :: setup_time = 0.0_dp
  real(kind=dp) :: memory_usage = 0.0_dp
  integer :: eigenvalue_solves = 0
  integer :: convergence_iterations = 0
  logical :: performance_targets_met = .false.
contains
  procedure :: validate_performance => metrics_validate
end type
```

### 6. Validation Contract

#### Validation Test Cases

```fortran
function get_validation_test_cases() result(test_cases)
  type(ValidationTestCase), allocatable :: test_cases(:)

  ! Allocate standard validation cases
  allocate(test_cases(6))

  ! Bulk GaAs conduction band
  test_cases(1) = ValidationTestCase( &
    "GaAs_bulk_conduction", &
    "GaAs", &
    -0.44_dp, &           ! analytical g-factor
    0.01_dp, &            ! tolerance
    1, &                  ! target state
    'z' &                 ! direction
  )

  ! Additional test cases...
end function
```

#### Validation Criteria

```fortran
function validate_calculation_results(g_tensor, config) result(validation_passed)
  type(GFactorTensor), intent(in) :: g_tensor
  type(PerturbationConfig), intent(in) :: config
  logical :: validation_passed

  ! Check convergence
  if (.not. all(g_tensor%converged)) then
    validation_passed = .false.
    return
  end if

  ! Check physical reasonableness
  if (any(abs(g_tensor%components) > 50.0_dp)) then
    validation_passed = .false.
    return
  end if

  ! Check uncertainty bounds
  if (any(g_tensor%uncertainties < 0.0_dp)) then
    validation_passed = .false.
    return
  end if

  validation_passed = .true.
end function
```

## Implementation Requirements

### 1. Thread Safety
All public procedures must be thread-safe for parallel execution of different g-factor calculations.

### 2. Memory Management
All allocated arrays must be properly deallocated. No memory leaks allowed.

### 3. Error Propagation
All procedures must propagate error codes up the call chain with appropriate error messages.

### 4. Backward Compatibility
Existing functionality must remain unchanged. New features are additive only.

### 5. Documentation
All public interfaces must have comprehensive documentation including:
- Purpose and usage
- Parameter descriptions and units
- Error conditions and recovery
- Performance characteristics
- Example usage

## Testing Requirements

### 1. Unit Tests
Test each public procedure with various input combinations including edge cases.

### 2. Integration Tests
Test integration with existing modules (hamiltonianConstructor, outputFunctions, utils).

### 3. Validation Tests
Validate against analytical solutions for standard test cases.

### 4. Performance Tests
Verify performance targets are met for standard problem sizes.

### 5. Regression Tests
Ensure new functionality doesn't break existing calculations.

This API contract provides a comprehensive specification for implementing the generic numerical perturbation method while maintaining compatibility with the existing k·p solver architecture.