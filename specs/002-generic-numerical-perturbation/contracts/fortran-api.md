# Fortran API Specification: Generic Numerical Perturbation

**Date**: 2025-11-02
**Version**: 1.0
**Module**: numerical_perturbation

## Public Interface

### Module: numerical_perturbation

#### High-Level API

##### Subroutine: calculate_gfactor_tensor

Calculates the complete g-tensor for a specified electronic state.

```fortran
subroutine calculate_gfactor_tensor(material_config, k_point, target_state, perturbation_config, g_tensor_result, info)
  ! Input parameters
  type(MaterialConfig), intent(in) :: material_config          ! Material system configuration
  real(kind=dp), intent(in) :: k_point(2)                       ! In-plane wavevector (kx, ky)
  integer, intent(in) :: target_state                           ! Target eigenstate index
  type(PerturbationConfig), intent(in) :: perturbation_config  ! Perturbation settings

  ! Output parameters
  type(GFactorTensor), intent(out) :: g_tensor_result          ! Complete g-tensor result
  integer, intent(out) :: info                                 ! Success/error status (0 = success)
end subroutine calculate_gfactor_tensor
```

**Parameters:**
- `material_config`: Material system definition including layers, compositions, and band parameters
- `k_point`: 2D in-plane wavevector in units of 2π/a (lattice constant)
- `target_state`: Index of target eigenstate (1-based indexing)
- `perturbation_config`: Numerical method configuration and tolerances
- `g_tensor_result`: Complete g-tensor with error estimates and convergence status
- `info`: Status code (0 = success, negative = error, positive = warning)

**Error Codes:**
- `0`: Success
- `-1`: Invalid material configuration
- `-2`: Target state out of range
- `-3`: Numerical convergence failure
- `-4`: Memory allocation error
- `-5`: Eigenvalue solver error

##### Subroutine: validate_perturbation_method

Validates the numerical perturbation implementation against analytical solutions.

```fortran
subroutine validate_perturbation_method(test_suite, validation_results, info)
  ! Input parameters
  character(len=*), intent(in) :: test_suite                     ! Test suite identifier

  ! Output parameters
  type(ValidationResult), allocatable, intent(out) :: validation_results(:)  ! Array of test results
  integer, intent(out) :: info                                 ! Overall status (0 = all passed)
end subroutine validate_perturbation_method
```

**Test Suites Available:**
- `"bulk_semiconductors"`: Bulk material validation (GaAs, InAs, AlSb)
- `"quantum_wells"`: Quantum well structure validation
- `"convergence"`: Step size and grid convergence testing
- `"performance"`: Performance benchmark validation

##### Function: check_convergence

Performs adaptive convergence testing for perturbation step size.

```fortran
function check_convergence(H0, direction, initial_delta, convergence_result) result(converged)
  ! Input parameters
  complex(kind=dp), intent(in) :: H0(:,:)                      ! Base Hamiltonian matrix
  character(len=1), intent(in) :: direction                   ! Perturbation direction ('x', 'y', 'z')
  real(kind=dp), intent(in) :: initial_delta                  ! Initial perturbation strength

  ! Output parameters
  type(ConvergenceResult), intent(out) :: convergence_result  ! Detailed convergence information

  ! Return value
  logical :: converged                                         ! Overall convergence status
end function check_convergence
```

#### Utility Functions

##### Subroutine: initialize_perturbation_workspace

Initializes workspace for repeated perturbation calculations.

```fortran
subroutine initialize_perturbation_workspace(matrix_size, max_states, info)
  ! Input parameters
  integer, intent(in) :: matrix_size                           ! Hamiltonian matrix dimension
  integer, intent(in) :: max_states                           ! Maximum number of states to track

  ! Output parameters
  integer, intent(out) :: info                                 ! Initialization status (0 = success)
end subroutine initialize_perturbation_workspace
```

##### Subroutine: cleanup_perturbation_workspace

Cleans up allocated workspace memory.

```fortran
subroutine cleanup_perturbation_workspace()
  ! No parameters - cleans up all allocated memory
end subroutine cleanup_perturbation_workspace
```

##### Function: estimate_calculation_time

Estimates calculation time based on problem size and configuration.

```fortran
function estimate_calculation_time(matrix_size, num_states, use_sparse) result(estimated_time)
  ! Input parameters
  integer, intent(in) :: matrix_size                           ! Hamiltonian matrix dimension
  integer, intent(in) :: num_states                            ! Number of states to calculate
  logical, intent(in) :: use_sparse                            ! Use sparse solver

  ! Return value
  real(kind=dp) :: estimated_time                               ! Estimated time in seconds
end function estimate_calculation_time
```

### Module: perturbation_solver

#### Low-Level Solver Interface

##### Subroutine: solve_hamiltonian_perturbation

Solves perturbed Hamiltonian eigenvalue problems.

```fortran
subroutine solve_hamiltonian_perturbation(H_base, perturbation_matrix, delta, num_eigenvalues, eigenvalues, eigenvectors, info)
  ! Input parameters
  complex(kind=dp), intent(in) :: H_base(:,:)                  ! Base Hamiltonian (B=0)
  complex(kind=dp), intent(in) :: perturbation_matrix(:,:)    ! Perturbation matrix
  real(kind=dp), intent(in) :: delta                          ! Perturbation strength
  integer, intent(in) :: num_eigenvalues                      ! Number of eigenvalues to compute

  ! Output parameters
  complex(kind=dp), intent(out) :: eigenvalues(num_eigenvalues)        ! Perturbed eigenvalues
  complex(kind=dp), intent(out) :: eigenvectors(size(H_base,1), num_eigenvalues)  ! Eigenvectors
  integer, intent(out) :: info                                 ! Solver status (0 = success)
end subroutine solve_hamiltonian_perturbation
```

##### Subroutine: apply_magnetic_perturbation

Applies magnetic field perturbation to Hamiltonian.

```fortran
subroutine apply_magnetic_perturbation(H0, B_field, direction, H_perturbed, info)
  ! Input parameters
  complex(kind=dp), intent(in) :: H0(:,:)                      ! Base Hamiltonian
  real(kind=dp), intent(in) :: B_field                        ! Magnetic field strength (Tesla)
  character(len=1), intent(in) :: direction                   ! Field direction

  ! Output parameters
  complex(kind=dp), intent(out) :: H_perturbed(size(H0,1), size(H0,2))  ! Perturbed Hamiltonian
  integer, intent(out) :: info                                 ! Status (0 = success)
end subroutine apply_magnetic_perturbation
```

### Module: validation_framework

#### Validation Interface

##### Subroutine: run_validation_test

Runs a single validation test case.

```fortran
subroutine run_validation_test(test_case, result, info)
  ! Input parameters
  type(ValidationTestCase), intent(in) :: test_case           ! Test case definition

  ! Output parameters
  type(ValidationResult), intent(out) :: result              ! Test result
  integer, intent(out) :: info                                 ! Status (0 = success)
end subroutine run_validation_test
```

##### Function: compare_with_analytical

Compares numerical result with analytical solution.

```fortran
function compare_with_analytical(numerical_value, analytical_value, tolerance) result(agree)
  ! Input parameters
  real(kind=dp), intent(in) :: numerical_value               ! Numerically computed value
  real(kind=dp), intent(in) :: analytical_value              ! Analytical solution
  real(kind=dp), intent(in) :: tolerance                     ! Acceptable tolerance

  ! Return value
  logical :: agree                                            ! Agreement status
end function compare_with_analytical
```

## Data Types

### Derived Types (Public)

```fortran
type :: PerturbationConfig
  character(len=32) :: method = 'central_difference'
  real(kind=dp) :: magnetic_field_strength = 0.1_dp
  real(kind=dp) :: convergence_tolerance = 1e-12_dp
  integer :: max_iterations = 10
  logical :: use_sparse_solver = .true.
  logical :: validate_analytical = .true.
  character(len=1) :: perturbation_direction = 'z'
end type

type :: GFactorTensor
  real(kind=dp) :: gx = 0.0_dp
  real(kind=dp) :: gy = 0.0_dp
  real(kind=dp) :: gz = 0.0_dp
  real(kind=dp) :: error_estimate(3) = 0.0_dp
  logical :: converged(3) = .false.
  integer :: target_state = 1
  real(kind=dp) :: target_energy = 0.0_dp
  character(len=32) :: calculation_method = ''
end type

type :: ConvergenceResult
  real(kind=dp) :: optimal_delta = 0.0_dp
  logical :: converged = .false.
  real(kind=dp) :: error_estimate = 0.0_dp
  real(kind=dp) :: delta_values(4) = 0.0_dp
  real(kind=dp) :: g_estimates(4) = 0.0_dp
  integer :: convergence_iterations = 0
  character(len=128) :: convergence_message = ''
end type

type :: ValidationTestCase
  character(len=64) :: test_name = ''
  character(len=32) :: material_system = ''
  real(kind=dp) :: analytical_g_factor = 0.0_dp
  real(kind=dp) :: tolerance = 1e-6_dp
  integer :: target_state = 1
  character(len=1) :: direction = 'z'
end type

type :: ValidationResult
  type(ValidationTestCase) :: test_case
  logical :: test_passed = .false.
  real(kind=dp) :: numerical_g_factor = 0.0_dp
  real(kind=dp) :: absolute_error = 0.0_dp
  real(kind=dp) :: relative_error = 0.0_dp
  character(len=256) :: error_message = ''
end type
```

## Constants and Parameters

### Physical Constants

```fortran
real(kind=dp), parameter :: MU_BOHR = 9.274009994e-24_dp     ! Bohr magneton (J/T)
real(kind=dp), parameter :: EV_TO_JOULE = 1.602176634e-19_dp ! eV to Joules conversion
real(kind=dp), parameter :: HBAR = 1.054571817e-34_dp       ! Reduced Planck constant (J·s)
```

### Default Configuration Values

```fortran
real(kind=dp), parameter :: DEFAULT_DELTA_B = 0.1_dp         ! Default magnetic field (Tesla)
real(kind=dp), parameter :: DEFAULT_TOLERANCE = 1e-12_dp     ! Default convergence tolerance
integer, parameter :: DEFAULT_MAX_ITERATIONS = 10            ! Default max refinement iterations
real(kind=dp), parameter :: MIN_DELTA_B = 1e-6_dp            ! Minimum perturbation strength
real(kind=dp), parameter :: MAX_DELTA_B = 1e-1_dp            ! Maximum perturbation strength
```

## Error Handling

### Status Codes

| Code | Meaning | Action |
|------|---------|--------|
| 0 | Success | Continue processing |
| >0 | Warning | Log warning, continue |
| <0 | Error | Abort operation, return error |

### Common Error Messages

```fortran
character(len=*), parameter :: ERR_INVALID_CONFIG = "Invalid perturbation configuration"
character(len=*), parameter :: ERR_CONVERGENCE_FAILED = "Numerical convergence failed"
character(len=*), parameter :: ERR_MEMORY_ALLOCATION = "Memory allocation failed"
character(len=*), parameter :: ERR_SOLVER_FAILED = "Eigenvalue solver failed"
character(len=*), parameter :: ERR_INVALID_STATE = "Invalid target state index"
```

## Usage Examples

### Basic G-Factor Calculation

```fortran
program example_gfactor_calc
  use numerical_perturbation
  use hamiltonianConstructor
  implicit none

  type(MaterialConfig) :: mat_config
  type(PerturbationConfig) :: pert_config
  type(GFactorTensor) :: g_result
  real(kind=dp) :: k_point(2)
  integer :: info

  ! Initialize material configuration
  mat_config = create_gaas_quantum_well()

  ! Set up perturbation configuration
  pert_config = PerturbationConfig( &
    method = 'central_difference', &
    magnetic_field_strength = 0.1_dp, &
    convergence_tolerance = 1e-12_dp &
  )

  ! Calculate g-factor for first conduction subband
  k_point = [0.0_dp, 0.0_dp]  ! Γ point
  call calculate_gfactor_tensor(mat_config, k_point, 1, pert_config, g_result, info)

  if (info == 0) then
    print *, "G-factor tensor:", g_result%gx, g_result%gy, g_result%gz
    print *, "Convergence:", g_result%converged
    print *, "Error estimates:", g_result%error_estimate
  else
    print *, "Calculation failed with error:", info
  end if
end program example_gfactor_calc
```

### Validation Testing

```fortran
program example_validation
  use numerical_perturbation
  implicit none

  type(ValidationResult), allocatable :: results(:)
  integer :: info, i

  ! Run bulk semiconductor validation
  call validate_perturbation_method('bulk_semiconductors', results, info)

  if (info == 0) then
    print *, "All validation tests passed!"
    do i = 1, size(results)
      print *, "Test:", trim(results(i)%test_case%test_name), &
               "Passed:", results(i)%test_passed, &
               "Error:", results(i)%relative_error
    end do
  else
    print *, "Validation failed with", info, "test(s) failing"
  end if
end program example_validation
```

## Integration Notes

### Dependencies
- `hamiltonianConstructor`: For Hamiltonian matrix construction
- `finitedifferences`: For spatial discretization
- `parameters`: For material parameter database
- `utils`: For matrix operations and utilities
- Intel MKL: For eigenvalue solving and sparse matrix operations

### Thread Safety
- All public subroutines are thread-safe when called with different data
- Workspace initialization should be called once per thread
- Global state is read-only after initialization

### Memory Management
- Users should call `initialize_perturbation_workspace` before calculations
- `cleanup_perturbation_workspace` should be called when finished
- Memory usage scales as O(N²) for dense matrices, O(N) for sparse

### Performance Considerations
- Use sparse solvers for matrices larger than 500×500
- Cache workspace for repeated calculations with similar matrix sizes
- Parallel execution available for multiple states or k-points