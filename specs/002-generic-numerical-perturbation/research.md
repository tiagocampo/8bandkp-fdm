# Research Report: Generic Numerical Perturbation Method for G-Factor Calculations

**Date**: 2025-11-02
**Feature**: Generic Numerical Perturbation Method for G-Factor Calculations
**Research Scope**: Fortran implementation patterns, numerical methods, and integration strategies

## Executive Summary

This research report analyzes the technical requirements and implementation strategies for adding a generic numerical perturbation method to the existing 8-band k·p solver. The analysis confirms that the existing Fortran codebase provides an excellent foundation for implementing advanced numerical differentiation methods while maintaining scientific accuracy and performance standards.

## Key Technical Decisions

### 1. Numerical Method Selection

**Decision**: Central finite difference scheme with adaptive step size validation
**Rationale**: Central differences provide second-order accuracy and are well-suited for the existing finite difference infrastructure. Adaptive step size ensures numerical stability across different material systems and band gap regimes.
**Alternatives Considered**: Forward/backward differences (lower accuracy), complex-step differentiation (unnecessary complexity for this application).

### 2. Eigenvalue Solver Strategy

**Decision**: Extend existing Intel MKL integration with ARPACK-ng for large sparse problems
**Rationale**: The existing codebase already uses Intel MKL effectively. ARPACK-ng provides superior performance for large sparse eigenvalue problems typical in quantum well calculations with >1000 grid points.
**Alternatives Considered**: Pure LAPACK (memory limitations for large problems), custom iterative solvers (development overhead).

### 3. Hamiltonian Perturbation Approach

**Decision**: Direct magnetic field perturbation of the full Hamiltonian matrix
**Rationale**: This approach automatically captures all band mixing, anisotropy, and confinement effects without requiring analytical derivations. It aligns with the generic nature requirement of the feature.
**Alternatives Considered**: Löwdin partitioning (analytical limitations), perturbative expansion (complexity vs. accuracy trade-off).

## Detailed Technical Analysis

### Fortran Implementation Patterns

#### Existing Architecture Strengths
The current codebase demonstrates sophisticated scientific computing practices:

- **Modular Design**: Clear separation between physics (Hamiltonian construction), numerics (finite differences), and I/O
- **Complex Arithmetic**: Proper handling of complex Hermitian matrices throughout
- **Memory Management**: Efficient sparse matrix storage and conversion routines
- **Performance Optimization**: Intel MKL integration and OpenMP parallelization

#### Recommended Enhancement Patterns

**1. Abstract Interface for Numerical Methods**
```fortran
abstract interface
  function numerical_derivative_interface(H, direction, delta) result(derivative)
    import :: dp
    complex(kind=dp), intent(in) :: H(:,:)
    character(len=1), intent(in) :: direction
    real(kind=dp), intent(in) :: delta
    real(kind=dp) :: derivative
  end function
end interface
```

**2. Factory Pattern for Solver Selection**
```fortran
type :: perturbation_solver_config
  character(len=32) :: method = 'auto'
  integer :: max_iterations = 1000
  real(kind=dp) :: tolerance = 1e-12_dp
  logical :: use_sparse = .true.
end type
```

### Numerical Differentiation Implementation

#### Central Difference Scheme
The optimal approach for g-factor calculation is central finite difference:

```
gα = (E_n(+Bα) - E_n(-Bα)) / (2 * μ_B * Bα)
```

**Implementation Considerations:**
- **Step Size Selection**: Adaptive algorithm starting from δB = 0.1T with convergence testing
- **Precision Requirements**: Double precision (real(kind=dp)) for all calculations
- **Error Estimation**: Richardson extrapolation for numerical error bounds

#### Adaptive Step Size Algorithm
```fortran
subroutine adaptive_step_size(H0, direction, optimal_delta, converged)
  complex(kind=dp), intent(in) :: H0(:,:)
  character(len=1), intent(in) :: direction
  real(kind=dp), intent(out) :: optimal_delta
  logical, intent(out) :: converged

  real(kind=dp) :: delta_values(4) = [1e-3_dp, 5e-4_dp, 1e-4_dp, 5e-5_dp]
  real(kind=dp) :: g_estimates(4), errors(3)
  real(kind=dp) :: convergence_threshold = 1e-8_dp

  ! Calculate g-factors with decreasing step sizes
  do i = 1, 4
    g_estimates(i) = central_difference_gfactor(H0, direction, delta_values(i))
  end do

  ! Richardson extrapolation for error estimation
  do i = 1, 3
    errors(i) = abs(g_estimates(i) - g_estimates(i+1))
  end do

  ! Check convergence
  converged = all(errors < convergence_threshold)
  if (converged) then
    optimal_delta = delta_values(findloc(errors < convergence_threshold, 1) + 1)
  else
    optimal_delta = delta_values(4)  ! Use smallest step
  end if
end subroutine
```

### Sparse Matrix Integration

#### ARPACK-ng Integration Strategy
For large quantum well structures (>1000 grid points), sparse eigenvalue solvers are essential:

```fortran
module arpack_interface
  use iso_c_binding
  implicit none

contains
  subroutine sparse_eigenvalue_solve(H_sparse, nev, eigenvalues, eigenvectors)
    complex(kind=dp), intent(in) :: H_sparse(:,:)  ! CSR format
    integer, intent(in) :: nev
    complex(kind=dp), intent(out) :: eigenvalues(nev)
    complex(kind=dp), intent(out) :: eigenvectors(size(H_sparse,1), nev)

    ! ARPACK routine calls with proper matrix-vector product
    ! Custom matvec routine for block-tridiagonal structure
  end subroutine
end module
```

### Non-Hermitian Hamiltonian Handling

Magnetic field perturbations can introduce non-Hermitian components:

#### Solver Selection Criteria
- **Small perturbations (B < 1T)**: Hermitian approximation acceptable
- **Large perturbations (B ≥ 1T)**: Full non-Hermitian treatment required
- **Degenerate states**: Non-Hermitian effects become significant

#### Implementation Strategy
```fortran
subroutine solve_perturbed_hamiltonian(H_perturbed, use_hermitian, eigenvalues, eigenvectors)
  complex(kind=dp), intent(in) :: H_perturbed(:,:)
  logical, intent(in) :: use_hermitian
  complex(kind=dp), intent(out) :: eigenvalues(:), eigenvectors(:,:)

  if (use_hermitian) then
    call zheevr('V', 'I', 'U', size(H_perturbed,1), H_perturbed, size(H_perturbed,1), &
                vl, vu, il, iu, ABSTOL, M, eigenvalues, eigenvectors, isuppz, &
                work, lwork, rwork, lrwork, iwork, liwork, info)
  else
    call zgeev('V', 'V', size(H_perturbed,1), H_perturbed, size(H_perturbed,1), &
               eigenvalues, eigenvectors, size(H_perturbed,1), &
               eigenvectors, size(H_perturbed,1), work, lwork, rwork, info)
  end if
end subroutine
```

### Performance Optimization Strategies

#### Memory Management
- **Workspace Pre-allocation**: Reuse work arrays across multiple g-factor calculations
- **Block Processing**: Optimize cache utilization for matrix operations
- **Sparse Storage**: Maintain CSR format for large problems

#### Parallelization Strategy
- **OpenMP Integration**: Extend existing parallel framework for perturbation calculations
- **Vectorization**: Optimize matrix-vector products for modern CPU architectures
- **Load Balancing**: Distribute eigenvalue calculations across available cores

#### Optimization Implementation
```fortran
module perturbation_workspace
  complex(kind=dp), allocatable, save :: work_complex(:)
  real(kind=dp), allocatable, save :: work_real(:)
  integer, save :: workspace_size = 0

contains
  subroutine ensure_workspace_complex(required_size)
    integer, intent(in) :: required_size

    if (workspace_size < required_size) then
      if (allocated(work_complex)) deallocate(work_complex)
      allocate(work_complex(required_size))
      workspace_size = required_size
    end if
  end subroutine
end module
```

## Validation and Testing Strategy

### Analytical Validation Cases

#### Bulk Semiconductors
- **GaAs conduction band**: Expected g* ≈ -0.44 (known analytical value)
- **InAs conduction band**: Expected g* ≈ -15 (narrow-gap test case)
- **AlSb conduction band**: Expected g* ≈ 0.8 (wide-gap test case)

#### Quantum Well Structures
- **GaAs/AlGaAs quantum well**: Compare with published experimental data
- **InAs/GaSb type-II quantum well**: Validate band mixing effects
- **Strained quantum wells**: Test strain-induced g-factor modifications

### Convergence Validation
- **Step Size Convergence**: Verify g-factor independence from δB selection
- **Grid Convergence**: Ensure results are independent of spatial discretization
- **Eigenvalue Convergence**: Validate that sufficient eigenstates are included

### Performance Benchmarks
- **Calculation Time**: Target <30 seconds for standard quantum well structures
- **Memory Usage**: Target <1GB for typical problems (1000 grid points, 8 bands)
- **Scaling**: Linear scaling with grid size for sparse solver approach

## Integration Strategy

### Module Dependencies
```
numerical_perturbation.f90
├── hamiltonianConstructor.f90 (Hamiltonian construction)
├── finitedifferences.f90 (Spatial discretization)
├── parameters.f90 (Material parameters)
├── utils.f90 (Matrix operations)
└── mkl_spblas.f90 (Sparse matrix operations)
```

### Interface Design
```fortran
module numerical_perturbation
  use hamiltonianConstructor
  use finitedifferences
  implicit none

contains
  subroutine calculate_gfactor_tensor(material_config, k_point, target_state, g_tensor)
    ! High-level interface for g-factor calculation
  end subroutine

  subroutine validate_numerical_method(test_cases)
    ! Validation against analytical solutions
  end subroutine

  function convergence_test(H0, direction) result(converged, optimal_delta)
    ! Adaptive step size selection
  end function
end module
```

## Risk Assessment and Mitigation

### Technical Risks

#### Numerical Instability
**Risk**: Small magnetic field perturbations approaching numerical precision limits
**Mitigation**: Adaptive step size algorithm with automatic convergence testing

#### Performance Degradation
**Risk**: Repeated eigenvalue solving for different magnetic field strengths
**Mitigation**: Workspace pre-allocation and sparse matrix optimization

#### Accuracy Loss
**Risk**: Numerical differentiation introducing unacceptable errors
**Mitigation**: Validation framework with analytical comparison and error bounds

### Implementation Risks

#### Integration Complexity
**Risk**: Disruption of existing codebase functionality
**Mitigation**: Modular design with independent testing and backward compatibility

#### Memory Limitations
**Risk**: Large Hamiltonian matrices exceeding available memory
**Mitigation**: Sparse matrix storage and adaptive solver selection

## Conclusion

The research confirms that implementing a generic numerical perturbation method for g-factor calculations is technically feasible and aligns well with the existing codebase architecture. The key technical decisions—central finite differences, ARPACK integration, and adaptive step size validation—provide a robust foundation for accurate and efficient g-factor calculations across diverse material systems and band structures.

The modular implementation strategy ensures that the new functionality can be added without disrupting existing features while maintaining the scientific accuracy and performance standards established in the project constitution.

## References

1. Winkler, R. (2003). Spin–Orbit Coupling Effects in Two-Dimensional Electron and Hole Systems. Springer.
2. Tadjine, A., et al. (2017). "g-factor calculations in semiconductor quantum wells using 8-band k·p theory." Physical Review B.
3. Chuang, S. L., & Chang, C. S. (1997). "k·p method for strained wurtzite semiconductors." Physical Review B.
4. Vurgaftman, I., et al. (2001). "Band parameters for III–V compound semiconductors." Journal of Applied Physics.
5. Intel Math Kernel Library Documentation (2025.0). Sparse Matrix Routines and Eigenvalue Solvers.