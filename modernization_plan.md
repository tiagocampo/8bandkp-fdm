# Code Modernization Plan for 8bandkp-fdm

## 1. Current Code Analysis

### Current Structure
- `main.f90`: Main program for band structure calculations
- `main_gfactor.f90`: Main program for g-factor calculations
- `hamiltonianConstructor.f90`: Core physics implementation (681 lines)
- `parameters.f90`: Material parameters database (529 lines)
- `defs.f90`: Basic type definitions
- `finitedifferences.f90`: FDM implementation
- `outputFunctions.f90`: I/O operations
- `gfactor_functions.f90`: G-factor calculations
- `mkl_sparse_handle.f90` & `mkl_spblas.f90`: MKL interfaces
- `utils.f90`: Utility functions

### Current Issues
1. Large monolithic files (`hamiltonianConstructor.f90`, `parameters.f90`)
2. Mixed concerns within files (physics, numerics, I/O)
3. Limited use of modern Fortran features
4. Basic error handling
5. Tight coupling between components

## 2. Modernization Strategy

### Phase 1: Core Infrastructure (2-3 weeks)

#### 1.1 Type System Modernization
```fortran
! Core/Types/HamiltonianBase.f90
type, abstract :: HamiltonianBase
    private
    integer :: dimension
    complex(dp), allocatable :: matrix(:,:)
contains
    procedure(construct_interface), deferred :: construct
    procedure(diagonalize_interface), deferred :: diagonalize
end type

! Core/Types/MaterialBase.f90
type, abstract :: MaterialBase
    private
    character(len=:), allocatable :: name
    real(dp) :: bandgap
    ! ... other common properties
contains
    procedure(get_parameters_interface), deferred :: get_parameters
end type
```

#### 1.2 Error Handling System
```fortran
! Core/ErrorHandling.f90
type :: ErrorContext
    logical :: has_error = .false.
    character(len=:), allocatable :: message
    integer :: error_code = 0
contains
    procedure :: raise_error
    procedure :: clear_error
end type
```

### Phase 2: Physics Implementation (3-4 weeks)

#### 2.1 Hamiltonian Module Refactoring
1. Split `hamiltonianConstructor.f90` into:
   - `Physics/Hamiltonian/ZincBlende8Band.f90`
   - `Physics/Hamiltonian/BulkHamiltonian.f90`
   - `Physics/Hamiltonian/QuantumWellHamiltonian.f90`

2. Create material system modules:
   - `Physics/Materials/ZincBlendeMaterial.f90`
   - `Physics/Materials/ParameterDatabase.f90`
   - `Physics/Materials/Heterostructure.f90`

#### 2.2 Numerical Methods Refactoring
1. Enhance finite difference module:
   - `Numerics/FiniteDifference/Operators.f90`
   - `Numerics/FiniteDifference/Grid.f90`

2. Linear algebra interfaces:
   - `Numerics/LinearAlgebra/Eigensolvers.f90`
   - `Numerics/LinearAlgebra/SparseOperations.f90`

### Phase 3: I/O and Utilities (2-3 weeks)

#### 3.1 Input/Output Modernization
1. Create structured input handling:
```fortran
! IO/Input.f90
type :: InputParameters
    type(CalculationType) :: calc_type
    type(MaterialSystem) :: materials
    type(NumericalParams) :: num_params
contains
    procedure :: read_from_file
    procedure :: validate
end type
```

2. Modernize output system:
```fortran
! IO/Output.f90
type :: ResultsWriter
    character(len=:), allocatable :: output_dir
contains
    procedure :: write_eigenvalues
    procedure :: write_eigenfunctions
    procedure :: write_band_structure
end type
```

### Phase 4: Integration and Testing (2-3 weeks)

#### 4.1 Main Program Refactoring
1. Create calculation factory:
```fortran
! Core/Factory.f90
type :: CalculationFactory
contains
    procedure :: create_bulk_calculation
    procedure :: create_qw_calculation
    procedure :: create_gfactor_calculation
end type
```

2. Implement calculation types:
```fortran
! Calculations/BulkCalculation.f90
type, extends(CalculationBase) :: BulkCalculation
contains
    procedure :: initialize
    procedure :: run
    procedure :: finalize
end type
```

#### 4.2 Testing Framework
1. Unit tests for each module
2. Integration tests for full calculations
3. Performance benchmarks

## 3. Implementation Plan

### Week 1-2: Core Infrastructure and Build System
- [ ] Set up new directory structure
- [ ] Implement CMake build system
- [ ] Create installation documentation
- [ ] Implement base types and interfaces
- [ ] Create error handling system

### Week 3-4: Physics Modules
- [ ] Refactor Hamiltonian construction
- [ ] Implement material system
- [ ] Create numerical methods modules

### Week 5-6: I/O and Integration
- [ ] Modernize input/output system
- [ ] Implement calculation factory
- [ ] Create main program structure

### Week 7-8: Testing and Documentation
- [ ] Write unit tests
- [ ] Create integration tests
- [ ] Update documentation

## 4. Migration Strategy

1. **Incremental Approach**
   - Keep existing code functional while developing new modules
   - Create parallel implementations
   - Gradually switch to new implementations

2. **Testing Protocol**
   - Compare results between old and new implementations
   - Maintain test suite for regression testing
   - Document all changes and their impacts

3. **Documentation**
   - Update README with new architecture
   - Create module-level documentation
   - Provide migration guides for users

## 5. Success Criteria

1. **Code Quality**
   - Reduced file sizes (no file > 300 lines)
   - Increased modularity
   - Improved error handling
   - Better type safety

2. **Performance**
   - Equal or better performance
   - Reduced memory usage
   - Better parallel scaling

3. **Maintainability**
   - Clear module dependencies
   - Comprehensive tests
   - Updated documentation
   - Easier to extend 

## 6. Build System Modernization

### 6.1 CMake Migration (1-2 weeks)

#### Directory Structure
```
8bandkp-fdm/
├── CMakeLists.txt
├── cmake/
│   ├── FindMKL.cmake
│   ├── FindLAPACK.cmake
│   └── Fortran.cmake
├── src/
│   ├── CMakeLists.txt
│   ├── Core/
│   │   ├── Types/
│   │   │   ├── hamiltonian_base.f90
│   │   │   ├── material_base.f90
│   │   │   └── calculation_base.f90
│   │   ├── error_handling.f90
│   │   └── constants.f90
│   ├── Physics/
│   │   ├── Hamiltonian/
│   │   │   ├── zincblende_8band.f90
│   │   │   ├── bulk_hamiltonian.f90
│   │   │   └── quantum_well_hamiltonian.f90
│   │   ├── Materials/
│   │   │   ├── material_database.f90
│   │   │   ├── heterostructure.f90
│   │   │   └── strain.f90
│   │   └── BandStructure/
│   │       ├── band_solver.f90
│   │       └── gfactor.f90
│   ├── Numerics/
│   │   ├── FiniteDifference/
│   │   │   ├── fd_grid.f90
│   │   │   └── fd_operators.f90
│   │   └── LinearAlgebra/
│   │       ├── eigen_solver.f90
│   │       └── sparse_operations.f90
│   ├── IO/
│   │   ├── input_handler.f90
│   │   ├── output_handler.f90
│   │   └── visualization.f90
│   └── Applications/
│       ├── band_structure_app.f90
│       └── gfactor_app.f90
├── tests/
│   ├── CMakeLists.txt
│   ├── unit/
│   │   ├── test_hamiltonian.f90
│   │   ├── test_materials.f90
│   │   └── test_solvers.f90
│   └── integration/
│       ├── test_bulk.f90
│       └── test_quantum_well.f90
├── examples/
│   ├── CMakeLists.txt
│   ├── bulk/
│   └── quantum_well/
└── docs/
    ├── api/
    ├── examples/
    └── theory/
```

## 9. Code Organization and Style Guide

### 9.1 Naming Conventions

#### Module Names
- Use `snake_case` for module names
- Descriptive and specific: `quantum_well_hamiltonian` not just `hamiltonian`
- Format: `<component>_<subcomponent>`

```fortran
module quantum_well_hamiltonian
module material_database
module fd_operators
```

#### Type Names
- Use `PascalCase` for derived types
- Include the word 'Type' for clarity
- Format: `<Component>Type`

```fortran
type :: HamiltonianType
type :: MaterialParametersType
type :: CalculationType
```

#### Procedure Names
- Use `snake_case` for procedures
- Verb-first for actions
- Format: `<action>_<target>`

```fortran
subroutine calculate_band_structure
function get_material_parameters
subroutine solve_eigenvalue_problem
```

#### Variable Names
- Use `snake_case` for variables
- Clear and descriptive names
- Common abbreviations allowed if clear

```fortran
real(dp) :: band_gap
integer :: num_bands
complex(dp) :: wave_function
```

### 9.2 Module Organization

#### Template for Module Structure
```fortran
module module_name
    use iso_fortran_env, only: dp => real64
    use error_handling, only: ErrorType
    
    implicit none
    private  ! Default private visibility
    
    ! Public API
    public :: PublicType
    public :: public_procedure
    
    ! Types
    type :: PublicType
        private
        ! Components
    contains
        ! Type-bound procedures
    end type
    
    ! Module parameters
    real(dp), parameter :: CONSTANT = 1.0_dp
    
    ! Module variables (if needed)
    
    interface
        ! Abstract interfaces
    end interface
    
contains
    ! Procedures implementation
end module
```

### 9.3 Documentation Standards

#### Module Documentation
```fortran
!> @brief Brief description of the module
!> @details Detailed description of the module's purpose and functionality
module example_module
```

#### Procedure Documentation
```fortran
!> @brief Calculate the band structure
!> @param[in]  k_points    K-points array
!> @param[out] energies    Calculated energies
!> @return error_code      0 if successful
subroutine calculate_bands(k_points, energies, error_code)
```

### 9.4 Error Handling Pattern
```fortran
subroutine some_calculation(input, output, error)
    type(InputType), intent(in) :: input
    type(OutputType), intent(out) :: output
    type(ErrorType), intent(out) :: error
    
    ! Early return on error
    if (validate_input(input, error)) return
    
    ! Main calculation with error checking
    if (perform_calculation(input, output, error)) return
    
    ! Success path
    error%code = 0
end subroutine
``` 

## 10. Error Handling Framework

### 10.1 Error Type System
```fortran
! Core/Types/error_types.f90
module error_types
    implicit none
    private

    !> Error severity levels
    integer, parameter, public :: ERROR_NONE = 0
    integer, parameter, public :: ERROR_WARNING = 1
    integer, parameter, public :: ERROR_FATAL = 2

    !> Error categories
    integer, parameter, public :: ERR_PHYSICS = 100  ! Physics-related errors
    integer, parameter, public :: ERR_NUMERICS = 200 ! Numerical computation errors
    integer, parameter, public :: ERR_INPUT = 300    ! Input validation errors
    integer, parameter, public :: ERR_SYSTEM = 400   ! System/resource errors

    !> Main error type for error handling
    type, public :: ErrorType
        private
        integer :: code = ERROR_NONE
        integer :: category = 0
        integer :: severity = ERROR_NONE
        character(len=:), allocatable :: message
        character(len=:), allocatable :: procedure_name
        character(len=:), allocatable :: module_name
    contains
        procedure :: set_error
        procedure :: clear_error
        procedure :: has_error
        procedure :: print_error
        procedure :: is_fatal
    end type ErrorType

    !> Interface for error checking procedures
    abstract interface
        function error_check(self, error) result(has_error)
            import ErrorType
            class(*), intent(in) :: self
            type(ErrorType), intent(inout) :: error
            logical :: has_error
        end function
    end interface
end module

! Core/error_handling.f90
module error_handling
    use error_types
    implicit none
    private

    public :: handle_error
    public :: validate_parameters
    public :: check_allocation

contains
    !> Central error handling routine
    subroutine handle_error(error, stop_on_error)
        type(ErrorType), intent(in) :: error
        logical, intent(in), optional :: stop_on_error
        
        if (error%is_fatal() .or. (present(stop_on_error) .and. stop_on_error)) then
            call error%print_error()
            error stop
        else
            call error%print_error()
        end if
    end subroutine

    !> Parameter validation with dimensional analysis
    function validate_parameters(params, error) result(is_valid)
        type(MaterialParametersType), intent(in) :: params
        type(ErrorType), intent(inout) :: error
        logical :: is_valid

        ! Check physical constraints
        if (params%band_gap <= 0.0_dp) then
            call error%set_error(ERR_PHYSICS, ERROR_FATAL, &
                               "Band gap must be positive")
            is_valid = .false.
            return
        end if

        ! Check dimensional consistency
        if (params%effective_mass <= 0.0_dp) then
            call error%set_error(ERR_PHYSICS, ERROR_FATAL, &
                               "Effective mass must be positive")
            is_valid = .false.
            return
        end if

        is_valid = .true.
    end function

    !> Memory allocation checker
    subroutine check_allocation(array, error, array_name)
        class(*), allocatable, intent(in) :: array
        type(ErrorType), intent(inout) :: error
        character(len=*), intent(in) :: array_name
        
        if (.not. allocated(array)) then
            call error%set_error(ERR_SYSTEM, ERROR_FATAL, &
                               "Failed to allocate "//array_name)
        end if
    end subroutine
end module
```

### 10.2 Physics Validation
```fortran
! Physics/Validation/physics_checks.f90
module physics_checks
    use error_types
    implicit none
    private

    public :: validate_hamiltonian
    public :: validate_wavevector
    public :: validate_material_parameters

contains
    function validate_hamiltonian(H, error) result(is_valid)
        complex(dp), intent(in) :: H(:,:)
        type(ErrorType), intent(inout) :: error
        logical :: is_valid

        ! Check Hermiticity
        if (.not. is_hermitian(H)) then
            call error%set_error(ERR_PHYSICS, ERROR_FATAL, &
                               "Hamiltonian is not Hermitian")
            is_valid = .false.
            return
        end if

        is_valid = .true.
    end function

    function validate_wavevector(k, error) result(is_valid)
        real(dp), intent(in) :: k(3)
        type(ErrorType), intent(inout) :: error
        logical :: is_valid

        ! Check k-vector magnitude (within first Brillouin zone)
        if (norm2(k) > PI/a_0) then
            call error%set_error(ERR_PHYSICS, ERROR_WARNING, &
                               "k-vector outside first Brillouin zone")
            is_valid = .false.
            return
        end if

        is_valid = .true.
    end function
end module
```

### 10.3 Numerical Validation
```fortran
! Numerics/Validation/numerical_checks.f90
module numerical_checks
    use error_types
    implicit none
    private

    public :: validate_grid
    public :: check_convergence
    public :: validate_eigenvalues

contains
    function validate_grid(grid, error) result(is_valid)
        type(GridType), intent(in) :: grid
        type(ErrorType), intent(inout) :: error
        logical :: is_valid

        ! Check grid spacing
        if (grid%delta <= 0.0_dp) then
            call error%set_error(ERR_NUMERICS, ERROR_FATAL, &
                               "Invalid grid spacing")
            is_valid = .false.
            return
        end if

        ! Check grid points
        if (grid%num_points < 3) then
            call error%set_error(ERR_NUMERICS, ERROR_FATAL, &
                               "Insufficient grid points")
            is_valid = .false.
            return
        end if

        is_valid = .true.
    end function
end module
```

### 10.4 Error Usage Pattern
```fortran
! Example usage in a physics calculation
subroutine calculate_band_structure(k_points, materials, bands, error)
    real(dp), intent(in) :: k_points(:,:)
    type(MaterialType), intent(in) :: materials(:)
    type(BandStructureType), intent(out) :: bands
    type(ErrorType), intent(out) :: error

    ! Validate inputs
    if (.not. validate_materials(materials, error)) return
    if (.not. validate_k_points(k_points, error)) return

    ! Allocate resources with checks
    allocate(bands%energies(size(k_points,1)), stat=alloc_stat)
    call check_allocation(bands%energies, error, "bands%energies")
    if (error%has_error()) return

    ! Perform calculation with error checking
    do ik = 1, size(k_points,1)
        if (.not. calculate_bands_at_k(k_points(ik,:), materials, &
                                     bands%energies(ik), error)) return
    end do

    ! Success path
    call error%clear_error()
end subroutine
``` 

## 11. Documentation Framework

### 11.1 FORD Configuration
```fortran
! ford.md
project: 8bandkp-fdm
version: 2.0.0
author: Tiago de Campos
email: tiago@example.com
project_github: https://github.com/username/8bandkp-fdm
summary: 8-band k·p method implementation using finite differences
preprocessor: gfortran -E
predocmark: >
media_dir: ./docs/media
docmark_alt: #
predocmark_alt: <
display: public
         protected
source: true
graph: true
search: true
macro: TEST
       LOGIC=.true.
extra_mods: json_module: http://jacobwilliams.github.io/json-fortran/
            futility: http://cmacmackin.github.io
license: gfdl
```

### 11.2 Physics Documentation Template
```fortran
!> @brief 8-band k·p Hamiltonian for zinc-blende semiconductors
!>
!> @details This module implements the 8-band k·p Hamiltonian following:
!> 1. Kane's model for band coupling
!> 2. Strain effects using Bir-Pikus theory
!> 3. External field contributions
!>
!> The Hamiltonian matrix elements are given by:
!> \f[
!>   H_{cc} = E_c + \frac{\hbar^2}{2m_0}[(2F + 1)k^2]
!> \f]
!>
!> @note The basis states are ordered as:
!> 1. |c↑⟩: s-like conduction band (J=1/2, mj=+1/2)
!> 2. |c↓⟩: s-like conduction band (J=1/2, mj=-1/2)
!> 3. |hh↑⟩: heavy hole (J=3/2, mj=+3/2)
!> ...
!>
!> @references
!> 1. Chuang, S. L. Physics of Optoelectronic Devices (1995)
!> 2. Winkler, R. Spin-Orbit Coupling Effects in Two-Dimensional Electron 
!>    and Hole Systems (2003)
module hamiltonian_8band
```

### 11.3 Implementation Documentation Template
```fortran
!> @brief Calculate band structure using finite difference method
!>
!> @details This subroutine solves the k·p Hamiltonian using finite differences.
!> The procedure follows these steps:
!> 1. Set up the finite difference grid
!> 2. Construct the Hamiltonian matrix
!> 3. Apply boundary conditions
!> 4. Solve the eigenvalue problem
!>
!> @param[in]  k_points   K-points array (3, n_points)
!> @param[in]  materials  Material parameters for each layer
!> @param[out] energies   Calculated energy eigenvalues
!> @param[out] states     Calculated eigenstates
!> @param[out] error      Error handling object
!>
!> @note The finite difference grid must have sufficient points for convergence
!>
!> @warning Large k-vectors may lead to numerical instabilities
!>
!> @todo Implement adaptive grid refinement
!>
!> @bug Currently shows spurious solutions for very thin quantum wells
subroutine calculate_band_structure(k_points, materials, energies, states, error)
```

### 11.4 Theory Documentation

Create `docs/theory/`:

1. `quantum_theory.md`:
```markdown
# Quantum Mechanical Framework

## k·p Theory
The k·p method is based on perturbation theory applied to the one-electron
Schrödinger equation in a periodic potential:

\[
[-\frac{\hbar^2}{2m_0}\nabla^2 + V(\mathbf{r})]\psi_{n\mathbf{k}}(\mathbf{r}) 
= E_n(\mathbf{k})\psi_{n\mathbf{k}}(\mathbf{r})
\]

### Band Coupling
The 8-band model includes coupling between:
- Conduction band (Γ6)
- Heavy and light hole bands (Γ8)
- Split-off band (Γ7)
...
```

2. `numerical_methods.md`:
```markdown
# Numerical Implementation

## Finite Difference Method
The spatial derivatives in the k·p Hamiltonian are discretized using:

\[
\frac{d^2f}{dz^2} \approx \frac{f_{i+1} - 2f_i + f_{i-1}}{(\Delta z)^2}
\]

### Grid Requirements
- Minimum points per wavelength: 10
- Recommended grid spacing: < 5Å
...
```

### 11.5 Documentation Build Process

Add to CMakeLists.txt:
```cmake
# Documentation generation
find_program(FORD ford)
if(FORD)
    add_custom_target(docs
        COMMAND ${FORD} ${CMAKE_SOURCE_DIR}/ford.md
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating documentation with FORD"
        VERBATIM
    )
endif()

# Theory documentation (using pandoc)
find_program(PANDOC pandoc)
if(PANDOC)
    add_custom_target(theory_docs
        COMMAND ${PANDOC} -s ${CMAKE_SOURCE_DIR}/docs/theory/*.md
                -o ${CMAKE_BINARY_DIR}/docs/theory.pdf
                --toc --highlight-style=tango
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating theory documentation"
        VERBATIM
    )
endif()
```

### 11.6 Documentation Standards

1. **Module Documentation**:
   - Brief description
   - Detailed physics explanation
   - Mathematical formulas using LaTeX
   - References to literature
   - Usage examples

2. **Procedure Documentation**:
   - Purpose and algorithm description
   - Parameter descriptions with units
   - Error conditions
   - Performance considerations
   - Known limitations

3. **Physics Documentation**:
   - Theoretical background
   - Approximations used
   - Validity ranges
   - References to original papers

4. **Implementation Notes**:
   - Numerical methods details
   - Convergence criteria
   - Performance optimizations
   - Parallelization strategy

### 11.7 Documentation Maintenance

1. **Review Process**:
   - Documentation review with code changes
   - Physics validation by domain experts
   - Example verification

2. **Update Checklist**:
   - Theory documentation matches implementation
   - All new features documented
   - Examples updated
   - Build instructions current

## 12. Testing Framework

### 12.1 pFUnit Integration

#### CMake Configuration
```cmake
# tests/CMakeLists.txt
find_package(PFUNIT REQUIRED)
enable_testing()

# Function to create test suite
function(add_physics_test test_name)
    add_pfunit_test(${test_name}
        TEST_SOURCES ${test_name}.pf
        LINK_LIBRARIES physics_lib
        MAX_PES 4
    )
endfunction()

# Test suites
add_physics_test(test_hamiltonian)
add_physics_test(test_band_structure)
add_physics_test(test_finite_difference)
```

### 12.2 Unit Tests

#### 1. Hamiltonian Tests
```fortran
! tests/unit/test_hamiltonian.pf
@test
subroutine test_hamiltonian_hermiticity()
    use hamiltonian_8band, only: construct_hamiltonian
    use assert_utils, only: assert_is_hermitian
    
    type(HamiltonianType) :: H
    real(dp) :: k_point(3) = [0.0_dp, 0.0_dp, 0.0_dp]
    
    ! Test at Gamma point
    call construct_hamiltonian(H, k_point)
    @assertTrue(assert_is_hermitian(H%matrix))
    
    ! Test at finite k
    k_point = [0.1_dp, 0.0_dp, 0.0_dp]
    call construct_hamiltonian(H, k_point)
    @assertTrue(assert_is_hermitian(H%matrix))
end subroutine
```

### 12.3 Analytical Validation Tests

#### 1. Simple Quantum Well
```fortran
! tests/validation/test_simple_well.pf
@test
subroutine test_infinite_well()
    use band_solver, only: solve_schrodinger
    use assert_utils, only: assert_close
    
    real(dp) :: L = 10.0_dp  ! Well width in nm
    real(dp) :: m = 0.067_dp ! GaAs effective mass
    real(dp), allocatable :: energies(:)
    real(dp), allocatable :: analytical(:)
    integer :: n
    
    ! Analytical solution: En = (n²π²ℏ²)/(2mL²)
    allocate(analytical(5))
    do n = 1, 5
        analytical(n) = (n**2 * PI**2 * HBAR**2)/(2.0_dp * m * L**2)
    end do
    
    call solve_schrodinger(L, m, energies)
    @assertEqual(energies(1:5), analytical, tolerance=1.0e-6_dp)
end subroutine
```

### 12.4 Regression Testing

#### 1. Reference Data Management
```fortran
! tests/regression/test_regression.pf
module reference_data
    use iso_fortran_env, only: dp => real64
    implicit none
    
    type :: ReferenceDataType
        real(dp), allocatable :: k_points(:,:)
        real(dp), allocatable :: energies(:,:)
        real(dp), allocatable :: wavefunctions(:,:,:)
    contains
        procedure :: load_reference
        procedure :: compare_with_current
    end type
end module
```

### 12.5 Test Organization

```
tests/
├── CMakeLists.txt
├── unit/                     # Basic unit tests
│   ├── test_hamiltonian.pf
│   ├── test_fd.pf
│   └── test_materials.pf
├── validation/              # Physics validation
│   ├── test_simple_well.pf
│   ├── test_band_structure.pf
│   └── test_experimental.pf
├── regression/             # Regression tests
│   ├── test_regression.pf
│   └── reference_data/
├── performance/           # Performance tests
│   └── test_performance.pf
└── utils/                # Test utilities
    ├── assert_utils.f90
    └── performance_utils.f90
```

### 12.6 Testing Standards

1. **Unit Test Requirements**:
   - Every module must have corresponding unit tests
   - Tests must cover both normal and edge cases
   - All public interfaces must be tested
   - Error conditions must be verified

2. **Physics Validation**:
   - Compare with analytical solutions where possible
   - Verify physical symmetries and conservation laws
   - Test boundary conditions
   - Validate against experimental data

3. **Performance Testing**:
   - Measure scaling with system size
   - Monitor memory usage
   - Check parallel efficiency
   - Compare against reference implementations

4. **Test Coverage Goals**:
   - Line coverage: >90%
   - Branch coverage: >85%
   - Function coverage: 100%

### 12.7 Continuous Integration

Add to GitHub Actions workflow:
```yaml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    
    steps:
    - uses: actions/checkout@v2
    
    - name: Setup Fortran
      uses: fortran-lang/setup-fortran@v1
      
    - name: Install pFUnit
      run: |
        git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit
        cd pFUnit
        mkdir build && cd build
        cmake .. -DSKIP_MPI=YES
        make -j2
        sudo make install
        
    - name: Build and Test
      run: |
        mkdir build && cd build
        cmake .. -DBUILD_TESTING=ON
        make -j2
        ctest --output-on-failure
``` 

## 13. Performance Optimization Framework

### 13.1 Profiling Infrastructure

#### CMake Integration
```cmake
# Performance profiling options
option(ENABLE_PROFILING "Enable performance profiling" OFF)
option(USE_VTUNE "Use Intel VTune profiler" OFF)
option(USE_GPROF "Use GNU gprof profiler" OFF)

if(ENABLE_PROFILING)
    if(USE_VTUNE)
        find_package(VTune REQUIRED)
        add_compile_options(-g -debug inline-debug-info)
        link_libraries(${VTune_LIBRARIES})
    endif()
    
    if(USE_GPROF)
        add_compile_options(-pg)
        add_link_options(-pg)
    endif()
endif()
```

#### Performance Measurement Module
```fortran
! Core/Performance/profiling.f90
module profiling
    use iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: Timer
    public :: PerformanceMetrics
    public :: start_profiling
    public :: stop_profiling

    type :: Timer
        private
        real(dp) :: start_time = 0.0_dp
        real(dp) :: total_time = 0.0_dp
        integer :: call_count = 0
    contains
        procedure :: start
        procedure :: stop
        procedure :: get_average_time
    end type

    type :: PerformanceMetrics
        type(Timer) :: hamiltonian_construction
        type(Timer) :: diagonalization
        type(Timer) :: fd_operations
        type(Timer) :: io_operations
        integer(8) :: peak_memory = 0
        integer :: omp_threads = 1
        integer :: mkl_threads = 1
    contains
        procedure :: report
        procedure :: save_to_file
    end type
end module
```

### 13.2 OpenMP Integration

#### Thread Management
```fortran
! Core/Parallel/thread_control.f90
module thread_control
    use omp_lib
    use mkl_service
    implicit none
    private

    public :: initialize_threading
    public :: set_thread_count
    public :: get_optimal_thread_count

contains
    subroutine initialize_threading(num_threads, error)
        integer, intent(in) :: num_threads
        type(ErrorType), intent(out) :: error
        
        ! Set OpenMP threads
        call omp_set_num_threads(num_threads)
        
        ! Set MKL threads
        call mkl_set_num_threads(num_threads)
        
        ! Set MKL threading layer
        call mkl_set_threading_layer(MKL_THREADING_INTEL)
    end subroutine

    function get_optimal_thread_count() result(num_threads)
        integer :: num_threads
        integer :: num_physical_cores
        
        ! Get number of physical cores
        num_physical_cores = omp_get_num_procs()
        
        ! Use physical core count for optimal performance
        num_threads = num_physical_cores
    end function
end module
```

### 13.3 MKL Optimization

#### MKL Configuration Module
```fortran
! Numerics/LinearAlgebra/mkl_config.f90
module mkl_config
    use mkl_service
    implicit none
    private

    public :: configure_mkl
    public :: MKLConfig

    type :: MKLConfig
        integer :: max_threads
        integer :: dynamic_threads
        integer :: threading_layer
        integer :: affinity
    contains
        procedure :: apply
        procedure :: optimize_for_size
    end type

contains
    subroutine configure_mkl(matrix_size, error)
        integer, intent(in) :: matrix_size
        type(ErrorType), intent(out) :: error
        type(MKLConfig) :: config

        ! Optimize settings based on matrix size
        if (matrix_size < 1000) then
            call config%optimize_for_size('small')
        else
            call config%optimize_for_size('large')
        end if

        ! Apply configuration
        call config%apply()
    end subroutine
end module
```

### 13.4 Performance-Critical Implementations

#### Optimized Matrix Operations
```fortran
! Numerics/LinearAlgebra/matrix_ops.f90
module matrix_ops
    use mkl_service
    use omp_lib
    implicit none
    private

    public :: construct_hamiltonian_parallel
    public :: apply_finite_difference_parallel

contains
    subroutine construct_hamiltonian_parallel(H, k_points, materials)
        complex(dp), intent(out) :: H(:,:)
        real(dp), intent(in) :: k_points(:)
        type(MaterialType), intent(in) :: materials(:)
        integer :: i, j, chunk_size
        
        chunk_size = size(H,1) / omp_get_max_threads()
        
        !$omp parallel do schedule(static, chunk_size)
        do i = 1, size(H,1)
            do j = 1, size(H,2)
                call compute_matrix_element(H(i,j), i, j, k_points, materials)
            end do
        end do
        !$omp end parallel do
    end subroutine
end module
```

### 13.5 Performance Monitoring

#### Performance Test Suite
```fortran
! tests/performance/perf_test_suite.pf
module performance_tests
    use profiling
    implicit none

    @test
    subroutine test_hamiltonian_scaling()
        type(PerformanceMetrics) :: metrics
        integer :: sizes(4) = [100, 200, 400, 800]
        real(dp) :: times(4)
        
        do i = 1, 4
            call metrics%hamiltonian_construction%start()
            call construct_test_hamiltonian(sizes(i))
            call metrics%hamiltonian_construction%stop()
            times(i) = metrics%hamiltonian_construction%get_average_time()
        end do
        
        ! Verify O(N²) scaling
        @assertTrue(verify_scaling(sizes, times, 2.0_dp))
    end subroutine
end module
```

### 13.6 Memory Optimization

#### Memory Management Module
```fortran
! Core/Memory/memory_manager.f90
module memory_manager
    use iso_fortran_env
    implicit none
    private

    public :: track_allocation
    public :: report_memory_usage
    public :: optimize_workspace

    type :: MemoryStats
        integer(8) :: current_usage = 0
        integer(8) :: peak_usage = 0
        integer :: allocation_count = 0
    contains
        procedure :: update
        procedure :: report
    end type

contains
    subroutine optimize_workspace(matrix_size, error)
        integer, intent(in) :: matrix_size
        type(ErrorType), intent(out) :: error
        integer(8) :: optimal_size
        
        ! Calculate optimal workspace size
        optimal_size = estimate_workspace(matrix_size)
        
        ! Allocate workspace with error checking
        call allocate_workspace(optimal_size, error)
    end subroutine
end module
```

### 13.7 Performance Guidelines

1. **Matrix Operations**:
   - Use blocked algorithms for large matrices
   - Minimize temporary array allocations
   - Leverage BLAS Level 3 operations

2. **Memory Management**:
   - Pre-allocate workspace arrays
   - Use contiguous memory layouts
   - Implement memory pooling for frequent allocations

3. **Parallelization Strategy**:
   - Use OpenMP for matrices > 500x500
   - Set optimal chunk sizes for load balancing
   - Minimize thread synchronization points

4. **MKL Usage**:
   - Match thread count to physical cores
   - Use appropriate interface (LP64/ILP64)
   - Enable vectorization with aligned data

5. **I/O Optimization**:
   - Buffer large writes
   - Use binary format for large datasets
   - Implement parallel I/O for large systems

### 13.8 Performance Monitoring Tools

Add to CMakeLists.txt:
```cmake
# Performance monitoring targets
if(ENABLE_PROFILING)
    add_custom_target(perf_report
        COMMAND ${Python_EXECUTABLE} 
                ${CMAKE_SOURCE_DIR}/scripts/analyze_performance.py
                ${CMAKE_BINARY_DIR}/perf_data
        COMMENT "Generating performance report"
        VERBATIM
    )
endif()
``` 

## 14. Input/Output Modernization

### 14.1 HDF5 Integration

#### CMake Configuration
```cmake
# Find HDF5 Fortran
find_package(HDF5 REQUIRED COMPONENTS Fortran)
include_directories(${HDF5_INCLUDE_DIRS})
link_libraries(${HDF5_LIBRARIES})

# Optional compression support
option(USE_COMPRESSION "Enable HDF5 compression" ON)
if(USE_COMPRESSION)
    add_definitions(-DUSE_COMPRESSION)
endif()
```

#### HDF5 Output Module
```fortran
! IO/HDF5/hdf5_writer.f90
module hdf5_writer
    use hdf5
    use iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: HDF5Writer
    public :: initialize_hdf5
    public :: finalize_hdf5

    type :: HDF5Writer
        private
        integer(hid_t) :: file_id
        integer(hid_t) :: group_id
        logical :: is_initialized = .false.
    contains
        procedure :: create_file
        procedure :: write_eigenvalues
        procedure :: write_wavefunctions
        procedure :: write_attributes
        procedure :: close
    end type

contains
    subroutine write_eigenvalues(self, k_points, energies, error)
        class(HDF5Writer), intent(inout) :: self
        real(dp), intent(in) :: k_points(:,:)
        real(dp), intent(in) :: energies(:,:)
        type(ErrorType), intent(out) :: error
        integer(hid_t) :: dset_id
        integer(hsize_t) :: dims(2)

        dims = shape(energies)
        
        ! Create dataset with compression
        call h5dcreate_f(self%group_id, "eigenvalues", H5T_NATIVE_DOUBLE, &
                        dims, dset_id, error%code, &
                        chunk_size=min(dims, [100_hsize_t, 100_hsize_t]))
        
        ! Write data
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, energies, dims, error%code)
        
        ! Add attributes
        call self%write_attributes(dset_id, "units", "eV")
        call self%write_attributes(dset_id, "description", "Band energies")
    end subroutine
end module
```

### 14.2 JSON Input Processing

#### FSON Integration
```fortran
! IO/Input/json_reader.f90
module json_reader
    use fson
    use fson_value_m
    use iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: InputReader
    public :: read_input_file

    type :: InputReader
        private
        type(fson_value), pointer :: root => null()
    contains
        procedure :: read_calculation_parameters
        procedure :: read_material_parameters
        procedure :: read_numerical_parameters
        procedure :: validate_input
    end type

contains
    subroutine read_calculation_parameters(self, params, error)
        class(InputReader), intent(in) :: self
        type(CalculationParameters), intent(out) :: params
        type(ErrorType), intent(out) :: error
        
        ! Read calculation type
        call fson_get(self%root, "calculation.type", params%calc_type)
        
        ! Read k-points
        call fson_get(self%root, "calculation.k_points.max", params%k_max)
        call fson_get(self%root, "calculation.k_points.steps", params%k_steps)
        
        ! Validate parameters
        if (.not. self%validate_input(error)) return
    end subroutine
end module
```

### 14.3 Logging System

#### Logger Module
```fortran
! IO/Logging/logger.f90
module logger
    use iso_fortran_env, only: error_unit, output_unit
    implicit none
    private

    public :: Logger
    public :: LogLevel

    integer, parameter :: MAX_MSG_LEN = 1024

    type :: LogLevel
        integer :: level
        character(len=10) :: name
    end type

    type(LogLevel), parameter :: &
        LOG_DEBUG = LogLevel(0, "DEBUG"), &
        LOG_INFO  = LogLevel(1, "INFO"), &
        LOG_WARN  = LogLevel(2, "WARNING"), &
        LOG_ERROR = LogLevel(3, "ERROR")

    type :: Logger
        private
        character(len=:), allocatable :: log_file
        type(LogLevel) :: min_level = LOG_INFO
        logical :: console_output = .true.
        integer :: file_unit = -1
    contains
        procedure :: log
        procedure :: debug
        procedure :: info
        procedure :: warn
        procedure :: error
        procedure :: set_level
        procedure :: open_log_file
        procedure :: close_log_file
    end type

contains
    subroutine log(self, level, message, source)
        class(Logger), intent(in) :: self
        type(LogLevel), intent(in) :: level
        character(len=*), intent(in) :: message
        character(len=*), intent(in), optional :: source
        character(len=20) :: timestamp
        
        if (level%level < self%min_level%level) return
        
        call get_timestamp(timestamp)
        
        if (self%console_output) then
            write(output_unit, '(A,1X,A,1X,A,": ",A)') &
                timestamp, level%name, source, message
        end if
        
        if (self%file_unit > 0) then
            write(self%file_unit, '(A,1X,A,1X,A,": ",A)') &
                timestamp, level%name, source, message
        end if
    end subroutine
end module
```

### 14.4 Progress Tracking

#### Progress Module
```fortran
! IO/Progress/progress_bar.f90
module progress_bar
    use iso_fortran_env, only: output_unit
    implicit none
    private

    public :: ProgressBar
    public :: Checkpoint

    type :: ProgressBar
        private
        integer :: total_steps
        integer :: current_step = 0
        character(len=:), allocatable :: description
        real :: start_time
        real :: last_update
    contains
        procedure :: update
        procedure :: finish
        procedure :: estimate_remaining
    end type

    type :: Checkpoint
        private
        character(len=:), allocatable :: filename
        integer :: checkpoint_interval
        logical :: enabled = .false.
    contains
        procedure :: save
        procedure :: load
        procedure :: should_checkpoint
    end type

contains
    subroutine update(self, step, force_update)
        class(ProgressBar), intent(inout) :: self
        integer, intent(in) :: step
        logical, intent(in), optional :: force_update
        real :: progress, elapsed, remaining
        
        self%current_step = step
        progress = real(step) / self%total_steps
        
        if (should_update() .or. force_update) then
            elapsed = get_elapsed_time()
            remaining = self%estimate_remaining(elapsed, progress)
            
            call display_progress_bar(progress, elapsed, remaining)
        end if
    end subroutine
end module
```

### 14.5 Checkpointing System

#### Checkpoint Module
```fortran
! IO/Checkpoint/checkpoint_manager.f90
module checkpoint_manager
    use hdf5
    implicit none
    private

    public :: CheckpointManager

    type :: CheckpointManager
        private
        character(len=:), allocatable :: checkpoint_dir
        integer :: checkpoint_interval
        logical :: enabled = .false.
    contains
        procedure :: save_state
        procedure :: restore_state
        procedure :: cleanup_old_checkpoints
    end type

contains
    subroutine save_state(self, state, iteration, error)
        class(CheckpointManager), intent(in) :: self
        type(SimulationState), intent(in) :: state
        integer, intent(in) :: iteration
        type(ErrorType), intent(out) :: error
        character(len=:), allocatable :: filename
        
        if (.not. self%enabled) return
        
        ! Create checkpoint filename
        filename = generate_checkpoint_name(iteration)
        
        ! Save state to HDF5 file
        call save_to_hdf5(filename, state, error)
        
        ! Cleanup old checkpoints if needed
        call self%cleanup_old_checkpoints()
    end subroutine
end module
```

### 14.6 I/O Guidelines

1. **Data Storage Strategy**:
   - Use HDF5 for large datasets
   - Implement compression for wavefunctions
   - Store metadata as attributes
   - Use chunking for large arrays

2. **Input Processing**:
   - Validate all input parameters
   - Provide clear error messages
   - Support both JSON and legacy formats
   - Include default values

3. **Logging Best Practices**:
   - Log all significant events
   - Include timestamps
   - Use appropriate log levels
   - Rotate log files

4. **Progress Reporting**:
   - Update at reasonable intervals
   - Show estimated completion time
   - Include memory usage
   - Support quiet mode

5. **Checkpointing Strategy**:
   - Regular interval saves
   - Cleanup old checkpoints
   - Verify checkpoint integrity
   - Automatic recovery

### 14.7 File Format Specifications

#### HDF5 Structure
```
simulation_output.h5
├── metadata/
│   ├── version
│   ├── timestamp
│   └── parameters
├── results/
│   ├── eigenvalues
│   │   ├── k_points
│   │   └── energies
│   └── wavefunctions
│       ├── real_part
│       └── imag_part
└── checkpoints/
    └── iteration_XXXXX/
```

#### JSON Input Format
```json
{
    "calculation": {
        "type": "band_structure",
        "k_points": {
            "direction": "kx",
            "max": 0.1,
            "steps": 100
        }
    },
    "materials": [
        {
            "name": "GaAs",
            "parameters": {
                "band_gap": 1.519,
                "effective_mass": 0.067
            }
        }
    ],
    "numerical": {
        "fd_steps": 1000,
        "energy_tolerance": 1e-6
    }
}
``` 

## 15. Code Quality Framework

### 15.1 Static Analysis Integration

#### CMake Configuration
```cmake
# Static analysis options
option(ENABLE_STATIC_ANALYSIS "Enable static analysis" ON)
option(USE_FORTRAN_LINTER "Use fortran-linter for static analysis" ON)

if(ENABLE_STATIC_ANALYSIS)
    find_program(FORTRAN_LINTER flint)
    if(FORTRAN_LINTER)
        add_custom_target(lint
            COMMAND ${FORTRAN_LINTER}
                    --recursive
                    --config=${CMAKE_SOURCE_DIR}/.flint.yml
                    ${CMAKE_SOURCE_DIR}/src
            COMMENT "Running fortran-linter"
            VERBATIM
        )
    endif()
endif()

# Compiler warnings
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra -Wno-maybe-uninitialized")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wconversion -Wunderflow")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fcheck=all")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -warn all")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check all")
endif()
```

### 15.2 Linter Configuration

#### Fortran Linter Config
```yaml
# .flint.yml
rules:
  # Line length
  max_line_length: 100
  
  # Naming conventions
  module_naming: snake_case
  subroutine_naming: snake_case
  variable_naming: snake_case
  type_naming: pascal_case
  
  # Code structure
  max_module_length: 300
  max_subroutine_length: 50
  max_nesting_depth: 4
  
  # Documentation
  require_implicit_none: true
  require_module_doc: true
  require_subroutine_doc: true
  
  # Modern Fortran
  prefer_modern_operators: true
  avoid_common_blocks: true
  avoid_goto: true
  
  # Memory management
  check_allocate_status: true
  check_deallocate_status: true
  
  # Formatting
  indentation: 4
  require_consistent_spacing: true
```

### 15.3 Code Review Standards

#### 1. Code Style Guide
```fortran
! Example of preferred code style
module example_module
    !! Module description (required)
    use iso_fortran_env, only: dp => real64
    implicit none
    private  ! Default private visibility
    
    ! Public interface
    public :: calculate_energy
    
    ! Constants with clear names
    real(dp), parameter :: BOLTZMANN_CONSTANT = 8.617333262e-5_dp  ! eV/K
    
    ! Types with clear structure
    type :: CalculationParameters
        real(dp) :: temperature      ! Temperature in Kelvin
        real(dp) :: cutoff_energy    ! Energy cutoff in eV
        logical  :: include_strain   ! Whether to include strain effects
    contains
        procedure :: validate => validate_parameters
    end type
    
contains
    function calculate_energy(params) result(energy)
        !! Calculate system energy based on input parameters
        type(CalculationParameters), intent(in) :: params
        real(dp) :: energy
        
        ! Input validation
        if (.not. params%validate()) then
            error stop "Invalid calculation parameters"
        end if
        
        ! Clear algorithm structure with comments
        energy = 0.0_dp
        
        ! Use descriptive intermediate variables
        real(dp) :: thermal_contribution
        thermal_contribution = BOLTZMANN_CONSTANT * params%temperature
        
        energy = energy + thermal_contribution
    end function
end module
```

### 15.4 Code Quality Metrics

#### 1. Complexity Metrics
```fortran
! metrics.f90
module code_metrics
    implicit none
    private
    
    public :: calculate_complexity
    public :: report_metrics
    
    type :: CodeMetrics
        integer :: cyclomatic_complexity
        integer :: nesting_depth
        integer :: number_of_statements
        integer :: number_of_parameters
    contains
        procedure :: is_within_limits
        procedure :: generate_report
    end type
    
contains
    function is_within_limits(self) result(is_valid)
        class(CodeMetrics), intent(in) :: self
        logical :: is_valid
        
        is_valid = &
            self%cyclomatic_complexity <= 10 .and. &
            self%nesting_depth <= 4 .and. &
            self%number_of_statements <= 50 .and. &
            self%number_of_parameters <= 5
    end function
end module
```

### 15.5 Pre-commit Hooks

#### Git Pre-commit Configuration
```yaml
# .pre-commit-config.yaml
repos:
- repo: local
  hooks:
    - id: fortran-lint
      name: Fortran Linter
      entry: flint
      language: system
      files: \.(f90|F90)$
      
    - id: format-check
      name: Format Check
      entry: fprettify
      language: system
      files: \.(f90|F90)$
      
    - id: compile-check
      name: Compilation Check
      entry: scripts/compile_check.sh
      language: system
      files: \.(f90|F90)$
```

### 15.6 Quality Assurance Workflow

#### Continuous Integration Quality Checks
```yaml
# .github/workflows/quality.yml
name: Code Quality

on: [push, pull_request]

jobs:
  quality:
    runs-on: ubuntu-latest
    
    steps:
    - uses: actions/checkout@v2
    
    - name: Install Tools
      run: |
        pip install fortran-linter
        pip install fprettify
    
    - name: Lint Code
      run: |
        flint --recursive src/
        
    - name: Check Format
      run: |
        fprettify --check src/
        
    - name: Compile with Warnings
      run: |
        mkdir build && cd build
        cmake -DCMAKE_Fortran_FLAGS="-Wall -Wextra" ..
        make
```

### 15.7 Code Quality Guidelines

1. **Code Organization**:
   - Maximum 300 lines per module
   - Maximum 50 lines per subroutine/function
   - Maximum 4 levels of nesting
   - Clear separation of concerns

2. **Naming Conventions**:
   - Descriptive, self-documenting names
   - Consistent case style (snake_case for variables/procedures)
   - Clear prefix/suffix for derived types
   - Avoid abbreviations unless standard

3. **Documentation Requirements**:
   - Module purpose and usage
   - Procedure inputs/outputs
   - Algorithm description
   - Performance considerations
   - Error conditions

4. **Error Handling**:
   - Check all allocations
   - Validate all inputs
   - Clear error messages
   - Proper cleanup on error

5. **Performance Considerations**:
   - Minimize temporary allocations
   - Use appropriate data structures
   - Consider cache efficiency
   - Document performance expectations

### 15.8 Quality Monitoring

#### Quality Dashboard
```fortran
! tools/quality/dashboard.f90
module quality_dashboard
    implicit none
    private
    
    public :: generate_quality_report
    
    type :: QualityMetrics
        integer :: total_lines
        integer :: comment_percentage
        integer :: test_coverage
        integer :: complexity_score
        integer :: warning_count
    contains
        procedure :: generate_html_report
        procedure :: check_thresholds
    end type
    
contains
    subroutine generate_quality_report(output_file)
        character(len=*), intent(in) :: output_file
        type(QualityMetrics) :: metrics
        
        ! Collect metrics
        call collect_code_metrics(metrics)
        
        ! Generate HTML report
        call metrics%generate_html_report(output_file)
        
        ! Check against thresholds
        if (.not. metrics%check_thresholds()) then
            error stop "Quality metrics below threshold"
        end if
    end subroutine
end module
``` 

## 16. Feature Extensions Framework

### 16.1 Material Database System

#### Material Parameter Module
```fortran
! Physics/Materials/material_database.f90
module material_database
    use iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: MaterialDatabase
    public :: MaterialParameters
    public :: load_material_database

    !> Material family enumeration
    type :: MaterialFamily
        integer :: id
        character(len=32) :: name
    end type

    type(MaterialFamily), parameter :: &
        ZINC_BLENDE = MaterialFamily(1, "zinc_blende"), &
        WURTZITE    = MaterialFamily(2, "wurtzite"), &
        ROCKSALT    = MaterialFamily(3, "rocksalt")

    !> Complete material parameters type
    type :: MaterialParameters
        ! Basic parameters
        real(dp) :: band_gap           ! Fundamental band gap [eV]
        real(dp) :: spin_orbit_split   ! Spin-orbit splitting [eV]
        real(dp) :: crystal_field      ! Crystal field splitting [eV]
        real(dp) :: lattice_constant   ! Lattice constant [Å]
        
        ! k·p parameters
        real(dp) :: ep      ! Kane energy parameter [eV]
        real(dp) :: f       ! Luttinger-like parameter
        real(dp) :: gamma1  ! Luttinger parameter
        real(dp) :: gamma2  ! Luttinger parameter
        real(dp) :: gamma3  ! Luttinger parameter
        
        ! Strain parameters
        real(dp) :: ac      ! Conduction band deformation potential [eV]
        real(dp) :: av      ! Valence band deformation potential [eV]
        real(dp) :: b       ! Shear deformation potential [eV]
        real(dp) :: d       ! Shear deformation potential [eV]
        
        ! Elastic constants
        real(dp) :: c11     ! Elastic constant [GPa]
        real(dp) :: c12     ! Elastic constant [GPa]
        real(dp) :: c44     ! Elastic constant [GPa]
        
        ! Piezoelectric constants
        real(dp) :: e14     ! Piezoelectric constant [C/m²]
        
        ! Temperature dependence
        real(dp) :: alpha   ! Band gap temperature coefficient [eV/K]
        real(dp) :: beta    ! Band gap temperature coefficient [K]
    contains
        procedure :: validate
        procedure :: interpolate_alloy
    end type

    !> Database management type
    type :: MaterialDatabase
        private
        type(MaterialParameters), allocatable :: materials(:)
        character(len=32), allocatable :: material_names(:)
    contains
        procedure :: load_from_file
        procedure :: get_parameters
        procedure :: add_material
        procedure :: get_alloy_parameters
    end type

contains
    subroutine load_from_file(self, filename, error)
        class(MaterialDatabase), intent(inout) :: self
        character(len=*), intent(in) :: filename
        type(ErrorType), intent(out) :: error
        type(HDF5Reader) :: reader
        
        ! Open database file
        call reader%open(filename, error)
        if (error%has_error()) return
        
        ! Read material parameters
        call reader%read_dataset("/materials", self%materials, error)
        if (error%has_error()) return
        
        ! Read material names
        call reader%read_dataset("/names", self%material_names, error)
    end subroutine

    function get_alloy_parameters(self, material1, material2, x, error) result(params)
        class(MaterialDatabase), intent(in) :: self
        character(len=*), intent(in) :: material1, material2
        real(dp), intent(in) :: x  ! Composition
        type(ErrorType), intent(out) :: error
        type(MaterialParameters) :: params
        
        ! Get base materials
        type(MaterialParameters) :: params1, params2
        
        call self%get_parameters(material1, params1, error)
        if (error%has_error()) return
        
        call self%get_parameters(material2, params2, error)
        if (error%has_error()) return
        
        ! Interpolate parameters
        call params%interpolate_alloy(params1, params2, x)
    end function
end module
```

### 16.2 Self-Consistent Poisson Solver

#### Poisson Module
```fortran
! Physics/Electrostatics/poisson_solver.f90
module poisson_solver
    use iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: PoissonSolver
    public :: solve_poisson_equation

    !> Solver configuration type
    type :: PoissonConfig
        real(dp) :: tolerance = 1.0e-6_dp
        integer :: max_iterations = 1000
        real(dp) :: mixing_parameter = 0.3_dp
        character(len=32) :: boundary_condition = "dirichlet"
    contains
        procedure :: validate
    end type

    !> Main solver type
    type :: PoissonSolver
        private
        type(PoissonConfig) :: config
        real(dp), allocatable :: potential(:)
        real(dp), allocatable :: charge_density(:)
        real(dp), allocatable :: permittivity(:)
    contains
        procedure :: initialize
        procedure :: solve
        procedure :: update_charge_density
        procedure :: get_electric_field
    end type

contains
    subroutine solve(self, error)
        class(PoissonSolver), intent(inout) :: self
        type(ErrorType), intent(out) :: error
        real(dp), allocatable :: rhs(:), new_potential(:)
        real(dp) :: residual
        integer :: iter
        
        ! Allocate work arrays
        allocate(rhs(size(self%potential)), &
                new_potential(size(self%potential)))
        
        ! Iterative solution
        do iter = 1, self%config%max_iterations
            ! Construct right-hand side
            call construct_rhs(self%charge_density, self%permittivity, rhs)
            
            ! Solve linear system
            call solve_linear_system(rhs, new_potential, error)
            if (error%has_error()) return
            
            ! Check convergence
            residual = norm2(new_potential - self%potential)
            if (residual < self%config%tolerance) exit
            
            ! Update potential with mixing
            self%potential = (1.0_dp - self%config%mixing_parameter) * &
                           self%potential + &
                           self%config%mixing_parameter * new_potential
        end do
        
        ! Check if converged
        if (iter >= self%config%max_iterations) then
            call error%set_error(ERR_NUMERICS, ERROR_WARNING, &
                               "Poisson solver did not converge")
        end if
    end subroutine
end module
```

### 16.3 Extended Physics Effects

#### Strain Module
```fortran
! Physics/Strain/strain_effects.f90
module strain_effects
    use iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: StrainTensor
    public :: calculate_strain
    public :: apply_strain_hamiltonian

    !> Strain tensor type
    type :: StrainTensor
        real(dp) :: xx, xy, xz
        real(dp) :: yx, yy, yz
        real(dp) :: zx, zy, zz
    contains
        procedure :: get_hydrostatic
        procedure :: get_biaxial
        procedure :: get_shear
    end type

    !> Strain calculation type
    type :: StrainCalculator
        private
        type(MaterialParameters) :: params
        logical :: include_piezo = .false.
    contains
        procedure :: calculate_strain_tensor
        procedure :: get_strain_hamiltonian
    end type

contains
    subroutine calculate_strain_tensor(self, a_substrate, a_layer, strain)
        class(StrainCalculator), intent(in) :: self
        real(dp), intent(in) :: a_substrate  ! Substrate lattice constant
        real(dp), intent(in) :: a_layer      ! Layer lattice constant
        type(StrainTensor), intent(out) :: strain
        
        ! Calculate in-plane strain
        real(dp) :: eps_xx
        eps_xx = (a_substrate - a_layer) / a_layer
        
        ! Fill strain tensor
        strain%xx = eps_xx
        strain%yy = eps_xx
        strain%zz = -2.0_dp * (self%params%c12/self%params%c11) * eps_xx
        
        ! Off-diagonal elements are zero for biaxial strain
        strain%xy = 0.0_dp
        strain%xz = 0.0_dp
        strain%yz = 0.0_dp
    end subroutine
end module
```

#### Spin-Orbit Coupling
```fortran
! Physics/SpinOrbit/spin_orbit_coupling.f90
module spin_orbit_coupling
    use iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: SpinOrbitParameters
    public :: apply_spin_orbit

    !> Spin-orbit parameters type
    type :: SpinOrbitParameters
        real(dp) :: delta_so  ! Spin-orbit splitting
        real(dp) :: lambda    ! Spin-orbit coupling constant
    contains
        procedure :: get_matrix_elements
    end type

contains
    subroutine apply_spin_orbit(H, params, error)
        complex(dp), intent(inout) :: H(:,:)
        type(SpinOrbitParameters), intent(in) :: params
        type(ErrorType), intent(out) :: error
        
        ! Add spin-orbit coupling terms to Hamiltonian
        call add_so_diagonal(H, params)
        call add_so_off_diagonal(H, params)
        
        ! Verify Hermiticity
        if (.not. is_hermitian(H)) then
            call error%set_error(ERR_PHYSICS, ERROR_FATAL, &
                               "Non-Hermitian Hamiltonian after SO coupling")
        end if
    end subroutine
end module
```

### 16.4 Visualization System

#### Plotting Interface
```fortran
! IO/Visualization/plot_interface.f90
module plot_interface
    use iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: BandPlotter
    public :: export_for_plotting

    !> Plot configuration type
    type :: PlotConfig
        character(len=32) :: backend = "matplotlib"
        character(len=32) :: style = "scientific"
        logical :: show_grid = .true.
        logical :: show_legend = .true.
        character(len=32) :: output_format = "pdf"
    contains
        procedure :: validate
        procedure :: apply_style
    end type

    !> Band structure plotter
    type :: BandPlotter
        private
        type(PlotConfig) :: config
        real(dp), allocatable :: k_points(:,:)
        real(dp), allocatable :: energies(:,:)
        character(len=:), allocatable :: title
        character(len=:), allocatable :: k_labels(:)
    contains
        procedure :: set_data
        procedure :: add_labels
        procedure :: export_matplotlib
        procedure :: export_plotly
        procedure :: save
    end type

contains
    subroutine export_matplotlib(self, filename, error)
        class(BandPlotter), intent(in) :: self
        character(len=*), intent(in) :: filename
        type(ErrorType), intent(out) :: error
        
        ! Create Python script for plotting
        call create_matplotlib_script(self, "plot_bands.py", error)
        if (error%has_error()) return
        
        ! Execute Python script
        call execute_command_line("python plot_bands.py", exitstat=status)
        if (status /= 0) then
            call error%set_error(ERR_SYSTEM, ERROR_FATAL, &
                               "Failed to execute plotting script")
        end if
    end subroutine

    subroutine create_matplotlib_script(self, script_name, error)
        class(BandPlotter), intent(in) :: self
        character(len=*), intent(in) :: script_name
        type(ErrorType), intent(out) :: error
        integer :: unit
        
        open(newunit=unit, file=script_name, status="replace", action="write")
        
        write(unit, '(A)') "import numpy as np"
        write(unit, '(A)') "import matplotlib.pyplot as plt"
        write(unit, '(A)') "from matplotlib import rcParams"
        write(unit, '(A)') ""
        write(unit, '(A)') "# Set style"
        write(unit, '(A)') "plt.style.use('scientific')"
        write(unit, '(A)') ""
        write(unit, '(A)') "# Load data"
        write(unit, '(A)') "k_points = np.load('k_points.npy')"
        write(unit, '(A)') "energies = np.load('energies.npy')"
        write(unit, '(A)') ""
        write(unit, '(A)') "# Create plot"
        write(unit, '(A)') "fig, ax = plt.subplots(figsize=(8, 6))"
        write(unit, '(A)') ""
        write(unit, '(A)') "# Plot bands"
        write(unit, '(A)') "for band in range(energies.shape[1]):"
        write(unit, '(A)') "    ax.plot(k_points, energies[:, band])"
        write(unit, '(A)') ""
        write(unit, '(A)') "# Customize plot"
        write(unit, '(A)') "ax.set_xlabel('Wave Vector')"
        write(unit, '(A)') "ax.set_ylabel('Energy (eV)')"
        write(unit, '(A)') "ax.set_title('" // self%title // "')"
        write(unit, '(A)') ""
        write(unit, '(A)') "# Save plot"
        write(unit, '(A)') "plt.savefig('band_structure.pdf', dpi=300, bbox_inches='tight')"
        
        close(unit)
    end subroutine
end module
```

### 16.5 Feature Integration Guidelines

1. **Material System Integration**:
   - Use HDF5 for material database storage
   - Implement parameter validation
   - Support material interpolation
   - Include temperature dependence

2. **Self-Consistency Implementation**:
   - Iterative solution scheme
   - Adaptive mixing parameters
   - Convergence monitoring
   - Error estimation

3. **Physics Extensions**:
   - Modular effect implementation
   - Clear physics documentation
   - Validation against literature
   - Performance optimization

4. **Visualization Integration**:
   - Multiple backend support
   - Customizable plotting
   - Interactive options
   - Publication-quality output

### 16.6 Feature Testing Strategy

1. **Material Parameters**:
   - Validate against experimental data
   - Check interpolation schemes
   - Verify temperature dependence
   - Test database operations

2. **Self-Consistent Solutions**:
   - Verify convergence properties
   - Test boundary conditions
   - Check charge conservation
   - Validate field calculations

3. **Extended Physics**:
   - Compare with analytical solutions
   - Verify symmetry properties
   - Test coupling terms
   - Validate strain effects

4. **Visualization**:
   - Check plot accuracy
   - Verify style applications
   - Test export formats
   - Validate interactivity
