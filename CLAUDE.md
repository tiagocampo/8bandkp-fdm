# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Development Commands

### Build System
```bash
# Build both executables (bandStructure and gfactorCalculation)
make all

# Build specific executables
make bandStructure
make gfactor

# Clean build artifacts
make clean_all

# Run calculations
./bandStructure < input.cfg
./gfactorCalculation < gfactor.example

# Validate results
./scripts/validate_results.sh
./scripts/verify_build.sh

# Generate plots
./scripts/plot_all.sh
```

### Testing
```bash
# Run validation tests
./scripts/validate_results.sh

# Test with example inputs
./bandStructure < examples/bulk_InAs60Sb40.example
./gfactorCalculation < examples/gfactor_quantum_well.example
```

## Architecture Overview

This is a Fortran-based scientific computing package that implements the 8-band k·p method for semiconductor band structure and g-factor calculations using finite difference discretization.

### Core Scientific Workflow
1. **Input parsing** → **Material parameter loading** → **Hamiltonian construction**
2. **Finite difference discretization** → **Sparse matrix diagonalization** → **Eigenvalue extraction**
3. **Output formatting** → **File generation** → **Visualization**

### Key Modules

**`src/defs.f90`** - Fundamental definitions and physical constants
**`src/parameters.f90`** - Material parameter database for III-V semiconductors
**`src/hamiltonianConstructor.f90`** - 8×8 k·p Hamiltonian matrix construction
**`src/finitedifferences.f90`** - Finite difference discretization matrices
**`src/gfactor_functions.f90`** - Landé g-factor calculations using Löwdin partitioning
**`src/main.f90`** - Band structure calculator primary workflow
**`src/main_gfactor.f90`** - G-factor specialized calculator

### Dual Executable Architecture
- **`bandStructure`**: Calculates band structures for bulk and quantum well systems
- **`gfactorCalculation`**: Specialized for Landé g-factor tensor calculations

### Data Dependencies
- **Intel MKL 2025.0+**: Required for sparse matrix operations and BLAS/LAPACK
- **FFTW3**: Fast Fourier transforms
- **Alternative fallback**: Standard BLAS/LAPACK + FFTW3 (comment out MKL in Makefile)

### Input Format
Human-readable configuration files with sections for:
- Wave vector parameters (direction, range, discretization)
- Confinement geometry (bulk vs quantum well)
- Material layers and band offsets
- Numerical discretization (FDstep, number of bands)
- External fields (electric/magnetic)

### Output Structure
```
outputs/<timestamp>/
├── eigenvalues.dat               # Band structure vs k-point
├── eigenfunctions_k_*.dat       # Wave function coefficients
├── parts.dat                    # Additional analysis data
└── run.meta                     # Run metadata
```

## Development Guidelines

### Scientific Accuracy Requirements
- All material parameters must be sourced from authoritative references (Vurgaftman et al., 2001)
- Mathematical formulations must be traceable to peer-reviewed literature
- Numerical methods must follow established finite difference formulations
- Results must be validated against analytical solutions or published benchmarks

### Code Standards
- Physical constants and approximations must be explicitly documented with citations
- Code must include inline comments explaining physical significance, not just implementation details
- Core physics modules must be independently testable with single, well-defined responsibilities
- New features must be added as separate modules without modifying existing core functionality

### Repository Structure
- `src/` - Core Fortran modules
- `examples/` - Input file examples for bulk, quantum well, and g-factor calculations
- `scripts/` - Build, validation, and plotting automation
- `docs/` - Complete user guides and API documentation
- `specs/` - Feature specifications using Speckit methodology
- `outputs/` - Calculation results (automatically timestamped)

### Performance Characteristics
- **Build Time**: < 5 minutes
- **Calculation Time**: < 30 seconds for test cases
- **Memory Usage**: < 1 GB for typical problems
- **Precision**: Supports single, double, and quadruple precision

### Key Physical Constraints
- 8-band k·p Hamiltonian includes spin-orbit coupling
- Supports bulk systems and quantum wells with up to 10 material layers
- Finite difference discretization with configurable grid density
- Material database includes III-V semiconductors and ternary alloys

## Agent Collaboration

This project uses Speckit for structured feature development:
- `/.specify/` - Speckit methodology templates and memory
- `/specs/` - Feature specifications with user stories and validation criteria
- Use `/speckit.plan` after specification clarification to create implementation plans
- Constitution (v1.1.0) mandates generic numerical perturbation methods over analytical approximations

## Validation Framework

Test cases include:
1. **Bulk InAs60Sb40**: Band gap validation (~0.11 eV expected)
2. **GaSb/InAs/AlSb Quantum Well**: Type-II band alignment verification
3. **G-factor calculations**: Landé g-factor tensor validation

Physical accuracy checks:
- Energy convergence within 1 meV
- Wave function normalization
- Band gap verification against known values
- Subband structure validation for quantum wells

## Active Technologies
- Fortran 2008+ (with gfortran 15.2.1+ or Intel ifort/ifx) + Intel MKL 2025.0+, ARPACK-ng, FFTW3, OpenMP (002-generic-numerical-perturbation)
- File-based I/O with structured outputs in timestamped directories (002-generic-numerical-perturbation)

## Recent Changes
- 002-generic-numerical-perturbation: Added Fortran 2008+ (with gfortran 15.2.1+ or Intel ifort/ifx) + Intel MKL 2025.0+, ARPACK-ng, FFTW3, OpenMP
