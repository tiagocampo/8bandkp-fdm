# Implementation Plan: Generic Numerical Perturbation Method for G-Factor Calculations

**Branch**: `002-generic-numerical-perturbation` | **Date**: 2025-11-02 | **Spec**: [link](spec.md)
**Input**: Feature specification from `/specs/002-generic-numerical-perturbation/spec.md`

## Summary

This feature implements a generic numerical perturbation method for Landé g-factor calculations in the 8-band k·p solver, replacing analytical Löwdin partitioning with robust numerical differentiation. The implementation will use central finite difference schemes to compute g-factors as numerical derivatives of eigenvalues with respect to magnetic field, enabling automatic capture of band mixing, anisotropy, and confinement effects for any band structure without requiring new analytical derivations.

## Technical Context

**Language/Version**: Fortran 2008+ (with gfortran 15.2.1+ or Intel ifort/ifx)
**Primary Dependencies**: Intel MKL 2025.0+, ARPACK-ng, FFTW3, OpenMP
**Storage**: File-based I/O with structured outputs in timestamped directories
**Testing**: Fortran unit tests with validation against analytical solutions and published experimental data
**Target Platform**: Linux HPC systems with Intel/AMD64 architecture
**Project Type**: Scientific computing library (extension to existing codebase)
**Performance Goals**: <30 seconds for standard quantum well g-factor calculations, <1GB memory usage for typical problems
**Constraints**: Must maintain backward compatibility with existing input/output formats, support up to 14-band k·p models, handle non-Hermitian Hamiltonians from magnetic field perturbations
**Scale/Scope**: Extension to existing Fortran codebase with ~2000 additional lines of code, supporting quantum well structures with up to 10 material layers and 1000 grid points

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

### Scientific Accuracy Gates
- [x] Mathematical formulations traceable to peer-reviewed literature (Winkler 2003, Tadjine et al. 2017, Chuang & Chang 1997)
- [x] Material parameters sourced from authoritative references (Vurgaftman et al. 2001 database already implemented)
- [x] Numerical methods follow established finite difference formulations (central difference schemes with proper convergence analysis)
- [x] Physical constants and approximations explicitly documented with citations (existing codebase includes comprehensive documentation)

### Reproducibility Gates
- [x] Input files self-contained with complete parameter specifications (existing input.cfg format maintained)
- [x] Build process deterministic and documented (existing Makefile with Intel MKL integration)
- [x] Output formats include sufficient metadata for verification (timestamped output directories with run metadata)
- [x] Results identical across different platforms/compiler versions (validated through existing framework)

### Modular Architecture Gates
- [x] Core physics modules independently testable (separate modules for Hamiltonian construction, finite differences, g-factor calculations)
- [x] Single, well-defined responsibility per module (clear separation of concerns in existing architecture)
- [x] Minimal and stable interfaces between modules (well-defined module interfaces in existing codebase)
- [x] New features added as separate modules without modifying core functionality (new numerical perturbation module will be independent)

### Validation Gates
- [x] Validation against analytical solutions or published benchmarks planned (requirement SC-002: 0.01% for bulk, 0.1% for quantum wells)
- [x] Test cases cover both bulk and quantum well configurations (existing examples and validation framework)
- [x] Performance benchmarks established (requirement SC-001: <30 seconds for standard structures)
- [x] Regression testing catches numerical accuracy degradation (existing validation scripts will be extended)

### Documentation Gates
- [x] Scientific methods documented with mathematical formulations and references (existing comprehensive documentation)
- [x] Code includes inline comments explaining physical significance (existing codebase follows this practice)
- [x] User documentation includes example input files and expected output ranges (existing examples and docs will be extended)
- [x] API documentation specifies input/output units and physical meaning (existing codebase includes proper documentation)

### Generic Numerical Perturbation Gates
- [x] Numerical methods preferred over analytical approximations for generality (core requirement of this feature)
- [x] Physical properties computed as numerical derivatives of eigenvalues (central finite difference implementation)
- [x] Hamiltonian construction separated from physics calculation (requirement FR-002 and existing modular design)
- [x] Band mixing, anisotropy, and confinement effects captured automatically (inherent in full Hamiltonian approach)
- [x] Central finite difference schemes with step size validation (requirement FR-001 and FR-008)
- [x] Results validated against analytical formulas in limiting cases (requirement FR-005 and existing validation framework)

**Phase 1 Re-evaluation**: All Constitution Check gates remain PASSED after design completion. The modular architecture, comprehensive validation framework, and generic numerical approach fully comply with project scientific standards and constitutional requirements.

## Project Structure

### Documentation (this feature)

```text
specs/002-generic-numerical-perturbation/
├── plan.md              # This file (/speckit.plan command output)
├── research.md          # Phase 0 output (/speckit.plan command)
├── data-model.md        # Phase 1 output (/speckit.plan command)
├── quickstart.md        # Phase 1 output (/speckit.plan command)
├── contracts/           # Phase 1 output (/speckit.plan command)
└── tasks.md             # Phase 2 output (/speckit.tasks command - NOT created by /speckit.plan)
```

### Source Code (repository root)

```text
src/
├── defs.f90                           # Fundamental definitions and constants
├── parameters.f90                     # Material parameter database
├── hamiltonianConstructor.f90         # Hamiltonian matrix construction
├── finitedifferences.f90              # Finite difference discretization
├── gfactor_functions.f90              # Current analytical g-factor calculations
├── numerical_perturbation.f90        # NEW: Generic numerical perturbation engine
├── perturbation_solver.f90            # NEW: Specialized eigenvalue solvers for perturbation
├── validation_framework.f90           # NEW: Validation against analytical solutions
├── outputFunctions.f90               # Results output and visualization
├── utils.f90                         # Utility functions
├── main.f90                          # Primary band structure calculator
└── main_gfactor.f90                  # G-factor specialized calculator

tests/
├── unit/
│   ├── test_numerical_perturbation.f90  # NEW: Unit tests for perturbation methods
│   ├── test_validation.f90              # NEW: Validation test suite
│   └── test_integration.f90             # NEW: Integration tests
└── examples/
    ├── bulk_gfactor_test.example        # NEW: Bulk material test cases
    ├── quantum_well_gfactor_test.example # NEW: Quantum well test cases
    └── convergence_test.example         # NEW: Step size validation tests

examples/                              # Extended with new test cases
├── numerical_perturbation_bulk.example     # NEW: Bulk semiconductor examples
├── numerical_perturbation_qw.example       # NEW: Quantum well examples
└── validation_comparison.example            # NEW: Validation against analytical results
```

**Structure Decision**: Single Fortran project with modular scientific computing architecture. The new numerical perturbation functionality will be implemented as independent modules (`numerical_perturbation.f90`, `perturbation_solver.f90`, `validation_framework.f90`) that integrate with the existing codebase through well-defined interfaces, maintaining the established modular architecture without modifying core functionality.

## Complexity Tracking

> **Fill ONLY if Constitution Check has violations that must be justified**

| Violation | Why Needed | Simpler Alternative Rejected Because |
|-----------|------------|-------------------------------------|
| [e.g., 4th project] | [current need] | [why 3 projects insufficient] |
| [e.g., Repository pattern] | [specific problem] | [why direct DB access insufficient] |
