# Implementation Plan: Build Verification and Plotting

**Branch**: `001-build-verification-plotting` | **Date**: 2025-01-27 | **Spec**: [spec.md](./spec.md)
**Input**: Feature specification from `/specs/001-build-verification-plotting/spec.md`

**Note**: This template is filled in by the `/speckit.plan` command. See `.specify/templates/commands/plan.md` for the execution workflow.

## Summary

Primary requirement: Verify build process works correctly, validate calculation results against known physics, and create standardized gnuplot scripts for result visualization. Technical approach: Document build dependencies, create validation test cases using GaSb/InAs/AlSb quantum well and bulk InAs60Sb40, and develop reusable plotting scripts for band structure, quantum well, and g-factor visualization.

## Technical Context

**Language/Version**: Fortran 90/95 (any compiler), gnuplot 5.0+  
**Primary Dependencies**: BLAS/LAPACK (any implementation), FFTW3, Make  
**Storage**: Plain text files (input configs, output data, plotting scripts)  
**Output Location**: All outputs under `outputs/<run-id>/` (isolated per run)  
**Testing**: Manual validation against analytical solutions, build verification scripts  
**Target Platform**: Linux/Unix systems with Fortran compiler and scientific libraries  
**Project Type**: Single scientific computing project  
**Performance Goals**: Build completes in <5 minutes, plots generate in <30 seconds  
**Constraints**: Must work with any Fortran compiler + BLAS/LAPACK combination  
**Scale/Scope**: Two executables (bandStructure, gfactorCalculation), three plotting script types

## Constitution Check (Post-Design)

*Re-evaluated after Phase 1 design completion*

### Scientific Accuracy Gates
- [x] Mathematical formulations traceable to peer-reviewed literature
- [x] Material parameters sourced from authoritative references  
- [x] Numerical methods follow established finite difference formulations
- [x] Physical constants and approximations explicitly documented with citations

### Reproducibility Gates  
- [x] Input files self-contained with complete parameter specifications
- [x] Build process deterministic and documented
- [x] Output formats include sufficient metadata for verification
- [x] Results identical across different platforms/compiler versions

### Modular Architecture Gates
- [x] Core physics modules independently testable
- [x] Single, well-defined responsibility per module
- [x] Minimal and stable interfaces between modules
- [x] New features added as separate modules without modifying core functionality

### Validation Gates
- [x] Validation against analytical solutions or published benchmarks planned
- [x] Test cases cover both bulk and quantum well configurations
- [x] Performance benchmarks established
- [x] Regression testing catches numerical accuracy degradation

### Documentation Gates
- [x] Scientific methods documented with mathematical formulations and references
- [x] Code includes inline comments explaining physical significance
- [x] User documentation includes example input files and expected output ranges
- [x] API documentation specifies input/output units and physical meaning

**Status**: ✅ All gates continue to pass after design phase

## Project Structure

### Documentation (this feature)

```text
specs/[###-feature]/
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
├── main.f90                    # Band structure calculation main
├── main_gfactor.f90           # G-factor calculation main
├── parameters.f90            # Material parameter definitions
├── hamiltonianConstructor.f90 # Hamiltonian matrix construction
├── finitedifferences.f90     # Finite difference methods
├── gfactor_functions.f90     # G-factor calculation functions
├── outputFunctions.f90       # Output formatting functions
├── mkl_sparse_handle.f90     # MKL sparse matrix handling
├── mkl_spblas.f90            # MKL sparse BLAS operations
└── utils.f90                 # Utility functions

scripts/
├── plot_band_structure.gp    # Band structure plotting script
├── plot_quantum_well.gp      # Quantum well visualization script
└── plot_gfactor.gp           # G-factor plotting script

docs/
├── BUILD.md                  # Build documentation
├── VALIDATION.md             # Result validation guide
└── PLOTTING.md               # Plotting script documentation

examples/
├── bulk.example              # Bulk calculation input
├── quantumwell.example       # Quantum well calculation input
└── gfactor.example           # G-factor calculation input
```

outputs/
└── <run-id>/                 # All run artifacts (e.g., eigenvalues.dat, eigenfunctions_*.dat, parts.dat)

**Input Handling Decision**: Executables accept an optional CLI filename argument for the configuration (default `input.cfg`), enabling deterministic selection of inputs in scripts and CI. Scripts MUST explicitly set the input file and output directory for reproducibility.

**Structure Decision**: Single Fortran project with modular source organization, dedicated plotting scripts, comprehensive documentation, and example input files for validation.

## Complexity Tracking

> **Fill ONLY if Constitution Check has violations that must be justified**

| Violation | Why Needed | Simpler Alternative Rejected Because |
|-----------|------------|-------------------------------------|
| [e.g., 4th project] | [current need] | [why 3 projects insufficient] |
| [e.g., Repository pattern] | [specific problem] | [why direct DB access insufficient] |
