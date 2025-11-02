# Feature Specification: Build Verification and Plotting

**Feature Branch**: `001-build-verification-plotting`  
**Created**: 2025-01-27  
**Status**: Draft  
**Input**: User description: "i want to review the code and make sure that the building process is working, well documented. I want to be able to run bulk and quantum well calculations and verify their results. I want to create some standardized and reusable gnuplot scripts to plot the results of the simulations to verify their accuracy"

## Clarifications

### Session 2025-01-27

- Q: Baseline toolchain targets for deterministic builds → A: Minimum: any Fortran + any BLAS/LAPACK; document example only
- Q: Verification baseline for "physically accurate" results → A: GaSb/InAs/AlSb quantum well + bulk InAs60Sb40
- Q: Output file formats for gnuplot compatibility → A: Plain text columns
- Q: Error handling scope for build failures → A: Clear error messages only
- Q: Documentation depth for build process → A: Step-by-step with examples

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Build Process Verification (Priority: P1)

A researcher needs to verify that the 8bandkp-fdm codebase builds correctly and produces reliable executables for both bulk and quantum well calculations.

**Why this priority**: Without a working build process, no calculations can be performed, making this the foundational requirement for all other functionality.

**Independent Test**: Can be fully tested by running `make all` and verifying that both `bandStructure` and `gfactorCalculation` executables are created successfully, then running basic calculations with provided example input files.

**Acceptance Scenarios**:

1. **Given** a clean repository with required dependencies installed, **When** running `make all`, **Then** both executables are created without errors
2. **Given** the created executables, **When** running `./bandStructure` with `bulk.example`, **Then** calculation completes successfully and produces output files
3. **Given** the created executables, **When** running `./gfactorCalculation` with `gfactor.example`, **Then** calculation completes successfully and produces output files

---

### User Story 2 - Result Verification Framework (Priority: P2)

A researcher needs to verify that bulk and quantum well calculations produce physically accurate results that match expected scientific behavior.

**Why this priority**: Scientific accuracy is non-negotiable for physics simulations - incorrect results could lead to invalid scientific conclusions.

**Independent Test**: Can be fully tested by running calculations with known material parameters and comparing results against analytical solutions or published benchmarks for simple cases (e.g., infinite quantum well).

**Acceptance Scenarios**:

1. **Given** bulk calculation input with InAs60Sb40 material parameters, **When** running the calculation, **Then** band structure results match expected energy levels within acceptable tolerance
2. **Given** quantum well calculation with GaSb/InAs/AlSb structure, **When** running the calculation, **Then** confinement effects produce quantized energy levels as expected
3. **Given** g-factor calculation input, **When** running the calculation, **Then** g-factor values fall within physically reasonable ranges

---

### User Story 3 - Standardized Plotting Scripts (Priority: P3)

A researcher needs reusable gnuplot scripts to visualize calculation results for verification and analysis purposes.

**Why this priority**: Visualization is essential for result verification and scientific communication, but depends on having working calculations first.

**Independent Test**: Can be fully tested by running calculations to generate data files, then using gnuplot scripts to create publication-quality plots that clearly show the expected physical behavior.

**Acceptance Scenarios**:

1. **Given** calculation output files, **When** running band structure plotting script, **Then** gnuplot generates clear band structure plots with proper labels and units
2. **Given** quantum well calculation results, **When** running quantum well plotting script, **Then** gnuplot generates plots showing wave functions and energy levels
3. **Given** g-factor calculation results, **When** running g-factor plotting script, **Then** gnuplot generates plots showing g-factor dependence on parameters

---

### Edge Cases

- What happens when build dependencies are missing or incompatible versions?
- How does the system handle malformed input configuration files (missing fields, wrong order, invalid labels)?
- What happens when calculation parameters result in numerical instabilities?
- How does the system handle insufficient memory for large quantum well calculations?
- What happens when gnuplot is not installed or has compatibility issues?
- What happens if output directory is not writable or absent?

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST compile successfully with any Fortran compiler and any BLAS/LAPACK implementation (documented example provided)
- **FR-002**: System MUST produce two working executables: `bandStructure` and `gfactorCalculation`
- **FR-003**: System MUST run bulk calculations using provided example input files
- **FR-004**: System MUST run quantum well calculations using provided example input files
- **FR-005**: System MUST validate calculation results against known physical behavior
- **FR-006**: System MUST provide clear error messages for common build and runtime issues
- **FR-007**: System MUST generate standardized gnuplot scripts for band structure visualization
- **FR-008**: System MUST generate standardized gnuplot scripts for quantum well visualization
- **FR-009**: System MUST generate standardized gnuplot scripts for g-factor visualization
- **FR-010**: System MUST include step-by-step build documentation with examples and dependency requirements
- **FR-011**: System MUST include example input files with documented parameter meanings
- **FR-012**: System MUST produce output files in plain text column format compatible with gnuplot scripts
 - **FR-013**: System MUST accept an input configuration filename as first CLI argument (defaulting to `input.cfg` when omitted) and MUST document this behavior consistently across docs and scripts
 - **FR-014**: System MUST write all generated output artifacts to `outputs/<run-id>/` (e.g., timestamped), never to repository root, and include minimal metadata (input filename, timestamp)

### Key Entities *(include if feature involves data)*

- **Build Configuration**: Makefile settings, compiler flags, library linking options, and dependency management
- **Input Parameters**: Material properties, geometric parameters, calculation settings, and boundary conditions
- **Calculation Results**: Energy eigenvalues, wave functions, band structures, and g-factor values
- **Visualization Scripts**: Gnuplot scripts with standardized formatting, labels, and output formats

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Build process completes successfully in under 5 minutes on standard development machine
- **SC-002**: Both bulk and quantum well calculations run without errors using provided examples
- **SC-003**: Calculation results match expected physical behavior within 5% tolerance for simple test cases
- **SC-004**: Gnuplot scripts generate publication-quality plots within 30 seconds of calculation completion
- **SC-005**: Build documentation enables new users to successfully compile and run calculations within 30 minutes
- **SC-006**: All example input files produce valid results that can be visualized with provided plotting scripts