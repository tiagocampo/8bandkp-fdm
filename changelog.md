# Changelog

All notable changes to this project will be documented in this file.

## [2024-12-28] - Continuous Integration Setup

### Added
- Implemented GitHub Actions workflow:
  - Added CI configuration in `.github/workflows/ci.yml`
  - Set up matrix builds for Release/Debug and MKL/non-MKL configurations
  - Added automated dependency installation
  - Configured test execution
  - Added compiler warning checks

### Changed
- Enhanced CMake configuration:
  - Added compile commands generation for static analysis
  - Enabled comprehensive compiler warnings
  - Added Intel compiler support with appropriate flags
  - Improved warning configuration for both Debug and Release builds

## [2024-12-28] - Directory Structure Reorganization

### Added
- Created modular directory structure:
  - `src/core/`: Core functionality (defs.f90, parameters.f90, utils.f90)
  - `src/math/`: Mathematical operations (mkl_spblas.f90, mkl_sparse_handle.f90, finitedifferences.f90)
  - `src/physics/`: Physics calculations (hamiltonianConstructor.f90, gfactor_functions.f90)
  - `src/io/`: Input/Output operations (outputFunctions.f90)
  - `src/apps/`: Main applications (main.f90, main_gfactor.f90)
  - `tests/`: Test directories (unit, integration, regression)
  - `docs/`: Documentation directories (api, user, examples)

### Changed
- Updated CMake configuration to reflect new directory structure:
  - Modified source file paths in src/CMakeLists.txt
  - Updated file dependencies to use new paths
  - Verified build system works with new structure

## [2024-12-28] - Build System Modernization

### Added
- Implemented CMake build system:
  - Created root CMakeLists.txt with project configuration and build options
  - Added cmake/FindMKL.cmake module for MKL dependency detection
  - Created src/CMakeLists.txt for source build configuration
  - Added proper dependency tracking between source files
  - Files affected: `CMakeLists.txt`, `cmake/FindMKL.cmake`, `src/CMakeLists.txt`

### Changed
- Improved build organization:
  - Separated common code into static library
  - Added proper compiler flags management
  - Improved MKL library detection and linking
  - Added installation targets
  - Added CPack configuration for packaging

## [2024-03-19]

### Changed
- Modified diagonalization method in `src/main.f90`:
  - Reverted from `zheev` back to `zheevx` to compute only specified eigenvalues
  - Added proper range selection for quantum well states (il to iuu)
  - Added checks to limit requested bands to available states (2*fdStep for CB, 6*fdStep for VB)
  - Improved error handling for diagonalization failures
  - Files affected: `src/main.f90`

- Updated band structure plotting in `scripts/plot_bands.gp`:
  - Simplified plotting script to show all energy bands
  - Removed manual energy range restrictions
  - Improved plot readability with consistent line styles
  - Files affected: `scripts/plot_bands.gp`

## [2024-12-28]

### Added
- Created visualization scripts in `scripts/` directory:
  - `plot_bands.gp`: Band structure visualization
  - `plot_eigenfunctions.gp`: Eigenfunction components visualization
  - `plot_potential.gp`: Potential profile visualization
  - `plot_gfactor.gp`: g-factor components visualization
  - `plot_all.sh`: Shell script to run all plots
  - Files affected: New files in `scripts/` directory

- Created `.gitignore` file:
  - Added rules for build artifacts (*.o, *.mod)
  - Added rules for executables (bandStructure, gfactorCalculation)
  - Added rules for data files (*.dat, fort.*, *.txt)
  - Added rules for IDE and system-specific files
  - Files affected: `.gitignore`

### Changed
- Implemented output directory organization:
  - Added 'output' directory for all generated files
  - Added automatic directory creation
  - Moved all data files to output directory
  - Improved file unit handling
  - Files affected: `src/outputFunctions.f90`, `src/main.f90`, `src/main_gfactor.f90`

- Modified Makefile to improve build organization and dependency management:
  - Added `build` directory for all compilation outputs
  - Moved `.o` (object) files to `build` directory
  - Moved `.mod` (module) files to `build` directory using gfortran's `-J` flag
  - Removed duplicate `SRC_D` variable definition
  - Updated dependencies to properly track module files
  - Switched from threaded MKL to sequential MKL for better compatibility
  - Added pattern rules for more maintainable compilation
  - Files affected: `Makefile`

### Fixed
- Fixed missing variable declarations:
  - Added `iounit` variable declaration in main programs
  - Fixed implicit typing errors
  - Files affected: `src/main.f90`, `src/main_gfactor.f90`

- Fixed bulk calculations in `main.f90`:
  - Added explicit temporary array (HTmp) for bulk diagonalization
  - Switched from zheevx to zheev for bulk calculations
  - Fixed array dimensions to match 8x8 bulk Hamiltonian
  - Eliminated runtime warnings from array sections
  - Files affected: `src/main.f90`

- Fixed eigenfunction output in `outputFunctions.f90`:
  - Added separate handling for bulk and quantum well cases
  - Preserved complex components for bulk eigenvectors
  - Fixed output format for complex components
  - Improved readability with one complex number per line
  - Added proper implementation of `get_unit` subroutine
  - Removed dependency on external `get_lun` function
  - Files affected: `src/outputFunctions.f90`

- Fixed eigenfunction writing in `main_gfactor.f90`:
  - Added missing `is_bulk` parameter
  - Set to `.false.` for quantum well calculations
  - Files affected: `src/main_gfactor.f90`

- Resolved out-of-bounds array access in `src/outputFunctions.f90` when calculating probability densities for quantum well structures.
  - Modified the loop structure in `writeEigenfunctions` to correctly iterate over eigenstates and bands.
  - Added a helper subroutine `get_eigenvector_component` to extract spatial components for each band.
  - Files affected: `src/outputFunctions.f90`

- Adjusted the calculation of `numcb` and `numvb` for quantum well calculations in `src/main.f90`. These values are now set based on `fdStep` to ensure the `eigv` array is allocated with sufficient size, resolving the array bound mismatch error.
  - Files affected: `src/main.f90`

- Corrected the argument passed to the `ZB8bandBulk` subroutine in `src/main.f90`. The entire `params` array is now passed instead of a scalar element, resolving the rank mismatch error.
  - Files affected: `src/main.f90` 