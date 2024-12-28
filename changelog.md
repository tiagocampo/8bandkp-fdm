# Changelog

All notable changes to this project will be documented in this file.

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