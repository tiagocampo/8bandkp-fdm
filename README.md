# 8bandkp-fdm
Fortran implementation of 8band ZB k·p-method using finite difference method, working for bulk and quantum-wells

## Disclaim:
 * This software is for my personal usage only. It comes with no warranty at all.
Author: Tiago de Campos

## Description:
 * Solves k·p Hamiltonian by finite differences method for zinc-blende semiconductors;
 * Implements a full 8-band model (6 valence bands + 2 conduction bands) using modern Fortran;
 * Based on Chuang, S. L. and Chang, C. S., Semiconductor Science and Technology, 12, 252, 1997;
 * For a description of the derivatives see https://www.wias-berlin.de/people/john/LEHRE/NUM_PDE_FUB/num_pde_fub_2.pdf;
 * Uses renormalization of the Kane interband momentum matrix element, P, suggested by Foreman B A 1997 Physical Review B 56 R12748;
 * g-factor calculation using second order Lowdin partitioning, following:
   - Winkler, R. Spin-orbit coupling effects in two-dimensional electron and hole systems; Physics and Astronomy online Library 191; Springer, 2003
   - Tadjine, A; Niquet, Y.-M.; Delerue, C. Universal behavior of electron g-factors in semiconductor nanostructures. Physical Review B 2017, 95, 235437
 * Supports external electric field (non-self-consistent calculation)
 * Parallel processing support through OpenMP
 * High-performance linear algebra operations using LAPACK/MKL
 * Selective eigenvalue computation for improved performance
 * Automatic band structure visualization with customizable plotting

## Project Structure

The project is organized into the following directories:

### Source Code (`src/`)
* `core/`: Core functionality and definitions
  - `defs.f90`: Type definitions and constants
  - `parameters.f90`: Material parameters and configuration
  - `utils.f90`: Utility functions and common operations

* `math/`: Mathematical operations and numerical methods
  - `mkl_spblas.f90`: MKL sparse BLAS interface
  - `mkl_sparse_handle.f90`: MKL sparse matrix handling
  - `finitedifferences.f90`: Finite difference method implementation

* `physics/`: Physics calculations and models
  - `hamiltonianConstructor.f90`: k·p Hamiltonian construction
  - `gfactor_functions.f90`: g-factor calculations

* `io/`: Input/Output operations
  - `outputFunctions.f90`: File I/O and data output

* `apps/`: Main applications
  - `main.f90`: Band structure calculation program
  - `main_gfactor.f90`: g-factor calculation program

### Tests (`tests/`)
* `unit/`: Unit tests for individual components
* `integration/`: Integration tests for combined functionality
* `regression/`: Regression tests for verified results

### Documentation (`docs/`)
* `api/`: API documentation and code reference
* `user/`: User guides and tutorials
* `examples/`: Example configurations and use cases

### Build System
* `CMakeLists.txt`: Main CMake configuration
* `src/CMakeLists.txt`: Source build configuration
* `cmake/`: CMake modules and find scripts

### Other
* `scripts/`: Utility and plotting scripts
* `output/`: Generated data and results
* `build/`: Build directory (created during compilation)

## Key Features:
 * Full 8-band k·p Hamiltonian implementation
 * Support for both bulk and quantum well calculations
 * Finite difference method for spatial discretization
 * Multiple material layer support for heterostructures
 * Band structure calculations along arbitrary k-vector directions
 * g-factor calculations for spin-related properties
 * External electric field effects
 * Efficient sparse matrix techniques
 * OpenMP parallelization for improved performance
 * Selective eigenvalue computation:
   - Bulk: Up to 8 bands (2 CB + 6 VB)
   - Quantum wells: Up to 2×fdStep CB and 6×fdStep VB states
 * Automated visualization tools for band structure and wavefunctions

## Use:

### General info

 * makefile provided with MKL and standard LAPACK/BLAS options
 * example input.cfg provided
 * Built-in material parameters for common semiconductors (check parameters.f90)
 * Additional material parameters can be found in Vurgaftman, I. Meyer, J. R. and Ram-Mohan, L. R., Journal of Applied Physics, 11, 5815, 2001

### Compilation steps

You can build the project using either Make or CMake:

#### Using Make (Legacy method)
 * modify Makefile to reflect either intel mkl library or standard lapack and blas by selecting appropriate LDFLAGS
 * make all

#### Using CMake (Recommended)
```bash
# Create and enter build directory
mkdir build && cd build

# Configure with CMake
cmake ..                    # Default Release build
# OR
cmake -DCMAKE_BUILD_TYPE=Debug ..  # Debug build with extra checks

# Build the project
make                       # Or 'make -j$(nproc)' for parallel build

# Optional: Install the executables
sudo make install         # Installs to default system location
# OR
cmake -DCMAKE_INSTALL_PREFIX=/custom/path .. && make install  # Custom install location
```

CMake build options:
* BUILD_TESTING: Enable/disable tests (default: OFF)
* ENABLE_OPENMP: Enable/disable OpenMP support (default: ON)
* USE_MKL: Use Intel MKL instead of standard LAPACK/BLAS (default: ON)

Example with all options:
```bash
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DENABLE_OPENMP=ON -DUSE_MKL=ON ..
```

This will generate two executable files:
 * bandStructure: for electronic band structure calculations
 * gfactorCalculation: for g-factor computations

### Input file structure

 * waveVector: reciprocal space direction (kx, ky, or kz). For confined systems, use kx or ky (kz is the confined direction)
 * waveVectorMax: percentage of Brillouin Zone to compute (keep close to zone center for k·p validity)
 * waveVectorStep: number of k-points between k=0 and k=kmax
 * confinement: 0 for bulk, 1 for quantum well
 * FDStep: number of spatial discretization points
 * numLayers: number of layers in the quantum well structure (centered at 0)
 * material definitions:
    * material1: host material (defines outer regions)
    * material2...N: well/barrier materials with positions and band offsets
 * numcb: number of conduction bands to compute
    * For bulk: maximum 2
    * For quantum wells: maximum 2 × fdstep
 * numvb: number of valence bands to compute
    * For bulk: maximum 6
    * For quantum wells: maximum 6 × fdstep
 * ExternalField: 0/1 and type (EF for electric field)
 * EFParams: field strength parameter

### Output files

The program generates several output files depending on the calculation type:

#### Band Structure Calculation (bandStructure):
1. `eigenvalues.dat`: Contains the energy eigenvalues for each k-point
   * Format: k-vector value followed by energy values
   * For bulk: outputs 8 bands
   * For quantum wells: outputs all requested bands (numcb + numvb)

2. `eigenfunctions_k_#####_ev_#####.dat`: Wavefunctions at k=0
   * Generated for each eigenstate
   * For bulk: only the 8 main bands
   * For quantum wells: includes position (z) and wavefunction components for all bands
   * Format: z-position followed by wavefunction components for each band

3. `parts.dat`: Band character analysis
   * For bulk: direct eigenvector components
   * For quantum wells: integrated probability densities for each band component
   * Format: 8 columns representing contribution from each band

4. `fort.101`: Potential profile (quantum wells only)
   * Contains the potential profile along the growth direction
   * Format: z-position and potential values for different bands

#### G-factor Calculation (gfactorCalculation):
* Similar output structure to band structure calculation
* Additional g-factor specific results for magnetic field effects
* Includes Zeeman splitting information

### Pre-requisits

 * Fortran compiler (gfortran recommended)
 * BLAS and LAPACK libraries
 * FFTW3 library
 * OpenMP support (included in most compilers)

#### Ubuntu installation:
```bash
sudo apt install gfortran gcc g++ liblapack-dev libblas-dev libfftw3-dev
```

For optimal performance, Intel MKL is recommended (see Intel's website for installation)

# Citation

This code has been used in several published works. If you use it, please cite:

Primary citation: 10.1021/acsaelm.0c00269

Additional relevant works:
 * 10.1063/1.5096970 
 * 10.1088/1361-648X/ab38a1
 * 10.1103/PhysRevB.97.245402
 * 10.1088/0268-1242/31/10/105002
 * 10.1103/PhysRevB.93.235204
 * 10.1063/1.4901209

# LICENSE
GNU General Public License v3.0
Please cite 10.1021/acsaelm.0c00269 or one of the other listed articles if it better suits your work.



