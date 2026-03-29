# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Fortran 90 code solving the **8-band zinc-blende k.p Hamiltonian** via finite differences. Computes electronic band structures for bulk semiconductors and quantum wells, plus Landau g-factors via second-order Lowdin partitioning. GPL v3.0, authored by Tiago de Campos.

## Build Commands

```bash
# Configure (first time or after clean)
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl

# Build
cmake --build build                # builds both executables into build/src/

# Or use the Make wrapper
make all            # Configure + build both executables
make run            # Build and run bandStructure
make clean          # Remove build/ directory
make clean_all      # Also remove output files
make test           # Configure with BUILD_TESTING=ON, build, run ctest
```

**Executables** are at `build/src/bandStructure` and `build/src/gfactorCalculation`.

**Prerequisites:** gfortran, Intel MKL (sequential, LP64 interface), FFTW3, CMake >= 3.15, Ninja (optional). The MKL defaults are set in `CMakeLists.txt` (`MKL_INTERFACE=lp64`, `MKL_THREADING=sequential`).

**Stale `.mod` files:** Old `.mod` files in the project root can shadow fresh ones in `build/`. If you get inexplicable type mismatch errors, run `rm -f *.mod` before building.

## Testing

```bash
# Configure with tests enabled
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl \
    -DBUILD_TESTING=ON -DpFUnit_DIR=/path/to/pfunit/install
cmake --build build

# Run tests
ctest --test-dir build                    # all tests
ctest --test-dir build -L unit            # pFUnit unit tests only
ctest --test-dir build -L regression      # regression/golden-output tests only
ctest --test-dir build -V                 # verbose output
```

**Prerequisites for unit tests:** pFUnit >= 4.9 (built from source, provides `pFUnitConfig.cmake`). Install via `git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit && cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local/pfunit .. && make install`.

**Regression tests** use shell scripts + Python (`tests/regression/compare_output.py`) to compare numerical output against reference data. No extra dependencies beyond Python 3.

## Running

```bash
./build/src/bandStructure       # Reads input.cfg, outputs to output/
./build/src/gfactorCalculation  # Reads input.cfg, outputs to output/ (requires k=0)
```

Both programs read from `input.cfg` in the working directory. Output goes to `output/` directory. Plot results with scripts in `scripts/`.

Copy an example config to get started: write `bulk.example` contents to `input.cfg` (or `quantumwell.example` / `gfactor.example`). Use the Write tool directly rather than `cp` to avoid interactive alias issues.

## Architecture

### Directory structure

```
src/
  core/       defs.f90, parameters.f90, utils.f90
  math/       mkl_spblas.f90, mkl_sparse_handle.f90, finitedifferences.f90
  io/         outputFunctions.f90, input_parser.f90
  physics/    hamiltonianConstructor.f90, gfactor_functions.f90
  apps/       main.f90, main_gfactor.f90
tests/
  unit/       pFUnit .pf test files (test_defs.pf, test_fd.pf, etc.)
  integration/  shell scripts for full-executable tests
  regression/   configs/, data/, compare_output.py
cmake/        FindFFTW3.cmake
build/        .o, .mod, executables (created by cmake)
scripts/      gnuplot plotting scripts
```

### Two executables from one source tree

| Executable | Entry Point | Purpose |
|---|---|---|
| `bandStructure` | `src/apps/main.f90` (`program kpfdm`) | Band structure vs. k-vector sweep |
| `gfactorCalculation` | `src/apps/main_gfactor.f90` (`program gfactor`) | g-factor at Gamma point (k=0) |

### Module dependency graph

```
defs.f90                      (kinds, constants, derived types — no deps)
  <- parameters.f90           (material parameter database for 25+ semiconductors)
  <- mkl_spblas.f90           (vendor: Intel MKL sparse BLAS interface)
       <- utils.f90           (dense-to-sparse conversion, Simpson integration)
            <- finitedifferences.f90  (FD stencils/orders 2-10, Toeplitz matrices)
                 <- hamiltonianConstructor.f90  (8x8 bulk & 8NxN QW Hamiltonian)
                      <- gfactor_functions.f90  (Lowdin partitioning, spin/momentum matrix elements)
  <- outputFunctions.f90      (eigenvalue/eigenfunction file I/O)
```

### Key design concepts

- **Basis ordering** (fixed throughout code): bands 1-4 = valence (HH, LH, LH, SO), bands 5-6 = split-off, bands 7-8 = conduction
- **Dual mode**: `confinement=0` gives bulk (8x8 Hamiltonian), `confinement=1` gives quantum well (8N x 8N, where N = FDstep grid points)
- **QW Hamiltonian**: 8x8 block matrix, each block is NxN. Blocks correspond to k.p terms (Q, R, S, T, P, band offsets)
- **Spatial profile**: `confinementInitialization` precomputes material parameters at each z-point into `kpterms(fdStep, fdStep, 10)`
- **Sparse path**: g-factor code optionally uses MKL SpBLAS for large QW systems; band structure code uses dense LAPACK (`zheevx`)
- **Foreman renormalization**: disabled by default (`renormalization = .False.` in defs.f90)

### Input file (`input.cfg`)

Label-value format with comments (`!`). Key fields: `waveVector` (direction: kx/ky/kz/k0), `waveVectorMax` (k-range), `waveVectorStep` (total k-points in sweep), `confinement` (0/1), `FDstep` (grid points), `numLayers`, `materialN` (name + z-range), `numcb`/`numvb` (bands to compute), `ExternalField` + `EFParams`.

## Code Conventions

- Precision kinds in `defs.f90`: `sp` (single), `dp` (double), `qp` (quad), `iknd` (int64)
- k-vectors and wavevector sweeps use `real(dp)` throughout
- File/function length guidelines from `.clinerules`: 300 lines/file, 50 lines/function
- Prioritize: Correctness > Maintainability > Readability > Type-Exactness > Performance > Minimal
