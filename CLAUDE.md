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

**Executables** at `build/src/bandStructure` and `build/src/gfactorCalculation`.

**Prerequisites:** gfortran, Intel MKL (sequential, LP64 interface), FFTW3, CMake >= 3.15, Ninja (optional). MKL defaults: `MKL_INTERFACE=lp64`, `MKL_THREADING=sequential` (set in root `CMakeLists.txt`).

**Gotcha — stale `.mod` files:** Old `.mod` in project root shadow fresh ones in `build/`. Run `rm -f *.mod` if you get type mismatch errors.

**Gotcha — `cp -i` alias:** The shell has `cp -i` alias. Use the Write tool (not `cp`) when creating `input.cfg`.

## Testing

```bash
# Configure with tests enabled (pFUnit installed at $HOME/.local/pfunit/PFUNIT-<ver>)
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl \
    -DBUILD_TESTING=ON -DPFUNIT_DIR=$HOME/.local/pfunit/PFUNIT-4.16/cmake
cmake --build build

# Run tests
ctest --test-dir build                    # all 9 tests (5 unit + 4 regression)
ctest --test-dir build -L unit            # pFUnit unit tests only
ctest --test-dir build -L regression      # regression/golden-output tests only
ctest --test-dir build -V                 # verbose output
```

pFUnit installs as **uppercase `PFUNIT`** — use `-DPFUNIT_DIR=.../PFUNIT-<ver>/cmake`. Built from source: `git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit`.

Regression tests use shell scripts + Python (`tests/regression/compare_output.py`) comparing numerical output against reference data. No extra dependencies beyond Python 3.

## Running

```bash
./build/src/bandStructure       # Reads input.cfg, outputs to output/
./build/src/gfactorCalculation  # Reads input.cfg, outputs to output/ (requires k=0)
```

Write example config (`bulk.example`, `quantumwell.example`, or `gfactor.example`) contents to `input.cfg`.

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
  unit/       pFUnit .pf test files (test_defs, test_finitedifferences, test_utils, test_parameters, test_hamiltonian)
  integration/  shell scripts for full-executable tests
  regression/   configs/, data/, compare_output.py
cmake/        FindFFTW3.cmake
build/        .o, .mod, executables (created by cmake)
scripts/      gnuplot plotting scripts
```

### Two executables

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
  <- input_parser.f90         (shared input.cfg parser → simulation_config type)
```

### Key design concepts

- **Basis ordering** (fixed throughout): bands 1-4 = valence (HH, LH, LH, SO), bands 5-6 = split-off, bands 7-8 = conduction
- **Dual mode**: `confinement=0` → bulk (8x8), `confinement=1` → quantum well (8N x 8N, N = FDstep)
- **QW Hamiltonian**: 8x8 block matrix, each block NxN. Blocks: k.p terms (Q, R, S, T, P) + band offsets
- **`simulation_config`** derived type in `defs.f90` holds all parsed input parameters; `input_parser.f90` populates it
- **`confinementInitialization`** precomputes material parameters at each z-point into `kpterms(fdStep, fdStep, 10)`
- **Sparse vs dense**: g-factor uses MKL SpBLAS for large QW; band structure uses dense LAPACK (`zheevx`)
- **Foreman renormalization**: disabled by default (`renormalization = .False.` in defs.f90)
- **W-variant materials** (GaAsW, InAsW, etc.): use Winkler's parameter set with InSb as EV reference. Non-W materials use Vurgaftman parameters. EC = EV + Eg convention for all materials with EV defined

### Input file (`input.cfg`)

Label-value format with comments (`!`). Key fields: `waveVector` (kx/ky/kz/k0), `waveVectorMax`, `waveVectorStep`, `confinement` (0/1), `FDstep`, `FDorder` (default 2), `numLayers`, `materialN` (name + z-range), `numcb`/`numvb`, `ExternalField` + `EFParams`.

g-factor additional fields: `whichBand` (0=CB, 1=VB) and `bandIdx` (subband index, default 1).

## Code Conventions

- Precision kinds in `defs.f90`: `sp` (single), `dp` (double), `qp` (quad), `iknd` (int64)
- k-vectors and wavevector sweeps use `real(dp)` throughout
- File/function length guidelines from `.clinerules`: 300 lines/file, 50 lines/function
- Prioritize: Correctness > Maintainability > Readability > Type-Exactness > Performance > Minimal

## Git Workflow

- **Branch naming:** `feature/<description>`, `fix/<description>`
- **Commits:** Conventional Commits (`feat:`, `fix:`, `refactor:`, `test:`, `docs:`)
- **PRs:** All tests must pass before merge. Use `gh pr create` from feature branches.

## Boundaries

- **NEVER** modify material parameters in `parameters.f90` without verifying against published references (Vurgaftman 2001, Winkler 2003)
- **NEVER** change the basis ordering (bands 1-4 valence, 5-6 split-off, 7-8 conduction) — it is hardcoded throughout
- **NEVER** commit `input.cfg` with personal test configs — use `tests/regression/configs/` for test configs
- **Require approval** for: changes to `defs.f90` derived types, Hamiltonian construction in `hamiltonianConstructor.f90`, FD stencil coefficients in `finitedifferences.f90`
