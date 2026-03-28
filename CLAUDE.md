# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Fortran 90 code solving the **8-band zinc-blende k.p Hamiltonian** via finite differences. Computes electronic band structures for bulk semiconductors and quantum wells, plus Landau g-factors via second-order Lowdin partitioning. GPL v3.0, authored by Tiago de Campos.

## Build Commands

```bash
make all            # Build both executables: bandStructure and gfactorCalculation
make run            # Build and run bandStructure
make clean          # Remove build/ directory and data files
make clean_all      # Also remove executables
```

**Prerequisites:** gfortran, Intel MKL (sequential), FFTW3. The Makefile LDFLAGS must be set for your linear algebra backend.

**No test framework exists.** Verification is done by running executables against known results.

**Stale `.mod` files:** Old `.mod` files in the project root can shadow fresh ones in `build/`. If you get inexplicable type mismatch errors, run `rm -f *.mod` before building.

## Running

```bash
./bandStructure       # Reads input.cfg, outputs to output/
./gfactorCalculation  # Reads input.cfg, outputs to output/ (requires k=0)
```

Both programs read from `input.cfg` in the working directory. Output goes to `output/` directory. Plot results with scripts in `scripts/`.

Copy an example config to get started: `cp bulk.example input.cfg` (or `quantumwell.example` / `gfactor.example`).

## Architecture

### Directory structure

```
src/
  core/       defs.f90, parameters.f90, utils.f90
  math/       mkl_spblas.f90, mkl_sparse_handle.f90, finitedifferences.f90
  io/         outputFunctions.f90
  physics/    hamiltonianConstructor.f90, gfactor_functions.f90
  apps/       main.f90, main_gfactor.f90
build/        .o and .mod files (created by make)
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
