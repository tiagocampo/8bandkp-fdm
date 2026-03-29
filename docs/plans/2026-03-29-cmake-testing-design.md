# CMake Migration + Test Suite Design

**Date:** 2026-03-29
**Status:** Implemented

## Context

The project used a 94-line Makefile with manually tracked Fortran module dependencies. There was no test framework. We needed regression protection for future implementations.

## Decisions

### Build system: CMake + Ninja

- CMake auto-detects Fortran module dependencies (no manual tracking)
- Ninja generator for fast parallel builds
- `MKLConfig.cmake` from Intel oneAPI for MKL finding (LP64 + sequential defaults)
- Custom `cmake/FindFFTW3.cmake` for FFTW3
- Minimal Make wrapper preserving `make all/run/clean/test` commands

### Unit tests: pFUnit

- pFUnit is the standard Fortran unit testing framework
- Tests use `.pf` preprocessor files with `@test` and `@assertEqual` directives
- Integrated via `add_pfunit_ctest()` CMake macro
- 5 test suites: defs, finitedifferences, utils, parameters, hamiltonian

### Regression tests: shell + Python

- Shell scripts run executables in temp directories with known configs
- `compare_output.py` compares numerical output against reference data (tolerance-aware)
- 4 regression tests: bulk GaAs, bulk InAs, QW AlSbW/GaSbW/InAsW, g-factor CB
- Reference data generated from current (verified) code

## File Inventory

| File | Action |
|---|---|
| `CMakeLists.txt` | Created (root, with MKL LP64/sequential defaults) |
| `src/CMakeLists.txt` | Rewritten (added input_parser.f90, removed manual deps) |
| `cmake/FindFFTW3.cmake` | Created |
| `Makefile` | Replaced with CMake wrapper |
| `tests/CMakeLists.txt` | Created |
| `tests/unit/test_defs.pf` | Created |
| `tests/unit/test_finitedifferences.pf` | Created |
| `tests/unit/test_utils.pf` | Created |
| `tests/unit/test_parameters.pf` | Created |
| `tests/unit/test_hamiltonian.pf` | Created |
| `tests/integration/test_bulk_bandstructure.sh` | Created |
| `tests/integration/test_qw_bandstructure.sh` | Created |
| `tests/integration/test_gfactor.sh` | Created |
| `tests/regression/compare_output.py` | Created |
| `tests/regression/configs/*.cfg` | Created (4 configs) |
| `tests/regression/data/*/` | Created (4 reference datasets) |

## Verification

1. `cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl` — configures
2. `cmake --build build` — builds both executables
3. `./build/src/bandStructure` with bulk GaAs — produces correct eigenvalues
4. `./build/src/gfactorCalculation` with gfactor config — produces correct g-values
5. `make all/run/clean` — wrapper commands work
6. pFUnit tests compile and pass (requires pFUnit installation)
7. Regression tests pass with `ctest -L regression`
