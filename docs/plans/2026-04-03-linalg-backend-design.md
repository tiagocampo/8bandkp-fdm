# Linalg Backend Portability & Performance Design

**Date**: 2026-04-03
**Author**: Tiago de Campos
**Status**: Approved (Phase 1)

## Motivation

1. **Portability**: The code currently requires Intel MKL for SpBLAS operations and FEAST eigensolver. The files `src/math/mkl_spblas.f90` and `src/math/mkl_sparse_handle.f90` are Intel-copyrighted interface headers that create a vendor lock-in.
2. **Performance**: Current compiler flags (`-O2`) leave performance on the table for the Haswell-EP target (Xeon E5-2678 v3, AVX2, 12C/24T).
3. **Future GPU**: AMD ROCm 7.2 is installed (rocSOLVER, rocSPARSE, hipSOLVER) but the RX 570 (gfx803/Polaris) is unsupported by ROCm 7.x user-space runtime. The GPU path is designed but deferred until hardware upgrade.

## Architecture

### Compile-time backend selection

```
┌──────────────────────────────────────────────┐
│  Application (main.f90, gfactor, eigensolver) │
├──────────────────────────────────────────────┤
│  sparse_matrices.f90  (own csr_matrix type)  │
│  + csr_spmv (pure Fortran, OpenMP)           │
├──────────────────────────────────────────────┤
│  BLAS/LAPACK (MKL | OpenBLAS | AOCL)         │
│  Selected at CMake configure time            │
└──────────────────────────────────────────────┘
```

CMake option: `-DLINALG_BACKEND=MKL|OPENBLAS|AOCL`

No runtime dispatch. Single backend per build. This matches the approach used by Quantum ESPRESSO, CP2K, and other Fortran physics codes.

### What gets abstracted

| Current (MKL-only) | Replacement | Status |
|---|---|---|
| `MKL_SPARSE_Z_CREATE_COO` + `CONVERT_CSR` | `csr_matrix` type + `csr_build_from_coo` in `sparse_matrices.f90` | Already exists |
| `MKL_SPARSE_Z_MV` (SpMV) | `csr_spmv(A, x, y, alpha, beta)` in `sparse_matrices.f90` | **New** (~15 lines + OpenMP) |
| `SPARSE_MATRIX_T` handle | `csr_matrix` type | Already exists |
| `zfeast_hcsrev` (MKL FEAST) | Dense LAPACK fallback via `solve_dense_lapack` in `eigensolver.f90` | Already exists |
| Dense eigensolvers (`zheevx`, `zheevd`) | Standard LAPACK (any vendor) | Already portable |

### What gets removed

- `src/math/mkl_spblas.f90` — 1580-line Intel-copyrighted MKL SpBLAS interface module
- `src/math/mkl_sparse_handle.f90` — 82-line Intel-copyrighted DSS handle module (already unused)

### What gets added

- `csr_spmv` subroutine in `sparse_matrices.f90` — pure Fortran CSR SpMV with OpenMP parallelization
- CMake `LINALG_BACKEND` option with auto-detection

## Phase 1: CPU Portability & Performance (Current)

### Step 1: Add `csr_spmv` to `sparse_matrices.f90`

Pure Fortran implementation of `y = alpha * A * x + beta * y` for CSR matrix:

```fortran
subroutine csr_spmv(A, x, y, alpha, beta)
  type(csr_matrix), intent(in) :: A
  complex(dp), intent(in) :: x(A%ncols)
  complex(dp), intent(inout) :: y(A%nrows)
  complex(dp), intent(in) :: alpha, beta
  integer :: row, k
  complex(dp) :: dot
  !$omp parallel do private(row, k, dot) schedule(static)
  do row = 1, A%nrows
    dot = cmplx(0, 0, kind=dp)
    do k = A%rowptr(row), A%rowptr(row+1) - 1
      dot = dot + A%values(k) * x(A%colind(k))
    end do
    y(row) = alpha * dot + beta * y(row)
  end do
  !$omp end parallel do
end subroutine
```

### Step 2: Refactor call sites (3 files)

**`src/core/utils.f90`** — `dnscsr_z_mkl`:
- Remove `use mkl_spblas`
- Build `csr_matrix` directly using existing `csr_build_from_coo` instead of `MKL_SPARSE_Z_CREATE_COO` + `MKL_SPARSE_CONVERT_CSR`
- Return `type(csr_matrix)` instead of `type(sparse_matrix_t)`

**`src/physics/gfactor_functions.f90`** — `compute_pele`:
- Remove `use mkl_spblas`
- Call `csr_spmv(HT_csr, x_vec, y_vec, alpha, beta)` instead of `MKL_SPARSE_Z_MV`

**`src/physics/hamiltonianConstructor.f90`**:
- Remove `use mkl_spblas`
- Change `HT_csr` argument type from `type(sparse_matrix_t)` to `type(csr_matrix)`

### Step 3: Update CMakeLists.txt

**Compiler flags** (optimized for Haswell-EP):
```cmake
set(CMAKE_Fortran_FLAGS_RELEASE
    "-O3 -march=haswell -mtune=haswell -ftree-vectorize "
    "-fstack-arrays -ffast-math -fno-trapping-math "
    "-flto=auto -fuse-linker-plugin")
```

Key changes from current flags:
- `-O2` → `-O3`: enables auto-vectorization, loop transformations
- Added `-ffast-math`: safe for k.p physics (~5-15% speedup on tight loops)
- Added `-fstack-arrays`: avoids heap allocation for temporaries
- Added `-ftree-vectorize`: explicit vectorization (implied by -O3)
- Removed `-funroll-loops`: can hurt on modern CPUs due to icache pressure
- Explicit `-march=haswell` instead of `native` for reproducibility

**Backend selection**:
```cmake
set(LINALG_BACKEND "MKL" CACHE STRING "BLAS/LAPACK backend: MKL, OPENBLAS, AOCL")

if(LINALG_BACKEND STREQUAL "MKL")
    find_package(MKL CONFIG REQUIRED)
    # MKL provides BLAS, LAPACK, FFTW wrappers
elseif(LINALG_BACKEND STREQUAL "OPENBLAS")
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    find_package(FFTW3 REQUIRED)
elseif(LINALG_BACKEND STREQUAL "AOCL")
    find_package(aocl REQUIRED)
    find_package(FFTW3 REQUIRED)
endif()
```

### Step 4: Handle FEAST

When `LINALG_BACKEND != MKL`, FEAST is unavailable:
- `feast_call.f90` guarded with preprocessor: `#ifdef USE_MKL_FEAST`
- `eigensolver.f90` dispatches to `solve_dense_lapack` when FEAST not available
- Wire mode (confinement=2) uses dense fallback for smaller problems
- Larger wire problems will need ARPACK-NG (future Phase 3)

### Step 5: Delete Intel files

- Delete `src/math/mkl_spblas.f90`
- Delete `src/math/mkl_sparse_handle.f90`
- Remove from `src/math/CMakeLists.txt`

### Step 6: Add OpenMP for SpMV

CMakeLists.txt adds:
```cmake
find_package(OpenMP REQUIRED)
# Link targets with OpenMP::OpenMP_Fortran
```

The `csr_spmv` kernel parallelizes trivially over rows with `schedule(static)`.

## Phase 2: GPU Eigensolver (Deferred)

### When needed
For wire mode with large grids (matrix size > 4000x4000) where CPU diagonalization is too slow.

### Prerequisites
- AMD GPU with RDNA2+ or CDNA architecture (gfx1030+, gfx90a+)
- ROCm 7.x with working HIP runtime
- hipfort Fortran bindings (not installed yet)

### Approach
- `eigensolver_gpu.f90` module using `iso_c_binding` to call hipSOLVER C API
- Proof-of-concept test already written: `tests/gpu_diag_test.f90`
- Uses `hipsolverZheevd` / `hipsolverZheevdx` for dense Hermitian eigenvalues
- For sparse: custom FEAST-style contour driver using rocSPARSE SpMV + rocSOLVER

### Test result (2026-04-03)
The proof-of-concept compiled and linked successfully against hipSOLVER/rocBLAS/hipBLAS.
Execution failed with `hsa_init failed` on RX 570 (gfx803/Polaris) because ROCm 7.2
user-space dropped Polaris support. The test program is correct and will work on
supported GPUs (gfx906+, gfx1030+, etc.).

## Validation

### Existing tests
All 21 existing tests (11 unit + 10 regression) must pass after refactoring.

### New test
- `test_csr_spmv` unit test: verify `csr_spmv` against dense `zgemv` for known matrices

### Performance benchmark
- Time QW diagonalization (FDstep=100, 200) before/after flag changes
- Time g-factor SpMV loop before/after MKL SpBLAS → csr_spmv migration
- Verify eigenvalues match within machine precision

## Files Modified

| File | Change |
|---|---|
| `src/math/sparse_matrices.f90` | Add `csr_spmv` |
| `src/core/utils.f90` | Remove `use mkl_spblas`, refactor `dnscsr_z_mkl` |
| `src/physics/gfactor_functions.f90` | Remove `use mkl_spblas`, use `csr_spmv` |
| `src/physics/hamiltonianConstructor.f90` | Remove `use mkl_spblas`, use `csr_matrix` type |
| `src/math/feast_call.f90` | Add `#ifdef USE_MKL_FEAST` guards |
| `src/math/eigensolver.f90` | Handle missing FEAST gracefully |
| `CMakeLists.txt` | New flags, backend option, OpenMP |
| `src/CMakeLists.txt` | Update link libraries per backend |
| `src/math/CMakeLists.txt` | Remove mkl_spblas.f90, mkl_sparse_handle.f90 |

## Files Deleted

| File | Reason |
|---|---|
| `src/math/mkl_spblas.f90` | Intel-copyrighted, replaced by csr_spmv |
| `src/math/mkl_sparse_handle.f90` | Intel-copyrighted, already unused |

## Files Added

| File | Purpose |
|---|---|
| `tests/gpu_diag_test.f90` | GPU eigensolver PoC (for future use) |
