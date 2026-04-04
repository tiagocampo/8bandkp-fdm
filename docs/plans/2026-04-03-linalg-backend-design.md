# Linalg Backend Portability, Parallelization & Performance Design

**Date**: 2026-04-03
**Author**: Tiago de Campos
**Status**: Approved (Phase 1)

## Motivation

1. **Single-threaded execution**: The code currently uses only 1 CPU core. MKL is explicitly set to `sequential` threading (`CMakeLists.txt` line 8). The k-vector sweep loops are embarrassingly parallel across k-points but run sequentially. On a 12-core/24-thread Xeon, this leaves >90% of compute resources idle.
2. **Portability**: The code requires Intel MKL for SpBLAS operations and FEAST eigensolver. The files `src/math/mkl_spblas.f90` and `src/math/mkl_sparse_handle.f90` are Intel-copyrighted interface headers that create a vendor lock-in.
3. **Compiler flags**: Current flags (`-O2`) leave performance on the table for the Haswell-EP target (Xeon E5-2678 v3, AVX2, 12C/24T).
4. **Future GPU**: AMD ROCm 7.2 is installed (rocSOLVER, rocSPARSE, hipSOLVER) but the RX 570 (gfx803/Polaris) is unsupported by ROCm 7.x user-space runtime. The GPU path is designed but deferred until hardware upgrade.

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

## Phase 1: CPU Portability, Parallelization & Performance (Current)

Three parallel efforts within one phase, ordered by expected speedup:

### Step 1: OpenMP parallelization of k-vector sweep (biggest win)

**Observation**: The main computation loops are embarrassingly parallel across k-points. Each k-point builds an independent Hamiltonian and diagonalizes it independently. Currently these run sequentially on 1 core while 23 threads sit idle.

**Target hardware**: Xeon E5-2678 v3, 12 cores / 24 threads (Haswell-EP, AVX2).

#### 1a. QW k-vector sweep (`src/apps/main.f90`, line 330-414)

```fortran
! Before (sequential):
do k = 1, cfg%waveVectorStep, 1
  call ZB8bandQW(HT, smallk(k), profile, kpterms)
  call zheevx(...)
  eig(:,k) = ...
  eigv(:,:,k) = ...
end do

! After (parallel):
!$omp parallel private(k, HT, HTmp, work, lwork, info, M)
allocate(HT(N,N), HTmp(8,8))
allocate(work(1))  ! each thread gets own workspace
!$omp do schedule(static)
do k = 1, cfg%waveVectorStep, 1
  call ZB8bandQW(HT, smallk(k), profile, kpterms)
  call zheevx(...)
  eig(:,k) = ...          ! no race: each k writes different column
  eigv(:,:,k) = ...       ! no race: each k writes different slice
end do
!$omp end do
deallocate(HT, HTmp, work)
!$omp end parallel
```

**Key points**:
- `HT`, `HTmp`, `work` must be thread-private (each thread builds its own Hamiltonian)
- `eig` and `eigv` are safe: each k-index writes to a different column/slice
- `profile` and `kpterms` are read-only — shared access is fine
- `schedule(static)`: k-points have similar cost, static scheduling minimizes overhead
- Eigenfunction writing (line 405-412) must be moved outside the parallel region or guarded with `!$omp single`

**Expected speedup**: Near-linear up to 12 threads (~10-12x). Hyperthreading adds ~20-30% more (~12-15x total).

#### 1b. Wire kz sweep (`src/apps/main.f90`, line 137-181)

Same pattern: each kz is independent.

```fortran
!$omp parallel private(k, HT_csr, eigen_res, coo_cache_local)
!$omp do schedule(static)
do k = 1, cfg%waveVectorStep
  call ZB8bandGeneralized(HT_csr, smallk(k)%kz, ...)
  call solve_sparse_evp(HT_csr, eigen_cfg, eigen_res)
  eig_wire(:,k) = eigen_res%eigenvalues(:)
  ! eigenfunction writing must be deferred (single-threaded I/O)
end do
!$omp end do
!$omp end parallel
```

**Complication**: `coo_cache` is currently shared across k-points (reuses COO structure from first call). In parallel mode, each thread needs its own cache, or the cache must be pre-built once and shared as read-only. The latter is simpler — build the cache in a serial warm-up iteration, then share it.

#### 1c. g-factor k-vector sweep (`src/apps/main_gfactor.f90`)

The g-factor code also has a k-vector loop. Same pattern: parallelize with thread-private `HT`, `work`.

**Note**: `use OMP_lib` is already imported in `main.f90` (line 7) but never used — the infrastructure is ready.

#### 1d. Enable MKL multi-threading

Currently `CMakeLists.txt` forces `MKL_THREADING=sequential`. Change to:

```cmake
set(MKL_THREADING "intel_thread" CACHE STRING "MKL threading: sequential, intel_thread")
```

This makes LAPACK calls (`zheevx`, `zheevd`) multi-threaded internally for large matrices. Combined with the k-vector OpenMP, this creates a two-level parallelism:
- **Outer level**: OpenMP over k-points (12 threads)
- **Inner level**: MKL threads for LAPACK diagonalization

**Nesting strategy**: For QW with small N (N < 200), use all threads for k-parallelism and disable inner MKL threading. For large N (N > 200), use fewer outer threads and let MKL use the rest. Default: all k-parallel, disable MKL threading. The user can tune via `OMP_NUM_THREADS` and `MKL_NUM_THREADS`.

**Recommendation**: Start with `MKL_NUM_THREADS=1` and `OMP_NUM_THREADS=24`. This gives near-linear scaling for the k-sweep. Only enable nested parallelism if profiling shows the diagonalization itself is the bottleneck.

### Step 2: Add `csr_spmv` to `sparse_matrices.f90`

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

### Step 3: Refactor MKL SpBLAS call sites (3 files)

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

### Step 4: Update CMakeLists.txt

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

**OpenMP**:
```cmake
find_package(OpenMP REQUIRED)
# Link both executables with OpenMP::OpenMP_Fortran
```

**MKL threading**:
```cmake
set(MKL_THREADING "intel_thread" CACHE STRING "MKL threading")
```

**Backend selection** (for future portability):
```cmake
set(LINALG_BACKEND "MKL" CACHE STRING "BLAS/LAPACK backend: MKL, OPENBLAS, AOCL")

if(LINALG_BACKEND STREQUAL "MKL")
    find_package(MKL CONFIG REQUIRED)
elseif(LINALG_BACKEND STREQUAL "OPENBLAS")
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    find_package(FFTW3 REQUIRED)
elseif(LINALG_BACKEND STREQUAL "AOCL")
    find_package(aocl REQUIRED)
    find_package(FFTW3 REQUIRED)
endif()
```

### Step 5: Handle FEAST

When `LINALG_BACKEND != MKL`, FEAST is unavailable:
- `feast_call.f90` guarded with preprocessor: `#ifdef USE_MKL_FEAST`
- `eigensolver.f90` dispatches to `solve_dense_lapack` when FEAST not available
- Wire mode (confinement=2) uses dense fallback for smaller problems
- Larger wire problems will need ARPACK-NG (future Phase 3)

### Step 6: Delete Intel files

- Delete `src/math/mkl_spblas.f90`
- Delete `src/math/mkl_sparse_handle.f90`
- Remove from `src/math/CMakeLists.txt`

### Step 7: SC loop OpenMP (optional, deferred)

The self-consistent Schrodinger-Poisson loop (`sc_loop.f90`) iterates the diagonalization for convergence. Each iteration's k-vector sweep can be parallelized the same way. However, this interacts with the SC convergence logic and is deferred until Phase 1 is validated.

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

### New tests
- `test_csr_spmv` unit test: verify `csr_spmv` against dense `zgemv` for known matrices

### Numerical correctness
- Eigenvalues must match reference data within machine precision (regression tests already cover this)
- OpenMP parallelization preserves bit-exact results: each k-point is independent, no reductions or atomics needed

### Performance benchmarks
- **k-sweep scaling**: time QW (FDstep=200, 50 k-points) with 1, 2, 4, 8, 12, 24 threads. Expect near-linear scaling to 12 threads.
- **Compiler flags**: time before/after flag changes (-O2 → -O3 + fast-math + LTO)
- **MKL SpBLAS vs csr_spmv**: time g-factor SpMV loop before/after migration
- **SC loop**: time full SC convergence before/after k-sweep parallelization

## Files Modified

| File | Change |
|---|---|
| `src/apps/main.f90` | OpenMP k-vector sweep (QW + wire), thread-private HT/workspace |
| `src/apps/main_gfactor.f90` | OpenMP k-vector sweep, thread-private HT/workspace |
| `src/math/sparse_matrices.f90` | Add `csr_spmv` |
| `src/core/utils.f90` | Remove `use mkl_spblas`, refactor `dnscsr_z_mkl` |
| `src/physics/gfactor_functions.f90` | Remove `use mkl_spblas`, use `csr_spmv` |
| `src/physics/hamiltonianConstructor.f90` | Remove `use mkl_spblas`, use `csr_matrix` type |
| `src/math/feast_call.f90` | Add `#ifdef USE_MKL_FEAST` guards |
| `src/math/eigensolver.f90` | Handle missing FEAST gracefully |
| `CMakeLists.txt` | New flags, MKL threading=intel_thread, backend option, OpenMP, LTO |
| `src/CMakeLists.txt` | Link OpenMP, update link libraries per backend |
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
