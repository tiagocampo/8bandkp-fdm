# Modern Fortran Migration Design

**Date:** 2026-04-26
**Status:** Implementation complete in `docs/plans/2026-04-26-modern-fortran-migration-plan.md`. Merged via PR on branch `refactor/modern-fortran-migration`.
**Approach:** Phased by standard version (bottom-up)

## Current State

The codebase is mostly Fortran 2003/2008: modules everywhere, `implicit none`, allocatable components, `block` constructs, generic interfaces, `pure` functions. Not Fortran 77. But several powerful standard features are unused, and some legacy patterns remain.

### Gaps

| Missing Feature | Standard | What it enables |
|---|---|---|
| Type-bound procedures | 2003 | OOP encapsulation, methods on 17+ derived types |
| `iso_c_binding` / `iso_fortran_env` | 2003/2008 | Clean MKL/LAPACK interfaces, portable kinds |
| `do concurrent` | 2008 | Replace `forall` (deprecated F2018), GPU offloading |
| `error stop` | 2008 | Better error handling than `goto` |
| `contiguous` attribute | 2008 | Optimization hints for assumed-shape arrays |
| `elemental` procedures | 95 | Vectorize pure functions |
| C descriptors (assumed-rank) | 2018 | Generic library wrappers |
| `do concurrent` reduce | 2023 | Cleaner parallel reductions |
| Conditional expressions | 2023 | Ternary-like syntax |

### Legacy patterns to clean up

- 30 `goto` statements in `input_parser.f90`
- 13 `forall` uses (deprecated in F2018)
- 5 `external` declarations (should use explicit interfaces)
- `dsqrt` calls (legacy, just use `sqrt`)
- ~14 modules without `private` default

### stdlib assessment

Fortran stdlib v0.8.0 provides high-level wrappers for BLAS/LAPACK, matrix operations (`eig`, `eigh`, `solve`, `qr`, `cholesky`, `svd`, `schur`, `det`, `norm`), constructors (`eye`, `diag`, `kronecker_product`), and an experimental sparse OOP API. Can link against MKL as backend via `STDLIB_EXTERNAL_BLAS_I64`.

**Relevant:** constructors, norms, string handling for input parser, dense eigenvalue wrappers.
**Not relevant:** No FEAST or ARPACK wrappers; sparse API is experimental and doesn't match MKL SpBLAS performance. Performance-critical paths (CSR SpMV, FEAST, PARDISO) keep direct MKL calls.

### Compiler support

- **gfortran 16** (just released): most F2008, significant F2018, some F2023 including coarray shared memory
- **Intel ifx**: full F2018, significant F2023, GPU offloading for `do concurrent`
- **LLVM Flang**: most F2008, some F2018

---

## Phase 1: F2008 Baseline (1-2 days)

Mechanical cleanup. No logic or architecture changes. Verifiable by running existing test suite.

### 1.1 Compiler standard flag

Add `Fortran_STANDARD 2008` (or `-std=f2008`) in `CMakeLists.txt`. The code already uses `block` constructs and allocatable components, so it's de facto F2008. This makes it explicit and catches violations.

### 1.2 Replace `forall` (13 occurrences)

`forall` is deprecated in F2018 and has surprising semantics. Replace with:
- `do concurrent` where genuinely parallelizable (diagonal fills in `finitedifferences.f90`, Simpson coefficients in `utils.f90`)
- Regular `do` loops where `forall` was just array syntax (most of `hamiltonianConstructor.f90`)

### 1.3 Replace `goto` in `input_parser.f90` (30 occurrences)

All follow the same pattern: sequential `read` with `iostat` checks jumping to labeled `continue` on error. Refactor each section into a helper subroutine returning a status code, or use `block` + `exit` from a named construct.

### 1.4 Minor cleanups

- `dsqrt(3.0_dp)` -> `sqrt(3.0_dp)` in `defs.f90` (3 occurrences)
- `external :: pardiso` -> explicit interface block via `iso_c_binding`
- `complex(kind=dp), external :: zdotc` -> proper `interface` block

---

## Phase 2: Aggressive Encapsulation (5-7 days)

OOP refactor focused on hot paths. Type-bound procedures where they enable algorithmic optimizations or polymorphic dispatch.

### 2.1 `private` default everywhere

Add `private` to the ~14 modules that don't use it. Enumerate `public` exports. Verifies no accidental coupling.

### 2.2 `elemental` on scalar pure functions

Upgrade `pure` -> `elemental` on: `kronij`, `fermi_dirac`, `grid_ngrid`, and any other scalar pure functions that could be called elementally.

### 2.3 `csr_matrix` — type-bound + hidden internals

```fortran
type :: csr_matrix
  private
  integer :: nrows = 0, ncols = 0, nnz = 0
  ! ... storage details hidden
contains
  procedure :: mv            ! y = A*x (SpMV)
  procedure :: mvh           ! y = A^H*x
  procedure :: cleanup
  procedure :: norm_frobenius
  final :: csr_finalize
end type
```

Performance payoff: swap internal representation (CSR -> blocked CSR -> sell-C) without touching call sites. MKL sparse handle becomes an internal detail — optimize handle reuse across multiple SpMV calls in velocity matrix construction and optical spectra loops.

### 2.4 `spatial_grid` — hot-path encapsulation

```fortran
type :: spatial_grid
  private
contains
  procedure :: nkpoints
  procedure :: position_at
  procedure :: material_at
  procedure :: kpterms_at   ! bounds-checked, contiguous access
  final :: grid_finalize
end type
```

Performance payoff: `contiguous` on internal `kpterms` array helps vectorization in hot Hamiltonian assembly loop.

### 2.5 Polymorphic eigensolver dispatch

```fortran
type, abstract :: eigensolver_base
contains
  procedure(solve_interface), deferred :: solve
end type

type, extends(eigensolver_base) :: dense_solver    ! zheevx
type, extends(eigensolver_base) :: feast_solver    ! MKL FEAST
type, extends(eigensolver_base) :: arpack_solver   ! ARPACK
```

Performance payoff: eliminates runtime `select case` in hot path. Each solver type caches its own workspace (FEAST `fpm`, ARPACK `workd/workl`) — avoids allocate/deallocate per k-point.

### 2.6 `simulation_config` — validate-once

Bind `init` method that validates all parameters once. Downstream code trusts the result with no per-call validation overhead.

### 2.7 Finalizers everywhere

Add `final` procedures to every type with allocatable components. Prevents memory leaks in k-point sweep and SC loop iterations.

---

## Phase 3: Performance + Ecosystem (3-5 days)

### 3.1 `contiguous` attribute on assumed-shape arrays

Add `contiguous` to ~30-40 procedure arguments that always receive contiguous data. Targets: dense Hamiltonian assembly, SpMV wrappers, Simpson integration, optical spectra accumulators.

### 3.2 `do concurrent` with locality clauses (F2018)

Identify loops for parallel expression:
- k-point sweep in `main.f90` / `main_gfactor.f90` (currently OpenMP)
- velocity matrix element-wise construction in `hamiltonianConstructor.f90`
- optical spectra energy-point loop in `optical_spectra.f90`

OpenMP directives stay — `do concurrent` is complementary, communicates intent to compiler, enables GPU offloading path.

### 3.3 `iso_c_binding` for MKL interfaces

Replace `external` declarations and hand-written interface blocks in `linalg.f90` with proper `iso_c_binding` interfaces. Compile-time type checking, portable across compilers.

### 3.4 stdlib adoption (selective)

Add stdlib via CMake (`find_package(fortran_stdlib)`). Use for:

| Current code | stdlib replacement |
|---|---|
| Manual identity matrix fills | `eye(n)` |
| `simps` in `utils.f90` | `stdlib_quadrature` (check API) |
| String handling in `input_parser.f90` | `stdlib_stringlist`, `stdlib_strings` |
| Hand-rolled `zheevx` wrapper | `stdlib_linalg:eigh` for dense path |

Keep MKL SpBLAS and FEAST as direct calls.

### 3.5 fpm as alternative build (optional)

Add `fpm.toml` alongside `CMakeLists.txt`. fpm v0.13.0 supports external modules and link settings for MKL/FFTW3. Provides `fpm test` and `fpm run` for development convenience. CMake stays primary for production.

---

## Phase 4: Advanced / Exploratory (F2018/2023)

Each item needs prototyping before committing.

### 4.1 GPU offloading via `do concurrent`

Most promising target: k-point sweep (independent eigenvalue problems). Blocked by LAPACK not running on GPU. Options: move to FEAST (GPU support via oneMKL), use sparse path for all sizes, or only offload Hamiltonian construction. Eigensolver interface from Phase 2 enables future swap-in.

**Hardware constraint:** Developer machine has a Polaris GPU (Radeon RX Ellesmere, gfx803). ROCm 7.x dropped official Polaris support — `rocminfo` returns HSA errors, rocBLAS/rocFFT are unavailable for this architecture. GPU offloading code should be written portably (OpenMP `!$omp target` directives, `do concurrent`) so it compiles and runs correctly on CPU, and activates on GPU only when a supported device is present (AMD Instinct gfx90a/gfx942, NVIDIA, or Intel). Local testing is CPU-only; GPU offload verification requires access to a cluster with supported hardware.

### 4.2 `iso_fortran_env` kinds

Replace `selected_real_kind` in `defs.f90`:
```fortran
use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64, qp => real128
```

Functionally identical on common platforms. More portable and self-documenting. Add compatibility shim for platforms without `real128`.

### 4.3 F2023 conditional expressions

gfortran 16 supports conditional expressions. Adopt as encountered, don't hunt for opportunities.

### 4.4 Coarrays for distributed k-point sweeps

k-point sweep is embarrassingly parallel. Coarrays with shared-memory (gfortran 16, `-fcoarray=lib`) express this without MPI. Worth a prototype on a single multi-core machine.

### 4.5 `simple` attribute (F2023)

`simple` procedures are stricter than `pure` — only access data through arguments. Useful for functions called inside `do concurrent`. Stronger compiler optimization guarantees.

---

## Effort Summary

| Phase | Standard | Effort | Risk | Payoff |
|---|---|---|---|---|
| 1. Baseline | F2008 | 1-2 days | Very low | Clean compiler warnings, no deprecated features |
| 2. Encapsulation | F2003 OOP | 5-7 days | Low | Maintainability, polymorphic eigensolver, finalizers |
| 3. Performance + Ecosystem | F2008/2018 | 3-5 days | Low-Medium | `contiguous`, stdlib, fpm, cleaner MKL interfaces |
| 4. Advanced | F2018/2023 | Exploratory | Medium | GPU path documented, coarray prototype |

Total committed effort for Phases 1-3: ~10-14 days.

---

## Documentation Requirements

Every phase must update the following files to reflect the changes made:

- **`README.md`** — update build instructions, compiler requirements, new dependencies (stdlib, fpm), Fortran standard version, and any new build commands
- **`CLAUDE.md`** — update: code conventions section with new standard features used, architecture section if module structure changes, build commands if dependencies change, and any new gotchas or conventions introduced by the migration
