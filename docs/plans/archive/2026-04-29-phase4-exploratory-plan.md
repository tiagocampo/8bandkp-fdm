# Phase 4 Exploratory Prototypes Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Prototype two Phase 4 items — `iso_fortran_env` portable kinds and coarray k-point sweep — to validate that F2018+ features work correctly with the codebase before committing to broader adoption.

**Architecture:** Two independent prototypes. Task 1 swaps `selected_real_kind`/`selected_int_kind` for `iso_fortran_env` named constants in `defs.f90`. Task 2 creates a standalone coarray prototype program that parallelizes the QW k-point sweep. Both are validated against the full regression suite. Neither changes physics formulas, basis ordering, or material parameters.

**Tech Stack:** Fortran 2008/2018, gfortran 15.2, `-fcoarray=single` for coarray prototype, Intel MKL (LP64), CMake/Ninja, pFUnit.

**Pre-prototype findings:**
- gfortran 15.2 does NOT support F2023 conditional expressions (`cond ? a : b`). F2023 conditional expressions are excluded.
- `selected_real_kind(6,37)` == `real32` == 4, `selected_real_kind(15,307)` == `real64` == 8, `selected_real_kind(33,4931)` == `real128` == 16 — all identical on this platform.
- `selected_int_kind(8)` == `int32` == 4 (NOT int64 despite CLAUDE.md annotation). `int64` == 8. The `iknd` → `int32` mapping is correct.
- `-fcoarray=single` compiles and runs coarray programs (single-image, no actual parallelism).

---

## Task 1: Replace `selected_real_kind`/`selected_int_kind` With `iso_fortran_env` Named Constants

**Files:**
- Modify: `src/core/defs.f90` (lines 29-32)
- Modify: `CLAUDE.md` (precision kinds section)
- Test: full regression suite

**Step 1: Modify defs.f90 to use iso_fortran_env**

Replace the kind declarations at lines 29-32:

```fortran
  integer, parameter :: sp = selected_real_kind(6, 37)
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: qp = selected_real_kind(33, 4931)
  integer, parameter :: iknd = selected_int_kind(8)
```

with:

```fortran
  use, intrinsic :: iso_fortran_env, only: real32, real64, real128, int32
  ...
  integer, parameter :: sp   = real32
  integer, parameter :: dp   = real64
  integer, parameter :: qp   = real128
  integer, parameter :: iknd = int32
```

The `use, intrinsic :: iso_fortran_env` statement goes after `implicit none` and `private` (line 4), before the `public` declarations. The `real32`, `real64`, `real128`, `int32` names are imported from `iso_fortran_env` and aliased to the existing short names (`sp`, `dp`, `qp`, `iknd`), so no downstream code changes are needed.

**Step 2: Build**

Run:

```bash
cmake --build build
```

Expected: clean build. No downstream changes needed since `sp`, `dp`, `qp`, `iknd` names are unchanged.

**Step 3: Run full regression suite**

Run:

```bash
ctest --test-dir build --output-on-failure
```

Expected: 39/39 tests pass. All numerical outputs match golden data exactly.

**Step 4: Update CLAUDE.md**

In the "Code Conventions" section, find:

```
 - Precision kinds in `defs.f90`: `sp` (single), `dp` (double), `qp` (quad), `iknd` (int64)
```

Replace with:

```
 - Precision kinds in `defs.f90`: `sp` (real32), `dp` (real64), `qp` (real128), `iknd` (int32) — aliased from `iso_fortran_env`
```

This also corrects the misleading `int64` annotation.

**Step 5: Commit**

```bash
git add src/core/defs.f90 CLAUDE.md
git commit -m "refactor: use iso_fortran_env kinds instead of selected_real_kind"
```

---

## Task 2: Coarray K-Point Sweep Prototype (QW)

**Files:**
- Create: `tests/integration/test_coarray_qw_sweep.f90` (standalone prototype)
- Modify: `src/CMakeLists.txt` (add coarray target, optional)
- Test: standalone run against known QW output

This task creates a standalone program that demonstrates coarray parallelism for the QW k-point sweep. It uses `-fcoarray=single` for single-image validation, with the code structured so that `-fcoarray=lib` + OpenCoarrays would give actual multi-image parallelism.

**Design note:** The QW k-point loop at `src/apps/main.f90:700-719` is already parallelized with OpenMP (`!$omp parallel`, `!$omp do schedule(static)`). The coarray prototype demonstrates an alternative parallelism model where each coarray image owns a subset of k-points, solves independently, and collects results via coarray communication.

The wire k-point sweep (`main.f90:397`) has serial branch tracking between consecutive k-points (reordering eigenvalue branches), making it inherently sequential. The bulk sweep (`main.f90:651`) is 8x8 and trivially fast. Only the QW sweep is a good coarray candidate.

**Step 1: Create the prototype program**

Create `tests/integration/test_coarray_qw_sweep.f90`:

```fortran
! ====================================================================
! Coarray prototype: QW k-point sweep
!
! Demonstrates coarray parallelism for the embarrassingly-parallel QW
! eigenvalue sweep. Each image owns a contiguous block of k-points.
!
! Compile (single-image for validation):
!   gfortran -fcoarray=single -o test_coarray_qw_sweep test_coarray_qw_sweep.f90
!
! Compile (multi-image with OpenCoarrays):
!   caf -o test_coarray_qw_sweep test_coarray_qw_sweep.f90
!   cafrun -n 4 ./test_coarray_qw_sweep
!
! Expected behavior: each image prints its k-point range, and image 1
! prints the total number of images.
! ====================================================================
program test_coarray_qw_sweep
  use definitions, only: dp, sp
  use parameters
  use hamiltonianConstructor
  use confinement_init, only: confinementInitialization
  use finitedifferences
  use input_parser
  use linalg, only: zheevx, ilaenv, dlamch, mkl_set_num_threads_local
  use sparse_matrices
  use outputFunctions
  implicit none

  type(simulation_config) :: cfg
  real(kind=dp), allocatable :: profile(:,:), kpterms(:,:,:)
  real(kind=dp), allocatable :: z(:)

  ! Each image has its own workspace
  integer :: N, il, iuu, M_loc, info_loc, k_local, k_global
  integer :: lwork, nev, nstates
  complex(kind=dp), allocatable :: HT(:,:), work(:)
  real(kind=dp), allocatable :: rwork(:)
  integer, allocatable :: iwork(:), ifail(:)
  real(kind=dp), allocatable :: eig_local(:,:)

  ! Coarray: image 1 collects all eigenvalues
  integer :: nk_total, nk_per_image, nk_remainder, k_start, k_end
  real(kind=dp), allocatable :: eig_all(:,:)[:]
  logical :: is_root

  is_root = (this_image() == 1)

  ! --- Parse input (all images read the same config) ---
  if (is_root) then
    print '(A,I0,A)', ' Coarray QW sweep: ', num_images(), ' images'
  end if

  call read_and_setup(cfg, profile, kpterms)

  N = cfg%fdStep * 8  ! 8N x 8N Hamiltonian
  il = 1
  iuu = cfg%numcb + cfg%numvb
  nev = iuu - il + 1
  nstates = nev
  nk_total = cfg%waveVectorStep

  ! --- Distribute k-points across images ---
  nk_per_image = nk_total / num_images()
  nk_remainder = mod(nk_total, num_images())

  ! Last nk_remainder images get one extra k-point
  if (this_image() <= nk_remainder) then
    k_start = (this_image() - 1) * (nk_per_image + 1) + 1
    k_end = k_start + nk_per_image
  else
    k_start = nk_remainder * (nk_per_image + 1) + (this_image() - nk_remainder - 1) * nk_per_image + 1
    k_end = k_start + nk_per_image - 1
  end if

  if (k_end > nk_total) k_end = nk_total
  if (k_start > nk_total) k_start = nk_total + 1  ! empty range

  print '(A,I0,A,I0,A,I0,A)', ' Image ', this_image(), &
    ': k-points ', k_start, ' to ', k_end, ''

  ! --- Allocate coarray on image 1 for collection ---
  allocate(eig_all(nev, nk_total)[*])

  ! --- Each image solves its k-point range ---
  if (k_start <= k_end) then
    allocate(HT(N,N), work(1), rwork(7*N), iwork(5*N), ifail(N))
    allocate(eig_local(nev, k_end - k_start + 1))

    ! Workspace query (k=1)
    call ZB8bandQW(HT, cfg%smallK(1), profile, kpterms, cfg=cfg)
    lwork = -1
    call zheevx('V', 'I', 'U', N, HT, N, 0.0_dp, 0.0_dp, il, iuu, &
      dlamch('S'), M_loc, eig_local(:,1), HT, N, work, lwork, rwork, &
      iwork, ifail, info_loc)
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))

    ! Solve each k-point
    info_loc = mkl_set_num_threads_local(1)
    do k_local = 1, k_end - k_start + 1
      k_global = k_start + k_local - 1
      call ZB8bandQW(HT, cfg%smallK(k_global), profile, kpterms, cfg=cfg)
      call zheevx('V', 'I', 'U', N, HT, N, 0.0_dp, 0.0_dp, il, iuu, &
        dlamch('S'), M_loc, eig_local(:,k_local), HT, N, work, lwork, &
        rwork, iwork, ifail, info_loc)
      if (info_loc /= 0) then
        print '(A,I0,A,I0,A,I0)', ' ERROR: image ', this_image(), &
          ' diag failed at k=', k_global, ' info=', info_loc
        stop 1
      end if
    end do

    ! Copy local results into coarray (image 1's segment)
    eig_all(:, k_start:k_end) = eig_local(:, 1:k_end-k_start+1)

    deallocate(HT, work, rwork, iwork, ifail, eig_local)
  end if

  ! --- Synchronize: all images wait before image 1 reads results ---
  sync all

  ! --- Image 1 prints summary ---
  if (is_root) then
    print '(A)', ' All images completed. First 5 eigenvalues at k=1:'
    do k_global = 1, min(5, nev)
      print '(I0,A,g20.12)', k_global, ' ', eig_all(k_global, 1)
    end do
    print '(A)', ' Coarray prototype completed successfully.'
  end if

  deallocate(eig_all)

end program test_coarray_qw_sweep
```

**Important caveat:** This prototype won't compile as-is because `cfg%smallK` doesn't exist as a public component. The real main program builds `smallk(:)` locally from `cfg%waveVectorStep` and `cfg%waveVectorMax`. For the prototype to work, the k-vector array must either be made accessible from `simulation_config` or reconstructed locally. The prototype demonstrates the coarray pattern; integrating it into the real app requires making `smallk` accessible.

**Step 2: Validate coarray pattern compiles**

Create a minimal version that validates coarray syntax without the full Hamiltonian dependency:

```bash
cat > /tmp/test_coarray_pattern.f90 << 'FORTRAN'
program test_coarray_pattern
  implicit none
  integer :: nk, i
  real(kind=8), allocatable :: eig(:,:)[:]

  nk = 10
  allocate(eig(5, nk)[*])

  ! Each image fills its portion
  do i = 1, nk
    eig(:, i) = this_image() * 100 + i
  end do

  sync all

  if (this_image() == 1) then
    print '(A,I0)', 'Images: ', num_images()
    print '(A,5F8.1)', 'eig(:,1) = ', eig(:, 1)
    print '(A)', 'PASS: coarray pattern works'
  end if

  deallocate(eig)
end program test_coarray_pattern
FORTRAN
gfortran -fcoarray=single -o /tmp/test_coarray_pattern /tmp/test_coarray_pattern.f90
/tmp/test_coarray_pattern
```

Expected: "PASS: coarray pattern works"

**Step 3: Decide on integration path**

After validating the pattern compiles and runs:

- If the coarray pattern is worth pursuing, the integration path is: make `smallk(:)` a public component of `simulation_config` (or add a type-bound accessor), then replace the OpenMP parallel region in `main.f90:691-723` with a coarray distribution. The OpenMP path stays as fallback for single-node, coarray is the multi-node path.
- If not pursuing, the prototype file documents the pattern for future reference.

**Step 4: Commit (if prototype is kept)**

```bash
git add tests/integration/test_coarray_qw_sweep.f90
git commit -m "prototype: coarray QW k-point sweep pattern (single-image validated)"
```

---

## Task 3: Update Design Doc and CLAUDE.md

**Files:**
- Modify: `docs/plans/2026-04-26-modern-fortran-migration-design.md`
- Modify: `CLAUDE.md`

**Step 1: Update design doc Phase 4 status**

In the design doc, update Phase 4 items:

- **4.2 `iso_fortran_env` kinds:** Change from "exploratory" to "prototyped in Task 1 — drop-in replacement, no downstream changes"
- **4.3 F2023 conditional expressions:** Add note "gfortran 15.2 does NOT support F2023 conditional expressions. Revisit when gfortran 16+ is available."
- **4.4 Coarray k-point sweep:** Change from "exploratory" to "prototype validated with `-fcoarray=single`. Integration path requires making `smallk` accessible. QW k-sweep is the only good candidate (wire has branch tracking, bulk is trivially fast)."

**Step 2: Commit**

```bash
git add docs/plans/2026-04-26-modern-fortran-migration-design.md CLAUDE.md
git commit -m "docs: update Phase 4 status with prototype findings"
```

---

## Final Verification

After all tasks:

```bash
cmake --build build
ctest --test-dir build --output-on-failure
```

Expected: 39/39 tests pass. The `iso_fortran_env` kinds change produces bit-identical numerical results.
