# Strain Code Unification Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Unify the Bir-Pikus strain code so formulas and block topology exist in exactly one place each.

**Architecture:** Extract a `bp_scalar` type and `compute_bp_scalar` pure function as the single source of truth for Bir-Pikus formulas. Add `add_bp_strain_dense` helper to encode the 8x8 block topology once. Move `bir_pikus_blocks` into `simulation_config` to eliminate parameter threading. Remove deprecated `apply_pikus_bir`.

**Tech Stack:** Fortran 90, pFUnit tests, CMake build

**Design doc:** `docs/plans/2026-04-07-strain-unification-design.md`

---

## Task 1: Add `bp_scalar` type and move `bir_pikus_blocks` to defs.f90

**Files:**
- Modify: `src/core/defs.f90:39-41` (after `SQR3o2` constant, before `wavevector` type)
- Modify: `src/core/defs.f90:217` (add field to `simulation_config`)
- Modify: `src/physics/strain_solver.f90:24-66` (remove type, update public/imports)

**Step 1: Add types to defs.f90**

Insert after `UM` constant (line 41), before `type wavevector` (line 44):

```fortran
  ! Scalar Bir-Pikus strain result for a single grid point.
  ! Used by compute_bp_scalar (strain_solver) as single source of truth.
  type :: bp_scalar
    real(kind=dp) :: delta_Ec    = 0.0_dp   ! CB: ac*Tr(eps)
    real(kind=dp) :: delta_EHH   = 0.0_dp   ! HH: P_eps + Q_eps
    real(kind=dp) :: delta_ELH   = 0.0_dp   ! LH: P_eps - Q_eps
    real(kind=dp) :: delta_ESO   = 0.0_dp   ! SO: P_eps
    complex(kind=dp) :: R_eps    = cmplx(0.0_dp, 0.0_dp, kind=dp)
    complex(kind=dp) :: S_eps    = cmplx(0.0_dp, 0.0_dp, kind=dp)
    real(kind=dp) :: QT2_eps     = 0.0_dp   ! Q_eps - T_eps = 2*Q_eps
  end type bp_scalar

  ! Per-grid-point Bir-Pikus strain Hamiltonian components.
  ! Arrays have length ngrid; unallocated when unstrained.
  type :: bir_pikus_blocks
    real(kind=dp), allocatable :: delta_Ec(:)    ! CB bands 7,8
    real(kind=dp), allocatable :: delta_EHH(:)   ! HH bands 1,4
    real(kind=dp), allocatable :: delta_ELH(:)   ! LH bands 2,3
    real(kind=dp), allocatable :: delta_ESO(:)   ! SO bands 5,6
    complex(kind=dp), allocatable :: R_eps(:)    ! VB-VB coupling
    complex(kind=dp), allocatable :: S_eps(:)    ! VB-VB coupling
    real(kind=dp), allocatable :: QT2_eps(:)     ! VB-SO coupling (2*Q_eps)
  end type bir_pikus_blocks
```

**Step 2: Add `strain_blocks` field to `simulation_config`**

Insert before `end type simulation_config` (line 218):

```fortran
    type(bir_pikus_blocks)     :: strain_blocks    ! precomputed strain data
```

**Step 3: Remove `bir_pikus_blocks` type from strain_solver.f90**

Remove the type definition at lines 45-66 (everything from the comment block starting `! Per-grid-point` through `end type bir_pikus_blocks`).

Remove `bir_pikus_blocks` from the public list at line 29.

Remove the `use definitions, only: ...` partial import that includes strain types -- since `bir_pikus_blocks` is now in `definitions` module, strain_solver.f90 already gets it via `use definitions`.

**Step 4: Build to verify compilation**

Run: `cmake --build build 2>&1 | tail -5`
Expected: Build succeeds. All types now in defs.f90.

**Step 5: Commit**

```
feat: add bp_scalar type, move bir_pikus_blocks to defs.f90
```

---

## Task 2: Add `compute_bp_scalar` pure function

**Files:**
- Modify: `src/physics/strain_solver.f90` (add function after `bir_pikus_blocks_free`)
- Modify: `tests/unit/test_strain_solver.pf` (add test)

**Step 1: Write the test**

Add to `tests/unit/test_strain_solver.pf` before `end module test_strain_solver`:

```fortran
  @test
  subroutine test_compute_bp_scalar_gaas_biaxial()
    ! Verify compute_bp_scalar matches the canonical Bir-Pikus formulas
    ! for [001] biaxial strain with GaAs parameters.

    type(paramStruct) :: params(1)
    character(len=255) :: materials(1)
    type(bp_scalar) :: s
    real(kind=dp) :: eps_xx, eps_zz, Tr_eps, P_eps, Q_eps

    materials(1) = "GaAs"
    call paramDatabase(materials, 1, params)

    ! Tensile biaxial strain (GaAs on InP-like substrate)
    eps_xx = 0.0382_dp
    eps_zz = -0.0354_dp
    Tr_eps = 2.0_dp * eps_xx + eps_zz

    s = compute_bp_scalar(params(1), eps_xx, eps_xx, eps_zz, &
      0.0_dp, 0.0_dp, 0.0_dp)

    P_eps = -params(1)%av * Tr_eps
    Q_eps = params(1)%b_dp * 0.5_dp * (eps_zz - eps_xx)

    @assertEqual(params(1)%ac * Tr_eps, s%delta_Ec, tolerance=1.0e-12_dp)
    @assertEqual(P_eps + Q_eps, s%delta_EHH, tolerance=1.0e-12_dp)
    @assertEqual(P_eps - Q_eps, s%delta_ELH, tolerance=1.0e-12_dp)
    @assertEqual(P_eps, s%delta_ESO, tolerance=1.0e-12_dp)
    @assertEqual(2.0_dp * Q_eps, s%QT2_eps, tolerance=1.0e-12_dp)

    ! For [001] biaxial (eps_xx=eps_yy, eps_xy=0): R_eps=0
    @assertTrue(abs(s%R_eps) < 1.0e-15_dp)
    @assertTrue(abs(s%S_eps) < 1.0e-15_dp)

  end subroutine test_compute_bp_scalar_gaas_biaxial
```

**Step 2: Build test to verify it fails**

Run: `cmake --build build 2>&1 | tail -5`
Expected: Compilation error -- `compute_bp_scalar` not defined yet.

**Step 3: Implement `compute_bp_scalar`**

Add to `strain_solver.f90` after `bir_pikus_blocks_free`, inside `contains`:

```fortran
  pure function compute_bp_scalar(params, eps_xx, eps_yy, eps_zz, &
      eps_xy, eps_xz, eps_yz) result(s)
    type(paramStruct), intent(in) :: params
    real(kind=dp), intent(in) :: eps_xx, eps_yy, eps_zz
    real(kind=dp), intent(in) :: eps_xy, eps_xz, eps_yz
    type(bp_scalar) :: s

    real(kind=dp) :: Tr_eps, P_eps, Q_eps, T_eps

    Tr_eps = eps_xx + eps_yy + eps_zz
    P_eps = -params%av * Tr_eps
    Q_eps = params%b_dp * 0.5_dp * (eps_zz - 0.5_dp * (eps_yy + eps_xx))
    T_eps = -Q_eps

    s%delta_Ec  = params%ac * Tr_eps
    s%delta_EHH = P_eps + Q_eps
    s%delta_ELH = P_eps - Q_eps
    s%delta_ESO = P_eps
    s%QT2_eps   = Q_eps - T_eps

    s%R_eps = -SQR3 * (params%b_dp * 0.5_dp * (eps_xx - eps_yy) &
                       - IU * params%d_dp * eps_xy)
    s%S_eps = IU * 2.0_dp * SQR3 * params%d_dp * &
              cmplx(eps_xz, -eps_yz, kind=dp)

  end function compute_bp_scalar
```

Also add `compute_bp_scalar` to the public list in strain_solver.f90.

**Step 4: Build and run the test**

Run: `cmake --build build && ctest --test-dir build -R test_strain_solver -V`
Expected: All strain solver tests pass including the new one.

**Step 5: Commit**

```
feat: add compute_bp_scalar pure function with test
```

---

## Task 3: Refactor `compute_bir_pikus_blocks` to use `compute_bp_scalar`

**Files:**
- Modify: `src/physics/strain_solver.f90:882-927` (the main loop body)

**Step 1: Replace the inline formula computation**

Replace the loop body (lines 910-927) with a call to `compute_bp_scalar`:

```fortran
    do ij = 1, ngrid
      mid = material_id(ij)
      if (mid < 1 .or. mid > size(params)) then
        if (.not. warned_invalid_mat) then
          print *, 'WARNING: invalid material_id=', mid, &
            ' at grid point', ij, '(range: 1..', size(params), '). Skipping.'
          warned_invalid_mat = .true.
        end if
        cycle
      end if

      eps_xx = strain_out%eps_xx(ij)
      eps_yy = strain_out%eps_yy(ij)
      eps_zz = strain_out%eps_zz(ij)
      eps_yz = strain_out%eps_yz(ij)
      eps_xy = 0.0_dp
      eps_xz = 0.0_dp

      ! Skip negligible strain
      Tr_eps = eps_xx + eps_yy + eps_zz
      if (abs(Tr_eps) < 1.0e-20_dp .and. abs(eps_yz) < 1.0e-20_dp) cycle

      ! Single source of truth for Bir-Pikus formulas
      associate(s => compute_bp_scalar(params(mid), eps_xx, eps_yy, eps_zz, &
                                        eps_xy, eps_xz, eps_yz))
        bp%delta_Ec(ij)  = s%delta_Ec
        bp%delta_EHH(ij) = s%delta_EHH
        bp%delta_ELH(ij) = s%delta_ELH
        bp%delta_ESO(ij) = s%delta_ESO
        bp%R_eps(ij)     = s%R_eps
        bp%S_eps(ij)     = s%S_eps
        bp%QT2_eps(ij)   = s%QT2_eps
      end associate
    end do
```

Remove the now-unused local variables `Tr_eps`, `P_eps`, `Q_eps`, `T_eps`, `av`, `ac`, `b_dp`, `d_dp` from the subroutine declaration.

**Step 2: Build and run tests**

Run: `cmake --build build && ctest --test-dir build -R test_strain_solver -V`
Expected: All strain tests pass (same behavior, thinner code).

**Step 3: Commit**

```
refactor: compute_bir_pikus_blocks delegates to compute_bp_scalar
```

---

## Task 4: Add `add_bp_strain_dense` helper

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90` (add new private subroutine)
- Remove: `use strain_solver, only: bir_pikus_blocks` from line 7 (type now in defs)

**Step 1: Add the helper subroutine**

Add a private subroutine `add_bp_strain_dense` in hamiltonianConstructor module, after the existing `insert_profile_diagonal`. This encodes the 8x8 Bir-Pikus block topology once:

```fortran
    subroutine add_bp_strain_dense(HT, ii, N, bp)
      complex(kind=dp), intent(inout) :: HT(:,:)
      integer, intent(in) :: ii, N
      type(bir_pikus_blocks), intent(in) :: bp

      complex(kind=dp) :: R_eps_c, S_eps_c

      R_eps_c = conjg(bp%R_eps(ii))
      S_eps_c = conjg(bp%S_eps(ii))

      ! === Diagonal per-band ===
      HT(      ii,      ii) = HT(      ii,      ii) + bp%delta_EHH(ii)
      HT(  N + ii,  N + ii) = HT(  N + ii,  N + ii) + bp%delta_ELH(ii)
      HT(2*N + ii,2*N + ii) = HT(2*N + ii,2*N + ii) + bp%delta_ELH(ii)
      HT(3*N + ii,3*N + ii) = HT(3*N + ii,3*N + ii) + bp%delta_EHH(ii)
      HT(4*N + ii,4*N + ii) = HT(4*N + ii,4*N + ii) + bp%delta_ESO(ii)
      HT(5*N + ii,5*N + ii) = HT(5*N + ii,5*N + ii) + bp%delta_ESO(ii)
      HT(6*N + ii,6*N + ii) = HT(6*N + ii,6*N + ii) + bp%delta_Ec(ii)
      HT(7*N + ii,7*N + ii) = HT(7*N + ii,7*N + ii) + bp%delta_Ec(ii)

      ! === Off-diagonal: S_eps (HH-LH) ===
      HT(      ii,  N + ii) = HT(      ii,  N + ii) + S_eps_c
      HT(  N + ii,      ii) = HT(  N + ii,      ii) + bp%S_eps(ii)
      HT(2*N + ii,3*N + ii) = HT(2*N + ii,3*N + ii) - S_eps_c
      HT(3*N + ii,2*N + ii) = HT(3*N + ii,2*N + ii) - bp%S_eps(ii)

      ! === Off-diagonal: R_eps (HH-LH) ===
      HT(      ii,2*N + ii) = HT(      ii,2*N + ii) + R_eps_c
      HT(2*N + ii,      ii) = HT(2*N + ii,      ii) + bp%R_eps(ii)
      HT(  N + ii,3*N + ii) = HT(  N + ii,3*N + ii) + R_eps_c
      HT(3*N + ii,  N + ii) = HT(3*N + ii,  N + ii) + bp%R_eps(ii)

      ! === Off-diagonal: VB-SO coupling ===
      HT(      ii,4*N + ii) = HT(      ii,4*N + ii) - IU * RQS2 * S_eps_c
      HT(4*N + ii,      ii) = HT(4*N + ii,      ii) + IU * RQS2 * bp%S_eps(ii)
      HT(      ii,5*N + ii) = HT(      ii,5*N + ii) + IU * SQR2 * R_eps_c
      HT(5*N + ii,      ii) = HT(5*N + ii,      ii) - IU * SQR2 * bp%R_eps(ii)

      HT(  N + ii,4*N + ii) = HT(  N + ii,4*N + ii) + IU * RQS2 * bp%QT2_eps(ii)
      HT(4*N + ii,  N + ii) = HT(4*N + ii,  N + ii) - IU * RQS2 * bp%QT2_eps(ii)
      HT(  N + ii,5*N + ii) = HT(  N + ii,5*N + ii) - IU * SQR3o2 * S_eps_c
      HT(5*N + ii,  N + ii) = HT(5*N + ii,  N + ii) + IU * SQR3o2 * bp%S_eps(ii)

      HT(2*N + ii,4*N + ii) = HT(2*N + ii,4*N + ii) + IU * SQR3o2 * bp%S_eps(ii)
      HT(4*N + ii,2*N + ii) = HT(4*N + ii,2*N + ii) - IU * SQR3o2 * S_eps_c
      HT(2*N + ii,5*N + ii) = HT(2*N + ii,5*N + ii) + IU * RQS2 * bp%QT2_eps(ii)
      HT(5*N + ii,2*N + ii) = HT(5*N + ii,2*N + ii) - IU * RQS2 * bp%QT2_eps(ii)

      HT(3*N + ii,4*N + ii) = HT(3*N + ii,4*N + ii) - IU * SQR2 * bp%R_eps(ii)
      HT(4*N + ii,3*N + ii) = HT(4*N + ii,3*N + ii) + IU * SQR2 * R_eps_c
      HT(3*N + ii,5*N + ii) = HT(3*N + ii,5*N + ii) + IU * RQS2 * bp%S_eps(ii)
      HT(5*N + ii,3*N + ii) = HT(5*N + ii,3*N + ii) - IU * RQS2 * S_eps_c
    end subroutine add_bp_strain_dense
```

**Step 2: Remove `use strain_solver` import from hamiltonianConstructor.f90**

Line 7: Remove `use strain_solver, only: bir_pikus_blocks` -- the type is now in `definitions` which is already imported.

**Step 3: Build to verify**

Run: `cmake --build build 2>&1 | tail -5`
Expected: Build succeeds (helper not yet called, just compiled).

**Step 4: Commit**

```
feat: add add_bp_strain_dense helper (single source of truth for block topology)
```

---

## Task 5: Refactor `ZB8bandBulk` to use `compute_bp_scalar` + `add_bp_strain_dense`

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90:1086-1179` (bulk strain block)
- Add `use strain_solver, only: compute_bp_scalar` to hamiltonianConstructor imports

**Step 1: Replace the inline strain block**

Replace the entire `if (params(1)%strainSubstrate > 0.0_dp) then ... end if` block (lines 1086-1179) with:

```fortran
      if (params(1)%strainSubstrate > 0.0_dp) then
        block
          real(kind=dp) :: a_film, eps_xx, eps_zz
          type(bir_pikus_blocks) :: bp_bulk

          a_film = params(1)%a0
          if (a_film > 0.0_dp) then
            ! Biaxial [001] strain tensor
            eps_xx = (params(1)%strainSubstrate - a_film) / a_film
            eps_zz = -2.0_dp * params(1)%C12 / params(1)%C11 * eps_xx

            ! Compute Bir-Pikus for single grid point (N=1)
            allocate(bp_bulk%delta_Ec(1), bp_bulk%delta_EHH(1), &
              bp_bulk%delta_ELH(1), bp_bulk%delta_ESO(1), &
              bp_bulk%R_eps(1), bp_bulk%S_eps(1), bp_bulk%QT2_eps(1))

            associate(s => compute_bp_scalar(params(1), eps_xx, eps_xx, eps_zz, &
                                              0.0_dp, 0.0_dp, 0.0_dp))
              bp_bulk%delta_Ec(1)  = s%delta_Ec
              bp_bulk%delta_EHH(1) = s%delta_EHH
              bp_bulk%delta_ELH(1) = s%delta_ELH
              bp_bulk%delta_ESO(1) = s%delta_ESO
              bp_bulk%R_eps(1)     = s%R_eps
              bp_bulk%S_eps(1)     = s%S_eps
              bp_bulk%QT2_eps(1)   = s%QT2_eps
            end associate

            call add_bp_strain_dense(HT, 1, 1, bp_bulk)
            call bir_pikus_blocks_free(bp_bulk)
          end if
        end block
      end if
```

Add to the module imports: `use strain_solver, only: compute_bp_scalar, bir_pikus_blocks_free`

**Step 2: Build and run regression tests for bulk**

Run: `cmake --build build && ctest --test-dir build -R regression_bulk`
Expected: Both bulk regression tests pass (eigenvalues unchanged).

**Step 3: Commit**

```
refactor: ZB8bandBulk uses compute_bp_scalar + add_bp_strain_dense
```

---

## Task 6: Refactor `ZB8bandQW` to use `add_bp_strain_dense` + `cfg%strain_blocks`

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90:623-635` (signature)
- Modify: `src/physics/hamiltonianConstructor.f90:808-875` (strain insertion)

**Step 1: Remove `strain_blocks` parameter from `ZB8bandQW`**

Change signature from:
```fortran
  subroutine ZB8bandQW(HT, wv, profile, kpterms, sparse, HT_csr, g, strain_blocks)
```
to:
```fortran
  subroutine ZB8bandQW(HT, wv, profile, kpterms, sparse, HT_csr, g)
```

Remove the declaration:
```fortran
  type(bir_pikus_blocks), intent(in), optional :: strain_blocks
```

**Step 2: Replace the 68-line strain insertion block (lines 808-875) with:**

```fortran
      ! Full Bir-Pikus strain (k-independent, not in g-mode)
      if (.not. present(g)) then
        if (allocated(cfg%strain_blocks%delta_Ec)) then
          do ii = 1, N
            call add_bp_strain_dense(HT, ii, N, cfg%strain_blocks)
          end do
        end if
      end if
```

Wait -- `ZB8bandQW` doesn't currently receive `cfg`. It receives `profile` and `kpterms` but not the full config. We need to add `cfg` as a parameter.

Actually, check: `ZB8bandQW` currently does NOT have `cfg`. The wire `ZB8bandGeneralized` does. So for QW, we have two options:
1. Add `cfg` parameter to `ZB8bandQW`
2. Pass `strain_blocks` directly from the caller

Since the design says to eliminate threading and use `cfg`, add `cfg` to the signature. This is the cleanest approach.

**Revised Step 1:** Change signature to add `cfg` and remove `strain_blocks`:

```fortran
  subroutine ZB8bandQW(HT, wv, profile, kpterms, cfg, sparse, HT_csr, g)
    ...
    type(simulation_config), intent(in), optional :: cfg
```

**Step 2:** Replace the 68-line insertion with the helper call.

**Step 3: Update all call sites of `ZB8bandQW`**

Search for `call ZB8bandQW(` and add `cfg=cfg` argument, remove `strain_blocks=bp` argument. Call sites are in:
- `src/apps/main.f90` (workspace query, k-sweep loop)
- `src/apps/main_gfactor.f90` (band loop)
- `src/physics/sc_loop.f90` (SC inner loop)

**Step 4: Build and run QW regression tests**

Run: `cmake --build build && ctest --test-dir build -R regression_qw`
Expected: QW regression passes (same eigenvalues).

**Step 5: Commit**

```
refactor: ZB8bandQW uses add_bp_strain_dense + cfg%strain_blocks
```

---

## Task 7: Refactor `ZB8bandGeneralized` to use `cfg%strain_blocks`

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90:1203` (signature)
- Modify: `src/physics/hamiltonianConstructor.f90:1212` (remove parameter)

**Step 1: Remove `strain_blocks` from signature**

Change line 1203 from:
```fortran
  subroutine ZB8bandGeneralized(HT_csr, kz, profile_2d, kpterms_2d, cfg, coo_cache, g, strain_blocks)
```
to:
```fortran
  subroutine ZB8bandGeneralized(HT_csr, kz, profile_2d, kpterms_2d, cfg, coo_cache, g)
```

Remove line 1212: `type(bir_pikus_blocks), intent(in), optional :: strain_blocks`

**Step 2: Replace `strain_blocks` references with `cfg%strain_blocks`**

In the nnz_est section (~line 1542):
```fortran
  ! Change from:
  if (present(strain_blocks)) then
    if (allocated(strain_blocks%delta_Ec)) then
      nnz_est = nnz_est + 32 * N
  ! To:
  if (allocated(cfg%strain_blocks%delta_Ec)) then
    nnz_est = nnz_est + 32 * N
```

In the COO insertion section (~line 1753):
```fortran
  ! Change from:
  if (present(strain_blocks)) then
    if (allocated(strain_blocks%delta_Ec)) then
      call insert_strain_coo(..., strain_blocks, N)
  ! To:
  if (allocated(cfg%strain_blocks%delta_Ec)) then
    call insert_strain_coo(..., cfg%strain_blocks, N)
```

**Step 3: Update all call sites**

Search for `call ZB8bandGeneralized(` and remove `strain_blocks=bp` argument. Call sites are in:
- `src/apps/main.f90` (wire k-sweep)
- `src/apps/main_gfactor.f90` (wire at kz=0)
- `src/physics/sc_loop.f90` (wire SC loop)

**Step 4: Build and run wire regression test**

Run: `cmake --build build && ctest --test-dir build -R regression_wire`
Expected: Wire regression passes.

**Step 5: Commit**

```
refactor: ZB8bandGeneralized uses cfg%strain_blocks directly
```

---

## Task 8: Refactor sc_loop to use `cfg%strain_blocks`

**Files:**
- Modify: `src/physics/sc_loop.f90:15` (remove import)
- Modify: `src/physics/sc_loop.f90:53-54` (self_consistent_loop signature)
- Modify: `src/physics/sc_loop.f90:559-561` (self_consistent_loop_wire signature)
- Modify: `src/physics/sc_loop.f90:214` (ZB8bandQW call site)
- Modify: `src/physics/sc_loop.f90:704-710` (ZB8bandGeneralized call sites)

**Step 1: Remove `use strain_solver, only: bir_pikus_blocks`** from line 15.

**Step 2: Remove `strain_blocks` parameter from both subroutines**

For `self_consistent_loop`:
```fortran
  ! Remove from signature:
  type(bir_pikus_blocks), intent(in), optional :: strain_blocks
```

For `self_consistent_loop_wire`:
```fortran
  ! Remove from signature:
  type(bir_pikus_blocks), intent(in), optional :: strain_blocks
```

**Step 3: Update ZB8bandQW call inside SC loop (~line 214)**

```fortran
  ! Change from:
  call ZB8bandQW(HT, wv, profile, kpterms, strain_blocks=strain_blocks)
  ! To:
  call ZB8bandQW(HT, wv, profile, kpterms, cfg=cfg)
```

**Step 4: Update ZB8bandGeneralized calls inside wire SC loop (~lines 704-710)**

Remove `strain_blocks=strain_blocks` from both calls.

**Step 5: Build and run SC regression tests**

Run: `cmake --build build && ctest --test-dir build -R regression_sc`
Expected: SC tests pass (same physics).

**Step 6: Commit**

```
refactor: sc_loop uses cfg%strain_blocks, no parameter threading
```

---

## Task 9: Update main.f90 and main_gfactor.f90

**Files:**
- Modify: `src/apps/main.f90`
- Modify: `src/apps/main_gfactor.f90`

**Step 1: Replace local `bp` variable with `cfg%strain_blocks`**

In `main.f90`:
- Remove `type(bir_pikus_blocks) :: bp` declaration
- Change `call compute_bir_pikus_blocks(..., bp)` to `call compute_bir_pikus_blocks(..., cfg%strain_blocks)`
- Remove `strain_blocks=bp` from all call sites (already done in Tasks 6-8)
- Remove `call bir_pikus_blocks_free(bp)` if present

Same changes in `main_gfactor.f90`.

**Step 2: Build and run ALL tests**

Run: `cmake --build build && ctest --test-dir build`
Expected: All 24 tests pass.

**Step 3: Commit**

```
refactor: main programs populate cfg%strain_blocks, remove local bp
```

---

## Task 10: Remove `apply_pikus_bir` and update tests

**Files:**
- Modify: `src/physics/strain_solver.f90` (remove subroutine + public)
- Modify: `tests/unit/test_strain_solver.pf` (remove/update tests)

**Step 1: Remove `apply_pikus_bir` from strain_solver.f90**

- Remove `public :: apply_pikus_bir` from line 28
- Remove the entire subroutine (currently marked DEPRECATED)

**Step 2: Update or remove tests that call `apply_pikus_bir`**

The tests `test_apply_pikus_bir_qw` and `test_apply_pikus_bir_wire` (approximately lines 100-300 in test_strain_solver.pf) call the removed function. Either:
- Remove them (the Bir-Pikus physics is already tested via `test_bir_pikus_blocks_diagonal` and `test_compute_bp_scalar_gaas_biaxial`)
- Or rewrite them to test `compute_bir_pikus_blocks` + Hamiltonian insertion

Prefer removing them since `compute_bir_pikus_blocks` tests already cover the physics.

**Step 3: Build and run all tests**

Run: `cmake --build build && ctest --test-dir build`
Expected: All tests pass.

**Step 4: Commit**

```
refactor: remove deprecated apply_pikus_bir and its tests
```

---

## Task 11: Final verification

**Step 1: Clean build from scratch**

```bash
rm -rf build
cmake -G Ninja -B build -DMKL_DIR=$MKL_DIR/lib/cmake/mkl -DBUILD_TESTING=ON -DPFUNIT_DIR=$HOME/.local/pfunit/PFUNIT-4.16/cmake
cmake --build build
```

**Step 2: Run all tests**

```bash
ctest --test-dir build
```

Expected: 100% pass (24/24 or more depending on test count).

**Step 3: Verify no regressions in physics output**

```bash
python tests/regression/compare_output.py output/eigenvalues.dat tests/regression/data/bulk_gaas_kx_ref.dat
```

Expected: No output (all values match within tolerance).
