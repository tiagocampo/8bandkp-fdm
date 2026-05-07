# Bir-Pikus Sign Fix and Buffer Overrun Fix

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix two bugs found in code review of PR #9: (1) Bir-Pikus HH/LH strain formulas use wrong sign for the hydrostatic VB shift, producing incorrect gap changes; (2) stack buffer overrun in gain calculation when numvb > 6.

**Architecture:** Two independent fixes. Bug 1 changes three sign operations in `compute_bp_scalar` from `P_eps` to `-P_eps`. Bug 2 changes a fixed-size local array to allocatable. Both fixes are localized to single functions with test updates.

**Tech Stack:** Fortran 90 (source fixes), pFUnit (test updates).

---

## Task 1: Fix Bir-Pikus sign in `compute_bp_scalar`

**Files:**
- Modify: `src/physics/strain_solver.f90:841-843`
- Modify: `src/physics/strain_solver.f90:714-715` (comment)

**Step 1: Fix the formula**

In `src/physics/strain_solver.f90`, change lines 841-843 from:

```fortran
    s%delta_EHH = P_eps + Q_eps
    s%delta_ELH = P_eps - Q_eps
    s%delta_ESO = P_eps
```

to:

```fortran
    s%delta_EHH = -P_eps + Q_eps
    s%delta_ELH = -P_eps - Q_eps
    s%delta_ESO = -P_eps
```

Note: line 844 (`s%QT2_eps = Q_eps - T_eps`) is correct as-is because `T_eps = -Q_eps` makes `QT2_eps = 2*Q_eps` regardless of the sign of the VB diagonal shifts.

**Step 2: Update the comment block**

In `src/physics/strain_solver.f90`, change lines 714-715 from:

```fortran
  !   HH (bands 1,4): P_eps + Q_eps
  !   LH (bands 2,3): P_eps - Q_eps
  !   SO (bands 5,6): P_eps
```

to:

```fortran
  !   HH (bands 1,4): -P_eps + Q_eps
  !   LH (bands 2,3): -P_eps - Q_eps
  !   SO (bands 5,6): -P_eps
```

**Step 3: Build and verify compilation**

```bash
cmake --build build 2>&1 | tail -5
```

Expected: build succeeds with no errors.

**Step 4: Commit**

```bash
git add src/physics/strain_solver.f90
git commit -m "fix: correct Bir-Pikus HH/LH sign convention in compute_bp_scalar"
```

---

## Task 2: Update strain solver tests for corrected convention

**Files:**
- Modify: `tests/unit/test_strain_solver.pf:540-543`

**Step 1: Fix test assertions**

In `tests/unit/test_strain_solver.pf`, change lines 540-543 from:

```fortran
    @assertEqual(P_eps + Q_eps, s%delta_EHH, tolerance=1.0e-12_dp)
    @assertEqual(P_eps - Q_eps, s%delta_ELH, tolerance=1.0e-12_dp)
    @assertEqual(P_eps, s%delta_ESO, tolerance=1.0e-12_dp)
    @assertEqual(2.0_dp * Q_eps, s%QT2_eps, tolerance=1.0e-12_dp)
```

to:

```fortran
    @assertEqual(-P_eps + Q_eps, s%delta_EHH, tolerance=1.0e-12_dp)
    @assertEqual(-P_eps - Q_eps, s%delta_ELH, tolerance=1.0e-12_dp)
    @assertEqual(-P_eps, s%delta_ESO, tolerance=1.0e-12_dp)
    @assertEqual(2.0_dp * Q_eps, s%QT2_eps, tolerance=1.0e-12_dp)
```

Note: `QT2_eps` assertion is unchanged since `QT2_eps = 2*Q_eps` regardless of the sign fix.

**Step 2: Build with tests and run strain tests**

```bash
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl \
    -DBUILD_TESTING=ON -DPFUNIT_DIR=$HOME/.local/pfunit/PFUNIT-4.16/cmake
cmake --build build
ctest --test-dir build -L unit -R strain -V
```

Expected: all strain tests PASS.

**Step 3: Commit**

```bash
git add tests/unit/test_strain_solver.pf
git commit -m "test: update Bir-Pikus assertions for corrected sign convention"
```

---

## Task 3: Fix buffer overrun in `find_quasi_fermi_holes`

**Files:**
- Modify: `src/physics/optical_spectra.f90:793,800-803`

**Step 1: Change fixed array to allocatable**

In `src/physics/optical_spectra.f90`, change line 793 from:

```fortran
    real(kind=dp) :: neg_Evb(6)  ! negated VB energies
```

to:

```fortran
    real(kind=dp), allocatable :: neg_Evb(:)
```

**Step 2: Add allocate/deallocate**

Change lines 800-803 from:

```fortran
    ! Negate VB eigenvalues for hole calculation
    do s = 1, numvb
      neg_Evb(s) = -eigvals(s)
    end do
```

to:

```fortran
    ! Negate VB eigenvalues for hole calculation
    allocate(neg_Evb(numvb))
    do s = 1, numvb
      neg_Evb(s) = -eigvals(s)
    end do
```

**Step 3: Add deallocate before return**

Find the end of the subroutine (the `end subroutine` line for `find_quasi_fermi_holes`). Add `deallocate(neg_Evb)` just before it:

```fortran
    deallocate(neg_Evb)
  end subroutine find_quasi_fermi_holes
```

**Step 4: Build and verify**

```bash
cmake --build build 2>&1 | tail -5
```

Expected: build succeeds.

**Step 5: Commit**

```bash
git add src/physics/optical_spectra.f90
git commit -m "fix: allocate neg_Evb dynamically to prevent buffer overrun for numvb > 6"
```

---

## Task 4: Run full test suite

**Step 1: Build with tests**

```bash
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl \
    -DBUILD_TESTING=ON -DPFUNIT_DIR=$HOME/.local/pfunit/PFUNIT-4.16/cmake
cmake --build build
```

**Step 2: Run all tests**

```bash
ctest --test-dir build -V
```

Expected: all 13+ tests PASS (8 unit + 5 regression, or more if new tests exist).

**Step 3: If tests pass, commit any remaining changes**

```bash
git status
git add -A
git commit -m "fix: strain sign convention and buffer overrun (code review PR #9)"
```
