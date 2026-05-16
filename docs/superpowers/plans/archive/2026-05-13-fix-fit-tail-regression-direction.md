# Fix fit_tail_exponential Regression Direction for Right-Edge States

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix `fit_tail_exponential` so the regression always runs over a decaying density region, even when the tail is found to the left of the peak (right-edge Majorana case).

**Architecture:** Add a `forward_tail` boolean to track search direction. Split the regression loop into two branches — forward (current behavior, left-edge) and backward (new, right-edge). The backward branch iterates from `tail_start` down to 1, using `x_start - positions(i)` as the distance so the slope is negative (decaying).

**Tech Stack:** Fortran 2018, pFUnit tests, ctest

---

### Task 1: Write failing test for right-edge exponential decay

**Files:**
- Modify: `tests/unit/test_edge_states.pf` (append new test before `end module`)

- [ ] **Step 1: Write the failing test**

Append this test to `tests/unit/test_edge_states.pf`, before the final `end module test_edge_states`:

```fortran
  @test
  subroutine test_fit_exponential_decay_right_edge()
    ! Right-edge density: peak near end, tail to the left.
    ! Old bug: backward search finds tail_start, but regression runs forward
    ! over increasing density, producing meaningless xi.
    real(kind=dp), parameter :: xi_true = 5.0_dp
    real(kind=dp), parameter :: dz = 1.0_dp
    real(kind=dp) :: density(100), xvals(100)
    real(kind=dp) :: xi
    logical :: success
    integer :: i

    do i = 1, 100
      xvals(i) = real(i - 1, kind=dp) * dz
    end do

    ! Peak at index 80, exponential tail to the LEFT
    density = 0.0_dp
    density(80) = 1.0_dp
    do i = 79, 21, -1
      density(i) = exp(-2.0_dp * real(80 - i, kind=dp) / xi_true)
    end do
    ! density(1:20) stays 0.0 (before the tail)

    call fit_exponential_decay(density, xvals, xi, success)

    @assertTrue(success, message="fit_exponential_decay should succeed for right-edge profile")
    @assertTrue(xi > 0.0_dp, message="xi should be positive for right-edge profile")
    ! Recovered xi should be within 25% of xi_true/2 (exp uses -2*x/xi)
    @assertTrue(abs(xi - xi_true * 0.5_dp) < 1.0_dp, &
      & message="right-edge xi should match left-edge xi (within 25%)")
  end subroutine test_fit_exponential_decay_right_edge
```

- [ ] **Step 2: Build and run the test to verify it fails**

Run:
```bash
cmake --build build 2>&1 | tail -5
```

Then:
```bash
ctest --test-dir build -R test_edge_states --output-on-failure 2>&1 | tail -20
```

Expected: `test_fit_exponential_decay_right_edge` FAILS — the right-edge xi will be far from `xi_true/2` because the regression runs over increasing density.

- [ ] **Step 3: Commit the failing test**

```bash
git add tests/unit/test_edge_states.pf
git commit -m "test: add failing test for right-edge Majorana xi regression"
```

---

### Task 2: Fix fit_tail_exponential regression direction

**Files:**
- Modify: `src/physics/topological_analysis.f90:467-555` (`fit_tail_exponential` subroutine)

- [ ] **Step 1: Add `forward_tail` variable and track search direction**

In the variable declarations (line 478), add `forward_tail`:

```fortran
    integer :: peak_idx, tail_start, i
    logical :: forward_tail
    real(kind=dp) :: rho_peak, x_start, domain_extent
```

After the forward search block (line 501), set `forward_tail`:

```fortran
    ! Search forward from peak for tail start
    tail_start = 0
    forward_tail = .true.
    do i = peak_idx + 1, n_points
      if (density(i) < threshold * rho_peak) then
        tail_start = i
        exit
      end if
    end do

    ! Search backward if no forward tail found
    if (tail_start == 0) then
      forward_tail = .false.
      do i = peak_idx - 1, 1, -1
        if (density(i) < threshold * rho_peak) then
          tail_start = i
          exit
        end if
      end do
    end if
```

- [ ] **Step 2: Fix the early-return guard for backward tails**

Replace the guard at line 513:

Old:
```fortran
    if (tail_start == 0 .or. tail_start >= n_points) return
```

New:
```fortran
    if (tail_start == 0) return
    if (forward_tail .and. tail_start >= n_points) return
    if (.not. forward_tail .and. tail_start <= 1) return
```

- [ ] **Step 3: Split the regression loop by direction**

Replace the single regression loop (lines 524-532) with two branches:

```fortran
    if (forward_tail) then
      ! Tail is to the right of peak: iterate forward
      do i = tail_start, n_points
        if (density(i) > 1.0e-14_dp) then
          n_fit_actual = n_fit_actual + 1
          sum_y_log = sum_y_log + log(density(i))
          sum_x = sum_x + (positions(i) - x_start)
          sum_x2 = sum_x2 + (positions(i) - x_start)**2
          sum_xy = sum_xy + (positions(i) - x_start) * log(density(i))
        end if
      end do
    else
      ! Tail is to the left of peak: iterate backward toward index 1
      do i = tail_start, 1, -1
        if (density(i) > 1.0e-14_dp) then
          n_fit_actual = n_fit_actual + 1
          sum_y_log = sum_y_log + log(density(i))
          sum_x = sum_x + (x_start - positions(i))
          sum_x2 = sum_x2 + (x_start - positions(i))**2
          sum_xy = sum_xy + (x_start - positions(i)) * log(density(i))
        end if
      end do
    end if
```

- [ ] **Step 4: Fix the domain_extent calculation for R11 convergence check**

Replace the domain_extent line (line 545):

Old:
```fortran
    domain_extent = abs(positions(n_points) - positions(tail_start))
```

New:
```fortran
    if (forward_tail) then
      domain_extent = abs(positions(n_points) - positions(tail_start))
    else
      domain_extent = abs(positions(tail_start) - positions(1))
    end if
```

- [ ] **Step 5: Build and run tests**

```bash
cmake --build build 2>&1 | tail -5
```

```bash
ctest --test-dir build -R test_edge_states --output-on-failure 2>&1 | tail -20
```

Expected: ALL tests pass, including `test_fit_exponential_decay_right_edge`.

- [ ] **Step 6: Run full topology test suite**

```bash
ctest --test-dir build -R "z2|topology|bdg|edge" --output-on-failure 2>&1 | tail -20
```

Expected: All pass.

- [ ] **Step 7: Commit**

```bash
git add src/physics/topological_analysis.f90
git commit -m "fix(topology): correct fit_tail_exponential regression direction for right-edge states"
```

---

### Task 3: Verify existing tests are unchanged

**Files:** None (verification only)

- [ ] **Step 1: Run full test suite**

```bash
ctest --test-dir build -j4 --output-on-failure 2>&1 | tail -10
```

Expected: All tests pass. The left-edge Majorana tests (`test_fit_exponential_decay`, `test_fit_exponential_decay_accuracy`) continue to pass unchanged because the forward-tail path is identical to the old code.
