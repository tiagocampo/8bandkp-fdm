# Issue 01: Fix thread-unsafe SAVE cache in block table functions

**Priority:** Critical (C1)
**Commit:** `fix(thread-safety): pre-initialize block table caches before OpenMP parallel regions`
**Files:** `src/physics/hamiltonian_blocks.f90`, `src/physics/strain_solver.f90`, `src/apps/main.f90`

## Problem

Module-level SAVE caches (`kp_table_cached`/`kp_table_cache` in `hamiltonian_blocks.f90:68-69`, `strain_table_cached`/`strain_table_cache` in `strain_solver.f90:56-57`, `zeeman_table_cached`/`zeeman_table_cache` in `strain_solver.f90:75-76`) use a lazy-init pattern that races under OpenMP. When `ZB8bandQW` is called inside `!$omp parallel` regions (`main.f90:616` QW, `main.f90:675` Landau), multiple threads race on these SAVE variables.

## Fix

- [ ] Add pre-initialization calls in `main.f90` before the QW `!$omp parallel` region (~line 614):
  ```fortran
  ! Pre-initialize caches before OpenMP fork
  call get_kp_block_table()  ! populate kp block cache
  call get_strain_table()    ! populate strain table cache
  call get_zeeman_table()    ! populate zeeman table cache
  ```
- [ ] Add same pre-initialization before the Landau `!$omp parallel` region (~line 673)
- [ ] Alternatively, consolidate into a single subroutine in `simulation_setup.f90` or `main.f90`

## Verification

- [ ] Existing QW/Landau regression tests pass with `OMP_NUM_THREADS=4`
- [ ] No behavioral change — caches are deterministic, pre-init just moves timing
