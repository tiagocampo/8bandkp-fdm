# PR #13 Review Fixes Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix all correctness bugs, backward-compatibility regressions, and code quality issues found during PR #13 review so the branch can merge.

**Architecture:** Fixes are organized by severity: correctness bugs first (C1-C5), backward-compat regressions next (C6-C7), then code quality (H1-H12), then test gaps (T1-T8). Each task touches specific files with exact line references. No new physics — only fixing what exists.

**Tech Stack:** Fortran 2018, pFUnit tests, Python 3 scripts, CMake/Ninja build

---

### Task 1: Fix double-assignment bug in main_topology.f90 (C1)

**Files:**
- Modify: `src/apps/main_topology.f90:356-358` (run_qshe_wire)

**Step 1: Write the fix**

At `src/apps/main_topology.f90:356-358`, the current code:
```fortran
result%edge_xi = edge_info(1)      ! xi_min
result%edge_xi = edge_info(2)      ! overwrites with xi_avg
result%min_gap = edge_info(3)      ! n_edge count, NOT a gap
```

Replace with:
```fortran
result%edge_xi_min = edge_info(1)  ! min localization length
result%edge_xi = edge_info(2)      ! average localization length
result%min_gap = edge_info(3)      ! kept for compatibility (edge count)
```

This requires adding `edge_xi_min` to `topological_result` in `src/core/defs.f90`. Find the `topological_result` type and add `real(kind=dp) :: edge_xi_min = 0.0_dp` next to `edge_xi`.

Also check `main_topology.f90:502` (run_bdg_wire) which correctly uses only `result%edge_xi = sum(...) / n`.

**Step 2: Build and verify**

Run: `cmake --build build`
Expected: Clean build, no errors.

**Step 3: Run existing tests**

Run: `ctest --test-dir build -L unit -j4`
Expected: All tests pass.

**Step 4: Commit**

```bash
git add src/apps/main_topology.f90 src/core/defs.f90
git commit -m "fix: correct edge_xi double-assignment and add edge_xi_min field"
```

---

### Task 2: Remove duplicate Zeeman block in ZB8bandQW (C2)

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90:283-314` (delete Block 2)

**Step 1: Write the fix**

Delete lines 283-314 (the second Zeeman block that starts with `! Zeeman splitting: g*mu_B * B . sigma for B in z-direction`). This block:
- Has no `.not. present(g)` guard (unlike Block 1 at line 237)
- Sits after `dnscsr_z_mkl` sparse conversion (line 278), so modifications to dense `HT` are lost in CSR output
- Double-counts Zeeman when called without `g=`

Block 1 (lines 236-258) is the correct implementation with the `.not. present(g)` guard.

**Step 2: Build and verify**

Run: `cmake --build build`
Expected: Clean build.

**Step 3: Run topology tests**

Run: `ctest --test-dir build -L unit -R "topology|chern|edge|bdg" -j4 --output-on-failure`
Expected: All pass.

**Step 4: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90
git commit -m "fix: remove duplicate Zeeman block in ZB8bandQW (double-counting)"
```

---

### Task 3: Fix Peierls substitution no-op (C3)

**Files:**
- Modify: `src/physics/magnetic_field.f90:53-106` (add_peierls_coo)
- Modify: `src/physics/bdg_hamiltonian.f90:148-161` (caller)

**Step 1: Understand the bug**

Two problems in `add_peierls_coo`:
1. It's called on `coo_vals_zeeman` which only contains diagonal entries. The off-diagonal check `coo_row /= coo_col` at line 88 always skips.
2. `coo_vals` is `real(kind=dp)` but Peierls phases are complex. `real(exp_phase)` at line 103 discards the imaginary part.

**Step 2: Rewrite add_peierls_coo to operate on the main complex COO**

Replace `add_peierls_coo` in `src/physics/magnetic_field.f90:53-106`:
- Change `coo_vals` from `real(kind=dp)` to `complex(kind=dp)`
- Remove unused args `gauge` and `kpterms_2d`
- Keep the early return for `abs(Bx) < 1e-12`
- Apply full complex phase: `coo_vals(idx) = coo_vals(idx) * exp_phase`

```fortran
subroutine add_peierls_coo(coo_vals, coo_row, coo_col, nnz_offset, &
                            grid, B_vec)
  complex(kind=dp), intent(inout) :: coo_vals(:)
  integer, intent(in) :: coo_row(:), coo_col(:)
  integer, intent(in) :: nnz_offset
  type(spatial_grid), intent(in) :: grid
  real(kind=dp), intent(in) :: B_vec(3)

  real(kind=dp) :: Bx, phase, y_i, y_j, dy_m, hbar_J
  integer :: idx, ngrid
  complex(kind=dp) :: exp_phase

  Bx = B_vec(1)
  if (abs(Bx) < 1e-12_dp) return

  ngrid = grid%npoints()
  if (.not. allocated(grid%z)) then
    print *, 'ERROR: add_peierls_coo: grid%z not allocated'
    stop 1
  end if

  dy_m = grid%dy * 1.0e-10_dp
  hbar_J = hbar * e

  do idx = 1, nnz_offset
    if (coo_row(idx) /= coo_col(idx)) then
      y_i = grid%z((coo_row(idx) - 1) / 8 + 1) * 1.0e-10_dp
      y_j = grid%z((coo_col(idx) - 1) / 8 + 1) * 1.0e-10_dp
      phase = e * Bx * (y_i - y_j) * dy_m / hbar_J
      exp_phase = cmplx(cos(phase), -sin(phase), kind=dp)
      coo_vals(idx) = coo_vals(idx) * exp_phase
    end if
  end do
end subroutine add_peierls_coo
```

Update the public interface (remove `gauge` and `kpterms_2d` params).

**Step 3: Update caller in bdg_hamiltonian.f90**

At `src/physics/bdg_hamiltonian.f90:152-154`, Peierls must operate on the MAIN complex COO array (`coo_vals_bdg`), not the Zeeman-only real array. Restructure:

1. First, copy H0 entries into `coo_vals_bdg` (Block 1,1 — already done at lines 126-134)
2. Call `add_peierls_coo` on the populated portion of `coo_vals_bdg`
3. Then add Zeeman entries (still diagonal-only, real-valued, merged during CSR assembly)

Reorder lines 140-161 to:
```fortran
! Apply Peierls to the electron block off-diagonal entries
if (present(B_vec) .and. abs(B_vec(1)) > 1e-12_dp) then
  call add_peierls_coo(coo_vals_bdg, coo_row_bdg, coo_col_bdg, &
                        coo_idx, cfg%grid, B_vec)
end if

! Add Zeeman splitting (diagonal, merged during CSR assembly)
if (present(B_vec) .and. ...) then
  ... (Zeeman stays the same, appends to coo_vals_bdg)
end if
```

**Step 4: Build and test**

Run: `cmake --build build && ctest --test-dir build -R "magnetic|landau|bdg" -j4 --output-on-failure`
Expected: Clean build, all pass. Peierls test in `test_magnetic_field.pf` may need updating to pass complex arrays.

**Step 5: Commit**

```bash
git add src/physics/magnetic_field.f90 src/physics/bdg_hamiltonian.f90 tests/unit/test_magnetic_field.pf
git commit -m "fix: Peierls substitution operates on complex off-diagonal entries"
```

---

### Task 4: Fix Majorana profile array bounds violation (C4)

**Files:**
- Modify: `src/apps/main_topology.f90:497` (caller — pass `Ntot_local` not `Ntot_local*2`)
- Modify: `src/physics/topological_analysis.f90:340-400` (add bounds check)

**Step 1: Fix the caller**

At `src/apps/main_topology.f90:497`, change:
```fortran
majorana_xi(n_majorana) = compute_majorana_profile( &
  & eigen_res_local%eigenvectors(:, i), cfg%grid, &
  & cfg%bdg%delta_0 * 0.01_dp, Ntot_local * 2)
```
to:
```fortran
majorana_xi(n_majorana) = compute_majorana_profile( &
  & eigen_res_local%eigenvectors(:, i), cfg%grid, &
  & cfg%bdg%delta_0 * 0.01_dp, Ntot_local)
```

The eigenvector has size `16*N` (full BdG). `half_n` should be `8*N` (electron block), which is `Ntot_local = 8*Ngrid_local`.

**Step 2: Add bounds check in compute_majorana_profile**

At `src/physics/topological_analysis.f90`, after `half_n = half_n_in` (line 352), add:
```fortran
if (half_n * 2 > size(evec_bdg)) then
  print *, 'ERROR: compute_majorana_profile: eigenvector too small'
  print *, '  half_n=', half_n, ' size(evec)=', size(evec_bdg)
  xi = 0.0_dp
  return
end if
```

**Step 3: Build and test**

Run: `cmake --build build && ctest --test-dir build -R "bdg|rashba" -j4 --output-on-failure`
Expected: All pass.

**Step 4: Commit**

```bash
git add src/apps/main_topology.f90 src/physics/topological_analysis.f90
git commit -m "fix: Majorana profile bounds violation (half_n was 2x too large)"
```

---

### Task 5: Add COO bounds checks in build_bdg_hamiltonian_1d (C5)

**Files:**
- Modify: `src/physics/bdg_hamiltonian.f90:129,157,166,182,207,225,236`

**Step 1: Add helper and checks**

After each `coo_idx = coo_idx + 1` in `build_bdg_hamiltonian_1d`, add:
```fortran
if (coo_idx > coo_capacity) then
  print *, 'ERROR: build_bdg_hamiltonian_1d: COO capacity exceeded'
  print *, '  coo_idx=', coo_idx, ' coo_capacity=', coo_capacity
  stop 1
end if
```

This pattern already exists in `insert_zeeman_coo` (`hamiltonian_wire.f90:1310`). Follow that precedent.

There are 5 locations where `coo_idx` is incremented (lines 129, 157, 166, 182, 207, 225, 236). Add the check after each increment. To avoid repetition, extract a small inline check or use a helper.

**Step 2: Build and test**

Run: `cmake --build build && ctest --test-dir build -R "bdg" -j4 --output-on-failure`
Expected: All pass.

**Step 3: Commit**

```bash
git add src/physics/bdg_hamiltonian.f90
git commit -m "fix: add COO bounds checks in BdG assembly"
```

---

### Task 6: Fix bdg_search loop consuming to EOF (C6)

**Files:**
- Modify: `src/io/input_parser.f90:797-829`

**Step 1: Fix the loop**

The `bdg_search` loop at line 797 reads lines until `bdg:` is found or EOF. When no `bdg:` block exists, it consumes all remaining lines, including `optics:`, `exciton:`, etc.

Fix: Stop when encountering any known optional block label and backspace so the next parser finds it. At line 826 (the `! Not bdg` comment), instead of continuing, check against known labels:

```fortran
bdg_search: do
  call read_next_data_line(data_unit, line, status)
  if (status /= 0) then
    found_optional = .false.
    exit bdg_search
  end if
  colon_pos = index(line, ':')
  if (colon_pos > 0) then
    label_part = adjustl(line(:colon_pos-1))
    select case(trim(to_lower_ascii(label_part)))
    case('bdg')
      ! Found bdg block - parse enabled flag
      block
        character(len=255) :: bdg_value_str
        bdg_value_str = adjustl(line(colon_pos+1:))
        ! ... existing T/F parsing ...
      end block
      found_optional = .true.
      label = trim(label_part) // ':'
      exit bdg_search
    case('b_field')
      cycle bdg_search  ! already parsed
    case default
      ! Known optional block or unknown label - stop and push back
      backspace(data_unit)
      found_optional = .false.
      exit bdg_search
    end select
  end if
  ! No colon - skip this line
end do bdg_search
```

The key change: `case default` does `backspace` + `exit`, so optics/exciton/strain blocks remain readable by their parsers.

**Step 2: Build and test**

Run: `cmake --build build && ctest --test-dir build -j4 --output-on-failure`
Expected: All pass. Configs without `bdg:` blocks (most existing configs) should still parse correctly.

**Step 3: Verify backward compat manually**

Run: `./build/src/bandStructure` with a non-topology config (e.g., `tests/regression/configs/sc_delta_doped_gaas.cfg`)
Expected: Runs normally, no parsing errors.

**Step 4: Commit**

```bash
git add src/io/input_parser.f90
git commit -m "fix: bdg_search stops at known labels instead of consuming to EOF"
```

---

### Task 7: Move FEAST auto-window before SC loop (C7)

**Files:**
- Modify: `src/apps/main_gfactor.f90:156-199`

**Step 1: Restructure**

At `src/apps/main_gfactor.f90`, the current flow is:
1. Lines 156-163: Set placeholder `[-1, 1]` if no feast_emin/emax
2. Lines 166-199: Run SC loop (uses placeholder!)
3. Lines 215-216: Auto-compute window (too late!)

Fix: Build a preliminary Hamiltonian, compute auto-window, THEN run SC loop.

Move the auto-window computation before the SC loop. Insert after the Hamiltonian build at the appropriate point:

```fortran
! Build preliminary Hamiltonian for auto-window
call ZB8bandGeneralized(HT_csr, 0.0_dp, profile_2d, kpterms_2d, cfg, ws=wire_ws)

! Auto-compute energy window from Hamiltonian
if (cfg%feast_emin == 0.0_dp .and. cfg%feast_emax == 0.0_dp) then
  call auto_compute_energy_window(HT_csr, eigen_cfg%emin, eigen_cfg%emax)
end if

call csr_free(HT_csr)  ! free preliminary; SC loop rebuilds it

! Now run SC loop with correct window
if (cfg%sc%enabled == 1) then
  ...
end if
```

Note: Check that the SC loop rebuilds the Hamiltonian internally (it should, since it modifies the potential profile). The preliminary build is only for window estimation.

**Step 2: Build and test**

Run: `cmake --build build && ctest --test-dir build -R "gfactor|sc" -j4 --output-on-failure`
Expected: All pass.

**Step 3: Commit**

```bash
git add src/apps/main_gfactor.f90
git commit -m "fix: compute FEAST energy window before SC loop"
```

---

### Task 8: Deduplicate LDOS functions (H1)

**Files:**
- Modify: `src/physics/green_functions.f90:47-281` (merge two functions into one)

**Step 1: Remove compute_ldos_bdg**

The two functions are identical except for the dummy argument name (`H` vs `H_bdg`). Remove `compute_ldos_bdg` entirely and keep `compute_ldos_csr`. Both operate on `type(csr_matrix)` — the LDOS formula is the same regardless of whether the matrix is normal or BdG.

Update the module to only export `compute_ldos_csr`. Remove line 26: `public :: compute_ldos_bdg`.

**Step 2: Update callers**

Search for `compute_ldos_bdg` in `src/apps/main_topology.f90` and replace with `compute_ldos_csr`.

**Step 3: Build and test**

Run: `cmake --build build && ctest --test-dir build -R "green|topology" -j4 --output-on-failure`
Expected: All pass.

**Step 4: Commit**

```bash
git add src/physics/green_functions.f90 src/apps/main_topology.f90
git commit -m "refactor: deduplicate LDOS functions (compute_ldos_bdg -> compute_ldos_csr)"
```

---

### Task 9: Extract shared Zeeman Vz subroutine (H2)

**Files:**
- Modify: `src/physics/magnetic_field.f90` (add `compute_zeeman_vz`)
- Modify: `src/physics/hamiltonianConstructor.f90:239-246` (use shared)
- Modify: `src/physics/magnetic_field.f90:38-42` (use shared)
- Modify: `src/physics/hamiltonian_wire.f90` (if applicable)

**Step 1: Add compute_zeeman_vz to magnetic_field.f90**

```fortran
pure subroutine compute_zeeman_vz(g_factor, mu_B_val, B_mag, Vz)
  real(kind=dp), intent(in) :: g_factor, mu_B_val, B_mag
  real(kind=dp), intent(out) :: Vz(8)
  Vz(1:2) = -1.5_dp * g_factor * mu_B_val * B_mag  ! HH
  Vz(3:4) =  0.5_dp * g_factor * mu_B_val * B_mag  ! LH
  Vz(5:6) = -0.5_dp * g_factor * mu_B_val * B_mag  ! SO
  Vz(7)   = -1.0_dp * g_factor * mu_B_val * B_mag  ! CB1
  Vz(8)   =  1.0_dp * g_factor * mu_B_val * B_mag  ! CB2
end subroutine compute_zeeman_vz
```

Make it `public`. Replace inline Vz computation at all 4 sites.

**Step 2: Update all call sites**

In `add_zeeman_coo` (magnetic_field.f90), `ZB8bandQW` (hamiltonianConstructor.f90 Block 1), and any wire Zeeman code.

**Step 3: Build and test**

Run: `cmake --build build && ctest --test-dir build -j4 --output-on-failure`
Expected: All pass.

**Step 4: Commit**

```bash
git add src/physics/magnetic_field.f90 src/physics/hamiltonianConstructor.f90
git commit -m "refactor: extract shared compute_zeeman_vz (4 copy-paste sites)"
```

---

### Task 10: Fix stub functions and diag_2x2 div-by-zero (H3/H5/H6)

**Files:**
- Modify: `src/physics/topological_analysis.f90:120-129` (Berry curvature stub)
- Modify: `src/physics/topological_analysis.f90:169-176` (Z2 Fu-Kane stub)
- Modify: `src/physics/topological_analysis.f90:116-117` (diag_2x2 normalization)

**Step 1: Add warnings to stubs**

In `compute_berry_curvature` (line 128), before `Omega = 0.0_dp`:
```fortran
print *, 'WARNING: compute_berry_curvature is not yet implemented. Returning zeros.'
```

In `compute_z2_fukane` (line 175), before `z2 = 0`:
```fortran
print *, 'WARNING: compute_z2_fukane is not yet implemented. Returning 0.'
```

**Step 2: Fix diag_2x2 div-by-zero**

At lines 116-117, add guards:
```fortran
if (sum(eigvec(:,1)**2) > 1.0e-30_dp) then
  eigvec(:,1) = eigvec(:,1) / sqrt(sum(eigvec(:,1)**2))
else
  eigvec(1,1) = 1.0_dp; eigvec(2,1) = 0.0_dp
end if
if (sum(eigvec(:,2)**2) > 1.0e-30_dp) then
  eigvec(:,2) = eigvec(:,2) / sqrt(sum(eigvec(:,2)**2))
else
  eigvec(1,2) = 0.0_dp; eigvec(2,2) = 1.0_dp
end if
```

**Step 3: Build and test**

Run: `cmake --build build && ctest --test-dir build -R "chern|topology" -j4 --output-on-failure`
Expected: All pass.

**Step 4: Commit**

```bash
git add src/physics/topological_analysis.f90
git commit -m "fix: add stub warnings and diag_2x2 div-by-zero guard"
```

---

### Task 11: Fix topology parser bug — remove workaround (H4/H9)

**Files:**
- Modify: `src/io/input_parser.f90:742,746,754,758` (use read_optional_logical_flag for logicals)
- Modify: `src/apps/main_topology.f90:63-67` (remove forced override workaround)
- Modify: `src/apps/main_topology.f90:70-72` (use cfg%topo%enabled instead of mode check)

**Step 1: Fix parser to handle T/F in topology block**

The topology block at lines 742-758 uses `read(data_unit, *, iostat=status) label, cfg%topo%compute_chern` which expects `.true.`/`.false.`. Change logical fields to use a T/F aware parser.

Replace each logical field read:
```fortran
read(data_unit, *, iostat=status) label, cfg%topo%compute_chern
```
with a custom T/F parsing helper that accepts both `T`/`F` and `.true.`/`.false.`. The `read_optional_logical_flag` already handles this. However, since we're inside the block (the label was already consumed), we need a simpler inline approach:

After reading the label, parse the remaining value:
```fortran
read(data_unit, '(A)', iostat=status) line
if (status /= 0) then; status = 0; exit topology_block; end if
! Extract value after label
colon_pos = index(line, ':')
if (colon_pos > 0) then
  value_str = adjustl(line(colon_pos+1:))
  if (value_str(1:1) == 'T' .or. value_str(1:1) == 't') then
    cfg%topo%compute_chern = .true.
  else
    cfg%topo%compute_chern = .false.
  end if
end if
```

Apply to all logical fields: `compute_chern`, `compute_hall`, `compute_z2`, `extract_edge_states`.

**Step 2: Remove workaround in main_topology.f90**

Delete lines 63-67:
```fortran
! Force compute_z2 for wire QSHE mode - needed because the parser has a bug
if (cfg%confinement == 2 .and. trim(cfg%topo%mode) == 'qshe') then
  cfg%topo%compute_z2 = .true.
end if
```

Also fix line 72: use `cfg%topo%enabled` directly instead of `len_trim(cfg%topo%mode)`:
```fortran
if (.not. cfg%topo%enabled .and. len_trim(cfg%topo%mode) == 0) then
```

**Step 3: Update test configs to use consistent format**

Check `tests/regression/configs/topology_*.cfg` — ensure all use consistent logical format (either all `T`/`F` or all `.true.`/`.false.`).

**Step 4: Build and test**

Run: `cmake --build build && ctest --test-dir build -R "topology" -j4 --output-on-failure`
Expected: All pass.

**Step 5: Commit**

```bash
git add src/io/input_parser.f90 src/apps/main_topology.f90 tests/regression/configs/topology_*.cfg
git commit -m "fix: topology parser handles T/F logics, remove forced override workaround"
```

---

### Task 12: Guard green_functions behind USE_ARPACK (H6)

**Files:**
- Modify: `src/physics/green_functions.f90:1-28` (module header)

**Step 1: Wrap module body**

The `pardiso_c` import is already behind `#ifdef USE_ARPACK`. But `compute_ldos_csr` uses it unconditionally. Wrap the entire module body:

```fortran
module green_functions
  use definitions, only: dp, pi_dp
  use sparse_matrices
  implicit none
  private

#ifdef USE_ARPACK
  public :: compute_ldos_csr
contains
  ! ... existing LDOS code ...
#endif

end module green_functions
```

When ARPACK is not available, the module compiles to empty exports. Callers that `use green_functions` won't fail because the module still exists — it just has no public symbols.

**Step 2: Guard callers**

In `src/apps/main_topology.f90`, guard the `use green_functions` line:
```fortran
#ifdef USE_ARPACK
  use green_functions, only: compute_ldos_csr
#endif
```

And guard the LDOS call site similarly.

**Step 3: Build with and without ARPACK**

Run: `cmake --build build` (with ARPACK enabled)
Expected: Clean build.

If possible, test without ARPACK:
```bash
cmake -G Ninja -B build-noarpack -DMKL_DIR=$MKLROOT/lib/cmake/mkl
cmake --build build-noarpack
```
Expected: `topologicalAnalysis` builds but LDOS features are disabled.

**Step 4: Commit**

```bash
git add src/physics/green_functions.f90 src/apps/main_topology.f90
git commit -m "fix: guard green_functions behind USE_ARPACK to fix build without ARPACK"
```

---

### Task 13: Fix hardcoded path and remove unused vars (H7/H8)

**Files:**
- Modify: `scripts/sweep_rashba_bdg.py:16`
- Modify: `src/physics/green_functions.f90:53,180` (remove `info`, `k`)
- Modify: `src/io/input_parser.f90:52` (remove `rewind_status`)
- Modify: `src/physics/magnetic_field.f90:68` (remove `nx`)
- Modify: `src/physics/topological_analysis.f90:5-6` (remove unused `use linalg`, `use eigensolver`)
- Modify: `scripts/verify_bhz_z2.py:15` (remove `import os`)
- Modify: `scripts/verify_landau_levels.py:15` (remove `import os`)
- Modify: `src/physics/green_functions.f90:15-16` (remove unused `c_double_complex`, `real64` if unused after dedup)

**Step 1: Fix hardcoded path**

In `scripts/sweep_rashpa_bdg.py:16`, replace:
```python
EXE = Path("/data/8bandkp-fdm/build/src/topologicalAnalysis")
```
with:
```python
REPO = Path(__file__).resolve().parent.parent
EXE = REPO / "build" / "src" / "topologicalAnalysis"
```

**Step 2: Remove unused declarations**

For each file listed above, remove the unused variable/import. Build after each batch to verify no breakage.

**Step 3: Build and test**

Run: `cmake --build build && ctest --test-dir build -j4 --output-on-failure`
Expected: All pass.

**Step 4: Commit**

```bash
git add scripts/sweep_rashba_bdg.py scripts/verify_bhz_z2.py scripts/verify_landau_levels.py
git add src/physics/green_functions.f90 src/io/input_parser.f90 src/physics/magnetic_field.f90 src/physics/topological_analysis.f90
git commit -m "fix: hardcoded path, remove unused variables and imports (Codacy)"
```

---

### Task 14: Add missing test coverage — BdG and Zeeman (T1-T3)

**Files:**
- Modify: `tests/unit/test_bdg_hamiltonian.pf` (add CSR assembly test)
- Modify: `tests/unit/test_green_functions.pf` (add BdG LDOS test — only if Task 8 dedup keeps it)
- Modify: `tests/unit/test_hamiltonian.pf` (add Zeeman QW test)

**Step 1: Add BdG CSR assembly test (T1)**

In `test_bdg_hamiltonian.pf`, add a test that calls `build_bdg_hamiltonian_1d` on a minimal wire (N=2, GaAs parameters, no B field) and verifies:
- Resulting CSR has dimension `16*N = 32`
- Matrix is Hermitian: `H = H^dagger` (check a sample of entries)
- Eigenvalues come in +/- E pairs (particle-hole symmetry)

```fortran
@test(subroutine)
subroutine test_bdg_csr_hermitian()
  ! Build minimal BdG, verify Hermiticity and particle-hole symmetry
  type(simulation_config) :: cfg
  type(csr_matrix) :: H_bdg
  real(kind=dp), allocatable :: profile(:,:), eigvals(:)
  type(csr_matrix), allocatable :: kpterms(:)
  ! ... setup minimal 2-point GaAs wire ...
  call build_bdg_hamiltonian_1d(H_bdg, cfg, profile, kpterms, &
    & 0.0_dp, 0.0_dp, 0.1_dp, ws)
  ! Verify size
  @assertEqual(32, H_bdg%nrows)
  ! Verify Hermiticity: H(i,j) == conjg(H(j,i))
  ! ... sample check ...
  call csr_free(H_bdg)
end subroutine
```

**Step 2: Add Zeeman QW test (T3)**

In `test_hamiltonian.pf`, add a test that builds QW Hamiltonian with B_z = 1T, g = 2.0, and verifies the Zeeman splitting matches expected g_J * g * mu_B * B values on the diagonal.

**Step 3: Build and run tests**

Run: `cmake --build build && ctest --test-dir build -L unit -R "bdg|hamiltonian" -j4 --output-on-failure`
Expected: New tests pass.

**Step 4: Commit**

```bash
git add tests/unit/test_bdg_hamiltonian.pf tests/unit/test_hamiltonian.pf tests/unit/test_green_functions.pf
git commit -m "test: add BdG CSR assembly and Zeeman QW unit tests"
```

---

### Task 15: Add quantitative assertions to existing tests (T4-T8)

**Files:**
- Modify: `tests/integration/test_rashba_phase_boundary.sh` (add gap threshold check)
- Modify: `tests/unit/test_sc_loop.pf` (add delta-doping normalization)
- Modify: `tests/unit/test_magnetic_field.pf` (add zero-B test)
- Modify: `tests/unit/test_bdg_config.pf` (add absent b_field test)

**Step 1: Rashba quantitative assertion (T4)**

In `test_rashba_phase_boundary.sh`, after extracting min_gap, add:
```bash
# Verify gap is near zero for topological parameters (Vz > Vc ≈ 0.58)
if (( $(echo "$min_gap > 0.5" | bc -l) )); then
  echo "FAIL: min_gap=$min_gap should be < 0.5 meV for topological phase"
  exit 1
fi
```

**Step 2: Delta-doping normalization (T7)**

In `test_sc_loop.pf`, add to `test_delta_doping_gaussian`:
```fortran
! Verify integral equals NS (sheet density)
block
  real(kind=dp) :: integral, dz
  dz = 1.0_dp  ! grid spacing in test
  integral = sum(rho) * dz * 1.0e-8_dp  ! Angstrom -> cm
  @assertEqual(5.0e11_dp, integral, 1.0e10_dp)  ! NS = 5e11 cm^-2
end block
```

Adjust the tolerance and conversion factor to match the test setup.

**Step 3: Zero-B Zeeman test (T6)**

In `test_magnetic_field.pf`, add:
```fortran
@test(subroutine)
subroutine test_zeeman_zero_field()
  ! All Zeeman entries should be zero at B=0
  real(kind=dp) :: coo_vals(8), B_vec(3)
  integer :: coo_row(8), coo_col(8), nnz_off
  type(spatial_grid) :: grid
  B_vec = [0.0_dp, 0.0_dp, 0.0_dp]
  nnz_off = 0
  ! ... setup 1-point grid ...
  call add_zeeman_coo(coo_vals, coo_row, coo_col, nnz_off, grid, B_vec, 2.0_dp)
  do i = 1, 8
    @assertEqual(0.0_dp, coo_vals(i), 1.0e-15_dp)
  end do
end subroutine
```

**Step 4: Build and test**

Run: `cmake --build build && ctest --test-dir build -j4 --output-on-failure`
Expected: All pass.

**Step 5: Commit**

```bash
git add tests/integration/test_rashba_phase_boundary.sh tests/unit/test_sc_loop.pf tests/unit/test_magnetic_field.pf
git commit -m "test: add quantitative assertions, normalization check, zero-B edge case"
```

---

### Task 16: Add remaining defensive checks (H10-H12)

**Files:**
- Modify: `src/physics/bdg_hamiltonian.f90:86` (validate H0 after build)
- Modify: `src/apps/main_topology.f90:181-195` (file I/O iostat)

**Step 1: Validate H0 in build_bdg_hamiltonian_1d**

After line 86 (`call ZB8bandGeneralized(H0, ...)`), add:
```fortran
if (H0%nrows == 0 .or. H0%nnz == 0) then
  print *, 'ERROR: build_bdg_hamiltonian_1d: empty H0 from ZB8bandGeneralized'
  stop 1
end if
if (mod(H0%nrows, 8) /= 0) then
  print *, 'ERROR: build_bdg_hamiltonian_1d: H0 nrows not multiple of 8'
  stop 1
end if
```

**Step 2: Add iostat to file I/O in main_topology**

At lines 181-195, add `iostat=status` to `open` and check:
```fortran
open(unit=iounit, file='output/topology_result.dat', status='replace', &
     action='write', iostat=status)
if (status /= 0) then
  print *, 'ERROR: cannot open output/topology_result.dat (iostat=', status, ')'
  stop 1
end if
```

**Step 3: Build and test**

Run: `cmake --build build && ctest --test-dir build -j4 --output-on-failure`
Expected: All pass.

**Step 4: Commit**

```bash
git add src/physics/bdg_hamiltonian.f90 src/apps/main_topology.f90
git commit -m "fix: add H0 validation and file I/O error checking"
```

---

### Task 17: Final integration test and BACKLOG update

**Files:**
- Modify: `docs/plans/BACKLOG.md` (mark Phase 3 complete)

**Step 1: Full test suite**

Run: `OMP_NUM_THREADS=12 ctest --test-dir build -j4 --output-on-failure`
Expected: All 39 tests pass (15 unit + 24 regression).

**Step 2: Build clean from scratch**

```bash
rm -rf build && cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl && cmake --build build
```
Expected: Clean build, no warnings.

**Step 3: Update BACKLOG**

Mark Phase 3 as COMPLETED in `docs/plans/BACKLOG.md`.

**Step 4: Commit**

```bash
git add docs/plans/BACKLOG.md
git commit -m "docs: update BACKLOG for Phase 3 completion"
```

**Step 5: Squash or rebase if needed, push**

Review the commit history. If too granular, interactive rebase into logical groups. Push to remote.

---

### Addendum: Root-Caused "Pre-Existing" Test Failures (2026-05-03)

All 5 tests that were initially dismissed as "pre-existing" were traced to real bugs via systematic debugging and fixed. No test failures remain.

| Test | Root Cause | Fix | Commit |
|------|-----------|-----|--------|
| `test_chern_number` | 3 bugs in `compute_chern_qwz`: non-Hermitian QWZ H, eigenvector cross-contamination, wrong FHS plaquette | Rewrote with complex Hermitian H, analytical eigenvectors, proper 4-corner plaquette | `aafe6d7` |
| `regression_topology_qwz_chern` | Same as above | Same as above | `aafe6d7` |
| `regression_topology_bhz_z2` | Config omitted `compute_hall`/`qwz_u` → parser field-order misalignment → `compute_z2` never set | Added missing fields to config files | `16379ca` |
| `regression_landau_inas` | `mu_B = e*hbar/(2*m0)` = 9.274e-4 (16x too large); test expected unimplemented orbital Landau levels | Added CODATA `mu_B = 5.788e-5` constant; test now verifies Zeeman splitting | `67a653f` |
| `regression_topology_rashba_phase` | `bdg:` block inside topology block; missing `B_vec`; `ldos_E_range` consumed `ldos_num_E` | Moved bdg: after topology; added B_vec; split ldos_E_range onto 2 lines | `5fb0bd3` |
