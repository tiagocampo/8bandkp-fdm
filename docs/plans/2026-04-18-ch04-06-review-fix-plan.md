# Ch04-06 Review & Fix Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix the wire Hamiltonian missing-const bug, the strain sign error, and update Ch04-06 documentation with correct data.

**Architecture:** Two independent code fixes (wire const factor, strain sign) followed by documentation updates. The wire fix is the largest: fold `const = hbar²/2m₀` into the material parameter profiles in `confinementInitialization_2d`. The strain fix is a one-line formula change in both QW and wire paths. Documentation updates regenerate figures and tables from corrected code.

**Tech Stack:** Fortran 90, pFUnit tests, CMake/Ninja build, Python plotting scripts

---

### Task 1: Fix wire Hamiltonian missing `const` factor

**Root cause:** `confinementInitialization_2d` builds kpterms_2d operators from material profiles (gamma1, gamma2, gamma3, A) and FD operators (Laplacian, gradient). The profiles are dimensionless but the FD operators are in Å⁻² (Laplacian) or Å⁻¹ (gradient). The product gives Å⁻² or Å⁻¹ — NOT eV. The code is missing `const = hbar²/(2m₀) = 3.80998 eV·Å²` that converts to eV. The QW code includes const explicitly (`kpterms(i,j,1) = -gamma1 * const * D2(j,j)`) but the wire code does not.

**Fix:** Scale all gamma and A material profiles by `const` before building kpterms_2d. Leave `P` unscaled (it already includes const: `P = sqrt(EP * const)`).

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90:496-533` (profile construction)
- Modify: `src/physics/hamiltonianConstructor.f90:542-616` (kpterms_2d building)
- Test: `tests/unit/test_hamiltonian_2d.pf`

**Step 1: Write the failing test**

Add a test to `tests/unit/test_hamiltonian_2d.pf` that verifies the wire Hamiltonian diagonal has the correct scale. For a uniform GaAs wire at kz=0, the CB block diagonal at an interior point should be:

```
EC + A*const*(2/dx² + 2/dy²)
= 0.719 + 14.93 * 3.81 * (2/9 + 2/9)
= 0.719 + 25.3
≈ 26.0 eV
```

```fortran
@test
subroutine test_wire_hamiltonian_diagonal_units()
  ! Verify that the wire Hamiltonian diagonal has the correct energy scale.
  ! For a uniform GaAs wire, the CB block diagonal should include
  ! const * A * Laplacian_diag + EC, where const = hbar^2/(2*m_0).
  ! Without const, the kinetic term is ~3.8x too small.

  integer, parameter :: nx = 5, ny = 5
  real(kind=dp), parameter :: dx = 3.0_dp, dy = 3.0_dp
  integer :: ngrid, Ntot, ij, mid, band7_diag_idx
  real(kind=dp) :: expected_cb_diag_min, actual_cb_diag, const_val

  type(spatial_grid) :: grid
  type(paramStruct) :: params(1)
  character(len=255) :: materials(1)
  type(region_spec) :: regions(1)
  real(kind=dp), allocatable :: profile_2d(:,:)
  type(csr_matrix), allocatable :: kpterms_2d(:)
  type(csr_matrix) :: HT_csr
  type(simulation_config) :: cfg
  type(wire_coo_cache) :: coo_cache

  const_val = hbar2O2m0  ! 3.80998 eV*Ang^2

  materials(1) = "GaAs"
  call paramDatabase(materials, 1, params)

  ! Set up 5x5 rectangular grid, all GaAs
  grid%ndim = 2
  grid%nx = nx; grid%ny = ny
  grid%dx = dx; grid%dy = dy
  allocate(grid%x(nx), grid%z(ny))
  do ij = 1, nx; grid%x(ij) = (ij-1)*dx; end do
  do ij = 1, ny; grid%z(ij) = (ij-1)*dy; end do
  ngrid = nx * ny
  allocate(grid%material_id(ngrid))
  grid%material_id = 1

  regions(1)%material = "GaAs"
  regions(1)%fraction = 100.0_dp

  ! Build kpterms_2d
  call confinementInitialization_2d(grid, params, regions, &
    profile_2d, kpterms_2d, FDorder=2)

  ! Build Hamiltonian at kz=0
  call init_simulation_config(cfg)
  cfg%grid = grid
  cfg%params = params
  cfg%regions = regions
  cfg%FDorder = 2

  call ZB8bandGeneralized(HT_csr, 0.0_dp, profile_2d, kpterms_2d, cfg, coo_cache)

  ! Check an interior CB diagonal element (band 7, center grid point)
  ! Center point index (3,3) in column-major: (3-1)*5 + 3 = 13
  ! Band 7 offset: 6 * ngrid = 6 * 25 = 150
  ! Global index: 150 + 13 = 163
  band7_diag_idx = 6 * ngrid + 13

  ! The CB diagonal should be EC + A*const*Laplacian_diag
  ! For interior point with 2nd order FD: Laplacian_diag = -2/dx^2 - 2/dy^2
  ! With const: A*const*|Laplacian_diag| = 14.93*3.81*0.444 ≈ 25.3
  ! Plus EC = 0.719 → total ≈ 26.0 eV
  expected_cb_diag_min = params(1)%EC + 10.0_dp  ! conservative lower bound

  ! Extract diagonal element from CSR
  actual_cb_diag = 0.0_dp
  do ij = HT_csr%rowptr(band7_diag_idx), HT_csr%rowptr(band7_diag_idx + 1) - 1
    if (HT_csr%colind(ij) == band7_diag_idx) then
      actual_cb_diag = real(HT_csr%values(ij), kind=dp)
      exit
    end if
  end do

  @assertTrue(actual_cb_diag > expected_cb_diag_min, &
    message="Wire CB diagonal includes const factor (got " // &
    trim(adjustl(str_from_dp(actual_cb_diag))) // " expected > " // &
    trim(adjustl(str_from_dp(expected_cb_diag_min))) // ")")

  call csr_free(HT_csr)
  call wire_coo_cache_free(coo_cache)
  ! cleanup kpterms_2d
  do ij = 1, size(kpterms_2d)
    call csr_free(kpterms_2d(ij))
  end do
  deallocate(kpterms_2d, profile_2d)
  deallocate(grid%x, grid%z, grid%material_id)
end subroutine test_wire_hamiltonian_diagonal_units
```

Note: The helper `str_from_dp` and `init_simulation_config` may not exist. If not, use a simpler assertion that just checks `actual_cb_diag > 10.0`.

**Step 2: Run test to verify it fails**

Run: `ctest --test-dir build -L unit -R test_hamiltonian_2d -V`
Expected: FAIL — the CB diagonal will be ~6.6 eV (without const) instead of ~26 eV.

**Step 3: Implement the fix**

In `src/physics/hamiltonianConstructor.f90`, modify the profile construction (around line 527-533) to multiply gamma and A profiles by `const`:

```fortran
! Material parameter profiles (scaled by const = hbar^2/(2m_0))
! The gamma and A parameters are dimensionless; multiplying by const
! converts the kpterms operators to energy units (eV).
! P is NOT scaled because it already includes const: P = sqrt(EP*const).
prof_gamma1(ij) = params(mid)%gamma1 * params(mid)%const
prof_gamma2(ij) = params(mid)%gamma2 * params(mid)%const
prof_gamma3(ij) = params(mid)%gamma3 * params(mid)%const
prof_P(ij)      = params(mid)%P               ! already includes const
prof_A(ij)      = params(mid)%A * params(mid)%const
prof_gm12g2(ij) = (params(mid)%gamma1 - 2.0_dp * params(mid)%gamma2) * params(mid)%const
prof_gp12g2(ij) = (params(mid)%gamma1 + 2.0_dp * params(mid)%gamma2) * params(mid)%const
```

This affects kpterms_2d(1-3, 5, 7-11, 14-15) — all terms that use gamma or A profiles. kpterms_2d(4, 6, 12-13) use P and are NOT changed (P already has const).

**Step 4: Run test to verify it passes**

Run: `ctest --test-dir build -L unit -R test_hamiltonian_2d -V`
Expected: PASS

**Step 5: Run all tests to check for regressions**

Run: `ctest --test-dir build -V`
Expected: All 13 tests pass (8 unit + 5 regression).

**Step 6: Run wire simulation and verify eigenvalues**

```bash
cp tests/regression/configs/wire_gaas_rectangle.cfg input.cfg
./build/src/bandStructure
head -3 output/eigenvalues.dat
```

Expected: CB1 eigenvalue at k=0 should be above EC = 0.719 eV (plus confinement). For a 63 Å box with 21×21 grid, CB1 ≈ 0.72-0.80 eV range.

**Step 7: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90 tests/unit/test_hamiltonian_2d.pf
git commit -m "fix: add missing hbar²/2m₀ factor in wire Hamiltonian kinetic terms"
```

---

### Task 2: Fix strain sign convention

**Root cause:** `compute_strain_qw` (line 148) and `compute_strain_wire` (line 259) use `eps = (a_mat - a_ref) / a_ref`, which gives the opposite sign from the standard pseudomorphic strain formula `eps = (a_sub - a_layer) / a_layer`. For InAs on GaAs (compressive), the code gives positive strain instead of negative.

**Fix:** Change the formula to `(a0_ref - a0_mat) / a0_mat`.

**Files:**
- Modify: `src/physics/strain_solver.f90:148` (QW path)
- Modify: `src/physics/strain_solver.f90:259` (wire path)
- Modify: `tests/unit/test_strain_solver.pf:58` (update expected value)
- Test: `tests/unit/test_strain_solver.pf`

**Step 1: Update the test expected value**

In `tests/unit/test_strain_solver.pf`, line 58, change:

```fortran
eps_0_expected = (params(2)%a0 - a0_ref) / a0_ref
```
to:
```fortran
eps_0_expected = (a0_ref - params(2)%a0) / params(2)%a0
```

This changes the expected value from `(6.0583 - 5.65325) / 5.65325 = +0.0716` (tensile, wrong for InAs on GaAs) to `(5.65325 - 6.0583) / 6.0583 = -0.0669` (compressive, correct).

**Step 2: Run test to verify it fails**

Run: `ctest --test-dir build -L unit -R test_strain_solver -V`
Expected: FAIL — the code still produces the old formula.

**Step 3: Fix the QW strain formula**

In `src/physics/strain_solver.f90`, line 148, change:

```fortran
eps_0 = (a0_mat - a0_ref) / a0_ref
```
to:
```fortran
eps_0 = (a0_ref - a0_mat) / a0_mat
```

**Step 4: Fix the wire strain formula**

In `src/physics/strain_solver.f90`, line 259, change:

```fortran
eps_xx_fixed(ij) = (params(mid)%a0 - a0_ref) / a0_ref
```
to:
```fortran
eps_xx_fixed(ij) = (a0_ref - params(mid)%a0) / params(mid)%a0
```

**Step 5: Run test to verify it passes**

Run: `ctest --test-dir build -L unit -R test_strain_solver -V`
Expected: PASS

**Step 6: Add a physics sign verification test**

Add to `tests/unit/test_strain_solver.pf`:

```fortran
@test
subroutine test_strain_sign_inas_on_gaas()
  ! InAs on GaAs substrate: InAs has LARGER lattice constant.
  ! Pseudomorphic strain should be COMPRESSIVE (negative eps_xx).
  ! eps_xx = (a_sub - a_layer) / a_layer = (5.65325 - 6.0583) / 6.0583 < 0
  type(paramStruct) :: params(2)
  character(len=255) :: materials(2)
  type(spatial_grid) :: grid
  type(strain_config) :: strain_cfg
  type(strain_result) :: strain_out
  real(kind=dp) :: a0_ref
  integer, parameter :: ny = 4

  materials(1) = "GaAs"
  materials(2) = "InAs"
  call paramDatabase(materials, 2, params)

  a0_ref = params(1)%a0  ! GaAs = 5.65325

  grid%ndim = 1; grid%nx = 1; grid%ny = ny
  grid%dx = 0.0_dp; grid%dy = 1.0_dp
  allocate(grid%z(ny), grid%material_id(ny))
  grid%material_id = [1, 2, 2, 1]
  strain_cfg%enabled = .true.

  call compute_strain(grid, params, grid%material_id, strain_cfg, a0_ref, strain_out)

  ! InAs points (2,3): eps_xx should be NEGATIVE (compressive)
  @assertTrue(strain_out%eps_xx(2) < 0.0_dp, message="InAs on GaAs: compressive eps_xx")
  @assertTrue(strain_out%eps_xx(3) < 0.0_dp, message="InAs on GaAs: compressive eps_xx")
  ! GaAs points (1,4): eps_xx should be zero (reference material)
  @assertEqual(0.0_dp, strain_out%eps_xx(1), tolerance=1.0e-15_dp)
  @assertEqual(0.0_dp, strain_out%eps_xx(4), tolerance=1.0e-15_dp)

  call strain_result_free(strain_out)
  deallocate(grid%z, grid%material_id)
end subroutine test_strain_sign_inas_on_gaas
```

**Step 7: Run all tests**

Run: `ctest --test-dir build -V`
Expected: All tests pass.

**Step 8: Commit**

```bash
git add src/physics/strain_solver.f90 tests/unit/test_strain_solver.pf
git commit -m "fix: correct strain sign to standard (a_sub-a_layer)/a_layer convention"
```

---

### Task 3: Update Ch04 strain documentation

**Files:**
- Modify: `docs/lecture/04-strain.md:69-80` (sign convention note)
- Modify: `docs/lecture/04-strain.md:72-74` (code excerpt)

**Step 1: Update the code excerpt and note**

At line 72, update the code to reflect the fixed formula:

```fortran
eps_0 = (a0_ref - a0_mat) / a0_mat
strain_out%eps_xx(ij) = eps_0
strain_out%eps_yy(ij) = eps_0
```

At lines 77-80, update the note:

```
Note the sign: if $a_0^{\text{mat}} > a_0^{\text{ref}}$ (e.g., InAs on GaAs substrate),
then $\varepsilon_0 < 0$ — the layer is compressed.  The formula
$(a_{\text{ref}} - a_0^{\text{mat}})/a_0^{\text{mat}}$ is equivalent to the
textbook expression $(a_{\text{sub}} - a_{\text{layer}})/a_{\text{layer}}$
since $a_0^{\text{ref}} = a_{\text{sub}}$ and $a_0^{\text{mat}} = a_{\text{layer}}$.
```

**Step 2: Verify Example A numbers still match**

The examples (lines 520-601) already use the textbook formula `(a_sub - a_layer)/a_layer`, so they remain correct. No numerical changes needed.

**Step 3: Commit**

```bash
git add docs/lecture/04-strain.md
git commit -m "docs: update Ch04 strain formula to match corrected code"
```

---

### Task 4: Generate real wire data and update Ch06

**Depends on:** Task 1 (wire Hamiltonian fix)

**Files:**
- Modify: `docs/lecture/06-optical-properties.md:308-433` (Section 6.6)

**Step 1: Run wire band structure**

```bash
cp tests/regression/configs/wire_gaas_rectangle.cfg input.cfg
./build/src/bandStructure
```

Extract eigenvalues at k=0 from `output/eigenvalues.dat`.

**Step 2: Run wire g-factor / optical transitions**

Create a config for the wire gfactor calculation (or use existing wire config with gfactorCalculation). Extract optical transition data from `output/optical_transitions.dat`.

**Step 3: Generate wire density figure**

Use the eigenfunction output to generate `wire_density_2d.png` with a plotting script.

**Step 4: Update Ch06 Section 6.6**

Replace the fabricated Table 6.3 data with actual computed values. Update:
- Section 6.6.1: wire parameters (grid size, materials) to match the config used
- Section 6.6.2: Table 6.3 with real oscillator strengths
- Section 6.6.3: Figure 6.2 caption with correct energy values
- Remove or flag the circular wire claim (only rectangular configs exist)

**Step 5: Commit**

```bash
git add docs/lecture/06-optical-properties.md docs/figures/wire_*.png
git commit -m "docs: replace fabricated wire data with actual simulation results"
```

---

### Task 5: Verify QW optical data in Ch06

**Depends on:** Nothing (independent of Tasks 1-2)

**Files:**
- Possibly modify: `docs/lecture/06-optical-properties.md:403-423` (QW optical section)

**Step 1: Run QW gfactor for optical transitions**

```bash
cp tests/regression/configs/qw_gaas_algaas_optics.cfg input.cfg
# Add gfactor-specific fields (whichBand, bandIdx, waveVector: k0, waveVectorStep: 0)
./build/src/gfactorCalculation
cat output/optical_transitions.dat
```

**Step 2: Compare with Ch06 values**

Ch06 claims CB1→VB2 at ~1.561 eV, |px|²≈|py|²≈51.6. Verify these match the output.

**Step 3: Update if needed**

If values don't match, update the documentation.

**Step 4: Commit (if changes needed)**

```bash
git add docs/lecture/06-optical-properties.md
git commit -m "docs: update Ch06 QW optical data from actual code output"
```

---

## Dependencies

- Task 1 (wire const fix) → Task 4 (wire data generation)
- Task 2 (strain sign) → Task 3 (Ch04 docs)
- Task 5 (QW optical verify) is independent

## Files Modified Summary

| File | Tasks | Changes |
|---|---|---|
| `src/physics/hamiltonianConstructor.f90` | 1 | Scale gamma/A profiles by const |
| `src/physics/strain_solver.f90` | 2 | Fix strain formula sign (2 locations) |
| `tests/unit/test_hamiltonian_2d.pf` | 1 | Add diagonal units test |
| `tests/unit/test_strain_solver.pf` | 2 | Fix expected value, add sign test |
| `docs/lecture/04-strain.md` | 3 | Update code excerpt and note |
| `docs/lecture/06-optical-properties.md` | 4, 5 | Replace wire data, verify QW data |
| `docs/figures/wire_*.png` | 4 | Regenerate from actual data |
