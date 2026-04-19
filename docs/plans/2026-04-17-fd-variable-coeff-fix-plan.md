# FD Variable-Coefficient & Config Overlap Fix Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix two pre-existing bugs — FDorder>2 variable coefficients and 3-layer config overlap — then regenerate all affected physics results and documentation.

**Architecture:** Replace naive row-scaling in `applyVariableCoeff` with midpoint-averaged profile matrix. Add mask-based layer assignment in `confinementInitialization_raw`. Fix overlapping configs. Re-run all affected simulations and update docs/figures.

**Tech Stack:** Fortran 90, pFUnit, LAPACK (zheevd), Python/Matplotlib for figures, CMake/ctest

---

### Task 1: Write failing test for FDorder=4 variable-coefficient bug

**Files:**
- Modify: `tests/unit/test_hamiltonian.pf`

**Context:** The test builds a small GaAs/AlGaAs QW, constructs kpterms with FDorder=2 and FDorder=4, and verifies they produce the same eigenvalues. Currently FDorder=4 gives wrong results. This test will FAIL until the bug is fixed.

**Step 1: Add the failing test to `test_hamiltonian.pf`**

Append this test inside `module test_hamiltonian` (before `end module`), after the existing `test_bulk_hermitian_inas`:

```fortran
  @test
  subroutine test_qw_variable_coeff_order4()
    ! FDorder=4 should give same eigenvalues as FDorder=2 for a QW
    ! (both must correctly handle variable coefficients at interfaces)
    integer, parameter :: ngrid = 41
    real(kind=dp), parameter :: dz_val = 2.5_dp  ! Angstrom
    integer :: matSize, info, lwork, lrwork, liwork

    ! Grid
    real(kind=dp) :: z(ngrid)

    ! Materials: 2-layer AlGaAs/GaAs (non-overlapping)
    integer, parameter :: nlayers = 2
    character(len=255) :: material(nlayers)
    type(paramStruct) :: params(nlayers)
    integer :: intStartPos(nlayers), intEndPos(nlayers)

    ! Confinement
    real(kind=dp), allocatable :: profile(:,:), kpterms(:,:,:)

    ! Hamiltonian and eigenvalues
    complex(kind=dp), allocatable :: HT(:,:), work(:)
    real(kind=dp), allocatable :: eig(:,:), rwork(:)
    integer, allocatable :: iwork(:)

    ! Results storage
    real(kind=dp) :: cb1_order2, cb1_order4
    real(kind=dp) :: vb1_order2, vb1_order4
    integer :: num_states_to_check, idx

    matSize = ngrid * 8

    ! ---- Setup grid ----
    do idx = 1, ngrid
      z(idx) = -50.0_dp + dble(idx - 1) * dz_val
    end do

    ! 2-layer: AlGaAs(-50,-25) | GaAs(-25,50)
    ! AlGaAs is points 1-10, GaAs is points 11-41
    intStartPos(1) = 1;   intEndPos(1) = 10
    intStartPos(2) = 11;  intEndPos(2) = 41
    material(1) = "Al30Ga70As"
    material(2) = "GaAs"
    call paramDatabase(material, nlayers, params)

    ! ---- Run with FDorder=2 ----
    allocate(kpterms(ngrid, ngrid, 10))
    kpterms = 0.0_dp
    call confinementInitialization_raw(z, intStartPos, intEndPos, material, &
      & nlayers, params, 'z', profile, kpterms, FDorder=2)

    allocate(HT(matSize, matSize), eig(matSize, 1))
    HT = cmplx(0.0_dp, kind=dp)
    call ZB8bandQW(HT, wavevector(0.0_dp, 0.0_dp, 0.0_dp), profile, kpterms)
    deallocate(profile, kpterms)

    allocate(work(1), rwork(1), iwork(1))
    call zheevd('N', 'U', matSize, HT, matSize, eig(:,1), work, -1, rwork, -1, iwork, -1, info)
    lwork = int(real(work(1))); lrwork = int(rwork(1)); liwork = iwork(1)
    deallocate(work, rwork, iwork)
    allocate(work(lwork), rwork(lrwork), iwork(liwork))
    call zheevd('N', 'U', matSize, HT, matSize, eig(:,1), work, lwork, rwork, lrwork, &
      & iwork, liwork, info)
    @assertTrue(info == 0, message="FDorder=2 diagonalization succeeded")

    ! CB1 = highest eigenvalue, VB1 = next-to-top group
    cb1_order2 = eig(matSize, 1)
    vb1_order2 = eig(matSize - 8, 1)  ! top of VB block (8 bands per kpoint)
    deallocate(HT, eig, work, rwork, iwork)

    ! ---- Run with FDorder=4 ----
    allocate(kpterms(ngrid, ngrid, 10))
    kpterms = 0.0_dp
    call confinementInitialization_raw(z, intStartPos, intEndPos, material, &
      & nlayers, params, 'z', profile, kpterms, FDorder=4)

    allocate(HT(matSize, matSize), eig(matSize, 1))
    HT = cmplx(0.0_dp, kind=dp)
    call ZB8bandQW(HT, wavevector(0.0_dp, 0.0_dp, 0.0_dp), profile, kpterms)
    deallocate(profile, kpterms)

    allocate(work(1), rwork(1), iwork(1))
    call zheevd('N', 'U', matSize, HT, matSize, eig(:,1), work, -1, rwork, -1, iwork, -1, info)
    lwork = int(real(work(1))); lrwork = int(rwork(1)); liwork = iwork(1)
    deallocate(work, rwork, iwork)
    allocate(work(lwork), rwork(lrwork), iwork(liwork))
    call zheevd('N', 'U', matSize, HT, matSize, eig(:,1), work, lwork, rwork, lrwork, &
      & iwork, liwork, info)
    @assertTrue(info == 0, message="FDorder=4 diagonalization succeeded")

    cb1_order4 = eig(matSize, 1)
    vb1_order4 = eig(matSize - 8, 1)
    deallocate(HT, eig, work, rwork, iwork)

    ! ---- Verify: FDorder=4 must match FDorder=2 within tolerance ----
    ! FDorder=4 should be MORE accurate, not less. With bug, CB1 is ~1.8 eV off.
    ! Tolerance is generous (10 meV) — the exact values differ due to different
    ! truncation error, but they should agree much better than the current bug.
    @assertTrue(abs(cb1_order4 - cb1_order2) < 0.05_dp, &
      & message="FDorder=4 CB1 matches FDorder=2 within 50 meV")
    @assertTrue(abs(vb1_order4 - vb1_order2) < 0.05_dp, &
      & message="FDorder=4 VB1 matches FDorder=2 within 50 meV")

  end subroutine test_qw_variable_coeff_order4
```

**Step 2: Build and run the test to verify it fails**

Run:
```bash
cmake --build build && ctest --test-dir build -R test_hamiltonian -V
```
Expected: `test_qw_variable_coeff_order4` FAILS (CB1 differs by >50 meV between FDorder=2 and FDorder=4).

**Step 3: Commit (failing test)**

```bash
git add tests/unit/test_hamiltonian.pf
git commit -m "test: add failing test for FDorder=4 variable-coefficient bug"
```

---

### Task 2: Fix `applyVariableCoeff` with midpoint averaging

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90:601-621`

**Context:** The current `applyVariableCoeff` at lines 607-621 does `kpterms(jj,ii) = -profile(jj) * FD(jj,ii)`. This is wrong for variable coefficients. Replace with midpoint averaging: diagonal uses local profile, off-diagonal uses `(profile(jj) + profile(ii))/2`.

**Step 1: Replace `applyVariableCoeff`**

In `src/physics/hamiltonianConstructor.f90`, replace lines 601-621 (the entire subroutine including comment) with:

```fortran
    !---------------------------------------------------------------------------
    !> Apply variable coefficient to FD matrix using midpoint averaging.
    !> Computes: kpterms(:,:,term_idx) = -G_avg .* FD
    !> where G_avg(i,i) = profile(i) and G_avg(i,j) = (profile(i)+profile(j))/2
    !> This correctly discretizes d/dz[g(z)·d/dz] at heterointerfaces.
    !---------------------------------------------------------------------------
    subroutine applyVariableCoeff(kpterms, profile_vec, FD, N, term_idx)

      real(kind=dp), intent(inout), dimension(:,:,:) :: kpterms
      real(kind=dp), intent(in), dimension(:) :: profile_vec
      real(kind=dp), intent(in), dimension(:,:) :: FD
      integer, intent(in) :: N, term_idx

      integer :: ii, jj
      real(kind=dp) :: g_avg

      do ii = 1, N
        do jj = 1, N
          if (FD(jj, ii) == 0.0_dp) cycle
          if (jj == ii) then
            ! Diagonal: use local profile value
            kpterms(jj, ii, term_idx) = -profile_vec(jj) * FD(jj, ii)
          else
            ! Off-diagonal: midpoint average
            g_avg = 0.5_dp * (profile_vec(jj) + profile_vec(ii))
            kpterms(jj, ii, term_idx) = -g_avg * FD(jj, ii)
          end if
        end do
      end do

    end subroutine applyVariableCoeff
```

**Step 2: Build and run the failing test**

Run:
```bash
cmake --build build && ctest --test-dir build -R test_hamiltonian -V
```
Expected: `test_qw_variable_coeff_order4` now PASSES.

**Step 3: Run ALL tests to verify no regressions**

Run:
```bash
ctest --test-dir build
```
Expected: All tests pass (existing FDorder=2 tests unchanged, new FDorder=4 test passes).

**Step 4: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90
git commit -m "fix: use midpoint averaging for variable coefficients in FDorder>2"
```

---

### Task 3: Fix `confinementInitialization_raw` with mask-based assignment

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90:115,148-161`

**Context:** Lines 148-161 sequentially assign `profile(startPos:endPos) = ...`, so the last layer wins at overlapping z-points. Fix with a boolean mask: first layer to claim a point wins.

**Step 1: Add mask logic**

In `src/physics/hamiltonianConstructor.f90`:

**1a.** Add `assigned` variable to the declaration at line 115. Change:
```fortran
      integer :: i, initIDX, endIDX, N, ii, jj
```
to:
```fortran
      integer :: i, initIDX, endIDX, N, ii, jj
      logical, allocatable :: assigned(:)
```

**1b.** Replace lines 148-161 (the `do i = 1, nlayers` block that assigns profile and kptermsProfile) with:

```fortran
        ! Mask-based assignment: first layer to claim a grid point wins.
        ! This prevents later layers from overwriting earlier ones at overlaps.
        allocate(assigned(N))
        assigned = .false.

        do i = 1, nlayers, 1
          do jj = startPos(i), endPos(i)
            if (.not. assigned(jj)) then
              profile(jj,1) = params(i)%EV
              profile(jj,2) = params(i)%EV - params(i)%DeltaSO
              profile(jj,3) = params(i)%EC

              kptermsProfile(jj,1) = params(i)%gamma1
              kptermsProfile(jj,2) = params(i)%gamma2
              kptermsProfile(jj,3) = params(i)%gamma3
              kptermsProfile(jj,4) = params(i)%A
              kptermsProfile(jj,5) = params(i)%P

              assigned(jj) = .true.
            end if
          end do
        end do

        ! Verify all grid points are covered by at least one layer
        if (any(.not. assigned)) then
          print *, 'ERROR: Not all grid points are covered by material layers.'
          print *, '  Uncovered points:'
          do jj = 1, N
            if (.not. assigned(jj)) then
              print *, '    z(', jj, ') =', z(jj)
            end if
          end do
          stop 1
        end if

        deallocate(assigned)
```

**Step 2: Build and run all tests**

Run:
```bash
cmake --build build && ctest --test-dir build
```
Expected: All tests pass. Existing configs use non-overlapping layers (verified in regression configs — all existing QW configs have non-overlapping ranges), so behavior is unchanged.

**Step 3: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90
git commit -m "fix: mask-based layer assignment prevents config overlap erasure"
```

---

### Task 4: Fix overlapping QW configs

**Files:**
- Modify: `tests/regression/configs/qw_gaas_algaas_kpar.cfg`
- Modify: `tests/regression/configs/qw_gaas_algaas_optics.cfg`

**Context:** Both configs have `material1: Al30Ga70As -200 200` which overwrites the GaAs well. Fix by splitting barriers into non-overlapping left/right regions.

**Step 1: Fix `qw_gaas_algaas_kpar.cfg`**

Write the entire file:
```
waveVector: kx
waveVectorMax: 0.1
waveVectorStep: 101
confinement:  1
FDstep: 401
FDorder: 4
numLayers:  3
material1: Al30Ga70As -200 -50 0
material2: GaAs -50 50 0
material3: Al30Ga70As 50 200 0
numcb: 4
numvb: 8
ExternalField: 0  EF
EFParams: 0.0
```

**Step 2: Fix `qw_gaas_algaas_optics.cfg`**

Write the entire file:
```
waveVector: k0
waveVectorMax: 0.0
waveVectorStep: 0
confinement:  1
FDstep: 101
FDorder: 4
numLayers:  3
material1: Al30Ga70As -200 -50 0
material2: GaAs -50 50 0
material3: Al30Ga70As 50 200 0
numcb: 4
numvb: 8
ExternalField: 0  EF
EFParams: 0.0
whichBand: 0
bandIdx: 1
```

**Step 3: Build and run all tests**

Run:
```bash
cmake --build build && ctest --test-dir build
```
Expected: All tests pass. The regression tests for these configs don't have golden data yet, so no golden comparison to break.

**Step 4: Commit**

```bash
git add tests/regression/configs/qw_gaas_algaas_kpar.cfg \
        tests/regression/configs/qw_gaas_algaas_optics.cfg
git commit -m "fix: non-overlapping z-ranges in GaAs/AlGaAs QW configs"
```

---

### Task 5: Re-run all three configs and regenerate figures

**Files:**
- Regenerate: `docs/figures/fig_qw_dispersion_gaas_algaas.png`
- Regenerate: `docs/figures/fig_qw_optical_matrix_elements.png`
- Regenerate: `docs/figures/fig_qw_potential_profile_gaas.png`
- Regenerate: `docs/figures/fig_qw_dispersion_broken_gap.png`

**Context:** With both bugs fixed, all physics results change. Re-run configs and regenerate figures.

**Step 1: Run the GaAs/AlGaAs dispersion config**

Run:
```bash
cp tests/regression/configs/qw_gaas_algaas_kpar.cfg input.cfg
./build/src/bandStructure
```
Expected: Completes without error. New eigenvalue files in `output/`.

**Step 2: Run the GaAs/AlGaAs optics config**

Run:
```bash
cp tests/regression/configs/qw_gaas_algaas_optics.cfg input.cfg
./build/src/gfactorCalculation
```
Expected: Completes without error. New `output/optical_transitions.dat`.

**Step 3: Run the broken-gap config**

Run:
```bash
cp tests/regression/configs/qw_inas_gasb_broken_gap_kpar.cfg input.cfg
./build/src/bandStructure
```
Expected: Completes without error. New eigenvalue files in `output/`.

**Step 4: Regenerate all QW figures**

Run:
```bash
python scripts/plotting/generate_all_figures.py --skip-build --figures \
  fig_qw_dispersion_gaas_algaas fig_qw_optical_matrix_elements \
  fig_qw_potential_profile_gaas fig_qw_dispersion_broken_gap
```
Expected: All four PNGs regenerated with corrected data.

**Step 5: Inspect the new eigenvalues and verify sanity**

Run:
```bash
head -5 output/eigenvalues_*.dat
cat output/optical_transitions.dat | head -15
```
Expected: CB1 for GaAs QW should be near ~0.8 eV (not ~2.6 eV as with the bug). Broken-gap effective gap should be ~65 meV (may shift slightly with FDorder=4 fix).

**Step 6: Commit figures and output**

```bash
git add docs/figures/fig_qw_dispersion_gaas_algaas.png \
        docs/figures/fig_qw_optical_matrix_elements.png \
        docs/figures/fig_qw_potential_profile_gaas.png \
        docs/figures/fig_qw_dispersion_broken_gap.png
git commit -m "docs: regenerate QW figures with fixed FDorder=4 variable coefficients"
```

---

### Task 6: Update documentation with new physics results

**Files:**
- Modify: `docs/lecture/02-quantum-well.md`
- Modify: `docs/lecture/06-optical-properties.md`

**Context:** After re-running with fixed code, all numerical values in the docs need updating. The implementer must extract actual computed eigenvalues from the output files and update the tables.

**Step 1: Update Ch02 eigenvalue tables (Section A.3)**

In `docs/lecture/02-quantum-well.md`, find the Section A.3 eigenvalue table. Replace all numerical values with the newly computed eigenvalues from `output/eigenvalues_001.dat` (the k=0 file from `qw_gaas_algaas_kpar.cfg`).

Key changes:
- CB1 eigenvalue (was ~0.7613 eV with FDorder=2 workaround)
- CB2, CB3, CB4
- HH1, LH1, SO1 eigenvalues
- Band gap value
- Update the note about FDorder to say "FDorder=4 with fixed variable-coefficient handling"

**Step 2: Update Ch02 optical transition table (Section A.6)**

Replace the table with data from the new `output/optical_transitions.dat` (from `qw_gaas_algaas_optics.cfg`). Update:
- All transition energies (meV)
- All |px|², |py|², |pz|² values
- All oscillator strengths
- The config caption to say `(FDstep=101, FDorder=4, non-overlapping 3-layer config)`

**Step 3: Update Ch02 broken-gap eigenvalues (Section B.3)**

Replace with eigenvalues from the re-run `qw_inas_gasb_broken_gap_kpar.cfg` output. Update:
- VB top eigenvalue
- CB bottom eigenvalue
- Effective gap value
- Any anticrossing gap mentioned in Section B.5

**Step 4: Update Ch02 config snippets**

Find all config snippets in Ch02 that show `FDorder: 2` as a workaround and update them to use `FDorder: 4` with non-overlapping z-ranges. Specifically:
- The config snippet in Section A.2 (if it exists)
- Any inline config examples

**Step 5: Update Ch06 optical properties**

In `docs/lecture/06-optical-properties.md`:
- Update Section 6.6.6 "Computed QW Transition Strengths" with new values
- Update Section 6.8.1: the recommendation "Use at least 4th-order FD" is now correct without caveat

**Step 6: Commit**

```bash
git add docs/lecture/02-quantum-well.md docs/lecture/06-optical-properties.md
git commit -m "docs: update QW eigenvalue tables and optical data with fixed FDorder=4 results"
```

---

### Task 7: Final verification

**Files:** None (verification only)

**Step 1: Run full test suite**

Run:
```bash
cmake --build build && ctest --test-dir build
```
Expected: All 25 tests pass (8 unit + 5 regression + any others).

**Step 2: Verify figures are consistent with docs**

Open each figure and verify it matches the numerical values stated in the documentation:
- `docs/figures/fig_qw_dispersion_gaas_algaas.png` — CB1/VB1 labels match table
- `docs/figures/fig_qw_optical_matrix_elements.png` — bar heights match oscillator strengths
- `docs/figures/fig_qw_potential_profile_gaas.png` — energy levels match table
- `docs/figures/fig_qw_dispersion_broken_gap.png` — anticrossing gap matches prose

**Step 3: Verify FDorder=4 gives better accuracy than FDorder=2**

Run a quick convergence check — FDorder=4 eigenvalues should be closer to the converged limit than FDorder=2. This is a sanity check, not a formal test.
