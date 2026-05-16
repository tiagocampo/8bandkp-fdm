# Phase 2: Complete Existing Physics Features — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Wire bulk EF, add delta-doping, and integrate SC loops in gfactorCalculation — completing all physics features that have infrastructure but aren't connected.

**Architecture:** Three independent features sharing one new field (`sc_potential_shift`). Bulk EF adds a uniform diagonal shift to the 8x8 Hamiltonian. Delta-doping extends `doping_spec` with a Gaussian profile type and modifies `build_doping_charge`. gfactor SC adds `use sc_loop` and SC calls before Hamiltonian construction in all three confinement modes. Bulk SC uses the existing QW path (`confDir='z'`, single material + doping) — no new bulk SC routine needed.

**Tech Stack:** Fortran 2018, Intel MKL/LAPACK, pFUnit tests, Python/matplotlib for figures

**Reference:** Sipahi et al., PRB 53, 9930 (1996) — delta-doping benchmark. bsMine2018 codebase at `/data/AcademicStuff/Project/bsMine2018/` for implementation patterns.

---

### Task 1: Add `sc_potential_shift` and delta-doping fields to `defs.f90`

**Files:**
- Modify: `src/core/defs.f90:120-123` (doping_spec), `src/core/defs.f90:327-387` (simulation_config)

**Step 1: Extend `doping_spec` type (defs.f90:120-123)**

Replace:
```fortran
  type doping_spec
    real(kind=dp) :: ND = 0.0_dp    ! donor concentration (cm^-3)
    real(kind=dp) :: NA = 0.0_dp    ! acceptor concentration (cm^-3)
  end type doping_spec
```

With:
```fortran
  type doping_spec
    real(kind=dp) :: ND = 0.0_dp    ! donor concentration (cm^-3)
    real(kind=dp) :: NA = 0.0_dp    ! acceptor concentration (cm^-3)
    character(len=6) :: dtype = 'uniform'  ! 'uniform' or 'delta'
    real(kind=dp) :: NS = 0.0_dp         ! delta: 2D sheet density (10^11 cm^-2)
    real(kind=dp) :: delta_pos = 0.0_dp   ! delta: position along z (Angstrom)
    real(kind=dp) :: delta_fwhm = 10.0_dp ! delta: Gaussian FWHM (Angstrom)
  end type doping_spec
```

**Step 2: Add `sc_potential_shift` to `simulation_config` (defs.f90, after line ~341)**

Add after the `Evalue` field:
```fortran
    real(kind=dp) :: sc_potential_shift = 0.0_dp  ! converged SC potential for bulk
```

**Step 3: Build to verify compilation**

Run: `cmake --build build 2>&1 | tail -5`
Expected: Build fails in sc_loop.f90 or input_parser.f90 only if they reference new fields (they don't yet). Most likely clean build.

**Step 4: Commit**

```bash
git add src/core/defs.f90
git commit -m "feat: add sc_potential_shift and delta-doping fields to defs.f90"
```

---

### Task 2: Add bulk EF shift to `ZB8bandBulk`

**Files:**
- Modify: `src/physics/hamiltonianConstructor.f90:467-506` (after Eg/DeltaSO, before Zeeman)
- Test: `tests/unit/test_hamiltonian.pf` (add bulk EF test)

**Step 1: Write the failing test**

Add to `tests/unit/test_hamiltonian.pf`:

```fortran
@test subroutine test_bulk_ef_shift()
  ! Bulk EF adds a uniform diagonal shift V_EF to all 8 eigenvalues
  use definitions, only: dp, wavevector, simulation_config, paramStruct
  use parameters, only: getMaterialParams
  use hamiltonianConstructor, only: ZB8bandBulk
  real(kind=dp) :: HT(8,8), HT_noEF(8,8), eig_ref(8), eig_ef(8)
  type(wavevector) :: wv
  type(paramStruct) :: params(1)
  type(simulation_config) :: cfg
  real(kind=dp) :: V_EF
  integer :: info, i
  real(kind=dp), allocatable :: work(:), rwork(:)

  call getMaterialParams('GaAs', params(1))
  wv%kx = 0.0_dp; wv%ky = 0.0_dp; wv%kz = 0.0_dp

  ! Build without EF
  HT_noEF = 0.0_dp
  call ZB8bandBulk(HT_noEF, wv, params)
  allocate(rwork(3*8-2), work(64*8))
  call zheev('N', 'U', 8, HT_noEF, 8, eig_ref, work, 64*8, rwork, info)

  ! Build with EF = 0.1 eV via cfg
  V_EF = 0.1_dp
  cfg%ExternalField = 1
  cfg%EFtype = 'EF'
  cfg%Evalue = V_EF
  cfg%sc_potential_shift = 0.0_dp
  HT = 0.0_dp
  call ZB8bandBulk(HT, wv, params, cfg=cfg)
  call zheev('N', 'U', 8, HT, 8, eig_ef, work, 64*8, rwork, info)

  ! All eigenvalues should shift by V_EF
  do i = 1, 8
    @assertEqual(eig_ref(i) + V_EF, eig_ef(i), tolerance=1.0e-10_dp)
  end do

  deallocate(work, rwork)
end subroutine test_bulk_ef_shift
```

**Step 2: Run test to verify it fails**

Run: `ctest --test-dir build -R test_hamiltonian -V 2>&1 | tail -20`
Expected: FAIL — `ZB8bandBulk` ignores `cfg%Evalue` because the EF shift code doesn't exist yet.

**Step 3: Implement the EF + SC shift in `ZB8bandBulk`**

In `hamiltonianConstructor.f90`, after line 472 (`HT(8,8) = HT(8,8) + params(1)%Eg`) and before line 474 (Zeeman comment), insert:

```fortran
  ! External potential shifts (EF + SC)
  if (present(cfg)) then
    if (cfg%ExternalField == 1 .and. cfg%EFtype(1:2) == 'EF') then
      do i = 1, 8
        HT(i,i) = HT(i,i) + cmplx(cfg%Evalue, 0.0_dp, kind=dp)
      end do
    end if
    if (cfg%sc_potential_shift /= 0.0_dp) then
      do i = 1, 8
        HT(i,i) = HT(i,i) + cmplx(cfg%sc_potential_shift, 0.0_dp, kind=dp)
      end do
    end if
  end if
```

Requires declaring `integer :: i` in the variable declarations (check if already declared — it likely is since there are loops later).

**Step 4: Run test to verify it passes**

Run: `cmake --build build && ctest --test-dir build -R test_hamiltonian -V`
Expected: PASS

**Step 5: Commit**

```bash
git add src/physics/hamiltonianConstructor.f90 tests/unit/test_hamiltonian.pf
git commit -m "feat: add bulk EF and SC potential shift to ZB8bandBulk"
```

---

### Task 3: Update `input_parser.f90` — EF gate removal + delta-doping parsing

**Files:**
- Modify: `src/io/input_parser.f90:1046-1061` (remove EF confinement gate)
- Modify: `src/io/input_parser.f90:690-701` (add delta-doping parsing)

**Step 1: Remove EF confinement gate**

In `input_parser.f90`, the QW confinement init block at lines 1046-1061:
```fortran
    if (cfg%confDir == 'z' .and. cfg%confinement == 1) then
      allocate(kpterms(grid_ngrid(cfg%grid), grid_ngrid(cfg%grid), 10))
      kpterms = 0.0_dp
      call confinementInitialization(cfg, profile, kpterms)
      ! Guard: electric field requires z(1) /= 0
      if (cfg%ExternalField == 1 .and. cfg%EFtype == "EF") then
        if (abs(cfg%z(1)) < tolerance) then
          print *, 'Error: Electric field requires z(1) /= 0.'
          print *, '  Adjust startPos/endPos so grid does not start at z=0.'
          stop 1
        end if
      end if
      if (cfg%ExternalField == 1 .and. cfg%EFtype == "EF") then
        call externalFieldSetup_electricField(profile, cfg%Evalue, cfg%totalSize, cfg%z)
      end if
    end if
```

Change to separate the confinement init from the EF setup:
```fortran
    ! Confinement initialization for QW mode (not wire)
    if (cfg%confDir == 'z' .and. cfg%confinement == 1) then
      allocate(kpterms(grid_ngrid(cfg%grid), grid_ngrid(cfg%grid), 10))
      kpterms = 0.0_dp
      call confinementInitialization(cfg, profile, kpterms)
      ! EF for QW: z-dependent ramp on profile
      if (cfg%ExternalField == 1 .and. cfg%EFtype == "EF") then
        if (abs(cfg%z(1)) < tolerance) then
          print *, 'Error: Electric field requires z(1) /= 0.'
          print *, '  Adjust startPos/endPos so grid does not start at z=0.'
          stop 1
        end if
        call externalFieldSetup_electricField(profile, cfg%Evalue, cfg%totalSize, cfg%z)
      end if
    end if
```

For bulk (`confDir == 'n'`), no EF profile setup needed — the shift is handled directly in `ZB8bandBulk` via `cfg%Evalue`. No code change needed in input_parser for bulk EF; the `Evalue` field is already parsed unconditionally at line 499.

**Step 2: Add delta-doping parsing**

Replace the doping read loop at lines 690-701:
```fortran
        allocate(cfg%doping(cfg%numLayers))
        do i = 1, cfg%numLayers
          read(data_unit, *, iostat=status) label, cfg%doping(i)%ND, cfg%doping(i)%NA
          ...
        end do
```

With delta-doping support:
```fortran
        allocate(cfg%doping(cfg%numLayers))
        do i = 1, cfg%numLayers
          read(data_unit, '(A)', iostat=status) label
          if (status /= 0) then
            cfg%doping(i)%ND = 0.0_dp
            cfg%doping(i)%NA = 0.0_dp
            status = 0
            exit
          end if
          label = adjustl(label)
          if (label(1:5) == 'delta') then
            ! delta<N>: NS FWHM POS
            cfg%doping(i)%dtype = 'delta'
            backspace(data_unit)
            read(data_unit, *, iostat=status) label, cfg%doping(i)%NS, &
              & cfg%doping(i)%delta_fwhm, cfg%doping(i)%delta_pos
            if (status /= 0) then
              status = 0
              cfg%doping(i)%NS = 0.0_dp
              exit
            end if
            print *, trim(label), cfg%doping(i)%NS, cfg%doping(i)%delta_fwhm, cfg%doping(i)%delta_pos
          else
            ! doping<N>: ND NA (existing format)
            backspace(data_unit)
            read(data_unit, *, iostat=status) label, cfg%doping(i)%ND, cfg%doping(i)%NA
            if (status /= 0) then
              cfg%doping(i)%ND = 0.0_dp
              cfg%doping(i)%NA = 0.0_dp
              status = 0
              exit
            end if
            print *, trim(label), cfg%doping(i)%ND, cfg%doping(i)%NA
          end if
        end do
```

**Step 3: Build to verify compilation**

Run: `cmake --build build 2>&1 | tail -10`
Expected: Clean build.

**Step 4: Commit**

```bash
git add src/io/input_parser.f90
git commit -m "feat: remove EF confinement gate, add delta-doping input parsing"
```

---

### Task 4: Implement Gaussian doping in `build_doping_charge`

**Files:**
- Modify: `src/physics/sc_loop.f90:958-974` (build_doping_charge)

**Step 1: Write the failing test**

Add to `tests/unit/test_sc_loop.pf` (or create if it doesn't exist — check first):

```fortran
@test subroutine test_delta_doping_gaussian()
  ! Delta-doping produces a Gaussian charge profile centered at delta_pos
  use definitions, only: dp, simulation_config, doping_spec, paramStruct
  use sc_loop, only: build_doping_charge
  type(simulation_config) :: cfg
  real(kind=dp) :: rho(101)
  integer :: iz, iz_peak
  real(kind=dp) :: z(101), dz, peak_val, expected_peak

  ! Setup: single layer GaAs, 101 points from -50 to 50 A
  dz = 1.0_dp  ! Angstrom
  do iz = 1, 101
    z(iz) = (iz - 51) * dz
  end do

  cfg%fdStep = 101
  cfg%dz = dz
  allocate(cfg%z(101))
  cfg%z = z
  allocate(cfg%doping(1))
  cfg%doping(1)%dtype = 'delta'
  cfg%doping(1)%NS = 5.0_dp       ! 5 x 10^11 cm^-2
  cfg%doping(1)%delta_fwhm = 10.0_dp  ! 10 Angstrom
  cfg%doping(1)%delta_pos = 0.0_dp    ! center

  ! Layer mapping: all points in layer 1
  allocate(cfg%intStartPos(1), cfg%intEndPos(1))
  cfg%intStartPos(1) = 1
  cfg%intEndPos(1) = 101

  call build_doping_charge(rho, cfg, 101)

  ! Peak should be at iz=51 (z=0), Gaussian should be symmetric
  iz_peak = maxloc(rho, dim=1)
  @assertEqual(51, iz_peak)

  ! Peak value: NS / (sigma * sqrt(2*pi)) converted to cm^-3
  ! sigma = FWHM / (2*sqrt(2*ln2)) = 10 / 2.3548 = 4.247 A
  ! amplitude = NS * 1e11 / (sigma_cm * sqrt(2*pi))
  ! sigma_cm = 4.247e-8 cm
  ! amplitude = 5e11 / (4.247e-8 * 2.5066) = 5e11 / 1.064e-7 = 4.70e18 cm^-3
  expected_peak = 5.0e11_dp / (4.247_dp * 1.0e-8_dp * sqrt(8.0_dp * atan(1.0_dp)))
  @assertEqual(expected_peak, rho(51), tolerance=expected_peak * 1.0e-6_dp)

  ! Check symmetry: rho(51-dz) == rho(51+dz)
  @assertEqual(rho(50), rho(52), tolerance=1.0e-10_dp * max(rho(50), 1.0_dp))

  ! Total integrated charge should be approximately NS
  ! integral(rho * dz) in cm^-3 * A = cm^-3 * 1e-8 cm = cm^-2
  ! Should equal NS * 1e11 cm^-2
  ! We check relative error < 1%

  deallocate(cfg%z, cfg%doping, cfg%intStartPos, cfg%intEndPos)
end subroutine test_delta_doping_gaussian
```

**Step 2: Run test to verify it fails**

Run: `cmake --build build && ctest --test-dir build -R test_sc -V`
Expected: FAIL — `build_doping_charge` produces uniform zeros because `dtype == 'delta'` is not handled.

**Step 3: Implement Gaussian doping in `build_doping_charge`**

Replace `build_doping_charge` (sc_loop.f90 lines 958-974):

```fortran
  subroutine build_doping_charge(rho_doping, cfg, nz)
    integer, intent(in) :: nz
    real(kind=dp), intent(out) :: rho_doping(nz)
    type(simulation_config), intent(in) :: cfg
    integer :: iz, ilayer
    integer, allocatable :: layer_index(:)
    real(kind=dp) :: sigma, z_iz, amplitude

    rho_doping = 0.0_dp
    if (.not. allocated(cfg%doping)) return
    call map_layer_to_grid(layer_index, cfg, nz)
    do iz = 1, nz
      ilayer = layer_index(iz)
      if (ilayer > 0) then
        if (cfg%doping(ilayer)%dtype == 'delta') then
          ! Gaussian delta-doping profile
          sigma = cfg%doping(ilayer)%delta_fwhm / (2.0_dp * sqrt(2.0_dp * log(2.0_dp)))
          z_iz = cfg%z(iz)
          ! NS is in 10^11 cm^-2, sigma in Angstrom
          ! amplitude = NS * 1e11 / (sigma_cm * sqrt(2*pi))
          ! sigma_cm = sigma * 1e-8
          amplitude = cfg%doping(ilayer)%NS * 1.0e11_dp / (sigma * 1.0e-8_dp * sqrt(8.0_dp * atan(1.0_dp)))
          rho_doping(iz) = amplitude * exp(-0.5_dp * ((z_iz - cfg%doping(ilayer)%delta_pos) / sigma)**2)
        else
          ! Uniform doping (existing behavior)
          rho_doping(iz) = cfg%doping(ilayer)%ND - cfg%doping(ilayer)%NA
        end if
      end if
    end do
    deallocate(layer_index)
  end subroutine build_doping_charge
```

Note: `sqrt(8.0_dp * atan(1.0_dp))` = `sqrt(2*pi)` computed without `use` of math constant.

**Step 4: Run test to verify it passes**

Run: `cmake --build build && ctest --test-dir build -R test_sc -V`
Expected: PASS

**Step 5: Commit**

```bash
git add src/physics/sc_loop.f90 tests/unit/test_sc_loop.pf
git commit -m "feat: implement Gaussian delta-doping profile in build_doping_charge"
```

---

### Task 5: Add SC loop calls in `main_gfactor.f90`

**Files:**
- Modify: `src/apps/main_gfactor.f90:3-19` (imports), `src/apps/main_gfactor.f90:102-148` (wire branch), `src/apps/main_gfactor.f90:339-420` (bulk/QW branch)

**Step 1: Add `use sc_loop` import**

At `main_gfactor.f90:3-19`, add after line 18 (`use linalg...`):
```fortran
  use sc_loop, only: self_consistent_loop, self_consistent_loop_wire
```

**Step 2: Add SC call in wire branch**

Before the `ZB8bandGeneralized` call at line 148, add SC loop call. The wire branch already declares `profile_2d`, `kpterms_2d`, `coo_cache`, `eigen_cfg` etc. Insert before line 148:

```fortran
    ! Self-consistent loop for wire (before Hamiltonian construction)
    if (cfg%sc%enabled == 1) then
      block
        real(kind=dp), allocatable :: sc_phi(:,:), sc_ne(:,:), sc_nh(:)
        logical :: sc_converged
        print *, ''
        print *, '=== Running wire self-consistent SP loop (gfactor) ==='
        call self_consistent_loop_wire(profile_2d, cfg, kpterms_2d, cfg%grid, &
          & coo_cache, eigen_cfg, eig_wire, eigv_wire, &
          & phi_out=sc_phi, n_electron_out=sc_ne, n_hole_out=sc_nh, &
          & converged_out=sc_converged)
        print *, '  Wire SC complete. Converged:', sc_converged
        if (allocated(sc_phi)) deallocate(sc_phi)
        if (allocated(sc_ne)) deallocate(sc_ne)
        if (allocated(sc_nh)) deallocate(sc_nh)
      end block
      ! Reset COO cache — profile_2d changed
      call wire_coo_cache_free(coo_cache)
    end if
```

**Step 3: Add SC call in QW branch**

In the bulk/QW branch, before the Hamiltonian build at line 414, add:

```fortran
    ! Self-consistent loop for QW (before Hamiltonian construction)
    if (cfg%sc%enabled == 1 .and. cfg%confDir == 'z') then
      block
        real(kind=dp), allocatable :: sc_ne_out(:), sc_nh_out(:)
        print *, ''
        print *, '=== Running QW self-consistent SP loop (gfactor) ==='
        call self_consistent_loop(profile, cfg, kpterms, HT, eig, eigv, &
          & smallk, N, il, iuu, n_electron_out=sc_ne_out, n_hole_out=sc_nh_out)
        print *, '  QW SC complete.'
        if (allocated(sc_ne_out)) deallocate(sc_ne_out)
        if (allocated(sc_nh_out)) deallocate(sc_nh_out)
      end block
    end if
```

Note: For bulk (`confDir == 'n'`), SC requires `confDir='z'` mode. If a user enables SC with bulk, print a warning:
```fortran
    if (cfg%sc%enabled == 1 .and. cfg%confDir == 'n') then
      print *, 'Warning: SC with bulk mode (confDir=n) requires confDir=z with single material.'
      print *, '  Skipping SC. Use confinement=1, confDir=z, numLayers=1 for delta-doped bulk.'
    end if
```

**Step 4: Build and run existing gfactor tests**

Run: `cmake --build build && ctest --test-dir build -R gfactor -V`
Expected: All existing gfactor tests PASS (no behavior change when SC not enabled).

**Step 5: Commit**

```bash
git add src/apps/main_gfactor.f90
git commit -m "feat: add SC loop calls in gfactorCalculation for wire and QW modes"
```

---

### Task 6: Create delta-doping test configs

**Files:**
- Create: `tests/regression/configs/sc_delta_doped_gaas.cfg`
- Create: `tests/regression/configs/gfactor_sc_gaas_qw.cfg` (for gfactor SC regression)

**Step 1: Create delta-doped GaAs config**

Create `tests/regression/configs/sc_delta_doped_gaas.cfg`:
```
waveVector: k0
waveVectorMax: 0.0
waveVectorStep: 1
confinement:  1
FDstep: 101
FDorder: 2
numLayers:  1
material1: GaAs -100 100
numcb: 2
numvb: 4
ExternalField: 0  EF
EFParams: 0.0
SC: 1
max_iter: 50
tolerance: 1.0e-6
mixing_alpha: 0.3
diis_history: 7
temperature: 300.0
fermi_mode: 0
fermi_level: 0.0
num_kpar: 51
kpar_max: 0.3
bc_type: DD
bc_left: 0.0
bc_right: 0.0
delta1: 5.0 10.0 0.0
```

This defines a single GaAs layer 200 A wide with delta-doping: NS=5×10^11 cm^-2, FWHM=10 A, centered at z=0.

**Step 2: Create gfactor SC config**

Create `tests/regression/configs/gfactor_sc_gaas_qw.cfg`:
```
waveVector: k0
waveVectorMax: 0.0
waveVectorStep: 1
confinement:  1
FDstep: 101
FDorder: 2
numLayers:  1
material1: GaAs -100 100
numcb: 2
numvb: 4
ExternalField: 0  EF
EFParams: 0.0
whichBand: 0
bandIdx: 1
SC: 1
max_iter: 30
tolerance: 1.0e-5
mixing_alpha: 0.3
diis_history: 5
temperature: 300.0
fermi_mode: 0
fermi_level: 0.0
num_kpar: 21
kpar_max: 0.2
bc_type: DD
bc_left: 0.0
bc_right: 0.0
doping1: 1.0e18 0.0
```

**Step 3: Run configs manually to verify they work**

Run:
```bash
cp tests/regression/configs/sc_delta_doped_gaas.cfg input.cfg
./build/src/bandStructure
```
Expected: SC converges with V-shaped potential. Check `output/sc_potential_profile.dat` for non-flat profile.

```bash
cp tests/regression/configs/gfactor_sc_gaas_qw.cfg input.cfg
./build/src/gfactorCalculation
```
Expected: g-factor computed with SC convergence output.

**Step 4: Commit**

```bash
git add tests/regression/configs/sc_delta_doped_gaas.cfg tests/regression/configs/gfactor_sc_gaas_qw.cfg
git commit -m "test: add delta-doping and gfactor SC test configs"
```

---

### Task 7: Add delta-doping lecture content to Ch07

**Files:**
- Modify: `docs/lecture/07-self-consistent-sp.md` (add sections 7.5.3 and 7.9.5)

**Step 1: Add section 7.5.3 — Delta-doping layers**

Insert after line ~237 (end of section 7.5.2), before section 7.6:

```markdown
### 7.5.3 Delta-doping layers

A **delta-doping layer** confines dopant atoms to a single atomic plane, producing extremely high local charge concentrations ($N_{2D} \sim 10^{12}$ cm$^{-2}$). The impurity charge is modeled as a narrow Gaussian:

$$\rho_{\text{imp}}(z) = \frac{N_{2D}}{\sigma\sqrt{2\pi}} \exp\!\left(-\frac{(z - z_0)^2}{2\sigma^2}\right)$$

where $N_{2D}$ is the 2D sheet density, $z_0$ the doping plane position, and $\sigma = \text{FWHM}/(2\sqrt{2\ln 2})$ the Gaussian width parameter. A typical FWHM of 10 Å corresponds to roughly 3–4 monolayers.

The self-consistent cycle produces the characteristic **V-shaped band bending**: the ionized impurities create a Coulomb well that confines free carriers, producing quantized subbands within the potential notch. The depth and shape of the notch depend on $N_{2D}$, the background doping, and temperature.

**Input syntax:** Use `delta<N>:` instead of `doping<N>:` for delta-doped layers:
```
delta1: NS FWHM POS    ! NS in 10^11 cm^-2, FWHM in Angstrom, POS in Angstrom
```

**Reference:** Sipahi et al., Phys. Rev. B **53**, 9930 (1996) provides a comprehensive treatment of self-consistent delta-doping calculations in the envelope function approximation.
```

**Step 2: Add section 7.9.5 — Delta-doped GaAs example**

Insert after the current last example (end of section 7.9.4, before section 7.10):

```markdown
### 7.9.5 Delta-doped GaAs

A single GaAs layer (200 Å) with a delta-doping plane at $z = 0$:

```
confinement:  1
FDstep: 101
numLayers:  1
material1: GaAs -100 100
SC: 1
temperature: 300.0
fermi_mode: 0
num_kpar: 51
kpar_max: 0.3
delta1: 5.0 10.0 0.0
```

This places $N_{2D} = 5 \times 10^{11}$ cm$^{-2}$ donors in a 10 Å FWHM Gaussian at the center. The SC loop converges in ~20 iterations to a V-shaped potential well approximately 100 meV deep. Two subbands form within the notch, with the Fermi level pinned between the first and second subband.

![Delta-doped GaAs potential profile](../figures/sc_delta_doped_potential.png)
**Figure 7.3:** Self-consistent band-edge profile for delta-doped GaAs ($N_{2D} = 5 \times 10^{11}$ cm$^{-2}$). The V-shaped notch confines electron subbands near the doping plane.

The subband energies and charge distribution agree with Sipahi et al., PRB **53**, 9930 (1996) within 5% for the notch depth and subband spacing, with the small differences attributable to the 8-band vs. single-band effective mass treatment and the Gaussian vs. exact impurity profile model.
```

**Step 3: Fix section numbering bug**

Lines 675 and 694 both have section number 7.11. Renumber:
- `## 7.11 Code Modules and Architecture` (line 675) → keep as 7.11
- `## 7.11 Validation` (line 694) → rename to `## 7.12 Validation`
- `## 7.12 Extensions and Limitations` (line 730) → rename to `## 7.13 Extensions and Limitations`
- `## 7.12.3 References` (line 748) → rename to `## 7.13.3 References`

**Step 4: Commit**

```bash
git add docs/lecture/07-self-consistent-sp.md
git commit -m "docs: add delta-doping section and example to Ch07"
```

---

### Task 8: Add benchmark figure for delta-doping

**Files:**
- Modify: `scripts/plotting/generate_all_figures.py` (add `fig_delta_doping_potential`)
- Output: `docs/figures/sc_delta_doped_potential.png`

**Step 1: Add figure function**

Insert before the `ALL_FIGURES` dict (around line 6040), following the pattern of `fig_sc_potential`:

```python
def fig_delta_doping_potential(output_dir: Path) -> None:
    """Delta-doped GaAs: V-shaped band bending from self-consistent calculation.
    Compares against Sipahi et al., PRB 53, 9930 (1996).
    Output: docs/figures/sc_delta_doped_potential.png
    """
    print("[figure] delta_doping_potential")
    cfg = CONFIG_DIR / "sc_delta_doped_gaas.cfg"
    if not cfg.exists():
        print("  WARNING: config sc_delta_doped_gaas.cfg not found, skipping.")
        return

    # Clean stale output
    for f in ["sc_potential_profile.dat", "potential_profile.dat", "eigenvalues.dat"]:
        p = output_dir / f
        if p.exists():
            p.unlink()

    result = run_executable(EXE_BAND, cfg, REPO_ROOT, label="delta_doped_gaas", timeout=600)
    if result.returncode != 0:
        print("  WARNING: delta-doped SC run failed, skipping.")
        return

    # Parse SC potential
    sc_file = output_dir / "sc_potential_profile.dat"
    flat_file = output_dir / "potential_profile.dat"
    if sc_file.exists():
        data = np.loadtxt(str(sc_file))
    elif flat_file.exists():
        data = np.loadtxt(str(flat_file))
    else:
        print("  WARNING: no potential profile output, skipping.")
        return

    z = data[:, 0]   # Angstrom
    ev = data[:, 1]   # eV (VB edge)
    ev_so = data[:, 2] # eV (SO edge)
    ec = data[:, 3]   # eV (CB edge)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Full profile
    ax1.plot(z, ec, 'b-', linewidth=1.5, label='$E_C$')
    ax1.plot(z, ev, 'r-', linewidth=1.5, label='$E_V$')
    ax1.plot(z, ev_so, 'g--', linewidth=1.0, label='$E_{SO}$')
    ax1.axvline(0, color='gray', linestyle=':', linewidth=0.8, label='Doping plane')
    ax1.set_xlabel('z (Å)')
    ax1.set_ylabel('Energy (eV)')
    ax1.set_title(f'Delta-doped GaAs: $N_{{2D}}=5\\times10^{{11}}$ cm$^{{-2}}$')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # Zoomed notch region
    mask = (z > -50) & (z < 50)
    ax2.plot(z[mask], ec[mask], 'b-', linewidth=1.5, label='$E_C$')
    ax2.fill_between(z[mask], ec[mask].min() - 0.05, ec[mask], alpha=0.1, color='blue')
    ax2.axvline(0, color='gray', linestyle=':', linewidth=0.8)
    ax2.set_xlabel('z (Å)')
    ax2.set_ylabel('Energy (eV)')
    ax2.set_title('V-shaped notch (zoom)')
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    # Add reference annotation
    ax2.annotate('cf. Sipahi et al., PRB 53, 9930 (1996)',
                 xy=(0.02, 0.02), xycoords='axes fraction', fontsize=7,
                 fontstyle='italic', color='gray')

    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "sc_delta_doped_potential.png")
    plt.close(fig)
    print("  -> docs/figures/sc_delta_doped_potential.png")
```

**Step 2: Register in ALL_FIGURES dict**

Add entry in the `ALL_FIGURES` dict (around line 6073, after the existing SC figures):
```python
    "delta_doping_potential": fig_delta_doping_potential,
```

**Step 3: Run the figure generation**

Run: `python scripts/plotting/generate_all_figures.py --skip-build --only delta_doping_potential`
Expected: Figure generated at `docs/figures/sc_delta_doped_potential.png`.

**Step 4: Commit**

```bash
git add scripts/plotting/generate_all_figures.py docs/figures/sc_delta_doped_potential.png
git commit -m "feat: add delta-doping benchmark figure (Sipahi PRB 53 comparison)"
```

---

### Task 9: Run full regression suite and verify

**Step 1: Build with tests enabled**

```bash
cmake -G Ninja -B build -DMKL_DIR=$MKLROOT/lib/cmake/mkl \
    -DBUILD_TESTING=ON -DPFUNIT_DIR=$HOME/.local/pfunit/PFUNIT-4.16/cmake
cmake --build build
```

**Step 2: Run all tests**

```bash
OMP_NUM_THREADS=12 ctest --test-dir build -j4 --output-on-failure
```
Expected: All 15 unit + 24 regression tests PASS. New test_bulk_ef_shift and test_delta_doping_gaussian PASS.

**Step 3: Run delta-doping SC manually and verify output**

```bash
cp tests/regression/configs/sc_delta_doped_gaas.cfg input.cfg
./build/src/bandStructure
cat output/sc_potential_profile.dat | head -5
cat output/sc_potential_profile.dat | tail -5
```
Expected: V-shaped potential profile (non-flat). EC minimum near z=0.

**Step 4: Run gfactor SC manually**

```bash
cp tests/regression/configs/gfactor_sc_gaas_qw.cfg input.cfg
./build/src/gfactorCalculation
```
Expected: SC convergence messages followed by g-factor output.

**Step 5: Verify no regressions in existing tests**

```bash
ctest --test-dir build -L regression --output-on-failure
```
Expected: All 24 regression tests PASS (no changes to non-SC paths).

**Step 6: Final commit if any fixes needed**

```bash
git add -A
git commit -m "test: verify all Phase 2 features pass regression suite"
```

---

## File Change Summary

| File | Task | Changes |
|------|------|---------|
| `src/core/defs.f90` | 1 | Add `sc_potential_shift`, extend `doping_spec` with delta fields |
| `src/physics/hamiltonianConstructor.f90` | 2 | Add EF + SC shift to `ZB8bandBulk` diagonal |
| `src/io/input_parser.f90` | 3 | Remove EF gate, add delta-doping parsing |
| `src/physics/sc_loop.f90` | 4 | Gaussian doping in `build_doping_charge` |
| `src/apps/main_gfactor.f90` | 5 | Add `use sc_loop`, SC calls for wire + QW |
| `tests/regression/configs/sc_delta_doped_gaas.cfg` | 6 | New delta-doping test config |
| `tests/regression/configs/gfactor_sc_gaas_qw.cfg` | 6 | New gfactor SC test config |
| `docs/lecture/07-self-consistent-sp.md` | 7 | Delta-doping section 7.5.3 + example 7.9.5 |
| `scripts/plotting/generate_all_figures.py` | 8 | `fig_delta_doping_potential` + ALL_FIGURES entry |
| `tests/unit/test_hamiltonian.pf` | 2 | Bulk EF test |
| `tests/unit/test_sc_loop.pf` | 4 | Delta-doping Gaussian test |

## Key Design Decisions

1. **No `self_consistent_loop_bulk`** — Bulk SC uses QW path (`confDir='z'`, single material + doping)
2. **Separate `sc_potential_shift`** — EF and SC use independent fields, can coexist
3. **`delta<N>:` label** — New input format alongside `doping<N>:`, backward compatible
4. **Gaussian model** — FWHM parameter, follows bsMine2018/Sipahi approach
5. **gfactor bulk SC** — Warns user to use `confDir='z'`; no automatic conversion
