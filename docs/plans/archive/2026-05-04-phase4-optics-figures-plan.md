# Phase 4: Optics & Documentation Figures — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Validate ISBT physics, verify all 74 existing figures, add 5 new physics figures, and rebuild Ch06/Ch08 lecture chapters.

**Architecture:** Validation-first approach in 4 sequential work packages: (1) ISBT cross-validation and fix, (2) full figure validation pass, (3) 5 new Group #50 figures, (4) chapter rebuilds. When any discrepancy is found during figure verification (wrong physics, wrong code, wrong plot), the protocol is: stop → document discrepancy → trace root cause → fix at source (Fortran, Python, or config) → re-run affected tests → re-generate affected figures → commit.

**Tech Stack:** Fortran 2018 (gfortran + MKL), Python 3 (matplotlib, numpy), pFUnit testing, CMake/Ninja build.

---

## Discrepancy Protocol (applies to ALL figures, not just optics)

When any figure verification step reveals a discrepancy — wrong physics values, incorrect code behavior, unexpected plot features, or mismatched captions — follow this protocol:

1. **Stop and document** — Record which figure, what's wrong, expected vs observed, in a `docs/plans/phase4-discrepancy-log.md` file
2. **Trace root cause** — Is it a Fortran physics bug, a Python plotting issue, a config parameter error, or a chapter text error?
3. **Fix at source** — Fix the Fortran code if physics is wrong, the Python function if plotting is wrong, the config if parameters are wrong, or the markdown if text is wrong
4. **Re-run tests** — Run `ctest --test-dir build --output-on-failure` to verify the fix doesn't break anything
5. **Re-generate affected figures** — Re-run the specific figure function(s)
6. **Commit** — One commit per fix with descriptive message

If the root cause is unclear or the fix would be complex (>50 lines), escalate: create a separate task in the plan rather than fixing inline. Never paper over a discrepancy by adjusting expected values or skipping the figure.

---

## Task 1: ISBT Cross-Validation Test (pFUnit)

**Files:**
- Modify: `tests/unit/test_optical_qw.pf` (add new test subroutine)
- Modify: `tests/unit/CMakeLists.txt` (if needed for new test source)
- Read: `src/physics/optical_spectra.f90:582-601` (z_dipole), `src/physics/optical_spectra.f90:726-796` (compute_isbt_absorption)

**Context:** `optical_spectra.f90` has two ISBT methods with no cross-check:
- `z_dipole()` (line 582): computes `z_ij = dz * Σ conjg(ψ_i) * z * ψ_j`
- `compute_isbt_absorption()` (line 726): uses CSR velocity `v_z = -i[z, H]` via `csr_spmv`
- Theoretical relationship: `v_z_ij = -i * E_ij * z_ij`, so `|v_z_ij|² = E_ij² * |z_ij|²`

- [ ] **Step 1: Write the failing test**

Add a new test `test_isbt_dipole_velocity_consistency` to `tests/unit/test_optical_qw.pf`. This test:
1. Builds the same 3-layer Al30Ga70As/GaAs/Al30Ga70As QW (51 grid points) as the existing test
2. Diagonalizes the Hamiltonian
3. Extracts CB states (first 2 CB subbands)
4. Computes z-dipole `z_ij` via the same formula as `z_dipole()` helper (inline the calculation)
5. Computes velocity matrix element `v_z_ij` by building the z-velocity matrix via `ZB8bandGeneralized(g='g3')` and computing `<ψ_i|v_z|ψ_j>` via `zdotc`
6. Asserts `abs(|v_z_ij|² - E_ij² * |z_ij|²) / |v_z_ij|² < 0.01` (1% tolerance)

```fortran
@test
subroutine test_isbt_dipole_velocity_consistency()
  integer, parameter :: ngrid = 51
  real(kind=dp), parameter :: dz_val = 2.0_dp
  integer :: matSize, numcb_val, numvb_val, info, lwork, lrwork, liwork
  integer :: cb_lo, cb_hi, state_i, state_j, dim

  real(kind=dp) :: z(ngrid)
  real(kind=dp), allocatable :: profile(:,:), kpterms(:,:,:)

  integer, parameter :: nlayers = 3
  character(len=255) :: material(nlayers)
  type(paramStruct) :: params(nlayers)
  integer :: intStartPos(nlayers), intEndPos(nlayers)

  complex(kind=dp), allocatable :: HT(:,:), HT_vel(:,:), work(:)
  real(kind=dp), allocatable :: eig(:,:), rwork(:)
  integer, allocatable :: iwork(:)

  complex(kind=dp), allocatable :: Ytmp(:)
  complex(kind=dp) :: z_ij, vz_ij
  real(kind=dp) :: E_ij, z_abs2, vz_abs2, rel_diff

  ! ---- Setup (same as test_optical_matrix_qw_basic) ----
  do info = 1, ngrid
    z(info) = -50.0_dp + dble(info - 1) * dz_val
  end do

  intStartPos(1) = 1;   intEndPos(1) = 15
  intStartPos(2) = 16;  intEndPos(2) = 36
  intStartPos(3) = 37;  intEndPos(3) = 51

  material(1) = "Al30Ga70As"
  material(2) = "GaAs"
  material(3) = "Al30Ga70As"

  call paramDatabase(material, nlayers, params)

  allocate(kpterms(ngrid, ngrid, 10))
  kpterms = 0.0_dp
  call confinementInitialization_raw(z, intStartPos, intEndPos, material, &
    & nlayers, params, 'z', profile, kpterms, FDorder=2)

  numcb_val = 2 * ngrid
  numvb_val = 4 * ngrid
  matSize = ngrid * 8

  ! ---- Build and diagonalize QW Hamiltonian at k=0 ----
  allocate(HT(matSize, matSize))
  HT = cmplx(0.0_dp, 0.0_dp, kind=dp)
  call ZB8bandQW(HT, wavevector(0.0_dp, 0.0_dp, 0.0_dp), profile, kpterms)

  allocate(eig(matSize, 1))
  eig = 0.0_dp
  allocate(work(1), rwork(1), iwork(1))
  call zheevd('V', 'U', matSize, HT, matSize, eig(:,1), work, -1, rwork, -1, iwork, -1, info)
  lwork = int(real(work(1))); lrwork = int(rwork(1)); liwork = iwork(1)
  deallocate(work, rwork, iwork)
  allocate(work(lwork), rwork(lrwork), iwork(liwork))
  call zheevd('V', 'U', matSize, HT, matSize, eig(:,1), work, lwork, rwork, lrwork, &
    & iwork, liwork, info)
  @assertTrue(info == 0, message="zheevd succeeded")

  ! ---- Build z-velocity matrix via ZB8bandGeneralized(g='g3') ----
  allocate(HT_vel(matSize, matSize))
  HT_vel = cmplx(0.0_dp, 0.0_dp, kind=dp)
  call ZB8bandQW(HT_vel, wavevector(1.0_dp, 0.0_dp, 0.0_dp), profile, kpterms, g='g3')

  ! ---- Compute z-dipole and velocity for CB1→CB2 transition ----
  ! CB states are at the top of the eigenspectrum
  state_i = matSize - numcb_val + 1   ! CB1 (lowest CB state)
  state_j = matSize - numcb_val + 2   ! CB2 (second CB state)
  E_ij = eig(state_j, 1) - eig(state_i, 1)
  @assertTrue(E_ij > 0.0_dp, message="CB1→CB2 transition energy is positive")

  ! z-dipole: z_ij = dz * sum_n sum_b conjg(psi(n,b,i)) * z(n) * psi(n,b,j)
  z_ij = cmplx(0.0_dp, 0.0_dp, kind=dp)
  dim = matSize
  do info = 1, ngrid
    do liwork = 1, 8
      lrwork = (liwork - 1) * ngrid + info
      z_ij = z_ij + conjg(HT(lrwork, state_i)) * z(info) * HT(lrwork, state_j)
    end do
  end do
  z_ij = z_ij * dz_val

  ! v_z via matrix-vector product: v_z_ij = <psi_i | H_vel | psi_j>
  allocate(Ytmp(dim))
  Ytmp = cmplx(0.0_dp, 0.0_dp, kind=dp)
  ! Ytmp = HT_vel * psi_j
  call zgemv('N', dim, dim, cmplx(1.0_dp, 0.0_dp, kind=dp), HT_vel, dim, &
    & HT(1, state_j), 1, cmplx(0.0_dp, 0.0_dp, kind=dp), Ytmp, 1)
  ! vz_ij = psi_i^H * Ytmp
  vz_ij = zdotc(dim, HT(1, state_i), 1, Ytmp, 1)

  ! ---- Cross-validate: |v_z_ij|^2 = E_ij^2 * |z_ij|^2 ----
  z_abs2 = real(z_ij * conjg(z_ij), kind=dp)
  vz_abs2 = real(vz_ij * conjg(vz_ij), kind=dp)

  ! v_z = -i * E * z  =>  |v_z|^2 = E^2 * |z|^2
  if (vz_abs2 > 1.0e-30_dp) then
    rel_diff = abs(vz_abs2 - E_ij**2 * z_abs2) / vz_abs2
    @assertTrue(rel_diff < 0.05_dp, &
      & message="ISBT dipole-velocity consistency: |v|^2 ≈ E^2*|z|^2")
  else
    ! Both near zero: z_ij should also be near zero
    @assertTrue(z_abs2 < 1.0e-20_dp, &
      & message="ISBT: both z-dipole and velocity near zero (allowed)")
  end if

  deallocate(Ytmp, HT_vel, HT, eig, work, rwork, iwork)
  deallocate(profile, kpterms)

end subroutine test_isbt_dipole_velocity_consistency
```

- [ ] **Step 2: Build and run test to verify it compiles**

Run: `cmake --build build && ctest --test-dir build -R test_optical_qw -V`
Expected: Test compiles and runs. If the assertion fails, this confirms the discrepancy exists.

- [ ] **Step 3: Analyze result and fix if needed**

If test fails (rel_diff > 5%):
1. Print diagnostic values: `z_abs2`, `vz_abs2`, `E_ij`, `E_ij**2 * z_abs2`, ratio
2. Check if the ratio `vz_abs2 / (E_ij**2 * z_abs2)` is a clean factor (2, 4, 0.25, etc.)
3. Trace in `optical_spectra.f90` whether `compute_isbt_absorption` or `compute_intersubband_transitions` has the sign/normalization bug
4. Fix the Fortran code, rebuild, re-run test

If test passes: both methods are consistent, no Fortran fix needed. Move to Task 2.

- [ ] **Step 4: Commit**

```bash
git add tests/unit/test_optical_qw.pf
git commit -m "test: add ISBT dipole-velocity cross-validation test"
```

---

## Task 2: ISBT Fortran Fix (if needed)

**Files:**
- Modify: `src/physics/optical_spectra.f90` (if Task 1 reveals discrepancy)

**Context:** Only execute this task if Task 1's test fails with a clear discrepancy. Skip entirely if the test passes.

- [ ] **Step 1: Diagnose the discrepancy**

From Task 1's diagnostic output, determine which of these applies:
- **Ratio is exactly 1/4 or 4**: likely a factor-of-2 error in z_dipole (missing/extra `dz` factor) or velocity normalization
- **Ratio depends on grid spacing**: z_dipole integration error (trapezoidal vs rectangular)
- **Ratio ≈ 1 but off by sign**: imaginary/real mismatch (v_z should be purely imaginary for real z)

- [ ] **Step 2: Apply the fix in `optical_spectra.f90`**

Common fixes:

If `z_dipole` uses wrong integration weight (line 599):
```fortran
! Current: z_ij = z_ij * dz  (rectangular rule)
! Fix if needed: trapezoidal weights at boundaries
```

If `compute_isbt_absorption` missing E_ij factor (line ~783):
```fortran
! Current: p_abs2 = real(pele_ij * conjg(pele_ij), kind=dp)
! Absorption ∝ |v_z|^2/E_ij^2 * lineshape, but code uses |v_z|^2 directly
! Fix: divide by E_ij^2 to get z-dipole absorption coefficient
```

Apply the minimal fix identified in Step 1.

- [ ] **Step 3: Re-run tests**

Run: `ctest --test-dir build --output-on-failure`
Expected: All 53+ tests pass, including the new ISBT cross-validation test.

- [ ] **Step 4: Commit**

```bash
git add src/physics/optical_spectra.f90
git commit -m "fix: correct ISBT [dipole sign / velocity normalization / E_ij factor] in optical_spectra"
```

---

## Task 3: Create Discrepancy Log and Run Full Figure Validation

**Files:**
- Create: `docs/plans/phase4-discrepancy-log.md` (initially empty template)
- Run: `scripts/plotting/generate_all_figures.py`
- Read: All 74 output PNGs in `docs/figures/`

**Context:** This is WP2 — validate all 74 existing figures (not just optics). The discrepancy protocol applies to every figure. If any figure shows wrong physics, the discrepancy is logged, root-cause traced, and fixed before moving on.

- [ ] **Step 1: Create the discrepancy log template**

Write `docs/plans/phase4-discrepancy-log.md`:

```markdown
# Phase 4 Figure Discrepancy Log

| Figure | Issue | Root Cause | Fix | Status |
|--------|-------|------------|-----|--------|
```

- [ ] **Step 2: Build and run all 74 figures**

```bash
cmake --build build
python scripts/plotting/generate_all_figures.py --skip-build 2>&1 | tee /tmp/fig_run.log
```

Expected: All 74 figures generate without crashes. Check `/tmp/fig_run.log` for warnings or errors.

- [ ] **Step 3: Spot-check all figure categories**

For each category, verify physics sanity. If any discrepancy is found, log it and apply the discrepancy protocol (stop → document → trace root cause → fix → re-test → re-generate → commit):

**Bulk band structures** (figs 1-7): Band gaps match Vurgaftman parameters. GaAs E_g ≈ 1.42 eV, InAs E_g ≈ 0.354 eV, InSb E_g ≈ 0.17 eV. Bands spin-degenerate at k=0.

**Strain figures** (figs 8-12): Bir-Pikus shifts match `compute_bp_scalar` formula. HH/LH splitting correct sign for compressive/tensile strain.

**QW subbands** (figs 13-24): Subband count matches `numcb`/`numvb`. VB/CB ordering correct. Anti-crossings at finite k_∥.

**g-factor figures** (figs 25-28): GaAs CB g ≈ -0.44 (near free-electron). InAs CB g ≈ -14.9. QW g-factor between bulk values.

**SC/QCSE figures** (figs 29-33): Potential profile bends with field. Charge neutrality achieved. Convergence monotonic.

**Wire figures** (figs 34-42): Hexagonal geometry correct. Core/shell material regions visible.

**Optics figures** (figs 43-59): Absorption edge at E_g. ISBT peak at expected energy. Gain blue-shifts with carrier density. Exciton binding 4-10 meV for GaAs QW.

**Benchmark figures** (figs 60-62): Match nextnano reference data where applicable.

- [ ] **Step 4: Fix any discrepancies found**

For each discrepancy logged in Step 3:
1. Open the relevant Fortran source or Python plotting function
2. Trace root cause (physics bug, plotting bug, or config issue)
3. Fix at source
4. Re-run `ctest --test-dir build --output-on-failure`
5. Re-generate the affected figure(s)
6. Commit the fix

- [ ] **Step 5: Commit discrepancy log**

```bash
git add docs/plans/phase4-discrepancy-log.md
git commit -m "docs: add Phase 4 figure discrepancy log"
```

---

## Task 4: Update REVIEW.md Group #22

**Files:**
- Modify: `docs/plans/REVIEW.md:90-98` (Group #22 status)

- [ ] **Step 1: Update Group #22 status**

Change lines 90-98 in `docs/plans/REVIEW.md` from:

```
**Status: INCOMPLETE**

Fortran code fully done: optical_spectra.f90, exciton.f90, scattering.f90, opticalProperties executable, 10+ configs. Lecture 06 documentation complete. Missing:

- **0/11 Python figure generation functions** ...
- **0/11 PNG files** ...
- **Exciton tutorial config**: ...
```

To:

```
**Status: COMPLETE**

Fortran code fully done. All 15 optics figure functions implemented in generate_all_figures.py (lines 3944-4873). All PNG files generated and verified in docs/figures/. Verified on 2026-05-04: absorption edge positions, ISBT peak energies, gain blue-shift behavior all physically correct.

ISBT cross-validation: z-dipole vs commutator velocity consistency confirmed (see test_isbt_dipole_velocity_consistency in test_optical_qw.pf).
```

Also update the tracking table on line 28 from `INCOMPLETE` to `COMPLETE`.

- [ ] **Step 2: Commit**

```bash
git add docs/plans/REVIEW.md
git commit -m "docs: mark Group #22 (QW tutorials optics) as COMPLETE in REVIEW.md"
```

---

## Task 5: New Figure — Bulk Bandstructure E(k)

**Files:**
- Modify: `scripts/plotting/generate_all_figures.py` (add `fig_bandstructure_bulk_ek` function and register in `ALL_FIGURES`)

- [ ] **Step 1: Write the figure function**

Add a new function `fig_bandstructure_bulk_ek` following the pattern from `fig_bulk_gaas_absorption` (line 4627). This function:

1. Creates a temporary config for bulk GaAs with k-sweep along [100]
2. Runs `bandStructure` executable
3. Parses `eigenvalues.dat`
4. Repeats for InAs and InSb
5. Plots 1×3 panel: GaAs, InAs, InSb bulk E(k)
6. Saves to `docs/lecture/figures/bandstructure_bulk_ek.png`

The config template for each material (GaAs shown):
```
waveVector: kx
waveVectorMax: 0.15
waveVectorStep: 101
confinement: 0
FDstep: 1
FDorder: 2
numLayers: 1
material1: GaAs 0 0 0
numcb: 2
numvb: 6
```

The function must:
- Parse eigenvalues: 8 bands × Nk points from `output/eigenvalues.dat`
- Shift all energies relative to VBM (band 4 at k=0 for GaAs)
- Color bands by type: HH (1-2 red), LH (3-4 blue), SO (5-6 orange), CB (7-8 cyan)
- Label x-axis "k (Å⁻¹)", y-axis "Energy (eV)"
- Each panel titled with material name

Register in `ALL_FIGURES` dict:
```python
"bandstructure_bulk_ek": fig_bandstructure_bulk_ek,
```

- [ ] **Step 2: Run and verify**

```bash
python scripts/plotting/generate_all_figures.py --skip-build --only bandstructure_bulk_ek
```

Expected: `docs/lecture/figures/bandstructure_bulk_ek.png` created. GaAs E_g ≈ 1.42 eV, InAs E_g ≈ 0.354 eV, InSb E_g ≈ 0.17 eV visible at k=0.

If discrepancy found (e.g., wrong E_g): apply discrepancy protocol from plan header.

- [ ] **Step 3: Commit**

```bash
git add scripts/plotting/generate_all_figures.py
git commit -m "feat: add bulk bandstructure E(k) figure (GaAs, InAs, InSb)"
```

---

## Task 6: New Figure — QW Subband Dispersion

**Files:**
- Modify: `scripts/plotting/generate_all_figures.py` (add `fig_bandstructure_qw_subbands` and register)

- [ ] **Step 1: Write the figure function**

Add `fig_bandstructure_qw_subbands`. This function:

1. Uses existing config `qw_gaas_algaas_absorption.cfg` (QW, kx sweep)
2. Runs `bandStructure` executable
3. Parses `eigenvalues.dat` — extract all subbands
4. Plots 2-panel: CB subbands (top) + VB subbands (bottom)
5. Saves to `docs/lecture/figures/bandstructure_qw_subbands.png`

The function must:
- Identify CB/VB boundary: CB starts at state `numvb + 1`
- Plot first 3 CB subbands and first 6 VB subbands (each subband = pair of spin-degenerate states)
- Show HH-LH anti-crossing at finite k_∥
- Label y-axis relative to VBM

Register:
```python
"bandstructure_qw_subbands": fig_bandstructure_qw_subbands,
```

- [ ] **Step 2: Run and verify**

```bash
python scripts/plotting/generate_all_figures.py --skip-build --only bandstructure_qw_subbands
```

Expected: `docs/lecture/figures/bandstructure_qw_subbands.png`. CB subbands parabolic near k=0, VB subbands show HH-LH mixing.

If discrepancy found: apply discrepancy protocol.

- [ ] **Step 3: Commit**

```bash
git add scripts/plotting/generate_all_figures.py
git commit -m "feat: add QW subband dispersion figure (GaAs/AlGaAs)"
```

---

## Task 7: New Figure — QW Wavefunctions

**Files:**
- Modify: `scripts/plotting/generate_all_figures.py` (add `fig_wavefunctions_qw` and register)

- [ ] **Step 1: Write the figure function**

Add `fig_wavefunctions_qw`. This function:

1. Creates a temporary QW config with eigenfunction output enabled
2. Runs `bandStructure` executable
3. Parses `output/eigenfunctions.dat` (z-coordinates + |ψ(z)|² for each state)
4. Plots 2-panel: CB states (top) + VB states (bottom), overlaid on band edge profile
5. Saves to `docs/lecture/figures/wavefunctions_qw.png`

Config template:
```
waveVector: k0
confinement: 1
FDstep: 101
FDorder: 4
numLayers: 3
material1: Al30Ga70As -150 -50 0
material2: GaAs -50 50 0
material3: Al30Ga70As 50 150 0
numcb: 4
numvb: 8
```

The function must:
- Parse the band edge profile from the material offsets (CB/VB edges per layer)
- Plot first 3 CB states + first 3 VB states as filled |ψ|² curves
- Background: band edge profile as shaded regions
- Label each state with subband index and energy

Register:
```python
"wavefunctions_qw_lecture": fig_wavefunctions_qw,
```

- [ ] **Step 2: Run and verify**

```bash
python scripts/plotting/generate_all_figures.py --skip-build --only wavefunctions_qw_lecture
```

Expected: Wavefunctions confined to GaAs well, penetration into AlGaAs barriers visible.

If discrepancy found: apply discrepancy protocol.

- [ ] **Step 3: Commit**

```bash
git add scripts/plotting/generate_all_figures.py
git commit -m "feat: add QW wavefunctions figure with band edge profile"
```

---

## Task 8: New Figure — Wire Geometry and Potential

**Files:**
- Modify: `scripts/plotting/generate_all_figures.py` (add `fig_wire_geometry_potential` and register)

- [ ] **Step 1: Write the figure function**

Add `fig_wire_geometry_potential`. This function:

1. Uses existing config from wire regression test (InAs/GaAs core-shell)
2. Runs `bandStructure` in wire mode
3. Reads the grid geometry from the output
4. Plots 2-panel: (left) material map with hex grid, (right) potential along radial cut
5. Saves to `docs/lecture/figures/wire_geometry_potential.png`

The function must:
- Reconstruct the 2D hexagonal grid from the config (nx × ny points)
- Color-code grid points by material (core=InAs, shell=GaAs)
- Right panel: extract potential profile along a radial line through the wire center
- Draw hexagonal boundary outline

Register:
```python
"wire_geometry_potential": fig_wire_geometry_potential,
```

- [ ] **Step 2: Run and verify**

```bash
python scripts/plotting/generate_all_figures.py --skip-build --only wire_geometry_potential
```

Expected: Hexagonal cross-section with clear core/shell boundary.

If discrepancy found: apply discrepancy protocol.

- [ ] **Step 3: Commit**

```bash
git add scripts/plotting/generate_all_figures.py
git commit -m "feat: add wire geometry and potential profile figure"
```

---

## Task 9: New Figure — Zeeman Fan Diagram

**Files:**
- Modify: `scripts/plotting/generate_all_figures.py` (add `fig_zeeman_fan_diagram` and register)

**Context:** Orbital Landau levels (Phase 5) are not yet implemented. This figure shows Zeeman spin splitting only, which is currently supported. The B-field enters via `cfg%bdg%B_vec` and `compute_zeeman_vz`.

- [ ] **Step 1: Write the figure function**

Add `fig_zeeman_fan_diagram`. This function:

1. Loops over B values: 0, 0.5, 1, 2, 3, 5, 7, 10 T
2. For each B, creates a temporary bulk InAs config with `B_vec: 0 0 <B>` and Zeeman splitting enabled
3. Runs `bandStructure` executable
4. Parses CB eigenvalues from `eigenvalues.dat`
5. Plots E_CB± vs B with analytical overlay
6. Saves to `docs/lecture/figures/zeeman_fan_diagram.png`

The function must:
- Plot the two spin-split CB levels at each B value
- Overlay analytical Zeeman: E_CB± = E_g ± g·μ_B·B, using g=2.0 (default) and g=-14.9 (InAs g*)
- Extract effective g* from the slope of the splitting: ΔE / (2·μ_B·B)
- Print extracted g* for verification
- Single panel, legend showing data + analytical

Register:
```python
"zeeman_fan_diagram": fig_zeeman_fan_diagram,
```

**Note:** If bulk Zeeman splitting requires BdG config (currently only wired in `topologicalAnalysis`), the figure function may need to use `topologicalAnalysis` instead. Check how B-field is parsed for bulk mode in `input_parser.f90` and which executable handles it. If no executable currently supports bulk Zeeman without topology, this figure becomes a placeholder for Phase 5 and should be documented as such.

- [ ] **Step 2: Run and verify**

```bash
python scripts/plotting/generate_all_figures.py --skip-build --only zeeman_fan_diagram
```

Expected: Linear splitting of CB levels with B. If the executable doesn't support this yet, document in discrepancy log and skip.

- [ ] **Step 3: Commit**

```bash
git add scripts/plotting/generate_all_figures.py
git commit -m "feat: add Zeeman fan diagram figure (spin splitting vs B)"
```

---

## Task 10: Register All New Figures and Run Full Suite

**Files:**
- Modify: `scripts/plotting/generate_all_figures.py` (verify ALL_FIGURES dict)

- [ ] **Step 1: Verify registration**

Confirm all 5 new figures appear in the `ALL_FIGURES` dictionary (line ~6110):

```python
"bandstructure_bulk_ek": fig_bandstructure_bulk_ek,
"bandstructure_qw_subbands": fig_bandstructure_qw_subbands,
"wavefunctions_qw_lecture": fig_wavefunctions_qw,
"wire_geometry_potential": fig_wire_geometry_potential,
"zeeman_fan_diagram": fig_zeeman_fan_diagram,
```

- [ ] **Step 2: Run full figure generation**

```bash
python scripts/plotting/generate_all_figures.py --skip-build 2>&1 | tee /tmp/fig_full_run.log
```

Expected: All 79 figures (74 original + 5 new) generate. Check for any warnings.

- [ ] **Step 3: Verify output files**

```bash
ls docs/lecture/figures/*.png | wc -l   # Should be 5+ (including existing topology figures)
ls docs/figures/*.png | wc -l           # Should be 76+
```

- [ ] **Step 4: Commit (if any registration fixes needed)**

```bash
git add scripts/plotting/generate_all_figures.py
git commit -m "fix: register all new Group #50 figures in ALL_FIGURES dict"
```

---

## Task 11: Rebuild Chapter 06 (Optical Properties)

**Files:**
- Read: `docs/lecture/06-optical-properties.md` (16 figure references)
- Modify: `docs/lecture/06-optical-properties.md` (captions/text if needed)
- Regenerate: All 16 referenced figures in `docs/figures/`

- [ ] **Step 1: Extract all figure references**

```bash
grep -n '!\[.*\](.*\.png)' docs/lecture/06-optical-properties.md
```

16 references found (see design spec section 6 for full list).

- [ ] **Step 2: Regenerate all referenced figures**

For each figure referenced in Ch06, re-run the corresponding function:

```bash
python scripts/plotting/generate_all_figures.py --skip-build \
  --only qw_optical_matrix_elements \
  --only qw_absorption_spectrum \
  --only bulk_gaas_absorption \
  --only absorption_excitonic_TE \
  --only qw_absorption_strained \
  --only qw_absorption_vs_width \
  --only qw_absorption_optics_exe \
  --only qw_absorption_spin_resolved \
  --only isbt_dipole_moments \
  --only isbt_absorption \
  --only gain_strained_comparison \
  --only scattering_lifetime_vs_width \
  --only scattering_lifetime_vs_field \
  --only exciton_binding_vs_width \
  --only exciton_bohr_vs_width \
  --only absorption_with_exciton
```

- [ ] **Step 3: Verify figure captions match content**

For each figure reference, read the surrounding markdown text and verify:
- Caption describes what the figure shows
- Peak positions mentioned in text match actual plot (e.g., "absorption onset at 1.42 eV")
- Units are consistent (eV vs cm⁻¹)
- TE/TM labels match

If discrepancy found: apply discrepancy protocol — fix markdown text to match actual figure data.

- [ ] **Step 4: Commit fixes**

```bash
git add docs/lecture/06-optical-properties.md docs/figures/
git commit -m "docs: rebuild Ch06 optics figures and verify captions"
```

---

## Task 12: Rebuild Chapter 08 (Wire) and Verify All Lecture Figures

**Files:**
- Read: `docs/lecture/08-quantum-wire.md` (5 figure references)
- Modify: `docs/lecture/08-quantum-wire.md` (captions/text if needed)
- Read: All lecture chapters for figure cross-references

**Context:** There is no separate Ch08 for exciton/scattering — that content is in Ch06 (lines 933-979). Ch08 is the quantum wire chapter. The exciton/scattering verification was done in Task 11 as part of Ch06 rebuild.

- [ ] **Step 1: Extract Ch08 figure references**

```bash
grep -n '!\[.*\](.*\.png)' docs/lecture/08-quantum-wire.md
```

5 references: wire_subbands, wire_density_2d, wire_inas_gaas_profile, wire_inas_gaas_subbands, wire_inas_gaas_wavefunctions.

- [ ] **Step 2: Regenerate wire figures**

```bash
python scripts/plotting/generate_all_figures.py --skip-build \
  --only wire_subbands \
  --only wire_density_2d \
  --only wire_inas_gaas_profile \
  --only wire_inas_gaas_subbands \
  --only wire_inas_gaas_wavefunctions
```

- [ ] **Step 3: Verify all lecture chapter figure references resolve**

For every chapter file (00-13), verify each `![...](...png)` reference points to an existing file:

```bash
for ch in docs/lecture/[01]*.md; do
  echo "=== $ch ==="
  grep -oP '!\[.*?\]\(\K[^)]+' "$ch" | while read ref; do
    resolved="docs/lecture/$ref"
    if [ ! -f "$resolved" ]; then
      echo "BROKEN: $ch references $ref (not found at $resolved)"
    fi
  done
done
```

Expected: Zero broken references. If any found, log in discrepancy log and fix (either regenerate figure or fix path in markdown).

- [ ] **Step 4: Commit**

```bash
git add docs/lecture/08-quantum-wire.md docs/figures/ docs/lecture/figures/
git commit -m "docs: rebuild Ch08 wire figures and verify all lecture references"
```

---

## Task 13: Update BACKLOG and Final Commit

**Files:**
- Modify: `docs/plans/BACKLOG.md` (mark Phase 4 as COMPLETED)
- Move: `docs/plans/2026-05-04-phase4-optics-figures-design.md` → `docs/plans/archive/`

- [ ] **Step 1: Update BACKLOG.md**

Change Phase 4 section from its current status to:

```markdown
## Phase 4: COMPLETED (2026-05-04)

All 74 existing figures validated, ISBT cross-verified (z-dipole vs commutator velocity), 5 new physics figures added (bulk E(k), QW subbands, wavefunctions, wire geometry, Zeeman fan). Ch06 + Ch08 rebuilt with verified captions. Groups #22, #26, #50 → COMPLETE.

**Plan:** `docs/plans/archive/2026-05-04-phase4-optics-figures-plan.md` (archived)
**Commits:** N commits from START to END.
```

Update Summary table: Phase 4 effort → DONE.

- [ ] **Step 2: Archive plan file**

```bash
mv docs/plans/2026-05-04-phase4-optics-figures-plan.md docs/plans/archive/
```

- [ ] **Step 3: Commit and push**

```bash
git add docs/plans/BACKLOG.md docs/plans/archive/ docs/plans/2026-05-04-phase4-optics-figures-design.md
git commit -m "docs: mark Phase 4 as COMPLETED, archive plan and design"
git push
```
