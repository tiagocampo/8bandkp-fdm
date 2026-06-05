# Issue 5: Wire strain quantitative validation

**Type**: AFK
**Blocked by**: None — can start immediately
**User stories**: 19, 20, 21

## What to build

Replace the 60% tolerance wire strain profile test with quantitative band-edge shift and HH-LH splitting checks against analytical Bir-Pikus predictions. Create a new wide-core wire configuration and Python verification script.

### Part 1: Wide-core wire strain TOML config

Create `tests/regression/configs/wire_inas_gaas_strain_wide.toml`:

```toml
confinement = "wire"
FDorder = 2
fd_step = 1
which_band = 0
band_idx = 1

[wave_vector]
mode = "kz"
max = 0.01
nsteps = 2

[bands]
num_cb = 4
num_vb = 8

[wire]
nx = 40
ny = 40
dx = 5.0
dy = 5.0

[wire.geometry]
shape = "rectangle"
width = 200.0
height = 200.0

[[region]]
material = "GaAs"
inner = 100.0
outer = 200.0

[[region]]
material = "InAs"
inner = 0.0
outer = 100.0

[strain]
reference = "GaAs"
solver = "pardiso"
piezoelectric = false

[feast]
emin = -1.5
emax = 2.0
m0 = -1
```

This gives a 100Å InAs core in a 200Å GaAs domain — the interior region (r < 50Å) should have near-biaxial strain.

### Part 2: Quantitative verification script

Create `tests/integration/verify_strain_wire_quantitative.py`. The script:

1. **Run strained and unstrained simulations** using `run_exe` from `star_helpers.py`:
   - Run `bandStructure` with the wide-core config (strained)
   - Run `bandStructure` with a modified config where `[strain]` section is removed (unstrained reference)

2. **Extract band-edge shifts** from the strained vs unstrained eigenvalue output:
   - Parse eigenvalues from `output/` files (follow existing pattern from `verify_strain_rung5_bulk.py`)
   - Identify CB bottom (lowest CB eigenvalue), VB top (highest VB eigenvalue)
   - Compute shifts: delta_CB = CB_strained - CB_unstrained, delta_VB = VB_strained - VB_unstrained
   - Compute gap shift: delta_gap = delta_CB - delta_VB

3. **Compute analytical Bir-Pikus reference** using `bir_pikus_biaxial_001()` from `star_helpers.py`:
   - InAs parameters: a0 = 6.0583 Å, C11 = 832.9 kbar, C12 = 452.6 kbar, ac = -5.08 eV, av = 1.00 eV, b_dp = -1.80 eV (from `parameters.f90` Vurgaftman table)
   - Substrate: GaAs a_sub = 5.6533 Å
   - This gives analytical delta_Ec, delta_EHH, delta_ELH, HH_LH_splitting

4. **Assert quantitative agreement:**
   - |delta_CB_measured - delta_Ec_analytical| / |delta_Ec_analytical| < 0.10 (10%)
   - |delta_VB_measured - delta_EHH_analytical| / |delta_EHH_analytical| < 0.10 (10%)
   - |HH_LH_splitting_measured - HH_LH_splitting_analytical| / |HH_LH_splitting_analytical| < 0.10 (10%)

   Use 10% tolerance (not 1% as in bulk rung5) because the wire strain field differs from ideal biaxial due to free-surface relaxation, even at 100Å core. The interior is "mostly biaxial" but not perfectly so.

5. **Interior strain profile check** (improved from current 60%):
   - Parse the strain output file (ε_zz, ε_yy, ε_xx at interior grid points)
   - Compute average interior Tr(ε) for grid points with r < core_radius/2
   - Compare against analytical Tr(ε) = 2ε_xx + ε_zz (biaxial)
   - Assert within 10%

6. **HH-LH splitting from eigenvalues:**
   - From the unstrained simulation, extract the top-of-VB degenerate level (HH = LH at k=0 for unstrained zincblende)
   - From the strained simulation, extract the highest VB and next VB eigenvalues
   - Their splitting should match the Bir-Pikus HH-LH splitting

### CTest registration

Register as `strain_validation_wire_quantitative` with label `strain-validation;verification`.

## Acceptance criteria

- [ ] `tests/regression/configs/wire_inas_gaas_strain_wide.toml` created
- [ ] `tests/integration/verify_strain_wire_quantitative.py` created
- [ ] CB band-edge shift within 10% of Bir-Pikus analytical prediction
- [ ] VB band-edge shift within 10% of Bir-Pikus analytical prediction
- [ ] HH-LH splitting within 10% of Bir-Pikus analytical prediction
- [ ] Interior strain Tr(ε) within 10% of biaxial reference
- [ ] CTest registered under `strain-validation;verification` label
- [ ] All existing tests pass

## Blocked by

None — can start immediately.
