# Issue 6: Verification ladder rung 7 — g-factor

**Type**: AFK
**Blocked by**: None — can start immediately
**User stories**: 14, 15

## What to build

Create `tests/integration/verify_8band_rung7_gfactor.py` — a verification ladder rung that tests the `gfactorCalculation` executable against the Roth analytical formula for multiple materials and geometries. This extends the existing 4-rung + 2 strain-rung ladder to cover g-factor observables.

### Test structure

Follow the pattern established by `verify_8band_rung1_bulk_k0.py` and `verify_8band_rung2_dispersion.py`: a standalone Python script that runs the Fortran executable, parses output, compares against analytical references, and prints PASS/FAIL per section.

### R7.1: Bulk CB g-factor vs Roth formula (4 materials)

For each material, run `gfactorCalculation` with a bulk config, extract g_z from `output/gfactor.dat`, compare against `roth_gfactor(Ep, Eg, DeltaSO)` from `star_helpers.py`.

| Material | Config | Ep (eV) | Eg (eV) | DeltaSO (eV) | Expected g* | Tolerance |
|----------|--------|---------|---------|-------------|-------------|-----------|
| GaAs | `gfactor_bulk_gaas_cb.toml` | 28.8 | 1.519 | 0.341 | ~-0.06 | 2% |
| InAs | `gfactor_bulk_inasw_cb.toml` | 21.5 | 0.354 | 0.38 | ~-14.9 | 2% |
| InSb | New config needed | 23.3 | 0.235 | 0.81 | ~-51 | 2% |
| GaSb | New config needed | 22.4 | 0.812 | 0.76 | ~-9.1 | 2% |

Materials parameters from Vurgaftman 2001 Table III (accessed via `star_helpers.py` constants or hardcoded).

The InSb config (`gfactor_bulk_insb_cb.toml`) and GaSb config (`gfactor_bulk_gasb_cb.toml`) need to be created. Use `gfactor_bulk_gaas_cb.toml` as template, change `confinement = "bulk"`, material name, `which_band = 0`, `band_idx = 1`.

**Tolerance:** 2% for bulk — the Roth formula is the 2-band limit of the 8-band model; the 8-band code should agree with it to ~1-2% for these materials (the deviation comes from remote band contributions not captured by the Roth formula).

### R7.2: QW g-factor vs Roth (1 material)

Run `gfactorCalculation` with `gfactor_qw_cb.toml` (GaAs/AlGaAs QW). Extract g_z. Compare against Roth formula with the QW effective band gap (Eg + confinement energies). The confinement shift reduces Eg and changes the Roth prediction.

The effective Eg for the QW is: Eg_eff = Eg(GaAs) + E_CB1 + |E_VB1|, where E_CB1 and E_VB1 are the confinement energies extracted from the eigenvalue output.

**Tolerance:** 5% for QW — the Roth formula is less accurate for confined structures due to the modified wavefunction and non-parabolicity effects.

### Parsing pattern

Follow the g-factor parsing pattern from `verify_star_inas_bulk.py`:
1. Run executable via `run_exe("gfactorCalculation", config_path, cwd=temp_dir)`
2. Parse `output/gfactor.dat` — contains gz value (column format, extract the numerical value)
3. Compare against analytical reference

### CTest registration

Register as `verification_rung7_gfactor` with label `verification;standard-star`.

## Acceptance criteria

- [ ] `tests/integration/verify_8band_rung7_gfactor.py` created
- [ ] `tests/regression/configs/gfactor_bulk_insb_cb.toml` created
- [ ] `tests/regression/configs/gfactor_bulk_gasb_cb.toml` created
- [ ] R7.1: Bulk CB g-factor within 2% of Roth for GaAs, InAs, InSb, GaSb
- [ ] R7.2: QW CB g-factor within 5% of Roth with confinement correction for GaAs/AlGaAs QW
- [ ] Script prints PASS/FAIL per section
- [ ] CTest registered under `verification;standard-star` label
- [ ] All existing tests pass

## Blocked by

None — can start immediately.
