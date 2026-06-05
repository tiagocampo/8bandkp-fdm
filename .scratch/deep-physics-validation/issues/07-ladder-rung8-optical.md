# Issue 7: Verification ladder rung 8 — optical observables

**Type**: AFK
**Blocked by**: None — can start immediately
**User stories**: 16, 17, 18

## What to build

Create `tests/integration/verify_8band_rung8_optical.py` — a verification ladder rung that tests the `opticalProperties` executable against analytical references for Kane matrix elements and absorption edge energies.

### Test structure

Follow the pattern from `verify_8band_rung1_bulk_k0.py`: standalone Python script, runs executable, parses output, compares against references, prints PASS/FAIL.

### R8.1: Bulk Kane matrix element Ep

The Kane interband matrix element Ep = (2m₀/ℏ²)|P|² is encoded in the k·p parameters. For GaAs, Vurgaftman 2001 gives Ep = 28.8 eV.

Verification approach:
1. Run `opticalProperties` with a bulk GaAs config that has `[optics]` section
2. Parse the absorption spectrum from `output/absorption_TE.dat` (or `absorption_TM.dat`)
3. The absorption spectrum near the band edge is proportional to the joint density of states × |P|²
4. Alternative (simpler): Run `bandStructure` with a bulk GaAs config at a small k value, extract the CB effective mass m* from a parabolic fit of E(k) near k=0. The 8-band Kane mass is: m*/m₀ = 1/(1 + Ep·(Eg+2ΔSO/3)/(Eg·(Eg+ΔSO))). Given m* and Eg and ΔSO, solve for Ep. Compare against Vurgaftman Ep.

Use the mass-based approach — it's simpler and doesn't require parsing optical spectra. The mass extraction is already done in `verify_8band_rung2_dispersion.py` (parabolic fit near k=0).

Actually, the simplest approach: compute Ep directly from the k·p parameters that are already in `star_helpers.py` or `parameters.f90`. The test is that the code's Ep (embedded in the Hamiltonian construction) produces the correct effective mass, which is already validated in rung 2. So for R8.1, we validate the Ep → m* → g* self-consistency chain:

1. Take Vurgaftman Ep for GaAs = 28.8 eV
2. Compute expected m* from Kane formula: 1/(1 + Ep·(Eg+2ΔSO/3)/(Eg·(Eg+ΔSO)))
3. Compute expected g* from Roth: g = 2 - 2·Ep·ΔSO/(3·Eg·(Eg+ΔSO))
4. Verify that the code's m* (from rung 2) and g* (from rung 7) are consistent with the same Ep

This is a consistency check between rungs 2 and 7, not a new executable run. It validates that the k·p parameter Ep is correctly embedded in both the dispersion (m*) and the magnetic (g*) response.

### R8.2: QW absorption edge

Run `opticalProperties` with a GaAs/AlGaAs QW config that has `[optics]` section (use existing `qw_gaas_algaas_optics.toml` or `qw_optics_commutator.toml`).

1. Parse `output/absorption_TE.dat` — two columns: energy (eV), absorption (cm⁻¹)
2. Find the absorption onset: the energy where absorption first exceeds a threshold (e.g., 10 cm⁻¹)
3. Compute expected onset from self-consistent ladder data: E_onset = Eg(GaAs) + E_CB1 + |E_VB1|
   - Eg = 1.519 eV (GaAs gap)
   - E_CB1 and E_VB1 from the QW subband energies validated in rung 3 (or compute them independently by running `bandStructure` with the same QW config)
4. Assert |E_onset_measured - E_onset_expected| < 0.010 eV (10 meV)

The 10 meV tolerance accounts for Lorentzian broadening (default linewidth = 30 meV FWHM in `optics_config%linewidth_lorentzian`), which smears the sharp absorption edge. The onset is shifted to lower energy by approximately one linewidth, but the exact shift depends on the lineshape.

### R8.3: TE/TM polarization ordering

A qualitative consistency check: for a QW, TE absorption (in-plane polarized, px+py) should onset before TM absorption (z-polarized, pz), because the QW confinement pushes the LH states down, making the HH-CB transition (TE) lower in energy than the LH-CB transition (TM).

1. Parse `output/absorption_TE.dat` and `output/absorption_TM.dat`
2. Find onset energies for both
3. Assert E_onset_TE < E_onset_TM (TE onset before TM)
4. Assert α_TE_peak > α_TM_peak near the edge (TE stronger)

This is already tested in `verify_qw_absorption_polarization.py`, but including it in the ladder gives diagnostic precision (if R8.2 passes but R8.3 fails, the bug is in the polarization decomposition, not the energy scale).

### Configs needed

- Bulk GaAs optics: use existing `bulk_gaas_optics.toml`
- QW GaAs/AlGaAs optics: use existing `qw_optics_commutator.toml`

No new configs needed for this issue.

### CTest registration

Register as `verification_rung8_optical` with label `verification`.

## Acceptance criteria

- [ ] `tests/integration/verify_8band_rung8_optical.py` created
- [ ] R8.1: Kane Ep self-consistency between m* and g* verified for GaAs
- [ ] R8.2: QW absorption onset within 10 meV of Eg + ECB1 + EVB1
- [ ] R8.3: TE onset < TM onset and TE peak > TM peak near edge
- [ ] Script prints PASS/FAIL per section
- [ ] CTest registered under `verification` label
- [ ] All existing tests pass

## Blocked by

None — can start immediately.
