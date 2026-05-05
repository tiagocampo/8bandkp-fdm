# Magnetic Field Integration — Landau Levels Implementation

## Context

Goal: Enable magnetic field effects (Landau quantization + Zeeman splitting) across all three confinement modes.

**What's already done:**
- Wire (confinement=2): `add_zeeman_coo` + `add_peierls_coo` applied in `bdg_hamiltonian.f90` ✓
- Zeeman in bulk `ZB8bandBulk`: lines 438-460 in `hamiltonianConstructor.f90` ✓
- `add_peierls_coo` implemented in `magnetic_field.f90` ✓

**What's needed:**
1. QW (confinement=1): Add Zeeman to `ZB8bandQW` for B⊥
2. Bulk (confinement=0): Post-processed analytical Landau markers on computed band structure
3. Wire: Verify full Peierls works (not just Zeeman)

## Phase 1: Verify Wire Peierls

The `add_peierls_coo` in `magnetic_field.f90` is called from `bdg_hamiltonian.f90:149`. For **wire with confinement=2**, this applies Peierls to the CSR wire Hamiltonian.

### Verification

Run `tests/integration/test_bhz_z2.sh` with topological config to verify Z2 is still correct after the BHZ wire fix (Task 7). The BHZ wire fix used M-sign heuristic, not Peierls — verify nothing regressed.

```bash
bash tests/integration/test_bhz_z2.sh \
  build/src/topologicalAnalysis \
  tests/regression/configs/topology_bhz_z2_trivial.cfg \
  tests/regression/configs/topology_bhz_z2_topological.cfg
```

Expected: Z2=0 for trivial (d=58Å), Z2=1 for topological (d=70Å)

## Phase 2: QW Zeeman (confinement=1)

### Where to Add

`ZB8bandQW` in `hamiltonianConstructor.f90` — add Zeeman diagonal after Hamiltonian construction (after line ~430). Pattern identical to `ZB8bandBulk` (lines 438-460).

### Implementation

At the end of `ZB8bandQW`, after all HT matrix elements are set and before the optional `g` derivative mode:

```fortran
! ---------------------------------------------------------------
! Zeeman splitting: g*mu_B * B . sigma for B in z-direction
! Applied to diagonal blocks when B_z /= 0
! ---------------------------------------------------------------
if (present(cfg)) then
  if (cfg%bdg%enabled) then
    block
      real(kind=dp) :: mu_B, B_z, Vz(8)
      B_z = cfg%bdg%B_vec(3)
      if (abs(B_z) < 1e-12_dp) exit block
      mu_B = e * hbar / (2.0_dp * m0)
      ! g_J eigenvalues: HH=-1.5, LH=+0.5, SO=-0.5, CB=+1.0
      Vz(1:2) = -1.5_dp * cfg%bdg%g_factor * mu_B * B_z  ! HH
      Vz(3:4) =  0.5_dp * cfg%bdg%g_factor * mu_B * B_z  ! LH
      Vz(5:6) = -0.5_dp * cfg%bdg%g_factor * mu_B * B_z  ! SO
      Vz(7)   = -1.0_dp * cfg%bdg%g_factor * mu_B * B_z  ! CB1
      Vz(8)   =  1.0_dp * cfg%bdg%g_factor * mu_B * B_z  ! CB2
      ! Add to each 8x8 block diagonal (bands 1-8, repeated for each z grid point)
      do i = 1, N
        HT((i-1)*8+1:(i-1)*8+8, (i-1)*8+1:(i-1)*8+8) = &
          HT((i-1)*8+1:(i-1)*8+8, (i-1)*8+1:(i-1)*8+8) + diag(Vz)
      end do
    end block
  end if
end if
```

Note: Only `B_z` is used (perpendicular B for QW). In-plane B components would require full Peierls (not implemented for QW due to translational invariance constraint).

### Physical Justification

For a QW with strong z-confinement, B⊥ (z-direction) primarily causes Zeeman spin splitting, not orbital Landau quantization in the z direction. The orbital effect (cyclotron motion in x-y) requires Peierls with y-discretization, which is expensive. Zeeman alone captures the dominant physics for typical QW experiments (g-factor splitting is directly measurable).

### Test

Config `landau_InAs.cfg` with `confinement=1`, `material1: InAs`, `ExternalField: 5.0 EF` — verify CB Zeeman splitting of ~g·μ_B·B ≈ 2·0.058·5 ≈ 0.58 meV for InAs.

## Phase 3: Bulk Analytical Landau Markers (material-generic)

For **bulk (confinement=0)**, post-process computed band structure with analytical Landau levels. Material-generic: uses `meff` from `paramStruct` if available, otherwise computes from `γ1`.

### Formula

For a 2D electron gas in perpendicular B:
```
ω_c = e·B / m*          (cyclotron frequency)
E_n = E_CB + (n + ½)·ħ·ω_c
```

### CB Effective Mass

Priority:
1. **If `meff` is set** in `paramStruct` → use it directly
2. **Else compute from γ1** (for zinc-blende CB):
   ```
   m*_CB = m₀ / γ1
   ```
   This is the Kane model approximation for the conduction band.

### ħω_c Computation

```python
def compute_cyclotron_energy(m_eff, B, hbar=6.582119e-16, e=1.602e-19):
    """
    Compute cyclotron energy ħ·ω_c = ħ·e·B/m*
    Returns ħω_c in meV.
    m_eff: effective mass in units of m₀ (unitless)
    B: magnetic field in Tesla
    """
    # ħ in eV·s, e in C, result in eV
    hbar_eV = hbar * 1e-3  # convert to eV·s
    omega_c = e * B / (m_eff * 9.109e-31)  # m* = m_eff * m₀
    # Actually use: ħω = ħ · e·B/m*  where m* = m_eff · m₀
    # In eV: ħ (eV·s) · e (C) · B (T) / (m_eff · m₀)
    # Convert: ħω (eV) = ħ(eV·s) × e(C) × B(T) / (m_eff × m₀(kg))
    hbar_J = 1.054571817e-34  # J·s
    m_kg = m_eff * 9.109e-31  # kg
    omega_c = e * B / m_kg  # rad/s
    hbar_omega = hbar_J * omega_c  # J
    return hbar_omega * 6.241509e18  # convert J to eV
```

### Landau Level Indexing

For each band index n = 0, 1, 2, 3:
```
E_n = CB_edge + (n + 0.5) * hbar_omega_c
```

### Material Database

```python
MATERIAL_DB = {
    'InAs': {'meff': 0.026, 'Eg': 0.417},
    'GaAs': {'meff': 0.067, 'Eg': 1.519},
    'InP':  {'meff': 0.077, 'Eg': 1.4236},
    # ... fallback to gamma1 if meff not in params
}
```

### Test

Run `scripts/verify_landau_levels.py` with various materials (InAs, GaAs) and verify:
- InAs at B=5T: ħω_c ≈ 22.26 meV → E_0 ≈ 11.13 meV above CB
- GaAs at B=5T: ħω_c ≈ 8.67 meV → E_0 ≈ 4.34 meV above CB

## Files to Modify

| File | Change |
|------|--------|
| `src/physics/hamiltonianConstructor.f90` | Add Zeeman to `ZB8bandQW` |
| `scripts/verify_landau_levels.py` | Analytical Landau overlay |

## Verification

```bash
# QW Zeeman test
cd build && ctest -R "test_qw\|test_hamiltonian" --output-on-failure

# Full suite
cd build && ctest --output-on-failure
```

## Summary of Physics Coverage

| Mode | Orbital (Peierls) | Zeeman | Landau Levels |
|------|-------------------|--------|---------------|
| Wire (2D grid) | ✓ Full | ✓ | N/A |
| QW (1D grid z) | Not implemented | ✓ NEW | Analytical |
| Bulk (k-space) | Not applicable | ✓ Existing | Analytical |

Wire gets full orbital + Zeeman. QW gets Zeeman only (orbital requires y-grid). Bulk gets analytical post-process.

---

## Phase 4: BdG Majorana Phase Diagram

### Background

The sweep script `scripts/sweep_rashra_bdg.py` was created but not yet executed to generate the figure. The BdG mode had a parser issue: `bdg: T` format doesn't set `cfg%bdg%enabled = .true.` due to label-aware parsing. Use `bdg T` (space, no colon).

### Verification Steps

1. **Quick stability test**: Run a single BdG config to verify it completes without segfault/hang
   ```bash
   cd /tmp && mkdir -p test_bdg && cd test_bdg && mkdir output
   # Use bdg T not bdg: T
   cat > input.cfg << 'EOF'
   waveVector: kz
   confinement: 2
   wire_nx: 21; wire_ny: 21; wire_width: 50.0
   material1: GaAs
   numcb: 4; numvb: 4
   topology: T
   mode: bdg
   compute_z2: F
   compute_ldos: F
   bdg T
   mu: 0.0005
   delta_0: 0.0003
   g_factor: 2.0
   EOF
   /data/8bandkp-fdm/build/src/topologicalAnalysis
   # Check exit code and output/topology_result.dat
   ```

2. **Sweep B field**: Run `scripts/sweep_rashra_bdg.py` to generate phase diagram
   - Sweep B from 0 to 10 T in 0.5 T steps
   - Parse min gap from `topology_result.dat` for each B
   - Plot min gap vs B showing gap closure at B_crit

### Analytical B_crit

For s-wave pairing + Zeeman:
```
B_crit = sqrt(mu^2 + Delta^2) / (g * mu_B)
mu_B = 0.05788 meV/T
```

With mu=0.5 meV, Delta=0.3 meV, g=2.0: B_crit ≈ 5.0 T

### Phase Diagram Figure

`docs/lecture/figures/rashra_majorana_phase_diagram.png`:
- x-axis: B (T) from 0 to 10
- y-axis: min gap (meV)
- Vertical dashed line at B_crit
- Computed data points (circles)
- Expected gap closure at B_crit annotated

### Files Generated

| File | Description |
|------|-------------|
| `docs/lecture/figures/rashra_majorana_phase_diagram.png` | Phase diagram |
| `docs/lecture/figures/rashra_majorana_phase_diagram.txt` | Raw sweep data |

---

## Phase 5: QW Z2 via Fu-Kane

### Current State

`compute_z2_fukane` stub in `topological_analysis.f90:169-176` returns 0 and is **never called**. QW mode with Z2 currently prints "Z2 = 0 (placeholder)".

### Fu-Kane Method for Z2

For 2D systems (QW or bulk), the Z2 invariant is computed via parity products at Time-Reversal Invariant Momenta (TRIM):

```
(-1)^Z2 = prod_i prod_m xi_m(Gamma_i)
```

where:
- `Gamma_i` ∈ {Γ, M1, M2, M3} are the 4 TRIM points in 2D BZ
- `xi_m = +/-1` is the parity eigenvalue of occupied band m at that TRIM
- For zinc-blende (with inversion symmetry): parities are determined by band character at Γ

### QW Z2 Implementation

The Fu-Kane method for QW requires:
1. Compute eigenvalues at 4 TRIM k-points (Γ, M1, M2, M3) using `ZB8bandQW`
2. Determine band parities (+/-1) at each TRIM — requires analyzing eigenvector symmetry or using k.p selection rules
3. Compute parity product over occupied bands at each TRIM
4. Z2 = 1 if product is negative (odd number of band inversions)

### Complexity

The full Fu-Kane implementation is non-trivial:
- Needs eigenvector analysis to determine parity at each TRIM
- Requires 4 full Hamiltonian diagonalizations (one per TRIM)
- QW mode's primary use is optical spectra, not topological analysis

### Decision

**For QW mode, defer Fu-Kane implementation.** QW is primarily used for optical properties. Z2 topological classification for QW structures is an advanced feature not required for the current lecture scope. The stub remains as documentation of the intended interface.

If QW Z2 becomes necessary in the future, the approach would be:
1. Add `compute_z2_fukane_qw` that takes kpterms and profile
2. Loop over 4 TRIM points calling `ZB8bandQW` at each
3. Use eigenvector symmetry to determine parities
4. This is a significant addition — not a quick fix

### Alternative: Gap-Sweep for QW

QW (confinement=1) could use gap-sweep Z2 similar to wire, but this is also not currently implemented. The gap-sweep would scan for edge states in the QW band gap, similar to how wire uses `compute_z2_gap`.

---

## Files Modified Summary

| File | Change | Phase |
|------|--------|-------|
| `src/physics/hamiltonianConstructor.f90` | Add Zeeman to `ZB8bandQW` | 2 |
| `scripts/verify_landau_levels.py` | Material-generic analytical Landau overlay | 3 |
| `scripts/sweep_rashra_bdg.py` | Execute to generate figure | 4 |
| `docs/lecture/figures/rashra_majorana_phase_diagram.png` | Generated figure | 4 |

## Commit

```
feat: add Zeeman splitting to QW mode for B_perp
feat: add analytical Landau level markers to landau verification figure
feat: add BdG Majorana phase diagram generation
docs: document QW Z2 Fu-Kane as deferred (QW not used for topology)
```
