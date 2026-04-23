# Wire Code Benchmarking: nextnano and Literature Comparison

Date: 2026-04-22

## 1. Systems Studied

| System | Code | Cross-section | Grid | Material | Confinement |
|--------|------|---------------|------|----------|-------------|
| GaAs rectangular wire | ours | 63 x 63 A (6.3 x 6.3 nm) | 21 x 21, dx=dy=3 A | GaAs (single) | Hard-wall Dirichlet |
| InAs/GaAs core-shell | ours | 150 x 150 A, 80 A InAs core | 30 x 30, dx=dy=5 A | InAs core / GaAs shell | Band offset + strain |
| InSb wire g-factor | ours | 55 x 55 A (5.5 x 5.5 nm) | 11 x 11, dx=dy=5 A | InSbW (Winkler params) | Hard-wall Dirichlet |
| InAs rectangular wire | nextnano | 10 x 10 nm | 0.5 nm spacing | InAs (single) | Hard-wall (infinite barrier GaAs) |
| ZB InSb cylindrical NW | Faria Jr. et al. | R = 15-50 nm | plane-wave expansion | InSb (14-band Kane) | Hard-wall cylindrical |

## 2. GaAs Rectangular Wire: Confinement Gap Analysis

### 2.1 Our results (63 x 63 A, 21 x 21 grid)

At kz = 0:

| Quantity | Value |
|----------|-------|
| VB top | -0.885 eV |
| CB1 (character-selected edge) | 1.150 eV |
| CB2 | 1.150 eV |
| CB3 | 1.705 eV |
| CB4 | 1.708 eV |
| CB5 | 1.715 eV |
| Effective near-gap split | 2.035 eV |
| Total VB states | 534 (267 Kramers pairs) |
| Total CB states | 28 (14 Kramers pairs) |
| CB1-CB2 spacing | 60.6 meV |
| CB2-CB3 spacing | 62.1 meV |

All eigenvalues are Kramers-degenerate (pairs), consistent with time-reversal symmetry at kz = 0.

### 2.2 Particle-in-a-box estimates

For a 2D square infinite well with side L = 6.3 nm:

| Level | (nx, ny) | CB (m* = 0.067 m0) | HH (m* = 0.35 m0) |
|-------|----------|--------------------|--------------------|
| 1 | (1,1) | 141 meV | 27 meV |
| 2 | (1,2)=(2,1) | 707 meV | 135 meV |
| 3 | (2,2) | 1132 meV | 216 meV |

The current clean figure run reports a near-gap split of 2.035 eV after the 2D boundary-operator fix. This replaces the stale 818 meV value from the earlier figure-generation pass. Because the wire state classification is now based on band-character weights and still lacks a direct external wire benchmark, this number should be re-benchmarked before being used as a quantitative nextnano comparison.

The earlier CB1-CB2 spacing number is no longer carried as a validated benchmark in the lecture figures. The current figures prioritize robust near-gap state selection and avoid connecting untracked eigenvalue branches across $k_z$.

### 2.3 Comparison with nextnano InAs wire (10 x 10 nm)

nextnano computes a 10 x 10 nm InAs wire with single-band effective mass and 6x6 k.p:

**nextnano single-band results (InAs, 10 nm, m* = 0.023 m0):**

| State | Energy (meV) | PIB Estimate |
|-------|-------------|--------------|
| 1st (1,1) | -18.3 | 164 |
| 2nd/3rd (1,2)=(2,1) | -45.8 | 818 |
| 4th (2,2) | -73.3 | 1309 |

The nextnano energies are measured from the InAs VB edge (set to 0), so these are hole confinement energies (negative = deeper into the valence band). The degeneracy of the 2nd and 3rd states is a consequence of the square symmetry, which we also observe in our GaAs wire.

**nextnano 6x6 k.p results (InAs, 10 nm, magnetic field):**

| State | Energy (meV) |
|-------|-------------|
| 1st/2nd (Kramers pair) | -34.2 |
| 3rd/4th (Kramers pair) | -36.7 |
| 5th/6th | -55.7 |
| 7th/8th | -58.1 |

The 6x6 k.p lifts the degeneracy of the (1,2) and (2,1) states through the LH-HH mixing encoded in the off-diagonal Luttinger parameters. Our 8-band code produces the same qualitative effect: nearly-degenerate Kramers pairs with small splittings from the off-diagonal k.p terms.

**Key comparison point:** The old 61 meV CB-spacing statement came from a stale figure-generation pass and is no longer used as a validated benchmark. The current figure set instead reports the selected GaAs near-gap split and leaves detailed CB-spacing validation as follow-up work.

## 3. InAs/GaAs Core-Shell Wire: Strained Heterostructure

### 3.1 Our results (150 x 150 A, 80 A InAs core, 30 x 30 grid)

The InAs core is under compressive biaxial strain (matched to GaAs lattice constant). Strain is computed with the continuum elasticity solver using GaAs as the reference substrate.

| Quantity | Value |
|----------|-------|
| InAs core EC | -0.173 eV |
| InAs core EV | -0.590 eV |
| Strained InAs gap (core center) | 0.417 eV |
| GaAs shell EC | 0.719 eV |
| GaAs shell EV | -0.800 eV |
| CB offset (Delta_EC) | 0.892 eV |
| VB offset (Delta_EV) | 0.210 eV |
| Effective near-edge split | ~0.017 eV |
| CB-like states | 24 (`CB` character >= 0.5 in the subband run) |
| Wavefunction localization | Interface-weighted mixed states on this coarse strained grid |

Band alignment: Type-I by band-edge profile, but the current near-edge eigenstates are strongly mixed and interface-weighted. Do not use this coarse strained run as a validated nextnano nanowire benchmark.

### 3.2 nextnano core-shell nanowires

nextnano provides tutorials for GaAs/AlGaAs circular and hexagonal core-shell nanowires:

- GaAs core (5 nm radius) / Al0.33Ga0.67As shell (15 nm radius)
- Single-band effective mass calculation
- States confined in the GaAs core have s-like ground states
- Higher states may exceed the barrier and become unconfined

Our InAs/GaAs core-shell shows similar qualitative behavior:
- Single CB1 state strongly localized in the InAs core
- Deep CB offset (0.892 eV) creates tight electron confinement
- Shallower VB offset (0.210 eV) allows hole states to spread more

### 3.3 Strain effects on band edges

The compressive strain in the InAs core has three effects captured by our Bir-Pikus solver:

1. **Hydrostatic shift**: Both CB and VB shift, modifying the effective gap
2. **HH-LH splitting**: The shear strain splits HH and LH bands. Under compressive strain, HH moves up relative to LH (toward the gap)
3. **Gap modification**: Strained InAs gap increases from 0.354 eV (unstrained) to 0.417 eV (strained core), a +63 meV increase

These strain-induced modifications are essential for accurate core-shell wire calculations. Without strain, the InAs/GaAs offset would be different and the subband structure would change qualitatively.

## 4. InSb Wire g-Factor: Detailed Analysis

### 4.1 Our results (55 x 55 A, 11 x 11 grid, Winkler parameters)

| Component | Value | Notes |
|-----------|-------|-------|
| gx | -49.94 | Transverse (x), band-character-selected CB edge |
| gy | -50.05 | Transverse (y), band-character-selected CB edge |
| gz | -49.97 | Axial (z, free direction) |
| |g| (RMS) | ~86.6 | sqrt(gx^2 + gy^2 + gz^2) |
| Bulk InSb |g*| | ~51 | All directions |

### 4.2 Faria Junior et al. (PRB 97, 245402, 2018)

This paper (authored by Tiago Campos, Paulo E. Faria Junior, Martin Gmitra, Guilherme M. Sipahi, and Jaroslav Fabian) uses a **14-band Kane model** with plane-wave expansion for ZB InSb cylindrical nanowires oriented along [001].

Key findings from the paper:

| Diameter (nm) | g-factor (CB1) | Notes |
|---------------|---------------|-------|
| 30 | ~-25 | Strong confinement |
| 40 | ~-30 | |
| 50 | ~-35 | |
| 60 | ~-38 | |
| 80 | ~-40 | |
| 100 | ~-42 | Approaching bulk |
| Bulk | ~-51 | Roth formula |

The g-factor is **isotropic in the transverse plane** for cylindrical wires (gx = gy) due to rotational symmetry, but differs from gz (axial). The paper uses cylindrical geometry, not rectangular.

### 4.3 Comparison with our results

| Aspect | Our code | Faria Junior et al. | Agreement |
|--------|---------|-------------------|-----------|
| Material | InSbW (Winkler params) | InSb (14-band Kane) | Different parameter sets |
| Geometry | Rectangular 55x55 A | Cylindrical R = 15-50 nm | Different shape |
| Grid | 11x11 (very coarse) | Plane-wave (converged) | Coarse grid in ours |
| Diameter | 5.5 nm | 30-100 nm | Much thinner wire |
| Model | 8-band k.p + Lowdin | 14-band Kane | Fewer bands in ours |
| |g| (transverse) | ~50 | ~25-42 (30-100 nm) | Bulk-scale; not converged |
| gz (axial) | ~50 | Not separately reported | - |
| Bulk limit | ~51 | ~51 | Agreement |

**Analysis of discrepancies:**

1. **Wire diameter**: Our wire is 5.5 nm, much smaller than the 30-100 nm range studied by Faria Junior. For thin wires, the g-factor should generally move toward $g_{\rm free}=2$ as confinement quenches orbital mixing. The current coarse-grid regression remains close to the bulk value, so it should be treated as a state-selection sanity check rather than a converged confinement trend.

2. **Coarse grid**: Our 11x11 grid is extremely coarse. The g-factor is sensitive to wavefunction accuracy near the boundary. Convergence studies with finer grids (41x41 or higher) are needed.

3. **Rectangular vs cylindrical**: The rectangular cross-section breaks rotational symmetry, while the paper's cylindrical wires have $g_x = g_y$. The present coarse 11x11 regression is nearly isotropic, so it does not yet validate the expected rectangular anisotropy.

4. **8-band vs 14-band**: Our 8-band model may miss contributions from higher bands that affect the g-factor, especially in narrow-gap materials like InSb. The 14-band model includes explicit Dresselhaus BIA terms from higher conduction bands.

### 4.4 Literature comparison: InSb nanowire g-factors

| Source | System | Diameter | g-factor | Method |
|--------|--------|----------|----------|--------|
| **Our code** | InSb rect. wire | 5.5 nm | gx≈gy≈gz≈-50 | 8-band k.p + Lowdin |
| Faria Jr. et al. 2018 | InSb cyl. NW | 30 nm | ~-25 | 14-band Kane |
| Faria Jr. et al. 2018 | InSb cyl. NW | 100 nm | ~-42 | 14-band Kane |
| Nilsson et al. 2009 (Nano Lett.) | InSb NW QD | ~100 nm | up to ~70 | Transport measurement |
| Quay et al. 2005 (PRB) | InSb 2DEG | - | ~-30 | Magnetotransport |
| Roth formula | InSb bulk | inf | -51 | Analytical |

Nilsson et al. (Nano Lett. 9, 3151, 2009) report **giant, level-dependent g-factors** in InSb nanowire quantum dots with values up to ~70. This exceeds the bulk value because of orbital contributions from the quantum dot confinement geometry. Level-dependent fluctuations arise from spin-orbit interaction. The large values are consistent with the strong orbital magnetic response in InSb.

### 4.5 Assessment

Our $g \approx -50$ for a 5.5 nm InSb wire is:
- **Bulk-scale**: The sign and magnitude are now consistent with InSb material scale after band-character selection.
- **Not a confinement benchmark**: A 5.5 nm wire should require a convergence/radius study before any trend claim.
- **Still not a Faria Junior reproduction**: The comparison is complicated by different geometry, a very coarse grid, and 8-band vs 14-band Hamiltonians.

The coarse 11x11 grid is the primary concern. A convergence study with progressively finer grids is essential before drawing quantitative conclusions.

## 5. Subband Energies and Ordering

### 5.1 GaAs wire subband classification

Our GaAs wire at kz = 0 produces subbands in the following energy ranges:

| Band type | Energy range (eV) | Count (pairs) | Character |
|-----------|-------------------|---------------|-----------|
| Deep VB | -1.498 to -0.385 | ~130 pairs | Predominantly HH/LH mixed |
| Upper VB | -0.385 to 0.799 | ~135 pairs | Strongly mixed HH/LH/SO |
| **Near-gap split** | **-0.885 to 1.150** | **--** | **2.035 eV** |
| CB1/CB2 | 1.150 | 1 near-degenerate pair | centered square-wire ground density |
| CB3/CB4 | 1.705-1.708 | 1 near-degenerate pair | higher confined modes |
| CB5/CB6 | 1.715-1.716 | 1 near-degenerate pair | higher confined modes |
| Higher CB-like states | >1.7 | several | window-dependent |

The lowest CB-like pair is nearly degenerate at 1.150 eV, consistent with Kramers symmetry at $k_z=0$. The detailed excited-state ordering remains a follow-up benchmark item because the current wire output writer pads missing states at later $k_z$ values and does not yet provide robust branch-resolved character at every point.

### 5.2 nextnano InAs wire ordering comparison

The nextnano 10 nm InAs wire shows:

| Level | Single-band | 6x6 k.p | Notes |
|-------|------------|---------|-------|
| 1st | -18.3 meV | -34.2 meV | HH-like ground state |
| 2nd/3rd | -45.8 meV (degenerate) | -36.7 / -55.7 meV (split) | LH mixing splits degeneracy |

The same qualitative pattern is seen in our code: the multi-band k.p model lifts degeneracies that exist in the single-band approximation, producing a richer subband structure with non-uniform spacings.

## 6. Wavefunction Localization

### 6.1 Our results

From the band decomposition (parts.dat) for the InAs/GaAs core-shell wire:

| State | HH1 | HH2 | LH1 | LH2 | SO1 | SO2 | CB1 | CB2 |
|-------|-----|-----|-----|-----|-----|-----|-----|-----|
| VB1 | ~0% | 16.3% | ~0% | 44.5% | 37.9% | ~0% | ~0% | 1.4% |
| VB2 | 41.3% | ~0% | 31.5% | ~0% | ~0% | 27.1% | ~0% | ~0% |
| VB3 | 28.3% | ~0% | 32.2% | ~0% | ~0% | 38.7% | ~0% | ~0% |

The top VB states have mixed HH/LH/SO character, confirming the strong band mixing in the 8-band model. The CB1 state is not shown in the VB decomposition (it has vanishing VB character), confirming clear CB-VB separation.

### 6.2 Expected wavefunction patterns

For a rectangular wire with Dirichlet boundaries, the wavefunctions should be:

- **CB1**: s-like (ground state, no nodes), peaked at center
- **CB2/CB3**: p-like (one nodal line), oriented along x or y
- **CB4**: d-like (two nodal lines), approximately (2,2) mode
- **VB states**: More complex due to HH-LH mixing and the anisotropic effective masses

The nextnano tutorial expectation for simple square/cylindrical wells is a center-peaked ground state with higher nodal modes. Our homogeneous GaAs hard-wall wire now follows that qualitative pattern in the regenerated density figure. The strained InAs/GaAs heterostructure output does not: its selected near-edge states are mixed and interface-weighted, so it remains a diagnostic rather than a validated nextnano reproduction.

Do not use the current InAs/GaAs wavefunction panel as evidence of a clean centered CB ground state.

## 7. Summary Comparison Table

| Quantity | Our Code | nextnano | Literature (Faria Jr.) | Assessment |
|----------|---------|----------|----------------------|------------|
| **GaAs wire near-gap split (6.3 nm)** | 2.035 eV | -- | -- | Recomputed after boundary fix; needs external benchmark |
| **GaAs CB edge** | 1.150 eV | -- | -- | Current clean figure run |
| **GaAs Kramers degeneracy** | Yes (all kz=0 states) | Yes | Yes | Correct |
| **InAs/GaAs near-edge split (8 nm core)** | ~0.017 eV | -- | -- | Mixed/interface-weighted qualitative output |
| **InAs/GaAs CB-like states** | 24 | -- | -- | Not yet a validated nextnano reproduction |
| **InSb gx (5.5 nm)** | -49.94 | -- | ~-25 (30 nm) | Needs convergence study |
| **InSb gz (5.5 nm)** | -49.97 | -- | ~-25 (30 nm) | Needs convergence study |
| **InSb anisotropy** | weak in coarse regression | -- | gx = gy (cylindrical) | Not yet validated |
| **VB band mixing** | Strong HH/LH/SO mixing | 6x6 k.p shows similar | 14-band shows similar | Correct qualitative behavior |
| **Subband ordering** | Kramers pairs, square symmetry | Same pattern | -- | Consistent |

## 8. Recommended Next Steps

1. **Grid convergence for g-factor**: Run InSb wire g-factor on 21x21, 31x31, 41x41 grids to assess convergence. The 11x11 grid is too coarse for quantitative g-factor work.

2. **Diameter sweep**: Systematically vary the InSb wire cross-section from 5 nm to 50 nm to reproduce the Faria Junior diameter-dependent g-factor trend. Use the same cylindrical geometry if possible.

3. **14-band comparison**: Investigate whether the 8-band model systematically underestimates or overestimates the g-factor relative to the 14-band Kane model for narrow-gap materials.

4. **CB subband spacing validation**: Compare GaAs CB subband spacings with an independent effective-mass calculation using the energy-dependent effective mass from the 8-band dispersion.

5. **nextnano wire reproduction**: Attempt to reproduce the nextnano InAs 10x10 nm wire tutorial results directly, using the same material parameters and grid resolution, to establish a precise numerical benchmark.

## 9. References

1. P. E. Faria Junior, T. Campos, M. Gmitra, G. M. Sipahi, J. Fabian, "Spin-orbit coupling effects in zinc-blende InSb and wurtzite InAs nanowires: Realistic calculations with multiband k.p method," Phys. Rev. B 97, 245402 (2018). arXiv:1802.06734

2. H. A. Nilsson, P. Caroff, C. Thelander, M. Larsson, J. B. Wagner, L.-E. Wernersson, L. Samuelson, H. Q. Xu, "Giant, Level-Dependent g Factors in InSb Nanowire Quantum Dots," Nano Lett. 9, 3151 (2009).

3. nextnano3 tutorial: "2D QWR magnetic - InAs quantum wire" (2Dwire_vbmassesfrom6x6kp_iso)

4. nextnano3 tutorial: "2D Core-Shell nanowire" (2DGaAs_AlGaAs_circle)

5. I. Vurgaftman, J. R. Meyer, L. R. Ram-Mohan, "Band parameters for III-V compound semiconductors and their alloys," J. Appl. Phys. 89, 5815 (2001).

6. O. V. Kibis, "Hole subbands in freestanding nanowires: six-band versus eight-band k.p modelling," Phys. Status Solidi RRL 6, 155 (2012).

7. T. Campos, "Spin-orbit coupling effects and g-factors in zinc-blende InSb and wurtzite InAs nanowires," PhD thesis, Instituto de Fisica de Sao Carlos, Universidade de Sao Paulo (2018).
