# Chapter 02: Quantum Well Band Structure

In Chapter 01 we solved the 8-band k.p Hamiltonian for bulk semiconductors -- a
single, infinite, translationally invariant crystal. Real devices are built from
heterostructures: thin layers of different semiconductors stacked along a growth
direction. The simplest and most important heterostructure is the **quantum well**
(QW): a narrow layer of one material sandwiched between barriers of another, confining
carriers in one spatial dimension.

This chapter develops the theory of the quantum well Hamiltonian, explains how it is
implemented in the code, and walks through two computed examples in detail: a classic
type-I GaAs/AlGaAs well and a type-III broken-gap AlSbW/GaSbW/InAsW well.

## 1. Theory

### 1.1 From bulk to confined: the envelope function approximation

In Chapter 01 we studied the 8-band k.p Hamiltonian for bulk semiconductors, where
the crystal has full translational symmetry in all three spatial directions. The
eigenstates are Bloch waves characterized by a continuous wavevector
$\mathbf{k} = (k_x, k_y, k_z)$, and the Hamiltonian is an $8 \times 8$ complex
Hermitian matrix for each $\mathbf{k}$-point.

A semiconductor heterostructure breaks this translational symmetry. When a thin layer
of one material (e.g., InAs) is sandwiched between layers of another (e.g., AlSb),
the potential becomes a function of the growth direction $z$. The key insight of the
envelope function approximation is that we can still expand the wavefunction in the
periodic part of the Bloch functions $|u_n\rangle$ of the constituent materials, but
now the expansion coefficients become $z$-dependent:

$$\Psi(\mathbf{r}) = \sum_{n=1}^{8} F_n(z) \, e^{i(k_x x + k_y y)} \, |u_n\rangle$$

where $F_n(z)$ are the envelope functions. The in-plane wavevector
$\mathbf{k}_\parallel = (k_x, k_y)$ remains a good quantum number because
translational symmetry is preserved in the $x$-$y$ plane. Along $z$, however, $k_z$
is no longer conserved -- it is replaced by a differential operator.

### 1.2 Replacing $k_z$ with the derivative operator

The central step in going from the bulk 8-band Hamiltonian to the quantum well
Hamiltonian is the substitution:

$$k_z \longrightarrow -i \frac{d}{dz}$$

This transforms every occurrence of $k_z$ and $k_z^2$ in the bulk Hamiltonian into
differential operators acting on the envelope functions $F_n(z)$. The k.p terms that
were simple scalars in the bulk case now become operators. Specifically:

- Terms proportional to $k_z^2$ become second-derivative operators $d^2/dz^2$
- Terms proportional to $k_z$ become first-derivative operators $d/dz$
- Terms proportional to $k_x$ and $k_y$ remain as scalar multipliers (the in-plane
  directions are free)

For the zinc-blende 8-band basis, the principal k.p terms and their $z$-dependence
are:

| Symbol | Bulk form | QW form |
|--------|-----------|---------|
| $Q$ | $-(\gamma_1 + \gamma_2)(k_x^2 + k_y^2) - (\gamma_1 - 2\gamma_2)k_z^2$ | $-(\gamma_1 + \gamma_2)k_\parallel^2 - (\gamma_1 - 2\gamma_2)\,d^2/dz^2$ |
| $T$ | $-(\gamma_1 - \gamma_2)(k_x^2 + k_y^2) - (\gamma_1 + 2\gamma_2)k_z^2$ | $-(\gamma_1 - \gamma_2)k_\parallel^2 - (\gamma_1 + 2\gamma_2)\,d^2/dz^2$ |
| $S$ | $i2\sqrt{3}\,\gamma_3\, k_- k_z$ | $2\sqrt{3}\,\gamma_3\, k_-\, d/dz$ |
| $R$ | $-\sqrt{3}\bigl(\gamma_2(k_x^2 - k_y^2) - 2i\gamma_3 k_x k_y\bigr)$ | unchanged (no $k_z$ dependence) |
| $A$ | $A \, k^2$ | $A\bigl(k_\parallel^2 + d^2/dz^2\bigr)$ |
| $P_z$ | $P \, k_z$ | $-i P\, d/dz$ |

where $k_\pm = k_x \pm i k_y$ and the Luttinger parameters $\gamma_1, \gamma_2,
\gamma_3$, the interband momentum matrix element $P$, and the remote-band parameter
$A$ all become position-dependent functions $\gamma_1(z), \gamma_2(z), \ldots$ that
take the values of whichever material is present at position $z$.

After discretization on $N$ grid points with spacing $\Delta z$, the derivative
operators become FD matrices. For second-order accuracy (`FDorder = 2`):

$$D_2 = \frac{1}{\Delta z^2}
\begin{pmatrix}
-2 & 1 & & \\
1 & -2 & 1 & \\
 & \ddots & \ddots & \ddots \\
 & & 1 & -2
\end{pmatrix}, \qquad
D_1 = \frac{1}{2\,\Delta z}
\begin{pmatrix}
0 & 1 & & \\
-1 & 0 & 1 & \\
 & \ddots & \ddots & \ddots \\
 & & -1 & 0
\end{pmatrix}$$

Variable coefficients $g(z)$ are applied as $G \cdot D_2$ where $G =
\mathrm{diag}(g_1, \ldots, g_N)$, yielding tridiagonal coupling between nearest-
neighbor grid points. Higher FD orders produce banded matrices with wider stencils
(see Section 4.1).

### 1.3 Block matrix structure: $8 \times 8$ blocks, each $N \times N$

After discretizing the $z$-direction on a grid of $N$ points (the `FDstep`
parameter), the quantum well Hamiltonian becomes a large $8N \times 8N$ complex
Hermitian matrix. Its natural structure is an $8 \times 8$ block matrix, where each
block $(\alpha, \beta)$ is itself an $N \times N$ matrix:

$$H_{\text{QW}} = \begin{pmatrix}
H_{11} & H_{12} & \cdots & H_{18} \\
H_{21} & H_{22} & \cdots & H_{28} \\
\vdots & & \ddots & \vdots \\
H_{81} & H_{82} & \cdots & H_{88}
\end{pmatrix}$$

The band index pairs $(\alpha, \beta) \in \{1,\ldots,8\}$ follow the standard basis
ordering:
- Bands 1--4: heavy-hole (HH), light-hole (LH), LH', split-off (SO) valence bands
- Bands 5--6: split-off doublet
- Bands 7--8: conduction band (CB) doublet

Each off-diagonal block $H_{\alpha\beta}$ contains the k.p coupling terms (Q, R, S,
T, P operators) connecting bands $\alpha$ and $\beta$, discretized on the $N$-point
grid. Each diagonal block $H_{\alpha\alpha}$ contains the self-energy of band
$\alpha$, which includes:

1. The kinetic energy (k.p terms involving $\gamma_i$ or $A$)
2. The **band offset** $V_\alpha(z)$ from the `profile` array

A schematic of the block layout (for $N$ grid points per band) is:

```
       HH    LH    LH'   SO    SO'   SO''  CB    CB'
      ┌─────┬─────┬─────┬─────┬─────┬─────┬─────┬─────┐
HH    │  Q  │ SC  │     │     │     │     │ iP+ │     │
      ├─────┼─────┼─────┼─────┼─────┼─────┼─────┼─────┤
LH    │     │  T  │  R  │ -S  │     │     │     │ iP+ │
      ├─────┼─────┼─────┼─────┼─────┼─────┼─────┼─────┤
LH'   │     │     │  T  │  C  │     │     │     │     │
      ├─────┼─────┼─────┼─────┼─────┼─────┼─────┼─────┤
SO    │     │     │     │  T  │     │     │ -Pz │     │
      ├─────┼─────┼─────┼─────┼─────┼─────┼─────┼─────┤
SO'   │     │     │     │     │ ... │     │     │     │
      ├─────┼─────┼─────┼─────┼─────┼─────┼─────┼─────┤
CB    │     │     │     │     │     │     │  A  │     │
      ├─────┼─────┼─────┼─────┼─────┼─────┼─────┼─────┤
CB'   │     │     │     │     │     │     │     │  A  │
      └─────┴─────┴─────┴─────┴─────┴─────┴─────┴─────┘
  Each cell = NxN matrix (N = FDstep grid points)
```

### 1.4 The kpterms array: precomputed position-dependent operators

The code precomputes a three-dimensional array `kpterms(N, N, 10)` during the
initialization phase (`confinementInitialization`). The ten terms encode:

| Index | Content | Derivative order |
|-------|---------|-----------------|
| 1 | $\gamma_1(z)$ (diagonal) | 0 |
| 2 | $\gamma_2(z)$ (diagonal) | 0 |
| 3 | $\gamma_3(z)$ (diagonal) | 0 |
| 4 | $P(z)$ (diagonal) | 0 |
| 5 | $A(z) \cdot d^2/dz^2$ | 2nd |
| 6 | $P(z) \cdot d/dz$ | 1st |
| 7 | $(\gamma_1 - 2\gamma_2)(z) \cdot d^2/dz^2$ (Q-term) | 2nd |
| 8 | $(\gamma_1 + 2\gamma_2)(z) \cdot d^2/dz^2$ (T-term) | 2nd |
| 9 | $\gamma_3(z) \cdot d/dz$ (S-term) | 1st |
| 10 | $A(z)$ (diagonal) | 0 |

Terms 1--4 and 10 are purely diagonal: they contain the material parameters evaluated
at each grid point. Terms 5--9 couple neighboring grid points through the finite
difference stencil. The key design principle is that the position dependence of the
material parameters is baked into these operators once, so that the Hamiltonian
assembly at each $\mathbf{k}_\parallel$ point is a fast linear combination of
precomputed matrices.

### 1.5 Band offsets and heterostructure alignment

The heterostructure potential is stored in the `profile(N, 3)` array:

- `profile(:, 1)` = $E_V(z)$: valence band edge (applied to bands 1--4)
- `profile(:, 2)` = $E_V(z) - \Delta_{\text{SO}}(z)$: split-off band edge (applied
  to bands 5--6)
- `profile(:, 3)` = $E_C(z)$: conduction band edge (applied to bands 7--8)

These profiles are built by assigning the band-edge parameters of each material to
the grid points within that material's spatial range. The parameters $E_V$ and $E_C$
are stored in the material database (`parameters.f90`) and define the relative
alignment of bands across the heterostructure.

For an electric field along $z$ with magnitude $\mathcal{E}$, the code adds a linear
potential:

$$V_{\text{elec}}(z) = -\mathcal{E} \cdot L \cdot \frac{z + z_{\min}}{2\,z_{\min}}$$

where $L$ is the total simulation size and $z_{\min}$ is the leftmost grid coordinate.
This tilts the band edges, breaking the structural inversion symmetry.

The relative position of conduction and valence band edges across a heterointerface
determines the **band alignment type**:

| Type | Name | Alignment | Confinement |
|------|------|-----------|-------------|
| I | Straddling gap | $E_g(\text{well}) < E_g(\text{barrier})$, both CB and VB of well inside barrier gap | Electrons and holes in same layer |
| II | Staggered | CB of one material below CB of other but above its VB | Electrons and holes in different layers |
| III | Broken gap | $E_C(\text{material A}) < E_V(\text{material B})$ | Strong interband coupling |

**Type-I (straddling gap):** The band gap of the well material sits entirely within
the gap of the barrier. Both electrons and holes are confined in the same spatial
region. Example: GaAs quantum well in AlGaAs barriers.

**Type-II (staggered):** The conduction band of one material is below that of the
other, but above the valence band. Electrons and holes are spatially separated.

**Type-III (broken gap):** The conduction band of one material lies *below* the
valence band of the other. The canonical example is the InAs/GaSb system:

$$E_C(\text{InAs}) < E_V(\text{GaSb})$$

In a GaSb/InAs/GaSb quantum well, electrons accumulate in the InAs layer while holes
reside in the GaSb layers, creating natural spatial separation. Under an applied
electric field, the hybridization of these electron and hole states can be tuned
continuously, giving rise to a topological phase transition.

### 1.6 Analytical reference solutions

Before trusting numerical results, it is useful to compare against simple analytical
formulas that provide order-of-magnitude estimates.

**Infinite square well.** For a particle of effective mass $m^*$ confined in an
infinite well of width $L_w$, the energy levels are:

$$E_n = \frac{\hbar^2 \pi^2 n^2}{2\,m^*\,L_w^2}, \qquad n = 1, 2, 3, \ldots$$

For an electron in GaAs ($m^* = 0.067\,m_0$) in a 100 A well:

$$E_1 = \frac{(1.055 \times 10^{-34})^2 \cdot \pi^2}{2 \times 0.067 \times 9.109
\times 10^{-31} \times (100 \times 10^{-10})^2} \approx 56 \text{ meV}$$

This sets the scale: confinement energies in typical semiconductor QWs are tens to
hundreds of meV.

**Bastard finite well formula.** For a finite well of depth $V_0$ and width $L_w$,
the ground-state energy $E_1$ satisfies the transcendental equation (Bastard 1981):

$$\sqrt{V_0 - E_1}\,\tan\!\left(\frac{L_w}{2}\sqrt{\frac{2m_w^* E_1}{\hbar^2}}\right)
= \sqrt{\frac{m_b^*}{m_w^*}} \,\sqrt{E_1}$$

where $m_w^*$ and $m_b^*$ are the effective masses in the well and barrier
respectively. This formula is a useful sanity check: for a GaAs/Al$_{0.3}$Ga$_{0.7}$As
well with $V_0 = 458$ meV (the CB offset) and $L_w = 100$ A, it predicts $E_1
\approx 30$--$35$ meV above the GaAs CB edge. The 8-band k.p calculation will
deviate from this because it includes band mixing, nonparabolicity, and the
correct multi-band coupling.

---

## 2. In the Code

### 2.1 Initialization: `confinementInitialization`

When `confinement = 1` (QW mode), the input parser triggers
`confinementInitialization` (in `hamiltonianConstructor.f90`). This routine:

1. **Builds the z-grid:** From `startPos` and `endPos` of the first layer, the total
   size $L$ is computed. The grid spacing is $\Delta z = L / (N - 1)$, and the
   coordinate array is `z(i) = startPos(1) + (i-1) * delta`.

2. **Fills the profile array:** For each layer $i$, grid points in the range
   `[intStartPos(i) : intEndPos(i)]` receive:
   - `profile(:, 1) = params(i)%EV`
   - `profile(:, 2) = params(i)%EV - params(i)%DeltaSO`
   - `profile(:, 3) = params(i)%EC`

3. **Fills the kpterms array:** Material parameters are extracted per grid point and
   combined with FD stencil matrices. For order 2, this uses a forward/central/
   backward stencil decomposition to build tridiagonal operators. For higher orders,
   the code calls `buildFD2ndDerivMatrix` and `buildFD1stDerivMatrix` from
   `finitedifferences.f90`, then applies `applyVariableCoeff` to multiply each FD
   matrix by the position-dependent parameter profile.

4. **Applies electric field:** If `ExternalField = 1` with type `EF`, the routine
   `externalFieldSetup_electricField` adds a linear tilt to the profile:
   `profile(i,:) -= (Evalue * totalSize) * (z(i) + z(1)) / (2 * z(1))`.

### 2.2 Hamiltonian assembly: `ZB8bandQW`

For each $\mathbf{k}_\parallel$ point in the wavevector sweep, the routine
`ZB8bandQW` (in `hamiltonianConstructor.f90`) assembles the full $8N \times 8N$
Hamiltonian:

1. **Compute k.p blocks:** The precomputed `kpterms` are combined with the in-plane
   wavevector components $(k_x, k_y)$ to build the ten $N \times N$ blocks: `Q`,
   `T`, `S`, `SC`, `R`, `RC`, `PP`, `PM`, `PZ`, `A`.

   For example, the Q block is:
   ```
   Q(ii,jj) = -((kpterms(ii,jj,1) + kpterms(ii,jj,2)) * k_par^2
                 + kpterms(ii,jj,7))
   ```
   where `kpterms(:,:,1) = gamma1`, `kpterms(:,:,2) = gamma2`, and
   `kpterms(:,:,7)` already contains the discretized
   $-(\gamma_1 - 2\gamma_2) \cdot d^2/dz^2$ operator.

2. **Populate the 8x8 block matrix:** The code fills the Hamiltonian using Fortran
   array sections:
   ```fortran
   HT(1 + 0*N : 1*N, 1 + 0*N : 1*N) = Q          ! (1,1)
   HT(1 + 0*N : 1*N, 1 + 1*N : 2*N) = SC         ! (1,2)
   HT(1 + 0*N : 1*N, 1 + 6*N : 7*N) = IU * PP    ! (1,7)
   ...
   ```
   Each block assignment inserts the full $N \times N$ matrix of k.p couplings.

3. **Add band offsets:** The profile is added to the diagonal:
   ```fortran
   HT(ii, ii)     = HT(ii, ii) + profile(ii, 1)   ! bands 1-4: EV
   HT(N+ii, N+ii) = HT(N+ii, N+ii) + profile(ii, 1)
   ...
   HT(6*N+ii, 6*N+ii) = HT(6*N+ii, 6*N+ii) + profile(ii, 3)  ! bands 7-8: EC
   ```

4. **Diagonalize:** The band structure executable calls LAPACK's `zheevx` to find
   the requested number of eigenvalues (`numcb` conduction + `numvb` valence) at
   each k-point.

### 2.3 Step-by-step kpterms construction for a 3-layer system

Consider a three-layer QW: barrier/well/barrier with $N = 7$ grid points for
illustration. Points 1--2 are barrier, 3--5 are well, 6--7 are barrier. The second
derivative operator with FD order 2 gives a tridiagonal stencil with $-2/\Delta
z^2$ on the diagonal and $1/\Delta z^2$ on the sub/super-diagonals. For the Q-term
(index 7), which carries the coefficient $(\gamma_1 - 2\gamma_2)(z)$:

```
gamma1 - 2*gamma2 at each grid point (illustrative):
  z:    1      2      3      4      5      6      7
  g:  4.00   4.00   2.88   2.88   2.88   4.00   4.00

kpterms(:,:,7) = diag(g) @ D2  (second-derivative stencil)

Result: each row i gets the stencil weighted by g(i)
  Row 3 (well): g(3)*[-1, +2, -1] / dz^2 = 2.88*[-1, +2, -1] / dz^2
  Row 5 (well): g(5)*[-1, +2, -1] / dz^2 = 2.88*[-1, +2, -1] / dz^2
  Row 2 (barrier/well interface): g(2)*[-1, +2, -1] / dz^2 = 4.00*[...]
```

This is the essence: the position-dependent material parameter is multiplied into
the stencil at each grid point, so that the kinetic energy operator automatically
reflects the material composition at each location. When `confinementInitialization`
has finished, all ten kpterms matrices are ready, and the Hamiltonian assembly loop
over $\mathbf{k}_\parallel$ values needs only to form the linear combinations of
these precomputed blocks.

### 2.4 Input parsing for QW mode

The input parser (`input_parser.f90`) handles QW setup when `confinement = 1`:

- Reads `numLayers` material specifications, each with a name, start position, and
  end position (in Angstroms)
- Computes integer grid indices `intStartPos(i)` and `intEndPos(i)` that map each
  material to a contiguous range of FD grid points
- Stores the grid in `cfg%z(:)` and copies it to `cfg%grid%z(:)` via
  `init_grid_from_config`
- Calls `paramDatabase` to fill `cfg%params(:)` with material parameters from the
  database

---

## 3. Computed Examples

We present two full examples that illustrate contrasting physics: a classic type-I
GaAs/AlGaAs quantum well and a type-III broken-gap AlSbW/GaSbW/InAsW system.

### Example A: GaAs/Al$_{0.3}$Ga$_{0.7}$As Type-I Quantum Well

#### A.1 Configuration

This example is taken from `docs/benchmarks/qw_gaas_algaas.cfg`:

```
waveVector: kx
waveVectorMax: 0.1
waveVectorStep: 21
confinement:  1
FDstep: 201
FDorder: 2
numLayers:  3
material1: Al30Ga70As -200 200 0
material2: GaAs -50 50 0
material3: Al30Ga70As -200 200 0
numcb: 4
numvb: 8
ExternalField: 0  EF
EFParams: 0.0
```

#### A.2 Structure walkthrough

This defines a symmetric type-I quantum well:

| Layer | Material | Range (A) | $E_V$ (eV) | $E_C$ (eV) | $E_g$ (eV) |
|-------|----------|-----------|------------|------------|------------|
| 1 | Al$_{0.3}$Ga$_{0.7}$As | $[-200, 200]$ | $-0.959$ | $+1.018$ | 1.977 |
| 2 | GaAs | $[-50, 50]$ | $-0.800$ | $+0.719$ | 1.519 |
| 3 | Al$_{0.3}$Ga$_{0.7}$As | $[-200, 200]$ | $-0.959$ | $+1.018$ | 1.977 |

- **Domain:** $z \in [-200, 200]$ A, total 400 A
- **Well width:** 100 A of GaAs
- **Grid:** $N = 201$ points, $\Delta z = 400/200 = 2.0$ A
- **Band offsets:** $\Delta E_C = 1.018 - 0.719 = 0.299$ eV (CB),
  $\Delta E_V = 0.959 - 0.800 = 0.159$ eV (VB)
- **Sweep:** 21 k-points from $k_x = 0$ to $0.1$ A$^{-1}$
- **Spectrum:** 4 CB states + 8 VB states requested

This is the textbook type-I alignment: both the conduction band and valence band
edges of GaAs lie inside the gap of AlGaAs, so electrons are confined in the GaAs
layer by the CB offset of 299 meV, and holes are confined by the VB offset of
159 meV. The barrier layers overlap the full domain, providing a uniform background
outside the well.

#### A.3 Numerical results

At $k_\parallel = 0$, the code computes the following eigenvalues:

**Valence subbands (8 requested):**

| State | Energy (eV) | Character |
|-------|-------------|-----------|
| VB-8 | $-0.961$ | HH1 (second well state) |
| VB-7 | $-0.961$ | HH1 (degenerate partner) |
| VB-6 | $-0.960$ | LH1 |
| VB-5 | $-0.960$ | LH1 (degenerate partner) |
| VB-4 | $-0.960$ | HH2 |
| VB-3 | $-0.960$ | HH2 (degenerate partner) |
| VB-2 | $-0.959$ | SO-related |
| VB-1 | $-0.959$ | SO-related |

The VB top lies at approximately $-0.959$ eV, just above the GaAs valence band edge
($-0.800$ eV) plus the confinement shift. The energies are close to the AlGaAs
barrier edge because 8 states are requested and the GaAs well is narrow -- the
higher VB states are weakly confined and approach the barrier edge.

**Conduction subbands (4 requested):**

| State | Energy (eV) | Character |
|-------|-------------|-----------|
| CB-1 | $1.021$ | CB1 (ground state) |
| CB-2 | $1.021$ | CB1 (degenerate partner) |
| CB-3 | $1.031$ | CB2 (first excited state) |
| CB-4 | $1.031$ | CB2 (degenerate partner) |

The CB ground state at $1.021$ eV lies $1.021 - 0.719 = 0.302$ eV above the GaAs CB
edge, corresponding to the confinement energy of the electron in the 100 A well.
The degeneracy of 2 reflects the Kramers degeneracy (time-reversal symmetry at
$k_\parallel = 0$).

#### A.4 Comparison with Bastard formula

For the GaAs/AlGaAs well:
- $V_0 = 299$ meV (CB offset), $L_w = 100$ A
- $m_w^* = 0.067\,m_0$ (GaAs), $m_b^* = 0.093\,m_0$ (Al$_{0.3}$Ga$_{0.7}$As)

The Bastard transcendental equation predicts $E_1 \approx 35$--$40$ meV above the
GaAs CB edge. The 8-band k.p result shows the CB1 state at 302 meV above the GaAs
edge. The difference arises because the "simple" offset picture is complicated by
band mixing in the 8-band model: the CB state includes admixtures of VB/SO states
through the k.p couplings, which renormalize the effective confinement. Additionally,
the AlGaAs barrier CB edge at 1.018 eV is very close to the computed CB1 energy of
1.021 eV, indicating that the ground state is actually quite weakly confined and
extends significantly into the barrier -- a regime where the Bastard formula's
single-band approximation breaks down.

![GaAs/AlGaAs subband dispersion](../figures/qw_gaas_algaas_subbands.png)

*Figure 1: Subband dispersion $E(k_\parallel)$ for the GaAs/Al$_{0.3}$Ga$_{0.7}$As
quantum well. The type-I alignment confines both electrons (upper set) and holes
(lower set) in the GaAs layer. The near-parabolic CB dispersion is characteristic
of the light GaAs electron mass. The VB subbands show strong nonparabolicity due to
HH-LH mixing at finite $k_\parallel$.*

---

### Example B: AlSbW/GaSbW/InAsW Type-III Broken-Gap Quantum Well

#### B.1 Configuration

This example is taken from the regression test config
`tests/regression/configs/qw_alsbw_gasbw_inasw.cfg`:

```
waveVector: kx
waveVectorMax: 0.1
waveVectorStep: 11
confinement:  1
FDstep: 101
FDorder: 2
numLayers: 3
material1: AlSbW -250  250 0
material2: GaSbW -135  135 0.2414
material3: InAsW  -35   35 -0.0914
numcb: 32
numvb: 32
ExternalField: 0  EF
EFParams: 0.0005
```

#### B.2 Structure walkthrough

This defines a type-III broken-gap quantum well with three materials using Winkler
parameter sets (the "W" suffix denotes Winkler 2003 parameters with InSb as the
valence band energy reference):

| Layer | Material | Range (A) | $E_V$ (eV) | $E_C$ (eV) | $E_g$ (eV) | $\Delta_{\text{SO}}$ (eV) |
|-------|----------|-----------|------------|------------|------------|---------------------------|
| 1 | AlSbW | $[-250, 250]$ | $-0.410$ | $+1.974$ | 2.384 | 0.673 |
| 2 | GaSbW | $[-135, 135]$ | $-0.030$ | $+0.782$ | 0.812 | 0.760 |
| 3 | InAsW | $[-35, 35]$ | $-0.590$ | $-0.172$ | 0.418 | 0.380 |

Key observations:

- **Domain:** $z \in [-250, 250]$ A, total 500 A
- **GaSbW layer:** width 270 A, centered. Its $E_V = -0.030$ eV is much higher than
  the AlSbW $E_V = -0.410$ eV, creating a **deep well for holes** with $\Delta E_V =
  0.380$ eV.
- **InAsW layer:** width 70 A, nested inside GaSbW. Its $E_C = -0.172$ eV is far
  below the AlSbW $E_C = 1.974$ eV, creating a **deep well for electrons** with
  $\Delta E_C = 2.146$ eV.
- **Broken gap:** $E_C(\text{InAsW}) = -0.172$ eV lies above
  $E_V(\text{GaSbW}) = -0.030$ eV by 0.142 eV. This is a **broken-gap** alignment:
  the InAs conduction band sits above the GaSb valence band, causing strong
  hybridization between InAs electron states and GaSb hole states.
- **Grid:** $N = 101$ points, $\Delta z = 500/100 = 5.0$ A
- **Spectrum:** 32 CB + 32 VB states -- the large number reflects the deep wells
  that host many bound states

![QW band-edge profile](../figures/qw_potential_profile.png)

*Figure 2: Band-edge profile of the AlSbW/GaSbW/InAsW heterostructure, showing
$E_V$, $E_{\Delta SO}$, and $E_C$ as functions of position $z$. The AlSbW barriers
provide a large gap (2.384 eV) for confinement. GaSbW creates a deep valence-band
well for holes, while InAsW creates a deep conduction-band well for electrons. The
broken-gap alignment between InAsW and GaSbW is visible where the InAsW $E_C$ sits
above the GaSbW $E_V$.*

#### B.3 Numerical results

At $k_\parallel = 0$, the diagonalization of the $8 \times 101 = 808$ dimensional
Hamiltonian yields a rich spectrum. The valence band top lies at approximately
$0.020$ eV and the conduction band bottom at approximately $0.290$ eV, giving an
effective gap of about 270 meV. This gap is determined by the hybridization of InAsW
electron states and GaSbW hole states, not by any single material's band gap.

The 32 VB and 32 CB states span a wide energy range:
- The deepest VB states are bound deep in the GaSbW/AlSbW valence-band wells
- The VB top states are derived from the GaSbW valence band edge, pushed up by
  confinement
- The CB bottom states are InAsW-derived, pulled up from $E_C = -0.172$ eV by
  confinement in the 70 A InAs well
- Higher CB states extend toward the AlSbW barrier edge at 1.974 eV

The large number of bound states (32+32) reflects the deep confinement potentials:
the InAsW CB well depth exceeds 2 eV, and the GaSbW VB well depth is about 0.38 eV.

![QW subband dispersion](../figures/qw_alsbw_gasbw_inasw_bands.png)

*Figure 3: Subband dispersion $E(k_\parallel)$ for the AlSbW/GaSbW/InAsW quantum
well. Red curves are valence subbands, cyan curves are conduction subbands. The
type-III alignment brings electron and hole subbands into close proximity near the
effective gap. The strong nonparabolicity and anticrossings in the VB subbands are
due to HH-LH mixing and the coupling to InAsW conduction-band states through the
interband k.p matrix element $P$.*

#### B.4 Band alignment diagram

The band alignment of this heterostructure can be summarized schematically:

```
Energy (eV)
  2.0 ─────────────────────────────────── E_C(AlSbW)
      |                                 |
  1.5 |          AlSbW barrier          |
      |                                 |
  1.0 |                                 |
      |              E_C(GaSbW)=0.782   |
  0.5 |     ┌───────────────────────┐   |
      |     │       GaSbW           │   |
  0.0 ──────│─── E_V(GaSbW)=-0.030 │───│──── effective gap ~0.27 eV
      |     │                       │   |
 -0.2 |     │    ┌─────────────┐    │   |
      |     │    │ E_C(InAsW)  │    │   |
 -0.4 |     │    │  =-0.172    │    │   |
      |     │    │   InAsW     │    │   |
 -0.6 |     │    │ E_V(InAsW)  │    │   |
      |     │    │  =-0.590    │    │   |
 -0.8 |     │    └─────────────┘    │   |
      |     │                       │   |
      |     └───────────────────────┘   |
 -1.0 ─────────────────────────────────── E_V(AlSbW)=-0.41
```

Electrons are confined in the narrow InAsW layer (70 A), while holes are spread
across the wider GaSbW layer (270 A). The spatial separation of carriers is the
hallmark of type-III alignment.

---

## 4. Discussion

### 4.1 Type-I versus Type-III: a comparison

| Property | Type-I (GaAs/AlGaAs) | Type-III (AlSbW/GaSbW/InAsW) |
|----------|----------------------|-------------------------------|
| Electron confinement | Same layer as holes | InAsW layer |
| Hole confinement | Same layer as electrons | GaSbW layer |
| Electron-hole overlap | Large | Small |
| Effective gap | Determined by well width and offset | Determined by hybridization |
| Optical transition strength | Strong (direct) | Weak (spatially indirect) |
| Topological properties | None | Can exhibit topological phase |
| Typical applications | Lasers, LEDs, HEMTs | Topological insulators, IR detectors |

The key physical difference is the spatial separation of electrons and holes in the
type-III system. This reduces the optical matrix element for interband transitions
but enables phenomena that are impossible in type-I structures, such as the electric-
field-tunable hybridization gap that can close and reopen -- the signature of a
topological phase transition.

### 4.2 Convergence considerations

The spatial discretization introduces two convergence parameters: the grid density
($N$ or `FDstep`) and the FD accuracy order (`FDorder`). For second-order FD, the
energy error scales as $O(\Delta z^2)$, so doubling $N$ reduces the error by a
factor of 4. Higher-order schemes ($O(\Delta z^4)$ and above) converge much faster
but each k.p term block has a bandwidth proportional to `FDorder`, increasing the
Hamiltonian fill-in.

Practical guidelines:

| System | Recommended `FDstep` | `FDorder` | Energy accuracy |
|--------|---------------------|-----------|-----------------|
| Shallow type-I QW (100 A) | 100--150 | 2 | sub-meV |
| Deep type-III QW (multi-layer) | 150--250 | 2 | ~1 meV |
| Critical (hybridization gaps) | 200--400 | 4 | sub-0.1 meV |
| Production + strain | 300--500 | 4--6 | < 0.01 meV |

The boundary conditions are hard-wall (Dirichlet): the wavefunction is forced to
zero at the grid boundaries. The simulation domain must extend far enough into the
barrier that the wavefunction decays to negligible amplitude. A rule of thumb is
3--5 decay lengths beyond the quantum well on each side, which typically means
50--100 A of barrier for shallow wells and 100--200 A for deep wells.

### 4.3 Connection to other chapters

The quantum well band structure computed in this chapter provides the foundation
for several subsequent topics:

- **Chapter 03 (Wavefunctions):** The eigenvectors of the QW Hamiltonian are
  $8N$-component vectors that encode the spatial profile and band composition of
  each subband. Visualizing the envelope functions $F_n(z)$ reveals how carriers are
  distributed across the heterostructure.

- **Chapter 04 (Strain):** Lattice-mismatched heterostructures (e.g., InGaAs/GaAs)
  develop strain that shifts band edges through deformation potentials. Strain adds
  terms to the diagonal blocks of the Hamiltonian, effectively modifying the
  `profile` array and the k.p coupling strengths.

- **Chapter 07 (Self-Consistent SP):** For doped structures, the charge density from
  occupied subbands creates an electrostatic potential that feeds back into the band
  profile. The Schrodinger-Poisson loop iterates between diagonalization and Poisson
  solving until self-consistency is reached, modifying the `profile` array at each
  iteration.

- **Chapter 08 (Quantum Wire):** The QW formalism extends to 2D confinement by
  discretizing both $k_x$ and $k_y$. The block structure of the Hamiltonian is
  preserved, but the blocks become sparse CSR matrices built from Kronecker products
  of 1D FD operators.

### 4.4 Limitations

The standard QW mode assumes:

- Growth along $z$ (enforced by `confDir = 'z'`)
- Uniform grid spacing across the entire domain
- No strain (though strain parameters are in the material database for future use)
- No self-consistent charge treatment (the SC loop is a separate module; see
  Chapter 07)
- Hard-wall boundary conditions (the wavefunction vanishes at the domain edges)

For strained quantum wells, the strain Hamiltonian adds additional terms to the
diagonal blocks, shifting band edges via deformation potentials ($a_c$, $a_v$, $b$,
$d$). For doped structures, the self-consistent Schrodinger-Poisson loop iterates
between the eigenvalue solve and a Poisson solve for the electrostatic potential,
modifying the `profile` array at each iteration.

### 4.5 Variable material parameters at interfaces

A subtlety of the position-dependent k.p approach is the treatment of interfaces.
When the Luttinger parameters ($\gamma_1, \gamma_2, \gamma_3$) or the interband
matrix element $P$ change abruptly at a heterointerface, the simple product
$\gamma(z) \cdot d^2/dz^2$ is not the correct Hermitian operator. The code uses the
"forward/backward" stencil approach for order 2 and the `applyVariableCoeff` routine
for higher orders, which effectively symmetrize the kinetic energy operator. This is
the standard approach in k.p finite difference codes and is accurate for slowly-
varying envelopes, though it may introduce small errors at abrupt interfaces. The
Foreman renormalization (disabled by default via `renormalization = .False.` in
`defs.f90`) provides a more rigorous treatment at the cost of additional complexity.

### 4.6 Comparison with bulk mode

The quantum well mode shares the same 8-band basis and block topology as the bulk
Hamiltonian (Chapter 01). The key differences are:

| Aspect | Bulk (`confinement=0`) | QW (`confinement=1`) |
|--------|----------------------|---------------------|
| Matrix size | $8 \times 8$ | $8N \times 8N$ |
| $k_z$ | scalar wavevector | $-i\,d/dz$ discretized |
| Material | single | multiple layers |
| Band offsets | constant | position-dependent profile |
| External field | not applicable | electric field tilt |
| Eigenvalue solver | all 8 eigenvalues | `zheevx` partial spectrum |
| Computation time per k-point | microseconds | milliseconds to seconds |
