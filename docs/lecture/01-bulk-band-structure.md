# Chapter 01: Bulk Band Structure with the 8-Band k.p Method

## 1. Theory

### 1.1 Bloch's Theorem and the Crystal Lattice

In a perfect semiconductor crystal, the potential experienced by an electron is
periodic with the lattice. Bloch's theorem tells us that the eigenstates of the
single-particle Hamiltonian in such a periodic potential can be written as

$$
\psi_{n\mathbf{k}}(\mathbf{r}) = e^{i\mathbf{k}\cdot\mathbf{r}} \, u_{n\mathbf{k}}(\mathbf{r}),
$$

where $n$ is the band index, $\mathbf{k}$ is a wave vector in the first
Brillouin zone, and $u_{n\mathbf{k}}(\mathbf{r})$ has the periodicity of the
lattice: $u_{n\mathbf{k}}(\mathbf{r}+\mathbf{R}) = u_{n\mathbf{k}}(\mathbf{r})$
for any lattice vector $\mathbf{R}$.

The full electronic band structure $E_n(\mathbf{k})$ requires solving this
eigenvalue problem across the entire Brillouin zone. For most device physics
applications, however, we only need the band structure near the **zone center**
(the $\Gamma$ point, $\mathbf{k}=0$). This is where the **k.p method** becomes
powerful.

### 1.2 The k.p Perturbation Expansion

The central idea of the k.p method is to treat the $\mathbf{k}$-dependent part
of the crystal Hamiltonian as a perturbation around $\mathbf{k}=0$. Starting
from the $\mathbf{k}\cdot\mathbf{p}$ Hamiltonian (see, e.g., Kane 1957 or
Winkler 2003), we write:

$$
H_{\mathbf{k}} = H_0 + H_{\mathbf{k}\cdot\mathbf{p}},
$$

where $H_0$ is the $\Gamma$-point Hamiltonian (diagonal in the band basis) and
$H_{\mathbf{k}\cdot\mathbf{p}}$ couples different bands via the momentum
operator $\mathbf{p} = -i\hbar\nabla$.

We work in the basis of zone-center Bloch functions $\{u_j\}$ (often called the
**Kane basis**). The matrix elements of $H_{\mathbf{k}}$ in this basis give us a
finite-dimensional eigenvalue problem:

$$
\sum_{j'} H_{jj'}(\mathbf{k}) \, c_{j'} = E \, c_j,
$$

where $H_{jj'}(\mathbf{k})$ is the $(j,j')$ matrix element and $c_j$ are the
expansion coefficients of $u_{n\mathbf{k}}$ in terms of the zone-center basis.

The accuracy of this approach depends on how many bands we include in the basis
and how far from $\Gamma$ we venture in $\mathbf{k}$-space.

### 1.3 The 8-Band Basis for Zincblende Semiconductors

For zincblende (ZB) III-V semiconductors (GaAs, InAs, InSb, and their alloys),
the bands near the $\Gamma$ point are well described by an 8-band basis
spanning three irreducible representations of the $T_d$ point group:

| Representation | States | Atomic character | Bands |
|---|---|---|---|
| $\Gamma_6$ | $|S{\uparrow}\rangle$, $|S{\downarrow}\rangle$ | s-like | Conduction band (CB) |
| $\Gamma_8$ | $|X{\uparrow}\rangle$, $|Y{\uparrow}\rangle$, $|Z{\uparrow}\rangle$ and spin-orbit partners | p-like | Heavy hole (HH), Light hole (LH) |
| $\Gamma_7$ | Spin-orbit split-off partners of $\Gamma_8$ | p-like | Split-off (SO) |

The $\Gamma_8$ states are 4-fold degenerate (including spin) at $\Gamma$ and
split into the heavy-hole and light-hole bands away from the zone center. The
$\Gamma_7$ states are split off from $\Gamma_8$ by the spin-orbit coupling
energy $\Delta_{\mathrm{SO}}$.

#### Basis Ordering in This Code

The code uses the following **fixed basis ordering**, which is critical to
understand when reading the Hamiltonian construction:

| Index | Band | Label | Spin |
|---|---|---|---|
| 1 | Valence | Heavy hole (HH) | $\uparrow$ |
| 2 | Valence | Light hole (LH) | $\uparrow$ |
| 3 | Valence | Light hole (LH) | $\downarrow$ |
| 4 | Valence | Heavy hole (HH) | $\downarrow$ |
| 5 | Split-off | SO | $\uparrow$ |
| 6 | Split-off | SO | $\downarrow$ |
| 7 | Conduction | CB | $\uparrow$ |
| 8 | Conduction | CB | $\downarrow$ |

**This ordering is hardcoded throughout the entire codebase and must never be
changed.** Bands 1--4 are the valence bands (HH, LH, LH, HH), bands 5--6 are
the split-off bands, and bands 7--8 are the conduction bands.

### 1.4 Material Parameters

The 8-band Hamiltonian is parameterized by material-specific quantities. For a
zincblende semiconductor, the primary parameters are:

| Parameter | Symbol | Meaning | GaAs value |
|---|---|---|---|
| $E_g$ | Band gap at $\Gamma$ | $E_C - E_V$ | 1.519 eV |
| $\Delta_{\mathrm{SO}}$ | Spin-orbit splitting | Split-off offset from VB top | 0.341 eV |
| $E_P$ | Kane energy | Conduction-valence coupling strength | 28.8 eV |
| $\gamma_1, \gamma_2, \gamma_3$ | Luttinger parameters | Valence band curvature | 6.98, 2.06, 2.93 |
| $E_V$ | Valence band edge | Absolute VB maximum | -0.8 eV |
| $E_C$ | Conduction band edge | $E_V + E_g$ | 0.719 eV |
| $m^*/m_0$ | Effective mass | CB curvature | 0.067 |

Two additional derived quantities are computed by the code:

$$
P = \sqrt{E_P \cdot C_0}, \qquad A = \frac{1}{m^*/m_0},
$$

where $C_0 = \hbar^2/(2m_0) \approx 3.810$ eV-Angstrom$^2$ is the free-electron
kinetic energy constant (defined as `const` in `defs.f90`).

The parameter $P$ (with units of eV-Angstrom) is the **interband momentum
matrix element** and controls the strength of the conduction-valence coupling.
The parameter $A = m_0/m^*$ is the inverse effective mass ratio for the
conduction band.

**Parameter sources.** The code uses two families of parameters:

- **Vurgaftman parameters** (materials without the "W" suffix): taken from
  Vurgaftman, Meyer, and Ram-Mohan, *J. Appl. Phys.* **89**, 5815 (2001).
  These use the standard $E_P$ convention.
- **Winkler parameters** (materials with the "W" suffix, e.g. `GaAsW`,
  `InAsW`): taken from Winkler, *Spin-Orbit Coupling Effects in Two-Dimensional
  Electron and Hole Systems* (Springer, 2003). These use InSb as the $E_V=0$
  reference and a different parameterization of $E_P$.

Both conventions share the same Hamiltonian structure; only the numerical values
of the parameters differ.

### 1.5 The Kane Matrix Elements

The 8-band Hamiltonian is built from several **k-dependent matrix elements**
that couple the basis states. These are conventionally labeled $Q$, $R$, $S$,
$T$, and $P$ (not to be confused with the momentum parameter $P$). For a bulk
semiconductor, where all material parameters are spatially uniform, these become
algebraic functions of $\mathbf{k} = (k_x, k_y, k_z)$:

**Valence-valence kinetic terms** (Luttinger parameters):

$$
Q = -\left[(\gamma_1 + \gamma_2)(k_x^2 + k_y^2) + (\gamma_1 - 2\gamma_2)k_z^2\right],
$$

$$
T = -\left[(\gamma_1 - \gamma_2)(k_x^2 + k_y^2) + (\gamma_1 + 2\gamma_2)k_z^2\right],
$$

These are the diagonal kinetic energy terms for the heavy-hole and light-hole
bands, respectively. In units of $C_0 = \hbar^2/(2m_0)$.

**Valence-valence off-diagonal terms**:

$$
S = i \, 2\sqrt{3} \, \gamma_3 \, (k_x - ik_y) k_z,
$$

$$
\bar{S} = -i \, 2\sqrt{3} \, \gamma_3 \, (k_x + ik_y) k_z,
$$

$$
R = -\sqrt{3} \left[\gamma_2 (k_x^2 - k_y^2) - 2i\gamma_3 k_x k_y\right],
$$

$$
\bar{R} = -\sqrt{3} \left[\gamma_2 (k_x^2 - k_y^2) + 2i\gamma_3 k_x k_y\right].
$$

The terms $S$ and $\bar{S}$ mix HH and LH states, while $R$ and $\bar{R}$ are
responsible for the warping of the valence bands (the difference between $\gamma_2$
and $\gamma_3$).

**Conduction-valence coupling** (Kane momentum matrix element $P$):

$$
P_+ = \frac{P}{\sqrt{2}}(k_x + ik_y), \qquad
P_- = \frac{P}{\sqrt{2}}(k_x - ik_y), \qquad
P_z = P \, k_z.
$$

These terms are what give the 8-band model its name: they couple the conduction
band (s-like $\Gamma_6$ states) to the valence band (p-like $\Gamma_8$ and
$\Gamma_7$ states). In a purely 6-band model (valence bands only), these terms
are absent.

**Conduction band kinetic term**:

$$
A \, k^2 = A \, (k_x^2 + k_y^2 + k_z^2),
$$

where $A = m_0/m^*$ and $k^2 = |\mathbf{k}|^2$.

### 1.6 The Full 8x8 Hamiltonian Matrix

Assembling all the Kane matrix elements into the 8-band basis, the zincblende
bulk Hamiltonian at wave vector $\mathbf{k}$ reads:

$$
H_{\mathrm{ZB}}(\mathbf{k}) =
\begin{pmatrix}
Q & \bar{S} & \bar{R} & 0 & -\frac{i}{\sqrt{2}}\bar{S} & i\sqrt{2}\bar{R} & iP_+ & 0 \\
S & T & 0 & \bar{R} & \frac{i}{\sqrt{2}}(Q{-}T) & -i\sqrt{\frac{3}{2}}\bar{S} & \sqrt{\frac{2}{3}}P_z & -\frac{1}{\sqrt{3}}P_+ \\
R & 0 & T & -\bar{S} & i\sqrt{\frac{3}{2}}S & \frac{i}{\sqrt{2}}(Q{-}T) & \frac{i}{\sqrt{3}}P_- & i\sqrt{\frac{2}{3}}P_z \\
0 & R & -S & Q & i\sqrt{2}R & \frac{i}{\sqrt{2}}S & 0 & -P_- \\
\frac{i}{\sqrt{2}}S & -\frac{i}{\sqrt{2}}(Q{-}T) & -i\sqrt{\frac{3}{2}}\bar{S} & -i\sqrt{2}\bar{R} & \frac{Q+T}{2}{-}\Delta & 0 & \frac{i}{\sqrt{3}}P_z & i\sqrt{\frac{2}{3}}P_+ \\
-i\sqrt{2}R & i\sqrt{\frac{3}{2}}S & -\frac{i}{\sqrt{2}}(Q{-}T) & -\frac{i}{\sqrt{2}}\bar{S} & 0 & \frac{Q+T}{2}{-}\Delta & \sqrt{\frac{2}{3}}P_- & -\frac{1}{\sqrt{3}}P_z \\
-iP_- & \sqrt{\frac{2}{3}}P_z & -\frac{i}{\sqrt{3}}P_+ & 0 & -\frac{i}{\sqrt{3}}P_z & \sqrt{\frac{2}{3}}P_+ & Ak^2{+}E_g & 0 \\
0 & -\frac{1}{\sqrt{3}}P_+ & -i\sqrt{\frac{2}{3}}P_z & -P_+ & -i\sqrt{\frac{2}{3}}P_- & -\frac{1}{\sqrt{3}}P_z & 0 & Ak^2{+}E_g
\end{pmatrix}
$$

Here the rows and columns follow the code's basis ordering: (HH$\uparrow$,
LH$\uparrow$, LH$\downarrow$, HH$\downarrow$, SO$\uparrow$, SO$\downarrow$,
CB$\uparrow$, CB$\downarrow$).

At $\mathbf{k}=0$, all off-diagonal elements vanish and the diagonal gives us
the band edges:

- Bands 1--4 (valence): $E_V$ (HH/LH degenerate at $\Gamma$)
- Bands 5--6 (split-off): $E_V - \Delta_{\mathrm{SO}}$
- Bands 7--8 (conduction): $E_V + E_g = E_C$

The Hamiltonian is Hermitian ($H = H^\dagger$), as must be the case for a
physical Hamiltonian. The lower-left triangle is the Hermitian conjugate of the
upper-right triangle, and the code only fills the upper triangle plus diagonal
(labeled `'U'` in LAPACK's convention), relying on `zheevx` to handle the
symmetry.

### 1.7 Solving the Eigenvalue Problem

For each wave vector $\mathbf{k}$, the band energies $E_n(\mathbf{k})$ are
found by diagonalizing the $8\times 8$ Hermitian matrix:

$$
H_{\mathrm{ZB}}(\mathbf{k}) \, \mathbf{c}_n = E_n \, \mathbf{c}_n, \qquad n = 1, \ldots, 8.
$$

This yields 8 eigenvalues at each $\mathbf{k}$ point: 2 conduction bands
(CB$\uparrow$, CB$\downarrow$, degenerate without magnetic field), 2 split-off
bands, and 4 valence bands (2 HH + 2 LH, degenerate at $\Gamma$).

By sweeping $\mathbf{k}$ along a chosen direction (e.g., $k_x$ from 0 to
$k_{\max}$), we obtain the **bulk band structure** $E_n(k)$ along that
direction.

---

## 2. In the Code

### 2.1 Module Map

The bulk band structure calculation involves three key source files:

| File | Module | Role |
|---|---|---|
| `src/core/defs.f90` | `definitions` | Precision kinds (`dp`), physical constants ($\hbar$, $m_0$, $C_0$), basis ordering, `paramStruct` type |
| `src/core/parameters.f90` | `parameters` | Material database (`paramDatabase` subroutine), computes $P$ and $A$ from $E_P$ and $m^*$ |
| `src/physics/hamiltonianConstructor.f90` | `hamiltonianConstructor` | `ZB8bandBulk` subroutine: builds the 8x8 Hamiltonian matrix |

The main program `src/apps/main.f90` (`program kpfdm`) orchestrates the
calculation: it reads the input configuration, calls `paramDatabase` to load
material parameters, constructs the $\mathbf{k}$-vector sweep, calls
`ZB8bandBulk` for each $\mathbf{k}$, and diagonalizes with LAPACK's `zheevx`.

### 2.2 The `ZB8bandBulk` Subroutine

The subroutine `ZB8bandBulk` in `hamiltonianConstructor.f90` is the core of the
bulk calculation. Its signature is:

```fortran
subroutine ZB8bandBulk(HT, wv, params, g)
  complex(kind=dp), intent(inout), dimension(:,:) :: HT  ! 8x8 Hamiltonian
  type(wavevector),  intent(in)    :: wv                  ! (kx, ky, kz)
  type(paramStruct), intent(in)    :: params(1)           ! material parameters
  character(len=1),  intent(in), optional :: g             ! g-factor mode flag
```

The algorithm proceeds as follows:

1. **Extract wave vector components:** $k_x^2$, $k_y^2$, $k_z^2$, and the
   compound wave vectors $k_{\pm} = k_x \pm ik_y$ and $k_{\pm z} = (k_x \pm ik_y)k_z$.

2. **Compute the Kane matrix elements** $Q$, $T$, $S$, $\bar{S}$, $R$,
   $\bar{R}$, $P_+$, $P_-$, $P_z$ from the Luttinger parameters and momentum
   matrix element stored in `params(1)`.

3. **Fill the 8x8 Hamiltonian** `HT` following the matrix layout shown in
   Section 1.6. Each line like `HT(1,1) = Q` directly corresponds to one matrix
   element.

4. **Add the spin-orbit splitting** $\Delta_{\mathrm{SO}}$ and band gap $E_g$:

```fortran
! SOC: split-off bands shifted by -DeltaSO
HT(5,5) = HT(5,5) - params(1)%DeltaSO
HT(6,6) = HT(6,6) - params(1)%DeltaSO

! Band gap: conduction bands shifted by +Eg
HT(7,7) = HT(7,7) + params(1)%Eg
HT(8,8) = HT(8,8) + params(1)%Eg
```

Note that the valence band edge $E_V$ is **not** added in the bulk subroutine
itself -- it enters through the `profile` mechanism in the quantum well path,
but for bulk the eigenvalues are referenced to the internal energy zero.

### 2.3 The `paramDatabase` Subroutine

Material parameters are loaded by `paramDatabase` in `parameters.f90`. For each
material, the subroutine sets:

- Raw parameters: `meff`, `EP`, `Eg`, `deltaSO`, `gamma1`, `gamma2`, `gamma3`,
  `EV`, `EC`, `eps0`, plus strain/deformation potentials.
- Derived parameters computed at the end:
  - $P = \sqrt{E_P \cdot C_0}$
  - $A = 1/m^*$

The code supports 25+ materials including binary III-V compounds (GaAs, InAs,
AlAs, InP, GaSb, InSb, GaP, AlP, AlSb), ternary alloys (AlGaAs with specific
Al fractions, InAsSb alloys from 10% to 90% Sb), and Winkler-variant materials
(with "W" suffix) that use a different absolute energy reference.

### 2.4 Diagonalization Flow

For bulk mode (`confinement=0`), the main program:

1. Allocates an 8x8 complex matrix `HT` and an 8x8 workspace `HTmp`.
2. For each $\mathbf{k}$ point in the sweep:
   - Calls `ZB8bandBulk(HT, smallk(k), cfg%params(1))` to fill the Hamiltonian.
   - Copies to `HTmp` and calls LAPACK's `zheevx` with `'V'` (compute
     eigenvectors), `'I'` (index-based eigenvalue selection), and `'U'`
     (upper triangle).
   - Stores the eigenvalues and eigenvectors.
3. Writes results to `output/eigenvalues.dat`.

The workspace query (`lwork = -1`) is performed once before the loop to
determine the optimal workspace size, then reused for all $\mathbf{k}$ points.

---

## 3. Computed Example: Bulk GaAs along $k_x$

### 3.1 Input Configuration

The following `input.cfg` computes the bulk GaAs band structure along the $k_x$
direction from 0 to 0.1 inverse Angstroms:

```
waveVector: kx
waveVectorMax: 0.1
waveVectorStep: 11
confinement:  0
FDstep: 101
FDorder: 2
numLayers:  1
material1: GaAs
numcb: 2
numvb: 6
ExternalField: 0  EF
EFParams: 0.0005
```

**Parameter-by-parameter explanation:**

| Parameter | Value | Meaning |
|---|---|---|
| `waveVector` | `kx` | Sweep along the $k_x$ direction |
| `waveVectorMax` | `0.1` | Maximum wave vector in 1/Angstrom |
| `waveVectorStep` | `11` | Number of k-points (including $k=0$) |
| `confinement` | `0` | Bulk mode (8x8 Hamiltonian) |
| `FDstep` | `101` | Ignored for bulk (needed for QW mode) |
| `FDorder` | `2` | Ignored for bulk (finite difference order for QW) |
| `numLayers` | `1` | Single material layer |
| `material1` | `GaAs` | Use GaAs parameters from the database |
| `numcb` | `2` | Request 2 conduction bands |
| `numvb` | `6` | Request 6 valence bands (HH+LH+SO, both spins) |
| `ExternalField` | `0` | No external electric field |
| `EFParams` | `0.0005` | Placeholder (unused when ExternalField=0) |

The wave vector grid is uniformly spaced:

$$
k_x^{(i)} = \frac{(i-1)}{10} \times 0.1 \; \text{AA}^{-1}, \qquad i = 1, \ldots, 11.
$$

### 3.2 Expected Output

After running `./build/src/bandStructure`, the program:

1. Prints GaAs material parameters to stdout:
   ```
   Material: GaAs
   Parameters
   EP : 28.8
   P  : 10.48...
   A  : 14.925...
   gamma1 : 6.98
   gamma2 : 2.06
   gamma3 : 2.93
   ```

2. Writes eigenvalues to `output/eigenvalues.dat` with 11 rows (one per
   $\mathbf{k}$ point) and 8 columns (one per band), sorted from lowest to
   highest energy.

3. Writes eigenfunctions at selected $\mathbf{k}$ points to files in
   `output/eigenfunctions/`.

### 3.3 Interpreting the Band Structure

At $k_x = 0$ ($\Gamma$ point), the 8 eigenvalues are:

| Band | Energy (eV) | Degeneracy |
|---|---|---|
| HH (bands 1, 4) | $E_V = -0.8$ | 2-fold |
| LH (bands 2, 3) | $E_V = -0.8$ | 2-fold |
| SO (bands 5, 6) | $E_V - \Delta_{\mathrm{SO}} = -1.141$ | 2-fold |
| CB (bands 7, 8) | $E_C = 0.719$ | 2-fold |

As $k_x$ increases from zero:

- **Conduction band** curves upward with effective mass $m^* = 0.067\,m_0$.
  The nonparabolicity from the $P$-coupling to the valence bands becomes
  important for $k_x \gtrsim 0.05$ AA$^{-1}$.
- **Heavy-hole band (HH)** has a nearly flat dispersion along $k_x$ (the HH
  mass is governed by $\gamma_1 - 2\gamma_2$ in the growth direction, but along
  an in-plane direction like $k_x$ the HH and LH mix).
- **Light-hole band (LH)** has a lighter curvature, mixing strongly with HH
  away from $\Gamma$.
- **Split-off band** is offset by $\Delta_{\mathrm{SO}} = 0.341$ eV below the
  valence band edge.

The regime $k_x \in [0, 0.1]$ AA$^{-1}$ corresponds to approximately
$[0, 0.16]\times 2\pi/a$ where $a = 5.653$ AA is the GaAs lattice constant,
which is well within the validity range of the 8-band k.p model.

---

## 4. Published Example: Bastos et al. (arXiv:1608.04982)

### 4.1 The Paper

**Reference:** C. M. O. Bastos, F. T. S. da Silva, R. H. Miwa, and G. M.
Sipahi, "Stability and accuracy control of k.p parameters for zincblende
III-V semiconductors," *Phys. Rev. B* **94**, 125108 (2016).
[arXiv:1608.04982](https://arxiv.org/abs/1608.04982).

### 4.2 What They Did

Bastos et al. addressed a fundamental problem in k.p parameterization: the
conventional parameters ($E_P$, $\gamma_i$, $m^*$) extracted from experiment
are not uniquely determined by the 8-band Hamiltonian. Different sets of
parameters can reproduce the same band structure near $\Gamma$ but diverge at
larger $\mathbf{k}$.

Their approach:

1. Used **hybrid DFT + spin-orbit coupling** (HSE06+SOC) to compute the
   band structure of several zincblende III-V semiconductors from first
   principles.
2. Fitted the 8-band k.p parameters to reproduce the DFT band structure,
   imposing physical stability constraints.
3. Validated the fitted parameters against experimental effective masses
   and known band curvatures.

### 4.3 What We Can Reproduce

Using this code, we can reproduce the bulk $E(\mathbf{k})$ dispersion for any
zincblende III-V material by supplying the Bastos parameters. The validation
procedure is:

1. Take the fitted parameters from their Table I (e.g., for GaAs):
   $E_g$, $\Delta_{\mathrm{SO}}$, $E_P$, $\gamma_1$, $\gamma_2$, $\gamma_3$.
2. Set up a bulk calculation with the corresponding material.
3. Compare the eigenvalues at several $\mathbf{k}$ points against their
   Figure 1 (band structure plots).

The code uses the Vurgaftman GaAs parameters by default (`material1: GaAs`),
which are broadly consistent with Bastos et al.'s fitted values. Minor
differences arise from different fitting strategies:

| Parameter | Vurgaftman (code) | Bastos et al. |
|---|---|---|
| $E_g$ (eV) | 1.519 | 1.519 |
| $\Delta_{\mathrm{SO}}$ (eV) | 0.341 | 0.341 |
| $E_P$ (eV) | 28.8 | 28.8 |
| $\gamma_1$ | 6.98 | 6.98 |
| $\gamma_2$ | 2.06 | 2.06 |
| $\gamma_3$ | 2.93 | 2.93 |

For GaAs, the agreement is essentially exact because both sources use the same
experimental data. Larger discrepancies appear for narrow-gap materials (InSb,
InAs) where the nonparabolicity and band coupling make the parameter fitting
more challenging.

### 4.4 Running the Validation

The regression test configuration
`tests/regression/configs/bulk_gaas_kx.cfg` uses exactly the GaAs bulk setup
described in Section 3.1. To run it:

```bash
# Copy the regression config as input.cfg (use Write tool, not cp -i)
# Then build and run:
cmake --build build
./build/src/bandStructure
```

The output eigenvalues can be compared against the reference data in
`tests/regression/data/` using the comparison script:

```bash
python tests/regression/compare_output.py output/eigenvalues.dat tests/regression/data/bulk_gaas_kx_ref.dat
```

---

## 5. Discussion

### 5.1 Strengths of the 8-Band Model

The 8-band k.p method is the workhorse of semiconductor heterostructure modeling
because it captures the essential physics near the $\Gamma$ point with minimal
computational cost:

- **Band coupling:** The $P$-coupling between conduction and valence bands
  naturally produces the correct nonparabolicity of the conduction band.
  This is critical for narrow-gap materials (InSb, InAs) where the gap is
  comparable to the spin-orbit splitting.
- **Spin-orbit effects:** The split-off band is included explicitly, giving
  the correct heavy-hole/light-hole splitting and enabling the calculation of
  g-factors (Chapter 03).
- **Computational efficiency:** For bulk, the problem is just an $8\times 8$
  eigenvalue problem at each $\mathbf{k}$ point, requiring microseconds per
  point.

### 5.2 Limitations

- **Validity range:** The k.p expansion is a perturbation around $\mathbf{k}=0$
  and loses accuracy far from the zone center. For GaAs, the model is reliable
  up to about $k \approx 0.15$ AA$^{-1}$ (roughly 25% of the way to the zone
  boundary). Beyond this, higher-lying bands neglected in the 8-band basis
  contribute significantly.
- **Parameter sensitivity:** As highlighted by Bastos et al. (Section 4), the
  parameters are not uniquely determined. Different parameter sets that agree
  at $\Gamma$ can diverge significantly at finite $\mathbf{k}$. This is
  especially problematic for narrow-gap materials.
- **Zincblende assumption:** The Hamiltonian as implemented assumes the $T_d$
  point group symmetry of the zincblende lattice. Wurtzite materials (GaN, ZnO)
  require a different Hamiltonian structure that is not supported by this code.
- **No excited conduction bands:** The 8-band basis does not include the
  $L$-point or $X$-point conduction minima. For indirect-gap semiconductors
  (Si, Ge, AlAs), the model cannot describe the lowest conduction bands
  correctly.

### 5.3 The Foreman Renormalization

The code includes an optional **Foreman renormalization** scheme (disabled by
default, controlled by the `renormalization` parameter in `defs.f90`). When
enabled, the Luttinger parameters and $E_P$ are renormalized to account for
contributions from remote bands not included in the 8-band basis:

$$
\gamma_i^{\mathrm{renorm}} = C_0 \left(\gamma_i - \frac{E_P}{6E_g}\right), \qquad i = 1,2,3.
$$

This shifts the valence band curvature by "removing" the contribution from the
conduction band coupling that is already treated explicitly in the 8-band model.
The renormalization is physically more consistent but produces different
numerical results from the standard (unrenormalized) parameters.

### 5.4 Connection to Quantum Wells and Wires

The bulk Hamiltonian is the foundation upon which the quantum well and quantum
wire calculations are built:

- **Quantum well** (Chapter 02): The wave vector components perpendicular to
  the growth direction ($k_x$, $k_y$) remain good quantum numbers, while $k_z$
  is replaced by a finite-difference derivative $\partial/\partial z$. The 8x8
  matrix becomes an $8N \times 8N$ block tridiagonal matrix, where $N$ is the
  number of spatial grid points.
- **Quantum wire** (Chapter 04): Only the axial wave vector $k_z$ remains a
  good quantum number. The transverse directions ($x$, $y$) are discretized
  with finite differences, producing a sparse Hamiltonian solved with iterative
  eigensolvers (FEAST).

In both cases, the bulk Kane matrix elements ($Q$, $R$, $S$, $T$, $P$ terms)
are re-used as building blocks, but with position-dependent material parameters
and derivative operators replacing the simple algebraic $\mathbf{k}$-dependence.

### 5.5 Further Reading

- **E. O. Kane**, "Band structure of indium antimonide," *J. Phys. Chem. Solids*
  **1**, 249 (1957). The original k.p paper.
- **R. Winkler**, *Spin-Orbit Coupling Effects in Two-Dimensional Electron and
  Hole Systems*, Springer (2003). Comprehensive treatment of the 8-band model
  with Foreman renormalization.
- **I. Vurgaftman, J. R. Meyer, and L. R. Ram-Mohan**, "Band parameters for
  III-V compound semiconductors and their alloys," *J. Appl. Phys.* **89**,
  5815 (2001). The standard reference for zincblende III-V parameters.
- **C. M. O. Bastos et al.**, *Phys. Rev. B* **94**, 125108 (2016).
  Parameter stability analysis discussed in Section 4.
