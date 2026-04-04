# Chapter 02: Quantum Well Band Structure

## 1. Theory

### 1.1 From bulk to confined: the envelope function approximation

In Chapter 01 we studied the 8-band k.p Hamiltonian for bulk semiconductors, where the crystal has full translational symmetry in all three spatial directions. The eigenstates are Bloch waves characterized by a continuous wavevector $\mathbf{k} = (k_x, k_y, k_z)$, and the Hamiltonian is an $8 \times 8$ complex Hermitian matrix for each $\mathbf{k}$-point.

A semiconductor heterostructure breaks this translational symmetry. When a thin layer of one material (e.g., InAs) is sandwiched between layers of another (e.g., AlSb), the potential becomes a function of the growth direction $z$. The key insight of the envelope function approximation is that we can still expand the wavefunction in the periodic part of the Bloch functions $|u_n\rangle$ of the constituent materials, but now the expansion coefficients become $z$-dependent:

$$\Psi(\mathbf{r}) = \sum_{n=1}^{8} F_n(z) \, e^{i(k_x x + k_y y)} \, |u_n\rangle$$

where $F_n(z)$ are the envelope functions. The in-plane wavevector $\mathbf{k}_\parallel = (k_x, k_y)$ remains a good quantum number because translational symmetry is preserved in the $x$-$y$ plane. Along $z$, however, $k_z$ is no longer conserved.

### 1.2 Replacing $k_z$ with the derivative operator

The central step in going from the bulk 8-band Hamiltonian to the quantum well Hamiltonian is the substitution:

$$k_z \longrightarrow -i \frac{d}{dz}$$

This transforms every occurrence of $k_z$ and $k_z^2$ in the bulk Hamiltonian into differential operators acting on the envelope functions $F_n(z)$. The k.p terms that were simple scalars in the bulk case now become operators. Specifically:

- Terms proportional to $k_z^2$ become second-derivative operators $\frac{d^2}{dz^2}$
- Terms proportional to $k_z$ become first-derivative operators $\frac{d}{dz}$
- Terms proportional to $k_x$ and $k_y$ remain as scalar multipliers (the in-plane directions are free)

For the zinc-blende 8-band basis, the principal k.p terms and their $z$-dependence are:

| Symbol | Bulk form | QW form |
|--------|-----------|---------|
| $Q$ | $-(\gamma_1 + \gamma_2)(k_x^2 + k_y^2) - (\gamma_1 - 2\gamma_2)k_z^2$ | $-(\gamma_1 + \gamma_2)k_\parallel^2 - (\gamma_1 - 2\gamma_2)\frac{d^2}{dz^2}$ |
| $T$ | $-(\gamma_1 - \gamma_2)(k_x^2 + k_y^2) - (\gamma_1 + 2\gamma_2)k_z^2$ | $-(\gamma_1 - \gamma_2)k_\parallel^2 - (\gamma_1 + 2\gamma_2)\frac{d^2}{dz^2}$ |
| $S$ | $i2\sqrt{3}\,\gamma_3\, k_- k_z$ | $2\sqrt{3}\,\gamma_3\, k_- \frac{d}{dz}$ |
| $R$ | $-\sqrt{3}\bigl(\gamma_2(k_x^2 - k_y^2) - 2i\gamma_3 k_x k_y\bigr)$ | unchanged (no $k_z$ dependence) |
| $A$ | $A \, k^2$ | $A\bigl(k_\parallel^2 + \frac{d^2}{dz^2}\bigr)$ |
| $P_z$ | $P \, k_z$ | $-i P \frac{d}{dz}$ |

where $k_\pm = k_x \pm i k_y$ and the Luttinger parameters $\gamma_1, \gamma_2, \gamma_3$, the interband momentum matrix element $P$, and the remote-band parameter $A$ all become position-dependent functions $\gamma_1(z), \gamma_2(z), \ldots$ that take the values of whichever material is present at position $z$.

### 1.3 Block matrix structure: $8 \times 8$ blocks, each $N \times N$

After discretizing the $z$-direction on a grid of $N$ points (the `FDstep` parameter), the quantum well Hamiltonian becomes a large $8N \times 8N$ complex Hermitian matrix. Its natural structure is an $8 \times 8$ block matrix, where each block $(\alpha, \beta)$ is itself an $N \times N$ matrix:

$$H_{\text{QW}} = \begin{pmatrix}
H_{11} & H_{12} & \cdots & H_{18} \\
H_{21} & H_{22} & \cdots & H_{28} \\
\vdots & & \ddots & \vdots \\
H_{81} & H_{82} & \cdots & H_{88}
\end{pmatrix}$$

Here, the band index pairs $(\alpha, \beta) \in \{1,\ldots,8\}$ follow the standard basis ordering:
- Bands 1--4: heavy-hole (HH), light-hole (LH), LH', split-off (SO) valence bands
- Bands 5--6: split-off doublet
- Bands 7--8: conduction band (CB) doublet

Each off-diagonal block $H_{\alpha\beta}$ contains the k.p coupling terms (Q, R, S, T, P operators) connecting bands $\alpha$ and $\beta$, discretized on the $N$-point grid. Each diagonal block $H_{\alpha\alpha}$ contains the self-energy of band $\alpha$, which includes:

1. The kinetic energy (k.p terms involving $\gamma_i$ or $A$)
2. The **band offset** $V_\alpha(z)$ from the `profile` array

### 1.4 The kpterms array: precomputed position-dependent operators

The code precomputes a three-dimensional array `kpterms(N, N, 10)` during the initialization phase (`confinementInitialization`). The ten terms encode:

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

Terms 1--4 and 10 are purely diagonal: they contain the material parameters evaluated at each grid point. Terms 5--9 couple neighboring grid points through the finite difference stencil. The key design principle is that the position dependence of the material parameters is baked into these operators once, so that the Hamiltonian assembly at each $\mathbf{k}_\parallel$ point is a fast linear combination of precomputed matrices.

### 1.5 Band offsets and the profile array

The heterostructure potential is stored in the `profile(N, 3)` array:

- `profile(:, 1)` = $E_V(z)$: valence band edge (applied to bands 1--4)
- `profile(:, 2)` = $E_V(z) - \Delta_{\text{SO}}(z)$: split-off band edge (applied to bands 5--6)
- `profile(:, 3)` = $E_C(z)$: conduction band edge (applied to bands 7--8)

These profiles are built by assigning the band-edge parameters of each material to the grid points within that material's spatial range. The parameters $E_V$ and $E_C$ are stored in the material database (`parameters.f90`) and define the relative alignment of bands across the heterostructure.

For an electric field along $z$ with magnitude $\mathcal{E}$, the code adds a linear potential:

$$V_{\text{elec}}(z) = -\mathcal{E} \cdot L \cdot \frac{z + z_{\min}}{2\,z_{\min}}$$

where $L$ is the total simulation size and $z_{\min}$ is the leftmost grid coordinate. This tilts the band edges, breaking the structural inversion symmetry.

### 1.6 Type-I versus type-II band alignment

The relative position of conduction and valence band edges across a heterointerface determines the band alignment type:

**Type-I (straddling gap):** The band gap of material A sits entirely within the gap of material B. Both electrons and holes are confined in the same spatial region. Example: GaAs quantum well in AlGaAs barriers.

**Type-II (staggered or broken gap):** The conduction band of one material lies below the valence band of the other. Electrons and holes are spatially separated. The canonical example is the **InAs/GaSb** system:

$$E_C(\text{InAs}) < E_V(\text{GaSb})$$

This is called a **broken-gap** alignment because the InAs conduction band edge sits approximately 150 meV below the GaSb valence band edge. In a GaSb/InAs/GaSb quantum well, electrons accumulate in the InAs layer while holes reside in the GaSb layers, creating a natural spatial separation. Under an applied electric field, the hybridization of these electron and hole states can be tuned continuously, giving rise to a topological phase transition. This is the physical mechanism exploited in the published examples discussed in Section 4.

### 1.7 Finite difference discretization

The $z$-direction is discretized on a uniform grid of $N$ points with spacing $\Delta z$. The code supports FD accuracy orders 2, 4, 6, 8, and 10 (controlled by the `FDorder` parameter, default 2).

For the standard second-order scheme (`FDorder = 2`), the second derivative is discretized as:

$$\frac{d^2 f}{dz^2}\bigg|_{z_i} \approx \frac{f_{i-1} - 2f_i + f_{i+1}}{(\Delta z)^2}$$

and the first derivative as:

$$\frac{df}{dz}\bigg|_{z_i} \approx \frac{f_{i+1} - f_{i-1}}{2\,\Delta z}$$

This yields tridiagonal matrices for both the first and second derivative operators. The variable-coefficient product $g(z) \cdot d^2f/dz^2$ is handled by computing the matrix-vector product of the stencil matrices with the material parameter profile, then inserting the result into the appropriate tridiagonal band of `kpterms`.

For higher orders (`FDorder >= 4`), the code uses banded FD matrices with wider stencils. The interior points use central stencils, while the boundary points use one-sided forward/backward stencils to maintain the requested accuracy order at the edges of the simulation domain. This avoids periodic boundary assumptions and properly handles hard-wall (zero-wavefunction) boundary conditions.

The grid must satisfy $N \geq \text{FDorder} + 1$ for the stencil to fit. In practice, $N$ (the `FDstep` parameter) is typically 50--500, providing ample resolution.

---

## 2. In the Code

### 2.1 Initialization: `confinementInitialization`

When `confinement = 1` (QW mode), the input parser triggers `confinementInitialization` (in `hamiltonianConstructor.f90`). This routine:

1. **Builds the z-grid:** From `startPos` and `endPos` of the first layer, the total size $L$ is computed. The grid spacing is $\Delta z = L / (N - 1)$, and the coordinate array is `z(i) = startPos(1) + (i-1) * delta`.

2. **Fills the profile array:** For each layer $i$, grid points in the range `[intStartPos(i) : intEndPos(i)]` receive:
   - `profile(:, 1) = params(i)%EV`
   - `profile(:, 2) = params(i)%EV - params(i)%DeltaSO`
   - `profile(:, 3) = params(i)%EC`

3. **Fills the kpterms array:** Material parameters are extracted per grid point and combined with FD stencil matrices. For order 2, this uses a forward/central/backward stencil decomposition to build tridiagonal operators. For higher orders, the code calls `buildFD2ndDerivMatrix` and `buildFD1stDerivMatrix` from `finitedifferences.f90`, then applies `applyVariableCoeff` to multiply each FD matrix by the position-dependent parameter profile.

4. **Applies electric field:** If `ExternalField = 1` with type `EF`, the routine `externalFieldSetup_electricField` adds a linear tilt to the profile: `profile(i,:) -= (Evalue * totalSize) * (z(i) + z(1)) / (2 * z(1))`.

### 2.2 Hamiltonian assembly: `ZB8bandQW`

For each $\mathbf{k}_\parallel$ point in the wavevector sweep, the routine `ZB8bandQW` (in `hamiltonianConstructor.f90`) assembles the full $8N \times 8N$ Hamiltonian:

1. **Compute k.p blocks:** The precomputed `kpterms` are combined with the in-plane wavevector components $(k_x, k_y)$ to build the ten $N \times N$ blocks: `Q`, `T`, `S`, `SC`, `R`, `RC`, `PP`, `PM`, `PZ`, `A`.

   For example, the Q block is:
   ```
   Q(ii,jj) = -((kpterms(ii,jj,1) + kpterms(ii,jj,2)) * k_par^2 + kpterms(ii,jj,7))
   ```
   where `kpterms(:,:,1) = gamma1`, `kpterms(:,:,2) = gamma2`, and `kpterms(:,:,7)` already contains the discretized $-(\gamma_1 - 2\gamma_2) \cdot d^2/dz^2$ operator.

2. **Populate the 8x8 block matrix:** The code fills the Hamiltonian using Fortran array sections:
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

4. **Diagonalize:** The band structure executable calls LAPACK's `zheevx` to find the requested number of eigenvalues (`numcb` conduction + `numvb` valence) at each k-point.

### 2.3 Input parsing for QW mode

The input parser (`input_parser.f90`) handles QW setup when `confinement = 1`:

- Reads `numLayers` material specifications, each with a name, start position, and end position (in Angstroms)
- Computes integer grid indices `intStartPos(i)` and `intEndPos(i)` that map each material to a contiguous range of FD grid points
- Stores the grid in `cfg%z(:)` and copies it to `cfg%grid%z(:)` via `init_grid_from_config`
- Calls `paramDatabase` to fill `cfg%params(:)` with material parameters from the database

### 2.4 Material database: band offsets

The material parameters in `parameters.f90` include $E_V$ and $E_C$ for each material. These define the absolute band edge positions used for the heterostructure profile. The convention is $E_C = E_V + E_g$ within each material, but the *absolute* values of $E_V$ (and hence $E_C$) differ between materials, defining the band offsets at the interface.

For the example config materials:

| Material | $E_V$ (eV) | $E_C$ (eV) | $E_g$ (eV) | $\Delta_{\text{SO}}$ (eV) |
|----------|-----------|-----------|------------|--------------------------|
| AlSb     | $-0.41$   | $+1.976$  | 2.386      | 0.676                    |
| GaSb     | $-0.80$   | $+0.012$  | 0.812      | 0.76                     |
| InAs     | $-0.59$   | $-0.173$  | 0.417      | 0.39                     |

Note that $E_C(\text{InAs}) = -0.173$ eV lies well above $E_V(\text{GaSb}) = -0.80$ eV -- the InAs conduction band is above the GaSb valence band by about 0.627 eV. This means the AlSb/GaSb/InAs system is a type-II staggered gap, not a broken gap. For the broken-gap configuration (InAs/GaSb with $E_C(\text{InAs}) < E_V(\text{GaSb})$), the band offsets would need different alignment values (e.g., using Winkler parameter sets `InAsW`, `GaSbW`).

---

## 3. Computed Example

### 3.1 Input configuration: AlSb/GaSb/InAs type-II quantum well

The following example is taken from the regression test config `tests/regression/configs/qw_alsb_gasb_inas.cfg`:

```
waveVector: kx
waveVectorMax: 0.1
waveVectorStep: 11
confinement:  1
FDstep: 101
FDorder: 2
numLayers:  3
material1: AlSb -250  250 0
material2: GaSb -135  135 0.2414
material3: InAs  -35   35 -0.0914
numcb: 32
numvb: 32
ExternalField: 0  EF
EFParams: 0.0005
```

### 3.2 Structure walkthrough

This input defines a three-layer heterostructure spanning $z \in [-250, 250]$ Angstroms:

- **Layer 1 (AlSb):** Full range $[-250, 250]$ A. This is the barrier material that defines the simulation domain. Its band edges are $E_V = -0.41$ eV and $E_C = +1.976$ eV, providing a large band gap (2.386 eV) for confinement.

- **Layer 2 (GaSb):** Range $[-135, 135]$ A, width 270 A. The GaSb layer has $E_V = -0.80$ eV, significantly lower than the AlSb valence band. This creates a deep potential well for holes.

- **Layer 3 (InAs):** Range $[-35, 35]$ A, width 70 A. The InAs layer has $E_C = -0.173$ eV, well below the AlSb conduction band, creating a confining well for electrons at the center.

The grid uses `FDstep = 101` points, giving $\Delta z = 500/100 = 5$ Angstroms. The wavevector sweep goes from $k_x = 0$ to $k_x = 0.1$ A$^{-1}$ in 11 steps. A small electric field of 0.5 kV/cm is applied (`ExternalField = 0` disables it in this config, but the `EFParams` line is present for compatibility).

### 3.3 Running the example

```bash
# Copy the config to input.cfg (use Write tool, not cp, due to cp -i alias)
# Then run:
./build/src/bandStructure
```

The program reads `input.cfg`, initializes the QW grid, builds `profile` and `kpterms`, then sweeps over 11 k-points from $k_x = 0$ to $0.1$ A$^{-1}$. At each k-point, it diagonalizes the $8 \times 101 = 808$ dimensional complex Hermitian Hamiltonian using `zheevx`, requesting the 32 lowest and 32 highest eigenvalues. Output files are written to the `output/` directory.

### 3.4 Expected physics

The band structure shows:
- Several subbands in the conduction band well (primarily InAs-derived)
- Heavy-hole and light-hole subbands in the GaSb valence band well
- The type-II alignment causes the InAs electron subbands and GaSb hole subbands to lie close in energy, with possible anticrossings at finite $k_\parallel$
- If an electric field were applied, the relative position of electron and hole states would shift, modifying the anticrossing gaps

---

## 4. Published Examples

### P1: Electrical tuning of helical edge states in topological multilayers

**Reference:** Campos et al., arXiv:1903.02687

This work studies GaSb/InAs/GaSb multilayer structures using the full 3D 8-band k.p method with electric field tuning. The key physics is the interplay between the type-II broken-gap alignment of InAs and GaSb and the structural inversion asymmetry induced by an external electric field.

**Physics:** In a symmetric GaSb/InAs/GaSb double quantum well, the InAs conduction band states hybridize with the GaSb valence band states. The resulting band inversion creates a topologically non-trivial phase. When an electric field is applied perpendicular to the layers, it breaks the inversion symmetry and tunes the hybridization gap. At a critical field, the gap closes and reopens, signaling a topological phase transition between a quantum spin Hall insulator and a trivial insulator.

**Reproducibility with this code:** The full band structure calculations in this paper can be reproduced using the QW mode (`confinement = 1`) with an appropriate `GaSb/InAs/GaSb` layer structure. The electric field is enabled via `ExternalField: 1 EF` and controlled with the `EFParams` value. The wavevector sweep along $k_x$ (or $k_y$) maps out the dispersion $E(k_\parallel)$ that reveals the hybridization gap and its evolution with field.

Key input parameters to reproduce:
- Three-layer structure: GaSb / InAs / GaSb (the Winkler parameter sets `GaSbW` and `InAsW` may be needed for the broken-gap alignment)
- `ExternalField: 1 EF` with varying `EFParams` to scan the electric field dependence
- Sufficient `FDstep` (150--300) for convergence of the subband energies
- `FDorder: 4` or higher for improved accuracy of the confined state energies

### P5: Superconducting proximity effect in InAsSb surface QWs

**Reference:** Mayer et al., arXiv:1909.12571

This paper investigates InAs$_{1-x}$Sb$_x$ surface quantum wells on AlAs$_{1-y}$Sb$_y$ barriers. The InAsSb alloy forms a narrow-gap quantum well near the surface, enabling strong coupling to a superconducting contact deposited on top.

**Physics:** The InAsSb quantum well hosts a high-mobility two-dimensional electron gas with strong spin-orbit coupling due to the structural inversion asymmetry of the surface well. The band structure calculated by 8-band k.p reveals the subband spacing, effective masses, and spin-orbit splitting that determine the proximity effect efficiency.

**Reproducibility with this code:** The QW band structure (subband dispersions, effective masses, spin-orbit couplings) can be computed using the QW mode with an appropriate InAsSb/AlAsSb layer structure. The InAsSb alloy parameters would need to be added to `parameters.f90` using Vegard's law interpolation between InAs and InSb (both already in the database). The superconducting proximity effect itself (induced gap, Majorana physics) is beyond the scope of the k.p calculation but the normal-state band structure is the essential input.

---

## 5. Discussion

### 5.1 Convergence considerations

The spatial discretization introduces two convergence parameters: the grid density ($N$ or `FDstep`) and the FD accuracy order (`FDorder`). For second-order FD, the energy error scales as $O(\Delta z^2)$, so doubling $N$ reduces the error by a factor of 4. Higher-order schemes ($O(\Delta z^4)$ and above) converge much faster, but each k.p term block has a bandwidth proportional to `FDorder`, increasing the Hamiltonian fill-in.

In practice, for quantum wells with typical widths of 50--200 A and barrier regions of 50--100 A on each side, `FDstep = 100--200` with `FDorder = 2` gives sub-meV accuracy for the lowest subbands. For critical applications (e.g., resolving small hybridization gaps in type-II systems), `FDorder = 4` with `FDstep = 150--300` is recommended.

### 5.2 Boundary conditions

The FD discretization implicitly imposes hard-wall (Dirichlet) boundary conditions: the wavefunction is zero at the grid boundaries. This is physically appropriate when the barriers are thick enough that the wavefunction decays to negligible amplitude before reaching the boundary. If the barrier is too thin, spurious reflections can occur. The rule of thumb is that the simulation domain should extend at least 3--5 decay lengths beyond the quantum well on each side. For typical semiconductor parameters, 50--100 A of barrier on each side is sufficient.

### 5.3 Variable material parameters at interfaces

A subtlety of the position-dependent k.p approach is the treatment of interfaces. When the Luttinger parameters ($\gamma_1, \gamma_2, \gamma_3$) or the interband matrix element $P$ change abruptly at a heterointerface, the simple product $\gamma(z) \cdot d^2/dz^2$ is not the correct Hermitian operator. The code uses the "forward/backward" stencil approach for order 2 and the `applyVariableCoeff` routine for higher orders, which effectively symmetrize the kinetic energy operator. This is the standard approach in k.p finite difference codes and is accurate for slowly-varying envelopes, though it may introduce small errors at abrupt interfaces. The Foreman renormalization (disabled by default via `renormalization = .False.` in `defs.f90`) provides a more rigorous treatment at the cost of additional complexity.

### 5.4 Comparison with bulk mode

The quantum well mode shares the same 8-band basis and block topology as the bulk Hamiltonian (Chapter 01). The key differences are:

| Aspect | Bulk (`confinement=0`) | QW (`confinement=1`) |
|--------|----------------------|---------------------|
| Matrix size | $8 \times 8$ | $8N \times 8N$ |
| $k_z$ | scalar wavevector | $-i\,d/dz$ discretized |
| Material | single | multiple layers |
| Band offsets | constant | position-dependent profile |
| External field | not applicable | electric field tilt |
| Eigenvalue solver | all 8 eigenvalues | `zheevx` partial spectrum |

### 5.5 Connection to wire mode

The quantum well mode is the 1D confinement ancestor of the 2D wire mode (`confinement = 2`). In wire mode, both $k_x$ and $k_y$ are discretized, and the k.p terms become sparse CSR matrices built from Kronecker products of 1D FD operators. The physical principles and the block topology of the $8 \times 8$ Hamiltonian are identical. The key difference is computational: the QW Hamiltonian is dense within each block (bandwidth $O(\text{FDorder})$ per k.p term) and solved with dense LAPACK, while the wire Hamiltonian is sparse and solved with MKL SpBLAS. See Chapter 03 for the wire extension.

### 5.6 Limitations and extensions

The standard QW mode assumes:
- Growth along $z$ (the code enforces `confDir = 'z'`)
- Uniform grid spacing within each layer
- No strain (though strain parameters are in the database for future use)
- No self-consistent treatment of charge (the SC loop is a separate module)

For strained quantum wells, the strain Hamiltonian adds additional terms to the diagonal blocks, shifting band edges via deformation potentials ($a_c$, $a_v$, $b$, $d$). For doped structures, the self-consistent Schrodinger-Poisson loop iterates between the eigenvalue solve and a Poisson solve for the electrostatic potential, modifying the `profile` array at each iteration.
