# Chapter 05: Landau g-Factor from Lowdin Partitioning

## 5.1 Theory

### 5.1.1 The Zeeman effect in semiconductors

When a magnetic field $\mathbf{B}$ is applied to a semiconductor, the electron states split according to their spin and orbital angular momentum. For a free electron, the Zeeman Hamiltonian is

$$
H_Z = \frac{\mu_B}{2} \, \mathbf{B} \cdot \left( g_0 \, \boldsymbol{\sigma} + 2\mathbf{L} \right),
$$

where $\mu_B = e\hbar/(2m_0)$ is the Bohr magneton, $g_0 \approx 2.00231$ is the free-electron g-factor (the code uses the CODATA value `g_free` from `defs.f90`), $\boldsymbol{\sigma} = (\sigma_x, \sigma_y, \sigma_z)$ are the Pauli matrices, and $\mathbf{L}$ is the orbital angular momentum operator. The first term describes the coupling of $\mathbf{B}$ to the electron spin; the second describes the coupling to orbital angular momentum.

In a crystal, the orbital angular momentum is no longer a good quantum number. The spin-orbit coupling inherited from the atomic $p$-like valence band states mixes spin and orbital character, producing an **effective g-factor** that deviates significantly from $g_0$. For conduction-band electrons in narrow-gap semiconductors like InSb or InAs, the g-factor can be $|g^*| \gg g_0$ and even change sign. The deviation $\Delta g = g^* - g_0$ is almost entirely due to inter-band coupling, which we compute via second-order perturbation theory.

### 5.1.2 The 8-band spin matrices

The code defines $8 \times 8$ spin matrices $\Sigma_x$, $\Sigma_y$, $\Sigma_z$ in the zincblende Kane basis. Their explicit form is initialized in `init_spin_matrices()` within `gfactor_functions.f90`. The basis ordering is:

| Index | State | Band |
|---|---|---|
| 1 | $\|3/2, +3/2\rangle$ | HH |
| 2 | $\|3/2, +1/2\rangle$ | LH |
| 3 | $\|3/2, -1/2\rangle$ | LH |
| 4 | $\|3/2, -3/2\rangle$ | HH |
| 5 | $\|1/2, +1/2\rangle$ | SO |
| 6 | $\|1/2, -1/2\rangle$ | SO |
| 7 | $\|S, +1/2\rangle$ | CB |
| 8 | $\|S, -1/2\rangle$ | CB |

The conduction-band block (rows/columns 7--8) yields the standard Pauli matrices:

$$
\Sigma_z^{(CB)} = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}, \quad
\Sigma_x^{(CB)} = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad
\Sigma_y^{(CB)} = \begin{pmatrix} 0 & i \\ -i & 0 \end{pmatrix}.
$$

The valence-band and cross-band blocks contain Clebsch--Gordan coefficients arising from the total angular momentum $J = 3/2$ (HH, LH) and $J = 1/2$ (SO) representations. For example,

$$
(\Sigma_z)_{11} = 1, \quad (\Sigma_z)_{22} = \tfrac{1}{3}, \quad
(\Sigma_z)_{33} = -\tfrac{1}{3}, \quad (\Sigma_z)_{44} = -1.
$$

The off-diagonal VB--SO elements such as $(\Sigma_z)_{25} = -\frac{2i\sqrt{2}}{3}$ arise from the coupling between the $J=3/2$ and $J=1/2$ manifolds. These off-diagonal elements are essential: they generate the inter-band contributions to $\Delta g$ when the VB and SO states act as intermediate states in the Lowdin partitioning.

![Zeeman splitting of the conduction band doublet](../figures/gfactor_zeeman.png)

### 5.1.3 Second-order Lowdin partitioning

#### Motivation

The 8-band k.p Hamiltonian at a given $\mathbf{k}$ point produces eigenstates that are linear combinations of all 8 basis states. For the g-factor, we are interested in the **Kramers doublet** of the conduction band (or a valence subband) near $\mathbf{k}=0$. Rather than working with the full $8N \times 8N$ problem (for QWs with $N$ finite-difference grid points), we can project the Zeeman interaction onto the $2 \times 2$ subspace of interest using **Lowdin partitioning**.

#### P and Q subspaces

We partition the Hilbert space into:

- **P subspace** (primary): the Kramers doublet of interest. For the CB, this is the pair of spin-split conduction subbands $\{|\psi_n^+\rangle, |\psi_n^-\rangle\}$ at the band index `bandIdx` (and its Kramers partner at `bandIdx+1`).

- **Q subspace** (remote): all other states -- the remaining conduction subbands, all valence subbands, and all split-off subbands.

The effective spin Hamiltonian in the P subspace, to second order in the magnetic perturbation, takes the form

$$
(H_{\mathrm{eff}})_{nm} = E_n \, \delta_{nm} + \sum_{l \in Q}
\frac{\langle n | V_\alpha | l \rangle \langle l | V_\beta | m \rangle
      - \langle n | V_\beta | l \rangle \langle l | V_\alpha | m \rangle}
     {E_n - E_l + E_m - E_l},
$$

where $V_\alpha$ and $V_\beta$ are the momentum perturbation operators corresponding to the two spatial directions transverse to the magnetic field component being computed, and the sum runs over all remote (Q-subspace) intermediate states $|l\rangle$.

#### Direction mapping

The cross-product structure of the angular momentum leads to a specific mapping between the tensor component $d$ (corresponding to $g_x$, $g_y$, $g_z$) and the pair of perturbation directions:

| Tensor component $d$ | Magnetic field direction | $\mathrm{mod}_1$ | $\mathrm{mod}_2$ |
|---|---|---|---|
| 1 ($g_x$) | $x$ | $y$ (2) | $z$ (3) |
| 2 ($g_y$) | $y$ | $z$ (3) | $x$ (1) |
| 3 ($g_z$) | $z$ | $x$ (1) | $y$ (2) |

This implements the relation $\mathbf{L} = \mathbf{r} \times \mathbf{p}$. For the $z$-component, $L_z = x p_y - y p_x$, so the perturbation operators are along $x$ (mod1=1) and $y$ (mod2=2). The antisymmetric combination $P_{\alpha\beta} = P_1 P_2 - P_2 P_1$ in the numerator captures the commutator structure of the angular momentum.

### 5.1.4 The g-tensor formula

The code computes the $2 \times 2$ g-tensor for each spatial direction $d$ as

$$
G_{ij}^{(d)} = -\frac{i}{\hbar^2/(2m_0)} \sum_{l \in Q}
\frac{P_1(n,l) \, P_2(l,m) - P_2(n,l) \, P_1(l,m)}{E_n + E_m - 2E_l}
\;-\; \frac{g_0}{2}\, \Sigma^{(d)}_{ij},
$$

where:

- $n, m \in \{n_0, n_0+1\}$ are the two states of the Kramers doublet (`bandIdx` and `bandIdx+1`),
- $l$ runs over all Q-subspace intermediate states,
- $P_1(n,l) = \langle \psi_n | H_{\mathrm{mod}_1} | \psi_l \rangle$ is the momentum matrix element in the first transverse direction,
- $P_2(l,m) = \langle \psi_l | H_{\mathrm{mod}_2} | \psi_m \rangle$ is the momentum matrix element in the second transverse direction,
- $\Sigma^{(d)}_{ij} = \langle \psi_{n_0+i-1} | \Sigma_d | \psi_{n_0+j-1} \rangle$ is the spin matrix element between the doublet states,
- $\hbar^2/(2m_0) \approx 3.81\;\mathrm{eV\,\AA^2}$ (the constant `hbar2O2m0` in `defs.f90`).

The effective g-factor along direction $d$ is extracted as the smallest eigenvalue of $2G^{(d)}$ (with sign). The code diagonalizes each $2 \times 2$ tensor slice with LAPACK's `zheev` and reports $g_d = 2\lambda_{\min}^{(d)}$.

### 5.1.5 Inter-band and intra-band contributions

The sum over intermediate states $l \in Q$ naturally separates into two physically distinct parts:

**Inter-band contributions** (CB $\leftrightarrow$ VB or VB $\leftrightarrow$ CB). For the CB g-factor, the intermediate states $|l\rangle$ run over all valence subbands. The energy denominator is

$$
E_{\mathrm{CB}}(n) + E_{\mathrm{CB}}(m) - 2E_{\mathrm{VB}}(l),
$$

which is dominated by the band gap $E_g$. For narrow-gap materials (small $E_g$), this denominator is small and the inter-band contribution is large. This is why InSb ($E_g \approx 0.17\;\mathrm{eV}$) has $|g^*| \approx 50$ while GaAs ($E_g \approx 1.52\;\mathrm{eV}$) has $g^* \approx -0.44$.

The dominant inter-band term can be estimated analytically. Keeping only the $P$ matrix element (Kane parameter) coupling the CB $s$-like states to the VB $p$-like states, Roth's formula gives

$$
\Delta g \approx -\frac{2E_P \, \Delta_{\mathrm{SO}}}{3E_g(E_g + \Delta_{\mathrm{SO}})},
$$

where $E_P = 2m_0 P^2/\hbar^2$ is the Kane energy and $\Delta_{\mathrm{SO}}$ is the spin-orbit splitting. This is the classic result: $|\Delta g|$ grows as $E_g$ shrinks.

**Intra-band contributions** (CB $\leftrightarrow$ CB or VB $\leftrightarrow$ VB). The intermediate states run over the other conduction subbands (excluding the doublet states themselves, i.e., $l \neq n, m$). The energy denominator is now

$$
E_{\mathrm{CB}}(n) + E_{\mathrm{CB}}(m) - 2E_{\mathrm{CB}}(l),
$$

which involves subband spacings rather than the band gap. In quantum wells and wires with strong confinement, these spacings can be comparable to $E_g$ for narrow structures, making intra-band corrections significant.

The code explicitly separates these two sums (see the two nested loops over `numvb` and `numcb` in `gfactorCalculation`), computing four momentum matrix elements per intermediate state:

$$
P_1(n,l),\; P_2(l,m),\; P_2(n,l),\; P_1(l,m).
$$

The antisymmetric combination $P_1 P_2 - P_2 P_1$ ensures that only the angular-momentum-like part of the coupling contributes (the symmetric part cancels). Contributions with near-zero energy denominators ($|\mathrm{denom}| < 10^{-7}\;\mathrm{eV}$) are skipped with a warning, preventing numerical divergence.

### 5.1.6 The momentum matrix element $P_\mathrm{ele}$

The quantity $P_\mathrm{ele} = \langle \psi_a | H_{\mathrm{pert}} | \psi_b \rangle$ is the matrix element of the **linear-k perturbation Hamiltonian** between two eigenstates. It is computed by building the Hamiltonian with a unit wave vector in a single direction and evaluating the overlap.

#### The `g` flag in Hamiltonian construction

When the optional argument `g='g'` is passed to `ZB8bandBulk` or `ZB8bandQW`, the Hamiltonian is constructed with **only the momentum (inter-band coupling) terms** active, while the intra-band Luttinger parameters are set to zero:

$$
\gamma_1 = 0, \quad \gamma_2 = 0, \quad \gamma_3 = 0, \quad A = 0.
$$

Only the Kane parameter $P$ is retained. This extracts the purely linear-$k$ part of the Hamiltonian, which is proportional to the crystal momentum operator. The perturbation Hamiltonian for direction $\alpha$ is

$$
H^{(\alpha)}_{\mathrm{pert}} = \left.\frac{\partial H}{\partial k_\alpha}\right|_{\mathbf{k}=0},
$$

evaluated by setting $k_\alpha = 1$ (unit perturbation) and all other components to zero.

In the 8-band bulk Hamiltonian, the linear-$k$ terms appear as:

$$
P_+ = \frac{P}{\sqrt{2}}(k_x + ik_y), \quad
P_- = \frac{P}{\sqrt{2}}(k_x - ik_y), \quad
P_z = P k_z.
$$

These couple the CB ($s$-like, bands 7--8) to the VB and SO bands ($p$-like, bands 1--6). With `g='g'` and $k_z = 1$, only $P_z$ survives, giving the $z$-direction perturbation. Similarly, $k_x = 1$ activates $P_+$ and $P_-$.

#### Wire mode: `g1`, `g2`, `g3` flags

For quantum wires (`confinement=2`), the perturbation Hamiltonian must account for the fact that two directions ($x$, $y$) involve spatial gradients ($d/dx$, $d/dy$) on a 2D finite-difference grid, while the third ($z$) is the free propagation direction with a simple $k_z$ perturbation. The code uses three separate flags:

| Flag | Direction | Mechanism |
|---|---|---|
| `g='g1'` | $x$ (confinement) | 2D gradient operator $d/dx$ applied to $P$-dependent blocks |
| `g='g2'` | $y$ (confinement) | 2D gradient operator $d/dy$ applied to $P$-dependent blocks |
| `g='g3'` | $z$ (free axis) | Unit $k_z$ perturbation (same topology as QW `g='g'`) |

All three produce a sparse CSR matrix that is applied to the eigenstate via `csr_spmv`, avoiding dense matrix construction for the large wire Hamiltonian.

### 5.1.7 Spin matrix elements

The spin matrix element between two eigenstates is

$$
\Sigma^{(d)}_{nm} = \langle \psi_n | \Sigma_d | \psi_m \rangle,
$$

where $\Sigma_d$ is the $8 \times 8$ spin matrix for direction $d$.

**Bulk** (single layer, `nlayers = 1`): The eigenstates are 8-component vectors. The matrix element is simply $\Sigma^{(d)}_{nm} = \psi_n^\dagger \Sigma_d \psi_m$, computed via `zgemv` followed by `zdotc`.

**Quantum well** (`nlayers > 1`): The eigenstates are $8N$-component vectors, where $N$ is the number of FD grid points (`fdStep`). At each grid point $i$, the 8-band components are extracted as $c_j^{(n)}(i) = \psi_n(i + (j-1)N)$ for $j = 1,\ldots,8$. The spin matrix element becomes a spatial integral:

$$
\Sigma^{(d)}_{nm} = \frac{1}{\Delta z}\int_{z_1}^{z_2}
\sum_{j,j'=1}^{8} \bar{c}_j^{(n)}(z) (\Sigma_d)_{jj'} c_{j'}^{(m)}(z) \, dz,
$$

evaluated numerically by Simpson's rule (`simpson(v, startz, endz) / dz`) over the FD grid.

**Wire** (`confinement=2`): The eigenstates are $8 N_\mathrm{grid}$-component vectors, where $N_\mathrm{grid} = n_x \times n_y$. At each grid point $(i_x, i_y)$ with flat index $f = (i_y - 1)n_x + i_x$, the 8-band components are $c_j^{(n)}(f) = \psi_n((j-1)N_\mathrm{grid} + f)$. The integration uses simple summation with uniform grid spacing:

$$
\Sigma^{(d)}_{nm} = \Delta x \, \Delta y \sum_{f=1}^{N_\mathrm{grid}}
\sum_{j,j'=1}^{8} \bar{c}_j^{(n)}(f) (\Sigma_d)_{jj'} c_{j'}^{(m)}(f).
$$

---

## 5.2 In the Code

### 5.2.1 Module structure

The g-factor computation is implemented in `src/physics/gfactor_functions.f90` (the `gfactorFunctions` module), with the executable entry point in `src/apps/main_gfactor.f90` (program `gfactor`).

| File | Role |
|---|---|
| `src/physics/gfactor_functions.f90` | Spin matrices, $\sigma$ elements, $P_\mathrm{ele}$ calculation, Lowdin sum, tensor assembly. Contains both QW and wire variants. |
| `src/apps/main_gfactor.f90` | Reads input, builds/diagonalizes Hamiltonian, calls g-factor routines, writes output. |

Key subroutines:

| Subroutine | Purpose |
|---|---|
| `init_spin_matrices()` | Initializes $8 \times 8$ $\Sigma_x, \Sigma_y, \Sigma_z$ on first call (cached). |
| `sigmaElem(state1, state2, dir, ...)` | QW/bulk: computes $\langle \psi_1 | \Sigma_\mathrm{dir} | \psi_2 \rangle$ with Simpson integration. |
| `sigmaElem_2d(state1, state2, dir, grid)` | Wire: same with uniform $dx \cdot dy$ summation over 2D grid. |
| `set_perturbation_direction(d, smallk)` | Sets a unit wave vector along direction $d \in \{1,2,3\}$. |
| `pMatrixEleCalc(Pele, d, ...)` | QW/bulk: computes $\langle a | H^{(d)}_\mathrm{pert} | b \rangle$ using dense or sparse Hamiltonian. |
| `pMatrixEleCalc_2d(Pele, d, ...)` | Wire: builds per-direction CSR perturbation with `g1`/`g2`/`g3` flags, uses `csr_spmv`. |
| `compute_pele(Pele, ...)` | Dispatcher that selects `pMatrixEleCalc` with correct arguments for bulk/QW. |
| `compute_pele_2d(Pele, ...)` | Dispatcher for wire mode, calling `pMatrixEleCalc_2d`. |
| `gfactorCalculation(tensor, ...)` | QW/bulk: full Lowdin partitioning, inter-band + intra-band sums, tensor assembly. |
| `gfactorCalculation_wire(tensor, ...)` | Wire: same algorithm with 2D sparse operators. |
| `compute_optical_matrix_wire(...)` | Wire: optical transition matrix elements for all CB-VB pairs. |

### 5.2.2 Program flow

The g-factor calculation is driven by program `gfactor`. The workflow is:

1. **Read input** via `read_and_setup` (shared `input_parser` module).
2. **Validate**: g-factor requires $k=0$ only (`waveVector: k0`, `waveVectorStep: 0`).
3. **Build Hamiltonian** at $\mathbf{k} = 0$ and diagonalize fully:
   - Bulk/QW: dense LAPACK (`zheev`/`zheevd`)
   - Wire: sparse FEAST eigensolver (`solve_sparse_evp`)
4. **Sort eigenstates**: VB states stored in descending energy order, CB states in ascending order.
5. **Select doublet**: the Kramers pair at `bandIdx` and `bandIdx+1`.
6. **Compute spin matrix** $\Sigma^{(d)}_{ij}$ between doublet states (3 directions).
7. **Loop over directions** $d = 1,2,3$: for each, select the perturbation pair $(\mathrm{mod}_1, \mathrm{mod}_2)$, build the perturbation Hamiltonians, and accumulate the momentum matrix element sums over all intermediate states.
8. **Assemble g-tensor** and diagonalize each $2 \times 2$ slice to extract $g_x, g_y, g_z$.
9. **Write output** to `output/gfactor.dat`.

### 5.2.3 CB vs VB g-factor

The `whichBand` flag selects which band edge to compute:

- `whichBand = 0` (CB): The P subspace is the CB Kramers doublet. Intermediate states are all VB subbands (inter-band) and all other CB subbands excluding the doublet (intra-band).

- `whichBand = 1` (VB): The P subspace is a VB Kramers doublet. Intermediate states are all CB subbands (inter-band) and all other VB subbands excluding the doublet (intra-band).

The logic is symmetric: the two code blocks in `gfactorCalculation` are structurally identical, with CB and VB arrays swapped. The `bandIdx` parameter selects which subband doublet to target (1 = ground state, 2 = first excited, etc.).

### 5.2.4 Sparse mode for quantum wells

For multi-layer QW structures, the Hamiltonian is large ($8N \times 8N$ with $N = \mathtt{fdStep}$). The code builds the perturbation Hamiltonians in CSR sparse format via `ZB8bandQW(..., sparse=.True., HT_csr=..., g='g')`. The momentum matrix elements are then computed with CSR sparse matrix-vector multiplication (`csr_spmv`), avoiding the need to store dense $8N \times 8N$ perturbation matrices. The two perturbation CSR matrices (`HT_csr_mod1` and `HT_csr_mod2`) are built once per direction $d$ and reused across all intermediate states.

### 5.2.5 Wire mode g-factor

For wire mode (`confinement=2`), the program follows a different branch in `main_gfactor.f90`:

1. Build the 2D sparse k.p terms via `confinementInitialization_2d`.
2. Optionally compute and apply strain via `compute_strain` + `apply_pikus_bir`.
3. Build the full wire Hamiltonian at $k_z = 0$ with `ZB8bandGeneralized`.
4. Solve with FEAST (`solve_sparse_evp`) using an auto-computed energy window.
5. Extract CB/VB states from sorted eigenvalues.
6. Call `gfactorCalculation_wire`, which uses `sigmaElem_2d` for spin integration (uniform $dx \cdot dy$ summation over the 2D grid) and `pMatrixEleCalc_2d` for momentum matrix elements (builds per-direction CSR perturbation matrices with `g1`/`g2`/`g3` flags).

The wire g-factor computes all three spatial directions ($g_x$, $g_y$, $g_z$) in a single run, reflecting the full anisotropy of the wire cross-section. For a cylindrical wire, $g_x = g_y$ by symmetry; for a rectangular wire, all three components can differ.

### 5.2.6 Tensor assembly

After the momentum matrix element sums and spin terms are computed, the code assembles the full g-tensor in two steps (lines 637--638 of `gfactor_functions.f90`):

```fortran
tensor(:,:,:) = -cmplx(0, 1, dp) * tensor(:,:,:) / hbar2O2m0
tensor(:,:,:) = tensor(:,:,:) - (g_free / 2) * sigma(:,:,:)
```

The first line multiplies the accumulated $\sum_l (P_1 P_2 - P_3 P_4)/\mathrm{denom}$ by $-i/(\hbar^2/2m_0)$, converting from the code's natural units to the physical g-tensor convention. The second line subtracts the free-electron spin contribution $(g_0/2)\Sigma$.

The final g-factor values are obtained by diagonalizing each $2 \times 2$ slice of the tensor in `main_gfactor.f90`:

```fortran
call zheev('N', 'U', 2, tensor(:,:,d), 2, gfac(:,d), work, lwork, rwork, info)
g_eff(d) = 2 * gfac(1, d)
```

The factor of 2 converts the eigenvalue convention (half the Zeeman splitting) to the standard g-factor definition. The three components are written to `output/gfactor.dat`.

### 5.2.7 Memory layout for eigenstates

Understanding the memory layout is essential for following the code:

**QW** (dense, column-major): The eigenstate vector for band $n$ has components $\psi_n(i + (j-1) \cdot N_\mathrm{fd})$ where $i = 1,\ldots,N_\mathrm{fd}$ is the spatial grid index and $j = 1,\ldots,8$ is the band index. This is the **band-interleaved** layout: all spatial points for band 1 come first, then all points for band 2, etc.

**Wire** (flat 2D, column-major): The eigenstate has components $\psi_n((j-1) \cdot N_\mathrm{grid} + f)$ where $f = (i_y - 1) n_x + i_x$ is the flat grid index and $j = 1,\ldots,8$ is the band index. The grid traverses $x$ fastest (inner loop), $y$ slowest.

---

## 5.3 Computed Example

### 5.3.1 Bulk GaAs conduction band g-factor

The simplest case is a bulk semiconductor. The input configuration (`gfactor_bulk_gaas_cb.cfg`) is:

```ini
waveVector: k0
waveVectorMax: 0.1
waveVectorStep: 0
confinement:  0
FDstep: 1
FDorder: 2
numLayers:  1
material1: GaAs
numcb: 2
numvb: 6
ExternalField: 0  EF
EFParams: 0.0005
whichBand: 0
bandIdx: 1
```

Key observations:
- `confinement: 0` selects bulk mode (8x8 Hamiltonian).
- `FDstep: 1` gives a single spatial point.
- `numcb: 2` and `numvb: 6` give the full 8-band basis.
- `whichBand: 0` selects the conduction band.
- `bandIdx: 1` selects the ground-state Kramers doublet (CB bands 7--8).

For bulk GaAs, the expected result is $g_{\text{eff}} \approx -0.44$ in all three directions (isotropic by cubic symmetry). The breakdown into contributions:

| Contribution | Approximate value |
|---|---|
| Free-electron $g_0$ | +2.002 |
| Inter-band $\Delta g$ (VB) | -2.44 |
| **Total** | **-0.44** |

The inter-band contribution dominates and is negative, reflecting the strong coupling to the valence band through the Kane parameter $P$.

### 5.3.2 QW conduction band g-factor: InAs/GaSb/AlSb type-II system

A more interesting case is the quantum well configuration (`gfactor_qw_cb.cfg`):

```ini
waveVector: k0
waveVectorMax: 0.1
waveVectorStep: 0
confinement:  1
FDstep: 101
FDorder: 2
numLayers:  3
material1: AlSbW -250  250 0
material2: GaSbW -135  135 0.2414
material3: InAsW  -35   35 -0.0914
numcb: 32
numvb: 32
ExternalField: 0  EF
EFParams: 0.0005
whichBand: 0
bandIdx: 1
```

This defines an InAs/GaSb broken-gap quantum well with AlSb barriers:
- 101 FD grid points, second-order finite differences.
- Three material layers with Winkler parameters (the `W` suffix selects the Winkler parameter set).
- The band offsets are specified as the third parameter on each `material` line.
- `numcb: 32` and `numvb: 32` are multiplied by `fdStep` internally to give 202 CB and 606 VB states for the full $808 \times 808$ QW Hamiltonian.

The g-factor output is written to `output/gfactor.dat` as three values:

```
  gx_value  gy_value  gz_value
```

For a QW grown along $z$, the growth-direction g-factor $g_z$ typically differs from the in-plane $g_x = g_y$ due to the reduced symmetry. The inter-band contribution is enhanced for the in-plane directions because the confinement-induced VB mixing couples more strongly to the in-plane momentum components.

The terminal output also shows the intermediate quantities:

1. **Sigma matrix** -- the spin Zeeman contribution ($2 \times 2$ for each direction).
2. **Tensor** -- the full orbital contribution before normalization ($2 \times 2$ complex for each direction).
3. **Eigenvalues** of each $2 \times 2$ block, multiplied by 2 to give $g_x$, $g_y$, $g_z$.

### 5.3.3 Wire g-factor: all three spatial directions

For a quantum wire with confinement in the $x$-$y$ plane and free propagation along $z$, the g-factor is generally anisotropic in all three directions. The wire computation uses `gfactorCalculation_wire`, which:

1. Solves the wire Hamiltonian at $k_z = 0$ using FEAST.
2. Computes `sigmaElem_2d` over the full 2D grid ($nx \times ny$ points).
3. Builds perturbation Hamiltonians for all three directions using `g='g1'`, `g='g2'`, `g='g3'`.
4. Sums the Lowdin contributions over VB and CB intermediates.

The three g-factor components $g_x$, $g_y$, $g_z$ reflect the broken rotational symmetry: $g_x \neq g_y$ for a rectangular wire, while $g_x = g_y$ is recovered for a circular cross-section. The $g_z$ component (along the free wire axis) is typically closest to the QW limit because the wave function is unconfined in this direction.

![g-factor components for different confinement geometries](../figures/gfactor_components.png)

---

## 5.4 Discussion

### 5.4.1 Physical interpretation

The effective g-factor is a probe of the **band mixing** in the target doublet. Several physical trends emerge from the Lowdin formula:

1. **Band gap dependence.** The dominant inter-band contribution scales as $\Delta g \propto -E_P / E_g$. Narrow-gap semiconductors (InSb, InAs) have large negative g-factors; wide-gap materials (GaN, ZnO) have g-factors close to $g_0 \approx 2$.

2. **Confinement effects.** In a quantum well, the quantized subband energies change the effective energy denominators in the Lowdin sum. The ground-state CB subband is pushed up in energy, increasing $E_n - E_l$ for VB intermediates and thus reducing $|\Delta g|$. This is why quantum-confined g-factors are generally **closer to $g_0$** than their bulk counterparts.

3. **Spin-orbit coupling.** The SO band (bands 5--6) contributes to $\Delta g$ through the split-off energy $\Delta_{\mathrm{SO}}$ in the denominator. Materials with small $\Delta_{\mathrm{SO}}$ (InSb) have enhanced SO contributions.

4. **Anisotropy.** Bulk zincblende has cubic symmetry, giving $g_x = g_y = g_z$. Quantum wells break this to $C_{\infty v}$ (or lower), splitting $g_\parallel = g_z$ from $g_\perp = g_x = g_y$. Quantum wires can further split $g_x \neq g_y$, producing a fully anisotropic g-tensor.

### 5.4.2 Limitations and accuracy

Several factors affect the accuracy of the computed g-factor:

1. **Band basis truncation.** The 8-band model includes only the CB, HH, LH, and SO bands. Remote bands (higher conduction bands, deeper valence bands) are absent. Their contribution to $\Delta g$ is estimated to be of order 0.1--0.5, which is significant for materials with $|g|$ near zero (e.g., GaAs).

2. **Converged subband count.** The Lowdin sum requires enough intermediate states to converge. Using too few `numcb` or `numvb` subbands truncates the sum and underestimates $|\Delta g|$. A convergence test -- increasing `numcb`/`numvb` until $g$ stabilizes -- is recommended.

3. **Finite-difference order.** The perturbation Hamiltonian inherits the FD discretization of the QW or wire Hamiltonian. Higher FD orders (4th, 6th, ...) improve the momentum matrix elements and thus the g-factor accuracy. Second-order ($FDorder=2$) is the default and is adequate for most applications.

4. **Energy denominator near zero.** When two states are nearly degenerate ($E_n \approx E_l$), the denominator $(E_n - E_l) + (E_m - E_l)$ approaches zero. The code skips such terms with a warning:

   ```
   WARNING: near-zero energy denominator, n=... m=... l=... denom=...
   ```

   This is correct behavior: nearly degenerate states should be in the P subspace (treated exactly), not in the Q subspace (perturbatively). The user should increase the doublet range or adjust `bandIdx` if many warnings appear.

5. **Spin matrix phase convention.** The spin matrices `SIGMA_X/Y/Z` in `gfactor_functions.f90` follow the phase convention of `ZB8bandBulk` in `hamiltonianConstructor.f90`. Off-diagonal VB-SO elements may differ from the Chuang & Chang or Winkler Table 2.3 conventions by factors of $\pm i$. This does not affect the final g-factor (which is gauge-invariant), but the intermediate sigma matrix printout may look different from textbook values.

### 5.4.3 Tips for accurate g-factor calculations

- **Always start from bulk.** Verify that your bulk g-factor matches published values before moving to confined structures. This validates the material parameters.
- **Use enough grid points.** The g-factor is sensitive to the wave function shape. For QW calculations, at least 80--100 FD points across the well is recommended.
- **Check convergence in `numcb`/`numvb`.** Double the number of subbands and verify that $g$ changes by less than 1%.
- **Monitor warnings.** More than a few "near-zero denominator" warnings suggests that the Lowdin partitioning is not well-separated and results may be unreliable.
- **For wires, use sufficient grid resolution** in both $x$ and $y$. The g-factor anisotropy $g_x - g_y$ is sensitive to the wire cross-section shape, which is resolved by the 2D grid.
- **The `g_free` constant** in `defs.f90` is set to the CODATA value 2.00231. If comparing with older literature that uses $g_0 = 2.0$, the difference of 0.0023 may be significant for materials with $|g| \approx 0$.

### 5.4.4 Connection to other chapters

- The momentum matrix elements used in the g-factor computation (this chapter) are the same operators that determine optical transition strengths (Chapter 6). The inter-band coupling that renormalizes $g$ also sets the oscillator strength.
- The self-consistent Schrodinger-Poisson solver (Chapter 7) modifies the band profiles and thus the wave functions, which in turn affect the g-factor. For doped structures, running the SP loop before the g-factor calculation is essential.
- Strain (Chapter 4) changes the band offsets and mixes the VB states, modifying the inter-band contributions to $\Delta g$.
