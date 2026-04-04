# Chapter 05: Landau g-Factor from Second-Order Lowdin Partitioning

## 1. Theory

### 1.1 The Zeeman Effect in Semiconductors

When a magnetic field $\mathbf{B}$ is applied to a semiconductor, the electron
states split according to their spin and orbital angular momentum. In the free
electron, this gives the familiar Zeeman Hamiltonian

$$
H_Z = \frac{\mu_B}{2} \, \mathbf{B} \cdot \left( g_0 \, \boldsymbol{\sigma} + 2\mathbf{L} \right),
$$

where $\mu_B = e\hbar/(2m_0)$ is the Bohr magneton, $g_0 \approx 2.00231$ is
the free-electron g-factor (the code uses the CODATA value `g_free` from
`defs.f90`), $\boldsymbol{\sigma} = (\sigma_x, \sigma_y, \sigma_z)$ are the
Pauli matrices, and $\mathbf{L}$ is the orbital angular momentum operator. The
first term describes the coupling of $\mathbf{B}$ to the electron spin; the
second describes the coupling to orbital angular momentum.

In a crystal, the orbital angular momentum is no longer a good quantum number.
The spin-orbit coupling inherited from the atomic $p$-like valence band states
mixes spin and orbital character, producing an **effective g-factor** that
deviates significantly from $g_0$. For conduction-band electrons in narrow-gap
semiconductors like InSb or InAs, the g-factor can be $|g^*| \gg g_0$ and even
change sign. The deviation $\Delta g = g^* - g_0$ is almost entirely due to
inter-band coupling, which we compute via second-order perturbation theory.

### 1.2 The 8-Band Spin Matrices

The code defines $8\times 8$ spin matrices $\Sigma_x$, $\Sigma_y$, $\Sigma_z$
in the zincblende Kane basis. Their explicit form is initialized in
`init_spin_matrices()` within `gfactor_functions.f90`. The basis ordering is:

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

The conduction-band block (rows/columns 7--8) yields the standard Pauli
matrices:

$$
\Sigma_z^{(CB)} = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}, \quad
\Sigma_x^{(CB)} = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad
\Sigma_y^{(CB)} = \begin{pmatrix} 0 & i \\ -i & 0 \end{pmatrix}.
$$

The valence-band and cross-band blocks contain Clebsch--Gordan coefficients
arising from the total angular momentum $J = 3/2$ (HH, LH) and $J = 1/2$ (SO)
representations. For example,

$$
(\Sigma_z)_{11} = 1, \quad (\Sigma_z)_{22} = \tfrac{1}{3}, \quad
(\Sigma_z)_{33} = -\tfrac{1}{3}, \quad (\Sigma_z)_{44} = -1.
$$

The off-diagonal VB--SO elements such as $(\Sigma_z)_{25} = -\frac{2i\sqrt{2}}{3}$
arise from the coupling between the $J=3/2$ and $J=1/2$ manifolds. These
off-diagonal elements are essential: they generate the inter-band contributions
to $\Delta g$ when the VB and SO states act as intermediate states in the
Lowdin partitioning.

### 1.3 Second-Order Lowdin Partitioning

#### Motivation

The 8-band k.p Hamiltonian at a given $\mathbf{k}$ point produces eigenstates
that are linear combinations of all 8 basis states. For the g-factor, we are
interested in the **Kramers doublet** of the conduction band (or a valence
subband) near $\mathbf{k}=0$. Rather than working with the full $8N\times 8N$
problem (for QWs with $N$ finite-difference grid points), we can project the
Zeeman interaction onto the $2\times 2$ subspace of interest using
**Lowdin partitioning**.

#### P and Q Subspaces

We partition the Hilbert space into:

- **P subspace** (primary): the Kramers doublet of interest. For the CB, this is
  the pair of spin-split conduction subbands $\{|\psi_n^+\rangle, |\psi_n^-\rangle\}$
  at the band index `bandIdx` (and its Kramers partner at `bandIdx+1`).

- **Q subspace** (remote): all other states --- the remaining conduction subbands,
  all valence subbands, and all split-off subbands.

The effective spin Hamiltonian in the P subspace, to second order in the
magnetic perturbation, takes the form

$$
(H_{\mathrm{eff}})_{nm} = E_n \, \delta_{nm} + \sum_{l \in Q}
\frac{\langle n | V_\alpha | l \rangle \langle l | V_\beta | m \rangle
      - \langle n | V_\beta | l \rangle \langle l | V_\alpha | m \rangle}
     {E_n - E_l + E_m - E_l},
$$

where $V_\alpha$ and $V_\beta$ are the momentum perturbation operators
corresponding to the two spatial directions transverse to the magnetic field
component being computed, and the sum runs over all remote (Q-subspace)
intermediate states $|l\rangle$.

#### Direction Mapping

The cross-product structure of the angular momentum leads to a specific mapping
between the tensor component $d$ (corresponding to $g_x$, $g_y$, $g_z$) and
the pair of perturbation directions:

| Tensor component $d$ | Magnetic field direction | $\mathrm{mod}_1$ | $\mathrm{mod}_2$ |
|---|---|---|---|
| 1 ($g_x$) | $x$ | $y$ (2) | $z$ (3) |
| 2 ($g_y$) | $y$ | $z$ (3) | $x$ (1) |
| 3 ($g_z$) | $z$ | $x$ (1) | $y$ (2) |

This implements the relation $\mathbf{L} = \mathbf{r} \times \mathbf{p}$. For
the $z$-component, $L_z = x p_y - y p_x$, so the perturbation operators are
along $x$ (mod1=1) and $y$ (mod2=2). The antisymmetric combination
$P_{\alpha\beta} = P_1 P_2 - P_2 P_1$ in the numerator captures the
commutator structure of the angular momentum.

### 1.4 The g-Tensor Formula

The code computes the $2\times 2$ g-tensor for each spatial direction $d$ as

$$
G_{ij}^{(d)} = -\frac{i}{\hbar^2/(2m_0)} \sum_{l \in Q}
\frac{P_1(n,l) \, P_2(l,m) - P_2(n,l) \, P_1(l,m)}{E_n + E_m - 2E_l}
\;-\; \frac{g_0}{2}\, \Sigma^{(d)}_{ij},
$$

where:

- $n, m \in \{n_0, n_0+1\}$ are the two states of the Kramers doublet
  (`bandIdx` and `bandIdx+1`),
- $l$ runs over all Q-subspace intermediate states,
- $P_1(n,l) = \langle \psi_n | H_{\mathrm{mod}_1} | \psi_l \rangle$ is the
  momentum matrix element in the first transverse direction,
- $P_2(l,m) = \langle \psi_l | H_{\mathrm{mod}_2} | \psi_m \rangle$ is the
  momentum matrix element in the second transverse direction,
- $\Sigma^{(d)}_{ij} = \langle \psi_{n_0+i-1} | \Sigma_d | \psi_{n_0+j-1} \rangle$
  is the spin matrix element between the doublet states,
- $\hbar^2/(2m_0) \approx 3.81\;\mathrm{eV\,\AA^2}$ (the constant `hbar2O2m0`
  in `defs.f90`).

The effective g-factor along direction $d$ is extracted as the smallest
eigenvalue of $2G^{(d)}$ (with sign). The code diagonalizes each $2\times 2$
tensor slice with LAPACK's `zheev` and reports $g_d = 2\lambda_{\min}^{(d)}$.

### 1.5 Inter-Band and Intra-Band Contributions

The sum over intermediate states $l \in Q$ naturally separates into two
physically distinct parts:

**Inter-band contributions** (CB$\leftrightarrow$VB or VB$\leftrightarrow$CB).
For the CB g-factor, the intermediate states $|l\rangle$ run over all valence
subbands. The energy denominator is

$$
E_{\mathrm{CB}}(n) + E_{\mathrm{CB}}(m) - 2E_{\mathrm{VB}}(l),
$$

which is dominated by the band gap $E_g$. For narrow-gap materials (small
$E_g$), this denominator is small and the inter-band contribution is large.
This is why InSb ($E_g \approx 0.17\;\mathrm{eV}$) has $|g^*| \approx 50$
while GaAs ($E_g \approx 1.52\;\mathrm{eV}$) has $g^* \approx -0.44$.

The dominant inter-band term can be estimated analytically. Keeping only the
$P$ matrix element (Kane parameter) coupling the CB $s$-like states to the VB
$p$-like states, Roth's formula gives

$$
\Delta g \approx -\frac{2E_P \, \Delta_{\mathrm{SO}}}{3E_g(E_g + \Delta_{\mathrm{SO}})},
$$

where $E_P = 2m_0 P^2/\hbar^2$ is the Kane energy and $\Delta_{\mathrm{SO}}$
is the spin-orbit splitting. This is the classic result: $|\Delta g|$ grows as
$E_g$ shrinks.

**Intra-band contributions** (CB$\leftrightarrow$CB or VB$\leftrightarrow$VB).
The intermediate states run over the other conduction subbands (excluding the
doublet states themselves, i.e., $l \neq n, m$). The energy denominator is now

$$
E_{\mathrm{CB}}(n) + E_{\mathrm{CB}}(m) - 2E_{\mathrm{CB}}(l),
$$

which involves subband spacings rather than the band gap. In quantum wells and
wires with strong confinement, these spacings can be comparable to $E_g$ for
narrow structures, making intra-band corrections significant.

The code explicitly separates these two sums (see the two nested loops over
`numvb` and `numcb` in `gfactorCalculation`), computing four momentum matrix
elements per intermediate state:

$$
P_1(n,l),\; P_2(l,m),\; P_2(n,l),\; P_1(l,m).
$$

The antisymmetric combination $P_1 P_2 - P_2 P_1$ ensures that only the
angular-momentum-like part of the coupling contributes (the symmetric part
cancels). Contributions with near-zero energy denominators ($|\mathrm{denom}| <
10^{-7}\;\mathrm{eV}$) are skipped with a warning, preventing numerical
divergence.

### 1.6 The Momentum Matrix Element $P_\mathrm{ele}$

The quantity $P_\mathrm{ele} = \langle \psi_a | H_{\mathrm{pert}} | \psi_b \rangle$
is the matrix element of the **linear-k perturbation Hamiltonian** between two
eigenstates. It is computed by building the Hamiltonian with a unit wave vector
in a single direction and evaluating the overlap.

#### The `g` Flag in Hamiltonian Construction

When the optional argument `g='g'` is passed to `ZB8bandBulk` or `ZB8bandQW`,
the Hamiltonian is constructed with **only the momentum (inter-band coupling)
terms** active, while the intra-band Luttinger parameters are set to zero:

$$
\gamma_1 = 0, \quad \gamma_2 = 0, \quad \gamma_3 = 0, \quad A = 0.
$$

Only the Kane parameter $P$ is retained. This extracts the purely linear-$k$
part of the Hamiltonian, which is proportional to the crystal momentum
operator. The perturbation Hamiltonian for direction $\alpha$ is

$$
H^{(\alpha)}_{\mathrm{pert}} = \left.\frac{\partial H}{\partial k_\alpha}\right|_{\mathbf{k}=0},
$$

evaluated by setting $k_\alpha = 1$ (unit perturbation) and all other
components to zero.

In the 8-band bulk Hamiltonian, the linear-$k$ terms appear as:

$$
P_+ = \frac{P}{\sqrt{2}}(k_x + ik_y), \quad
P_- = \frac{P}{\sqrt{2}}(k_x - ik_y), \quad
P_z = P k_z.
$$

These couple the CB ($s$-like, bands 7--8) to the VB and SO bands ($p$-like,
bands 1--6). With `g='g'` and $k_z = 1$, only $P_z$ survives, giving the
$z$-direction perturbation. Similarly, $k_x = 1$ activates $P_+$ and $P_-$.

#### Wire Mode: `g1`, `g2`, `g3` Flags

For quantum wires (`confinement=2`), the perturbation Hamiltonian must account
for the fact that two directions ($x$, $y$) involve spatial gradients
($d/dx$, $d/dy$) on a 2D finite-difference grid, while the third ($z$) is the
free propagation direction with a simple $k_z$ perturbation. The code uses
three separate flags:

| Flag | Direction | Mechanism |
|---|---|---|
| `g='g1'` | $x$ (confinement) | Uses `kpterms_2d(12)` = $P \cdot d/dx$ and `kpterms_2d(14)` = $\gamma_3 \cdot d/dx$ |
| `g='g2'` | $y$ (confinement) | Uses `kpterms_2d(13)` = $P \cdot d/dy$ and `kpterms_2d(15)` = $\gamma_3 \cdot d/dy$ |
| `g='g3'` | $z$ (free axis) | Uses $k_z = 1$ with PP, PM, S, SC blocks (same topology as QW `g='g'`) |

The `g1` and `g2` modes build a CSR sparse perturbation Hamiltonian that
applies the precomputed $d/dx$ or $d/dy$ finite-difference operators to the
$P_z$ and $S$/$S^c$ blocks. The `g3` mode treats $k_z$ as the perturbation
variable, analogous to the QW case. All three produce a sparse matrix that is
applied to the eigenstate via CSR sparse matrix-vector multiplication
(`csr_spmv`), avoiding dense matrix construction for the large wire Hamiltonian.

### 1.7 Spin Matrix Elements: The $\sigma$ Term

The spin matrix element between two eigenstates is

$$
\Sigma^{(d)}_{nm} = \langle \psi_n | \Sigma_d | \psi_m \rangle,
$$

where $\Sigma_d$ is the $8\times 8$ spin matrix for direction $d$.

**Bulk** (single layer, `nlayers = 1`): The eigenstates are 8-component vectors.
The matrix element is simply $\Sigma^{(d)}_{nm} = \sum_{ij} \bar{c}_i^{(n)}
(\Sigma_d)_{ij} c_j^{(m)}$, computed via `zgemv` followed by `zdotc`.

**Quantum well** (`nlayers > 1`): The eigenstates are $8N$-component vectors,
where $N$ is the number of FD grid points (`fdStep`). At each grid point $i$,
the 8-band components are extracted as $c_j^{(n)}(i) = \psi_n(i + (j-1)N)$ for
$j = 1,\ldots,8$. The spin matrix element becomes a spatial integral:

$$
\Sigma^{(d)}_{nm} = \frac{1}{\Delta z}\int_{z_1}^{z_2}
\sum_{j,j'=1}^{8} \bar{c}_j^{(n)}(z) (\Sigma_d)_{jj'} c_{j'}^{(m)}(z) \, dz,
$$

evaluated numerically by Simpson's rule (`simpson(v, startz, endz) / dz`)
over the FD grid. The QW spin matrices are printed to stdout for verification.

**Wire** (`confinement=2`): The eigenstates are $8 N_\mathrm{grid}$-component
vectors, where $N_\mathrm{grid} = n_x \times n_y$. At each grid point
$(i_x, i_y)$ with flat index $f = (i_y - 1)n_x + i_x$, the 8-band components
are $c_j^{(n)}(f) = \psi_n((j-1)N_\mathrm{grid} + f)$. The integration uses
simple summation with uniform grid spacing:

$$
\Sigma^{(d)}_{nm} = \Delta x \, \Delta y \sum_{f=1}^{N_\mathrm{grid}}
\sum_{j,j'=1}^{8} \bar{c}_j^{(n)}(f) (\Sigma_d)_{jj'} c_{j'}^{(m)}(f).
$$

### 1.8 Assembling the g-Tensor

After the momentum matrix element sums and spin terms are computed, the code
assembles the full g-tensor in two steps (see lines 637--638 of
`gfactor_functions.f90`):

```fortran
tensor(:,:,:) = -cmplx(0,1,dp) * tensor(:,:,:) / hbar2O2m0
tensor(:,:,:) = tensor(:,:,:) - (g_free/2) * sigma(:,:,:)
```

The first line multiplies the accumulated
$\sum_l (P_1 P_2 - P_3 P_4)/\mathrm{denom}$ by $-i/(\hbar^2/2m_0)$,
converting from the code's natural units to the physical g-tensor convention.
The second line subtracts the free-electron spin contribution $(g_0/2)\Sigma$.

The final g-factor values are obtained by diagonalizing each $2\times 2$ slice
of the tensor. In `main_gfactor.f90`, the code uses `zheev` on each of the
three slices:

$$
g_d = 2\lambda_{\min}^{(d)}, \quad d \in \{x, y, z\},
$$

where $\lambda_{\min}^{(d)}$ is the smallest eigenvalue of the $2\times 2$
matrix $G^{(d)}$. For an isotropic bulk semiconductor, all three components
should be equal; for QWs, $g_z$ (growth direction) typically differs from
$g_x = g_y$ (in-plane); for wires, all three can differ.

The three components are written to `output/gfactor.dat`.

---

## 2. Implementation in the Code

### 2.1 Program Flow

The g-factor calculation is driven by the `gfactorCalculation` program
(`src/apps/main_gfactor.f90`). The workflow is:

1. **Read input** via `read_and_setup` (shared `input_parser` module).
2. **Validate**: g-factor requires $k=0$ only (`waveVector: k0`,
   `waveVectorStep: 0`).
3. **Build Hamiltonian** at $\mathbf{k} = 0$ and diagonalize fully
   (LAPACK `zheev`/`zheevd` for bulk/QW; FEAST for wire).
4. **Sort eigenstates**: VB states stored in descending energy order,
   CB states in ascending order.
5. **Select doublet**: the Kramers pair at `bandIdx` and `bandIdx+1`.
6. **Compute spin matrix** $\Sigma^{(d)}_{ij}$ between doublet states
   (3 directions).
7. **Loop over directions** $d = 1,2,3$: for each, select the perturbation
   pair $(\mathrm{mod}_1, \mathrm{mod}_2)$, build the perturbation
   Hamiltonians, and accumulate the momentum matrix element sums over all
   intermediate states.
8. **Assemble g-tensor** and diagonalize to extract $g_x, g_y, g_z$.

### 2.2 Conduction Band vs. Valence Band

The `whichBand` flag selects which band edge to compute:

- `whichBand = 0` (CB): The P subspace is the CB Kramers doublet.
  Intermediate states are all VB subbands (inter-band) and all other CB
  subbands excluding the doublet (intra-band).

- `whichBand = 1` (VB): The P subspace is a VB Kramers doublet.
  Intermediate states are all CB subbands (inter-band) and all other VB
  subbands excluding the doublet (intra-band).

The logic is symmetric: the two code blocks in `gfactorCalculation` are
structurally identical, with CB and VB arrays swapped. The `bandIdx` parameter
selects which subband doublet to target (1 = ground state, 2 = first excited,
etc.).

### 2.3 Sparse Mode for Quantum Wells

For multi-layer QW structures, the Hamiltonian is large ($8N \times 8N$ with
$N = \mathtt{fdStep}$). The code builds the perturbation Hamiltonians in
CSR sparse format via `ZB8bandQW(..., sparse=.True., HT_csr=..., g='g')`.
The momentum matrix elements are then computed with CSR sparse matrix-vector
multiplication (`csr_spmv`), avoiding the need to store dense $8N \times 8N$
perturbation matrices. The two perturbation CSR matrices (`HT_csr_mod1` and
`HT_csr_mod2`) are built once per direction $d$ and reused across all
intermediate states.

### 2.4 Wire Mode g-Factor

For wire mode (`confinement=2`), the program follows a different branch in
`main_gfactor.f90`:

1. Build the 2D sparse k.p terms via `confinementInitialization_2d`.
2. Optionally compute and apply strain via `compute_strain` +
   `apply_pikus_bir`.
3. Build the full wire Hamiltonian at $k_z = 0$ with `ZB8bandGeneralized`.
4. Solve with FEAST (`solve_sparse_evp`) using an auto-computed energy window.
5. Extract CB/VB states from sorted eigenvalues.
6. Call `gfactorCalculation_wire`, which uses `sigmaElem_2d` for spin
   integration (uniform $dx \cdot dy$ summation over the 2D grid) and
   `pMatrixEleCalc_2d` for momentum matrix elements (builds per-direction
   CSR perturbation matrices with `g1`/`g2`/`g3` flags).

The wire g-factor computes all three spatial directions ($g_x$, $g_y$, $g_z$)
in a single run, reflecting the full anisotropy of the wire cross-section.
For a cylindrical wire, $g_x = g_y$ by symmetry.

### 2.5 Optical Transitions (Wire Mode)

As a by-product of the g-factor calculation in wire mode, the code also
computes **optical transition matrix elements** between all CB-VB subband
pairs via `compute_optical_matrix_wire`. For each pair $(i, j)$:

$$
|\langle \psi_i^{\mathrm{CB}} | p_\alpha | \psi_j^{\mathrm{VB}} \rangle|^2,
\quad \alpha \in \{x, y, z\},
$$

and the oscillator strength

$$
f_{ij} = \frac{1}{\hbar^2/(2m_0)} \,
\frac{|\langle \psi_i | H_{p_x} | \psi_j \rangle|^2
    + |\langle \psi_i | H_{p_y} | \psi_j \rangle|^2
    + |\langle \psi_i | H_{p_z} | \psi_j \rangle|^2}
     {\Delta E_{ij}}.
$$

These are written to `output/optical_transitions.dat`.

---

## 3. Practical Example

### 3.1 Published Result: InSb Nanowire g-Factor

In Faria Junior, Campos *et al.*, Phys. Rev. B **97**, 245402 (2018)
([arXiv:1802.06734](https://arxiv.org/abs/1802.06734)), the authors computed
the Landau g-factor for cylindrical InSb nanowires grown along [111], showing
strong diameter-dependent g-factor renormalization due to quantum confinement.
The g-factor decreases from the bulk value ($\approx -51$) toward zero as the
wire diameter shrinks below about 30 nm.

The current code can reproduce the **quantum well** g-factor with the existing
QW mode. The wire g-factor requires `confinement=2` mode with appropriate
2D geometry input. The Winkler parameter set (materials ending in `W`, e.g.,
`InAsW`, `GaSbW`, `AlSbW`) uses InSb as the valence-band energy reference,
which is essential for consistent band offsets in heterostructures.

### 3.2 Example: InAs/GaSb/AlSb Quantum Well

The following `input.cfg` computes the conduction-band g-factor for an
InAs/GaSb/AlSb broken-gap quantum well:

```
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

Key points:

- `waveVector: k0` with `waveVectorStep: 0` --- g-factor is computed at
  $\Gamma$ only.
- `FDstep: 101` --- 101 FD grid points across the structure.
- Three materials with Winkler parameters: AlSbW barriers, GaSbW inner well,
  InAsW core. The third column in `materialN` is the valence-band offset
  relative to the reference (InSb).
- `whichBand: 0` --- compute the conduction-band g-factor.
- `bandIdx: 1` --- ground-state CB Kramers doublet.
- `ExternalField: 0 EF` with `EFParams: 0.0005` --- a small electric field
  (0.5 kV/cm) applied to break structural inversion symmetry (optional).
- `numcb: 32` and `numvb: 32` --- these are multiplied by `fdStep` internally
  (`cfg%numcb = NUM_CB_STATES * cfg%fdStep`) to give 64 CB and 606 VB states
  for the full $808 \times 808$ QW Hamiltonian.

### 3.3 Interpreting the Output

The program prints the $2\times 2$ spin matrices $\Sigma^{(x,y,z)}$ for
verification. For a pure CB doublet with dominant $s$-like character, these
should be close to the Pauli matrices. Deviations indicate VB/CB mixing.

The final g-factor values are printed and written to `output/gfactor.dat`:

```
 gx
 ... (eigenvalue calculation)
 gy
 ...
 gz
 ...
```

Each line reports $g_d = 2\lambda_{\min}^{(d)}$. For a symmetric QW grown
along $z$, expect $g_x \approx g_y$ (in-plane) and $g_z$ (out-of-plane)
to differ. The sign convention is: negative g-factor means the Zeeman
splitting has the same sign as the spin-orbit-induced correction.

### 3.4 Convergence Considerations

Several numerical parameters affect the accuracy of the g-factor:

- **FD grid density** (`FDstep`): The momentum matrix elements involve spatial
  integrals. Too few grid points can miss rapid oscillations in the envelope
  functions. Test convergence by doubling `FDstep`.

- **Number of intermediate states** (`numcb`, `numvb`): All subbands up to
  the specified counts participate in the Lowdin sum. If too few states are
  included, the high-energy contributions are truncated. For narrow-gap
  materials, the dominant contribution comes from the nearest VB states, so
  moderate numbers often suffice.

- **FD order** (`FDorder`): Higher-order stencils (4, 6, 8, 10) give more
  accurate gradient operators, improving the momentum matrix elements.
  Second-order ($FDorder=2$) is the default and is adequate for most
  applications.

- **Energy denominator warnings**: The code prints a warning when
  $|\mathrm{denom}| < 10^{-7}\;\mathrm{eV}$ and skips that contribution.
  Occasional skips are normal (degenerate states), but many warnings suggest
  insufficient state separation or numerical issues.

---

## 4. Code Architecture

### 4.1 Module Structure

The g-factor code lives in two files:

| File | Role |
|---|---|
| `src/physics/gfactor_functions.f90` | `gfactorFunctions` module: spin matrices, $\sigma$ elements, $P_\mathrm{ele}$ calculation, Lowdin sum, tensor assembly. Contains both QW and wire variants. |
| `src/apps/main_gfactor.f90` | `gfactor` program: reads input, builds/diagonalizes Hamiltonian, calls g-factor routines, writes output. |

### 4.2 Key Subroutines

| Subroutine | Purpose |
|---|---|
| `init_spin_matrices()` | Initializes $8\times 8$ $\Sigma_x, \Sigma_y, \Sigma_z$ on first call (cached). |
| `sigmaElem(state1, state2, dir, ...)` | QW: computes $\langle \psi_1 | \Sigma_\mathrm{dir} | \psi_2 \rangle$ with Simpson integration. |
| `sigmaElem_2d(state1, state2, dir, grid)` | Wire: same as `sigmaElem` but with uniform $dx \cdot dy$ summation over 2D grid. |
| `set_perturbation_direction(d, smallk)` | Sets a unit wave vector along direction $d \in \{1,2,3\}$. |
| `pMatrixEleCalc(Pele, d, ...)` | QW/bulk: computes $\langle a | H^{(d)}_\mathrm{pert} | b \rangle$ using dense or sparse Hamiltonian. |
| `pMatrixEleCalc_2d(Pele, d, ...)` | Wire: builds per-direction CSR perturbation with `g1`/`g2`/`g3` flags, uses `csr_spmv`. |
| `gfactorCalculation(tensor, ...)` | QW/bulk: full Lowdin partitioning, inter-band + intra-band sums, tensor assembly. |
| `gfactorCalculation_wire(tensor, ...)` | Wire: same algorithm with 2D sparse operators. |
| `compute_optical_matrix_wire(...)` | Wire: optical transition matrix elements for all CB-VB pairs. |

### 4.3 Data Flow Diagram

```
input.cfg
    |
    v
read_and_setup() --> simulation_config, profile, kpterms
    |
    v
ZB8bandBulk / ZB8bandQW / ZB8bandGeneralized (k=0)
    |
    v
zheev / zheevd / FEAST --> eigenvalues, eigenvectors
    |
    v
Sort: cb_state(N, numcb), vb_state(N, numvb)
       cb_value(numcb),   vb_value(numvb)
    |
    v
For each direction d = 1,2,3:
    |
    +-- sigmaElem(cb(:,n), cb(:,m), dir) --> sigma(2,2,3)
    |
    +-- Build perturbation Hamiltonians:
    |       ZB8bandQW/Bulk(g='g') with mod1, mod2
    |       or ZB8bandGeneralized(g='g1'/'g2'/'g3') for wire
    |
    +-- Sum over intermediate states l:
    |       compute_pele(P1, mod1, n, l) * compute_pele(P2, mod2, l, m)
    |     - compute_pele(P3, mod2, n, l) * compute_pele(P4, mod1, l, m)
    |     / (En + Em - 2*El)
    |
    v
tensor(2,2,d) = -i * sum / hbar2O2m0 - (g_free/2) * sigma(2,2,d)
    |
    v
zheev(tensor(:,:,d)) --> g_d = 2 * lambda_min
    |
    v
output/gfactor.dat: gx, gy, gz
```

### 4.4 Memory Layout for Eigenstates

Understanding the memory layout is essential for following the code:

**QW** (dense, column-major): The eigenstate vector for band $n$ has components
$\psi_n(i + (j-1) \cdot N_\mathrm{fd})$ where $i = 1,\ldots,N_\mathrm{fd}$
is the spatial grid index and $j = 1,\ldots,8$ is the band index. This is the
**band-interleaved** layout: all spatial points for band 1 come first, then
all points for band 2, etc. The Hamiltonian is organized the same way.

**Wire** (flat 2D, column-major): The eigenstate has components
$\psi_n((j-1) \cdot N_\mathrm{grid} + f)$ where $f = (i_y - 1) n_x + i_x$
is the flat grid index and $j = 1,\ldots,8$ is the band index. The grid
traverses $x$ fastest (inner loop), $y$ slowest.

---

## 5. Summary

The g-factor calculation in this code implements the second-order Lowdin
partitioning of the 8-band k.p Hamiltonian to project the Zeeman interaction
onto a Kramers doublet of interest. The key physical ingredients are:

1. **Spin matrices** $\Sigma_d$ in the 8-band basis, encoding the
   spin-orbit-coupled angular momentum of each band.
2. **Momentum matrix elements** $\langle a | H^{(\alpha)}_\mathrm{pert} | b \rangle$,
   computed by building the linear-$k$ perturbation Hamiltonian with the
   `g`/`g1`/`g2`/`g3` flags.
3. **Inter-band coupling** (dominant for narrow-gap materials) through VB
   intermediate states, weighted by inverse energy denominators involving $E_g$.
4. **Intra-band coupling** through other subband states, important for
   strongly confined nanostructures.

The code supports all three geometries --- bulk (8-band), quantum well (1D
confinement, dense LAPACK), and quantum wire (2D confinement, sparse CSR with
FEAST) --- with the same underlying Lowdin formula. The wire mode additionally
provides optical transition matrix elements as a useful by-product.
