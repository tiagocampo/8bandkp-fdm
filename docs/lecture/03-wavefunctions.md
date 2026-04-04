# Chapter 3: Eigenstates and Wavefunctions

## 3.1 From Eigenvalues to Eigenvectors

In Chapter 2 we constructed the Hamiltonian matrix $H(\mathbf{k})$ and diagonalized it to obtain the energy eigenvalues $E_n(\mathbf{k})$. The diagonalization routine (LAPACK's `zheevx` for dense matrices, or MKL sparse eigensolvers for large sparse problems) returns not only the eigenvalues but also the corresponding eigenvectors. These eigenvectors are the wavefunctions of the system, and they carry far richer physical information than the eigenvalues alone.

For a quantum well with $N$ finite-difference grid points, the Hamiltonian is an $8N \times 8N$ Hermitian matrix. Its eigenvectors are column vectors of dimension $8N$ with complex entries. Each eigenvector $\mathbf{c}^{(n)}$ describes one electronic state:

$$
H(\mathbf{k})\, \mathbf{c}^{(n)} = E_n(\mathbf{k})\, \mathbf{c}^{(n)}, \qquad n = 1, 2, \ldots, 8N.
$$

The question we address in this chapter is: **what do these $8N$ complex numbers mean physically, and how do we extract spatial and band-resolved information from them?**

## 3.2 The Block Structure of Eigenvectors

### 3.2.1 The 8-band basis

Recall that the 8-band Kane basis for zinc-blende semiconductors comprises:

| Band index | Label | Character |
|---|---|---|
| 1 | $|J=3/2,\, m_J=+3/2\rangle$ | Heavy hole (HH1) |
| 2 | $|J=3/2,\, m_J=+1/2\rangle$ | Light hole (LH1) |
| 3 | $|J=3/2,\, m_J=-1/2\rangle$ | Light hole (LH2) |
| 4 | $|J=3/2,\, m_J=-3/2\rangle$ | Heavy hole (HH2) |
| 5 | $|J=1/2,\, m_J=+1/2\rangle$ | Split-off (SO1) |
| 6 | $|J=1/2,\, m_J=-1/2\rangle$ | Split-off (SO2) |
| 7 | $|\Gamma_6,\, \uparrow\rangle$ | Conduction band (CB1) |
| 8 | $|\Gamma_6,\, \downarrow\rangle$ | Conduction band (CB2) |

This ordering is fixed throughout the code (bands 1--4 are valence, 5--6 are split-off, 7--8 are conduction) and must never be changed.

### 3.2.2 Spatial decomposition: the band-block structure

The $8N$-dimensional eigenvector is organized into 8 contiguous blocks of $N$ entries each. Entry $i$ of the $b$-th block corresponds to the amplitude of band $b$ at the $i$-th grid point $z_i$:

$$
\mathbf{c}^{(n)} = \begin{pmatrix} \psi_1^{(n)}(z_1) \\ \vdots \\ \psi_1^{(n)}(z_N) \\ \psi_2^{(n)}(z_1) \\ \vdots \\ \psi_2^{(n)}(z_N) \\ \vdots \\ \psi_8^{(n)}(z_1) \\ \vdots \\ \psi_8^{(n)}(z_N) \end{pmatrix}
$$

where $\psi_b^{(n)}(z_i)$ is the complex amplitude of band $b$ in eigenstate $n$ at grid point $z_i$. The flat index for band $b$ at spatial point $i$ is:

$$
\text{flat\_idx}(b, i) = (b - 1) \times N + i, \qquad b = 1, \ldots, 8, \quad i = 1, \ldots, N.
$$

This mapping is implemented directly in the code's `get_eigenvector_component` subroutine, which extracts the $N$-element spatial profile for a given band from the full eigenvector:

```fortran
base_idx = (band_index - 1) * fdstep + i
component(i) = abs(eigenvectors(base_idx, eigenstate_index))
```

The full multi-component wavefunction is thus a sum over bands:

$$
\Psi^{(n)}(\mathbf{r}) = \sum_{b=1}^{8} \psi_b^{(n)}(z)\, |b\rangle,
$$

where $|b\rangle$ denotes the basis kets listed in the table above, and $\psi_b^{(n)}(z)$ is obtained by interpolating the discrete values $\psi_b^{(n)}(z_i)$.

### 3.2.3 Bulk eigenvectors (no spatial dependence)

For bulk semiconductors ($N = 1$, no confinement), the Hamiltonian is only $8 \times 8$, and each eigenvector has exactly 8 complex components. There is no spatial profile -- the wavefunction is a plane wave extended throughout the crystal. The output code handles this case by writing only the 8 complex amplitudes:

$$
\mathbf{c}^{(n)}_{\text{bulk}} = \begin{pmatrix} c_1^{(n)} \\ c_2^{(n)} \\ \vdots \\ c_8^{(n)} \end{pmatrix}.
$$

The "parts" (band-resolved probability) are computed trivially as $|c_b^{(n)}|^2$.

## 3.3 Band-Resolved Probability Density

### 3.3.1 Per-band spatial density

For a quantum well eigenstate $n$, the probability density contributed by band $b$ at position $z$ is:

$$
\rho_b^{(n)}(z) = |\psi_b^{(n)}(z)|^2.
$$

The code computes and writes this quantity directly. In `writeEigenfunctions`, for each eigenstate and each of the 8 bands, the absolute value of the complex amplitude is extracted at every grid point:

```fortran
do m = 1, 8
  call get_eigenvector_component(A, j, m, fdstep, N, component)
  eigv_abs(:,m) = abs(component)
end do
```

The output file `eigenfunctions_k_XXXXX_ev_YYYYY.dat` contains $N$ rows (one per grid point), each with 9 columns: the $z$-coordinate followed by $|\psi_b^{(n)}(z_i)|$ for $b = 1, \ldots, 8$.

### 3.3.2 Total probability density

The total probability density at position $z$ is the sum over all bands:

$$
\rho^{(n)}(z) = \sum_{b=1}^{8} |\psi_b^{(n)}(z)|^2.
$$

Normalization of the eigenvector (as returned by LAPACK) means:

$$
\sum_{b=1}^{8} \int |\psi_b^{(n)}(z)|^2\, dz = \sum_{b=1}^{8} \sum_{i=1}^{N} |\psi_b^{(n)}(z_i)|^2 \approx 1.
$$

This is not the same as saying each band individually integrates to unity -- quite the contrary. The distribution across bands is precisely the information that tells us about the state's character.

## 3.4 Band Character: The Parts File

### 3.4.1 Definition

The **integrated band probability** (called "parts" in the code) quantifies how much of eigenstate $n$ resides in each of the 8 bands:

$$
P_b^{(n)} = \int |\psi_b^{(n)}(z)|^2\, dz \approx \Delta z \sum_{i=1}^{N} |\psi_b^{(n)}(z_i)|^2,
$$

where $\Delta z$ is the uniform grid spacing. The code implements this as a simple rectangular quadrature:

```fortran
do m = 1, 8
  call get_eigenvector_component(A, j, m, fdstep, N, component)
  eigv_abs(:,m) = abs(component)
  parts(j,m) = sum(eigv_abs(:,m)**2) * (z(2) - z(1))
end do
```

The output file `parts.dat` contains `evnum` rows (one per eigenstate) with 8 columns, giving $P_b^{(n)}$ for each band. By construction:

$$
\sum_{b=1}^{8} P_b^{(n)} = 1 \qquad \text{(normalization)}.
$$

### 3.4.2 Physical interpretation

The parts vector $(P_1^{(n)}, P_2^{(n)}, \ldots, P_8^{(n)})$ tells you the **character** of eigenstate $n$:

- **Conduction band state**: $P_7^{(n)} + P_8^{(n)} \approx 1$, with bands 1--6 having negligible weight.
- **Heavy-hole state**: $P_1^{(n)} + P_4^{(n)} \approx 1$ (the two HH bands with $m_J = \pm 3/2$).
- **Light-hole state**: $P_2^{(n)} + P_3^{(n)} \approx 1$ (the two LH bands with $m_J = \pm 1/2$).
- **Split-off state**: $P_5^{(n)} + P_6^{(n)} \approx 1$.

In practice, at finite in-plane wavevector $\mathbf{k}_\parallel = (k_x, k_y)$, the bands couple and no eigenstate is purely of one character. The parts vector provides a quantitative measure of this **band mixing**.

## 3.5 Heavy-Hole vs Light-Hole Mixing

### 3.5.1 Mixing at finite in-plane wavevector

At the zone center ($\mathbf{k}_\parallel = 0$), the quantum well eigenstates have definite angular momentum character. The HH bands ($b = 1, 4$) and LH bands ($b = 2, 3$) decouple in the diagonal blocks $Q$ and $T$ of the Hamiltonian. However, the off-diagonal terms $S$, $R$, and their conjugates couple HH and LH at finite $k_\parallel$.

Consider the kp-term $S$, which appears in the Hamiltonian block structure as:

$$
S = 2\sqrt{3}\, k_{-}\, \gamma_3(z)\, \frac{d}{dz}, \qquad k_{-} = k_x - ik_y.
$$

This term couples HH1 (band 1) to LH1 (band 2) and appears in the off-diagonal block $H_{1,2}$. Similarly, $R$ couples HH to LH through the $k_x^2 - k_y^2$ dependence:

$$
R = -\sqrt{3}\left[\gamma_2(k_x^2 - k_y^2) - 2i\gamma_3 k_x k_y\right].
$$

The consequence is that at any finite $k_\parallel$, a nominally "heavy-hole" eigenstate acquires light-hole admixture, and vice versa. The degree of mixing depends on:

1. **The magnitude of $k_\parallel$**: larger $k$ means stronger coupling.
2. **The well width**: narrow wells push HH and LH subbands further apart (quantum confinement energy scales as $1/L^2$), reducing mixing.
3. **The Luttinger parameters** $\gamma_2$ and $\gamma_3$: materials with larger $\gamma_2/\gamma_3$ anisotropy show stronger mixing.
4. **The direction of $\mathbf{k}_\parallel$**: the mixing is anisotropic because $R$ depends on the angle $\phi = \arctan(k_y/k_x)$ in the $k_x$-$k_y$ plane.

### 3.5.2 Quantifying the mixing

The parts vector provides the simplest quantification. Define the **HH fraction** and **LH fraction** of eigenstate $n$:

$$
f_{\text{HH}}^{(n)} = P_1^{(n)} + P_4^{(n)}, \qquad f_{\text{LH}}^{(n)} = P_2^{(n)} + P_3^{(n)}.
$$

At $k_\parallel = 0$, a ground-state heavy hole will have $f_{\text{HH}} \approx 1$ and $f_{\text{LH}} \approx 0$. At finite $k_\parallel$, $f_{\text{HH}}$ decreases and $f_{\text{LH}}$ increases. The crossover point where $f_{\text{HH}} = f_{\text{LH}}$ is an important feature in the valence band dispersion and is directly related to the change in effective mass and optical polarization selection rules.

The conduction band states (bands 7--8) also acquire valence band character at large $k$ through the interband coupling term $P$ (the Kane momentum matrix element). This non-parabolicity effect is encoded in the parts as small but nonzero $P_1^{(n)}, \ldots, P_6^{(n)}$ for nominally CB states.

## 3.6 Example: Quantum Well Eigenstates

### 3.6.1 Ground-state conduction band electron

Consider a GaAs/Al$_{0.3}$Ga$_{0.7}$As quantum well of width $L = 100$ A with $N = 200$ grid points. The ground-state conduction band electron ($n = 1$, the lowest state with energy above the CB edge) has a parts vector approximately:

$$
(P_1, \ldots, P_8) \approx (0.001,\; 0.001,\; 0.001,\; 0.001,\; 0.002,\; 0.002,\; 0.496,\; 0.496).
$$

The state is overwhelmingly conduction-band in character ($P_7 + P_8 \approx 0.992$), with nearly equal spin-up and spin-down components ($P_7 \approx P_8$), as expected for a state at $k_\parallel = 0$ where spin is degenerate. The tiny valence band admixture ($\sim 0.8\%$) comes from the interband coupling via the Kane parameter $P$.

The spatial profile $|\psi_7^{(1)}(z)|$ and $|\psi_8^{(1)}(z)|$ are identical (spin degeneracy) and show a single half-sine-wave-like envelope localized within the well. The probability density $\rho^{(1)}(z)$ peaks at the well center and decays evanescently into the barriers.

### 3.6.2 Heavy-hole ground state

The highest-energy valence state (HH1) at $k_\parallel = 0$ typically has:

$$
(P_1, \ldots, P_8) \approx (0.98,\; 0.00,\; 0.00,\; 0.00,\; 0.01,\; 0.01,\; 0.00,\; 0.00).
$$

The state is almost purely $|3/2, +3/2\rangle$ (band 1). The other HH band (band 4, $|3/2, -3/2\rangle$) carries essentially zero weight because at $k_\parallel = 0$ the two HH bands are decoupled. The small split-off admixture ($\sim 2\%$) arises from the $\Delta_{\text{SO}}$ coupling between the $J = 3/2$ and $J = 1/2$ manifolds, which is never truly zero.

At $k_\parallel = 0.05\;\text{A}^{-1}$ along $k_x$, the same state might show:

$$
(P_1, \ldots, P_8) \approx (0.82,\; 0.08,\; 0.03,\; 0.01,\; 0.03,\; 0.02,\; 0.01,\; 0.00).
$$

Now $f_{\text{HH}} = 0.83$ and $f_{\text{LH}} = 0.11$: the state has acquired 11% light-hole character. This mixing modifies the effective mass, the optical matrix elements, and the spin properties of the state.

### 3.6.3 Excited states and nodal structure

Higher eigenstates exhibit nodal structure in their spatial wavefunctions. The $m$-th excited subband in a given band manifold has $m$ nodes in the envelope function. For example, the second conduction subband (CB2) has one node at the well center, with $|\psi_7^{(n)}(z)|$ peaking near the well edges and passing through zero at $z = L/2$.

These nodes are visible directly in the eigenfunction output files. They are analogous to the nodes of a particle-in-a-box wavefunction, but distorted by the spatially varying material parameters (effective mass, band offsets) and the band-mixing terms in the k.p Hamiltonian.

## 3.7 Implementation: How the Code Writes Wavefunctions

### 3.7.1 The writeEigenfunctions subroutine

The output pipeline is:

1. **Loop over eigenstates** ($j = 1$ to `evnum`): For each state, open a file named `eigenfunctions_k_XXXXX_ev_YYYYY.dat`.

2. **Extract band components**: For each of the 8 bands, call `get_eigenvector_component` to obtain the $N$-element spatial profile using the flat-index mapping $\text{flat\_idx} = (b-1) \times N + i$.

3. **Write spatial data**: For bulk, write 8 absolute amplitudes. For QW, write 9 columns: $z$ followed by $|\psi_b(z_i)|$ for $b = 1, \ldots, 8$.

4. **Compute parts**: Integrate $|\psi_b(z)|^2$ over $z$ using the rectangle rule with spacing $\Delta z = z_2 - z_1$, and write all eigenstates to a single `parts.dat` file.

### 3.7.2 The 2D wire case

For quantum wires (confinement mode 2), the code uses `writeEigenfunctions2d`. The eigenvector layout is analogous but now the spatial grid is 2D:

$$
\text{flat\_idx}(b, i_x, i_y) = (b - 1) \times N_{\text{grid}} + (i_y - 1) \times n_x + i_x,
$$

where $N_{\text{grid}} = n_x \times n_y$ is the total number of spatial grid points. The probability density at each grid point is summed over all 8 bands:

$$
\rho^{(n)}(x, y) = \sum_{b=1}^{8} |c^{(n)}_{\text{flat\_idx}(b, i_x, i_y)}|^2.
$$

The output format is three columns ($x$, $y$, $\rho$) with blank lines separating $y$-rows, directly plottable with `gnuplot splot`. The band-resolved parts are also computed, integrating over the 2D area element $dA = \Delta x \times \Delta y$.

### 3.7.3 Spin degeneracy and Kramers theorem

For systems without magnetic fields or structural inversion asymmetry, every eigenstate has a Kramers partner: a state at the same energy with opposite spin. In the output, this manifests as pairs of states with nearly identical parts vectors (e.g., $P_7 \approx P_8$ for CB states, or $P_1 \approx P_4$ for HH states at $k_\parallel = 0$). At finite $k_\parallel$ in asymmetric structures (e.g., under an external electric field), this degeneracy can be lifted and the parts vectors of the two spin partners may differ.

## 3.8 Practical Guide: Reading and Plotting the Output

### 3.8.1 Eigenfunction files

Each file `eigenfunctions_k_XXXXX_ev_YYYYY.dat` has the format:

```
   z_1    |psi_1(z_1)|    |psi_2(z_1)|    ...    |psi_8(z_1)|
   z_2    |psi_1(z_2)|    |psi_2(z_2)|    ...    |psi_8(z_2)|
   ...
   z_N    |psi_1(z_N)|    |psi_2(z_N)|    ...    |psi_8(z_N)|
```

To plot the total probability density of eigenstate 3 at k-step 1:

```gnuplot
plot 'output/eigenfunctions_k_00001_ev_00003.dat' using 1:(($2**2)+($3**2)+($4**2)+($5**2)+($6**2)+($7**2)+($8**2)+($9**2)) with lines title '|psi|^2'
```

To plot the conduction band component (columns 8 and 9 correspond to bands 7 and 8):

```gnuplot
plot 'output/eigenfunctions_k_00001_ev_00003.dat' using 1:($8**2+$9**2) with lines title 'CB character'
```

### 3.8.2 Parts file

The file `parts.dat` has `evnum` rows and 8 columns:

```
P_1^(1)  P_2^(1)  P_3^(1)  P_4^(1)  P_5^(1)  P_6^(1)  P_7^(1)  P_8^(1)
P_1^(2)  P_2^(2)  P_3^(2)  P_4^(2)  P_5^(2)  P_6^(2)  P_7^(2)  P_8^(2)
...
```

A quick way to identify which eigenstates are conduction-band-like: look for rows where columns 7 and 8 sum to near unity. Similarly, valence band states have large entries in columns 1--4, and split-off states in columns 5--6.

### 3.8.3 Eigenvalue file

The file `eigenvalues.dat` has `waveVectorStep` rows. The first column is $|\mathbf{k}|$ and the remaining `evnum` columns are the sorted eigenvalues in eV. By cross-referencing with the parts file, you can identify the band character of each subband and track how it evolves with $k$.

## 3.9 Summary

The eigenvectors of the 8-band k.p Hamiltonian encode the full multi-component wavefunction of each electronic state. By decomposing the $8N$-component vector into 8 spatial profiles, we obtain:

- **Band-resolved probability densities** $|\psi_b(z)|^2$ showing how the wavefunction distributes across the HH, LH, SO, and CB manifolds.
- **Integrated parts** $P_b^{(n)}$ quantifying the band character of each eigenstate.
- **Heavy-hole/light-hole mixing** that evolves with in-plane wavevector and determines the optical and spin properties of the system.

These quantities are computed directly by the code's `writeEigenfunctions` subroutine using the flat-index mapping $\text{flat\_idx} = (b-1) \times N + i$, and are written to per-eigenstate data files and the aggregated `parts.dat` file. Understanding and visualizing these outputs is essential for interpreting the physics of quantum well and quantum wire structures simulated by the code.
