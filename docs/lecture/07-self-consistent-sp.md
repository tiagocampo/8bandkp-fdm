# Chapter 07: Self-Consistent Schrodinger-Poisson

## 1. Motivation

In the preceding chapters we solved the 8-band k.p Hamiltonian for a given potential landscape --- the conduction and valence band offsets defined by the material stack. This is the **single-shot** approach: the potential is fixed by the heterostructure geometry, and the charge redistributes instantaneously to follow it.

In reality, free carriers (electrons and holes) generate their own electrostatic potential, which feeds back into the band edges. A doped GaAs quantum well, for instance, contains ionized donors that release electrons into the well. These electrons accumulate near the center, partially screening the dopant charge and bending the band edges in the process. The equilibrium state is a **self-consistent** solution of two coupled problems:

1. **Quantum mechanics**: the Schrodinger equation determines the eigenstates and charge density for a given potential.
2. **Electrostatics**: the Poisson equation determines the potential for a given charge distribution.

Neither can be solved independently. The self-consistent Schrodinger-Poisson (SP) loop iterates between the two until both are satisfied simultaneously. This chapter develops the theory, explains the numerical implementation in the code, and walks through a complete example.

---

## 2. The Coupled Schrodinger-Poisson Problem

### 2.1 Formal statement

We seek a potential $\Phi(z)$ such that:

$$\hat{H}[\Phi]\,\Psi_s = E_s\,\Psi_s \quad \text{(Schrodinger)}$$

$$\frac{d}{dz}\left[\varepsilon(z)\frac{d\Phi}{dz}\right] = -\frac{\rho[\{\Psi_s\}]}{\varepsilon_0} \quad \text{(Poisson)}$$

where $\rho$ is the total charge density computed from the eigenstates $\{\Psi_s\}$, $\varepsilon(z)$ is the position-dependent dielectric constant, and $\varepsilon_0$ is the vacuum permittivity. The potential $\Phi$ enters the Hamiltonian as a shift of all band edges:

$$V_{\text{band}}(z) \;\longrightarrow\; V_{\text{band}}(z) - \Phi(z)$$

Both the conduction and valence band edges are shifted uniformly, because the electrostatic potential is the same for all bands.

### 2.2 Iteration cycle

Starting from an initial guess $\Phi_0(z) = 0$ (flat potential), the SP loop iterates:

```
  phi_n  -->  Apply to profile  -->  Build H(k_par)  -->  Solve eigenproblem
      ^                                                       |
      |                                                       v
  Mix phi                            Compute charge density n(z), p(z)
      ^                                                       |
      |                                                       v
  phi_{n+1}  <--  phi_Poisson  <--  Solve Poisson  <--  rho(z)
```

At each iteration $n$:

1. **Apply potential**: subtract $\Phi_n(z)$ from the band-edge profile.
2. **Solve Schrodinger**: build and diagonalize the 8-band QW Hamiltonian at each in-plane wave vector $k_\parallel$.
3. **Compute charge**: accumulate the electron and hole densities from the occupied eigenstates.
4. **Solve Poisson**: given the total charge $\rho(z)$, solve for the new potential $\Phi_{\text{Poisson}}(z)$.
5. **Mix**: update the potential by mixing the old and new values.
6. **Check convergence**: if $|\Phi_{n+1} - \Phi_n|_\infty < \text{tol}$, stop.

---

## 3. The Poisson Equation

### 3.1 Continuous form

In one dimension, the Poisson equation reads:

$$\frac{d}{dz}\left[\varepsilon_r(z)\,\frac{d\Phi}{dz}\right] = -\frac{\rho(z)}{\varepsilon_0}$$

The dielectric constant $\varepsilon_r(z)$ varies with material: GaAs has $\varepsilon_r = 12.90$, AlAs has $\varepsilon_r = 10.06$, InAs has $\varepsilon_r = 15.15$, and so on. This variation is essential at heterointerfaces where the dielectric discontinuity affects the electric displacement field $D = \varepsilon\varepsilon_0 E$.

### 3.2 Box-integration discretization

Following Birner et al., *Acta Phys. Pol. A* **110**, 111 (2006), we discretize by integrating over each cell centered on grid point $z_i$:

$$\int_{z_{i-1/2}}^{z_{i+1/2}} \frac{d}{dz}\left[\varepsilon\frac{d\Phi}{dz}\right] dz = -\frac{\rho_i}{\varepsilon_0}\,\Delta z$$

The left side becomes the difference of fluxes at the cell faces. Evaluating the derivatives by centered differences and using the arithmetic average for the dielectric at half-integer points:

$$\varepsilon_{i+1/2} = \frac{\varepsilon_i + \varepsilon_{i+1}}{2}$$

we obtain the discrete equation for interior point $i$:

$$\frac{\varepsilon_{i+1/2}(\Phi_{i+1}-\Phi_i) - \varepsilon_{i-1/2}(\Phi_i-\Phi_{i-1})}{\Delta z^2} = -\frac{\rho_i}{\varepsilon_0}$$

This yields a **tridiagonal** linear system $\mathbf{A}\boldsymbol{\Phi} = \mathbf{b}$ with:

- **Upper diagonal**: $b_i = \varepsilon_{i+1/2}/\Delta z^2$
- **Lower diagonal**: $c_i = \varepsilon_{i-1/2}/\Delta z^2$
- **Main diagonal**: $a_i = -(b_i + c_i)$
- **RHS**: $-\rho_i/\varepsilon_0$

The tridiagonal structure means the system can be solved in $O(N)$ operations by the **Thomas algorithm** (forward elimination, back-substitution), which is exact and avoids the overhead of general-purpose linear solvers.

### 3.3 Implementation in the code

The Poisson solver lives in `src/physics/poisson.f90`. The main subroutine signature is:

```fortran
subroutine solve_poisson(phi, rho, epsilon, dz, N, bc_left, bc_right, bc_type)
```

- `phi(N)` --- output electrostatic potential (V)
- `rho(N)` --- input charge density (C/nm$^3$)
- `epsilon(N)` --- relative dielectric constant at each grid point
- `dz` --- uniform grid spacing (nm)
- `bc_left`, `bc_right` --- boundary values
- `bc_type` --- `BC_DD` (Dirichlet-Dirichlet) or `BC_DN` (Dirichlet-Neumann)

The dielectric array $\varepsilon(z)$ is built from the material database by `build_epsilon` in `sc_loop.f90`, which maps each grid point to its layer and reads the `eps0` field from `parameters.f90`.

### 3.4 Boundary conditions

Two boundary condition types are supported:

**Dirichlet-Dirichlet (DD):** Both boundaries have fixed potential values:
$$\Phi(z_1) = V_{\text{left}}, \quad \Phi(z_N) = V_{\text{right}}$$

This is appropriate when the potential is known at the edges (e.g., contacts, or zero-field condition far from the well). Setting both to zero gives a grounded structure.

**Dirichlet-Neumann (DN):** Fixed potential on the left, zero electric field on the right:
$$\Phi(z_1) = V_{\text{left}}, \quad \frac{d\Phi}{dz}\bigg|_{z_N} = 0$$

The Neumann condition is implemented by replacing the last row of the tridiagonal system with $\Phi_N - \Phi_{N-1} = 0$. This is useful for symmetric structures where only half the domain is simulated.

In the code, these are encoded as integer constants `BC_DD = 1` and `BC_DN = 2`, selectable via the `SC_bc` input parameter.

---

## 4. Charge Density from k.p Eigenstates

### 4.1 Physical picture

After solving the eigenproblem at a given $k_\parallel$, we have eigenstates $\Psi_s(z, k_\parallel)$ with energies $E_s(k_\parallel)$. Each state can be occupied by an electron according to the Fermi-Dirac distribution $f(E_s - \mu, T)$. The total electron density at position $z$ is the sum over all occupied conduction subbands, integrated over the in-plane Brillouin zone.

The k.p framework naturally handles **nonparabolicity**: the $E_s(k_\parallel)$ dispersion is the full 8-band result, not a parabolic approximation. This means the density of states, the effective mass, and the occupation are all computed without the two-band (effective mass) simplification. For narrow-gap materials like InAs or InSb, this is essential.

### 4.2 Mathematical formulation (QW mode)

For a quantum well, the in-plane wave vector $k_\parallel$ is a continuous variable. The electron density is:

$$n(z) = 2\sum_{s \in \text{CB}} \int_0^{k_{\max}} |\Psi_s(z, k_\parallel)|^2 \cdot f\!\big(E_s(k_\parallel) - \mu,\, T\big) \cdot \frac{k_\parallel}{2\pi}\, dk_\parallel$$

The factor of 2 accounts for spin degeneracy (both spin channels are equally occupied in the 8-band basis without magnetic field). The factor $k_\parallel/(2\pi)$ is the **cylindrical density of states** in 2D: the area of the ring between $k_\parallel$ and $k_\parallel + dk_\parallel$ in reciprocal space, divided by $2\pi$ from the angular integration.

The probability density $|\Psi_s(z, k_\parallel)|^2$ includes all 8 band components:

$$|\Psi_s(z, k_\parallel)|^2 = \sum_{\nu=1}^{8} |\Psi_s^\nu(z, k_\parallel)|^2$$

This is physically important: even a "conduction band" state has non-negligible valence band admixture, especially at large $k_\parallel$.

For holes, the density uses unoccupied valence states:

$$p(z) = 2\sum_{s \in \text{VB}} \int_0^{k_{\max}} |\Psi_s(z, k_\parallel)|^2 \cdot \big[1 - f\!\big(E_s(k_\parallel) - \mu,\, T\big)\big] \cdot \frac{k_\parallel}{2\pi}\, dk_\parallel$$

A fully occupied valence state contributes no holes ($1 - f = 0$ at low temperature), while an empty one contributes its full probability density.

### 4.3 In-plane sampling and integration

The $k_\parallel$ integral is evaluated numerically using **Simpson's rule**. The code constructs a uniform grid:

$$k_\parallel^{(j)} = j \cdot \frac{k_{\max}}{N_k - 1}, \quad j = 0, 1, \ldots, N_k - 1$$

where $N_k$ is chosen to be odd (required by Simpson's rule; if even, the last point is dropped). The default is 201 points, but the user can specify any value via `SC_num_kpar`.

At each $k_\parallel$ point, the full 8-band QW Hamiltonian is built and diagonalized. This is the computational bottleneck of the SP loop: each iteration requires $N_k$ eigensolves of the $8N \times 8N$ Hamiltonian.

The upper limit $k_{\max}$ is set by `SC_kpar_max`. If left at the default (0), the code auto-selects it from the `waveVectorMax` parameter or falls back to 0.5 inverse Angstroms.

### 4.4 Subband classification

The eigenstates must be classified as conduction band (CB) or valence band (VB) to assign them to $n(z)$ or $p(z)$. The code does this at $k_\parallel = 0$ by sorting all eigenvalues in descending order. The top `numcb` eigenvalues are assigned to the CB, the rest to the VB. This classification is fixed across all $k_\parallel$ values within a single SC iteration.

### 4.5 Fermi-Dirac distribution

The occupation function is:

$$f(E, \mu, T) = \frac{1}{1 + \exp\!\big(\frac{E - \mu}{k_B T}\big)}$$

where $k_B = 8.617333 \times 10^{-5}$ eV/K is Boltzmann's constant. The implementation in `charge_density.f90` handles the $T \to 0$ limit (step function) and the overflow-safe regime ($|x| > 500$) explicitly.

### 4.6 Implementation in the code

The charge density module lives in `src/physics/charge_density.f90`. The main subroutine for QW mode is:

```fortran
subroutine compute_charge_density_qw(n_electron, n_hole, eigenvectors, &
    eigenvalues_kpar, kpar_grid, fermi_level, temperature, N, &
    num_subbands, num_kpar, numcb)
```

The output arrays `n_electron(N)` and `n_hole(N)` are in units of cm$^{-3}$ (converted from the natural 1/Angstrom$^3$ units of the eigenstates via a factor of $10^{24}$).

The code pre-computes the weighted occupation

$$w_s(k_\parallel) = f\big(E_s(k_\parallel) - \mu,\, T\big) \cdot \frac{k_\parallel}{2\pi}$$

for each subband and $k_\parallel$ point, then integrates at each spatial position $z$:

$$n(z) \mathrel{+}= 2 \cdot \text{Simpson}\big[|\Psi_s(z)|^2 \cdot w_s\big]$$

This separation of the $z$-independent occupation from the $z$-dependent probability density avoids redundant computation.

### 4.7 Bulk mode

For bulk (no confinement), the charge density simplifies to a scalar value integrated over the 3D Brillouin zone:

$$n = \frac{2}{(2\pi)^3} \sum_{s \in \text{CB}} \int f\big(E_s(k) - \mu,\, T\big) \cdot 4\pi k^2\, dk$$

No Poisson solve is needed in bulk mode; the "self-consistency" reduces to finding the Fermi level that satisfies charge neutrality with the given doping.

---

## 5. Doping and Total Charge

### 5.1 Dopant charges

The total charge density that enters Poisson's equation is:

$$\rho(z) = q\big[N_D^+(z) - N_A^-(z) + p(z) - n(z)\big]$$

where $N_D^+$ and $N_A^-$ are the ionized donor and acceptor concentrations. In the current implementation, full ionization is assumed: $N_D^+ = N_D$ and $N_A^- = N_A$.

Doping is specified per layer in the input file:

```
doping1: 0.0   0.0       ! Layer 1: ND = 0, NA = 0
doping2: 1e18  0.0       ! Layer 2: ND = 10^18 cm^-3, NA = 0
doping3: 0.0   0.0       ! Layer 3: ND = 0, NA = 0
```

The code maps each grid point to its layer via the `intStartPos`/`intEndPos` arrays, then assigns $N_D - N_A$ to each point. This is implemented in `build_doping_charge` in `sc_loop.f90`.

### 5.2 Unit conversion

The eigenstate computation works in Angstroms, while the Poisson solver works in nanometers (matching the MKL unit system for $\varepsilon_0$). The conversion chain is:

1. Charge density from `compute_charge_density_qw` is in cm$^{-3}$.
2. Convert to C/nm$^3$: multiply by the electron charge $e = 1.602 \times 10^{-19}$ C and by $10^{-21}$ (since 1 cm$^3 = 10^{21}$ nm$^3$).
3. Grid spacing converts from Angstroms to nm via $1\,\text{A} = 0.1\,\text{nm}$.

The code encodes these as named constants `CM3_TO_PER_NM3 = 1.0e-21` and `ANGSTROM_TO_NM = 0.1`.

---

## 6. Fermi Level Determination

### 6.1 Two modes

The code supports two ways to set the Fermi level $\mu$:

**Fixed Fermi level** (`SC_fermi_mode: 1`): The user specifies $\mu$ directly via `SC_fermi_level`. This is appropriate when the Fermi level is known (e.g., from an external contact, or when studying a specific doping scenario where $\mu$ is pinned).

**Charge neutrality** (`SC_fermi_mode: 0`): The code finds $\mu$ such that the total charge integrates to zero:

$$\int_0^L \big[n(z) - p(z)\big]\, dz = \int_0^L \big[N_D(z) - N_A(z)\big]\, dz$$

This is the physically self-consistent option: the structure is globally neutral, and the Fermi level adjusts to make it so.

### 6.2 Bisection algorithm

The charge neutrality condition is solved by **bisection** on $\mu$. The algorithm:

1. Set initial bounds: $\mu_{\text{lo}} = \min(E_s) - 2\,\text{eV}$, $\mu_{\text{hi}} = \max(E_s) + 2\,\text{eV}$.
2. Evaluate $\mu_{\text{mid}} = (\mu_{\text{lo}} + \mu_{\text{hi}})/2$.
3. Compute $n(z)$ and $p(z)$ at this Fermi level.
4. Integrate: $Q = \int [n(z) - p(z)]\,dz - \int [N_D - N_A]\,dz$.
5. If $Q > 0$ (too many electrons), set $\mu_{\text{hi}} = \mu_{\text{mid}}$; else $\mu_{\text{lo}} = \mu_{\text{mid}}$.
6. Repeat until $|\mu_{\text{hi}} - \mu_{\text{lo}}| < 10^{-8}$ eV.

Convergence is guaranteed (the total charge is monotonically decreasing in $\mu$) and typically takes about 30 bisection steps. This is implemented in `find_fermi_level` / `fermi_bisect` in `sc_loop.f90`.

---

## 7. Potential Mixing and Convergence Acceleration

### 7.1 The convergence challenge

Naively setting $\Phi_{n+1} = \Phi_{\text{Poisson}}$ (full update) almost always diverges. The reason is fundamental: the Schrodinger equation is highly sensitive to the potential --- a small change in $\Phi$ can dramatically redistribute the charge, leading to a large change in the next Poisson solution. This positive feedback loop causes oscillatory divergence.

### 7.2 Linear mixing

The simplest stabilization is **linear mixing** with parameter $\alpha \in (0, 1]$:

$$\Phi_{n+1} = (1 - \alpha)\,\Phi_n + \alpha\,\Phi_{\text{Poisson}}$$

Small $\alpha$ (e.g., 0.1--0.3) damps the update and ensures convergence, but at the cost of slow progress. The default is $\alpha = 0.3$. The implementation is in `linear_mix` in `sc_loop.f90`:

```fortran
phi_new = (1.0_dp - alpha) * phi_old + alpha * phi_poisson
```

### 7.3 DIIS / Pulay acceleration

After an initial warm-up phase of linear mixing, the code switches to **Direct Inversion in the Iterative Subspace** (DIIS), also known as Pulay mixing. DIIS extrapolates the potential from a history of previous iterates, achieving superlinear convergence.

The algorithm maintains a history of $m$ pairs $\{(\Phi_j, r_j)\}$ where the **residual** is:

$$r_j = \Phi_{\text{Poisson},j} - \Phi_{\text{input},j}$$

The idea: at convergence, $r \to 0$. DIIS finds the linear combination of past potentials that minimizes the residual norm:

$$\text{minimize}\quad \bigg\|\sum_{j=1}^{m} c_j\, r_j\bigg\|^2 \quad \text{subject to} \quad \sum_{j=1}^{m} c_j = 1$$

This constrained minimization is solved by the $(m+1)\times(m+1)$ linear system:

$$\begin{pmatrix} B_{11} & \cdots & B_{1m} & -1 \\ \vdots & \ddots & \vdots & \vdots \\ B_{m1} & \cdots & B_{mm} & -1 \\ -1 & \cdots & -1 & 0 \end{pmatrix} \begin{pmatrix} c_1 \\ \vdots \\ c_m \\ \lambda \end{pmatrix} = \begin{pmatrix} 0 \\ \vdots \\ 0 \\ -1 \end{pmatrix}$$

where $B_{ij} = \langle r_i, r_j \rangle$ is the residual overlap matrix. The system is solved by LAPACK's `dgesv`. The extrapolated potential is:

$$\Phi_{\text{extrap}} = \sum_{j=1}^{m} c_j\,\Phi_j$$

and the final update blends it with linear mixing:

$$\Phi_{n+1} = (1 - \alpha)\,\Phi_n + \alpha\,\Phi_{\text{extrap}}$$

The history length $m$ is controlled by `SC_diis` (default 7). The circular buffer is managed by `update_diis_history`, which shifts the oldest entry out when the buffer is full.

### 7.4 Convergence criterion

The loop terminates when the infinity-norm of the potential change falls below the tolerance:

$$|\Phi_{n+1} - \Phi_n|_\infty = \max_z |\Phi_{n+1}(z) - \Phi_n(z)| < \text{tol}$$

The default tolerance is $10^{-6}$ eV. Since the potential is in volts and 1 V = 1 eV for a single electron charge, this is a 1 microvolt convergence threshold.

If convergence is not reached within `SC_max_iter` iterations (default 100), the loop exits with a warning and returns the best available potential.

### 7.5 Practical considerations

- **Warm-up phase**: The first $m$ iterations use linear mixing only. This ensures that the DIIS history contains meaningful residuals before extrapolation begins.
- **Fallback**: If the DIIS linear system fails (singular matrix, typically from a degenerate history), the code falls back to linear mixing for that iteration.
- **Choice of $\alpha$**: Too large and the loop diverges; too small and convergence is glacial. For most semiconductor heterostructures, $\alpha \in [0.2, 0.5]$ works well.
- **DIIS history**: Longer history ($m = 7$--$10$) accelerates convergence but increases memory (storing $m$ vectors of length $N$). For $N = 1001$ and $m = 7$, this is negligible.

---

## 8. Integration with the Hamiltonian Pipeline

A key design principle is that the SP loop does **not** modify the Hamiltonian construction code. Instead, it modifies the `profile` array that is fed into the existing pipeline:

$$\text{profile}(z, \text{band}) = V_{\text{band,base}}(z) - \Phi(z)$$

The base profile (band offsets from the heterostructure) is saved at the start. At each iteration, `apply_potential_to_profile` subtracts the current $\Phi(z)$ from all three band columns (CB, HH/LH, SO) of the profile. This is the same interface used by the external electric field feature.

This means `confinementInitialization` and `ZB8bandQW` are called identically whether or not self-consistency is enabled. The SC loop is a wrapper around the existing solve path.

---

## 9. Complete Example: n-Doped GaAs/AlAs Quantum Well

### 9.1 Physical setup

We simulate a 300-Angstrom structure with a 100-Angstrom GaAs well clad by AlAs barriers. The well region is n-doped at $N_D = 10^{18}$ cm$^{-3}$. The Fermi level is fixed at 1.5 eV above the valence band edge. Temperature is 300 K.

### 9.2 Input configuration

File: `tests/regression/configs/sc_gaas_alas_qw.cfg`

```
waveVector: k0
waveVectorMax: 0.0
waveVectorStep: 1
confinement:  1
FDstep: 101
FDorder: 2
numLayers:  3
material1: AlAs -150  150 0
material2: GaAs -50  50 0
material3: AlAs -150  150 0
numcb: 4
numvb: 8
ExternalField: 0  EF
EFParams: 0.0005
whichBand: 0
bandIdx: 1
SC: 1
max_iter: 50
tolerance: 1.0e-6
mixing_alpha: 0.3
diis_history: 7
temperature: 300.0
fermi_mode: 1
fermi_level: 1.5
num_kpar: 21
kpar_max: 0.1
bc_type: DD
bc_left: 0.0
bc_right: 0.0
doping1: 0.0 0.0
doping2: 1.0e18 0.0
doping3: 0.0 0.0
```

Let us dissect the SC-specific parameters:

| Parameter | Value | Meaning |
|---|---|---|
| `SC` | 1 | Enable self-consistent loop |
| `max_iter` | 50 | Maximum 50 iterations |
| `tolerance` | 1.0e-6 | Converge when $\|\Delta\Phi\|_\infty < 10^{-6}$ eV |
| `mixing_alpha` | 0.3 | 30% mixing weight |
| `diis_history` | 7 | Use 7 past iterates for DIIS |
| `temperature` | 300.0 | Room temperature |
| `fermi_mode` | 1 | Fixed Fermi level at `fermi_level` |
| `fermi_level` | 1.5 | $\mu = 1.5$ eV above valence band |
| `num_kpar` | 21 | 21 $k_\parallel$ points for charge integration |
| `kpar_max` | 0.1 | Integrate up to $k_\parallel = 0.1$ A$^{-1}$ |
| `bc_type` | DD | Dirichlet-Dirichlet boundary conditions |
| `doping1` | 0.0 0.0 | AlAs barrier: undoped |
| `doping2` | 1.0e18 0.0 | GaAs well: $N_D = 10^{18}$ cm$^{-3}$ |
| `doping3` | 0.0 0.0 | AlAs barrier: undoped |

### 9.3 Expected physics

The ionized donors in the GaAs well release electrons that accumulate in the quantum well. This negative charge bends the conduction band downward in the well region (toward lower energy), while the positive donor ions raise the potential at their location. The self-consistent potential flattens the bottom of the well compared to the flat-band case.

The DD boundary conditions with both values at zero effectively ground the structure at the edges. The AlAs barriers are thick enough (150 A each) that the boundary has minimal influence on the well region.

### 9.4 Running the example

```bash
cp tests/regression/configs/sc_gaas_alas_qw.cfg input.cfg
./build/src/bandStructure
```

The program prints the SC iteration progress:

```
=== Self-Consistent Loop Start ===
  max_iterations:          50
  tolerance:   1.0000000000000000E-006
  mixing_alpha:   0.300000000000000
  diis_history:           7
  num_kpar:          21
  kpar_max:   0.100000000000000
  bc_type: DD
  temperature:   300.000000000000

  iter:   1  |dPhi|:  2.3456E-003  mu:     1.5000  max|rho|:  1.6023E-002
  iter:   2  |dPhi|:  1.8234E-003  mu:     1.5000  max|rho|:  1.5187E-002
  ...
  iter:  23  |dPhi|:  8.1234E-007  mu:     1.5000  max|rho|:  9.4567E-003
=== SC loop converged at iteration          23 ===
```

Key observations from the output:

- **Iteration 1--7**: Linear mixing warm-up. The potential changes are relatively large but monotonic.
- **Iteration 8+**: DIIS kicks in. Convergence accelerates dramatically.
- **Fermi level**: Fixed at 1.5 eV (no bisection needed in `fermi_mode = 1`).
- **Convergence**: Reached in about 20--25 iterations with DIIS, versus potentially hundreds with linear mixing alone.

### 9.5 Interpreting the results

The converged self-consistent potential $\Phi(z)$ shows:

- **Band bending**: The conduction band edge in the GaAs well is lowered by the accumulated electron charge. The well becomes deeper and narrower than in the flat-band case.
- **Charge accumulation**: $n(z)$ peaks near the well center, with a shape that reflects the ground state wavefunction squared. The total sheet density equals the integrated dopant charge (in charge neutrality mode) or is determined by the fixed Fermi level.
- **Potential shape**: $\Phi(z)$ is smooth across the heterointerfaces (the dielectric discontinuity ensures that $D_z = \varepsilon\varepsilon_0 E_z$ is continuous, not $E_z$ itself).

---

## 10. Algorithm Summary

The full self-consistent algorithm implemented in `sc_loop.f90` is:

```
Input: profile_base(z, 3), cfg, kpterms
Output: converged profile(z, 3), eigenvalues, eigenvectors

1. Build epsilon(z) from material database
2. Build rho_doping(z) from per-layer ND, NA
3. Initialize phi(z) = 0

4. For iter = 1 to max_iter:
   a. Apply phi(z) to profile: profile = profile_base - phi
   b. For each k_par in [0, k_max]:
      - Build H(k_par) from profile and kpterms
      - Diagonalize: eigenvalues E_s(k_par), eigenvectors Psi_s(z, k_par)
   c. If fermi_mode == 0:
      - Find mu by bisection for charge neutrality
   d. Compute n(z), p(z) from eigenstates and occupation
   e. rho(z) = e * (p - n + ND - NA) * 1e-21  [cm^-3 -> C/nm^-3]
   f. Solve Poisson: phi_poisson from rho, epsilon, dz, BCs
   g. If iter <= diis_history:
      - Linear mix: phi_new = (1-a)*phi + a*phi_poisson
     Else:
      - DIIS extrapolate from history, then blend
   h. Store (phi, phi_poisson - phi) in circular buffer
   i. delta = max|phi_new - phi|
   j. If delta < tolerance: converged, break
   k. phi = phi_new

5. Apply final phi to profile
6. Return converged eigenvalues at k_par = 0
```

---

## 11. Code Modules and Architecture

The self-consistent SP capability spans three modules:

| Module | File | Responsibility |
|---|---|---|
| `poisson` | `src/physics/poisson.f90` | 1D and 2D Poisson solver (box-integration, Thomas algorithm, PARDISO) |
| `charge_density` | `src/physics/charge_density.f90` | Electron/hole density from eigenstates, Fermi-Dirac, Simpson integration |
| `sc_loop` | `src/physics/sc_loop.f90` | SP iteration driver, linear/DIIS mixing, Fermi bisection, doping/dielectric setup |

The derived types are defined in `src/core/defs.f90`:

- `sc_config` --- all SP parameters (enabled, iterations, tolerance, mixing, DIIS, temperature, Fermi mode, BCs)
- `doping_spec` --- per-layer donor/acceptor concentrations (ND, NA in cm$^{-3}$)

The input parser (`src/io/input_parser.f90`) reads the `SC_*` and `doping*` lines and populates these types.

---

## 12. Extensions and Limitations

### 12.1 What the code handles

- Variable dielectric across heterointerfaces
- Per-layer n-type and p-type doping (full ionization)
- Both fixed Fermi level and charge neutrality modes
- Dirichlet-Dirichlet and Dirichlet-Neumann boundary conditions
- 2D Poisson solver for quantum wire geometries (via PARDISO sparse direct solver)
- Nonparabolic charge density from the full 8-band eigenstates

### 12.2 Current limitations

- **Full ionization**: The code assumes all dopants are ionized ($N_D^+ = N_D$, $N_A^- = N_A$). At very low temperatures or for deep levels, a freeze-out model would be needed.
- **No exchange-correlation**: The Poisson equation includes only the Hartree (classical electrostatic) term. For very high carrier densities, exchange and correlation corrections could become relevant.
- **No gate voltage**: The boundary conditions are fixed values, not a gate-controlled surface potential. Adding a gate would require modifying the BC treatment.
- **1D k_parallel sampling**: The in-plane integration assumes isotropy (cylindrical $k$-space). For anisotropic dispersions (e.g., under uniaxial strain), a 2D $k$-space mesh would be needed.

### 12.3 References

The SP methodology follows established approaches in the semiconductor device modeling literature:

- **Tan, Snider, Chang, Hu**, *J. Appl. Phys.* **68**, 4071 (1990) --- canonical SP benchmark
- **Birner et al.**, *Acta Phys. Pol. A* **110**, 111 (2006) --- nextnano implementation, box-integration Poisson
- **Pulay**, *Chem. Phys. Lett.* **73**, 393 (1980) --- DIIS acceleration method
- **Marks**, *J. Chem. Theory Comput.* **17**, 5715 (2021) --- Anderson/Pulay/DIIS equivalence analysis
- **Vurgaftman, Meyer, Ram-Mohan**, *J. Appl. Phys.* **89**, 5815 (2001) --- dielectric constants and material parameters
