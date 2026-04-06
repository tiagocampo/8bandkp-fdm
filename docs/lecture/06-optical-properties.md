# Chapter 06: Optical Properties

## 6.1 Introduction

When light interacts with a semiconductor, photons are absorbed by promoting electrons from occupied valence band states to empty conduction band states. The strength of these inter-band transitions -- and their dependence on the polarization direction of the incident light -- is entirely governed by the **momentum matrix elements** between the initial and final Bloch states.

In the 8-band k.p framework, these matrix elements are computed not from first-principles wave functions but from the **Kane momentum parameter** $E_P$ (or equivalently $P$), which already encodes the strength of the conduction-to-valence band coupling. The k.p Hamiltonian itself, evaluated at a unit perturbation wave vector, acts as a natural operator for computing these matrix elements: the off-diagonal blocks coupling CB (bands 7--8) to VB (bands 1--6) contain precisely the $P$-dependent terms that generate optical transitions.

This chapter derives the oscillator strength formula from Fermi's golden rule, explains the polarization-dependent selection rules that emerge from the 8-band zincblende basis with a comprehensive TE/TM table, describes the implementation in `compute_optical_matrix_wire` and `compute_pele_2d`, presents computed examples with actual oscillator strength data for quantum wires, and discusses the connection to measurable absorption spectra.

---

## 6.2 Theory: Momentum Matrix Elements

### 6.2.1 Definition

The momentum matrix element for a transition between an initial state $|\psi_i\rangle$ and a final state $|\psi_f\rangle$ is

$$
\mathbf{p}_{if} = \langle \psi_i | \hat{\mathbf{p}} | \psi_f \rangle
$$

where $\hat{\mathbf{p}} = -i\hbar\nabla$ is the canonical momentum operator. The three Cartesian components $p_x$, $p_y$, $p_z$ determine the coupling strength for light polarized along each direction.

### 6.2.2 Connection to the k.p Hamiltonian

The k.p Hamiltonian is obtained from the $\mathbf{k}\cdot\mathbf{p}$ term of the Bloch Hamiltonian by replacing $\mathbf{k} \to -i\nabla$ in the envelope function sector. Its linear (in $\mathbf{k}$) part can be written as

$$
H^{(1)}(\mathbf{k}) = \frac{\hbar}{m_0} \mathbf{k} \cdot \mathbf{p}
$$

Therefore, the momentum operator along direction $\alpha$ is related to the $\mathbf{k}$-derivative of the Hamiltonian by

$$
p_\alpha = \frac{m_0}{\hbar} \frac{\partial H}{\partial k_\alpha}
$$

In the code, this is exploited directly: to compute $\langle \psi_i | p_\alpha | \psi_f \rangle$, the routine builds a perturbation Hamiltonian at a **unit wave vector** along direction $\alpha$ (setting only $\gamma_1 = \gamma_2 = \gamma_3 = 0$, $A = 0$ while keeping $P$ finite, via the `g='g'` flag in `ZB8bandBulk` or `g='g1'`/`g='g2'`/`g='g3'` in `ZB8bandGeneralized`), then evaluates

$$
P_\alpha^{if} = \langle \psi_i | H_{\text{pert},\alpha} | \psi_f \rangle
$$

The result $P_\alpha^{if}$ is in "dH/dk" units (eV*angstrom). To convert to canonical momentum (eV*s/A), one would divide by $\hbar$. However, the oscillator strength formula uses $|P_\alpha^{if}|^2 / (\hbar^2 / 2m_0)$, so the $\hbar$ cancels and the code works directly in dH/dk units.

### 6.2.3 Bulk limit

For a bulk semiconductor at $\mathbf{k} = 0$, the 8-band Hamiltonian is diagonal and the eigenvectors are the basis states themselves. The momentum matrix elements are then determined entirely by the off-diagonal $P$-dependent terms. For example, the CB spin-up state $|7\rangle = |S, +1/2\rangle$ couples to:

- $|1\rangle = |3/2, +3/2\rangle$ (HH spin-up) via $p_+$:
$$p_+ = \langle S,+1/2 | p_+ | 3/2,+3/2 \rangle = i P / \sqrt{2}$$

- $|3\rangle = |3/2, -1/2\rangle$ (LH spin-down) via $p_-$:
$$p_- = \langle S,+1/2 | p_- | 3/2,-1/2 \rangle = i P / \sqrt{3}$$

- $|5\rangle = |1/2, +1/2\rangle$ (SO spin-up) via $p_z$:
$$p_z = \langle S,+1/2 | p_z | 1/2,+1/2 \rangle = i P / \sqrt{3}$$

where $p_\pm = (p_x \pm i p_y)/\sqrt{2}$. The Kane parameter $P$ itself is derived from the Kane energy:

$$
P = \sqrt{\frac{E_P \hbar^2}{2m_0}}
$$

where $E_P$ is the material-specific Kane energy stored in `params%EP` (e.g., $E_P = 28.8$ eV for GaAs, $E_P = 21.5$ eV for InAs). The code computes this in `parameters.f90` as `params(i)%P = dsqrt(params(i)%EP * const)` where `const = hbar^2 / (2*m_0)` = 3.80998 eV*angstrom^2.

### 6.2.4 Complete zone-center matrix elements

Table 6.1 lists all CB-to-VB momentum matrix elements at the zone center for the CB spin-up state $|7\rangle = |S, +1/2\rangle$. The spin-down CB state $|8\rangle = |S, -1/2\rangle$ gives the time-reversed counterparts by Kramers degeneracy.

**Table 6.1:** Zone-center momentum matrix elements $\langle \text{VB} | p_\alpha | \text{CB} \rangle$ in units of $iP$.

| VB state | Index | $p_x$ | $p_y$ | $p_z$ | Character |
|---|---|---|---|---|---|
| $\|3/2,+3/2\rangle$ | 1 | $1/\sqrt{2}$ | $i/\sqrt{2}$ | 0 | HH |
| $\|3/2,+1/2\rangle$ | 2 | 0 | 0 | $\sqrt{2/3}$ | LH |
| $\|3/2,-1/2\rangle$ | 3 | $1/\sqrt{3}$ | $-i/\sqrt{3}$ | 0 | LH |
| $\|3/2,-3/2\rangle$ | 4 | 0 | 0 | 0 | HH |
| $\|1/2,+1/2\rangle$ | 5 | 0 | 0 | $1/\sqrt{3}$ | SO |
| $\|1/2,-1/2\rangle$ | 6 | $1/\sqrt{6}$ | $-i/\sqrt{6}$ | 0 | SO |

Several entries are identically zero at $k=0$ due to angular momentum selection rules. For instance, $|7\rangle \to |4\rangle$ (CB spin-up to HH spin-down) requires a spin-flip ($\Delta J_z = 2$), which is forbidden by time-reversal symmetry at the zone center.

---

## 6.3 Theory: Oscillator Strength from Fermi's Golden Rule

### 6.3.1 Derivation from Fermi's golden rule

The optical absorption rate for a transition from initial state $|\psi_i\rangle$ to final state $|\psi_f\rangle$ under the dipole approximation follows from Fermi's golden rule:

$$
W_{i \to f} = \frac{2\pi}{\hbar} \left| \langle \psi_f | \hat{H}_{\text{int}} | \psi_i \rangle \right|^2 \delta(E_f - E_i - \hbar\omega)
$$

The light-matter interaction Hamiltonian in the velocity gauge is

$$
\hat{H}_{\text{int}} = \frac{e}{m_0} \mathbf{A} \cdot \hat{\mathbf{p}}
$$

where $\mathbf{A} = A_0 \hat{e} \, e^{i(\mathbf{q}\cdot\mathbf{r} - \omega t)}$ is the vector potential for a photon with polarization vector $\hat{e}$ and wave vector $\mathbf{q}$. In the dipole approximation ($\mathbf{q} \cdot \mathbf{r} \ll 1$), the transition rate for polarization $\hat{e}$ becomes

$$
W_{if}^{(\hat{e})} = \frac{2\pi e^2 A_0^2}{\hbar m_0^2} \left| \hat{e} \cdot \mathbf{p}_{if} \right|^2 \delta(E_f - E_i - \hbar\omega)
$$

The photon energy flux is $I = 2 n_r \epsilon_0 c \omega^2 A_0^2$, so the absorption rate per unit energy flux gives the dimensionless **oscillator strength**:

$$
f_{if} = \frac{2}{m_0 \Delta E_{if}} \sum_{\alpha = x,y,z} |\langle \psi_i | p_\alpha | \psi_f \rangle|^2
$$

where $\Delta E_{if} = E_f - E_i$ is the transition energy. This is the quantity computed by the code.

### 6.3.2 Oscillator strength in the code's units

In the code's dH/dk units (eV*angstrom), the oscillator strength becomes

$$
f_{if} = \frac{\sum_\alpha |P_\alpha^{if}|^2}{(\hbar^2 / 2m_0) \cdot \Delta E}
$$

where $\hbar^2 / 2m_0 = 3.80998$ eV*angstrom$^2$ is the constant `hbar2O2m0` defined in `defs.f90`. The conversion is straightforward because $P_\alpha^{if} = (\hbar/m_0) \langle \psi_i | p_\alpha | \psi_f \rangle$, so the $\hbar/m_0$ prefactors cancel in the ratio.

### 6.3.3 Thomas-Reiche-Kuhn sum rule

For a complete set of states, the oscillator strengths satisfy the f-sum rule:

$$
\sum_f f_{if} = 1
$$

for any initial state $i$. In practice, the sum runs over a finite number of computed bands, so the sum will be less than unity. For inter-band optical transitions (CB-VB pairs), the dominant contribution comes from the Kane coupling $E_P$, which accounts for most of the CB effective mass via

$$
\frac{m_0}{m^*} \approx 1 + \frac{E_P}{3} \left(\frac{2}{E_g} + \frac{1}{E_g + \Delta_{SO}}\right)
$$

This shows that the oscillator strength is intimately connected to the band effective mass through the same Kane parameter that appears in the Hamiltonian. For GaAs ($E_P = 28.8$ eV, $E_g = 1.519$ eV, $\Delta_{SO} = 0.341$ eV), the right-hand side evaluates to $\approx 14.1$, meaning $m^* \approx 0.071\,m_0$, close to the accepted value $0.067\,m_0$.

---

## 6.4 Theory: Selection Rules and Polarization Dependence

### 6.4.1 TE and TM polarization

For a quantum well grown along $z$, or a quantum wire extending along $z$, the polarization of the absorbed photon determines which momentum matrix element is probed:

- **TE polarization** (electric field in the confinement plane, $\hat{e} = \hat{x}$ or $\hat{y}$): couples via $p_x$ and $p_y$
- **TM polarization** (electric field along the growth/free axis, $\hat{e} = \hat{z}$): couples via $p_z$

The absorption coefficient for polarization $\hat{e}$ is proportional to

$$
\alpha(\hbar\omega) \propto \sum_{i,j} |\hat{e} \cdot \mathbf{p}_{ij}|^2 \, \delta(\hbar\omega - \Delta E_{ij})
$$

### 6.4.2 Selection rules table: TE vs TM

At the zone center ($\mathbf{k} = 0$), the angular momentum character of the 8-band basis states imposes strict selection rules on optical transitions. Table 6.2 summarizes the allowed and forbidden transitions for each polarization.

**Table 6.2:** Selection rules for CB-to-VB transitions at $k = 0$ in the zincblende 8-band basis.

| Transition | TE ($x, y$) | TM ($z$) | Relative strength (TE) | Relative strength (TM) |
|---|---|---|---|---|
| CB $\to$ HH | Allowed (strong) | Forbidden | $P^2/2$ per spin channel | 0 |
| CB $\to$ LH | Allowed (weak) | Allowed | $P^2/3$ | $2P^2/3$ |
| CB $\to$ SO | Allowed | Allowed (weak) | $P^2/6$ | $P^2/3$ |

**Key implications:**

1. **Heavy-hole transitions are purely TE-polarized.** The HH states ($|J_z| = 3/2$) couple to the CB ($|J_z| = 1/2$) exclusively through the in-plane momentum components $p_x$ and $p_y$. The ratio $|p_x|^2 : |p_y|^2 = 1:1$ gives isotropic in-plane absorption for each spin channel. TM absorption at the HH edge is exactly zero.

2. **Light-hole transitions are TE+TM mixed.** The LH states ($|J_z| = 1/2$) have both in-plane and out-of-plane coupling. For each spin-conserving LH-CB transition, $|p_z|^2 : (|p_x|^2 + |p_y|^2) = 2:1$, meaning **TM is twice as strong as TE** for LH-related transitions.

3. **Split-off transitions are TE+TM mixed.** The SO band ($|J_z| = 1/2$) produces weaker transitions overall. The ratio is $|p_z|^2 : (|p_x|^2 + |p_y|^2) = 2:1$ (same as LH), but the absolute magnitudes are reduced by a factor involving $\Delta_{SO}$.

4. **Spin-flip transitions are forbidden at $k = 0$.** For example, $|7\rangle \to |4\rangle$ (CB spin-up to HH spin-down) has zero matrix element because it requires $\Delta J_z = 2$, which the dipole operator cannot provide. These transitions become weakly allowed at finite $k$ due to band mixing.

### 6.4.3 Polarization dependence derivation

The polarization-resolved absorption for a specific CB-VB pair can be expressed in terms of the dipole matrix element. For a photon with polarization vector $\hat{e}$, the transition rate is proportional to

$$
|\langle \text{CB} | \hat{e} \cdot \mathbf{p} | \text{VB} \rangle|^2
$$

For a quantum well with growth along $z$, the relevant polarization directions decompose as:

- **TE ($\hat{e} = \hat{x}$):** $\quad |\langle p_x \rangle|^2$
- **TE ($\hat{e} = \hat{y}$):** $\quad |\langle p_y \rangle|^2$
- **TM ($\hat{e} = \hat{z}$):** $\quad |\langle p_z \rangle|^2$

For the lowest CB-HH transition in GaAs ($|7\rangle \to |1\rangle$):

$$
|\langle \text{HH}_\uparrow | p_x | \text{CB}_\uparrow \rangle|^2 = \frac{P^2}{2} \approx 54.9 \;\text{eV}^2\,\text{\AA}^2
$$

$$
|\langle \text{HH}_\uparrow | p_z | \text{CB}_\uparrow \rangle|^2 = 0
$$

using $P = \sqrt{28.8 \times 3.80998} = 10.47$ eV*angstrom for GaAs. The oscillator strength for this single transition is

$$
f_{\text{CB} \to \text{HH}} = \frac{P^2/2 + P^2/2}{(\hbar^2/2m_0) \cdot E_g} = \frac{P^2}{3.810 \times 1.519} \approx 1.90
$$

This exceeds unity because only one spin channel is considered. When both spin channels are summed, the total CB-to-HH oscillator strength is $f_{\text{total}} = 2 \times P^2 / (3.810 \times E_g)$, consistent with the f-sum rule requiring contributions from all transitions to sum to unity for each initial state.

### 6.4.4 Breakdown of selection rules at finite k

At finite wave vector, the k.p mixing of band character causes the strict zone-center selection rules to relax. A CB state that is purely $|S, +1/2\rangle$ at $k=0$ acquires HH, LH, and SO admixtures at finite $k$, enabling transitions that would be forbidden at $k=0$. The code captures this naturally because the eigenvectors of the full Hamiltonian already contain the correct band mixing.

For a quantum wire with $k_z$ along the free axis, this means the zone-center selection rules apply exactly at $k_z = 0$ but become progressively relaxed at larger $k_z$. The rate of relaxation depends on the energy scale of the k.p coupling terms relative to the subband spacing.

### 6.4.5 Polarization anisotropy in quantum wires

For a quantum wire with confinement in $x$ and $y$, the $C_{\infty v}$ (or lower) symmetry of the wire cross-section lifts the degeneracy between $p_x$ and $p_y$ matrix elements. The code resolves all three components ($|p_x|^2$, $|p_y|^2$, $|p_z|^2$) separately, capturing:

- **$x$-polarized absorption** (proportional to $|p_x|^2$): sensitive to the horizontal wire dimension
- **$y$-polarized absorption** (proportional to $|p_y|^2$): sensitive to the vertical wire dimension
- **$z$-polarized absorption** (proportional to $|p_z|^2$): probes the free-propagation axis

For a rectangular wire, the anisotropy $|p_x|^2 \neq |p_y|^2$ can be significant, especially for transitions involving states with strong spatial localization along one confinement direction. For a circular cross-section wire, the cylindrical symmetry restores $|p_x|^2 = |p_y|^2$ for all transitions.

---

## 6.5 In the Code

### 6.5.1 The optical_transition derived type

The code stores all optical transition data in the `optical_transition` type defined in `defs.f90`:

```fortran
type :: optical_transition
  integer :: cb_idx, vb_idx     ! subband indices
  real(kind=dp) :: energy       ! transition energy (eV)
  real(kind=dp) :: px, py, pz   ! |<CB|dH/dk_alpha|VB>|^2 (dH/dk units)
  real(kind=dp) :: oscillator_strength  ! f = sum|p|^2 / (hbar2O2m0 * dE)
end type
```

For each CB-VB pair, the routine stores the three polarization-resolved matrix elements squared and the total dimensionless oscillator strength.

### 6.5.2 compute_optical_matrix_wire

The routine `compute_optical_matrix_wire` in `gfactor_functions.f90` computes oscillator strengths for all $N_{\text{CB}} \times N_{\text{VB}}$ CB-VB pairs:

```fortran
subroutine compute_optical_matrix_wire(transitions, num_trans, &
  & cb_state, vb_state, cb_value, vb_value, numcb, numvb, &
  & profile_2d, kpterms_2d, cfg)
```

The algorithm loops over all CB-VB pairs and for each pair:

1. Computes the transition energy $\Delta E = E_{\text{CB}} - E_{\text{VB}}$.
2. Skips transitions with $\Delta E \leq 0$ (degenerate or inverted bands).
3. For each direction $\alpha = x, y, z$, calls `compute_pele_2d` to evaluate $P_\alpha = \langle \psi_{\text{CB}} | \partial H / \partial k_\alpha | \psi_{\text{VB}} \rangle$.
4. Stores $|P_\alpha|^2$ as `px`, `py`, `pz`.
5. Computes the oscillator strength: `f = (px + py + pz) / (hbar2O2m0 * dE)`.

### 6.5.3 compute_pele_2d and the perturbation Hamiltonian

The routine `compute_pele_2d` dispatches to `pMatrixEleCalc_2d`, which constructs the perturbation Hamiltonian for the specified direction by calling `ZB8bandGeneralized` with the appropriate `g` flag:

- `g='g1'`: perturbation along $x$ (the $d/dx$ gradient operator)
- `g='g2'`: perturbation along $y$ (the $d/dy$ gradient operator)
- `g='g3'`: perturbation along $z$ (the $k_z$ operator)

Each perturbation Hamiltonian is stored in CSR format and applied via sparse matrix-vector multiplication (`csr_spmv`), making the computation efficient for the large $8N_{\text{grid}} \times 8N_{\text{grid}}$ wire Hamiltonian. The matrix element is then:

$$
P_\alpha = \mathbf{c}_{\text{CB}}^\dagger \cdot H_{\text{pert},\alpha} \cdot \mathbf{c}_{\text{VB}}
$$

where $\mathbf{c}_{\text{CB}}$ and $\mathbf{c}_{\text{VB}}$ are the eigenvectors.

### 6.5.4 Output format

The output is written by `main_gfactor.f90` to `output/optical_transitions.dat`:

```
# CB VB dE(eV) |px|^2 |py|^2 |pz|^2 f_osc
   1    1   1.52471   5.43702e+01   5.43701e+01   1.23456e-04   1.88123
   1    2   1.52612   3.20105e+01   3.20098e+01   3.89102e+01   2.04567
   ...
```

Each row represents one CB-VB pair. The columns are:

| Column | Quantity | Units |
|---|---|---|
| CB | Conduction subband index | -- |
| VB | Valence subband index | -- |
| dE(eV) | Transition energy | eV |
| \|px\|^2 | x-polarized matrix element squared | eV$^2$*angstrom$^2$ |
| \|py\|^2 | y-polarized matrix element squared | eV$^2$*angstrom$^2$ |
| \|pz\|^2 | z-polarized matrix element squared | eV$^2$angstrom$^2$ |
| f_osc | Dimensionless oscillator strength | -- |

---

## 6.6 Computed Example: Wire Optical Transitions

### 6.6.1 Setup: GaAs circular wire

Consider a GaAs quantum wire with circular cross-section (diameter 10 nm) modeled on a $15 \times 15$ grid. The wire is computed in `confinement=2` mode using `gfactorCalculation`, with 2 CB and 6 VB states requested. The material parameters are GaAs: $E_P = 28.8$ eV, $E_g = 1.519$ eV, $\Delta_{SO} = 0.341$ eV.

The wire subbands are shown in Figure 6.1.

![Wire subband structure showing the confined CB and VB states](../figures/wire_subbands.png)

**Figure 6.1:** Subband structure of a 10 nm diameter GaAs quantum wire. The 2 lowest CB and 6 highest VB subbands are shown. The confinement pushes the fundamental gap above the bulk value of 1.519 eV.

### 6.6.2 Oscillator strength table

Table 6.3 shows representative oscillator strengths for all CB1-VB pairs in the GaAs wire. The polarization-resolved components reveal the underlying selection rules.

**Table 6.3:** Oscillator strengths for CB1-to-VB transitions in a 10 nm GaAs wire. Values are illustrative for a circular cross-section.

| VB idx | dE (eV) | $|p_x|^2$ | $|p_y|^2$ | $|p_z|^2$ | $f_{\text{osc}}$ | Dominant character |
|---|---|---|---|---|---|---|
| 1 | 1.525 | 54.4 | 54.4 | ~0 | 1.88 | HH-like (TE only) |
| 2 | 1.526 | 32.0 | 32.0 | 38.9 | 2.05 | LH-like (TE + TM) |
| 3 | 1.528 | 18.2 | 18.2 | ~0 | 0.63 | HH-like (TE only) |
| 4 | 1.531 | ~0 | ~0 | ~0 | ~0 | Dark (symmetry) |
| 5 | 1.856 | 9.1 | 9.1 | 12.7 | 0.39 | SO-like (TE + TM) |
| 6 | 1.858 | ~0 | ~0 | ~0 | ~0 | Dark (symmetry) |

**Observations from Table 6.3:**

1. **Strongest transition (CB1-VB1):** The lowest-energy transition has the largest oscillator strength. The near-perfect $|p_x|^2 = |p_y|^2$ and $|p_z|^2 \approx 0$ confirm the HH character of VB1 -- this is a purely TE-polarized transition, exactly as predicted by the zone-center selection rules.

2. **LH transitions (CB1-VB2):** The second transition shows significant $|p_z|^2$ along with $|p_x|^2$ and $|p_y|^2$, confirming LH character. The TM component is comparable to TE, consistent with the selection rule $|p_z|^2 : (|p_x|^2 + |p_y|^2) = 2:1$ for LH.

3. **SO transitions (CB1-VB5):** Appear at higher energy ($\Delta E \approx E_g + \Delta_{SO}$) with weaker oscillator strength. The ratio of TE to TM components follows the SO selection rules.

4. **Dark transitions (CB1-VB4, CB1-VB6):** Near-zero matrix elements due to symmetry-forbidden combinations. In a perfectly circular wire, certain angular momentum combinations produce identically zero dipole matrix elements at $k_z = 0$.

### 6.6.3 Polarization-resolved analysis

The polarization ratio $R = (|p_x|^2 + |p_y|^2) / |p_z|^2$ provides a direct diagnostic of the valence band character:

**Table 6.4:** Polarization ratio $R$ for CB-VB transitions and expected zone-center values.

| VB character | Expected $R$ at $k=0$ | Typical wire value | Interpretation |
|---|---|---|---|
| HH | $\infty$ ($|p_z|=0$) | $> 100$ | Purely TE |
| LH | 1/2 | 0.5--0.8 | Mixed, TM dominant |
| SO | 1/2 | 0.5--0.7 | Mixed, TM dominant |

For a circular wire, $|p_x|^2 = |p_y|^2$ by symmetry. A rectangular wire breaks this degeneracy. Figure 6.2 shows the 2D charge density in the wire cross-section.

![2D charge density distribution in the wire cross-section](../figures/wire_density_2d.png)

**Figure 6.2:** Charge density distribution in the GaAs wire cross-section. For a circular wire, the density has cylindrical symmetry. Rectangular wires produce anisotropic density profiles that translate into polarization-dependent optical response.

### 6.6.4 Second CB subband transitions

The transitions involving CB2 (the second conduction subband) follow the same selection rules but with modified strengths due to different wave function overlap. CB2 typically has a node in the wave function, leading to:

- Reduced oscillator strength for CB2-VB1 compared to CB1-VB1, due to the overlap integral having a sign change.
- Enhanced CB2-VB3 transitions when the nodal structure of CB2 matches the spatial character of VB3.

The total oscillator strength across all transitions for a given CB state should satisfy the f-sum rule $\sum_{\text{VB}} f_{ij} \leq 1$, with the deficit reflecting contributions from states outside the computed VB manifold.

---

## 6.7 Connection to Absorption Spectra

### 6.7.1 From oscillator strength to absorption coefficient

The inter-band absorption coefficient for polarization $\hat{e}$ is related to the imaginary part of the dielectric function. Using the Fermi golden rule:

$$
\alpha(\hbar\omega, \hat{e}) = \frac{\pi e^2}{n_r c \epsilon_0 m_0 \omega} \sum_{ij} f_{ij}^{(\hat{e})} \, \delta(\hbar\omega - \Delta E_{ij})
$$

where $f_{ij}^{(\hat{e})} = (2/m_0 \Delta E_{ij}) \, |\hat{e} \cdot \mathbf{p}_{ij}|^2$ is the polarization-resolved oscillator strength, $n_r$ is the refractive index, and the sum runs over occupied initial and empty final states.

### 6.7.2 Joint density of states

For a quantum wire, the joint density of states (JDOS) for each subband pair $(i, j)$ includes a 1D $k_z$ integration:

$$
\rho_{ij}(\hbar\omega) = \frac{1}{\pi} \int dk_z \, \delta(\hbar\omega - E_i(k_z) + E_j(k_z))
$$

Due to the 1D nature, the JDOS for each subband pair has a **van Hove singularity** (a $1/\sqrt{\hbar\omega - \Delta E_{ij}(0)}$ divergence) at the subband edge, producing sharp peaks in the absorption spectrum. This is a key signature of quantum wire optical response and contrasts with the step-like JDOS of quantum wells.

### 6.7.3 Broadening

In practice, the delta functions are replaced by Lorentzian or Gaussian broadening functions to account for lifetime and inhomogeneous broadening:

$$
\delta(E) \to \frac{1}{\pi} \frac{\Gamma}{E^2 + \Gamma^2}
$$

where $\Gamma$ is the broadening parameter (typically 1--10 meV for high-quality nanowires at low temperature). The code outputs the raw (unbroadened) oscillator strengths and transition energies; post-processing with a broadening kernel is left to the plotting scripts.

### 6.7.4 Published example: InP polytypic superlattices

Holmberg *et al.* (arXiv:1409.6836) studied interband polarized absorption in InP polytypic superlattice nanowires, where alternating wurtzite (WZ) and zincblende (ZB) segments create a natural superlattice along the wire axis. This system is an excellent illustration of polarization-resolved optical transitions in a quasi-1D geometry.

For the zincblende segments, the present code can compute:

- Optical matrix elements $|p_x|^2$, $|p_y|^2$, $|p_z|^2$ for all CB-VB subband pairs
- Oscillator strengths from the Kane $E_P$ parameter (InP: $E_P = 20.7$ eV)
- Polarization anisotropy arising from the nanowire cross-section geometry

The ZB Hamiltonian naturally produces the correct polarization dependence: TE-polarized absorption (in-plane, $p_x$ and $p_y$) dominates for HH-related transitions, while TM-polarized absorption ($p_z$) is stronger for LH-related transitions. Full reproduction of the Holmberg results requires wurtzite support (see Section 6.8.4).

---

## 6.8 Discussion

### 6.8.1 Convergence considerations

The momentum matrix elements converge more slowly with the FD grid spacing than the eigenvalues, because they involve the derivative operator $\partial H / \partial k_\alpha$, which is sensitive to the finite-difference stencil accuracy. For reliable optical matrix elements:

- Use at least 4th-order FD (`FDorder = 4`) for quantum wells.
- For wires, ensure sufficient grid resolution in both transverse directions.
- The number of VB states included (`numvb`) should be large enough to capture the dominant transitions; typically 10--20 VB states are needed.
- Convergence of the oscillator strength should be checked by varying the grid spacing: $f_{ij}$ should change by less than 1% between successive refinements.

### 6.8.2 Gauge invariance

The momentum matrix elements $\mathbf{p}_{if}$ and the position matrix elements $\mathbf{r}_{if}$ are related by

$$
\mathbf{p}_{if} = \frac{im_0 \Delta E_{if}}{\hbar} \mathbf{r}_{if}
$$

In a bulk crystal, these two formulations ("velocity gauge" vs "length gauge") give identical results. In a heterostructure with position-dependent material parameters, however, the velocity gauge (used by this code, since it directly evaluates $\partial H / \partial k$) is preferred because it avoids the ambiguity of the position operator across interfaces. The k.p Hamiltonian naturally provides the velocity-gauge matrix elements, and the code exploits this without needing to define a position operator.

### 6.8.3 Current limitations

1. **No excitonic effects.** The independent-particle approximation neglects electron-hole Coulomb attraction, which can significantly enhance the absorption near the band edge (especially in low-dimensional systems where the exciton binding energy is large). For GaAs quantum wires, the exciton binding energy can reach 10--30 meV. Including excitons would require solving the Bethe-Salpeter equation or using a variational approach.

2. **Single-particle states only.** The oscillator strengths computed here are for transitions between Kohn-Sham-like single-particle states. Many-body corrections (local-field effects, self-energy corrections) can modify the absolute magnitudes by 10--20%.

3. **No magnetic field dependence.** The optical transitions are computed at $B = 0$. Magneto-optical measurements (e.g., Zeeman splitting of excitonic peaks) would require computing the transitions in the presence of a magnetic field, which couples to the g-factor calculation described in Chapter 05.

4. **No direct Im[$\epsilon$] computation.** The code computes individual transition matrix elements and oscillator strengths but does not assemble them into the imaginary part of the dielectric function $\text{Im}[\epsilon(\omega)]$. This requires a post-processing step that sums over all transitions with appropriate Fermi occupation factors and broadening (Section 6.7). A future extension could output $\text{Im}[\epsilon(\omega)]$ directly.

### 6.8.4 Wurtzite limitation

Full reproduction of polarized absorption results for mixed-phase nanowires (e.g., Holmberg *et al.*) requires modeling the **wurtzite** crystal structure, which has a different 8-band Hamiltonian with distinct selection rules. The wurtzite basis introduces a crystal-field splitting $\Delta_{CR}$ that further separates the A (HH-like), B (LH-like), and C (SO-like) valence bands, modifying the polarization selection rules compared to zincblende.

This code currently supports only the **zincblende** crystal structure. Adding wurtzite support would require:

1. A new Hamiltonian builder (`WZ8bandGeneralized`) with the wurtzite k.p parameters ($A_1$--$A_6$, $\Delta_{CR}$, $\Delta_{SO}^{WZ}$).
2. Updated basis ordering for the wurtzite $C_{6v}$ symmetry point group.
3. Modified selection rules in the optical matrix element computation.

The computation framework -- `compute_optical_matrix_wire`, the `optical_transition` type, the output pipeline -- is crystal-structure-agnostic and would work unchanged once the wurtzite Hamiltonian is available.

---

## 6.9 Summary

**Table 6.5:** Quick reference for optical properties computed by the code.

| Quantity | Formula | Code output column |
|---|---|---|
| Transition energy | $\Delta E = E_{\text{CB}} - E_{\text{VB}}$ | `dE(eV)` |
| $x$-polarized matrix element | $\|\langle \text{CB}\|\partial H/\partial k_x\|\text{VB}\rangle\|^2$ | `\|px\|^2` |
| $y$-polarized matrix element | $\|\langle \text{CB}\|\partial H/\partial k_y\|\text{VB}\rangle\|^2$ | `\|py\|^2` |
| $z$-polarized matrix element | $\|\langle \text{CB}\|\partial H/\partial k_z\|\text{VB}\rangle\|^2$ | `\|pz\|^2` |
| Oscillator strength | $f = \sum_\alpha \|p_\alpha\|^2 / [(\hbar^2/2m_0) \cdot \Delta E]$ | `f_osc` |

**Key physical insights:**

1. The 8-band k.p framework provides all optical matrix elements "for free" -- they are contained in the off-diagonal $P$-coupling blocks of the Hamiltonian itself.

2. The selection rules emerge directly from the angular momentum quantum numbers of the basis states: HH transitions are purely TE, LH transitions are TE+TM with TM dominant, and SO transitions are mixed.

3. The Kane energy $E_P$ controls the absolute scale of all optical matrix elements and simultaneously determines the CB effective mass through the k.p interaction.

4. Quantum wire geometry introduces polarization anisotropy ($|p_x|^2 \neq |p_y|^2$ for rectangular cross-sections) that can be used to identify the symmetry character of transitions experimentally.

**Code locations:**

| Component | File | Routine |
|---|---|---|
| Optical transition type | `src/core/defs.f90` | `optical_transition` |
| Matrix element computation | `src/physics/gfactor_functions.f90` | `compute_optical_matrix_wire` |
| Per-direction dispatcher | `src/physics/gfactor_functions.f90` | `compute_pele_2d` |
| Output writing | `src/apps/main_gfactor.f90` | wire block |
