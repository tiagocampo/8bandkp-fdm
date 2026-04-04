# Chapter 06: Optical Properties

## 6.1 Introduction

When light interacts with a semiconductor, photons are absorbed by promoting electrons from occupied valence band states to empty conduction band states. The strength of these inter-band transitions -- and their dependence on the polarization direction of the incident light -- is entirely governed by the **momentum matrix elements** between the initial and final Bloch states.

In the 8-band k.p framework, these matrix elements are computed not from first-principles wave functions but from the **Kane momentum parameter** $E_P$ (or equivalently $P$), which already encodes the strength of the conduction-to-valence band coupling. The k.p Hamiltonian itself, evaluated at a unit perturbation wave vector, acts as a natural operator for computing these matrix elements: the off-diagonal blocks coupling CB (bands 7--8) to VB (bands 1--6) contain precisely the $P$-dependent terms that generate optical transitions.

This chapter derives the oscillator strength formula, explains the polarization-dependent selection rules that emerge from the 8-band zincblende basis, describes the implementation in `compute_optical_matrix_wire` and `compute_pele_2d`, and discusses the connection to measurable absorption spectra.

---

## 6.2 Momentum Matrix Elements

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

where $E_P$ is the material-specific Kane energy stored in `params%EP` (e.g., $E_P = 28.8$ eV for GaAs, $E_P = 21.5$ eV for InAs). The code computes this in `parameters.f90` as `params(i)%P = dsqrt(params(i)%EP * const)` where `const = hbar^2 / (2*m_0)` = 3.80998 eV*AA^2.

---

## 6.3 Oscillator Strength

### 6.3.1 Definition

The dimensionless **oscillator strength** for the inter-band transition $i \to f$ with transition energy $\Delta E = E_f - E_i$ is

$$
f_{if} = \frac{2}{m_0 \Delta E} \sum_{\alpha = x,y,z} |\langle \psi_i | p_\alpha | \psi_f \rangle|^2
$$

In the code's dH/dk units, this becomes

$$
f_{if} = \frac{\sum_\alpha |P_\alpha^{if}|^2}{(\hbar^2 / 2m_0) \cdot \Delta E}
$$

where $\hbar^2 / 2m_0$ is the constant `hbar2O2m0` defined in `defs.f90`.

### 6.3.2 Thomas-Reiche-Kuhn sum rule

For a complete set of states, the oscillator strengths satisfy the f-sum rule:

$$
\sum_f f_{if} = 1
$$

for any initial state $i$. In practice, the sum runs over a finite number of computed bands, so the sum will be less than unity. For inter-band optical transitions (CB-VB pairs), the dominant contribution comes from the Kane coupling $E_P$, which accounts for most of the CB effective mass via

$$
\frac{m_0}{m^*} \approx 1 + \frac{E_P}{3} \left(\frac{2}{E_g} + \frac{1}{E_g + \Delta_{SO}}\right)
$$

This shows that the oscillator strength is intimately connected to the band effective mass through the same Kane parameter that appears in the Hamiltonian.

### 6.3.3 Computation in the code

The routine `compute_optical_matrix_wire` in `gfactor_functions.f90` computes oscillator strengths for all CB-VB pairs:

```fortran
! For each CB-VB pair (i, j):
dE = cb_value(i) - vb_value(j)

! Compute |<CB_i|dH/dk_alpha|VB_j>|^2 for alpha = x, y, z
do dir = 1, 3
  call compute_pele_2d(Pele, dir, cb_state(:,i), vb_state(:,j), ...)
  select case(dir)
  case(1); transitions(idx)%px = real(Pele * conjg(Pele))
  case(2); transitions(idx)%py = real(Pele * conjg(Pele))
  case(3); transitions(idx)%pz = real(Pele * conjg(Pele))
  end select
end do

! Oscillator strength
sum_p = px + py + pz
f = sum_p / (hbar2O2m0 * dE)
```

The output is written to `output/optical_transitions.dat` with columns:

```
# CB VB dE(eV) |px|^2 |py|^2 |pz|^2 f_osc
```

---

## 6.4 Polarization Dependence and Selection Rules

### 6.4.1 TE and TM polarization

For a quantum well grown along $z$, the polarization of the absorbed photon determines which momentum matrix element is probed:

- **TE polarization** (electric field in the plane, $\hat{e} = \hat{x}$ or $\hat{y}$): couples via $p_x$ and $p_y$
- **TM polarization** (electric field along $z$, $\hat{e} = \hat{z}$): couples via $p_z$

The absorption coefficient for polarization $\hat{e}$ is proportional to

$$
\alpha(\hbar\omega) \propto \sum_{i,j} |\hat{e} \cdot \mathbf{p}_{ij}|^2 \, \delta(\hbar\omega - \Delta E_{ij})
$$

### 6.4.2 Selection rules from the 8-band basis

At the zone center ($\mathbf{k} = 0$), the angular momentum character of the basis states imposes strict selection rules. In the zincblende basis used by this code (ordering: HH, LH, LH, HH, SO, SO, CB, CB), the dominant CB-to-VB optical transitions obey:

| Transition | $p_x$ | $p_y$ | $p_z$ | Polarization |
|---|---|---|---|---|
| CB $\uparrow$ $\to$ HH $\uparrow$ ($|7\rangle \to |1\rangle$) | $P/\sqrt{2}$ | $iP/\sqrt{2}$ | 0 | TE |
| CB $\uparrow$ $\to$ LH $\uparrow$ ($|7\rangle \to |2\rangle$) | 0 | 0 | $\sqrt{2/3}\,P$ | TM |
| CB $\uparrow$ $\to$ LH $\downarrow$ ($|7\rangle \to |3\rangle$) | $P/\sqrt{3}$ | $-iP/\sqrt{3}$ | 0 | TE |
| CB $\uparrow$ $\to$ SO $\uparrow$ ($|7\rangle \to |5\rangle$) | 0 | 0 | $P/\sqrt{3}$ | TM |
| CB $\uparrow$ $\to$ SO $\downarrow$ ($|7\rangle \to |6\rangle$) | $P/\sqrt{6}$ | $-iP/\sqrt{6}$ | 0 | TE |

Key observations:

1. **HH-CB transitions** are purely TE-polarized ($p_x$ and $p_y$ only, no $p_z$). This is because the HH states ($|J_z| = 3/2$) couple to the CB ($|J_z| = 1/2$) exclusively through the in-plane momentum components. The ratio $|p_x|^2 : |p_y|^2 = 1:1$ gives isotropic in-plane absorption.

2. **LH-CB transitions** have both TE and TM components. The LH states ($|J_z| = 1/2$) mix in-plane and out-of-plane coupling. Specifically, $|p_z|^2 : (|p_x|^2 + |p_y|^2) = 2:1$ for LH, meaning TM is twice as strong as TE.

3. **SO-CB transitions** similarly mix TE and TM, with the ratio depending on the specific $J_z$ quantum numbers involved.

4. **Spin-flip transitions** (where the spin quantum number changes) are forbidden at $\mathbf{k} = 0$ by time-reversal symmetry. For example, $|7\rangle \to |4\rangle$ (CB $\uparrow \to$ HH $\downarrow$) has zero matrix element at $k = 0$.

### 6.4.3 Breakdown of selection rules at finite k

At finite wave vector, the k.p mixing of band character causes the strict zone-center selection rules to relax. A CB state that is purely $|S, +1/2\rangle$ at $k=0$ acquires HH, LH, and SO admixtures at finite $k$, enabling transitions that would be forbidden at $k=0$. The code captures this naturally because the eigenvectors of the full Hamiltonian already contain the correct band mixing.

### 6.4.4 Polarization anisotropy in quantum wires

For a quantum wire with confinement in $x$ and $y$, the $C_{\infty v}$ (or lower) symmetry of the wire cross-section lifts the degeneracy between $p_x$ and $p_y$ matrix elements. The code resolves all three components ($|p_x|^2$, $|p_y|^2$, $|p_z|^2$) separately, capturing:

- **$x$-polarized absorption** (proportional to $|p_x|^2$): sensitive to the horizontal wire dimension
- **$y$-polarized absorption** (proportional to $|p_y|^2$): sensitive to the vertical wire dimension
- **$z$-polarized absorption** (proportional to $|p_z|^2$): probes the free-propagation axis

For a rectangular wire, the anisotropy $|p_x|^2 \neq |p_y|^2$ can be significant, especially for transitions involving states with strong spatial localization along one confinement direction.

---

## 6.5 Wire Optical Transitions

### 6.5.1 Inter-subband transitions in quantum wires

In a quantum wire, the 2D confinement quantizes the transverse motion into discrete subbands, while the wave vector $k_z$ along the wire axis remains a good quantum number. The optical transitions are computed between **all CB-VB pairs** at a fixed $k_z$ (typically $k_z = 0$):

$$
\text{transitions} = \{(\text{CB}_i, \text{VB}_j) \mid i = 1,\ldots,N_{\text{CB}};\; j = 1,\ldots,N_{\text{VB}}\}
$$

Each transition is characterized by the `optical_transition` derived type:

```fortran
type :: optical_transition
  integer :: cb_idx, vb_idx         ! band indices
  real(dp) :: energy                ! transition energy (eV)
  real(dp) :: px, py, pz           ! |<i|dH/dk_alpha|j>|^2
  real(dp) :: oscillator_strength  ! f = sum(|p|^2) / (hbar2O2m0 * dE)
end type
```

### 6.5.2 Implementation details

The routine `compute_optical_matrix_wire` uses `compute_pele_2d` to evaluate each momentum matrix element. The perturbation Hamiltonian for direction $\alpha$ is built by `pMatrixEleCalc_2d`, which calls `ZB8bandGeneralized` with the appropriate `g` flag:

- `g='g1'`: perturbation along $x$ (the $d/dx$ gradient operator)
- `g='g2'`: perturbation along $y$ (the $d/dy$ gradient operator)
- `g='g3'`: perturbation along $z$ (the $k_z$ operator)

Each perturbation Hamiltonian is stored in CSR format and applied via sparse matrix-vector multiplication (`csr_spmv`), making the computation efficient for the large $8N_{\text{grid}} \times 8N_{\text{grid}}$ wire Hamiltonian.

### 6.5.3 Interpreting the output

The file `output/optical_transitions.dat` lists all $N_{\text{CB}} \times N_{\text{VB}}$ transitions. Typical observations for a GaAs-based quantum wire:

1. **Dominant transitions**: The lowest CB subband to the topmost VB subband has the largest oscillator strength, because these states have the strongest spatial overlap and the smallest transition energy.

2. **Weak transitions**: Higher subbands (CB2, CB3, ...) coupling to deep VB states often have negligible oscillator strength due to poor wave function overlap.

3. **Dark transitions**: Some transitions have $f \approx 0$ due to symmetry-forbidden selection rules (e.g., odd-to-odd parity combinations in symmetric wires).

4. **Polarization signature**: The ratio $|p_x|^2 / |p_y|^2$ deviates from unity for non-circular wire cross-sections, providing a polarization-dependent absorption spectrum.

---

## 6.6 Connection to Absorption Coefficient

### 6.6.1 From oscillator strength to absorption

The inter-band absorption coefficient for polarization $\hat{e}$ is related to the imaginary part of the dielectric function. Using the Fermi golden rule:

$$
\alpha(\hbar\omega, \hat{e}) = \frac{\pi e^2}{n_r c \epsilon_0 m_0 \omega} \sum_{ij} f_{ij}^{(\hat{e})} \, \delta(\hbar\omega - \Delta E_{ij})
$$

where $f_{ij}^{(\hat{e})} = (2/m_0 \Delta E_{ij}) \, |\hat{e} \cdot \mathbf{p}_{ij}|^2$ is the polarization-resolved oscillator strength, $n_r$ is the refractive index, and the sum runs over occupied initial and empty final states.

### 6.6.2 Joint density of states

For a quantum wire, the joint density of states (JDOS) for each subband pair $(i, j)$ includes a 1D $k_z$ integration:

$$
\rho_{ij}(\hbar\omega) = \frac{1}{\pi} \int dk_z \, \delta(\hbar\omega - E_i(k_z) + E_j(k_z))
$$

Due to the 1D nature, the JDOS for each subband pair has a **van Hove singularity** (a $1/\sqrt{\hbar\omega - \Delta E_{ij}(0)}$ divergence) at the subband edge, producing sharp peaks in the absorption spectrum. This is a key signature of quantum wire optical response.

### 6.6.3 Broadening

In practice, the delta functions are replaced by Lorentzian or Gaussian broadening functions to account for lifetime and inhomogeneous broadening:

$$
\delta(E) \to \frac{1}{\pi} \frac{\Gamma}{E^2 + \Gamma^2}
$$

where $\Gamma$ is the broadening parameter. The code outputs the raw (unbroadened) oscillator strengths and transition energies; post-processing with a broadening kernel is left to the plotting scripts.

---

## 6.7 Published Example: Interband Polarized Absorption in InP Polytypic Superlattices

### 6.7.1 Context

Holmberg *et al.* (arXiv:1409.6836) studied interband polarized absorption in InP polytypic superlattice nanowires, where alternating wurtzite (WZ) and zincblende (ZB) segments create a natural superlattice along the wire axis. This system is an excellent illustration of polarization-resolved optical transitions in a quasi-1D geometry.

### 6.7.2 What the code can compute

For the zincblende segments, the present code can compute:

- **Optical matrix elements** $|p_x|^2$, $|p_y|^2$, $|p_z|^2$ for all CB-VB subband pairs
- **Oscillator strengths** from the Kane $E_P$ parameter (InP: $E_P = 20.7$ eV)
- **Polarization anisotropy** arising from the nanowire cross-section geometry

The ZB Hamiltonian naturally produces the correct polarization dependence: TE-polarized absorption (in-plane, $p_x$ and $p_y$) dominates for HH-related transitions, while TM-polarized absorption ($p_z$) is stronger for LH-related transitions.

### 6.7.3 Limitation: wurtzite crystal structure

Full reproduction of the polarized absorption results in Holmberg *et al.* requires modeling the **wurtzite** crystal structure, which has a different 8-band Hamiltonian with distinct selection rules. The wurtzite basis introduces a crystal-field splitting $\Delta_{CR}$ that further separates the A (HH-like), B (LH-like), and C (SO-like) valence bands, modifying the polarization selection rules compared to zincblende.

This code currently supports only the **zincblende** crystal structure. Adding wurtzite support would require:

1. A new Hamiltonian builder (`WZ8bandGeneralized`) with the wurtzite k.p parameters ($A_1$--$A_6$, $\Delta_{CR}$, $\Delta_{SO}^{WZ}$)
2. Updated basis ordering for the wurtzene $C_{6v}$ symmetry point group
3. Modified selection rules in the optical matrix element computation

The computation framework -- `compute_optical_matrix_wire`, the `optical_transition` type, the output pipeline -- is crystal-structure-agnostic and would work unchanged once the wurtzite Hamiltonian is available. This is noted as a TODO for future development.

---

## 6.8 Discussion

### 6.8.1 Convergence considerations

The momentum matrix elements converge more slowly with the FD grid spacing than the eigenvalues, because they involve the derivative operator $\partial H / \partial k_\alpha$, which is sensitive to the finite-difference stencil accuracy. For reliable optical matrix elements:

- Use at least 4th-order FD (`FDorder = 4`) for quantum wells
- For wires, ensure sufficient grid resolution in both transverse directions
- The number of VB states included (`numvb`) should be large enough to capture the dominant transitions; typically 10--20 VB states are needed

### 6.8.2 Gauge invariance

The momentum matrix elements $\mathbf{p}_{if}$ and the position matrix elements $\mathbf{r}_{if}$ are related by

$$
\mathbf{p}_{if} = \frac{im_0 \Delta E_{if}}{\hbar} \mathbf{r}_{if}
$$

In a bulk crystal, these two formulations ("velocity gauge" vs "length gauge") give identical results. In a heterostructure with position-dependent material parameters, however, the velocity gauge (used by this code, since it directly evaluates $\partial H / \partial k$) is preferred because it avoids the ambiguity of the position operator across interfaces.

### 6.8.3 Limitations

1. **No excitonic effects**: The independent-particle approximation neglects electron-hole Coulomb attraction, which can significantly enhance the absorption near the band edge (especially in low-dimensional systems). Including excitons would require solving the Bethe-Salpeter equation or using a variational approach.

2. **Single-particle states only**: The oscillator strengths computed here are for transitions between Kohn-Sham-like single-particle states. Many-body corrections (local-field effects, self-energy corrections) can modify the absolute magnitudes.

3. **No magnetic field dependence**: The optical transitions are computed at $B = 0$. Magneto-optical measurements (e.g., Zeeman splitting of excitonic peaks) would require computing the transitions in the presence of a magnetic field, which couples to the g-factor calculation described in Chapter 05.

---

## 6.9 Summary

| Quantity | Formula | Code output column |
|---|---|---|
| Transition energy | $\Delta E = E_{\text{CB}} - E_{\text{VB}}$ | `dE(eV)` |
| $x$-polarized matrix element | $\|\langle \text{CB}|\partial H/\partial k_x|\text{VB}\rangle\|^2$ | `|px|^2` |
| $y$-polarized matrix element | $\|\langle \text{CB}|\partial H/\partial k_y|\text{VB}\rangle\|^2$ | `|py|^2` |
| $z$-polarized matrix element | $\|\langle \text{CB}|\partial H/\partial k_z|\text{VB}\rangle\|^2$ | `|pz|^2` |
| Oscillator strength | $f = \sum_\alpha |p_\alpha|^2 / [(\hbar^2/2m_0) \cdot \Delta E]$ | `f_osc` |

The key physical insight is that the 8-band k.p framework provides all optical matrix elements "for free" -- they are contained in the off-diagonal $P$-coupling blocks of the Hamiltonian itself. The code extracts them by evaluating the Hamiltonian at a unit perturbation wave vector, making the computation both exact (within the k.p model) and efficient (using the same sparse infrastructure as the eigenvalue solver).
