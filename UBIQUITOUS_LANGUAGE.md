# Ubiquitous Language

The shared vocabulary for the **8-band zinc-blende k·p FDM solver** — the words a
developer and a semiconductor physicist use to talk about this system, with the
specific meanings they carry *here*.

> **Relationship to `CONTEXT.md`.** `CONTEXT.md` is the companion glossary for
> **resolved ambiguities** — terms overloaded across the four confinement modes
> (`block`, `mode`, `geometry`, `confinement direction`, `velocity operator`,
> eigensolver dispatch, the `AUTO`/band-window contracts). This file is the
> **domain-language backbone** those ambiguities sit on top of: the bands, the
> confinement structures, the k·p model, and the observables. When a term is
> overloaded, this file links to `CONTEXT.md` rather than restating the
> disambiguation.

## Bands & basis

| Term | Definition | Aliases to avoid |
| --- | --- | --- |
| **Conduction band (CB)** | The upper band pair (basis bands 7–8) from which conduction electrons originate | Upper band |
| **Valence band (VB)** | The lower band quartet (basis bands 1–4: HH, LH, LH, HH) | Lower band |
| **Split-off band (SO)** | Bands 5–6, separated from the valence band by the spin–orbit coupling Δ_SO | |
| **Heavy hole (HH)** | The \|J=3/2, m_J=±3/2⟩ valence state (basis bands 1, 4) | |
| **Light hole (LH)** | The \|J=3/2, m_J=±1/2⟩ valence state (basis bands 2, 3) | |
| **Band gap (Eg)** | The conduction–valence energy separation at a given k-point | Energy gap |
| **Basis ordering** | The fixed mapping band 1–4 = valence, 5–6 = SO, 7–8 = CB; hardcoded throughout | Band ordering |

## Confinement structures

| Term | Definition | Aliases to avoid |
| --- | --- | --- |
| **Confinement** | The spatial restriction imposed by the structure; selects one of four modes (bulk, QW, wire, Landau) | Geometry (reserved for the wire cross-section — see `CONTEXT.md`) |
| **Bulk** | No spatial confinement — the homogeneous crystal, solved as an 8×8 Hamiltonian at a single point | |
| **Quantum well (QW)** | A 1D-confined heterostructure, discretized along the growth axis z (8N×8N Hamiltonian) | Well (ambiguous — see flagged ambiguities) |
| **Quantum wire (wire)** | A 2D-confined channel in the x–y cross-section, free along its propagation axis z | Nanowire, channel |
| **Landau** | Magnetic-confinement mode producing discrete Landau levels, discretized along x | |
| **Growth direction** | The epitaxial axis (z) along which the QW is deposited; also strain's reference axis | |
| **Propagation direction** | The wire's translationally free axis (z) — the opposite role from the QW's confinement axis | Free direction |

## The k·p model

| Term | Definition | Aliases to avoid |
| --- | --- | --- |
| **k·p Hamiltonian** | The effective-mass envelope-function Hamiltonian expanded in k about a band extremum — the model this code solves | |
| **Kane P (P)** | The interband momentum-matrix coupling between conduction and valence bands; enters the off-diagonal k·p terms | Kane element |
| **k·p term** | One labeled off-diagonal block descriptor (Q, R, S, T, P) coupling a band pair in the block table | Block (overloaded — see `CONTEXT.md`); block-formula |
| **Foreman renormalization** | An optional correction to the Luttinger-like valence terms; off by default | |

## States & transitions

| Term | Definition | Aliases to avoid |
| --- | --- | --- |
| **Eigenstate (ψ)** | A stationary envelope-function solution of the Hamiltonian at a given k | Eigenfunction, wavefunction |
| **Energy eigenvalue (E)** | The energy of an eigenstate — the scalar the eigensolver returns | Eigenvalue (solver-internal sense), level |
| **Subband** | A discrete quantized energy level arising from confinement within a single band manifold (QW/wire) | Band (a subband is a confined slice of a band) |
| **Intersubband transition (ISBT)** | An optical transition between two subbands within the same band manifold | |
| **Gamma point (k=0)** | The Brillouin-zone center; band-edge and g-factor quantities are evaluated here | Γ-point |
| **Band structure** | Energy vs. k-vector dispersion of the eigenstates | Dispersion |

## Magnetic & strain response

| Term | Definition | Aliases to avoid |
| --- | --- | --- |
| **g-factor (Landé g)** | The dimensionless factor linking a state's magnetic moment to its angular momentum; computed by second-order perturbation here | |
| **Löwdin partitioning** | The second-order perturbative scheme that folds remote-band contributions into the g-factor | |
| **Zeeman splitting** | The linear-in-B energy splitting of degenerate band states | |
| **Landau level (LL)** | A discrete, magnetic-field-quantized energy level of a free carrier | |
| **Bir–Pikus strain** | The k·p-compatible deformation-potential correction to the Hamiltonian from lattice mismatch | Strain Hamiltonian (loose) |
| **Compressive / tensile strain** | The sign of the lattice mismatch relative to the substrate; compressive raises HH and lowers LH | |

## Optical response

| Term | Definition | Aliases to avoid |
| --- | --- | --- |
| **Optical absorption** | The photon-energy-dependent rate at which the material takes up photons via inter/subband transitions | |
| **Gain** | Negative absorption, occurring under population inversion | |
| **Spontaneous emission** | The photon-emission rate independent of an incident field | |
| **Oscillator strength (f)** | The dimensionless transition weight f = (2/m₀)\|p\|²/ΔE | |
| **Velocity operator / momentum matrix** | The interband transition matrix v_α = −i[r_α,H]/ħ (≡ p_α/m₀); central to both optics and g-factor — equivalence and units are subtle, see `CONTEXT.md` | Current operator (not used here) |

## Self-consistent loop

| Term | Definition | Aliases to avoid |
| --- | --- | --- |
| **Self-consistent Schrödinger–Poisson (SP)** | The iterative loop coupling the quantum eigenstates to the electrostatic potential they produce | SC loop |
| **DIIS / Pulay mixing** | The direct-inversion-in-iterative-subspace acceleration scheme used to converge the SP loop | |
| **Charge density** | The carrier density n(z)/p(z) built from the k·p eigenstates and k∥ sampling | |
| **Doping** | Intentional donor/acceptor charge that pins the Fermi level and bends the bands | Donor/acceptor (those are the species, not the amount) |
| **Fermi level** | The chemical potential setting state occupation | Fermi energy |

## Topology & superconductivity

| Term | Definition | Aliases to avoid |
| --- | --- | --- |
| **Chern number** | The integer topological invariant of a 2D quantum-Hall system (computed via Fukui–Hatsugai–Suzuki here) | |
| **Z2 invariant** | The binary time-reversal topological invariant of a quantum-spin-Hall system | |
| **BdG (Bogoliubov–de Gennes)** | The mean-field superconducting Hamiltonian written in Nambu (particle–hole) space | |
| **Nambu space** | The doubled particle–hole basis in which BdG lives (16N×16N, four 8N×8N quarters) | |
| **Majorana mode** | A self-conjugate zero-energy edge state of a topological superconductor | |

## Discretization

| Term | Definition | Aliases to avoid |
| --- | --- | --- |
| **Finite difference (FDM)** | The spatial-discretization scheme that turns the continuous k·p Hamiltonian into a matrix | |
| **FD order** | The stencil order (2–10) controlling discretization accuracy | |
| **Grid** | The discrete set of spatial points on which the Hamiltonian is built (`grid%npoints()` total) | Mesh |

## Numerics

- **Gershgorin bound**: the upper bound on |eigenvalue| of a matrix derived from the sum of absolute values in each row. For a CSR matrix with diagonal + off-diagonal pattern, the bound scales with the FD stencil coefficient and the grid spacing (typically ±tens of eV for the 8-band k.p Hamiltonian). When fed to FEAST as an energy window, the Gershgorin scale samples the FD-Nyquist tail (see below) and returns spurious states. Cross-reference: `docs/solutions/best-practices/2026-06-21-fd-nyquist-spurious-tail.md`.

- **sentinel gap value**: the value `min_gap = -1.0_dp` (or, in error-stop variants, a non-zero exit) returned by `run_bdg_wire` / `eval_wire_bdg_gap` when the FEAST eigensolve finds 0 eigenvalues in the configured physics window. The sentinel means "no BdG modes found in this window" — typically caused by μ in the band gap or a mis-sized `[solver] emin/emax`. Consumers (lecture scripts, phase-diagram tools) MUST treat the sentinel as a configuration error, NOT as a real negative gap.

- **FD-Nyquist tail**: spurious high-energy eigenstates produced when a wide energy window (e.g., the Gershgorin scale) is fed to FEAST on an FD-discretized Hamiltonian. The FD stencil introduces non-physical dispersion at high k that FEAST samples indiscriminately. Detection: `min_gap` near the window boundary rather than near the pairing gap. Mitigation: use a physics-sized window (≤ ~50 meV for BdG). Cross-reference: `docs/solutions/best-practices/2026-06-21-fd-nyquist-spurious-tail.md`.

- **PHS operator (Xi = tau_x K)**: the particle-hole symmetry operator used in BdG cross-checks is `Xi = tau_x K = [[0, I_n], [I_n, 0]] K`, where K is complex conjugation. This is the correct operator for a wire-builder BdG of the form `H = [[H_ee, Delta], [Delta^*, -H_ee^T]]` — the swap-and-conjugate maps the electron block to the negated hole block. The simpler `C = diag(I_n, -I_n) K` does NOT satisfy PHS for this builder form. Verified numerically by `test_bdg_phs_at_finite_bx` (C8, commit `214acf4`).

## Relationships

- A **Confinement** structure is one of **bulk**, **quantum well**, **quantum wire**, or **Landau**; it fixes the matrix size (8×8 for bulk, 8N×8N otherwise).
- A **k·p term** (Q, R, S, T, P) couples one band pair and is assembled into a **spatial block** of the built Hamiltonian (see `CONTEXT.md`).
- Confinement splits each **band** into a ladder of **subbands**; an **ISBT** connects two subbands within one band manifold.
- A **velocity operator / momentum matrix** drives both **oscillator strength** (→ optical absorption/gain/emission) and the **g-factor** via **Löwdin partitioning**.
- **Bir–Pikus strain** and **Zeeman splitting** are perturbative corrections added to the same k·p Hamiltonian.
- The **self-consistent SP** loop wraps the eigensolve: **doping** sets the **Fermi level**, eigenstates give the **charge density**, which feeds Poisson, which bends the bands.
- **BdG** extends the Hamiltonian into **Nambu space**; its edge spectrum hosts **Majorana modes**, classified by **Chern**/**Z2** invariants.

## Example dialogue

> **Dev:** "To compute the **g-factor** for this **quantum wire**, do I diagonalize at the **Gamma point** like the **QW**?"
>
> **Domain expert:** "No. The **QW** is 1D-confined along its **growth direction**, so its **subbands** are read off at k∥=0. The **wire** is 2D-confined in the cross-section and free along its **propagation direction**, so you sweep kz — there's no single Gamma point for it."
>
> **Dev:** "Same **Löwdin partitioning** folding in the remote **split-off band**?"
>
> **Domain expert:** "Same second-order **Löwdin**, yes. But the **velocity operator** for the confined directions comes from the commutator −i[r,H] on the CSR matrix, while the free-z direction uses the k-derivative ∂H/∂k. Same operator, two sources."
>
> **Dev:** "Does the strained InAs/GaAs case change the **band gap** I target?"
>
> **Domain expert:** "**Bir–Pikus strain** shifts **HH** and **LH** differently — compressive strain raises HH and lowers LH — so the effective **band gap** you read off depends on which valence **subband** you track."

## Flagged ambiguities

- **Two glossaries now exist — keep the split clean.** `CONTEXT.md` = resolved *ambiguities* (overloaded across modes); this file = the *domain-language backbone*. Do not migrate resolved ambiguities here, and do not restate them — link. If a term earns an entry in both, the backbone entry stays a one-liner and defers to `CONTEXT.md`.
- **"Well" is overloaded.** It means both the **quantum well** structure and the narrow **well layer** of material inside it. Prefer **quantum well** for the structure and **well layer** for the material slab; never bare "the well" where either could be meant.
- **"Eigenvalue" / "energy" / "level"** are three words for the same scalar. Use **energy eigenvalue** (or just **energy**) for the physical quantity; reserve **eigenvalue** for the solver-internal, returned-by-the-eigensolver sense (which may be clamped — see `CONTEXT.md`'s band-window entry).
- **"Band structure" vs "dispersion"** name the same thing. Prefer **band structure**; reserve **dispersion** when contrasting against a single subband's curvature.
- **"Bandedge" (one word) vs "band edge" (two words)** appear in different registers — test names use the compound, prose uses the two-word form. Prefer **band edge**; it means the k=0 energy of a band, not the whole edge region.
- **"Doping" vs "donor/acceptor".** Doping is the act/amount (a `doping_spec`); donor/acceptor are the species. Don't say "increase the donors" when you mean "raise the doping."
- The mode-level overloads (**block**, **mode**, **geometry**, **confinement direction**) are already disambiguated in `CONTEXT.md` and are not repeated here.
