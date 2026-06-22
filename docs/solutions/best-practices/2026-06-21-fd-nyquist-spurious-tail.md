---
title: "FD-Nyquist spurious-state tail in finite-difference k.p spectra (and why it is not a bug)"
date: "2026-06-21"
category: best-practices
module: src/physics/hamiltonian_wire.f90
problem_type: best_practice
component: physics_engine
severity: low
applies_when:
  - "A wire or QW k.p dense/FEAST spectrum shows eigenvalues spread over +/- tens of eV and you suspect the Hamiltonian build is broken"
  - "Choosing an eigensolver energy window for a confined (wire/QW) BdG or normal-state solve"
  - "An auto/wide energy window returns unphysical eigenvalues (e.g. -29 eV for GaAs, which lives in [-0.8, +0.719] eV)"
  - "Vetting a near-zero BdG eigenstate as a real Majorana vs a grid-scale spurious mode"
tags: [finite-differences, kp-model, spurious-solutions, nyquist, eigensolver-window, wire, bdg, convergence, kdotpy]
related_components:
  - src/physics/hamiltonianConstructor.f90
  - src/physics/confinement_init.f90
  - src/apps/main_topology.f90
  - src/math/finitedifferences.f90
---

# FD-Nyquist spurious-state tail in finite-difference k.p spectra

## Context

A dense or FEAST diagonalization of an 8-band k.p wire (or QW) Hamiltonian routinely
returns a spectrum spread over **±(20…80) eV**, with a wide empty region around the
physical band gap. On first sight this looks catastrophic — "the wire build must be
completely wrong" — and it triggers an understandable panic. It is **not** a bug. It is
the expected **Nyquist spurious-solutions tail** of any finite-difference k.p
discretization, and it is harmless to low-energy physics provided you interrogate the
spectrum correctly.

This doc records the mechanism, a quantitative estimate that closes the loop, the
failure modes where it *does* bite, and the mitigations — so the next person who sees a
±80 eV wire spectrum diagnoses it in five minutes instead of two hours.

## The mechanism

A finite-difference grid of spacing `dx` can only represent momenta up to the Nyquist
wavevector `k_N = π/dx`. The 2nd-order FD second derivative on a plane wave is not
`-k²` but

```
d²/dx²  →  -(4/dx²) sin²(k dx / 2)
```

which equals `-k²` only for `k·dx ≪ 1` and **saturates at `-4/dx²`** at `k = π/dx`.
The k.p kinetic terms (`const·A·k²` for conduction via `A = 1/m*`, `const·γ·k²` for
valence via the Luttinger `γ`, with `const = ℏ²/2m₀ = 3.81 eV·Å²`) therefore produce a
batch of eigenstates at the Nyquist scale whose wavefunctions are grid-scale
oscillations `c_i ∝ (-1)ⁱ`. Those are the spurious solutions.

### Quantitative check (GaAs wire, dx = 3 Å, both transverse dirs at Nyquist)

With `k_x² + k_y² → 8/dx²` and `dx² = 9 Å²`:

| band family | coefficient | Nyquist energy `const·(coef)·8/dx²` |
|---|---|---|
| valence (γ₁ ≈ 7) | `const·γ₁` | 3.81 · 7 · 8/9 ≈ **24 eV** |
| conduction (A = 1/0.067 = 14.9) | `const·A` | 3.81 · 14.9 · 8/9 ≈ **50 eV** |

Total spread ≈ ±(24…50) eV → ~80 eV range, matching the measured ~84 eV exactly. **The
estimate closing the loop is the proof of what the tail is.**

## How to recognize it (diagnostic signatures)

- **Constant across wire *size*** at fixed `dx` (it scales as `1/dx²`, not with `L`).
  A size scan (widen the wire) leaves the spread unchanged while the band gap shrinks
  toward bulk Eg — the opposite of what a broken build would do.
- **Scales as `1/dx²`** when `dx` is refined.
- **Wavefunctions are alternating-sign grid oscillations**, not smooth physical states.
- The **low-energy bands converge** to bulk Eg as the wire widens (physical confinement
  ∝ 1/L²), confirming the build is correct underneath the tail.

Empirical confirmation on this tree (2026-06-21, GaAs wire, dense LAPACK): a size scan
gave gap 3.97 / 2.44 / 2.12 eV at 3.3 / 4.5 / 6.3 nm (→ bulk Eg = 1.519), with the
±84 eV spread identical at every size. FEAST and dense agreed to 5 decimals on the
same matrix, so the tail is a property of the *matrix*, not the solver.

## Does it interfere with the result? — conditionally

**No, for low-energy physics.** The tail sits at |E| ~ tens of eV; band edges, subbands,
effective masses, and the BdG minigap at E ≈ 0 live within ±~0.1 eV of the gap. The two
manifolds are separated by 2–3 orders of magnitude and do not hybridize (spurious
wavefunctions are `(-1)ⁱ`, projecting ≈ 0 onto smooth states). Anything computed *from
the low-energy manifold* is trustworthy.

**Yes, in three failure modes — all about how you interrogate the spectrum:**

1. **Wide eigensolver windows sample the tail.** This is the live failure mode. The
   ±69 eV Gershgorin *auto-window* in `run_bdg_wire` returns unphysical eigenvalues
   (e.g. −29 eV for GaAs) and misses the ±E particle-hole partners, because it
   diagonalizes inside the spurious tail. **Never use the auto-window for BdG.** Always
   solve in a tight window around the physics (`±50 meV` for BdG, `±5δ₀` for the
   minigap sweep). `min_gap = 2·min|E|` over the found set is then physical.
2. **Spectrum-wide observables** (sum rules, full-spectrum partition functions,
   `minval(|E|)` over a wide solve) are contaminated — restrict them to the low-energy
   subspace explicitly.
3. **A spurious state near E ≈ 0** (rare; narrow-gap materials, strong mixing, coarse
   `dx`) can masquerade as a near-zero mode and be **mis-identified as a Majorana**.
   Vet near-zero candidates by wavefunction: a real MZM is smooth and edge-localized; a
   spurious state alternates sign site-to-site. Sticlet polarization (plan unit U6) is a
   complementary discriminator.

## Mitigations (best → worst for this codebase)

1. **Tight eigensolver windows** — costs nothing, fully effective for low-energy work.
   Keep the auto-window disabled on the BdG path (already implied by ADR 0005 / KTD6).
2. **Raise `FDorder` (2 → 4 or 6).** The code supports orders 2–10. A higher-order
   `d²/dx²` stencil keeps `k²_FD ≈ k²` accurate closer to `k_N`, pushing the spurious
   onset to higher energy and cleaning the band-edge region at fixed `dx`. Highest-leverage
   numerical knob available.
3. **Refine `dx`.** Raises the Nyquist energy, pushing the tail further from any window
   you might need. Cost: larger matrices.
4. **Wavefunction-character filtering** — drop `(-1)ⁱ`-dominant eigenstates before
   computing observables; a defensive layer if a wide window is ever unavoidable.
5. **Plane-wave / Fourier basis (kdotpy).** No FD lattice → no FD-Nyquist; it has a
   reciprocal cutoff `G_max` managed differently. This is why kdotpy is the clean
   reference and the deferred `hzy` wire comparison is the gold-standard check.
6. **Foreman renormalization** addresses a *different* spurious class (continuum
   k.p-construction valence divergence), **not** FD-Nyquist — see the dormant-units-bug
   backlog item before touching it.

## What this is NOT

- Not a wire-build bug — the size-scan vindicates the construction (gap → bulk Eg).
- Not a solver bug — FEAST and dense agree exactly.
- Not fixed by adding a core/shell — a proper barrier adds confinement energy (bands
  push apart, as it should) but does **not** remove the tail; the tail is intrinsic to
  FD k.p at the chosen `dx`, shell or no shell (verified: GaAs/AlAs core/shell has the
  same ±45 eV spread as shell-less).
- Not removed by Foreman renormalization (wrong spurious class).

## Related

- Backlog: dormant Foreman-renormalization `const` units bug (`parameters.f90:798-800,814`
  vs consumer `hamiltonianConstructor.f90:473`) — a *different* latent issue.
- Memory: `project_bdg_allzero_mu_misparam` — the μ-in-gap all-zero is orthogonal to
  this tail; the tail only bit via the auto-window failure mode.
- ADR 0005 / KTD6 — one stable, physics-sized window per sweep via `apply_solver_window`,
  never the auto-window.
