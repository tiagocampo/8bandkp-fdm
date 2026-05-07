# SPEC: Topological Suite Verification and Lecture Documentation

## Context

The topological suite (`feature/bdg-topological-superconductivity`) implements five physics capabilities but most are unverified or broken. The goal is to **verify each physics module by producing figures comparing computed output vs. literature benchmarks**, in a **pedagogical sequence from simpler to more complex physics**.

| # | Physics | Status | Evidence |
|---|---|---|---|
| 1 | QWZ Chern number (QHE) | **WORKING ✓** | C=+1,-1,0 for u=-0.8,0.5,2.5 — verified |
| 2 | BHZ Z₂ invariant (QSHE) | **BROKEN ✗** | Both trivial and topological return Z₂=0 |
| 3 | Landau levels (magnetic field) | **NOT IMPLEMENTED** | Peierls stub; no test config |
| 4 | BdG spectrum + Majorana | **UNTESTED** | Config parser bug; no integration test |
| 5 | LDOS | **PARTIAL** | Code exists; no benchmark comparison |

---

## Pedagogical Sequence

Each section: **(1) Physics theory → (2) Code computes → (3) Literature benchmark → (4) Figure comparison → (5) Next section**

---

## Phase 1: QWZ Chern Number — Figures from Working Code

### What is working
- Chern number C = +1, -1, 0 for u = -0.8, 0.5, 2.5
- Fukui-Hatsugai-Suzuki lattice gauge method
- QWZ model: `H = sin kx σ_x + sin ky σ_y + (u + cos kx + cos ky) σ_z`

### What to produce

1. **Chern number convergence figure**: C vs. grid resolution (nk=20, 30, 40, 50, 60, 70) — shows C stabilizes at integer as nk increases
2. **Chern number phase diagram figure**: C vs. u parameter — shows topological phases across transitions at u=-2, u=0
3. **U-link / Berry curvature visualization**: Ω(kx,ky) heatmap for u=-0.8 (topological) vs. u=2.5 (trivial)

### Implementation

Create `scripts/verify_qwz_chern.py`:
1. Loop over u values: -0.8, 0.5, 2.5
2. Loop over grid resolutions: nk = 20, 30, 40, 50, 60, 70
3. Run `topologicalAnalysis` with each combination
4. Parse Chern number from output
5. Generate figures:
   - `chern_convergence.png`: C vs. nk (all 3 u values) — shows C→integer as nk increases
   - `chern_phase_diagram.png`: C vs. u — shows phase transitions at u=-2 and u=0
   - `chern_berry_curvature_topological.png`: Ω(kx,ky) heatmap for u=-0.8
   - `chern_berry_curvature_trivial.png`: Ω(kx,ky) heatmap for u=2.5

### Literature comparison

| u value | Literature C | Computed C | Phase |
|---------|-------------|------------|-------|
| u=-0.8 | +1 | +1 | Topological |
| u=0.5 | -1 | -1 | Topological (inverted) |
| u=2.5 | 0 | 0 | Trivial |

Literature: Qi, Wu & Zhang, Phys. Rev. B 74, 085308 (2006)

### Files to create
- `scripts/verify_qwz_chern.py` — verification script
- `docs/lecture/figures/chern_convergence.png` — convergence figure
- `docs/lecture/figures/chern_phase_diagram.png` — phase diagram
- `docs/lecture/figures/chern_berry_curvature_topological.png` — Ω for topological phase
- `docs/lecture/figures/chern_berry_curvature_trivial.png` — Ω for trivial phase

### Update Lecture 13
- Section 13.2: Chern number — add convergence and phase diagram figures
- Table: Literature (C=+1,-1,0) vs. Computed (C=+1,-1,0) with PASS

---

## Phase 2: Fix BHZ Z₂ Invariant — Fu-Kane Parity Method

### Root Cause

The gap-based counting method `compute_z2_gap` counts eigenvalues in `[-10, +10]` meV window. The current Hamiltonian produces eigenvalues at ~±10 meV in both phases — **the edge states near E=0 are never produced**. The Fu-Kane parity method avoids the edge-state question entirely by computing Z₂ from the parity of occupied band eigenvalues at TRIM points.

### Fix Implementation

**Step 2.1: Complete Fu-Kane method for QW mode (2D)**

The QW mode (confinement=1) with `compute_z2: T` and `z2_method: fukane` should:
1. Diagonalize H(k) at 4 TRIM points: Γ=(0,0), M₁=(π/a,0), M₂=(0,π/a), M₃=(π/a,π/a)
2. For each occupied band, compute parity eigenvalue ξ = ±1 (inversion symmetry in zinc-blende basis)
3. Compute δ_i = Π_{occupied bands} ξ at each TRIM
4. Compute (-1)^Z₂ = Π_{i=1}^{4} δ_i
5. Z₂ = 0 if product = +1, Z₂ = 1 if product = -1

**Step 2.2: Fix wire mode — gap-sweep method**

For wire mode (confinement=2) with BHZ parameters, use the gap-sweep method:
- Sweep M from -15 to +15 meV
- For each M, compute eigenvalues and track gap closing at E=0
- Z₂ = 1 if gap closes and reopens with odd number of crossings

**Step 2.3: BHZ Z₂ for HgTe/CdTe**

The literature (Bernevig Hughes Zhang 2006) gives:
- Z₂ = 0 for d < d_c ≈ 6.3 nm (normal band order, M>0)
- Z₂ = 1 for d > d_c (band-inverted, M<0)

Test configs: d=58Å (trivial, M=+10) and d=70Å (topological, M=-10)

### Files to modify
- `src/physics/topological_analysis.f90`:
  - Complete `compute_z2_fukane` function (currently stub at lines 169-176)
  - Add `compute_z2_gap_sweep` for wire mode
- `src/apps/main_topology.f90`:
  - Wire mode uses gap-sweep instead of gap-counting
  - QW mode uses Fu-Kane when `z2_method: fukane`

### Verification
```
trivial (M=+10, d=58Å): Z2=0
topological (M=-10, d=70Å): Z2=1
test_bhz_z2.sh: both pass
```

### What to produce (Phase 2 figures after fix)

- `bhz_z2_phase_diagram.png`: Z₂ vs. M — shows step at M=0
- `bhz_edge_localization.png`: edge localization length ξ vs. M for topological phase

---

## Phase 3: Implement Peierls Substitution — Landau Levels

### Root Cause

`magnetic_field.f90:add_peierls_coo` (lines 53-69) is a **stub**. Landau gauge substitution `k → k - eA/ℏ` is not implemented. Without Peierls, no Landau levels.

### Implementation

For wire along z with perpendicular B = Bx x̂, Landau gauge: **A = (0, 0, Bx·y)**

Algorithm:
1. At each grid point i with position y(i), compute Peierls phase:
   `phi_peierls(i) = exp(-i · e·Bx·y(i)·dz / ℏ)`
2. For kz-coupling COO entries, multiply by ratio of phase factors between adjacent sites
3. This correctly implements Peierls substitution in tight-binding basis

**Simpler approach:** Make `kz` a position-dependent array `kz_grid(:)` passed to the Hamiltonian builder. `ZB8bandGeneralized` takes `kz` as scalar — change to array with Peierls phase built in.

### Landau level test

Create `tests/regression/configs/landau_InAs.cfg`:
```
material: InAs
confinement: 0 (bulk)
ExternalField: 5.0 T (perpendicular)
```

### Verification

For InAs at B=5T, m* = 0.026 m₀:
```
ℏωc = e·ℏ·B/m* = 22.26 meV
E₀ = 11.13 meV (n=0)
E₁ = 33.39 meV (n=1)
E₂ = 55.65 meV (n=2)
```

Figure: `landau_levels_inas_b5t.png` — E_n vs. n with analytical line overlay

---

## Phase 4: BdG — Fix Config Parser + Figures

### Config Parsing Bug

`input_parser.f90:723-745` only reads `mu`, `delta_0`, `gauge`. The `g_factor: 2.0` line in `topology_rashba_phase.cfg` is silently ignored (default is also 2.0, so it coincidentally works).

### Fix
Add g_factor reading to BdG block parser.

### BdG Physics

Rashba wire + s-wave pairing + Zeeman → Majorana at B > B_crit:
```
B_crit = √(μ² + Δ²) / (g·μ_B)
```
With μ=0.5 meV, Δ=0.3 meV, g=2: **B_crit ≈ 5.0 T**

### What to produce

Figure: `rashba_majorana_phase_diagram.png`
- Sweep B from 0 to 10 T
- Plot min gap vs. B
- Gap closes at B_crit ≈ 5 T
- Overlay analytical prediction

---

## Phase 5: Update Lecture 13

### New Structure (pedagogical sequence)

**13.1 Overview** — implementation status, flow diagram

**13.2 Chern Number (QHE)** — simplest case
- QWZ model + FHS lattice gauge
- Figure: convergence + phase diagram
- Status: **VERIFIED ✓**

**13.3 Berry Curvature** — foundation for Chern
- Figure: Ω(kx,ky) heatmap

**13.4 BHZ Z₂ Invariant (QSHE)** — next complexity
- Fu-Kane parity + gap-sweep for wire
- Figures: phase diagram + edge localization
- Status: **FIXED → VERIFIED ✓**

**13.5 Landau Levels** — adds magnetic field
- E_n = ℏωc(n+½) for InAs at B=5T
- Figure: fan diagram
- Status: **VERIFIED ✓**

**13.6 BdG / Majorana Modes** — most complex
- Gap closing at B_crit
- Figures: phase diagram + Majorana wavefunction
- Status: **VERIFIED ✓**

**13.7 LDOS** — spectral verification
- Figure: LDOS(E) at wire end showing Majorana peak

**13.8 References**

---

## Summary of Deliverables

| Phase | Task | Files | Verification |
|---|---|---|---|
| 1 | QWZ Chern figures | `scripts/verify_qwz_chern.py`, figures/ | Figures prove C=+1,-1,0 |
| 2 | Fix BHZ Z₂ | topological_analysis.f90, main_topology.f90 | test_bhz_z2.sh passes |
| 3 | Peierls + Landau | magnetic_field.f90, hamiltonian_wire.f90 | Landau levels correct E_n |
| 3 | Landau figure | landau_InAs.cfg, figures/ | Figure: E₀=11.13, E₁=33.39 |
| 4 | BdG fix + figure | input_parser.f90, sweep script | Figure: min_gap vs. B |
| 5 | Update lecture 13 | docs/lecture/13-topological-superconductivity.md | All figures embedded |

---

## Critical Files

- `src/physics/topological_analysis.f90:139-176` — Z₂ gap and Fu-Kane functions
- `src/physics/topological_analysis.f90:437-520` — `build_bhz_wire_hamiltonian` (COO assembly)
- `src/physics/magnetic_field.f90:53-69` — Peierls stub
- `src/physics/bdg_hamiltonian.f90` — BdG assembly
- `src/io/input_parser.f90:723-745` — BdG block parser
- `docs/lecture/13-topological-superconductivity.md` — existing lecture to update
