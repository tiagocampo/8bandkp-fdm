# Ch01-03 Review & Fix Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix all 9 findings from the review of lecture chapters 01, 02, and 03 — stale data, consistency issues, content gaps, and minor fixes.

**Architecture:** Three chapters get a consistency pass. Ch03 gets a major expansion (new type-I example + k-evolution data). Ch02 Sec 2.3 gets rewritten to match the staggered-grid FD fix. All three chapters get normalized section numbering. Simulations are re-run first to provide fresh data.

**Tech Stack:** Fortran 90 codebase (bandStructure executable), gnuplot, Python for data extraction, markdown editing.

---

### Task 1: Re-run broken-gap simulation with fine grid (F1 data source)

This task produces the fresh simulation data that Tasks 5, 6, and 7 depend on.

**Files:**
- Create: `input.cfg` (temporary, do not commit)
- Read: `output/eigenvalues.dat`, `output/parts.dat`, `output/eigenfunctions_k_00001_ev_*.dat`

**Step 1: Create fine-grid broken-gap input.cfg**

Write `input.cfg` with the broken-gap config upgraded to FDstep=401, FDorder=4:

```
waveVector: kx
waveVectorMax: 0.1
waveVectorStep: 51
confinement:  1
FDstep: 401
FDorder: 4
numLayers:  3
material1: AlSbW -250  250 0
material2: GaSbW -135  135 0.2414
material3: InAsW  -35   35 -0.0914
numcb: 10
numvb: 10
ExternalField: 0  EF
EFParams: 0.0005
```

Note: Request fewer states (10+10 instead of 32+32) since we only need the states near the gap for Ch03. Use 51 k-steps to get k-evolution data (F6). Do NOT commit this file.

**Step 2: Build and run**

```bash
cmake --build build
./build/src/bandStructure
```

Expected: completes in ~1-2 minutes. Output in `output/`.

**Step 3: Extract key data points**

From `output/eigenvalues.dat`:
- Row 1 (k=0): extract the 10 lowest eigenvalues for VB and 10 highest for CB
- Identify CB1 energy and VB-1 energy
- Verify effective gap matches Ch02's ~65 meV (should be very close with FDstep=401)

From `output/parts.dat`:
- Row for the CB1 state: extract all 8 parts values
- Row for the VB-1 state: extract all 8 parts values

From `output/eigenfunctions_k_00001_ev_*.dat`:
- Read the CB1 ground-state eigenfunction file
- Extract values at z=-250, -150, -40, -35, -10, 0, +10, +35, +150, +250

From `output/eigenvalues.dat`:
- Extract CB1 energy at k=0, k=0.05, k=0.10 for k-evolution table (F6)

Record all extracted values for use in Tasks 5, 6, 7.

**Step 4: Verify data consistency**

Check that CB1 energy from this run matches Ch02 Section B.3 values (+0.0319 eV) within ~1 meV. If values differ significantly, investigate before proceeding — the fine-grid config should produce consistent results.

---

### Task 2: Re-run type-I GaAs/AlGaAs simulation for Ch03 expansion (F5 data source)

This task produces data for the new type-I wavefunction example in Ch03.

**Files:**
- Create: `input.cfg` (temporary, do not commit)
- Read: `output/eigenvalues.dat`, `output/parts.dat`, `output/eigenfunctions_k_00001_ev_*.dat`

**Step 1: Create type-I config for wavefunction extraction**

Write `input.cfg`:

```
waveVector: kx
waveVectorMax: 0.1
waveVectorStep: 21
confinement:  1
FDstep: 401
FDorder: 4
numLayers:  2
material1: Al30Ga70As -200 200 0
material2: GaAs -50 50 0
numcb: 4
numvb: 8
ExternalField: 0  EF
EFParams: 0.0
```

This matches the Ch02 A.3 config exactly. Do NOT commit.

**Step 2: Build and run**

```bash
./build/src/bandStructure
```

Expected: completes in seconds. Output in `output/`.

**Step 3: Extract key data points**

From `output/parts.dat` at k=0:
- Extract parts for CB1 (state with highest energy), CB2
- Extract parts for VB-1 (HH1 ground state), VB-5 (LH1)

From `output/eigenfunctions_k_00001_ev_*.dat`:
- Extract CB1 spatial profile at z=-200, -100, -60, -50, -25, 0, +25, +50, +60, +100, +200
- Extract HH1 spatial profile at same positions

From `output/eigenvalues.dat`:
- Verify CB1 = +0.7613 eV and VB-8 = -0.8262 eV match Ch02 A.3

Record all values for Task 5.

---

### Task 3: Rewrite Ch02 Section 2.3 with staggered-grid description (F2)

**Files:**
- Modify: `docs/lecture/02-quantum-well.md:319-345` (Section 2.3)

**Step 1: Rewrite the section**

Replace lines 319-345 with a new Section 2.3 that describes the staggered-grid conservative form. The new content should:

1. Introduce the same 7-point example (barrier/well/barrier, points 1-2, 3-5, 6-7)
2. Explain the **conservative form**: $d/dz[g(z) \cdot d/dz] \approx D_{\text{outer}} \cdot \text{diag}(g_{1/2}) \cdot D_{\text{inner}}$
3. Define $D_{\text{inner}}$ as the (N-1)×N forward half-point difference: `[-1/dz, +1/dz]` per row
4. Define $D_{\text{outer}} = -D_{\text{inner}}^T$ (ensures self-adjointness)
5. Show the half-point coefficient $g_{j+1/2}$:
   - FDorder=2: $(g_j + g_{j+1})/2$ (simple average)
   - FDorder>=4: $(-g_{j-1} + 9g_j + 9g_{j+1} - g_{j+2})/16$ (4th-order Lagrange)
6. Show the worked example with illustrative g values:

```
g:  4.00   4.00   2.88   2.88   2.88   4.00   4.00
g_half (2nd-order):  4.00  3.44  2.88  2.88  3.44  4.00
                    (avg of adjacent g values)

kpterms(:,:,7) = D_outer @ diag(g_half) @ D_inner

Row 3 (well interior):
  g_half(2)*D_inner(row 2) + g_half(3)*D_inner(row 3)
  = 3.44*[-1/dz, +1/dz, 0, ...] + 2.88*[0, -1/dz, +1/dz, ...]
  → diagonal element: (3.44 + 2.88)/dz^2 = 6.32/dz^2
  → off-diag left: -3.44/dz^2, off-diag right: -2.88/dz^2

Row 2 (barrier/well interface):
  g_half(1)*D_inner(row 1) + g_half(2)*D_inner(row 2)
  = 4.00*[-1/dz, +1/dz, 0, ...] + 3.44*[0, -1/dz, +1/dz, ...]
  → diagonal: (4.00 + 3.44)/dz^2 = 7.44/dz^2
  → off-diag left: -4.00/dz^2, off-diag right: -3.44/dz^2
```

7. Point to Section 4.6 for further discussion and to Ch09 for the full FD stencil details

**Step 2: Update Section 4.6 for consistency**

Read `docs/lecture/02-quantum-well.md:940-955`. Update to mention the staggered-grid approach by name and reference the new Sec 2.3. The current text says "midpoint averaging" and "applyVariableCoeff" — update to say "staggered-grid conservative form" and "applyVariableCoeffStaggered" to match the actual code.

**Step 3: Commit**

```bash
git add docs/lecture/02-quantum-well.md
git commit -m "docs: rewrite Ch02 Sec 2.3 with staggered-grid variable-coeff description"
```

---

### Task 4: Add energy reference note to Ch02 (F3)

**Files:**
- Modify: `docs/lecture/02-quantum-well.md:169` (start of Section 1.5)

**Step 1: Insert energy convention note**

After the opening paragraph of Section 1.5 (line 171), add a note:

```
> **Energy convention change.** In bulk mode (Chapter 01), eigenvalues are
> referenced to an internal energy zero where $E_V = 0$. In quantum well mode,
> the material database band offsets ($E_V$, $E_C$) are added to the diagonal
> of the Hamiltonian, shifting all eigenvalues to an **absolute energy scale**.
> For GaAs, this means the valence band edge appears at $E_V = -0.800$ eV and
> the conduction band at $E_C = +0.719$ eV, rather than at 0 and $E_g = 1.519$
> eV respectively.
```

**Step 2: Commit**

```bash
git add docs/lecture/02-quantum-well.md
git commit -m "docs: add energy convention note at bulk-to-QW transition in Ch02"
```

---

### Task 5: Fix Ch01 Figure 7 caption (F8)

**Files:**
- Modify: `docs/lecture/01-bulk-band-structure.md:654-656` (Figure 7 caption)

**Step 1: Update caption**

Find the Figure 7 caption (around line 654-656):

```
*Figure 7: Evolution of band character from pure states at $\Gamma$ to mixed
states at finite $k$, along [110]. The top-left panels (states 1-4) show the
valence and SO states mixing as $k$ increases. States 7-8 (CB) remain nearly
pure at small $k$ but acquire increasing valence admixture from the $P$-coupling.*
```

Replace with:

```
*Figure 7: Evolution of band character from pure states at $\Gamma$ to mixed
states at finite $k$, along [110]. Eigenstates are ordered by ascending energy:
states 1--2 are split-off, states 3--6 are valence (LH then HH), states 7--8
are conduction. The top-left panels show the SO and LH states mixing as $k$
increases. States 7--8 (CB) remain nearly pure at small $k$ but acquire
increasing valence admixture from the $P$-coupling.*
```

**Step 2: Commit**

```bash
git add docs/lecture/01-bulk-band-structure.md
git commit -m "docs: clarify Ch01 Figure 7 eigenvalue ordering in caption"
```

---

### Task 6: Normalize section numbering across Ch01, Ch02, Ch03 (F4)

**Files:**
- Modify: `docs/lecture/01-bulk-band-structure.md` — headings
- Modify: `docs/lecture/03-wavefunctions.md` — headings

**Current state:**
- Ch01: `## 1. Theory` / `### 1.1 ...` / `#### Basis Function Table` — clean, no change needed
- Ch02: `## 1. Theory` / `### 1.1 ...` / `### Example A` / `#### A.1 ...` — clean, no change needed
- Ch03: `## 3.1 From Eigenvalues...` — **problem**: uses chapter-internal numbering starting at 3.x instead of 1.x like the other chapters

**Step 1: Renumber Ch03 sections**

Ch03 currently uses `## 3.N` through `## 3.9` (chapter-prefixed). Change to `## 1.` through `## 5.` to match Ch01 and Ch02's convention of starting each chapter at Section 1:

- `## 3.1 From Eigenvalues...` → `## 1. From Eigenvalues...`
- `## 3.2 The Block Structure...` → `## 2. The Block Structure...`
- `### 3.2.1 ...` → `### 2.1 ...`
- `### 3.2.2 ...` → `### 2.2 ...`
- `### 3.2.3 ...` → `### 2.3 ...`
- `## 3.3 Band-Resolved...` → `## 3. Band-Resolved...`
- `### 3.3.1 ...` → `### 3.1 ...`
- `### 3.3.2 ...` → `### 3.2 ...`
- `## 3.4 Band Character...` → `## 4. Band Character...`
- `### 3.4.1 ...` → `### 4.1 ...`
- `### 3.4.2 ...` → `### 4.2 ...`
- `## 3.5 Heavy-Hole...` → `## 5. Heavy-Hole...`
- `### 3.5.1 ...` → `### 5.1 ...`
- `### 3.5.2 ...` → `### 5.2 ...`
- `## 3.6 Implementation...` → `## 6. Implementation...`
- `### 3.6.1 ...` → `### 6.1 ...`
- `### 3.6.2 ...` → `### 6.2 ...`
- `### 3.6.3 ...` → `### 6.3 ...`
- `## 3.7 Computed Example...` → `## 7. Computed Example...` (will become 7a/7b after Task 7 adds type-I)
- `### 3.7.1 ...` → `### 7.1 ...` through `### 7.6 ...`
- `## 3.8 Discussion` → `## 8. Discussion`
- `### 3.8.1 ...` → `### 8.1 ...`
- `### 3.8.2 ...` → `### 8.2 ...`
- `### 3.8.3 ...` → `### 8.3 ...`
- `## 3.9 Summary` → `## 9. Summary`

Also update any internal cross-references (e.g., "Section 3.4" → "Section 4").

**Step 2: Verify Ch01 and Ch02 numbering is consistent**

Ch01: Uses `## N.` / `### N.M` / `#### subheading` — already consistent. No changes needed.

Ch02: Uses `## N.` / `### N.M` / `#### A.M` for examples — already consistent. No changes needed.

**Step 3: Commit**

```bash
git add docs/lecture/03-wavefunctions.md
git commit -m "docs: normalize Ch03 section numbering to match Ch01/Ch02 convention"
```

---

### Task 7: Update Ch03 with fresh type-III data (F1, F9)

This task uses simulation output from Task 1.

**Files:**
- Modify: `docs/lecture/03-wavefunctions.md` — Sections 7.1 through 7.6 (renumbered from 3.7.x)

**Step 1: Update Section 7.1 (input configuration)**

Update the config snippet to show FDstep=401, FDorder=4 (matching the fine-grid run from Task 1). Add a note: "This configuration uses a finer grid (FDstep=401, FDorder=4) than the reference regression test (FDstep=101, FDorder=2) for improved accuracy."

**Step 2: Update Section 7.2 (eigenfunction file excerpt)**

Replace the eigenfunction file excerpt (lines 249-254) with fresh data from Task 1's output. The new data will have 401 rows instead of 101 (different z-spacing). Show the first 5 lines as before.

**Step 3: Update Section 7.3 (CB1 wavefunction table)**

Replace the wavefunction table (lines 277-288) with fresh data from Task 1. Use actual values from the eigenfunction file at positions near z=-250, -150, -40, -35, -10, 0, +10, +35, +150, +250. Note: z-spacing changes from 5A to 1.25A with FDstep=401, so exact positions may differ slightly.

**Step 4: Update Section 7.5 (parts table)**

Replace the parts table (lines 315-325) with fresh data from Task 1's `parts.dat`. Use 4 decimal places for all values. Replace mechanical `0.000` with the actual small values (e.g., `0.0001`). Keep the bold formatting on key rows.

**Step 5: Update eigenvalue references in prose**

Update all inline eigenvalue references (e.g., "eigenvalue $E_{33} = +0.0205$ eV" on line 247) with the new values from Task 1. The state indices may change with the different numcb/numvb (10+10 vs 32+32) — identify the correct state index for CB1 in the new output.

**Step 6: Commit**

```bash
git add docs/lecture/03-wavefunctions.md
git commit -m "docs: update Ch03 type-III data from FDstep=401 FDorder=4 run"
```

---

### Task 8: Add type-I wavefunction example to Ch03 (F5)

This task uses simulation output from Task 2. It adds a new computed example section to Ch03.

**Files:**
- Modify: `docs/lecture/03-wavefunctions.md` — insert new section after the type-III example (after current Section 7.6)

**Step 1: Add new section header**

After Section 7.6 (Spatial profiles), insert:

```markdown
## 7b. Computed Example: GaAs/AlGaAs Type-I Quantum Well

For comparison with the broken-gap system above, we examine the wavefunctions
of the GaAs/Al$_{0.3}$Ga$_{0.7}$As quantum well from Chapter 02, Example A.

### 7b.1 Configuration

[Config snippet from Task 2 — same as Ch02 A.1 config]

### 7b.2 CB1 ground-state wavefunction

[Spatial profile table from Task 2 data]

Key difference from the broken-gap CB1: the type-I CB1 has $P_7 + P_8 > 0.99$ —
essentially pure conduction-band character. There is no significant valence-band
admixture because the band gap (1.519 eV) is much larger than the k.p coupling
strength. The wavefunction is a half-sine envelope localized in the GaAs layer,
with exponential tails penetrating into the AlGaAs barriers.

### 7b.3 HH1 ground-state wavefunction

[Spatial profile table from Task 2 data for the deepest VB state]

The HH1 state is confined in the GaAs layer with $P_1 + P_4 > 0.99$ — essentially
pure heavy-hole character at $k_\parallel = 0$. The confinement energy of
~26 meV below the GaAs VB edge is consistent with the 159 meV VB offset and the
100 A well width.

### 7b.4 Band-resolved parts

[Parts table from Task 2 — show CB1, CB2, HH1, HH2, LH1, LH2]

| State | $E$ (eV) | $P_{HH}$ | $P_{LH}$ | $P_{SO}$ | $P_{CB}$ | Character |
| [data from Task 2] |

Comparison with the broken-gap system:

| Property | Type-I (GaAs/AlGaAs) | Type-III (AlSbW/GaSbW/InAsW) |
|---|---|---|
| CB1 purity ($P_{CB}$) | >0.99 | ~0.66 |
| CB1 VB admixture | <1% | ~34% (LH2) |
| Wavefunction localization | GaAs layer | InAs layer |
| HH1 localization | GaAs layer | GaSb layer |
```

**Step 2: Commit**

```bash
git add docs/lecture/03-wavefunctions.md
git commit -m "docs: add type-I GaAs/AlGaAs wavefunction example to Ch03"
```

---

### Task 9: Add k-evolution data to Ch03 (F6)

This task uses simulation output from Task 1 (multiple k-points).

**Files:**
- Modify: `docs/lecture/03-wavefunctions.md` — insert into the Discussion section (Section 8) or add a new subsection to the type-III example

**Step 1: Add k-evolution subsection**

Add a new subsection (either within Section 7 or in Section 5 "Heavy-Hole vs Light-Hole Mixing") showing how the band character evolves with $k_\parallel$ for the broken-gap QW.

Use data from Task 1's 51 k-point run. Extract parts at k=0, k=0.03, k=0.05, k=0.08, k=0.10 for the CB1 state and the top VB state.

Create a table:

```
### 7.7 Band character evolution with $k_\parallel$

The integrated band parts change with in-plane wavevector as the off-diagonal
k.p couplings activate HH-LH mixing. For the broken-gap CB1 state:

| $k_\parallel$ (A$^{-1}$) | $P_{HH}$ | $P_{LH}$ | $P_{SO}$ | $P_{CB}$ |
|---|---|---|---|---|
| 0.00 | [from Task 1] | ... | ... | ... |
| 0.03 | ... | ... | ... | ... |
| 0.05 | ... | ... | ... | ... |
| 0.08 | ... | ... | ... | ... |
| 0.10 | ... | ... | ... | ... |

At $k=0$, the CB1 state has [~66% CB, ~33% LH] character. As $k_\parallel$
increases, [describe trend from data]. This contrasts with the type-I GaAs/AlGaAs
system where the CB1 state remains >99% conduction-band at all k values.
```

**Step 2: Commit**

```bash
git add docs/lecture/03-wavefunctions.md
git commit -m "docs: add k-evolution band character data to Ch03"
```

---

### Task 10: Update Section 4.6 in Ch02 for consistency (part of F2)

This is a follow-up to Task 3's Sec 2.3 rewrite. Ensures the discussion section matches.

**Files:**
- Modify: `docs/lecture/02-quantum-well.md:940-955`

**Step 1: Rewrite Section 4.6**

Replace the current text (which mentions "midpoint averaging" and `applyVariableCoeff`) with a description that matches the actual code:

```
### 4.6 Variable material parameters at interfaces

When the Luttinger parameters ($\gamma_1, \gamma_2, \gamma_3$) change abruptly
at a heterointerface, the product $\gamma(z) \cdot d^2/dz^2$ is not a Hermitian
operator. The code uses a **staggered-grid conservative form** (see Section 2.3
for a worked example):

$$d/dz[g(z) \cdot d/dz] \approx D_{\text{outer}} \cdot \text{diag}(g_{1/2}) \cdot D_{\text{inner}}$$

where $D_{\text{inner}}$ is the half-point forward difference, $D_{\text{outer}} = -D_{\text{inner}}^T$,
and $g_{1/2}$ is the coefficient interpolated to half-grid points. For FDorder=2,
this reduces to the standard 3-point stencil with $(g_j + g_{j+1})/2$ averaging
(implicit in the `dgemv`-based construction). For FDorder>=4, a 4th-order Lagrange
interpolation $(-g_{j-1} + 9g_j + 9g_{j+1} - g_{j+2})/16$ provides the half-point
values. This approach ensures the resulting operator is exactly self-adjoint and
correctly handles sharp interfaces at all FD orders. The Foreman renormalization
(disabled by default via `renormalization = .False.` in `defs.f90`) provides an
alternative treatment at the cost of additional complexity.
```

**Step 2: Commit**

```bash
git add docs/lecture/02-quantum-well.md
git commit -m "docs: update Ch02 Sec 4.6 to match staggered-grid implementation"
```

---

### Task 11: Final verification

**Step 1: Build and run all tests**

```bash
cmake --build build
ctest --test-dir build
```

Expected: all 15 tests pass.

**Step 2: Check for broken cross-references**

Grep all three chapter files for:
- References to "Section N.M" — verify they point to the correct (potentially renumbered) sections
- `![Figure](../figures/...)` — verify all referenced figures exist in `docs/figures/`

**Step 3: Verify eigenvalue consistency across chapters**

Check that eigenvalues mentioned in Ch02 and Ch03 are consistent:
- Ch02 A.3 (type-I): CB1=+0.7613 eV
- Ch02 B.3 (type-III): CB1=+0.0319 eV
- Ch03 type-I example: matches Ch02 A.3
- Ch03 type-III example: matches Ch02 B.3 (with fine-grid tolerance)

**Step 4: Clean up**

Remove the temporary `input.cfg` (do not commit). Ensure no simulation output files are staged.

---

## Summary of commits

1. `docs: rewrite Ch02 Sec 2.3 with staggered-grid variable-coeff description`
2. `docs: add energy convention note at bulk-to-QW transition in Ch02`
3. `docs: clarify Ch01 Figure 7 eigenvalue ordering in caption`
4. `docs: normalize Ch03 section numbering to match Ch01/Ch02 convention`
5. `docs: update Ch03 type-III data from FDstep=401 FDorder=4 run`
6. `docs: add type-I GaAs/AlGaAs wavefunction example to Ch03`
7. `docs: add k-evolution band character data to Ch03`
8. `docs: update Ch02 Sec 4.6 to match staggered-grid implementation`
