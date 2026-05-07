# Ch01-03 Review & Fix Design

**Date:** 2026-04-18
**Branch:** `feature/docs-overhaul`

## Scope

Fix all 9 findings from the review of chapters 01, 02, and 03 in a single pass.

## Findings

### Critical (accuracy/stale data)

**F1: Ch03 stale numerical data.** Ch03's computed example (Sec 3.7) uses the old `qw_alsbw_gasbw_inasw.cfg` with FDstep=101, FDorder=2. CB1=+0.021 eV vs Ch02's CB1=+0.0319 eV (FDstep=401, FDorder=4). Wavefunction tables (Sec 3.7.3), parts table (Sec 3.7.5), and eigenfunction file excerpts (Sec 3.7.2) all come from the old run.

**Fix:** Re-run with FDstep=401, FDorder=4 (same as Ch02's fine-grid run). Update all numerical tables and excerpts with new data.

### Significant (consistency)

**F2: Ch02 Sec 2.3 describes old buggy approach.** Lines 332-345 say "each row i gets the stencil weighted by g(i)" — the naive `diag(g) @ D2` row-scaling that was the FDorder>2 bug. This contradicts the corrected Sec 4.6 which describes the staggered-grid conservative form.

**Fix:** Rewrite Sec 2.3 to describe the staggered-grid approach: `D_outer @ diag(g_half) @ D_inner` where D_inner is the half-point forward difference and D_outer = -D_inner^T. Include a worked example for the 3-layer system with the corrected operator.

**F3: Energy reference convention not documented at transition.** Ch01 Sec 2.2 states bulk uses E_V=0, but the transition to absolute energies in Ch02 is not explained. A reader going from Ch01 to Ch02 would be confused by the sudden appearance of E_V = -0.8 eV.

**Fix:** Add a boxed note at the start of Ch02 Section 1.5 (band offsets) explaining: "In bulk mode (Ch01), eigenvalues are referenced to E_V=0. In quantum well mode, the band offsets from the material database (E_V, E_C) are added to the profile, shifting all eigenvalues to an absolute energy scale."

**F4: Section numbering format inconsistent.** Ch01/Ch03 use `### N.M`. Ch02 uses `#### A.M`/`#### B.M` for examples. Ch03 uses `### N.M.L` for subsections.

**Fix:** Adopt a consistent convention. Recommended:
- `## N.` for top-level sections (theory, code, examples, discussion)
- `### N.M` for subsections
- `### N.M.L` for sub-subsections only when needed
- Keep Ch02's `A.`/`B.` example labeling within the `### N.M` structure but as `#### A.M` under a `### 3. Computed Examples` parent. This matches current structure; normalize Ch01 and Ch03 to the same pattern.

### Content gaps

**F5: Ch03 missing type-I wavefunction example.** No GaAs/AlGaAs wavefunction data for comparison with the broken-gap system. Ch02 has both type-I and type-III; Ch03 only type-III.

**Fix:** Add a new section (3.7-equivalent) for the GaAs/AlGaAs QW wavefunctions. Include:
- CB1 ground-state spatial profile (localized in GaAs, no VB admixture)
- HH1 ground-state spatial profile
- Parts table showing P_7+P_8 ≈ 1.0 for CB states (contrast with broken-gap 33% LH admixture)
- Brief comparison with Ch03's type-III example

**F6: Ch03 no k-evolution data.** Sec 3.5 discusses HH/LH mixing at finite k. Sec 3.8.3 mentions "tracking band mixing with k" as a tip. But no figure or table shows this evolution.

**Fix:** Add a brief section or figure reference showing parts evolution for CB1 vs k_parallel. Can reference Ch01's Figure 7 (bulk parts vs k) as precedent. For QW, add a short table or figure showing P_HH and P_LH vs k for the first few VB states.

**F7: Ch03 too short.** ~400 lines vs ~800-970 for Ch01/Ch02. F5 and F6 will expand it to comparable length.

### Minor

**F8: Ch01 Figure 7 caption.** "states 1-4" refers to ascending eigenvalue ordering (SO, LH, HH), not band index. Could confuse readers.

**Fix:** Rephrase caption to "The top-left panels show the valence and SO eigenstates (ordered by ascending energy)..." or similar.

**F9: Ch03 parts table over-rounded.** Values like P_1=1.000 and 0.000 suggest exact values where they're approximate. Some 0.000 entries are likely just small.

**Fix:** After re-running (F1), use the actual output values with appropriate precision (4 decimal places). Replace mechanical 0.000 with "~0" or show 4 digits.

## Execution Plan

### Task 1: Re-run simulations for Ch03 data (F1)
- Run `qw_alsbw_gasbw_inasw.cfg` equivalent with FDstep=401, FDorder=4
- Extract eigenfunction data for CB1 (state 33 equivalent)
- Extract parts.dat at k=0
- Extract eigenfunction at a few k-points for k-evolution data (F6)

### Task 2: Rewrite Ch02 Sec 2.3 (F2)
- Replace diag(g)@D2 description with staggered-grid conservative form
- Include worked example with corrected operator for 3-layer system
- Ensure consistency with Sec 4.6

### Task 3: Add energy reference note (F3)
- Add a note box at Ch02 Sec 1.5 explaining the bulk-to-QW energy convention change

### Task 4: Normalize section numbering (F4)
- Audit all three chapters, adopt consistent numbering scheme
- Minimal changes — mainly ensuring Ch01 and Ch03 match Ch02's pattern

### Task 5: Expand Ch03 with type-I example (F5)
- Add new section with GaAs/AlGaAs wavefunction data
- Run the type-I config, extract wavefunctions and parts
- Write comparison prose

### Task 6: Add k-evolution data to Ch03 (F6)
- Use simulation data from Task 1 to show parts vs k_parallel
- Add table or figure reference

### Task 7: Fix minor issues (F8, F9)
- Fix Ch01 Figure 7 caption
- Fix Ch03 parts table precision after Task 1 re-run

## Dependencies

- Task 1 must complete before Tasks 5, 6, 7 (they need new simulation data)
- Task 2 is independent (pure writing)
- Tasks 3, 4 are independent (pure writing)
- Task 5 needs a type-I simulation run (separate from Task 1)

## Files to Modify

- `docs/lecture/01-bulk-band-structure.md` — F4 (numbering), F8 (figure caption)
- `docs/lecture/02-quantum-well.md` — F2 (Sec 2.3), F3 (energy note), F4 (numbering)
- `docs/lecture/03-wavefunctions.md` — F1 (stale data), F4 (numbering), F5 (type-I example), F6 (k-evolution), F7 (expansion), F9 (precision)
