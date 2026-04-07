# Design: Bulk Chapter Overhaul

Date: 2026-04-06

## Problem Statement

Five issues with `docs/lecture/01-bulk-band-structure.md` and supporting code:

1. **CRITICAL: `parts.dat` is overwritten** — only the last k-point's band decomposition survives, showing mixed states instead of pure Gamma-point character
2. **LaTeX formulas don't render** in VS Code markdown preview
3. **Missing directional sweeps** — no [110], [101], [011] directions to show valence band warping
4. **No bulk strain** — Pikus-Bir strain exists for QW/wire but not bulk mode
5. **Documentation needs rewriting** to cover all new features with nextnano tutorial validation

## Task 1: Fix `parts.dat` Overwrite Bug

### Root Cause

In `src/io/outputFunctions.f90:84`, `writeEigenfunctions` opens `parts.dat` with
`status='replace'` every call. The main program (`src/apps/main.f90:645-650`) calls
it at k=1, k=mid, k=end — only the LAST call's data survives.

At k=0 (Gamma), eigenvectors ARE pure one-hot vectors:
- State 1: 100% basis 5 (SO↑), State 2: 100% basis 6 (SO↓)
- States 3-6: 100% valence (HH/LH)
- States 7-8: 100% basis 7/8 (CB)

But the saved file shows k=0.1 data where bands mix.

### Fix

Change `parts.dat` to **multi-block gnuplot format** with k-point headers:

```fortran
! In writeEigenfunctions, for bulk:
if (is_bulk) then
  if (k == 1) then
    open(unit=iounit2, file=OUTPUT_DIR//'/parts.dat', status='replace', action='write')
  else
    open(unit=iounit2, file=OUTPUT_DIR//'/parts.dat', position='append', action='write')
  end if
  write(iounit2, '("# k = ", g14.6)') k_value  ! gnuplot index separator
  ! ... write parts rows ...
  close(iounit2)
end if
```

The gnuplot format uses blank lines + comment headers between blocks.
`gnuplot> plot 'parts.dat' index 0 ...` selects the Gamma-point block.

**Files modified:**
- `src/io/outputFunctions.f90` — writeEigenfunctions subroutine
- `scripts/plotting/generate_all_figures.py` — `fig_bulk_gaas_parts` to read first block

### Validation

Run bulk GaAs, verify parts.dat block 1 (k=0):
- Rows 1-2: 100% in columns 5,6 (SO)
- Rows 3-6: 100% in columns 1-4 (valence)
- Rows 7-8: 100% in columns 7,8 (CB)

## Task 2: VS Code LaTeX Rendering

Add `.vscode/settings.json` with math preview enabled:

```json
{
  "markdown.math.enabled": true
}
```

This enables built-in KaTeX rendering in VS Code 1.75+. No extensions needed.
The existing `$...$` and `$$...$$` syntax in docs will render immediately.

**Files created:**
- `.vscode/settings.json`

## Task 3: Diagonal Direction Sweeps (kxky, kxkz, kykz)

### Current State

`src/apps/main.f90:62-78` supports `kx`, `ky`, `kz`, `k0` — single-axis sweeps only.

### Implementation

Add three new waveVector modes in the select block:

```fortran
case ("kxky")
  ! [110] direction: kx = ky = k, kz = 0
  do i = 1, cfg%waveVectorStep
    k = (i-1) * cfg%waveVectorMax / (cfg%waveVectorStep - 1)
    smallk(i)%kx = k
    smallk(i)%ky = k
    smallk(i)%kz = 0.0_dp
  end do

case ("kxkz")
  ! [101] direction: kx = kz = k, ky = 0
  do i = 1, cfg%waveVectorStep
    k = (i-1) * cfg%waveVectorMax / (cfg%waveVectorStep - 1)
    smallk(i)%kx = k
    smallk(i)%kz = k
    smallk(i)%ky = 0.0_dp
  end do

case ("kykz")
  ! [011] direction: ky = kz = k, kx = 0
  do i = 1, cfg%waveVectorStep
    k = (i-1) * cfg%waveVectorMax / (cfg%waveVectorStep - 1)
    smallk(i)%ky = k
    smallk(i)%kz = k
    smallk(i)%kx = 0.0_dp
  end do
```

The eigenvalue output format already stores `|k|` (magnitude), so no output changes needed.
The plotting script should label the direction axis as `|k| [110]` etc.

### Physics: Valence Band Warping

Along [100]: HH effective mass is governed by γ₁ - 2γ₂
Along [110]: HH effective mass mixes with γ₃ → different curvature

The difference between [100] and [110] valence dispersion directly reveals
the cubic anisotropy parameterized by γ₂ ≠ γ₃. For GaAs: γ₂ = 2.06, γ₃ = 2.93.

### Config Files

Add regression configs:
- `tests/regression/configs/bulk_gaas_kxky.cfg`
- `tests/regression/configs/bulk_gaas_kxkz.cfg`

### Figures

New figures in `generate_all_figures.py`:
- `bulk_gaas_bands_110.png` — E(k) along [110] with comparison to [100]
- `bulk_gaas_warping.png` — overlay [100], [110], [111] dispersions (HH/LH only)

**Files modified:**
- `src/apps/main.f90` — add three cases to waveVector select
- `scripts/plotting/generate_all_figures.py` — new figure functions

## Task 4: Bulk Strain (Pikus-Bir)

### Current State

Strain is fully implemented for QW (`compute_strain_qw`) and wire (`compute_strain_wire`)
in `src/physics/strain_solver.f90`. The `apply_pikus_bir` subroutine applies deformation
potentials to the profile array. ALL material parameters are populated (C11, C12, C44, ac,
av, b_dp, d_dp, a0).

But bulk mode (ndim=0) returns immediately from `compute_strain` — no strain applied.

### Implementation: Uniform Biaxial Strain

For bulk, add strain directly to the Hamiltonian diagonal in `ZB8bandBulk`.
No need for the spatial strain solver — bulk is spatially uniform.

Add a new input parameter `strainSubstrate` (lattice constant in Angstrom):
```
strainSubstrate: 5.869   ! InP substrate lattice constant
```

When set (nonzero), compute uniform biaxial strain:
```
eps_parallel = (a_substrate - a_film) / a_film
eps_perp = -2 * C12/C11 * eps_parallel
Tr(eps) = 2*eps_parallel + eps_perp
```

Apply Pikus-Bir shifts to the 8x8 Hamiltonian diagonal:
```
delta_Ec = ac * Tr(eps)                        ! CB shift
delta_HH = av * Tr(eps) + b/2 * (eps_perp - eps_parallel)  ! HH shift
delta_LH = av * Tr(eps) - b/2 * (eps_perp - eps_parallel)  ! LH shift
delta_SO = av * Tr(eps)                                    ! SO shift (no shear term)
```

Note: `av` uses positive sign convention (already in code). The deformation
potential `b_dp` is the Bir shear parameter.

### Hamiltonian Modification

In `ZB8bandBulk`, after building the standard Hamiltonian, add strain shifts:

```fortran
! Strain shifts (if strainSubstrate > 0)
if (params(1)%strainSubstrate > 0.0_dp) then
  a_film = params(1)%a0
  eps_par = (params(1)%strainSubstrate - a_film) / a_film
  eps_per = -2.0_dp * params(1)%C12 / params(1)%C11 * eps_par
  Tr_eps = 2.0_dp * eps_par + eps_per

  ! CB (bands 7,8): ac * Tr(eps)
  HT(7,7) = HT(7,7) + params(1)%ac * Tr_eps
  HT(8,8) = HT(8,8) + params(1)%ac * Tr_eps

  ! HH (bands 1,4): av*Tr + b/2*(eps_perp - eps_par)
  HH_shift = params(1)%av * Tr_eps + params(1)%b_dp * 0.5_dp * (eps_per - eps_par)
  HT(1,1) = HT(1,1) + HH_shift
  HT(4,4) = HT(4,4) + HH_shift

  ! LH (bands 2,3): av*Tr - b/2*(eps_perp - eps_par)
  LH_shift = params(1)%av * Tr_eps - params(1)%b_dp * 0.5_dp * (eps_per - eps_par)
  HT(2,2) = HT(2,2) + LH_shift
  HT(3,3) = HT(3,3) + LH_shift

  ! SO (bands 5,6): av*Tr
  SO_shift = params(1)%av * Tr_eps
  HT(5,5) = HT(5,5) + SO_shift
  HT(6,6) = HT(6,6) + SO_shift
end if
```

### Config

```ini
waveVector: kx
waveVectorMax: 0.1
waveVectorStep: 11
confinement:  0
FDstep: 1
FDorder: 2
numLayers:  1
material1: GaAs
numcb: 2
numvb: 6
ExternalField: 0  EF
EFParams: 0.0005
strainSubstrate: 5.869
```

### Validation (nextnano tutorial)

GaAs pseudomorphically strained to InP (a₀ = 5.869 Å):
- GaAs a₀ = 5.653 Å → eps_∥ = (5.869 - 5.653)/5.653 = +3.82% (tensile)
- eps_⊥ = -2 × 566/1221 × 0.0382 = -0.0354
- Tr(eps) = 2(0.0382) + (-0.0354) = 0.041
- δE_C = -7.17 × 0.041 = -0.294 eV
- δE_HH = 1.16 × 0.041 + (-2.0)/2 × (-0.0354 - 0.0382) = 0.048 + 0.074 = 0.121 eV
- δE_LH = 1.16 × 0.041 - (-2.0)/2 × (-0.0354 - 0.0382) = 0.048 - 0.074 = -0.027 eV

Result: HH/LH splitting at Gamma = 0.148 eV. CB shifts down by 0.294 eV.

### Figures

- `bulk_gaas_strained_bands.png` — E(k) for strained GaAs showing HH/LH splitting
- `bulk_gaas_strain_comparison.png` — overlay unstrained vs strained

**Files modified:**
- `src/physics/hamiltonianConstructor.f90` — add strain to ZB8bandBulk
- `src/core/defs.f90` — add strainSubstrate to paramStruct
- `src/core/parameters.f90` — initialize strainSubstrate = 0
- `src/io/input_parser.f90` — parse strainSubstrate
- `scripts/plotting/generate_all_figures.py` — new figure functions

## Task 5: Config Files

New regression configs:
- `tests/regression/configs/bulk_gaas_kxky.cfg` — [110] sweep
- `tests/regression/configs/bulk_gaas_strained.cfg` — strained to InP

## Task 6: Figure Updates

Update `scripts/plotting/generate_all_figures.py`:

1. **Fix `fig_bulk_gaas_parts`** — read only the first block (k=0) from multi-block parts.dat
2. **Add `fig_bulk_gaas_bands_110`** — E(k) along [110]
3. **Add `fig_bulk_gaas_warping`** — overlay [100] vs [110] valence bands
4. **Add `fig_bulk_gaas_strained_bands`** — strained GaAs E(k)
5. **Add `fig_bulk_gaas_strain_comparison`** — unstrained vs strained overlay

## Task 7: Documentation Rewrite

Rewrite `docs/lecture/01-bulk-band-structure.md`:

- Section 1 (Theory): Keep existing content, fix any LaTeX issues
- Section 2 (In the Code): Update module map, add strain params
- Section 3 (Computed Examples):
  - 3.1 Unstrained GaAs along [100] (existing)
  - 3.2 **NEW** GaAs along [110] — show valence band warping
  - 3.3 **NEW** Directional comparison — overlay [100] vs [110] vs [011]
  - 3.4 **NEW** Strained GaAs (GaAs on InP) — HH/LH splitting
  - 3.5 InAs comparison (existing, keep)
  - 3.6 Effective mass extraction (existing)
  - 3.7 Band decomposition at Gamma (fix with correct data)
- Section 4 (nextnano validation): Keep, add strained comparison
- Section 5 (Discussion): Update with warping and strain discussion

## Implementation Order

1. Fix parts.dat overwrite (Task 1) — immediate, unblocks correct figure
2. Add .vscode/settings.json (Task 2) — trivial
3. Add kxky/kxkz/kykz sweeps (Task 3) — small Fortran change
4. Add bulk strain (Task 4) — medium Fortran change
5. Update Python figures (Task 6) — depends on 1, 3, 4
6. Add regression configs (Task 5) — depends on 3, 4
7. Rewrite docs (Task 7) — depends on all above

## Constraints

- **NEVER** change basis ordering (CLAUDE.md rule)
- **NEVER** modify material parameters without reference verification
- All new code follows existing patterns (select case for waveVector, etc.)
- Strain in bulk uses same sign conventions as QW/wire (positive av convention)
