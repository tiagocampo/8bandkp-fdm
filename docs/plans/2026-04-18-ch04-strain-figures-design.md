# Ch04 Strain Figures Design

**Date:** 2026-04-18
**Branch:** `feature/docs-overhaul`

## Scope

Add 7 figures to `docs/lecture/04-strain.md` — 2 existing (just add references), 2 new simulation-based, 3 new schematics.

## Figure List

| # | Filename | Type | Section | Source |
|---|---|---|---|---|
| 1 | `bulk_gaas_strained_bands.png` | Sim | 3.1 | Existing — add `![...]` reference |
| 2 | `bulk_gaas_strain_comparison.png` | Sim | 3.3 | Existing — add `![...]` reference |
| 3 | `qw_strained_band_edges.png` | Sim | 3.2 / 7.1 | New — AlSb/GaSb/InAs strained QW potential profile |
| 4 | `wire_strain_2d.png` | Sim | 6.6 | New — wire 2D eps_xx colormap |
| 5 | `strain_lattice_mismatch.png` | Schematic | 1.1 | New — lattice mismatch cartoon |
| 6 | `strain_biaxial_tensor.png` | Schematic | 2.4 | New — biaxial strain tensor diagram |
| 7 | `bir_pikus_band_shifts.png` | Schematic | 3.2 | New — BP energy level diagram (compressive + tensile) |

## Figure Details

### Figure 1: Strained bulk GaAs bands (existing)
- File: `docs/figures/bulk_gaas_strained_bands.png`
- Placement: Sec 3.1, after the Bir-Pikus Hamiltonian introduction
- Caption: "GaAs band structure under 0.65% biaxial strain on InP substrate. The HH/LH degeneracy at Gamma is lifted."

### Figure 2: Strain comparison (existing)
- File: `docs/figures/bulk_gaas_strain_comparison.png`
- Placement: Sec 3.3 (HH/LH splitting discussion)
- Dual panel: full bands + VB zoom showing dashed (unstrained) vs solid (strained)
- Caption: "Unstrained (dashed) vs strained (solid) GaAs valence bands. The HH/LH splitting at Gamma is clearly visible."

### Figure 3: QW strained band edges (new sim)
- File: `docs/figures/qw_strained_band_edges.png`
- Placement: Sec 7.1 (Example A), after Step 4 summary table
- Dual panel: unstrained (left) vs strained (right) potential profile
- Left panel: flat HH=LH=SO edges per material, CB offset
- Right panel: HH/LH/SO split apart, CB shifted — all with BP formulas applied
- Run config: AlSb/GaSb/InAs broken-gap QW with strain enabled
- Data source: `output/potential_profile.dat` from strained and unstrained runs
- Caption: "Band edge profile for the AlSb/GaSb/InAs QW: unstrained (left) vs strained (right)."

### Figure 4: Wire 2D strain map (new sim)
- File: `docs/figures/wire_strain_2d.png`
- Placement: Sec 6.6 (wire strain result), after the strain result description
- Single panel: colormap of eps_xx(x,y) on wire cross-section
- Shows strain concentration at core corners, decay in matrix
- Run config: InAs square core in GaAs matrix (from existing strain test parameters)
- Data source: `output/strain.dat` from wire run
- Caption: "Biaxial strain component eps_xx in an InAs/GaAs wire cross-section."

### Figure 5: Lattice mismatch cartoon (schematic)
- File: `docs/figures/strain_lattice_mismatch.png`
- Placement: Sec 1.1 (substrate vs layer)
- Two panels: "Free-standing" (different grid spacings) and "Pseudomorphic" (forced to match, with arrows)
- Matplotlib drawing: simple 2D grids with different spacing
- Caption: "Lattice mismatch and pseudomorphic strain. The epitaxial layer is forced to adopt the substrate lattice constant."

### Figure 6: Biaxial strain tensor diagram (schematic)
- File: `docs/figures/strain_biaxial_tensor.png`
- Placement: Sec 2.4 (complete tensor)
- 2D isometric view of a rectangular block with labeled arrows
- In-plane (x,y): blue double arrows labeled eps_xx = eps_yy
- Out-of-plane (z): red arrow labeled eps_zz = -2C12/C11 * eps_xx
- Shear: cross (zero) labeled eps_yz = 0
- Caption: "Biaxial strain geometry for a (001)-oriented quantum well."

### Figure 7: Bir-Pikus band edge shifts (schematic)
- File: `docs/figures/bir_pikus_band_shifts.png`
- Placement: Sec 3.2 (band edge shifts)
- Dual panel: compressive (left) + tensile (right)
- Each panel: vertical energy axis, horizontal lines for CB, HH, LH, SO
- Unstrained positions as dashed, strained as solid
- Arrows showing delta_Ec, delta_EHH, delta_ELH, delta_ESO
- Labels with formulas
- Caption: "Bir-Pikus strain-induced band edge shifts under compressive (left) and tensile (right) biaxial strain."

## Implementation Notes

- All new figure functions go into `scripts/plotting/generate_all_figures.py`
- Schematic figures are pure matplotlib (no external dependencies)
- Simulation figures require running Fortran executables (use existing `run_executable` helper)
- For Figure 3, need to run both strained and unstrained versions of the QW config
- For Figure 4, may need to write a config file for an InAs/GaAs wire with strain
- After generating figures, add `![...]` references in the appropriate sections of `04-strain.md`
- Also update the code reference table (Sec 9) since `apply_pikus_bir` was renamed to `compute_bir_pikus_blocks`

## Files to Modify

- `docs/lecture/04-strain.md` — add figure references
- `scripts/plotting/generate_all_figures.py` — add 5 new figure functions
- `docs/figures/` — 5 new PNG files
