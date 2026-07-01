# Issue 009: Optics + exciton + scattering parsers

**Type:** AFK
**Blocked by:** Issue 002 (bulk tracer bullet)
**User stories:** US 15

## Parent

PRD: `.scratch/input-parser-toml/PRD.md`

## What to build

Add `[optics]`, `[exciton]`, and `[scattering]` section parsing for optical spectra, exciton, and phonon scattering calculations.

### `[optics]` section
- All optical parameters: `linewidth_lorentzian`, `linewidth_gaussian`, `refractive_index`, `E_range` (3-element: min, max, npts), `temperature`, `carrier_density`, `gain_enabled`, `gain_carrier_density`, `ISBT`, `spontaneous`, `spin_resolved`

### `[exciton]` section
- `Ry` (exciton binding energy), `a0` (exciton Bohr radius)

### `[scattering]` section
- `temperature`, `deformation_potential`, `phonon_energy`

Convert all optics/exciton/scattering configs to TOML.

## Acceptance criteria

- [ ] `[optics]` section parses all parameters correctly
- [ ] `[exciton]` and `[scattering]` sections parse correctly
- [ ] All optics/exciton/scattering regression configs converted to TOML
- [ ] Optical spectra simulation outputs match pre-migration golden data
- [ ] Parser test has optics-specific test case

## Blocked by

- Issue 002 (bulk tracer bullet)
