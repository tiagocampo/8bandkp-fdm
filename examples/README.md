# Example Files

This directory contains working example input files for the 8-band k·p solver.

## Bandstructure Calculations

### Bulk Materials
- `bulk_GaAs.example` - GaAs bulk bandstructure
- `validation_bulk_*.example` - Additional bulk material validations

### Quantum Wells
- `qw_gaas_algaas.example` - GaAs/Al₀.₃Ga₀.₇As quantum well
- `qw_ingaas_gaas_strained.example` - Strained In₀.₂Ga₀.₈As/GaAs QW
- `qw_ingaas_inp.example` - InGaAs/InP quantum well
- `qw_inas_alsb.example` - InAs/AlSb quantum well
- `qw_inas_gasb_alsb_clean.example` - InAs/GaSb/AlSb type-II QW

## G-Factor Calculations

### Analytical Method (Bulk & QW)
- `gfactor_analytical_bulk_GaAs.example` - Bulk GaAs analytical g-factor
- `gfactor_bulk_*_analytical.example` - Other bulk materials
- `gfactor_qw_*_analytical.example` - Quantum well analytical g-factors

**Note**: Analytical method works for both bulk and quantum wells.

### Numerical Method (Bulk Only)
- `gfactor_bulk_GaAs_numerical.example` - Bulk GaAs numerical g-factor
- `gfactor_bulk_InAs_numerical.example` - Bulk InAs numerical g-factor
- `gfactor_bulk_*_numerical.example` - Other bulk numerical examples

**Note**: Numerical method currently works reliably for bulk materials only. For quantum wells, use the analytical method.

## Archived Examples

The `archive/` subdirectory contains experimental or non-working examples:
- Quantum well numerical g-factor examples (under development)
- Legacy example files
- Test configurations

## Usage

**Bandstructure**:
```bash
./bandStructure examples/bulk_GaAs.example
```

**G-Factor**:
```bash
./gfactorCalculation < examples/gfactor_analytical_bulk_GaAs.example
```

See `docs/` for detailed documentation.
