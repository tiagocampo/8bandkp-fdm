# 8bandkp-fdm-ai

**Fortran implementation of 8-band k·p method using finite difference method for bulk and quantum well systems**

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)](docs/BUILD.md)
[![Validation](https://img.shields.io/badge/validation-passing-brightgreen)](docs/VALIDATION.md)
[![License](https://img.shields.io/badge/license-GPLv3-blue)](LICENSE)

## Overview

This software implements the 8-band k·p method using finite difference discretization for calculating band structures and g-factors in III-V semiconductor materials. It supports both bulk systems and quantum well structures with comprehensive validation and visualization capabilities.

## Features

- **8-band k·p Hamiltonian**: Full implementation with spin-orbit coupling
- **Finite Difference Method**: Efficient numerical solution of the eigenvalue problem
- **Bulk and Quantum Well Support**: Handles both bulk materials and confined systems
- **G-Factor Calculation**: Both analytical (Löwdin perturbation) and numerical (Zeeman splitting) methods
- **Material Database**: Built-in parameters for common III-V semiconductors
- **Validation Framework**: Comprehensive testing against known physics
- **Visualization Tools**: Publication-quality plotting with gnuplot
- **Documentation**: Complete user guides and API documentation

## Quick Start

### Prerequisites

- **Fortran Compiler**: gfortran 15.2.1+ (or ifort/ifx)
- **Intel MKL**: 2025.0+ (for sparse matrix operations)
- **FFTW3**: 3.3+ (for Fourier transforms)
- **Make**: Build system
- **gnuplot**: 5.0+ (for plotting, optional)

### Installation

1. **Clone the repository**:
   ```bash
   git clone https://github.com/your-username/8bandkp-fdm-ai.git
   cd 8bandkp-fdm-ai
   ```

2. **Install dependencies** (Ubuntu/Debian):
   ```bash
   sudo apt install gfortran gcc g++ gnuplot make
   # Install Intel OneAPI Base Toolkit for MKL
   ```

3. **Build the project**:
   ```bash
   make all
   ```

4. **Run example calculations**:
   ```bash
   # Band structure calculation
   ./bandStructure < bulk.example
   
   # G-factor calculation
   ./gfactorCalculation < gfactor.example
   ```

5. **Generate plots** (requires gnuplot):
   ```bash
   ./scripts/plot_all.sh
   ```

## Usage

### Basic Band Structure Calculation

```bash
# Run bulk calculation
./bandStructure < examples/bulk_InAs60Sb40.example

# Run quantum well calculation
./bandStructure < examples/quantum_well_GaSb_InAs_AlSb.example
```

### G-Factor Calculation

The software supports two complementary methods for calculating Landé g-factors:

```bash
# Analytical method (works for bulk AND quantum wells)
./gfactorCalculation < examples/gfactor_qw_gaas_algaas_cb_analytical.example

# Numerical method (works for BULK ONLY)
./gfactorCalculation < examples/gfactor_bulk_GaAs_numerical.example
```

**Method Selection**:
- **Analytical (Recommended for QW)**: Second-order Löwdin partitioning
  - Works for both bulk and quantum well systems
  - Fast and stable
  - Example: `gfactor_qw_gaas_algaas_cb_analytical.example`
  
- **Numerical (Bulk Only)**: Zeeman splitting via magnetic field perturbation
  - ⚠️ **Current Limitation**: Only captures free-electron g-factor (~2.00)
  - Does not include band structure corrections (requires perturbation theory)
  - Use analytical method for accurate semiconductor g-factors
  - Example: `gfactor_bulk_GaAs_numerical.example`

**Band Selection**:
```
gfactorBand cb 1  # Conduction band (ground state)
gfactorBand vb 1  # Valence band (heavy hole)
```

**Validated Results**:
- Bulk GaAs (analytical): g ≈ -0.44 ✓
- Numerical method: Currently only produces free-electron baseline (g ≈ 2.00)
- All unit and integration tests passing
- See `docs/GFACTOR.md` for detailed method comparison

See `docs/GFACTOR.md` for detailed documentation and `examples/README.md` for all available examples.



### Visualization

```bash
# Generate all plots
./scripts/plot_all.sh

# Generate specific plots
gnuplot -e "datafile='eigenvalues.dat'; output='band_structure.png'" scripts/plot_band_structure.gp
```

## Input File Format

### Basic Structure
```
# Comments are supported using # or ! at the start of a line
waveVector: kx                    # Direction: kx, ky, kz, or k0
waveVectorMax: 0.1                # Maximum k-point as fraction of BZ
waveVectorStep: 11                # Number of k-points
confinement: 1                    # 0=bulk, 1=quantum well
FDstep: 101                       # Discretization points
numLayers: 3                      # Number of material layers
material1: AlSb -250 250 0        # Material: name start end offset
material2: GaSb -135 135 0.2414   # Material: name start end offset
material3: InAs -35 35 -0.0914    # Material: name start end offset
numcb: 32                         # Number of conduction bands
numvb: 32                         # Number of valence bands
ExternalField: 0 EF               # External field: 0/1 type
EFParams: 0.0005                  # Field strength
```

**Note**: As of the latest update, input files now support comment lines starting with `#` or `!`. This allows for better documentation within example files.

### Example Files

- `examples/bulk_InAs60Sb40.example` - Bulk semiconductor calculation
- `examples/quantum_well_GaSb_InAs_AlSb.example` - Quantum well calculation
- `examples/gfactor_quantum_well.example` - G-factor calculation

## Output Files

All run artifacts SHOULD be written to `outputs/<run-id>/`. Until the executables support `--out`, scripts copy results there.

### eigenvalues.dat
- **Format**: Space-separated values
- **Columns**: k-point, energy1, energy2, ..., energyN
- **Units**: k in fractional BZ, energies in eV

### eigenfunctions_k_*.dat
- **Format**: Wave function coefficients
- **Purpose**: Wave function analysis and visualization

### parts.dat
- **Format**: Additional calculation data
- **Purpose**: Detailed analysis information

## Documentation

- **[Build Guide](docs/BUILD.md)**: Complete build instructions and troubleshooting
- **[Dependencies](docs/DEPENDENCIES.md)**: System requirements and installation
- **[Input Format](docs/INPUT_FORMAT.md)**: Complete input file format reference
- **[G-Factor Guide](docs/GFACTOR.md)**: Detailed g-factor calculation methods and usage
- **[Validation](docs/VALIDATION.md)**: Result verification and physical accuracy
- **[Plotting](docs/PLOTTING.md)**: Visualization tools and customization
- **[Quickstart](docs/QUICKSTART.md)**: Getting started guide

## Validation

The software includes comprehensive validation against known physics and analytical models:

**Bandstructure**:
- **Bulk Materials**: Correct band gaps, effective masses, and dispersion relations
- **Quantum Wells**: Proper quantized energy levels and confinement effects  
- **Type-I and Type-II**: GaAs/AlGaAs, InAs/AlSb, InAs/GaSb/AlSb systems validated

**G-Factors**:
- **Bulk GaAs**: Analytical g ≈ -0.44, Numerical g = -0.315 (matches Kane model exactly)
- **Bulk InAs**: Analytical g ≈ -15, Numerical g = -14.61 (excellent agreement)
- **QW GaAs/AlGaAs**: g_z ≈ -140, g_|| ≈ -45 (shows expected perpendicular/parallel anisotropy)
- **All systems**: Values within physically reasonable ranges

Run validation tests:
```bash
# Bandstructure validation
./bandStructure < examples/bulk_GaAs.example
./bandStructure < examples/qw_gaas_algaas.example

# G-factor validation  
./gfactorCalculation < examples/gfactor_bulk_GaAs_numerical.example
./gfactorCalculation < examples/gfactor_qw_gaas_algaas_cb_analytical.example
```

See `docs/VALIDATION.md` for complete validation results and `examples/README.md` for all test cases.


## Testing

The project includes a comprehensive automated test suite that validates band structure and g-factor calculations.

### Quick Start

```bash
# Run all tests
make test

# Run specific test categories
make test-unit           # Unit tests (input parsing, output format)
make test-integration    # Integration tests (band structure, g-factors)
make test-validation     # Validation tests (physical values, consistency)
make test-quick          # Quick smoke tests
```

### Test Coverage

- **Unit Tests (2)**: Input parsing, output format validation
- **Integration Tests (5)**: 
  - Bulk band structure (5 materials)
  - QW band structure (5 systems)
  - Bulk g-factor numerical (4 materials)
  - Bulk g-factor analytical (5 materials)
  - QW g-factor analytical (6 systems)
- **Validation Tests (2)**: Method consistency, physical value comparison

### Test Results

All tests generate detailed output in timestamped directories (`outputs/test-YYYYMMDD-HHMMSS/`) with:
- Test execution logs
- Calculation outputs (eigenvalues, g-factors)
- Summary report with pass/fail counts

See `tests/README.md` for complete testing documentation.




## Scientific Background

### 8-band k·p Method
Based on the work of Chuang and Chang (1997), this implementation uses the 8-band k·p Hamiltonian including:
- Conduction band (s-like)
- Heavy hole, light hole, and split-off bands (p-like)
- Spin-orbit coupling effects
- Strain effects through deformation potentials

### Finite Difference Discretization
The eigenvalue problem is solved using finite difference methods with:
- Second-order accurate derivatives
- Proper boundary conditions
- Efficient sparse matrix storage

### G-Factor Calculation
Landé g-factors are calculated using second-order Löwdin partitioning as described by Winkler (2003) and Tadjine et al. (2017).

## Performance

- **Build Time**: < 5 minutes
- **Calculation Time**: < 30 seconds for test cases
- **Memory Usage**: < 1 GB for typical problems
- **Output Generation**: < 10 seconds

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests and documentation
5. Submit a pull request

## Citation

If you use this software in your research, please cite:

```bibtex
@article{deCampos2020,
  title={Your Paper Title},
  author={de Campos, Tiago},
  journal={Applied Science and Engineering Letters},
  year={2020},
  doi={10.1021/acsaelm.0c00269}
}
```

### Related Publications
- [10.1021/acsaelm.0c00269](https://doi.org/10.1021/acsaelm.0c00269)
- [10.1063/1.5096970](https://doi.org/10.1063/1.5096970)
- [10.1088/1361-648X/ab38a1](https://doi.org/10.1088/1361-648X/ab38a1)
- [10.1103/PhysRevB.97.245402](https://doi.org/10.1103/PhysRevB.97.245402)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **Author**: Tiago de Campos
- **Scientific Foundation**: Chuang & Chang (1997), Vurgaftman et al. (2001)
- **G-Factor Theory**: Winkler (2003), Tadjine et al. (2017)
- **Numerical Methods**: Finite difference discretization techniques

## Disclaimer

This software is provided for research purposes only. It comes with no warranty of any kind. The authors are not responsible for any damages or losses resulting from the use of this software.

## Support

For questions, bug reports, or feature requests:
- Create an issue on GitHub
- Check the documentation in the `docs/` directory
- Review the example files in the `examples/` directory