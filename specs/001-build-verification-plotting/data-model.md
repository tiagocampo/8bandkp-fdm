# Data Model: Build Verification and Plotting

**Feature**: Build Verification and Plotting  
**Date**: 2025-01-27  
**Status**: Complete

## Entity Definitions

### Build Configuration
**Purpose**: Manages compilation settings and dependency management  
**Fields**:
- `compiler`: Fortran compiler executable path
- `compiler_flags`: Compilation flags (optimization, debugging)
- `blas_lib`: BLAS library path and linking flags
- `lapack_lib`: LAPACK library path and linking flags
- `fftw_lib`: FFTW3 library path and linking flags
- `makefile_targets`: Available make targets (all, clean, bandStructure, gfactorCalculation)

**Validation Rules**:
- Compiler must be executable and support Fortran 90/95
- Libraries must be linkable and compatible
- Makefile must define all required targets

### Input Parameters
**Purpose**: Configuration for calculation runs  
**Fields**:
- `wave_vector`: Reciprocal space direction (kx, ky, kz, k0)
- `wave_vector_max`: Percentage of Brillouin Zone to compute
- `wave_vector_step`: Number of k-points
- `confinement`: System type (0=bulk, 1=quantum well)
- `fd_step`: Number of discretization points
- `num_layers`: Number of quantum well layers
- `materials`: Array of material definitions with positions and band offsets
- `num_cb`: Number of conduction bands
- `num_vb`: Number of valence bands
- `external_field`: Field type and parameters

**Validation Rules**:
- Wave vector must be valid direction
- Material positions must be physically reasonable
- Band counts must not exceed system limits
- Field parameters must be within numerical stability range

### Calculation Results
**Purpose**: Output data from simulations  
**Fields**:
- `eigenvalues`: Energy eigenvalues array
- `wave_functions`: Wave function coefficients
- `band_structure`: Energy vs k-point data
- `g_factors`: G-factor values and dependencies
- `metadata`: Calculation parameters, timestamps, convergence info

**Validation Rules**:
- Eigenvalues must be real and finite
- Wave functions must be normalized
- Band structure must show expected physical behavior
- G-factors must be within physically reasonable ranges

### Visualization Scripts
**Purpose**: Gnuplot scripts for result plotting  
**Fields**:
- `script_type`: Type of plot (band_structure, quantum_well, g_factor)
- `input_format`: Expected data file format
- `output_format`: Generated plot format (PNG, PDF, SVG)
- `plot_settings`: Axis labels, ranges, styling
- `dependencies`: Required gnuplot version and features

**Validation Rules**:
- Scripts must be valid gnuplot syntax
- Input format must match calculation output
- Output must be publication-quality
- Dependencies must be clearly documented

## State Transitions

### Build Process States
1. **Clean** → **Configured**: Dependencies detected and Makefile configured
2. **Configured** → **Compiling**: Make process initiated
3. **Compiling** → **Built**: Executables created successfully
4. **Compiling** → **Failed**: Build errors encountered

### Calculation States
1. **Input Validated** → **Running**: Calculation initiated
2. **Running** → **Completed**: Results generated successfully
3. **Running** → **Failed**: Numerical or convergence errors
4. **Completed** → **Validated**: Results verified against benchmarks

### Plotting States
1. **Data Available** → **Script Executing**: Gnuplot script running
2. **Script Executing** → **Plot Generated**: Visualization created
3. **Script Executing** → **Failed**: Script errors or missing data

## Relationships

- **Build Configuration** → **Input Parameters**: Build must support input parameter ranges
- **Input Parameters** → **Calculation Results**: Parameters determine result structure
- **Calculation Results** → **Visualization Scripts**: Results provide data for plotting
- **Build Configuration** → **Visualization Scripts**: Build must ensure compatible output formats
