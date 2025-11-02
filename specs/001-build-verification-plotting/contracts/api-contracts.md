# API Contracts: Build Verification and Plotting

**Feature**: Build Verification and Plotting  
**Date**: 2025-01-27  
**Status**: Complete

## Build System Contract

### Make Targets
```
make all                    # Build both executables
make bandStructure         # Build band structure calculator
make gfactorCalculation    # Build g-factor calculator
make clean                 # Remove build artifacts
make help                  # Show available targets
```

### Build Output Contract
- **bandStructure**: Executable for band structure calculations
- **gfactorCalculation**: Executable for g-factor calculations
- **Object files**: .o files for each source module
- **Module files**: .mod files for Fortran modules

## Executable Interface Contracts

### bandStructure Executable
**Invocation**: `./bandStructure [CONFIG_FILE] [--out outputs/<run-id>]`  
**Input**: Configuration file path (e.g., `examples/bulk_InAs60Sb40.example`). If omitted, defaults to `input.cfg` in CWD.  
**Output**: Plain text files with band structure data  
**Error Handling**: Clear error messages to stderr  
**Exit Codes**: 0=success, 1=input error, 2=calculation error, 3=output error

### gfactorCalculation Executable
**Invocation**: `./gfactorCalculation [CONFIG_FILE] [--out outputs/<run-id>]`  
**Input**: Configuration file path (e.g., `examples/gfactor_quantum_well.example`). If omitted, defaults to `input.cfg` in CWD.  
**Output**: Plain text files with g-factor data  
**Error Handling**: Clear error messages to stderr  
**Exit Codes**: 0=success, 1=input error, 2=calculation error, 3=output error

## Data Format Contracts

### Input Configuration Format
```
waveVector: <direction>           # kx, ky, kz, or k0 (labels case-sensitive; colon optional if separated by whitespace)
waveVectorMax: <float>            # Percentage of BZ
waveVectorStep: <integer>         # Number of k-points
confinement: <0|1>               # 0=bulk, 1=quantum well
FDstep: <integer>                # Discretization points
numLayers: <integer>             # Number of layers
material<N>: <name> <start> <end> [offset]  # Material definitions (offset optional if not used)
numcb: <integer>                 # Conduction bands
numvb: <integer>                 # Valence bands
ExternalField: <0|1> <type>      # Field configuration
EFParams: <float>                # Field parameters
```

### Output Data Format
**Band Structure Output (written to outputs/<run-id>/)**:
```
# k-point energy1 energy2 energy3 ...
0.0000  -0.1234  -0.0567  0.0123
0.0001  -0.1233  -0.0566  0.0124
...
```

**G-Factor Output (written to outputs/<run-id>/)**:
```
# parameter g_factor_x g_factor_y g_factor_z
0.0000  2.1234  2.1234  2.1234
0.0001  2.1235  2.1235  2.1235
...
```

## Plotting Script Contracts

### plot_band_structure.gp
**Input**: Band structure data file  
**Output**: Publication-quality band structure plot  
**Parameters**: Data file path, output file path, plot title  
**Dependencies**: gnuplot 5.0+

### plot_quantum_well.gp
**Input**: Quantum well calculation results  
**Output**: Wave function and energy level plots  
**Parameters**: Data file path, output file path, well parameters  
**Dependencies**: gnuplot 5.0+

### plot_gfactor.gp
**Input**: G-factor calculation results  
**Output**: G-factor dependence plots  
**Parameters**: Data file path, output file path, parameter range  
**Dependencies**: gnuplot 5.0+

## Output Location Contract

- All executables MUST write outputs under a dedicated run directory: `outputs/<run-id>/` where `<run-id>` is a timestamp or unique identifier.
- The directory MUST contain a small metadata file (e.g., `run.meta`) indicating input filename, timestamp, and executable name.

## Validation Contracts

### Build Validation
- All dependencies must be available and linkable
- Compilation must complete without errors
- Executables must be runnable and accept an optional CONFIG_FILE argument

### Calculation Validation
- Input files must parse without errors
- Calculations must converge within reasonable iterations
- Output files must contain valid numerical data
- Results must match expected physical behavior within tolerance

### Plotting Validation
- Scripts must execute without gnuplot errors
- Output plots must be generated successfully
- Plots must be publication-quality with proper labels and formatting
