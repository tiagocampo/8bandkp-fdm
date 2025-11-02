# Quickstart Guide: Build Verification and Plotting

**Feature**: Build Verification and Plotting  
**Date**: 2025-01-27  
**Status**: Complete

## Prerequisites

### Required Software
- Fortran compiler (gfortran, ifort, or ifx)
- BLAS/LAPACK library (OpenBLAS, Intel MKL, or system BLAS)
- FFTW3 library
- Make build system
- gnuplot 5.0+ (for plotting)

### Installation Examples

**Ubuntu/Debian**:
```bash
sudo apt install gfortran gcc g++ liblapack-dev libblas-dev libfftw3-dev gnuplot make
```

**CentOS/RHEL**:
```bash
sudo yum install gcc-gfortran lapack-devel blas-devel fftw3-devel gnuplot make
```

**Intel OneAPI** (optional, for MKL):
```bash
# Install Intel OneAPI Base Toolkit
# Source environment: source /opt/intel/oneapi/setvars.sh
```

## Build Process

### Step 1: Configure Build
```bash
# Clone and navigate to project
cd /path/to/8bandkp-fdm-ai

# Check dependencies
make check-deps  # If available, or manually verify
```

### Step 2: Compile
```bash
# Build both executables
make all

# Verify executables created
ls -la bandStructure gfactorCalculation
```

### Step 3: Test Build
```bash
# Test band structure calculation
./bandStructure < bulk.example

# Test g-factor calculation  
./gfactorCalculation < gfactor.example
```

## Running Calculations

### Bulk Calculation
```bash
# Run bulk band structure calculation
./bandStructure < bulk.example

# Expected output: eigenvalues.dat, band structure data
```

### Quantum Well Calculation
```bash
# Run quantum well calculation
./bandStructure < quantumwell.example

# Expected output: eigenvalues.dat, wave function data
```

### G-Factor Calculation
```bash
# Run g-factor calculation
./gfactorCalculation < gfactor.example

# Expected output: g-factor data files
```

## Plotting Results

### Band Structure Plot
```bash
# Generate band structure plot
gnuplot -e "datafile='eigenvalues.dat'; output='band_structure.png'" scripts/plot_band_structure.gp
```

### Quantum Well Plot
```bash
# Generate quantum well visualization
gnuplot -e "datafile='quantum_well_data.dat'; output='quantum_well.png'" scripts/plot_quantum_well.gp
```

### G-Factor Plot
```bash
# Generate g-factor plot
gnuplot -e "datafile='gfactor_data.dat'; output='gfactor.png'" scripts/plot_gfactor.gp
```

## Validation

### Expected Results
- **Bulk InAs60Sb40**: Band structure shows expected energy levels
- **GaSb/InAs/AlSb Quantum Well**: Clear confinement effects with quantized levels
- **G-Factors**: Values within physically reasonable ranges (typically 1.5-3.0)

### Troubleshooting
- **Build fails**: Check compiler and library versions
- **Calculation fails**: Verify input file format and parameters
- **Plotting fails**: Ensure gnuplot is installed and data files exist

## Next Steps
- Modify input parameters for different materials
- Create custom plotting scripts for specific analyses
- Validate results against published literature
- Contribute improvements back to the project
