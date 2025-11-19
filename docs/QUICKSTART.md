# Quickstart Guide: 8bandkp-fdm-ai

**Get up and running in 30 minutes!**

## Prerequisites

### Required Software
- **Fortran Compiler**: gfortran 15.2.1+ (or ifort/ifx)
- **Intel MKL**: 2025.0+ (for sparse matrix operations)
- **FFTW3**: 3.3+ (for Fourier transforms)
- **Make**: Build system
- **gnuplot**: 5.0+ (for plotting, optional)

### Installation Examples

**Ubuntu/Debian**:
```bash
sudo apt install gfortran gcc g++ gnuplot make
# Install Intel OneAPI Base Toolkit for MKL
```

**CentOS/RHEL**:
```bash
sudo yum install gcc-gfortran gnuplot make
# Install Intel OneAPI Base Toolkit for MKL
```

**Arch/Manjaro**:
```bash
sudo pacman -S gcc-fortran gnuplot make
# Install Intel OneAPI Base Toolkit for MKL
```

## Step 1: Build the Project (5 minutes)

```bash
# Clone and navigate to project
git clone https://github.com/your-username/8bandkp-fdm-ai.git
cd 8bandkp-fdm-ai

# Build the project
make all

# Verify executables were created
ls -la bandStructure gfactorCalculation
```

**Expected Output**:
```
-rwxr-xr-x 1 user user 343k Jan 27 23:59 bandStructure
-rwxr-xr-x 1 user user 497k Jan 27 23:59 gfactorCalculation
```

## Step 2: Run Your First Calculation (5 minutes)

### Bulk Band Structure
```bash
# Run bulk calculation
./bandStructure examples/bulk.example --out run/my_bulk

# Check results
head -5 run/my_bulk/eigenvalues.dat
```

**Expected Output**:
```
#k, values
0.00000    -0.153886    -0.150896    0.843198    0.175542
0.100000E-01  -0.166952  -0.151226    0.172399    0.116615
0.200000E-01  -0.188651  -0.153064    0.464922    0.176623
0.300000E-01  -0.193247  -0.162278    0.921608    0.244021
0.400000E-01  -0.196300  -0.181261    0.153665    0.315436
```

### G-Factor Calculation
```bash
# Analytical method (quantum well)
./gfactorCalculation examples/gfactor.example --out run/my_gfactor

# Numerical method (bulk GaAs conduction band)
./gfactorCalculation examples/gfactor_bulk_GaAs_numerical.example --out run/gfactor_cb

# Check results
cat run/gfactor_cb/gfactor.dat
```

**Expected Output**:
```
-0.31731403134144992  -0.31731403134144992  -0.31731403134144992
```

For more details on g-factor calculations, see [`docs/GFACTOR.md`](GFACTOR.md).


## Step 3: Visualize Results (10 minutes)

### Install gnuplot (if not already installed)
```bash
# Ubuntu/Debian
sudo apt install gnuplot

# Arch/Manjaro
sudo pacman -S gnuplot

# macOS
brew install gnuplot
```

### Generate Plots
```bash
# Navigate to output directory
cd run/my_gfactor

# Generate individual plots
gnuplot ../../scripts/plot_gfactor.gp
gnuplot ../../scripts/plot_quantum_well.gp

# For band structure
cd ../my_bulk
gnuplot ../../scripts/plot_band_structure.gp
```

**Expected Output**:
```
G-factor plot saved to gfactor.png
Quantum well wavefunction plot saved to quantum_well.png
Band structure plot saved to band_structure.png
```

## Step 4: Validate Results (5 minutes)

```bash
# Verify the calculations produced expected output
ls run/my_bulk/eigenvalues.dat
ls run/my_gfactor/gfactor.dat
ls run/my_gfactor/eigenfunctions_*.dat

# Check g-factor values are reasonable
cat run/my_gfactor/gfactor.dat
```

**Expected Output**:
```
run/my_bulk/eigenvalues.dat
run/my_gfactor/gfactor.dat
run/my_gfactor/eigenfunctions_k_00001_ev_00001.dat
...
-1.57e-13  -2.03e-13  -1666.07  # gx, gy, gz values
```

## Step 5: Explore Examples (5 minutes)

### Available Examples
```bash
ls examples/
```

**Output**:
```
bulk_InAs60Sb40.example
gfactor_quantum_well.example
quantum_well_GaSb_InAs_AlSb.example
validation_bulk_InAs60Sb40.example
validation_quantum_well_GaSb_InAs_AlSb.example
```

### Run Different Calculations
```bash
# Try different examples
./bandStructure examples/bulk.example --out run/test_bulk
./gfactorCalculation examples/gfactor.example --out run/test_gfactor
```

## Troubleshooting

### Common Issues

**1. Build Fails**
```bash
# Check dependencies
gfortran --version
pkg-config --exists fftw3 && echo "FFTW3 found"

# Clean and rebuild
make clean_all
make all
```

**2. Calculation Fails**
```bash
# Check input file format
cat examples/gfactor.example

# Run with output to check for errors
./gfactorCalculation examples/gfactor.example --out run/debug_output
```

**3. Plotting Fails**
```bash
# Check gnuplot installation
gnuplot --version

# Navigate to output directory and test
cd run/my_gfactor
gnuplot ../../scripts/plot_gfactor.gp
```

**4. Permission Denied**
```bash
# Make scripts executable
chmod +x scripts/*.sh
chmod +x scripts/*.gp
```

## Next Steps

### 1. Read the Documentation
- [Build Guide](BUILD.md) - Detailed build instructions
- [Dependencies](DEPENDENCIES.md) - System requirements
- [Validation](VALIDATION.md) - Result verification
- [Plotting](PLOTTING.md) - Visualization tools

### 2. Create Your Own Calculations
```bash
# Copy example file
cp examples/gfactor.example my_calculation.example

# Edit parameters
nano my_calculation.example

# Run calculation
./gfactorCalculation my_calculation.example --out run/my_custom
```

### 3. Customize Plots
```bash
# Edit plotting script
nano scripts/plot_band_structure.gp

# Generate custom plot
gnuplot -e "datafile='eigenvalues.dat'; output='my_plot.png'" scripts/plot_band_structure.gp
```

### 4. Explore Advanced Features
- Electric field effects
- Different material systems
- Custom material parameters
- Batch processing

## Performance Tips

### For Large Calculations
- Increase `FDstep` for higher accuracy
- Use more `numcb` and `numvb` for detailed band structure
- Consider memory usage for very large systems

### For Fast Testing
- Use small `FDstep` values
- Limit `waveVectorStep` for quick tests
- Use fewer bands for initial exploration

## Getting Help

1. **Check the documentation** in the `docs/` directory
2. **Review example files** in the `examples/` directory
3. **Run validation scripts** to verify your setup
4. **Create an issue** on GitHub for bugs or questions

## Success Criteria

You've successfully completed the quickstart if you can:
- ✅ Build the project without errors
- ✅ Run at least one calculation successfully
- ✅ Generate at least one plot
- ✅ Validate results using the provided scripts
- ✅ Understand the basic input file format

**Congratulations! You're ready to use 8bandkp-fdm-ai for your research!**
