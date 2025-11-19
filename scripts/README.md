# Plotting Scripts for 8-Band k·p FDM Simulator

This directory contains gnuplot scripts for visualizing simulation results.

## Prerequisites

- Gnuplot installed (`sudo dnf install gnuplot` or `sudo apt install gnuplot`)
- Simulation results in an output directory

## Quick Start

### 1. Run a Simulation

From the project root directory:

```bash
# For bulk band structure
./bandStructure examples/bulk.example --out run/bulk_output

# For quantum well g-factor calculation
./gfactorCalculation examples/gfactor.example --out run/gfactor_output
```

### 2. Generate Plots

Navigate to the output directory and run the appropriate plotting script:

```bash
cd run/bulk_output
gnuplot ../../scripts/plot_band_structure.gp

cd ../gfactor_output
gnuplot ../../scripts/plot_gfactor.gp
gnuplot ../../scripts/plot_quantum_well.gp
```

## Available Plotting Scripts

### `plot_band_structure.gp`

**Purpose**: Visualizes the electronic band structure  
**Input**: `eigenvalues.dat` (k-point, 8 eigenvalues)  
**Output**: `band_structure.png`  
**Shows**: Conduction and valence bands vs k-vector

**Example**:
```bash
cd run/bulk_output
gnuplot ../../scripts/plot_band_structure.gp
```

### `plot_gfactor.gp`

**Purpose**: Visualizes g-factor components  
**Input**: `gfactor.dat` (gx, gy, gz values)  
**Output**: `gfactor.png`  
**Shows**: Bar chart of g-factor anisotropy (gx, gy, gz)

**Example**:
```bash
cd run/gfactor_output
gnuplot ../../scripts/plot_gfactor.gp
```

### `plot_quantum_well.gp`

**Purpose**: Visualizes quantum well wavefunctions  
**Input**: `eigenfunctions_k_00001_ev_XXXXX.dat` files  
**Output**: `quantum_well.png`  
**Shows**: Probability density |ψ|² for selected states vs position

**Example**:
```bash
cd run/gfactor_output
gnuplot ../../scripts/plot_quantum_well.gp
```

## Output File Formats

### eigenvalues.dat
```
# k_value  eig1  eig2  eig3  eig4  eig5  eig6  eig7  eig8
0.00000   -0.27 -0.27 -0.00 -0.00 -0.00 -0.00  0.13  0.13
```

### gfactor.dat
```
gx  gy  gz
-0.00 -0.00 -1666.07
```

### eigenfunctions_k_00001_ev_00001.dat
```
# z_position  real(psi)  imag(psi)
-250.0  0.0001  0.0000
-245.0  0.0002  0.0001
...
```

## Customization

All scripts can be customized by editing the `.gp` files:

- **Image size**: `set terminal pngcairo size WIDTH,HEIGHT`
- **Colors**: Modify `set style line` definitions
- **Labels**: Edit `set xlabel`, `set ylabel`, `set title`
- **States plotted**: Change file numbers in plot commands

## Troubleshooting

**Problem**: "gnuplot: command not found"  
**Solution**: Install gnuplot using your package manager

**Problem**: "No such file or directory"  
**Solution**: Make sure you're in the correct output directory

**Problem**: Empty or missing data files  
**Solution**: Check that the simulation completed successfully

## Examples

### Complete Workflow Example

```bash
# 1. Clean build
cd /data/8bandkp-fdm-ai
make clean && make

# 2. Run bulk calculation
./bandStructure examples/bulk.example --out run/bulk_test

# 3. Generate band structure plot
cd run/bulk_test
gnuplot ../../scripts/plot_band_structure.gp
ls -l *.png  # View generated plot

# 4. Run quantum well g-factor calculation
cd ../..
./gfactorCalculation examples/gfactor.example --out run/gfactor_test

# 5. Generate all plots for quantum well
cd run/gfactor_test
gnuplot ../../scripts/plot_band_structure.gp
gnuplot ../../scripts/plot_gfactor.gp
gnuplot ../../scripts/plot_quantum_well.gp
ls -l *.png  # View all plots
```

## Notes

- All plots are saved as PNG files in the current directory
- Scripts print confirmation messages when plots are generated
- The quantum well script plots states 1, 5, 10, 15, 20 by default
- Edit the scripts to plot different states or change visualization settings
