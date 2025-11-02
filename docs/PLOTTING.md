# Plotting Documentation: 8bandkp-fdm-ai

**Date**: 2025-01-27  
**Status**: Complete

## Overview

This document provides comprehensive guidance for creating publication-quality plots from 8bandkp-fdm-ai calculation results using gnuplot. The plotting system includes automated scripts and customizable templates for various visualization needs.

## Prerequisites

### Required Software
- **gnuplot**: 5.0+ (for plotting)
- **ImageMagick**: Optional (for image processing)
- **LaTeX**: Optional (for enhanced text rendering)

### Installation

**Ubuntu/Debian**:
```bash
sudo apt install gnuplot imagemagick
```

**CentOS/RHEL**:
```bash
sudo yum install gnuplot ImageMagick
```

**Arch/Manjaro**:
```bash
sudo pacman -S gnuplot imagemagick
```

**macOS**:
```bash
brew install gnuplot imagemagick
```

## Available Plotting Scripts

### 1. Band Structure Plotting (`plot_band_structure.gp`)

**Purpose**: Create band structure plots for bulk and quantum well systems  
**Usage**:
```bash
# After running a calculation to outputs/<run-id>/bulk/eigenvalues.dat
gnuplot -e "datafile='outputs/<run-id>/bulk/eigenvalues.dat'; output='band_structure.png'" scripts/plot_band_structure.gp
```

**Features**:
- Multiple band visualization
- Customizable colors and styles
- Grid and axis labels
- Publication-quality output

**Parameters**:
- `datafile`: Input data file (default: eigenvalues.dat)
- `output`: Output image file (default: band_structure.png)
- `title`: Plot title (default: "Band Structure")

### 2. Quantum Well Plotting (`plot_quantum_well.gp`)

**Purpose**: Create quantum well subband structure plots  
**Usage**:
```bash
gnuplot -e "datafile='outputs/<run-id>/qw/eigenvalues.dat'; output='quantum_well.png'" scripts/plot_quantum_well.gp
```

**Features**:
- Subband structure visualization
- Multiple subband plotting
- Quantum confinement effects
- Enhanced styling for quantum wells

**Parameters**:
- `datafile`: Input data file (default: eigenvalues.dat)
- `output`: Output image file (default: quantum_well.png)
- `title`: Plot title (default: "Quantum Well Band Structure")

### 3. G-Factor Plotting (`plot_gfactor.gp`)

**Purpose**: Create g-factor dependence plots  
**Usage**:
```bash
gnuplot -e "datafile='gfactor_data.dat'; output='gfactor.png'" scripts/plot_gfactor.gp
```

**Features**:
- G-factor component visualization
- Anisotropy analysis
- Reference lines (g=0, g=2)
- Customizable styling

**Parameters**:
- `datafile`: Input data file (default: gfactor_data.dat)
- `output`: Output image file (default: gfactor.png)
- `title`: Plot title (default: "G-Factor Dependence")

### 4. Automated Plotting (`plot_all.sh`)

**Purpose**: Generate all plots automatically  
**Usage**:
```bash
./scripts/plot_all.sh
```

**Features**:
- Automatic plot generation
- Multiple plot types
- Publication-quality output
- Combined analysis plots

## Data Format Requirements

### eigenvalues.dat Format
```
#k, values
0.00000    -0.153886    -0.150896    0.843198    0.175542
0.100000E-01  -0.166952  -0.151226    0.172399    0.116615
...
```

**Columns**:
- Column 1: k-point values
- Columns 2+: Energy eigenvalues for each band

### gfactor_data.dat Format
```
# parameter gx gy gz
0.0  2.1  2.1  2.1
0.1  2.05  2.15  2.0
...
```

**Columns**:
- Column 1: Parameter values
- Column 2: G-factor x-component
- Column 3: G-factor y-component
- Column 4: G-factor z-component

## Customization Options

### Color Schemes

**Default Colors**:
- Band 1: Blue (#1f77b4)
- Band 2: Orange (#ff7f0e)
- Band 3: Green (#2ca02c)
- Band 4: Red (#d62728)
- Band 5: Purple (#9467bd)
- Band 6: Brown (#8c564b)

**Custom Colors**:
```gnuplot
set style line 1 lc rgb '#your_color' lw 2 pt 7 ps 0.5
```

### Line Styles

**Available Styles**:
- `lw`: Line width (1-5)
- `pt`: Point type (1-20)
- `ps`: Point size (0.1-2.0)
- `dt`: Dash type (1=solid, 2=dashed, 3=dotted)

### Font Settings

**Font Configuration**:
```gnuplot
set terminal pngcairo enhanced color size 800,600 font "Arial,12"
set title "Your Title" font "Arial,16"
set xlabel "X Label" font "Arial,14"
set ylabel "Y Label" font "Arial,14"
```

## Advanced Plotting

### Multiplot Layouts

**2x2 Grid**:
```gnuplot
set multiplot layout 2,2
# Plot 1
plot datafile using 1:2 with lines
# Plot 2
plot datafile using 1:3 with lines
# Plot 3
plot datafile using 1:4 with lines
# Plot 4
plot datafile using 1:5 with lines
unset multiplot
```

### 3D Plots

**Surface Plot**:
```gnuplot
set terminal pngcairo enhanced color size 800,600
set output "3d_plot.png"
set dgrid3d 50,50
set hidden3d
splot datafile using 1:2:3 with lines
```

### Animation

**GIF Animation**:
```gnuplot
set terminal gif animate delay 10
set output "animation.gif"
do for [i=1:100] {
    plot datafile using 1:2 with lines
}
```

## Publication Guidelines

### High-Resolution Output

**For Publications**:
```gnuplot
set terminal pngcairo enhanced color size 1600,1200 font "Arial,16"
set output "publication_plot.png"
```

**For Presentations**:
```gnuplot
set terminal pngcairo enhanced color size 1200,900 font "Arial,14"
set output "presentation_plot.png"
```

### LaTeX Integration

**EPS Output**:
```gnuplot
set terminal postscript eps enhanced color size 8,6
set output "plot.eps"
```

**PDF Output**:
```gnuplot
set terminal pdf enhanced color size 8,6
set output "plot.pdf"
```

## Troubleshooting

### Common Issues

**1. Gnuplot Not Found**
```
Error: gnuplot not found
Solution: Install gnuplot using package manager
```

**2. Data File Not Found**
```
Error: Cannot open data file
Solution: Check file path and ensure data file exists
```

**3. Empty Plots**
```
Error: No data points plotted
Solution: Check data file format and column numbers
```

**4. Font Issues**
```
Error: Font not found
Solution: Use system fonts or install required fonts
```

### Performance Optimization

**Large Data Sets**:
- Use `set samples` to limit data points
- Use `set dgrid3d` for 3D plots
- Consider data preprocessing

**Memory Issues**:
- Reduce plot resolution
- Use fewer data points
- Close unnecessary terminals

## Examples

### Basic Band Structure Plot
```bash
# Generate band structure plot (replace <run-id>)
gnuplot -e "datafile='outputs/<run-id>/bulk/eigenvalues.dat'; output='band_structure.png'" scripts/plot_band_structure.gp
```

### Custom Quantum Well Plot
```bash
# Generate custom quantum well plot
gnuplot -e "datafile='outputs/<run-id>/qw/eigenvalues.dat'; output='qw_custom.png'; title='Custom QW Plot'" scripts/plot_quantum_well.gp
```

### Batch Processing
```bash
# Generate all plots
./scripts/plot_all.sh

# Generate specific plots
gnuplot scripts/plot_band_structure.gp
gnuplot scripts/plot_quantum_well.gp
gnuplot scripts/plot_gfactor.gp
```

## Next Steps

1. Install gnuplot and required dependencies
2. Run example calculations to generate data files
3. Use plotting scripts to visualize results
4. Customize plots for specific needs
5. Integrate plots into reports and publications
6. Develop additional plotting scripts as needed
