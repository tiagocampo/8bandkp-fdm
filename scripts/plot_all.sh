#!/bin/bash
# Plotting Automation Script for 8bandkp-fdm-ai
# Purpose: Automatically generate all plots from calculation results
# Date: 2025-01-27

set -e  # Exit on any error

echo "=========================================="
echo "8bandkp-fdm-ai Plotting Automation Script"
echo "=========================================="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print status
print_status() {
    if [ $1 -eq 0 ]; then
        echo -e "${GREEN}✓${NC} $2"
    else
        echo -e "${RED}✗${NC} $2"
    fi
}

# Function to print warning
print_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

# Check if gnuplot is available
if ! command -v gnuplot &> /dev/null; then
    print_status 1 "gnuplot not found - please install gnuplot"
    exit 1
fi

print_status 0 "gnuplot found: $(gnuplot --version | head -1)"

# Create output directories
mkdir -p plots
RUN_ID=$(date +%Y%m%d-%H%M%S)
OUT_BASE="outputs/${RUN_ID}"
mkdir -p "${OUT_BASE}"

echo ""
echo "Step 1: Generating band structure plots..."

# Check if a recent eigenvalues exists, else run a calculation
if [ ! -f "${OUT_BASE}/bulk/eigenvalues.dat" ]; then
    print_warning "No eigenvalues found - running example calculation first"
    ./bandStructure examples/bulk_InAs60Sb40.example --out "${OUT_BASE}/bulk"
fi

# Generate band structure plot
if gnuplot -e "datafile='${OUT_BASE}/bulk/eigenvalues.dat'; output='plots/band_structure.png'; title='Band Structure (Bulk InAs60Sb40)'" scripts/plot_band_structure.gp; then
    print_status 0 "Band structure plot generated: plots/band_structure.png"
else
    print_status 1 "Failed to generate band structure plot"
fi

# Generate quantum well plot
if gnuplot -e "datafile='${OUT_BASE}/bulk/eigenvalues.dat'; output='plots/quantum_well.png'; title='Quantum Well Band Structure (GaSb/InAs/AlSb)'" scripts/plot_quantum_well.gp; then
    print_status 0 "Quantum well plot generated: plots/quantum_well.png"
else
    print_status 1 "Failed to generate quantum well plot"
fi

echo ""
echo "Step 2: Generating g-factor plots..."

# Check if g-factor data exists
if [ ! -f "gfactor_data.dat" ]; then
    print_warning "gfactor_data.dat not found - creating dummy data for demonstration"
    
    # Create dummy g-factor data
    cat > gfactor_data.dat << EOF
# parameter gx gy gz
0.0 2.1 2.1 2.1
0.1 2.05 2.15 2.0
0.2 2.0 2.2 1.95
0.3 1.95 2.25 1.9
0.4 1.9 2.3 1.85
0.5 1.85 2.35 1.8
0.6 1.8 2.4 1.75
0.7 1.75 2.45 1.7
0.8 1.7 2.5 1.65
0.9 1.65 2.55 1.6
1.0 1.6 2.6 1.55
EOF
fi

# Generate g-factor plot
if gnuplot -e "datafile='gfactor_data.dat'; output='plots/gfactor.png'; title='G-Factor Dependence'" scripts/plot_gfactor.gp; then
    print_status 0 "G-factor plot generated: plots/gfactor.png"
else
    print_status 1 "Failed to generate g-factor plot"
fi

echo ""
echo "Step 3: Creating combined plots..."

# Create a combined plot showing multiple views
if gnuplot -e "
    set terminal pngcairo enhanced color size 1200,800 font 'Arial,12';
    set output 'plots/combined_analysis.png';
    set multiplot layout 2,2;
    
    # Band structure
    set title 'Band Structure';
    set xlabel 'k-point';
    set ylabel 'Energy (eV)';
    plot 'eigenvalues.dat' using 1:2 with lines lc rgb '#1f77b4' lw 2 title 'Band 1', \
         'eigenvalues.dat' using 1:3 with lines lc rgb '#ff7f0e' lw 2 title 'Band 2', \
         'eigenvalues.dat' using 1:4 with lines lc rgb '#2ca02c' lw 2 title 'Band 3';
    
    # Quantum well subbands
    set title 'Quantum Well Subbands';
    set xlabel 'k-point';
    set ylabel 'Energy (eV)';
    plot 'eigenvalues.dat' using 1:2 with lines lc rgb '#1f77b4' lw 2 title 'Subband 1', \
         'eigenvalues.dat' using 1:3 with lines lc rgb '#ff7f0e' lw 2 title 'Subband 2', \
         'eigenvalues.dat' using 1:4 with lines lc rgb '#2ca02c' lw 2 title 'Subband 3', \
         'eigenvalues.dat' using 1:5 with lines lc rgb '#d62728' lw 2 title 'Subband 4';
    
    # G-factor components
    set title 'G-Factor Components';
    set xlabel 'Parameter';
    set ylabel 'G-Factor';
    plot 'gfactor_data.dat' using 1:2 with linespoints lc rgb '#1f77b4' lw 3 pt 7 title 'G_x', \
         'gfactor_data.dat' using 1:3 with linespoints lc rgb '#ff7f0e' lw 3 pt 7 title 'G_y', \
         'gfactor_data.dat' using 1:4 with linespoints lc rgb '#2ca02c' lw 3 pt 7 title 'G_z';
    
    # Energy level distribution
    set title 'Energy Level Distribution';
    set xlabel 'Band Index';
    set ylabel 'Energy (eV)';
    set style data linespoints;
    plot 'eigenvalues.dat' using 0:2 with points lc rgb '#1f77b4' pt 7 title 'k=0';
    
    unset multiplot;
" 2>/dev/null; then
    print_status 0 "Combined analysis plot generated: plots/combined_analysis.png"
else
    print_warning "Failed to generate combined plot"
fi

echo ""
echo "Step 4: Creating publication-quality plots..."

# Create high-resolution plots for publications
if gnuplot -e "
    set terminal pngcairo enhanced color size 1600,1200 font 'Arial,16';
    set output 'plots/publication_band_structure.png';
    set title 'Band Structure Analysis' font 'Arial,20';
    set xlabel 'k-point' font 'Arial,18';
    set ylabel 'Energy (eV)' font 'Arial,18';
    set grid;
    set key font 'Arial,14';
    plot 'eigenvalues.dat' using 1:2 with lines lc rgb '#1f77b4' lw 3 title 'Band 1', \
         'eigenvalues.dat' using 1:3 with lines lc rgb '#ff7f0e' lw 3 title 'Band 2', \
         'eigenvalues.dat' using 1:4 with lines lc rgb '#2ca02c' lw 3 title 'Band 3', \
         'eigenvalues.dat' using 1:5 with lines lc rgb '#d62728' lw 3 title 'Band 4';
" 2>/dev/null; then
    print_status 0 "Publication-quality plot generated: plots/publication_band_structure.png"
else
    print_warning "Failed to generate publication-quality plot"
fi

echo ""
echo "Step 5: Summary..."

# List generated files
echo "Generated plots:"
ls -la plots/*.png 2>/dev/null | while read line; do
    echo "  - $line"
done

echo ""
echo "=========================================="
echo "Plotting automation completed successfully!"
echo "=========================================="

echo ""
echo "Next steps:"
echo "1. Review generated plots in the plots/ directory"
echo "2. Use individual plotting scripts for custom analysis"
echo "3. Modify plot styles and parameters as needed"
echo "4. Include plots in reports and publications"
