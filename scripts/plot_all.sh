#!/bin/bash

# Create output directory if it doesn't exist
mkdir -p ../output

# Run all gnuplot scripts
echo "Plotting band structure..."
gnuplot plot_bands.gp

echo "Plotting eigenfunctions..."
gnuplot plot_eigenfunctions.gp

# Skip potential profile for bulk calculations
if [ -f "../output/potential_profile.dat" ]; then
    echo "Plotting potential profile..."
    gnuplot plot_potential.gp
else
    echo "Skipping potential profile plot (bulk calculation)"
fi

# Skip g-factor if not calculated
if [ -f "../output/gfactor.dat" ]; then
    echo "Plotting g-factor..."
    gnuplot plot_gfactor.gp
else
    echo "Skipping g-factor plot (not calculated yet)"
fi

# Plot band composition if available
if [ -f "../output/parts.dat" ]; then
    echo "Plotting band composition..."
    gnuplot plot_parts.gp
else
    echo "Skipping band composition plot (parts.dat not found)"
fi

echo "All plots have been generated in the output directory." 