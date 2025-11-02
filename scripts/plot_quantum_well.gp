#!/usr/bin/gnuplot
# Quantum Well Plotting Script for 8bandkp-fdm-ai
# Purpose: Create publication-quality quantum well visualization plots
# Usage: gnuplot -e "datafile='eigenvalues.dat'; output='quantum_well.png'" plot_quantum_well.gp
# Date: 2025-01-27

# Set default values if not provided
if (!exists("datafile")) datafile = "eigenvalues.dat"
if (!exists("output")) output = "quantum_well.png"
if (!exists("title")) title = "Quantum Well Band Structure"

# Set terminal and output
set terminal pngcairo enhanced color size 1000,700 font "Arial,12"
set output output

# Set style
set style line 1 lc rgb '#1f77b4' lw 2 pt 7 ps 0.5
set style line 2 lc rgb '#ff7f0e' lw 2 pt 7 ps 0.5
set style line 3 lc rgb '#2ca02c' lw 2 pt 7 ps 0.5
set style line 4 lc rgb '#d62728' lw 2 pt 7 ps 0.5
set style line 5 lc rgb '#9467bd' lw 2 pt 7 ps 0.5
set style line 6 lc rgb '#8c564b' lw 2 pt 7 ps 0.5

# Set grid
set grid xtics ytics mxtics mytics
set grid ls 0 lc rgb '#808080' lw 0.5

# Set labels and title
set title title font "Arial,16" offset 0,-0.5
set xlabel "k-point" font "Arial,14"
set ylabel "Energy (eV)" font "Arial,14"

# Set key
set key top right font "Arial,10"
set key spacing 1.2

# Set margins
set lmargin 12
set rmargin 5
set tmargin 3
set bmargin 5

# Set x-axis
set xrange [*:*]
set xtics nomirror
set mxtics 5

# Set y-axis
set yrange [*:*]
set ytics nomirror
set mytics 5

# Plot quantum well band structure
# Plot first 10 bands to show subband structure
plot datafile using 1:2 with lines ls 1 title "Subband 1", \
     datafile using 1:3 with lines ls 2 title "Subband 2", \
     datafile using 1:4 with lines ls 3 title "Subband 3", \
     datafile using 1:5 with lines ls 4 title "Subband 4", \
     datafile using 1:6 with lines ls 5 title "Subband 5", \
     datafile using 1:7 with lines ls 6 title "Subband 6", \
     datafile using 1:8 with lines ls 1 dt 2 title "Subband 7", \
     datafile using 1:9 with lines ls 2 dt 2 title "Subband 8", \
     datafile using 1:10 with lines ls 3 dt 2 title "Subband 9", \
     datafile using 1:11 with lines ls 4 dt 2 title "Subband 10"

# Add horizontal line at E=0 (Fermi level)
set arrow from graph 0,0 to graph 1,0 nohead lc rgb '#000000' lw 1 dt 2

# Add vertical lines to indicate k-points
set arrow from 0,graph 0 to 0,graph 1 nohead lc rgb '#808080' lw 1 dt 3

# Print information
print "Quantum well plot saved as: " . output
print "Data file: " . datafile
print "Number of k-points: " . (system("wc -l < " . datafile) - 1)
print "Number of bands: " . (system("head -1 " . datafile . " | wc -w") - 1)



