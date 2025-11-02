#!/usr/bin/gnuplot
# G-Factor Plotting Script for 8bandkp-fdm-ai
# Purpose: Create publication-quality g-factor plots
# Usage: gnuplot -e "datafile='gfactor_data.dat'; output='gfactor.png'" plot_gfactor.gp
# Date: 2025-01-27

# Set default values if not provided
if (!exists("datafile")) datafile = "gfactor_data.dat"
if (!exists("output")) output = "gfactor.png"
if (!exists("title")) title = "G-Factor Dependence"

# Set terminal and output
set terminal pngcairo enhanced color size 800,600 font "Arial,12"
set output output

# Set style
set style line 1 lc rgb '#1f77b4' lw 3 pt 7 ps 1.0
set style line 2 lc rgb '#ff7f0e' lw 3 pt 7 ps 1.0
set style line 3 lc rgb '#2ca02c' lw 3 pt 7 ps 1.0
set style line 4 lc rgb '#d62728' lw 3 pt 7 ps 1.0

# Set grid
set grid xtics ytics mxtics mytics
set grid ls 0 lc rgb '#808080' lw 0.5

# Set labels and title
set title title font "Arial,16" offset 0,-0.5
set xlabel "Parameter" font "Arial,14"
set ylabel "G-Factor" font "Arial,14"

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

# Plot g-factors
# Assuming data format: parameter gx gy gz
plot datafile using 1:2 with linespoints ls 1 title "G_x", \
     datafile using 1:3 with linespoints ls 2 title "G_y", \
     datafile using 1:4 with linespoints ls 3 title "G_z"

# Add horizontal line at g=2 (free electron value)
set arrow from graph 0,2 to graph 1,2 nohead lc rgb '#000000' lw 1 dt 2

# Add horizontal line at g=0
set arrow from graph 0,0 to graph 1,0 nohead lc rgb '#808080' lw 1 dt 3

# Print information
print "G-factor plot saved as: " . output
print "Data file: " . datafile
print "Number of data points: " . (system("wc -l < " . datafile) - 1)



