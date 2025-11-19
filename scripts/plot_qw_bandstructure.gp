# Enhanced gnuplot script for quantum well band structure with more bands
# Usage: cd <output_directory> && gnuplot ../../scripts/plot_qw_bandstructure.gp

set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'bandstructure.png'

set multiplot layout 1,1

set title "Quantum Well Band Structure" font 'Arial,16'
set xlabel "k_x (Å^{-1})" font 'Arial,14'
set ylabel "Energy (eV)" font 'Arial,14'
set grid

# Style settings - conduction bands in blue shades, valence bands in red shades
set style line 1 lc rgb '#8B0000' lt 1 lw 1.5 # Dark red for VB
set style line 2 lc rgb '#DC143C' lt 1 lw 1.5 # Crimson for VB
set style line 3 lc rgb '#FF6347' lt 1 lw 1.5 # Tomato for VB
set style line 4 lc rgb '#00008B' lt 1 lw 2   # Dark blue for CB
set style line 5 lc rgb '#0000FF' lt 1 lw 2   # Blue for CB
set style line 6 lc rgb '#4169E1' lt 1 lw 1.5 # Royal blue for CB

# The eigenvalues.dat file has columns: k_index eigenvalue1 eigenvalue2 ... eigenvalue_N
# We'll plot the first 32 valence bands and first 16 conduction bands

# Find zero energy line
set arrow from graph 0, first 0 to graph 1, first 0 nohead lc rgb 'black' lw 0.5 dt 2

# Plot valence bands (columns 2-33) in red, conduction bands (34-49) in blue
plot for [i=2:17] 'eigenvalues.dat' using 1:i with lines ls 1 notitle, \
     for [i=18:33] 'eigenvalues.dat' using 1:i with lines ls 2 notitle, \
     for [i=34:41] 'eigenvalues.dat' using 1:i with lines ls 4 notitle, \
     for [i=42:49] 'eigenvalues.dat' using 1:i with lines ls 5 notitle, \
     for [i=50:65] 'eigenvalues.dat' using 1:i with lines ls 6 notitle

unset multiplot
print "Band structure plot saved to bandstructure.png"
