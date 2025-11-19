# Enhanced gnuplot script for quantum well envelope functions (CB and VB states)
# Usage: cd <output_directory> && gnuplot ../../scripts/plot_qw_envelopes.gp

set terminal pngcairo size 1400,1000 enhanced font 'Arial,11'
set output 'envelope_functions.png'

set multiplot layout 2,3 title "Quantum Well Envelope Functions" font 'Arial,16'

# Style settings
set style line 1 lc rgb '#e41a1c' lt 1 lw 2
set style line 2 lc rgb '#377eb8' lt 1 lw 2
set style line 3 lc rgb '#4daf4a' lt 1 lw 2
set style line 4 lc rgb '#984ea3' lt 1 lw 2
set style line 5 lc rgb '#ff7f00' lt 1 lw 2
set style line 6 lc rgb '#a65628' lt 1 lw 2

# Plot 1: First conduction band state (lowest CB)
set title "CB State 1 (Ground)" font 'Arial,12'
set xlabel "Position z (Å)"
set ylabel "|ψ|^2 (arb. units)"
set grid
plot 'eigenfunctions_k_00001_ev_00033.dat' using 1:($2**2 + $3**2) with lines ls 2 notitle

# Plot 2: Second conduction band state
set title "CB State 2" font 'Arial,12'
plot 'eigenfunctions_k_00001_ev_00034.dat' using 1:($2**2 + $3**2) with lines ls 2 notitle

# Plot 3: Third conduction band state
set title "CB State 3" font 'Arial,12'
plot 'eigenfunctions_k_00001_ev_00035.dat' using 1:($2**2 + $3**2) with lines ls 2 notitle

# Plot 4: First valence band state (highest VB / heavy hole)
set title "VB State 1 (Heavy Hole)" font 'Arial,12'
set xlabel "Position z (Å)"
set ylabel "|ψ|^2 (arb. units)"
plot 'eigenfunctions_k_00001_ev_00032.dat' using 1:($2**2 + $3**2) with lines ls 1 notitle

# Plot 5: Second valence band state
set title "VB State 2 (Light Hole)" font 'Arial,12'
plot 'eigenfunctions_k_00001_ev_00031.dat' using 1:($2**2 + $3**2) with lines ls 1 notitle

# Plot 6: Third valence band state
set title "VB State 3" font 'Arial,12'
plot 'eigenfunctions_k_00001_ev_00030.dat' using 1:($2**2 + $3**2) with lines ls 1 notitle

unset multiplot
print "Envelope functions plot saved to envelope_functions.png"
