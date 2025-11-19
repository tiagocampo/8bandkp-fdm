# Gnuplot script for quantum well wavefunction visualization
# Usage: cd <output_directory> && gnuplot ../../scripts/plot_quantum_well.gp

set terminal pngcairo size 1000,700 enhanced font 'Arial,12'
set output 'quantum_well.png'

set title "Quantum Well Wavefunctions" font 'Arial,14'
set xlabel "Position z (Å)"
set ylabel "Probability Density |ψ|^2 (arb. units)"
set grid

# Style settings
set style line 1 lc rgb '#e41a1c' lt 1 lw 2
set style line 2 lc rgb '#377eb8' lt 1 lw 2
set style line 3 lc rgb '#4daf4a' lt 1 lw 2
set style line 4 lc rgb '#984ea3' lt 1 lw 2
set style line 5 lc rgb '#ff7f00' lt 1 lw 2

# Eigenfunction files have format: z real(psi) imag(psi)
# We compute probability density as real^2 + imag^2

# Plot first few states (adjust based on available files)
# File naming: eigenfunctions_k_00001_ev_XXXXX.dat

plot 'eigenfunctions_k_00001_ev_00001.dat' using 1:($2**2 + $3**2) with lines ls 1 title "State 1", \
     'eigenfunctions_k_00001_ev_00005.dat' using 1:($2**2 + $3**2) with lines ls 2 title "State 5", \
     'eigenfunctions_k_00001_ev_00010.dat' using 1:($2**2 + $3**2) with lines ls 3 title "State 10", \
     'eigenfunctions_k_00001_ev_00015.dat' using 1:($2**2 + $3**2) with lines ls 4 title "State 15", \
     'eigenfunctions_k_00001_ev_00020.dat' using 1:($2**2 + $3**2) with lines ls 5 title "State 20"

print "Quantum well wavefunction plot saved to quantum_well.png"
