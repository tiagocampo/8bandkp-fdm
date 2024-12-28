# Eigenfunction plotting script
set terminal pdfcairo enhanced color font 'Arial,12' size 8,6

# Function to generate filename
filename(k,ev) = sprintf("../output/eigenfunctions_k_%05d_ev_%05d.dat", k, ev)

# Plot settings
set xlabel 'Component'
set ylabel 'Amplitude'
set grid
set key outside right

# Style settings
set style line 1 lt 1 lw 2 pt 7 ps 0.5 lc rgb '#0060ad'  # Blue (Real part)
set style line 2 lt 1 lw 2 pt 7 ps 0.5 lc rgb '#dd181f'  # Red (Imaginary part)

# For bulk calculations (8 components per state)
set title 'Eigenfunction Components at k-point 1'
set output '../output/eigenfunctions_k1.pdf'
set multiplot layout 4,2 title 'Eigenfunctions at k-point 1'
do for [i=1:8] {
    set title sprintf('State %d', i)
    plot filename(1,i) using 0:1 title 'Re' w lp ls 1, \
         '' using 0:2 title 'Im' w lp ls 2
}
unset multiplot

# For k-point in the middle of the Brillouin zone
set output '../output/eigenfunctions_k_mid.pdf'
set multiplot layout 4,2 title 'Eigenfunctions at middle k-point'
do for [i=1:8] {
    set title sprintf('State %d', i)
    plot filename(6,i) using 0:1 title 'Re' w lp ls 1, \
         '' using 0:2 title 'Im' w lp ls 2
}
unset multiplot

# For last k-point
set output '../output/eigenfunctions_k_end.pdf'
set multiplot layout 4,2 title 'Eigenfunctions at last k-point'
do for [i=1:8] {
    set title sprintf('State %d', i)
    plot filename(11,i) using 0:1 title 'Re' w lp ls 1, \
         '' using 0:2 title 'Im' w lp ls 2
}
unset multiplot

# Reset terminal
set output
set terminal pop 