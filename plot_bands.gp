set title "Band Structure"
set xlabel "kx"
set ylabel "Eigenvalue"
set grid

set terminal pngcairo
set output "band_structure.png"

plot for [i=2:*] 'eigenvalues.dat' using 1:i with lines title columnheader(i)