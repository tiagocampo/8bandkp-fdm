# Band structure plotting script
set terminal pdfcairo enhanced color font 'Arial,12' size 8,6

# Output file
set output '../output/band_structure.pdf'

# Plot settings
set title 'Band Structure'
set xlabel 'k (1/Å)'
set ylabel 'Energy (eV)'
set grid

# Style settings
set style data lines
set style line 1 lw 2 lc rgb '#0060ad'  # Blue

# Plot data
plot '../output/eigenvalues.dat' using 1:2 title 'Band 1' ls 1, \
     for [i=3:*] '../output/eigenvalues.dat' using 1:i title sprintf('Band %d',i-1) ls 1

# Reset terminal
set output
set terminal pop 