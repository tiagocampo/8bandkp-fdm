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
set style line 1 lt 1 lw 2 pt 7 ps 0.5 lc rgb '#0060ad'  # Blue (CB)
set style line 2 lt 1 lw 2 pt 7 ps 0.5 lc rgb '#dd181f'  # Red (HH)
set style line 3 lt 1 lw 2 pt 7 ps 0.5 lc rgb '#00A000'  # Green (LH)
set style line 4 lt 1 lw 2 pt 7 ps 0.5 lc rgb '#9400D3'  # Purple (SO)

# Plot data
plot '../output/eigenvalues.dat' using 1:2 title 'VB1' w l ls 2, \
     '' using 1:3 title 'VB2' w l ls 2, \
     '' using 1:4 title 'VB3' w l ls 3, \
     '' using 1:5 title 'VB4' w l ls 3, \
     '' using 1:6 title 'VB5' w l ls 4, \
     '' using 1:7 title 'VB6' w l ls 4, \
     '' using 1:8 title 'CB1' w l ls 1, \
     '' using 1:9 title 'CB2' w l ls 1

# Reset terminal
set output
set terminal pop 