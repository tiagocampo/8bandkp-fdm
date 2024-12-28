# Potential profile plotting script
set terminal pdfcairo enhanced color font 'Arial,12' size 8,6

# Output file
set output '../output/potential_profile.pdf'

# Plot settings
set title 'Potential Profile'
set xlabel 'Position (Å)'
set ylabel 'Energy (eV)'
set grid

# Style settings
set style line 1 lt 1 lw 2 pt 7 ps 0.5 lc rgb '#0060ad'  # Blue
set style line 2 lt 1 lw 2 pt 7 ps 0.5 lc rgb '#dd181f'  # Red
set style line 3 lt 1 lw 2 pt 7 ps 0.5 lc rgb '#00A000'  # Green

# Plot data
plot '../output/potential_profile.dat' using 1:2 title 'CB' w l ls 1, \
     '' using 1:3 title 'HH/LH' w l ls 2, \
     '' using 1:4 title 'SO' w l ls 3

# Reset terminal
set output
set terminal pop 