# G-factor plotting script
set terminal pdfcairo enhanced color font 'Arial,12' size 8,6

# Output file
set output '../output/gfactor.pdf'

# Plot settings
set title 'g-factor Components'
set style data histogram
set style histogram cluster gap 1
set style fill solid
set boxwidth 0.9
set xtics format ""
set grid y
set ylabel 'g-factor'

# Style settings
set style line 1 lc rgb '#0060ad'  # Blue
set style line 2 lc rgb '#dd181f'  # Red
set style line 3 lc rgb '#00A000'  # Green

# Create temporary data file with labels
set table 'temp.dat'
plot '../output/gfactor.dat' using 1
unset table

# Plot data
plot '../output/gfactor.dat' using 1:xtic("g_x") title 'g_x' ls 1, \
     '' using 2:xtic("g_y") title 'g_y' ls 2, \
     '' using 3:xtic("g_z") title 'g_z' ls 3

# Reset terminal
set output
set terminal pop 