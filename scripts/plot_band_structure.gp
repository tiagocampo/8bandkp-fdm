# Gnuplot script for band structure visualization
# Usage: cd <output_directory> && gnuplot ../../scripts/plot_band_structure.gp

set terminal pngcairo size 1000,700 enhanced font 'Arial,12'
set output 'band_structure.png'

set title "8-Band k·p Band Structure" font 'Arial,14'
set xlabel "k-point index"
set ylabel "Energy (eV)"
set grid

# Style settings
set style line 1 lc rgb '#0060ad' lt 1 lw 2 # CB
set style line 2 lc rgb '#dd181f' lt 1 lw 2 # VB

# The eigenvalues.dat file has format: k_val eig1 eig2 ... eig8
# Columns 2-9 are the 8 eigenvalues
# Typically: columns 2-7 are VB, columns 8-9 are CB

plot 'eigenvalues.dat' using 1:2 with lines ls 2 title "VB1", \
     '' using 1:3 with lines ls 2 title "VB2", \
     '' using 1:4 with lines ls 2 title "VB3", \
     '' using 1:5 with lines ls 2 title "VB4", \
     '' using 1:6 with lines ls 2 title "VB5", \
     '' using 1:7 with lines ls 2 title "VB6", \
     '' using 1:8 with lines ls 1 title "CB1", \
     '' using 1:9 with lines ls 1 title "CB2"

print "Band structure plot saved to band_structure.png"
