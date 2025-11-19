# Gnuplot script for g-factor visualization
# Usage: cd <output_directory> && gnuplot ../../scripts/plot_gfactor.gp

set terminal pngcairo size 900,600 enhanced font 'Arial,12'
set output 'gfactor.png'

set title "G-Factor Components" font 'Arial,14'
set ylabel "g-factor value"
set grid ytics
set style fill solid 0.8 border -1
set boxwidth 0.6
set xtics ("g_x" 1, "g_y" 2, "g_z" 3)
set xrange [0.5:3.5]

# Read gfactor values from file
set arrow from 0,0 to 4,0 nohead lc rgb 'black' lw 1

# Plot the data directly - gfactor.dat has format: gx gy gz (single line)
plot "< awk '{print 1, $1; print 2, $2; print 3, $3}' gfactor.dat" using 1:2 with boxes lc rgb '#3366cc' title sprintf("Values")

print "G-factor plot saved to gfactor.png"
