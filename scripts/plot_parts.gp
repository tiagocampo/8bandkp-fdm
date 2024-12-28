# Band composition plotting script
set terminal pdfcairo enhanced color font 'Arial,12' size 12,8

# Output file
set output '../output/band_composition.pdf'

# Plot settings
set title 'Band Composition Analysis'
set xlabel 'State Index'
set ylabel 'Component Contribution'
set grid

# Style settings
set style line 1 lt 1 lw 2 lc rgb '#0060ad'  # Blue
set style line 2 lt 1 lw 2 lc rgb '#dd181f'  # Red
set style line 3 lt 1 lw 2 lc rgb '#00A000'  # Green
set style line 4 lt 1 lw 2 lc rgb '#9400D3'  # Purple
set style line 5 lt 1 lw 2 lc rgb '#FF8C00'  # Dark Orange
set style line 6 lt 1 lw 2 lc rgb '#4B0082'  # Indigo
set style line 7 lt 1 lw 2 lc rgb '#FF1493'  # Deep Pink
set style line 8 lt 1 lw 2 lc rgb '#00CED1'  # Dark Turquoise

# Plot data
plot '../output/parts.dat' using 0:1 title 'Component 1' w lp ls 1, \
     '' using 0:2 title 'Component 2' w lp ls 2, \
     '' using 0:3 title 'Component 3' w lp ls 3, \
     '' using 0:4 title 'Component 4' w lp ls 4, \
     '' using 0:5 title 'Component 5' w lp ls 5, \
     '' using 0:6 title 'Component 6' w lp ls 6, \
     '' using 0:7 title 'Component 7' w lp ls 7, \
     '' using 0:8 title 'Component 8' w lp ls 8

# Reset terminal
set output
set terminal pop 