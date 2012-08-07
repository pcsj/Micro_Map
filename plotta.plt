#!/gnuplot
FILE="risultatiFINALI.txt"
#set terminal postscript eps  enhanced colour solid rounded "Helvetica" 25
#set output "plot.eps"
set terminal png enhanced 15
set output "plot.png"
#unset key
set xlabel "s (m)"
#set ylabel "" 
#set xrange[0.5:60]
#set yrange[0:1.2]
plot FILE u 1:2  w lines  lt 1  lc rgb "blue" lw 2   t   'Beta x',\
FILE u 1:3   w lines  lt 1  lc rgb "brown" lw 2   t   'Beta y',\
FILE u 1:4   w lines  lt 1  lc rgb "dark-green" lw 2   t   'Alpha x',\
FILE u 1:5   w lines  lt 1  lc rgb "red" lw 2   t   'Alpha y'
#! epstopdf Fig_sp.eps && rm Fig_sp.eps
