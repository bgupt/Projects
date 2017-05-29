set term postscript eps color enhanced "Helvetica" 22 #fontfile 'cmmi10.pfb'
set xtics font "Helvetica, 23"
set ytics font "Helvetica, 23"

set output "scalarPower_1p17.eps"

#set log xy
#set log x
set xrange [0.01:3]
#set yrange [1:2]

set ytics 

set title "Scalar power spectrum"
set xlabel 'k/k_*'
set ylabel 'P_R'

p 'poweriter_1_1p17_scalar.dat' u ($2/17.716):($2**3*($3**2+$4**2)*4*3.14/(1-$9*$7/$8**2)/(2*3.14**2)) w l lw 2 ti 'Iteration=1', 'poweriter_10_1p17_scalar.dat' u ($2/2.528):($2**3*($3**2+$4**2)*4*3.14/(1-$9*$7/$8**2)/(2*3.14**2)) w l ti 'Iteration=10',ZZ\
2.474e-9*x**(0.9684-1) w l ti 'BD'

