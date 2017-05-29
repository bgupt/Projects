set term postscript eps color enhanced "Helvetica" 22 #fontfile 'cmmi10.pfb'
set xtics font "Helvetica, 23"
set ytics font "Helvetica, 23"

set output "tensorPower_1p17.eps"

#set log xy
#set log x
set xrange [0.01:3]
#set yrange [1:2]

set ytics 

set title "Tensor power spectrum"
set xlabel 'k/k_*'
set ylabel 'P_T'

p 'poweriter_1_1p17_v2.dat' u ($2/17.716):($2**3*($3**2+$4**2)*32*pi/(2*pi**2)) w l lw 2 ti 'Iteration=1',\
  'poweriter_10_1p17_v2.dat' u ($2/1.93288):($2**3*($3**2+$4**2)*32*pi/(2*pi**2)) w l ti 'Iteration=10'

set output "scalefactor_1p17.eps"

set log y
set xrange [:1.5e6]

set title "scale factor"
set xlabel "t"

p 'data_phib117_tensor/background_iter_001.dat' u 1:2 w l lw 1.5 ti "iter=1", 'data_phib117_tensor/background_iter_002.dat' u 1:2 w l ti "iter=2", 'data_phib117_tensor/background_iter_003.dat' u 1:2 w l ti "iter=3", 'data_phib117_tensor/background_iter_004.dat' u 1:2 w l ti "iter=4", 'data_phib117_tensor/background_iter_006.dat' u 1:2 w l ti "iter=6", 'data_phib117_tensor/background_iter_008.dat' u 1:2 w l ti "iter=8", 'data_phib117_tensor/background_iter_010.dat' u 1:2 w l lt 0 lw 2 lc rgb 'black' ti "iter=10"


set output "Hubble_1p17.eps"

unset log xy
set xrange [:1.5e6]

set title "Hubble rate"
set xlabel "t"

p 'data_phib117_tensor/background_iter_001.dat' u 1:($3/$2) w l lw 1.5 ti "iter=1", 'data_phib117_tensor/background_iter_002.dat' u 1:($3/$2) w l ti "iter=2", 'data_phib117_tensor/background_iter_003.dat' u 1:($3/$2) w l ti "iter=3", 'data_phib117_tensor/background_iter_004.dat' u 1:($3/$2) w l ti "iter=4", 'data_phib117_tensor/background_iter_006.dat' u 1:($3/$2) w l ti "iter=6", 'data_phib117_tensor/background_iter_008.dat' u 1:($3/$2) w l ti "iter=8", 'data_phib117_tensor/background_iter_010.dat' u 1:($3/$2) w l lt 0 lw 2 lc rgb 'black' ti "iter=10", 7.83e-6


set output "Hubble_1p17_zoom.eps"

unset log xy
set xrange [:1.5e6]
set yrange [7e-6:9e-6]

set title "Hubble rate"
set xlabel "t"

p 'data_phib117_tensor/background_iter_001.dat' u 1:($3/$2) w l lw 1.5 ti "iter=1", 'data_phib117_tensor/background_iter_002.dat' u 1:($3/$2) w l ti "iter=2", 'data_phib117_tensor/background_iter_003.dat' u 1:($3/$2) w l ti "iter=3", 'data_phib117_tensor/background_iter_004.dat' u 1:($3/$2) w l ti "iter=4", 'data_phib117_tensor/background_iter_006.dat' u 1:($3/$2) w l ti "iter=6", 'data_phib117_tensor/background_iter_008.dat' u 1:($3/$2) w l ti "iter=8", 'data_phib117_tensor/background_iter_010.dat' u 1:($3/$2) w l lt 0 lw 2 lc rgb 'black' ti "iter=10", 7.83e-6
