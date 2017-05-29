set term postscript eps color enhanced "Helvetica" 22 #fontfile 'cmmi10.pfb'
set xtics font "Helvetica, 23"
set ytics font "Helvetica, 23"

set output "beta_degen.eps"

set log xy
#set log x
#set xrange [0.01:3]
#set yrange [1:2]

set ytics 

set title "|Beta| wrt k/kstar for different phi_i showing degeneracy for UV modes"
set xlabel 'k/k_*'
set ylabel '|Beta|'

kstar=1.17583

set arrow from 0.1,2.8e-11 to 0.1,3.4e-11 nohead

set arrow from 9,2.8e-11 to 9,3.4e-11 nohead

#p 'poweriter_1_1p171_scalar_mid.dat' u 2:(sqrt($20**2+$21**2)) w lp, 'poweriter_2_1p171_scalar_mid.dat' u 2:(sqrt($20**2+$21**2)) w lp, 'poweriter_3_1p171_scalar_mid.dat'  u 2:(sqrt($20**2+$21**2)) w l, 'poweriter_4_1p171_scalar_mid.dat' u 2:(sqrt($20**2+$21**2)) w l, 'poweriter_5_1p171_scalar_mid.dat' u 2:(sqrt($20**2+$21**2)) w l, 'poweriter_6_1p171_scalar_mid.dat' u 2:(sqrt($20**2+$21**2)) w l, 'poweriter_7_1p171_scalar_mid.dat' u 2:(sqrt($20**2+$21**2)) w l

p 'poweriter_1_1p028_scalar_mid.dat' u ($2/kstar):(sqrt($20**2+$21**2)) w lp ti 'iter=1 phib =1.1028', \
'../Backreaction_tini30000_phib1171_GaussBeta/poweriter_10_1p171_scalar_mid.dat' u ($2/kstar):(sqrt($20**2+$21**2)) w l ti 'iter=10 phib=1.171' 
