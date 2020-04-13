set term "png" enhanced
set output "Scaling.png"
set border 31 lt -1 lw  3.000
set origin 0.02,0.02
set size 0.98,0.98
set key font "Times-Roman,15"
set key spacing 1.2
set key center at 0.06,-0.56
set yrange [-0.85:-0.5]
set xtics  font "Times-Roman,20"
set title "E=-0.670076" font "Times-Roman,20"
set ytics  font "Times-Roman,20"
set ylabel "E" font "Times-Roman,22"
set xlabel "{/Symbol \326}(L_1*L_2)" #font "Times-Roman,22"
p "fixData" u 1:2 w l lt 1 lw 4 t "","input.dat" u (1/$1):2 w p lt 6 ps 3 lw 4 t "",-0.669275 t "DMRG",-0.6669 t "iPEPS",-0.66191 t "VQMC","fixData1" u 1:2 w l lt 1 lw 4 t "","input.dat1" u (1/$1):2 w p lt 6 ps 3 lw 4 t ""

