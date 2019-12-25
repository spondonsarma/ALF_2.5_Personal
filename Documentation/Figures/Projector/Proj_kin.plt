set terminal epslatex color standalone
set output "Proj_kin.tex"
set key center top
#set key at graph 0.7,0.99
set multiplot
#set origin 0.2,0.5
set size 0.9,0.9
#set key at graph 0.5,-0.5
set pointsize 1.5
#set ytics -2.95,0.01, -2.85
#set yrange [-80.5:-79.5]
set xrange [0.0:0.21] 
set ylabel "$ \\langle \\hat{H}_{t}  \\rangle $"
set title "$L=6,  U/t=4$ "
set xlabel " $ \\frac{1}{\\beta t},  \\frac{1}{2 \\theta t + \\beta t} $ "
#set x2label " $ \\frac{1}{2 \\theta t + \\beta t} $ "
plot "Finite_T/Kin.dat" u (1/$1):2:3 w e lc rgb "black" lt 1 pt 8 t "Finite T", \
     "Proj_tt1/Kin.dat" u (1/($1*2 + 1.0)):2:3 w e lc rgb "black" lt 1 pt 7 t "t-t'", \
     "Proj_Kekule/Kin.dat" u (1/($1*2 + 1.0)):2:3 w e lc rgb "black" lt 1 pt 6 t "Kekule"
#plot "Langevin_kin.dat" u 1:2:3 w e lc rgb "black" lt 1 pt 7  lw 2 t "", f(x)  w l lt 1 lc rgb "blue"  t ""
unset multiplot
set output 
!pdflatex Proj_kin.tex
!open Proj_kin.pdf
