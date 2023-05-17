# to set filename, call gnuplot with
# -e "filename='yourfilename.dat'"
# before the name of this script.
if (!exists("filename")) filename='profile_mpi.dat'

set term qt font "arial,14"
set term qt size 600,450
set pointsize 2
f(x) = a * (x ** 2) * log(x) + c
c = 1e-5
fit f(x) filename using 1:8 via a,c
fstr = sprintf('x^2 log(x)')
# g(x) = b * (x ** 3) + d * (x ** 2) + e * x + f
# fit g(x) filename using 1:8 via b,d,e,f
g(x) = b * (x ** 3) + d
fit g(x) filename using 1:8 via b,d
gstr = sprintf('x^3')
h(x) = g * (x ** h)
h = 2.5
g = b
fit h(x) filename using 1:8 via g,h
hstr = sprintf('%.2g x^{%0.3g}', g, h)
set key top left
set xlabel 'domain size'
set ylabel 'Clustering wall-time [s]'
plot filename using 1:8 title "Serial" lw 1.5, f(x) title fstr lw 1.5, g(x) title gstr lw 1.5, h(x) title hstr lw 1.5
pause -1 "Hit any key to continue"