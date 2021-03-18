# f(x) = a1 * (x ** e1) * log( x ** e2 )
# a1 = 1e-5
# e1 = 3
# e2 = 2
# fit f(x) 'profile.dat' using 1:5 via a1,e1,e2
# g(x) = a2 * (x ** e3)
# a2 = 1e-5
# e3 = 5
# fit g(x) 'profile.dat' using 1:5 via a2,e3
set term qt font "arial,14"
set term qt size 600,450
set pointsize 2
# f(x) = a * (x ** b)
f(x) = a * (x ** 2) * log(x)
a = 1e-7
fit f(x) 'profile.dat' using 1:5 via a
# fstr = sprintf('x^{%.2g}', b)
fstr = "x^2 log(x)"
g(x) = c * (x ** d)
fit g(x) '../serial/profile.dat' using 1:5 via c,d
gstr = sprintf('x^{%.2g}', d)
# h(x) = e * (x ** f)
h(x) = e * (x ** 2) * log(x)
fit h(x) 'profile_serial.dat' using 1:5 via e
# hstr = sprintf('x^{%.2g}', f)
hstr = sprintf('x^2 log(x)')
i(x) = m * (x ** 2) * log(x)
fit i(x) 'profile_mpi.dat' using 1:5 via m
# hstr = sprintf('x^{%.2g}', f)
istr = sprintf('x^2 log(x)')
set key top left
# set yrange [0:2000]
set xlabel 'domain size'
set ylabel 'Clustering wall-time [s]'
plot 'profile.dat' using 1:5 title "Colin's OMP" lw 1.5, f(x) title fstr lw 1.5, '../serial/profile.dat' using 1:5 title "Serial" lw 1.5, g(x) title gstr lw 1.5, 'profile_serial.dat' using 1:5 title "Colin's Serial" lw 1.5, h(x) title hstr lw 1.5, 'profile_mpi.dat' using 1:5 title "MPI" lw 1.5, i(x) title istr lw 1.5
print sprintf(' OMP speedup = %.2g', e / a )
pause -1 "Hit any key to continue"