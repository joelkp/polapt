Plotting of improved (error-reduced) Taylor polynomials
=======================================================

Iterates and compares modifications of a Taylor polynomial
with libm sinf(x) and picks the best approximation arrived at.

Tries to minimize maximum error for -PI/2 <= x <= PI/2,
ignoring results outside this domain.

`
cc plot.c -lm
`

`
./a.out
`

`
gnuplot
plot 'plot.txt' w l
exit
`
