Polynomial approximation optimizer
==================================

This program compares iterated modifications
of a polynomial approximation of a function, with the function,
selecting the approximation which best minimizes maximum error.
Results are good but not perfect for the kind of approximations
it was designed for.

The default function is libm sinf(x) and the starting approximation
a Taylor polynomial, of degree 7 (4 terms). By default, the program
tries to minimize maximum error for -PI/2 <= x <= PI/2,
focusing only on that input range. (Taylor polynomials can be made
much more accurate for a given computational cost by limiting input
to the smallest range practical to reduce to.)

Each time a search for modifying coefficients for the approximation
completes, they are printed (for easy use in a program), along with
the maximum and end-point error when using them. A few searches are
made, each using the best pick of the previous as a starting point.
When searches are complete, a full picture of the error (difference
curve) for the final resulting approximation is written as a table,
to a text file for plotting with e.g. gnuplot.

If you want to test something else, currently
you'll need to edit the polapt.c file.

See the article "[Modifying Taylor polynomials for better accuracy](https://joelkp.frama.io/blog/modified-taylor.html)"
for background and the results for a Taylor polynomial of degree 7 (4 terms).

Plotting of error for improved (error-reduced) polynomial
---------------------------------------------------------

```
make
```

```
./polapt
```

```
gnuplot
plot 'plot.txt' w l
exit
```
