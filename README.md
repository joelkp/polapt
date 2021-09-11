Polynomial approximation optimizer
==================================

This program compares iterated modifications
of a polynomial approximation of a function, with the function,
selecting the approximation which best minimizes maximum error.
It's only practically useful when the number of coefficients to
find optimal values for is very low, as in at most 4. Something
else (e.g. the Remez algorithm) is recommended for more general
uses. This program does not search for the form of a polynomial
and requires it to be pre-entered by editing `polapt.c`.

When the starting point has the correct form and the search for
better coefficients doesn't take too long, it's possible to get
the optimal minimax result (except for floating point precision
limits). Simple code can also be used to transform that result,
e.g. to instead minimize end-point error primarily, afterwards.

The program presents more than one solution when asked for more
than one coefficient value: one per coefficient. That's because
the program does a search first for one, then two, etc. values;
the search time increases in polynomial time with the number to
find values for, so the earlier searches take a negligible time
compared to the last.

Each time a search for modifying coefficients for the approximation
completes, they are printed (for easy use in a program), along with
the maximum and end-point error when using them.

When searches are complete, a full picture of the error (difference
curve) for the final resulting approximation is written as a table,
to a text file for plotting with e.g. gnuplot.

See the article "[Modifying Taylor polynomials for better accuracy](https://joelkp.frama.io/blog/modified-taylor.html)"
for background and results for a Taylor polynomial of degree 7 (4 terms).

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

Configuration
-------------

Various commits change which function is tested, and the input
range (i.e. function domain) and other parameters along with that.

Taylor polynomials and other approximations can be made much more
accurate for a given computational cost by limiting input to the
smallest range practical to reduce to.

For approximating sin() or sinf() (single-precision version),
as starting points mainly Taylor polynomials are provided. Changing
the coefficients of the degree 5 (3 coeffs) or degree 7 (4 coeffs)
version can produce a minimax polynomial. By default, the domain
(input range) to try to minimize maximum error for
is -PI/2 <= x <= PI/2 for sin().

Some other functions are also provided.
