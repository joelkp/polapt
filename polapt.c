/* Polynomial approximation optimizer.
 * Copyright (c) 2021 Joel K. Pettersson
 * <joelkpettersson@gmail.com>.
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include <stdint.h>
#include <stdio.h>
#include <math.h>
#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

/*
 * Some polynomials for testing purposes.
 */

/* max error: 1.388788e-05
   (tweaked manually using graph plot program) */
static inline float SGS_sinf_t7_old(float x) {
	const float scale7 = -1.f/5040 * 29.f/30;
	float x2 = x*x;
	return x + x*x2*(-1.f/6 + x2*(1.f/120 + x2*scale7));
}

/* max error: 3.576279e-07
   (tweaked manually using graph plot program) */
static inline float SGS_sinf_t9_old(float x) {
	const float scale9 = 1.f/362880 * 44.f/45;
	float x2 = x*x;
	return x + x*x2*(-1.f/6 + x2*(1.f/120 + x2*(-1.f/5040 + x2*scale9)));
}

/* max error: 1.788139e-07
   https://web.archive.org/web/20200628195036/http://mooooo.ooo/chebyshev-sine-approximation/ */
static inline float moo_sine(float x) {
	const float coeffs[] = {
		-0.10132118f,          // x
		+0.0066208798f,        // x^3
		-0.00017350505f,       // x^5
		+0.0000025222919f,     // x^7
		-0.000000023317787f,   // x^9
		+0.00000000013291342f, // x^11
	};
	float pi_major = 3.1415927f;
	float pi_minor = -0.00000008742278f;
	float x2 = x*x;
	float p11 = coeffs[5];
	float p9  = p11*x2 + coeffs[4];
	float p7  = p9*x2  + coeffs[3];
	float p5  = p7*x2  + coeffs[2];
	float p3  = p5*x2  + coeffs[1];
	float p1  = p3*x2  + coeffs[0];
	return (x - pi_major - pi_minor) * (x + pi_major + pi_minor) * p1 * x;
}

/*
 * The program.
 */

/* What to run? */
#define TEST_Y test_sin
#define GOOD_Y sinf
#define TEST_T float

/* Test more than starting point? */
#define RUN_TESTS       1

/* Produce file suitable for gnuplot? */
#define WRITE_PLOT_FILE 1

static inline TEST_T test_sin(TEST_T x, double scale_adj[]) {
	const TEST_T scale[] = {
		-1.f/6 * scale_adj[0],
		+1.f/120 * scale_adj[1],
		-1.f/5040 * scale_adj[2],
	};
	TEST_T x2 = x*x;
	return x + x*x2*(scale[0] + x2*(scale[1] + x2*scale[2]));
}

#define TAB_LEN 1000 //64 //16 //128 //1024
#define SUB_LEN 10
#define MAX_ERR 1.f  // large-enough start value to accept any contender
#define EPSILON (1.f/SUB_LEN)
#define PDIM 3

TEST_T good_y[TAB_LEN];
double tryerr_y[TAB_LEN];
double selerr_y[TAB_LEN];
double minmaxerr_y = MAX_ERR;
double trymaxerr_y = MAX_ERR;

double scale_adj[PDIM] = {1.f, 1.f, 1.f};
double tryscale_adj[PDIM];
double selscale_adj[PDIM];
double trypos[PDIM];
double selpos[PDIM];

/*
 * For minimizing error at end of curve only.
 */
static int compare_enderr(double minerr) {
	uint32_t i = TAB_LEN - 1;
	double x = (1.f - 0.5f);
	double err = TEST_Y(x * M_PI, tryscale_adj) - good_y[i];
	double abserr = fabs(err);
	minerr = fabs(minerr);
	tryerr_y[i] = err;
	if (abserr > trymaxerr_y)
		trymaxerr_y = abserr;
	if (abserr >= minerr)
		return 0 - (abserr > minerr);
	return 1;
}

/*
 * For minimizing error all over curve.
 */
static int compare_maxerr(double minerr) {
	minerr = fabs(minerr);
	for (uint32_t i = 0, end = TAB_LEN - 1; i <= end; ++i) {
		double x = (i * 1.f/end - 0.5f);
		double err = TEST_Y(x * M_PI, tryscale_adj) - good_y[i];
		double abserr = fabs(err);
		tryerr_y[i] = err;
		if (abserr > trymaxerr_y)
			trymaxerr_y = abserr;
		if (abserr > minerr)
			return -1;
	}
	return 1;
}

/*
 * For minimizing error all over curve.
 *
 * Uses compare_enderr() as tie-breaker
 * when the threshold is reached but not exceeded.
 */
static int compare_maxerr_enderr(double minerr) {
	minerr = fabs(minerr);
	for (uint32_t i = 0, end = TAB_LEN - 1; i <= end; ++i) {
		double x = (i * 1.f/end - 0.5f);
		double err = TEST_Y(x * M_PI, tryscale_adj) - good_y[i];
		double abserr = fabs(err);
		tryerr_y[i] = err;
		if (abserr > trymaxerr_y)
			trymaxerr_y = abserr;
		if (abserr >= minerr) {
			if (abserr == minerr) {
				int cmp = compare_enderr(selerr_y[end]);
				if (cmp > 0)
					continue;
				return cmp;
			}
			return 0 - (abserr > minerr);
		}
	}
	return 1;
}

/*
 * For minimizing error at the end of the curve primarily,
 * over the rest of the curve secondarily.
 */
static int compare_enderr_maxerr(double minerr) {
	if (compare_enderr(selerr_y[TAB_LEN - 1]) < 0)
		return -1;
	return compare_maxerr(minerr);
}

static inline void set_candidate(uint32_t n, double pos) {
	trypos[n] = pos;
}

static int try_candidate(int (*compare)(double minerr), double minerr) {
	trymaxerr_y = 0.f;
	for (uint32_t j = 0; j < PDIM; ++j) {
		double q = trypos[j];
		tryscale_adj[j] = scale_adj[j];
		if (q > 0)
			tryscale_adj[j] *= (q - 1.f) / q;
	}
	return compare(minerr);
}

static void select_candidate(void) {
	for (uint32_t i = 0, end = TAB_LEN; i < end; ++i) {
		selerr_y[i] = tryerr_y[i];
	}
	if (trymaxerr_y < minmaxerr_y)  {
		minmaxerr_y = trymaxerr_y;
	}
	for (uint32_t j = 0; j < PDIM; ++j) {
		selpos[j] = trypos[j];
		selscale_adj[j] = tryscale_adj[j];
	}
}

static void print_report(void) {
	printf("Max.err. %.11e\tEnd.err. %.11e\n",
			minmaxerr_y, selerr_y[TAB_LEN - 1]);
	double selscale_report[PDIM][2];
	for (uint32_t j = 0; j < PDIM; ++j) {
		double q = selpos[j];
		if (q <= 0) {
			selscale_report[j][0] = 0.f;
			selscale_report[j][1] = 0.f;
		} else {
			selscale_report[j][0] = q - 1.f;
			selscale_report[j][1] = q;
		}
	}
	for (uint32_t j = 0; j < PDIM; ++j) {
		char label = 'A' + j;
		printf("%c==%.11e\t(%.11e * (%.1f/%.1f))\n",
				label, selscale_adj[j], scale_adj[j],
				selscale_report[j][0], selscale_report[j][1]);
	}
}

/*
 * Makes result of previous pass apply to next pass.
 */
static void apply_selected(void) {
	printf("(Applying selected coefficients.)\n");
	for (uint32_t j = 0; j < PDIM; ++j) {
		scale_adj[j] = selscale_adj[j];
	}
}

static const uint32_t loop_limits[PDIM] = {
	100000, //100000
	1000,   //10000
	100,    //1000
};

uint32_t bench_count;
static inline double probe_posi(uint32_t n, uint32_t pos, double minerr) {
	double err;
	++bench_count;
	set_candidate(n, pos);
	try_candidate(compare_maxerr, minerr);
	err = trymaxerr_y;
	return err;
}
static inline double probe_at(uint32_t n, uint32_t pos) {
	double err;
	++bench_count;
	set_candidate(n, pos);
	try_candidate(compare_maxerr, minmaxerr_y);
	err = trymaxerr_y;
	return err;
}
static inline double probe_above(uint32_t n, uint32_t pos) {
	double err;
	++bench_count;
	set_candidate(n, pos * 2.f);
	try_candidate(compare_maxerr, minmaxerr_y);
	err = trymaxerr_y;
	return err;
}
static inline double probe_below(uint32_t n, uint32_t pos) {
	double err;
	++bench_count;
	set_candidate(n, pos * .5f);
	try_candidate(compare_maxerr, minmaxerr_y);
	err = trymaxerr_y;
	return err;
}

static inline int probe_posf(uint32_t n, double pos, double minerr) {
	int found = 0;
	for (uint32_t i = 0, end = SUB_LEN - 1; i <= end; ++i) {
		double subpos = pos - (i * 1.f/end - 0.5f);
		set_candidate(n, subpos);
		if (try_candidate(compare_maxerr, minerr) > 0) {
			minerr = trymaxerr_y;
			select_candidate();
			found = 1;
		}
	}
	if (found)
		trymaxerr_y = minerr;
	return found;
}

static int test_linear(uint32_t n, uint32_t lbound, uint32_t ubound) {
	int found = 0;
	uint32_t middle = ((lbound + 1) >> 1) + (ubound >> 1);
	/*
	double le = probe_posi(n, lbound, minmaxerr_y);
	double me = probe_posi(n, middle, minmaxerr_y);
	double ue = probe_posi(n, ubound, minmaxerr_y);
	printf("le %e, me %e, ue %e | %e, %e\n",
			le, me, ue,
			(me - le), (ue - me)
	      );
	      */
	//printf("lbound==%u, ubound==%u\n", lbound, ubound);
#if 0
	for (uint32_t i = lbound - 1; i <= ubound; i += 2) {
		for (uint32_t j = 0, end = SUB_LEN << 1; j < end; ++j) {
			double pos = i + (j * 1.f/(end - 1) - 1.f);
			set_candidate(n, pos);
			if (try_candidate(compare_maxerr_enderr,
						minmaxerr_y) > 0) {
				select_candidate();
				found = 1;
			}
		}
	}
#else
	for (uint32_t i = lbound; i <= ubound; ++i) {
		double pos = i;
		set_candidate(n, pos);
		if (try_candidate(compare_maxerr_enderr,
					minmaxerr_y) > 0) {
			select_candidate();
			found = 1;
		}
	}
#endif
	if (found)
		printf("*[%f], l==%u, m==%u, u==%u\n", selpos[n], lbound, (lbound + ubound) >> 1, ubound);
	return found;
}

//static uint32_t find_middle

static int run_linear(uint32_t n) {
	uint32_t lbound, mpivot, ubound;
	double lerr, merr, uerr;
	double mdecerr, mincerr;
	/*
	 * Find upper and lower bound for search.
	 * (0 has the effect of infinite number.)
	 *
	 * First narrow it down to powers of two.
	 */
	uint8_t lpow = 0, mpow = 16, upow = 32; /* treat 1<<32 as above 1<<31 */
	lbound = 1; /* treated as smallest value */
	mpivot = 1<<16;
	ubound = 0; /* treated as UINT32_MAX + 1 */
	mincerr = probe_above(n, mpivot);
	merr = probe_at(n, mpivot);
	mdecerr = probe_below(n, mpivot);
	for (;;) {
		if (mdecerr < mincerr) {
			if (mdecerr < merr) {
				/* - - 0 */
				lerr = probe_at(n, lbound);
			GO_DOWN:
				upow = mpow;
				ubound = mpivot;
				mpow = (lpow + mpow + 1) >> 1;
				mpivot = 1<<mpow;
				if (mpivot == lbound) {
					merr = lerr;
					break;
				}
				mdecerr = probe_below(n, mpivot);
				merr = probe_at(n, mpivot);
				if (mdecerr <= merr) goto GO_DOWN;
				mincerr = probe_above(n, mpivot);
//				printf("-- l1<<%u==%u, m1<<%u==%u\n", lpow, lbound, mpow, mpivot);
			} else {
//				if (mdecerr == merr) puts("???");
				/* - + 0 */
//				printf("? {- + 0}\t(%u +/- 1)\n", mpivot);
				break; /* a middle found */
			}
		} else {
			if (merr <= mincerr) {
				/* + - 0 */
//				printf("! {+ - 0}\t(%u +/- 1)\n", mpivot);
				break; /* a middle found */
			} else {
//				if (merr == mincerr) puts("!!!");
				/* + + 0 */
				uerr = probe_at(n, ubound);
			GO_UP:
				lpow = mpow;
				lbound = mpivot;
				mpow = (upow + mpow + 1) >> 1;
				mpivot = (mpow == 32) ? 0 : 1<<mpow;
				if (mpivot == ubound) {
					merr = uerr;
					break;
				}
				mincerr = probe_above(n, mpivot);
				merr = probe_at(n, mpivot);
				if (mincerr <= merr) goto GO_UP;
				mdecerr = probe_below(n, mpivot);
//				printf("++ u1<<%u==%u, m1<<%u==%u\n", upow, ubound, mpow, mpivot);
			}
		}
	RETEST:
		// printf(", {l==%u, m==%u, u==%u}\n", lbound, mpivot, ubound);
		if ((ubound - 1) <= mpivot) break;
		if ((lbound + 1) >= mpivot) break;
	}
//	printf(". {l==%u, m==%u, u==%u}\n", lbound, mpivot, ubound);
	if (lbound > mpivot) /* happens if there's only up-going */
		lbound = mpivot;
	if (ubound > mpivot && ubound < (mpivot << 1)) /* ??? */
		ubound = mpivot << 1;
#if 0 //.
		if (probe_posi(n, mpivot - 1, minmaxerr_y) < merr)

		if (uerr < lerr) {
			if (uerr < merr) {
				if (ubound == 0) {
					/* sloppy but almost always right */
					lbound = ubound;
					goto SEARCH;
				}
				lbound = mpivot;
			} else {
				goto TO_M;
			}
		} else {
			if (merr < lerr) {
				goto TO_M;
			} else {
				if (lbound == 1) {
					/* sloppy but almost always right */
					ubound = lbound;
					goto SEARCH;
				}
				ubound = mpivot;
			}
		}
		mpivot = ((lbound + 1) >> 1) + ((ubound - 1) >> 1) + 1;
		if (ubound - mpivot < 1) break;
		if (lbound - mpivot < 1) break;
//		printf("2{l==%u, m==%u, u==%u}\n", lbound, mpivot, ubound);
		continue;
	TO_M:
		lbound = ((lbound + 1) >> 1) + (mpivot >> 1);
		ubound = ((ubound - 1) >> 1) + ((mpivot + 1) >> 1) + 1;
//		printf("{l==%u, m==%u, u==%u}\n", lbound, mpivot, ubound);
		if (ubound - mpivot < 1) break;
		if (mpivot - lbound < 1) break;
	}
#endif //.
#if 0
	/*
	 * Find upper and lower bound for search.
	 */
	uint32_t lbound, mpivot, ubound;
	double lerr, merr, uerr;
	lbound = 1; /* treated as smallest value */
	ubound = 0; /* treated as UINT32_MAX + 1 */
	/* move lbound up by powers of two until
	   the error curve gets increasing slope */
	do {
		lerr = probe_posi(n, lbound, minmaxerr_y);
		if (lerr <= probe_posi(n, lbound + 1, minmaxerr_y)) {
			if (lerr == probe_posi(n, lbound - 1, minmaxerr_y))
				lbound >>= 1;
			break;
		}
		lbound <<= 1;
	} while (lbound != ubound);
	if (lbound <= 2) {
		ubound = (lerr == trymaxerr_y) ? lbound + 1 : lbound;
		goto SEARCH;
	}
	lbound >>= 1;
	/* move ubound down powers of two, until
	   the error curve gets decreasing slope */
	do {
		uerr = probe_posi(n, ubound, minmaxerr_y);
		if (probe_posi(n, ubound - 1, minmaxerr_y) > uerr) {
			if (uerr == probe_posi(n, ubound + 1, minmaxerr_y))
				ubound <<= 1;
			break;
		}
		ubound = ((ubound - 1) >> 1) + 1;
	} while (ubound > lbound);
	if (ubound < lbound) {
		if (uerr <= lerr) {
			lbound = ubound;
			goto SEARCH;
		} else {
			puts("AAAA");
		//	ubound = lbound << 3;
			goto SEARCH;
		}
	}
	ubound <<= 1;
	if (ubound < lbound)
		printf("lbound==%u, ubound==%u\n", lbound, ubound);
	/*
	if (ubound == lbound)
		goto SEARCH;
	lerr = probe_posi(n, lbound, minmaxerr_y);
	for (;;) {
		mpivot = ((lbound + 1) >> 1) + ((ubound - 1) >> 1) + 1;
		merr = probe_posi(n, mpivot, minmaxerr_y);
	}
	*/
#elif 0
	uint32_t i, j, lbound, ubound;
	double ierr;
	/*
	 * Find upper and lower bound for search.
	 *
	 * Move down by powers of two from a very
	 * large start range until slope changes.
	 * (0 has the effect of infinite number.)
	 */
	i = lbound = ubound = 0; /* treat 0 as UINT32_MAX + 1 */
	do {
		ierr = probe_posi(n, i, minmaxerr_y);
		if (ierr < probe_posi(n, i - 1, minmaxerr_y)) {
			lbound = i;
			break;
		}
		ubound = i;
		i = ((i - 1) >> 1) + 1;
		if (ierr == trymaxerr_y)
			ubound -= i >> 1; /* optimization: don't double */
	} while (i > 1);
	if (i == 1) /* can't handle usefully */
		return 0;
	if (ierr == probe_posi(n, i + 1, minmaxerr_y))
		ubound += i; /* optimization: don't double */
	/*
	 * Move up by powers of two from lower to
	 * higher bound, until the slope changes.
	 */
	j = 2;
	i = lbound + j;
	do {
		ierr = probe_posi(n, i - 1, minmaxerr_y);
		if (ierr < probe_posi(n, i, minmaxerr_y)) {
			ubound = i;
			break;
		}
		j <<= 1;
		i = lbound + j;
	} while (i < ubound);
#endif
SEARCH:
//	printf("L==%u, U==%u\n", lbound, ubound);
	return test_linear(n, lbound, ubound);
}

static int recurse_linear(uint32_t j, uint32_t n) {
	int found = 0;
	if (j == 0) {
		found = (run_linear(PDIM - n) > 0);
	} else {
		const uint32_t limit = loop_limits[PDIM - n + j];
		for (uint32_t i = 0; i <= limit; ++i) {
			set_candidate(PDIM - n + j, i);
			if (recurse_linear(j - 1, n)) found = 1;
		}
	}
	return found;
}

/*
 * Run n-dimensional pass, from 0 to PDIM;
 * 0-pass tests the unmodified polynomial.
 */
static int run_pass(uint32_t n) {
	int found = 0;
	for (uint32_t j = 0; j < PDIM; ++j)
		set_candidate(j, 0);
	if (n == 0) {
		if (try_candidate(compare_maxerr, minmaxerr_y) > 0) {
			select_candidate();
			found = 1;
		}
	} else {
		found = (recurse_linear(n - 1, n) > 0);
	}
	if (found)
		print_report();
	return found;
}

int main(void) {
#if WRITE_PLOT_FILE
	FILE *f = fopen("plot.txt", "w");
#endif
	for (uint32_t i = 0, end = TAB_LEN - 1; i <= end; ++i) {
		double x = (i * 1.f/end - 0.5f);
		good_y[i] = GOOD_Y(x * M_PI);
		selerr_y[i] = MAX_ERR;
	}
	run_pass(0); /* also print stats for unmodified polynomial */
#if RUN_TESTS
	for (uint32_t n = 1; n <= PDIM; ++n) {
		run_pass(n);
		if (n < PDIM)
			apply_selected();
	}
#endif
#if WRITE_PLOT_FILE
	for (uint32_t i = 0, end = TAB_LEN - 1; i <= end; ++i) {
		double x = (i * 1.f/end - 0.5f);
		fprintf(f, "%.11f\t%.11f\n", x, selerr_y[i]);
	}
	fclose(f);
#endif
	printf("BENCH %u\n", bench_count);
	return 0;
}
