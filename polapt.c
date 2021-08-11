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
#define EPSILON 1.e-20
#define PDIM 3

TEST_T good_y[TAB_LEN];
double tryerr_y[TAB_LEN];
double selerr_y[TAB_LEN];
double minmaxerr_y = MAX_ERR;
double trymaxerr_y = MAX_ERR;
double stageminmaxerr_y[PDIM];
int stageresult[PDIM];

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

static inline void begin_stage(uint32_t n) {
	stageminmaxerr_y[n] = MAX_ERR;
	stageresult[n] = 0;
}

static inline void set_candidate(uint32_t n, double pos) {
	trypos[n] = pos;//sqrt(pos);
}

static inline void set_candidate_int(uint32_t n, uint32_t pos) {
	trypos[n] = (pos != 0) ? ((double) pos - 1 ) / ((double) pos) : 0.f;
}

static int try_candidate(int (*compare)(double minerr), double minerr) {
	trymaxerr_y = 0.f;
	for (uint32_t j = 0; j < PDIM; ++j) {
		tryscale_adj[j] = scale_adj[j] * trypos[j];
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
	for (uint32_t j = 0; j < PDIM; ++j) {
		char label = 'A' + j;
		printf("%c==%.20f\t(%.14f * %.14f)\n",
				label, selscale_adj[j],
				scale_adj[j], selpos[j]);
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

uint32_t bench_count, sub_bench_count;
static double run_one(uint32_t n, double pos) {
	++bench_count;
	++sub_bench_count;
	set_candidate(n, pos);
	if (try_candidate(compare_maxerr_enderr, stageminmaxerr_y[n]) > 0) {
		select_candidate();
		stageresult[n] = 1;
	}
	if (trymaxerr_y < stageminmaxerr_y[n])
		stageminmaxerr_y[n] = trymaxerr_y;
	return trymaxerr_y;
}

/*
 * Subdivision testing algorithm. Looks for number between 0.0 and 1.0,
 * halving the size of steps to take within the range until a number is
 * chosen. May make roughly a hundred tests, before picking the number.
 */
static double run_subdivide(uint32_t n) {
	double lpos, mpos, upos;
	double lerr, merr, uerr;
	double mlpos, mupos;
	double mlerr, muerr;
	sub_bench_count = 0;
	lerr = run_one(n, (lpos = 0.0625f));
	mlerr = lerr; mlpos = lpos;
	merr = mlerr; mpos = mlpos;
	muerr = merr; mupos = mpos;
	uerr = run_one(n, (upos = 1.f));
	/*
	 * Go up by big steps to find nice intial middle.
	 */
	while (uerr < muerr) {
		double old_mupos = mupos;
		mupos = sqrt(mupos);
		if (mupos - old_mupos < EPSILON) break;
		muerr = run_one(n, mupos);
		lerr = mlerr; lpos = mlpos;
		mlerr = merr; mlpos = mpos;
		merr = muerr; mpos = mupos; /* current best */
	}
	do {
		lerr = mlerr; lpos = mlpos;
		mlerr = merr; mlpos = mpos;
		merr = muerr; mpos = mupos; /* current best */
		double old_mupos = mupos;
		mupos = sqrt(mupos);
		if (mupos - old_mupos < EPSILON) break;
		muerr = run_one(n, mupos);
	} while (muerr < merr);
	if (mlpos != mpos) {
		lpos = mlpos; lerr = mlerr;
	} else {
		mlpos = lpos; mlerr = lerr;
	}
	uerr = muerr; upos = mupos;
SEARCH_SIDES:
	/*
	 * Search on the lower side of the middle.
	 */
	for (;;) {
		double old_mlpos = mlpos;
		mlpos = (lpos + mpos) * 0.5f;
		if (fabs(old_mlpos - mlpos) < EPSILON) break;
		mlerr = run_one(n, mlpos);
		if (merr > mlerr) {
			uerr = merr; upos = mpos;
			merr = mlerr; mpos = mlpos; /* current best */
			goto SEARCH_SIDES;
		} else {
			if (mlerr <= lerr) {
				lerr = mlerr; lpos = mlpos;
			} else {
				break;
			}
		}
	}
	/*
	 * Search on the upper side of the middle.
	 */
	for (;;) {
		double old_mupos = mupos;
		mupos = (upos + mpos) * 0.5f;
		if (fabs(old_mupos - mupos) < EPSILON) break;
		muerr = run_one(n, mupos);
		if (merr > muerr) {
			uerr = merr; upos = mpos;
			merr = muerr; mpos = mupos; /* current best */
			goto SEARCH_SIDES;
		} else {
			if (muerr <= uerr) {
				uerr = muerr; upos = mupos;
			} else {
				break;
			}
		}
	}
//	printf("iP: %e; %e; %e\n", lpos, mpos, upos);
//	printf(" E: %e; %e; %e\n", lerr, merr, uerr);
	return stageminmaxerr_y[n];
}

static double recurse_subdivide(uint32_t m, uint32_t n);

static double recurse_one(uint32_t m, uint32_t n, double pos) {
	uint32_t j = PDIM - n + m;
	double err;
	set_candidate(j, pos);
	err = recurse_subdivide(m - 1, n);
	if (stageresult[j - 1])
		stageresult[j] = 1;
	if (err < stageminmaxerr_y[j])
		stageminmaxerr_y[j] = err;
//	printf("newer subdivide\n\tmin err so far\t%e\n",
//			stageminmaxerr_y[j]);
	return err;
}

static double recurse_subdivide(uint32_t m, uint32_t n) {
	uint32_t j = PDIM - n + m;
	begin_stage(j);
	if (m == 0) {
		run_subdivide(j);
		goto DONE;
	}
#if 0
	const uint32_t limit = loop_limits[j];
	for (uint32_t i = 0; i <= limit; ++i) {
		set_candidate_int(j, i);
		double err = recurse_subdivide(m - 1, n);
		if (stageresult[j - 1])
			stageresult[j] = 1;
		if (err < stageminmaxerr_y[j])
			stageminmaxerr_y[j] = err;
		//printf("older loop\n\tmin err so far\t%e\n",
		//		stageminmaxerr_y[j]);
	}
	//printf("\tdone\n");
#else
	/*
	 * Subdivision testing algorithm, as in innermost
	 * search except adapted for this recursive step.
	 */
	double lpos, mpos, upos;
	double lerr, merr, uerr;
	double mlpos, mupos;
	double mlerr, muerr;
	sub_bench_count = 0;
	lerr = recurse_one(m, n, (lpos = 0.0625f));
	mlerr = lerr; mlpos = lpos;
	merr = mlerr; mpos = mlpos;
	muerr = merr; mupos = mpos;
	uerr = recurse_one(m, n, (upos = 1.f));
	/*
	 * Go up by big steps to find nice intial middle.
	 */
	while (uerr < muerr) {
		double old_mupos = mupos;
		mupos = sqrt(mupos);
		if (mupos - old_mupos < EPSILON) break;
		muerr = recurse_one(m, n, mupos);
		lerr = mlerr; lpos = mlpos;
		mlerr = merr; mlpos = mpos;
		merr = muerr; mpos = mupos; /* current best */
	}
	do {
		lerr = mlerr; lpos = mlpos;
		mlerr = merr; mlpos = mpos;
		merr = muerr; mpos = mupos; /* current best */
		double old_mupos = mupos;
		mupos = sqrt(mupos);
		if (mupos - old_mupos < EPSILON) break;
		muerr = recurse_one(m, n, mupos);
	} while (muerr < merr);
	if (mlpos != mpos) {
		lpos = mlpos; lerr = mlerr;
	} else {
		mlpos = lpos; mlerr = lerr;
	}
	uerr = muerr; upos = mupos;
SEARCH_SIDES:
	/*
	 * Search on the lower side of the middle.
	 */
	for (;;) {
		double old_mlpos = mlpos;
		mlpos = (lpos + mpos) * 0.5f;
		if (fabs(old_mlpos - mlpos) < EPSILON) break;
		mlerr = recurse_one(m, n, mlpos);
		if (merr > mlerr) {
			uerr = merr; upos = mpos;
			merr = mlerr; mpos = mlpos; /* current best */
			goto SEARCH_SIDES;
		} else {
			if (mlerr <= lerr) {
				lerr = mlerr; lpos = mlpos;
			} else {
				break;
			}
		}
	}
	/*
	 * Search on the upper side of the middle.
	 */
	for (;;) {
		double old_mupos = mupos;
		mupos = (upos + mpos) * 0.5f;
		if (fabs(old_mupos - mupos) < EPSILON) break;
		muerr = recurse_one(m, n, mupos);
		if (merr > muerr) {
			uerr = merr; upos = mpos;
			merr = muerr; mpos = mupos; /* current best */
			goto SEARCH_SIDES;
		} else {
			if (muerr <= uerr) {
				uerr = muerr; upos = mupos;
			} else {
				break;
			}
		}
	}
#endif
DONE:
	return stageminmaxerr_y[j];
}

/*
 * Run n-dimensional pass, from 0 to PDIM;
 * 0-pass tests the unmodified polynomial.
 */
static int run_pass(uint32_t n) {
	int found = 0;
	for (uint32_t j = 0; j < PDIM; ++j)
		set_candidate(j, 1.f);
	if (n == 0) {
		if (try_candidate(compare_maxerr, minmaxerr_y) >= 0) {
			select_candidate();
			found = 1;
		}
	} else {
		recurse_subdivide(n - 1, n);
		found = stageresult[PDIM - n + (n - 1)];
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
