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
#define PI 3.14159265358979323846
#define END_SHIFT(x) 1.f/(1.f + (x))

#define lin(x) x
#define linf(x) ((float) x)
#define cosramp(x) \
	(1.f - (cos((x) * PI) + 1.f)*0.5f)
#define cosrampf(x) \
	(1.f - (cosf((x) * PI) + 1.f)*0.5f)

static inline double sqrtp1(double x) {
	return sqrt(x + 1.0);
}

static inline float sqrtp1f(float x) {
	return sqrtf(x + 1.0);
}

static inline double srs(double x) {
	x = sin(x);
	return x >= 0.f ? sqrt(x) : -sqrt(-x);
}

static inline float srsf(float x) {
	x = sinf(x);
	return x >= 0.f ? sqrtf(x) : -sqrtf(-x);
}

#include "polapt-good_y.h"

/*
 * The program.
 */

/* What to run? */
#define TEST_X(x) ((x - 0.5f) * PI)
#define TEST_Y test_sin_t5
#define GOOD_Y sin
#define TEST_T long double
#define TEST_C compare_maxerr
#define ENDS_Y END_SHIFT(0.0f) // compensation for extra offset goes here

#define TAB_LEN 1000 //64 //16 //128 //1024
#define MAX_ERR 1.f  // large-enough start value to accept any contender
#define ERR_BIAS 1.f // value between 0 and 1 to give weighed preference
#define EPSILON 1.e-19
#define PDIM 3       // this also limits which test functions are usable

/* Test more than starting point? */
#define RUN_TESTS       1

/* Adjust result to primarily minimize error at endpoints? */
#define FIT_ENDPOINTS   1

/* Produce file suitable for gnuplot? */
#define WRITE_PLOT_FILE 1

#include "polapt-test_y.h"

/* -0.5, *PI */
static inline TEST_T test_sin_t5(TEST_T x, long double scale_adj[]) {
	const TEST_T scale[] = {
		+1.f     * scale_adj[0],
		-1.f/6   * scale_adj[1],
		+1.f/120 * scale_adj[2],
	};
	TEST_T x2 = x*x;
	return x*(scale[0] + x2*(scale[1] + x2*scale[2]));
}

/* -0.5, *PI */
static inline TEST_T test_sin_t7(TEST_T x, long double scale_adj[]) {
	const TEST_T scale[] = {
		+1.f      * scale_adj[0],
		-1.f/6    * scale_adj[1],
		+1.f/120  * scale_adj[2],
		-1.f/5040 * scale_adj[3],
	};
	TEST_T x2 = x*x;
	return x*(scale[0] + x2*(scale[1] + x2*(scale[2] + x2*scale[3])));
}

/* -0.5, *2.0 */
static inline TEST_T test_sqrtp1_t4(TEST_T x, long double scale_adj[]) {
	const TEST_T scale[] = {
		+1.f/2   * scale_adj[0],
		-1.f/8   * scale_adj[1],
		+1.f/16  * scale_adj[2],
		-5.f/128 * scale_adj[3],
	};
	return 1. + x*(scale[0] + x*(scale[1] + x*(scale[2] + x*scale[3])));
}

TEST_T good_y[TAB_LEN];
TEST_T tryerr_y[TAB_LEN];
TEST_T selerr_y[TAB_LEN];
TEST_T minmaxerr_y = MAX_ERR;
TEST_T trymaxerr_y = MAX_ERR;
TEST_T stageminmaxerr_y[PDIM];
int stageresult[PDIM];

long double scale_adj[PDIM];
long double tryscale_adj[PDIM];
long double selscale_adj[PDIM];
long double trypos[PDIM];
long double selpos[PDIM];

static inline int compare_point(TEST_T minerr, uint32_t i) {
	TEST_T x = i * 1.f/(TAB_LEN - 1);
	TEST_T err = TEST_Y(TEST_X(x), tryscale_adj) - good_y[i];
	TEST_T abserr = fabsl(err);
	tryerr_y[i] = err;
	if (abserr > trymaxerr_y)
		trymaxerr_y = abserr;
	if (abserr >= minerr)
		return 0 - (abserr > minerr);
	return 1;
}

/*
 * For minimizing error at the ends of the curve only.
 */
static int compare_enderr(TEST_T minerr) {
	int cmp0 = compare_point(minerr, 0);
	int cmp1 = compare_point(minerr, TAB_LEN - 1);
	int ret = 1;
	if (cmp0 <= 0 || cmp1 <= 0) {
		ret = (cmp0 < 0 || cmp1 < 0) ? -1 : 0;
	}
	return ret;
}

/*
 * For minimizing error all over curve.
 */
static int compare_maxerr(TEST_T minerr) {
	int ret = 1;
	for (uint32_t i = 0, end = TAB_LEN - 1; i <= end; ++i) {
		int cmp = compare_point(minerr, i);
		if (cmp <= 0) {
			if (cmp < 0)
				return -1;
			ret = 0;
		}
	}
	return ret;
}

/*
 * For minimizing error at the end of the curve primarily,
 * over the rest of the curve secondarily.
 */
static int compare_enderr_maxerr(TEST_T minerr) {
	int end_cmp = compare_enderr(minerr * ERR_BIAS);
	if (end_cmp < 0)
		return -1;
	return compare_maxerr(minerr);
}

static inline void begin_stage(uint32_t n) {
	stageminmaxerr_y[n] = MAX_ERR;
	stageresult[n] = 0;
}

static inline void set_candidate(uint32_t n, long double pos) {
	trypos[n] = pos;
}

static int try_candidate(int (*compare)(TEST_T minerr), TEST_T minerr) {
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
	printf("Max.err. %.11Le\tEnd.err. %.11Le\n",
			minmaxerr_y, selerr_y[TAB_LEN - 1]);
	for (uint32_t j = 0; j < PDIM; ++j) {
		char label = 'A' + j;
		printf("%c==%.20Lf\t(%.14Lf * %.14Lf)\n",
				label, selscale_adj[j],
				scale_adj[j], selpos[j]);
	}
}

/*
 * Makes result of previous pass apply to next pass.
 */
static void apply_selected(int msg) {
	if (msg)
		printf("(Applying selected coefficients.)\n");
	for (uint32_t j = 0; j < PDIM; ++j) {
		scale_adj[j] = selscale_adj[j];
	}
}

static inline long double moved_pos(long double from, long double to) {
	return (from + to) * 0.5;
}

uint32_t bench_count, sub_bench_count;
static TEST_T run_one(uint32_t n, long double pos) {
	++bench_count;
	++sub_bench_count;
	set_candidate(n, pos);
	try_candidate(TEST_C, stageminmaxerr_y[n]);
	if (trymaxerr_y < stageminmaxerr_y[n])
		stageminmaxerr_y[n] = trymaxerr_y;
	return trymaxerr_y;
}

/*
 * Subdivision testing algorithm. Looks for number between 0.0 and 1.0,
 * halving the size of steps to take within the range until a number is
 * chosen. May make roughly a hundred tests, before picking the number.
 */
static TEST_T run_subdivide(uint32_t n) {
	long double lpos, mpos, upos;
	long double lerr, merr, uerr;
	long double mlpos, mupos;
	long double mlerr, muerr;
	sub_bench_count = 0;
	lpos = .5f;
	mpos = 1.f;
	upos = 2.f;
	lerr = run_one(n, lpos);
	merr = run_one(n, mpos);
	uerr = run_one(n, upos);
	/* Find lower bound by going down through squaring. */
	while (lerr < merr) {
		if (merr < uerr) {
			upos = mpos;
			uerr = merr;
		}
		mpos = lpos;
		merr = lerr;
		lpos *= lpos;
		lerr = run_one(n, lpos);
	}
	/* Find upper bound by going up through squaring. */
	while (uerr < merr) {
		if (merr < lerr) {
			lpos = mpos;
			lerr = merr;
		}
		mpos = upos;
		merr = uerr;
		upos *= upos;
		uerr = run_one(n, upos);
	}
	/* Move towards bottom of error curve valley or slope. */
	for (;;) {
		mlpos = moved_pos(lpos, mpos);
		mlerr = run_one(n, mlpos);
		mupos = moved_pos(upos, mpos);
		muerr = run_one(n, mupos);
		if (mlerr < muerr) {
			/* ? ? + */
			if (mlerr < merr) {
				/* - 0 +  (Go down...) */
				upos = mpos;
				uerr = merr;
				mpos = mlpos;
				merr = mlerr;
			} else {
				/* 0 - +  (Go in...) */
				if (muerr < uerr) {
					upos = mupos;
					uerr = muerr;
				}
				lpos = mlpos;
				lerr = mlerr;
			}
			if (mpos <= lpos + EPSILON) break;
		} else {
			/* + ? ? */
			if (merr <= muerr) {
				/* + - 0  (Go in...) */
				if (mlerr <= lerr) {
					lpos = mlpos;
					lerr = mlerr;
				}
				upos = mupos;
				uerr = muerr;
			} else {
				/* + 0 -  (Go up...) */
				lpos = mpos;
				lerr = merr;
				mpos = mupos;
				merr = muerr;
			}
			if (mpos >= upos - EPSILON) break;
		}
		if ((mpos >= upos - EPSILON) && (mpos <= lpos + EPSILON))
			break;
	}
	set_candidate(n, mpos);
	if (try_candidate(TEST_C, minmaxerr_y) >= 0) {
		select_candidate();
		stageresult[n] = 1;
	}
	//printf("(SUB BENCH %u)\n", sub_bench_count);
	return stageminmaxerr_y[n];
}

static TEST_T recurse_subdivide(uint32_t m, uint32_t n);

static TEST_T recurse_one(uint32_t m, uint32_t n, long double pos) {
	uint32_t j = PDIM - n + m;
	TEST_T err;
	set_candidate(j, pos);
	err = recurse_subdivide(m - 1, n);
	if (stageresult[j - 1])
		stageresult[j] = 1;
	if (err < stageminmaxerr_y[j])
		stageminmaxerr_y[j] = err;
	return err;
}

static TEST_T recurse_subdivide(uint32_t m, uint32_t n) {
	uint32_t j = PDIM - n + m;
	begin_stage(j);
	if (m == 0) {
		run_subdivide(j);
		goto DONE;
	}
	/*
	 * Subdivision testing algorithm, as in innermost
	 * search except adapted for this recursive step.
	 */
	long double lpos, mpos, upos;
	long double lerr, merr, uerr;
	long double mlpos, mupos;
	long double mlerr, muerr;
	lpos = .5f;
	mpos = 1.f;
	upos = 2.f;
	lerr = recurse_one(m, n, lpos);
	merr = recurse_one(m, n, mpos);
	uerr = recurse_one(m, n, upos);
	/* Find lower bound by going down through squaring. */
	while (lerr < merr) {
		if (merr < uerr) {
			upos = mpos;
			uerr = merr;
		}
		mpos = lpos;
		merr = lerr;
		lpos *= lpos;
		lerr = recurse_one(m, n, lpos);
	}
	/* Find upper bound by going up through squaring. */
	while (uerr < merr) {
		if (merr < lerr) {
			lpos = mpos;
			lerr = merr;
		}
		mpos = upos;
		merr = uerr;
		upos *= upos;
		uerr = recurse_one(m, n, upos);
	}
	/* Move towards bottom of error curve valley or slope. */
	for (;;) {
		mlpos = moved_pos(lpos, mpos);
		mlerr = recurse_one(m, n, mlpos);
		mupos = moved_pos(upos, mpos);
		muerr = recurse_one(m, n, mupos);
		if (mlerr < muerr) {
			/* ? ? + */
			if (mlerr < merr) {
				/* - 0 +  (Go down...) */
				upos = mpos;
				uerr = merr;
				mpos = mlpos;
				merr = mlerr;
			} else {
				/* 0 - +  (Go in...) */
				if (muerr < uerr) {
					upos = mupos;
					uerr = muerr;
				}
				lpos = mlpos;
				lerr = mlerr;
			}
			if (mpos <= lpos + EPSILON) break;
		} else {
			/* + ? ? */
			if (merr <= muerr) {
				/* + - 0  (Go in...) */
				if (mlerr <= lerr) {
					lpos = mlpos;
					lerr = mlerr;
				}
				upos = mupos;
				uerr = muerr;
			} else {
				/* + 0 -  (Go up...) */
				lpos = mpos;
				lerr = merr;
				mpos = mupos;
				merr = muerr;
			}
			if (mpos >= upos - EPSILON) break;
		}
		if ((mpos >= upos - EPSILON) && (mpos <= lpos + EPSILON))
			break;
	}
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
		if (try_candidate(TEST_C, MAX_ERR) >= 0) {
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

/*
 * Uses end point error measure for current pick
 * to minimize end error at the cost of overall accuracy.
 */
static void generate_endfitted(void) {
	//if (selerr_y[0] > EPSILON || selerr_y[TAB_LEN - 1] > EPSILON)
	// FIXME: deal with both end points...
	TEST_T end1 = selerr_y[TAB_LEN - 1], abs_end1 = fabsl(end1);
	if (abs_end1 < EPSILON || abs_end1 >= 1.f)
		return;
	printf("(Creating variation for fitted endpoints.)\n");
	apply_selected(0);
	long double scale1 = END_SHIFT(end1 * ENDS_Y);
	for (uint32_t j = 0; j < PDIM; ++j) {
		scale_adj[j] *= scale1;
	}
	minmaxerr_y = MAX_ERR;
	run_pass(0);
}

int main(void) {
#if WRITE_PLOT_FILE
	FILE *f = fopen("plot.txt", "w");
#endif
	for (uint32_t j = 0; j < PDIM; ++j)
		scale_adj[j] = 1.f;
	for (uint32_t i = 0, end = TAB_LEN - 1; i <= end; ++i) {
		TEST_T x = i * 1.f/end;
		good_y[i] = GOOD_Y(TEST_X(x));
		selerr_y[i] = MAX_ERR;
	}
	run_pass(0); /* also print stats for unmodified polynomial */
#if RUN_TESTS
	for (uint32_t n = 1; n <= PDIM; ++n) {
		run_pass(n);
		if (n < PDIM)
			apply_selected(1);
	}
# if FIT_ENDPOINTS
	generate_endfitted();
# endif
#endif
#if WRITE_PLOT_FILE
	for (uint32_t i = 0, end = TAB_LEN - 1; i <= end; ++i) {
		TEST_T x = i * 1.f/end;
		fprintf(f, "%.11Lf\t%.11Lf\n", TEST_X(x), selerr_y[i]);
	}
	fclose(f);
#endif
	printf("BENCH %u\n", bench_count);
	return 0;
}
