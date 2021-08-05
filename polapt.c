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

#define WRITE_PLOT_FILE 1

static inline float test_sin(float x, double scale_adj[]) {
	const float scale[] = {
		-1.f/6 * scale_adj[0],
		+1.f/120 * scale_adj[1],
		-1.f/5040 * scale_adj[2],
	};
	float x2 = x*x;
	return x + x*x2*(scale[0] + x2*(scale[1] + x2*scale[2]));
}

#define LENGTH 1000 //64 //16 //128 //1024
#define PDIM 3
#define MAXERR 1.f  // large-enough start value to accept any contender

float good_sinf[LENGTH];
float tryerr_sinf[LENGTH];
float selerr_sinf[LENGTH];
float minmaxerr_sinf = MAXERR;
float trymaxerr_sinf = MAXERR;

double scale_adj[PDIM] = {1.f, 1.f, 1.f};
double tryscale_adj[PDIM];
double selscale_adj[PDIM];
uint32_t selscale_report[PDIM][2];

/*
 * For minimizing error at end of curve only.
 */
static int compare_enderr(float minerr) {
	uint32_t i = LENGTH - 1;
	float x = (1.f - 0.5f);
	float err = test_sin(x * M_PI, tryscale_adj) - good_sinf[i];
	float abserr = fabs(err);
	minerr = fabs(minerr);
	tryerr_sinf[i] = err;
	if (abserr > trymaxerr_sinf)
		trymaxerr_sinf = abserr;
	if (abserr >= minerr)
		return 0 - (abserr > minerr);
	return 1;
}

/*
 * For minimizing error all over curve.
 *
 * Uses compare_enderr() as tie-breaker
 * when the threshold is reached but not exceeded.
 */
static int compare_maxerr(float minerr) {
	minerr = fabs(minerr);
	for (uint32_t i = 0, end = LENGTH - 1; i <= end; ++i) {
		float x = (i * 1.f/end - 0.5f);
		float err = test_sin(x * M_PI, tryscale_adj) - good_sinf[i];
		float abserr = fabs(err);
		tryerr_sinf[i] = err;
		if (abserr > trymaxerr_sinf)
			trymaxerr_sinf = abserr;
		if (abserr >= minerr) {
			if (abserr == minerr) {
				int cmp = compare_enderr(selerr_sinf[end]);
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
static int compare_endfirstmaxerr(float minerr) {
	if (compare_enderr(selerr_sinf[LENGTH - 1]) < 0)
		return -1;
	return compare_maxerr(minerr);
}

static int try_candidate(const uint32_t pcoeffs[PDIM], float err_threshold) {
	trymaxerr_sinf = 0.f;
	for (uint32_t j = 0; j < PDIM; ++j) {
		uint32_t q = pcoeffs[j];
		tryscale_adj[j] = scale_adj[j];
		if (q > 0)
			tryscale_adj[j] *= ((double)(q - 1) / (double) q);
	}
//	return compare_endfirstmaxerr(err_threshold);
	return compare_maxerr(err_threshold);
}

static void select_candidate(const uint32_t pcoeffs[PDIM]) {
	for (uint32_t i = 0, end = LENGTH; i < end; ++i) {
		selerr_sinf[i] = tryerr_sinf[i];
	}
	if (trymaxerr_sinf < minmaxerr_sinf)  {
		minmaxerr_sinf = trymaxerr_sinf;
	}
	for (uint32_t j = 0; j < PDIM; ++j) {
		uint32_t q = pcoeffs[j];
		if (q == 0) {
			selscale_report[j][0] = 0;
			selscale_report[j][1] = 0;
		} else {
			selscale_report[j][0] = q - 1;
			selscale_report[j][1] = q;
		}
		selscale_adj[j] = tryscale_adj[j];
	}
}

static void print_report(void) {
	printf("Max.err. %.11e\tEnd.err. %.11e\n",
			minmaxerr_sinf, selerr_sinf[LENGTH - 1]);
	for (uint32_t j = 0; j < PDIM; ++j) {
		char label = 'A' + j;
		printf("%c==%.11e\t(%.11e * (%d.0/%d.0))\n",
				label, selscale_adj[j], scale_adj[j],
				selscale_report[j][0], selscale_report[j][1]);
	}
}

/*
 * Makes result of previous pass apply to next pass.
 */
static void apply_selected(void) {
	printf("(Applying selected coefficents.)\n");
	for (uint32_t j = 0; j < PDIM; ++j) {
		scale_adj[j] = selscale_adj[j];
	}
}

static const uint32_t loop_limits[PDIM] = {
	100000, //100000
	1000,   //10000
	100,    //1000
};

static int test_linear(uint32_t pcoeffs[PDIM], uint32_t n,
		uint32_t from, uint32_t to) {
	uint32_t hits = 0;
	for (uint32_t i = from; i <= to; ++i) {
		pcoeffs[n] = i;
		if (try_candidate(pcoeffs, minmaxerr_sinf) > 0) {
			select_candidate(pcoeffs);
			++hits;
		}
	}
	return hits;
}

static int recurse_linear(uint32_t pcoeffs[PDIM], uint32_t n) {
	const uint32_t limit = loop_limits[n];
	int found = 0;
	if (n == PDIM - 1) {
		found = (test_linear(pcoeffs, n, 0, limit) > 0);
	} else {
		for (uint32_t i = 0; i <= limit; ++i) {
			pcoeffs[n] = i;
			if (recurse_linear(pcoeffs, n + 1)) found = 1;
		}
	}
	return found;
}

/*
 * Run n-dimensional pass, from 0 to PDIM;
 * 0-pass tests the unmodified polynomial.
 */
static int run_pass(uint32_t n) {
	uint32_t pcoeffs[PDIM] = {0, 0, 0};
	int found = 0;
	if (n == 0) {
		if (try_candidate(pcoeffs, minmaxerr_sinf) > 0) {
			select_candidate(pcoeffs);
			found = 1;
		}
	} else {
		found = (recurse_linear(pcoeffs, PDIM - n) > 0);
	}
	if (found)
		print_report();
	return found;
}

int main(void) {
#if WRITE_PLOT_FILE
	FILE *f = fopen("plot.txt", "w");
#endif
	for (uint32_t i = 0, end = LENGTH - 1; i <= end; ++i) {
		float x = (i * 1.f/end - 0.5f);
		good_sinf[i] = sinf(x * M_PI);
		selerr_sinf[i] = MAXERR;
	}
	run_pass(0); /* also print stats for unmodified polynomial */
	for (uint32_t n = 1; n <= PDIM; ++n) {
		run_pass(n);
		if (n < PDIM)
			apply_selected();
	}
#if WRITE_PLOT_FILE
	for (uint32_t i = 0, end = LENGTH - 1; i <= end; ++i) {
		float x = (i * 1.f/end - 0.5f);
		fprintf(f, "%.11f\t%.11f\n", x, selerr_sinf[i]);
	}
	fclose(f);
#endif
	return 0;
}
