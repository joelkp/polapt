/* Plotting program for iterating improved Taylor polynomials
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

/* max error: 1.388788e-05 */
static inline float SGS_sinf_t7_old(float x) {
	const float scale7 = -1.f/5040 * 29.f/30;
	float x2 = x*x;
	return x + x*x2*(-1.f/6 + x2*(1.f/120 + x2*scale7));
}

/* max error: 3.576279e-07 */
static inline float SGS_sinf_t9_old(float x) {
	const float scale9 = 1.f/362880 * 44.f/45;
	float x2 = x*x;
	return x + x*x2*(-1.f/6 + x2*(1.f/120 + x2*(-1.f/5040 + x2*scale9)));
}

/* max error: 1.788139e-07 */
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

static inline float test_sin(float x, double scale_adj[]) {
	const float scale[] = {
		-1.f/6 * scale_adj[0],
		+1.f/120 * scale_adj[1],
		-1.f/5040 * scale_adj[2],
	};
	float x2 = x*x;
	return x + x*x2*(scale[0] + x2*(scale[1] + x2*scale[2]));
}

#define LENGTH 1000
#define PDIM 3
#define MAXERR 1.f  // large-enough start value to accept any contender

float good_sinf[LENGTH];
float tryerr_sinf[LENGTH];
float selerr_sinf[LENGTH];
float maxerr_sinf = MAXERR;

double scale_adj[PDIM] = {1.f, 1.f, 1.f};
double tryscale_adj[PDIM];
uint32_t scale_report[PDIM][2];

static float try_candidate(const uint32_t pcoeffs[PDIM], float err_threshold) {
	float try_maxerr = 0.f;
	for (uint32_t j = 0; j < PDIM; ++j) {
		uint32_t q = pcoeffs[j];
		if (q == 0) {
			tryscale_adj[j] = 1.f;
		} else {
			tryscale_adj[j] = ((double)(q - 1) / (double) q);
		}
	}
	for (uint32_t i = 0, end = LENGTH; i < end; ++i) {
		float x = (i * 1.f/(end - 1) - 0.5f);
		float err = good_sinf[i] - test_sin(x * M_PI, tryscale_adj);
		float abserr = fabs(err);
		if (abserr > err_threshold)
			return abserr;
		if (abserr > try_maxerr)
			try_maxerr = abserr;
		tryerr_sinf[i] = err;
	}
	return try_maxerr;
}

static void select_candidate(const uint32_t pcoeffs[PDIM], float try_maxerr) {
	/*
	 * New selection made...
	 */
	for (uint32_t i = 0, end = LENGTH; i < end; ++i) {
		selerr_sinf[i] = tryerr_sinf[i];
	}
	maxerr_sinf = try_maxerr;
	for (uint32_t j = 0; j < PDIM; ++j) {
		uint32_t q = pcoeffs[j];
		if (q == 0) {
			scale_report[j][0] = 0;
			scale_report[j][1] = 0;
		} else {
			scale_report[j][0] = q - 1;
			scale_report[j][1] = q;
		}
		scale_adj[j] = tryscale_adj[j];
	}
}

static void print_report(void) {
	printf("Max.err. %e\n", maxerr_sinf);
	for (uint32_t j = 0; j < PDIM; ++j) {
		char label = 'A' + j;
		printf("%c==%.11f\t(%d, %d)\n", label, scale_adj[j],
				scale_report[j][0], scale_report[j][1]);
	}
}

#define A_TRY 100000 //100000
#define B_TRY 1000   //10000
#define C_TRY 100    //1000

static int test_C(uint32_t pcoeffs[PDIM], uint32_t n) {
	uint32_t hits = 0;
	float tryerr, trymaxerr = MAXERR;
	for (uint32_t i = 0; i <= C_TRY; ++i) {
		pcoeffs[n] = i;
		tryerr = try_candidate(pcoeffs, trymaxerr);
		if (tryerr < trymaxerr) {
			trymaxerr = tryerr;
			if (tryerr < maxerr_sinf) {
				select_candidate(pcoeffs, tryerr);
				++hits;
			}
		}
	}
	return hits;
}

int main(void) {
	FILE *f = fopen("plot.txt", "w");
	int new_report = 0;
	for (uint32_t i = 0, end = LENGTH; i < end; ++i) {
		float x = (i * 1.f/(end - 1) - 0.5f);
		good_sinf[i] = sinf(x * M_PI);
	}
	uint32_t pcoeffs[PDIM] = {0, 0, 0};
	for (uint32_t A = 0; A <= A_TRY + 1; ++A) {
		pcoeffs[0] = A;
		for (uint32_t B = 0; B <= B_TRY + 1; ++B) {
			pcoeffs[1] = B;
			if (test_C(pcoeffs, 2) > 0)
				new_report = 1;
			if (B == 0 && new_report) {
				new_report = 0;
				print_report();
			}
		}
		if (A == 0 && new_report) {
			new_report = 0;
			print_report();
		}
	}
	if (new_report) {
		new_report = 0;
		print_report();
	}
	for (uint32_t i = 0, end = LENGTH; i < end; ++i) {
		float x = (i * 1.f/(end - 1) - 0.5f);
		fprintf(f, "%.11f\t%.11f\n", x, selerr_sinf[i]);
	}
	fclose(f);
	return 0;
}
