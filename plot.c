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

#include <stdio.h>
#include <math.h>

static inline float SGS_sinf_t7(float x) {
	const float scale7 = -1.f/5040 * 29.f/30;
	float x2 = x*x;
	return x + x*x2*(-1.f/6 + x2*(1.f/120 + x2*scale7));
}

static inline float SGS_sinf_t9(float x) {
	const float scale9 = 1.f/362880 * 44.f/45;
	float x2 = x*x;
	return x + x*x2*(-1.f/6 + x2*(1.f/120 + x2*(-1.f/5040 + x2*scale9)));
}

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

/* two varied */
//	const float scale5 = 1.f/120 * (869.f/870);
//	const float scale7 = -1.f/5040 * 18.f/19;
/* three varied */
//	const float scale3 = -1.f/6 * (15000.f/15001);
//	const float scale5 = 1.f/120 * (359.f/360);
//	const float scale7 = -1.f/5040 * 13.f/14;
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
float good_sinf[LENGTH];
float selerr_sinf[LENGTH];
float curerr_sinf[LENGTH];
float maxerr_sinf = 1.f; // large-enough start value to accept any contender

#define A_TRY 100000 //10000 //100000
#define B_TRY 1000  //10000
#define C_TRY 100   //1000

int main() {
	FILE *f = fopen("plot.txt", "w");
	double scale_adj[3] = {1.f, 1.f, 29.f/30};
	int scale_report[3][2];
	double try_scale_adj[3];
	int new_report = 0, want_new_report = 1;
	for (int i = 0, end = LENGTH; i < end; ++i) {
		float x = (i * 1.f/(end - 1) - 0.5f);
		good_sinf[i] = sinf(x * M_PI);
	}
	for (int A = 0; A <= A_TRY + 1; ++A) {
		if (A == 0) {
			try_scale_adj[0] = 1.f;
		} else {
			try_scale_adj[0] = ((double)(A - 1) / (double) A);
		}
		for (int B = 0; B <= B_TRY + 1; ++B) {
			if (B == 0) {
				try_scale_adj[1] = 1.f;
			} else {
				try_scale_adj[1] =
					((double)(B - 1) / (double) B);
			}
			for (int C = 0; C <= C_TRY + 1; ++C) {
				if (C == 0) {
					try_scale_adj[2] = 1.f;
				} else {
					if (C == C_TRY && (!B || B == B_TRY))
						want_new_report = 1;
					try_scale_adj[2] =
						((double)(C - 1) / (double) C);
				}
				for (int i = 0, end = LENGTH; i < end; ++i) {
					float x = (i * 1.f/(end - 1) - 0.5f);
					float err = good_sinf[i] -
						test_sin(x * M_PI,
								try_scale_adj);
					if (fabs(err) > maxerr_sinf)
						goto SKIP;
					curerr_sinf[i] = err;
				}
				/*
				 * New selection made...
				 */
				float selmaxerr = 0.f;
				for (int i = 0, end = LENGTH; i < end; ++i) {
					if (selmaxerr < curerr_sinf[i])
						selmaxerr = curerr_sinf[i];
					selerr_sinf[i] = curerr_sinf[i];
				}
				maxerr_sinf = selmaxerr;
				new_report = 1;
				scale_adj[0] = try_scale_adj[0];
				scale_adj[1] = try_scale_adj[1];
				scale_adj[2] = try_scale_adj[2];
				if (A == 0) {
					scale_report[0][0] = 0;
					scale_report[0][1] = 0;
				} else {
					scale_report[0][0] = A - 1;
					scale_report[0][1] = A;
				}
				if (B == 0) {
					scale_report[1][0] = 0;
					scale_report[1][1] = 0;
				} else {
					scale_report[1][0] = B - 1;
					scale_report[1][1] = B;
				}
				if (C == 0) {
					scale_report[2][0] = 0;
					scale_report[2][1] = 0;
				} else {
					scale_report[2][0] = C - 1;
					scale_report[2][1] = C;
				}
			SKIP:
				if (new_report && want_new_report) {
					new_report = 0;
					want_new_report = 0;
					printf(
"Max.err. %.11f\nA==%.11f\t(%d, %d)\nB==%.11f\t(%d, %d)\nC==%.11f\t(%d, %d)\n",
maxerr_sinf,
scale_adj[0], scale_report[0][0], scale_report[0][1],
scale_adj[1], scale_report[1][0], scale_report[1][1],
scale_adj[2], scale_report[2][0], scale_report[2][1]
						);
				}
				continue;
			}
		}
	}
	for (int i = 0, end = LENGTH; i < end; ++i) {
		float x = (i * 1.f/(end - 1) - 0.5f);
		fprintf(f, "%.11f\t%.11f\n", x, selerr_sinf[i]);
	}
	fclose(f);
//	printf("A==%.11f\t(%d, %d)\nB==%.11f\t(%d, %d)\nC==%.11f\t(%d, %d)\n",
//			scale_adj[0], scale_report[0][0], scale_report[0][1],
//			scale_adj[1], scale_report[1][0], scale_report[1][1],
//			scale_adj[2], scale_report[2][0], scale_report[2][1]
//			);
	return 0;
}
