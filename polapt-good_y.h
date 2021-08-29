/* Polynomial approximation optimizer: Extra reference functions.
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

/*
 * Some extra polynomials that can be interesting
 * to use as "good" functions for
 * comparing with the test functions.
 *
 * (This is in addition to sin() and other near-"ideal" functions.)
 */

/* - sin(x) approximations - */

/* max error: 1.388788e-05
   (tweaked manually using graph plot program to be close near +/- PI/2) */
static inline float SGS_sinf_t7_old(float x) {
	const float scale7 = -1.f/5040 * 29.f/30;
	float x2 = x*x;
	return x + x*x2*(-1.f/6 + x2*(1.f/120 + x2*scale7));
}

/* max error: 3.576279e-07
   (tweaked manually using graph plot program to be close near +/- PI/2) */
static inline float SGS_sinf_t9_old(float x) {
	const float scale9 = 1.f/362880 * 44.f/45;
	float x2 = x*x;
	return x + x*x2*(-1.f/6 + x2*(1.f/120 + x2*(-1.f/5040 + x2*scale9)));
}

/*
 * A higher-order Chebyshev approximation by Colin Wallace that's
 * much faster than glibc sin(x) but a little heavy-duty for some purposes.
 *
 * https://web.archive.org/web/20200628195036/http://mooooo.ooo/chebyshev-sine-approximation/
 *
 * max error: 1.788139e-07
 */
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

/* - other functions - */

/*
 * Simplistic sine approximation I made in 2010,
 * with added scale factor to bring output close
 * to +/- 1.0 near +/- PI/2. A distorting touch.
 *
 * max error: 5.267441e-03
 */
static inline float sin_simple2010(float x) {
	const float a = 0.99535533335807127063;
	float v = x * 0.5f;
	float v2 = v * v, v3 = v2 * v, va = fabs(v);
	return (x - v3 - v3*va + v3*v2) * a;
}

/*
 * Simplistic approximation of square root of sine
 * (mirrored for the negative half), scaled so end (+/- PI/2)
 * of input range has near exactly +/- 1.0 output.
 *
 * max error: 9.937980e-02
 */
static inline float simple_srsf(float x) {
	const float scale[] = {
		+21.0/6  * 1.00855031486950772112,
		-44.f/7  * 1.00855031486950772112,
		+74.f/15 * 1.00855031486950772112,
		-4.f/3   * 1.00855031486950772112,
	};
	float xa = fabs(x);
	return x*(scale[0] + xa*(scale[1] + xa*(scale[2] + xa*scale[3])));
}
