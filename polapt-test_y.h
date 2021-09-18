/* Polynomial approximation optimizer: Extra iterative test functions.
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

/* -0.0, *1.0 */
static inline TEST_T test_sqrt_r1_d4(TEST_T x, long double scale_adj[]) {
	const TEST_T scale[] = {
		+24.344885f/6 * scale_adj[0],
		-58.344885f/6 * scale_adj[1],
		+67.344885f/6 * scale_adj[2],
		-27.344885f/6 * scale_adj[3],
	};
	TEST_T x2 = x*x;
	TEST_T xa = fabsl(x);
	return x*(scale[0] + xa*(scale[1] + xa*(scale[2] + xa*scale[3])));
}

/* 0.0, *1.0 */
static inline TEST_T test_cosramp_jkp(TEST_T x, long double scale_adj[]) {
	const TEST_T scale[] = {
		+118.f/75 * scale_adj[0],
		-196.f/75 * scale_adj[1],
		+88.0f/75 * scale_adj[2],
	};
	x -= 0.5f;
	TEST_T x2 = x*x;
	return 0.5f + x*(scale[0] + x2*(scale[1] + x2*scale[2]));
#if 0 /* version using a cos-style polynomial, worse for this use */
	const TEST_T scale[] = {
		+2.f    * scale_adj[0],
		-15.f/8 * scale_adj[1],
		+1.f/2  * scale_adj[2],
	};
	TEST_T x2 = x*x;
	return x2*(scale[0] + x2*(scale[1] + x2*scale[2]));
#endif
}

/* -0.5, *2.0 */
static inline TEST_T test_fabs_d16(TEST_T x, long double scale_adj[]) {
	const TEST_T scale[] = {
		+1.f * scale_adj[0],
		-1.f * scale_adj[1],
		+1.f * scale_adj[2],
		-1.f * scale_adj[3],
	};
	TEST_T xp = (x+x - x*fabsl(x));
	xp *= xp;
	return xp*(scale[0] + xp*(scale[1] + xp*(scale[2] + xp*scale[3])));
}

/* -0.5, *2.0 */
static inline TEST_T test_fabs_d12(TEST_T x, long double scale_adj[]) {
	const TEST_T scale[] = {
		+1.f * scale_adj[0],
		-1.f * scale_adj[1],
		+1.f * scale_adj[2],
		-1.f * scale_adj[3],
	};
	TEST_T xp = (3 - 2*fabsl(x))*x*x;
	return xp*(scale[0] + xp*(scale[1] + xp*(scale[2] + xp*scale[3])));
}
