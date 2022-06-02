/* Polynomial approximation optimizer: Extra iterative test functions.
 * Copyright (c) 2021-2022 Joel K. Pettersson
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

/* Get array of per-dimension TEST_T offset by number of dimensions skipped. */
#define ARR_DIMOFFSET(arr, used) \
	((arr) + sizeof(arr)/sizeof(TEST_T) - (used))

/* -0.0, *1.0 */
static inline TEST_T test_sqrt_r1_d4(TEST_T x, long double scale_adj[]) {
	TEST_T scale[] = {
		+24.344885f/6,
		-58.344885f/6,
		+67.344885f/6,
		-27.344885f/6,
	};
	for (size_t i = 0; i < PDIM; ++i)
		ARR_DIMOFFSET(scale, PDIM)[i] *= scale_adj[i];
	TEST_T x2 = x*x;
	TEST_T xa = fabsl(x);
	return x*(scale[0] + xa*(scale[1] + xa*(scale[2] + xa*scale[3])));
}

/* 0.0, *1.0 */
static inline TEST_T test_sinramp_jkp(TEST_T x, long double scale_adj[]) {
	TEST_T scale[] = {
		+118.f/75,
		-196.f/75,
		+88.0f/75,
	};
	for (size_t i = 0; i < PDIM; ++i)
		ARR_DIMOFFSET(scale, PDIM)[i] *= scale_adj[i];
	x -= 0.5f;
	TEST_T x2 = x*x;
	return 0.5f + x*(scale[0] + x2*(scale[1] + x2*scale[2]));
}

/* -0.5, *2.0 */
static inline TEST_T test_fabs_d16(TEST_T x, long double scale_adj[]) {
	TEST_T scale[] = {
		+1.f,
		-1.f,
		+1.f,
		-1.f,
	};
	for (size_t i = 0; i < PDIM; ++i)
		ARR_DIMOFFSET(scale, PDIM)[i] *= scale_adj[i];
	TEST_T xp = (x+x - x*fabsl(x));
	xp *= xp;
	return xp*(scale[0] + xp*(scale[1] + xp*(scale[2] + xp*scale[3])));
}

/* -0.5, *2.0 */
static inline TEST_T test_fabs_d12(TEST_T x, long double scale_adj[]) {
	TEST_T scale[] = {
		+1.f,
		-1.f,
		+1.f,
		-1.f,
	};
	for (size_t i = 0; i < PDIM; ++i)
		ARR_DIMOFFSET(scale, PDIM)[i] *= scale_adj[i];
	TEST_T xp = (3 - 2*fabsl(x))*x*x;
	return xp*(scale[0] + xp*(scale[1] + xp*(scale[2] + xp*scale[3])));
}
