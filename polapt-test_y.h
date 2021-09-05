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
static inline TEST_T test_sqrt_r1_d4(TEST_T x, double scale_adj[]) {
	const TEST_T scale[] = {
		+24.344885f/6 * scale_adj[0],
		-58.344885f/6 * scale_adj[1],
		+67.344885f/6 * scale_adj[2],
		-27.344885f/6 * scale_adj[3],
	};
	TEST_T x2 = x*x;
	TEST_T xa = fabs(x);
	return x*(scale[0] + xa*(scale[1] + xa*(scale[2] + xa*scale[3])));
}

/* -0.5, *2.0 */
static inline TEST_T test_fabs_d16(TEST_T x, double scale_adj[]) {
	const TEST_T scale[] = {
		+1.f * scale_adj[0],
		-1.f * scale_adj[1],
		+1.f * scale_adj[2],
		-1.f * scale_adj[3],
	};
	TEST_T xp = (x+x - x*fabs(x));
	xp *= xp;
	return xp*(scale[0] + xp*(scale[1] + xp*(scale[2] + xp*scale[3])));
}
