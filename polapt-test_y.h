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

/* 0.0, *1.0 */
static inline TEST_T test_cosline_jkp(TEST_T x, double scale_adj[]) {
	const TEST_T scale[] = {
		+0.49998307401098157632 * scale_adj[PDIM - 4],
		+1.57026686101023282838 * scale_adj[PDIM - 3],
		-2.56836572869766729953 * scale_adj[PDIM - 2],
		+1.14973477027553139252 * scale_adj[PDIM - 1],
	};
	x -= 0.5f;
	TEST_T x2 = x*x;
	return scale[0] + x*(scale[1] + x2*(scale[2] + x2*scale[3]));
	//const TEST_T scale[] = {
	//	+2.f    * scale_adj[PDIM - 3],
	//	-15.f/8 * scale_adj[PDIM - 2],
	//	+1.f/2  * scale_adj[PDIM - 1],
	//};
	//TEST_T x2 = x*x;
	//return x2*(scale[0] + x2*(scale[1] + x2*scale[2]));
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

/* -0.5, *2.0 */
static inline TEST_T test_fabs_d12(TEST_T x, double scale_adj[]) {
	const TEST_T scale[] = {
		+1.f * scale_adj[0],
		-1.f * scale_adj[1],
		+1.f * scale_adj[2],
		-1.f * scale_adj[3],
	};
	TEST_T xp = (3 - 2*fabs(x))*x*x;
	return xp*(scale[0] + xp*(scale[1] + xp*(scale[2] + xp*scale[3])));
}
