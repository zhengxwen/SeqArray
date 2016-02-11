// ===========================================================
//
// vectorization.h: compiler optimization with vectorization
//
// Copyright (C) 2016    Xiuwen Zheng
//
// This file is part of CoreArray.
//
// CoreArray is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License Version 3 as
// published by the Free Software Foundation.
//
// CoreArray is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with CoreArray.
// If not, see <http://www.gnu.org/licenses/>.

/**
 *	\file     vectorization.c
 *	\author   Xiuwen Zheng [zhengx@u.washington.edu]
 *	\version  1.0
 *	\date     2016
 *	\brief    compiler optimization with vectorization
 *	\details
**/

#include "vectorization.h"



/// get the number of non-zero
size_t vec_byte_count(const uint8_t *p, size_t n)
{
	size_t ans = 0;

#ifdef COREARRAY_SIMD_SSE2

	const __m128i ZERO = { 0LL, 0LL };
	const __m128i ONES = { 0x0101010101010101LL, 0x0101010101010101LL };
	const __m128i ONE  = { 1LL, 1LL };

	// header 1, 16-byte aligned
	size_t h = (16 - ((size_t)p & 0x0F)) & 0x0F;
	for (; (n > 0) && (h > 0); n--, h--)
		ans += (*p++) ? 1 : 0;

#   ifdef COREARRAY_SIMD_AVX2

	// header 2, 32-byte aligned
	if ((n >= 16) && ((size_t)p & 0x10))
	{
		__m128i c = _mm_cmpeq_epi8(_mm_load_si128((__m128i const*)p), ZERO);
		__m128i bit = _mm_and_si128(c, ONES);
		p += 16; n -= 16;

		uint64_t array[2] __attribute__((aligned(16)));
		*((__m128i*)array) = bit;
		ans += 16 - POPCNT_U64(array[0]) - POPCNT_U64(array[1]);
	}

	const __m256i ZERO2 = { 0LL, 0LL, 0LL, 0LL };
	const __m256i ONES2 = { 0x0101010101010101LL, 0x0101010101010101LL,
							0x0101010101010101LL, 0x0101010101010101LL };

	// body, AVX2
	for (; n >= 256; n -= 256)
	{
		__m256i c = _mm256_cmpeq_epi8(_mm256_load_si256((__m256i const*)p), ZERO2);
		__m256i bit = _mm256_and_si256(c, ONES2);
		p += 32;

		c = _mm256_cmpeq_epi8(_mm256_load_si256((__m256i const*)p), ZERO2);
		bit = _mm256_or_si256(_mm256_sll_epi64(bit, ONE), _mm256_and_si256(c, ONES2));
		p += 32;

		c = _mm256_cmpeq_epi8(_mm256_load_si256((__m256i const*)p), ZERO2);
		bit = _mm256_or_si256(_mm256_sll_epi64(bit, ONE), _mm256_and_si256(c, ONES2));
		p += 32;

		c = _mm256_cmpeq_epi8(_mm256_load_si256((__m256i const*)p), ZERO2);
		bit = _mm256_or_si256(_mm256_sll_epi64(bit, ONE), _mm256_and_si256(c, ONES2));
		p += 32;

		c = _mm256_cmpeq_epi8(_mm256_load_si256((__m256i const*)p), ZERO2);
		bit = _mm256_or_si256(_mm256_sll_epi64(bit, ONE), _mm256_and_si256(c, ONES2));
		p += 32;

		c = _mm256_cmpeq_epi8(_mm256_load_si256((__m256i const*)p), ZERO2);
		bit = _mm256_or_si256(_mm256_sll_epi64(bit, ONE), _mm256_and_si256(c, ONES2));
		p += 32;

		c = _mm256_cmpeq_epi8(_mm256_load_si256((__m256i const*)p), ZERO2);
		bit = _mm256_or_si256(_mm256_sll_epi64(bit, ONE), _mm256_and_si256(c, ONES2));
		p += 32;

		c = _mm256_cmpeq_epi8(_mm256_load_si256((__m256i const*)p), ZERO2);
		bit = _mm256_or_si256(_mm256_sll_epi64(bit, ONE), _mm256_and_si256(c, ONES2));
		p += 32;

		uint64_t array[4] __attribute__((aligned(32)));
		*((__m256i*)array) = bit;
		ans += 256 - POPCNT_U64(array[0]) - POPCNT_U64(array[1]) -
			 POPCNT_U64(array[2]) - POPCNT_U64(array[3]);
	}

#   endif

	// body, SSE2
	for (; n >= 128; n -= 128)
	{
		__m128i c = _mm_cmpeq_epi8(_mm_load_si128((__m128i const*)p), ZERO);
		__m128i bit = _mm_and_si128(c, ONES);
		p += 16;

		c = _mm_cmpeq_epi8(_mm_load_si128((__m128i const*)p), ZERO);
		bit = _mm_or_si128(_mm_sll_epi64(bit, ONE), _mm_and_si128(c, ONES));
		p += 16;

		c = _mm_cmpeq_epi8(_mm_load_si128((__m128i const*)p), ZERO);
		bit = _mm_or_si128(_mm_sll_epi64(bit, ONE), _mm_and_si128(c, ONES));
		p += 16;

		c = _mm_cmpeq_epi8(_mm_load_si128((__m128i const*)p), ZERO);
		bit = _mm_or_si128(_mm_sll_epi64(bit, ONE), _mm_and_si128(c, ONES));
		p += 16;

		c = _mm_cmpeq_epi8(_mm_load_si128((__m128i const*)p), ZERO);
		bit = _mm_or_si128(_mm_sll_epi64(bit, ONE), _mm_and_si128(c, ONES));
		p += 16;

		c = _mm_cmpeq_epi8(_mm_load_si128((__m128i const*)p), ZERO);
		bit = _mm_or_si128(_mm_sll_epi64(bit, ONE), _mm_and_si128(c, ONES));
		p += 16;

		c = _mm_cmpeq_epi8(_mm_load_si128((__m128i const*)p), ZERO);
		bit = _mm_or_si128(_mm_sll_epi64(bit, ONE), _mm_and_si128(c, ONES));
		p += 16;

		c = _mm_cmpeq_epi8(_mm_load_si128((__m128i const*)p), ZERO);
		bit = _mm_or_si128(_mm_sll_epi64(bit, ONE), _mm_and_si128(c, ONES));
		p += 16;

		uint64_t array[2] __attribute__((aligned(16)));
		*((__m128i*)array) = bit;
		ans += 128 - POPCNT_U64(array[0]) - POPCNT_U64(array[1]);
	}

	for (; n >= 16; n -= 16)
	{
		__m128i c = _mm_cmpeq_epi8(_mm_load_si128((__m128i const*)p), ZERO);
		__m128i bit = _mm_and_si128(c, ONES);
		p += 16;
		uint64_t array[2] __attribute__((aligned(16)));
		*((__m128i*)array) = bit;
		ans += 16 - POPCNT_U64(array[0]) - POPCNT_U64(array[1]);
	}

#else

	// header, 8-byte aligned
	size_t h = (8 - ((size_t)p & 0x07)) & 0x07;
	for (; (n > 0) && (h > 0); n--, h--)
		ans += (*p++) ? 1 : 0;
	// body, unroll
	for (; n >= 8; n -= 8)
	{
		ans += (p[0] ? 1 : 0) + (p[1] ? 1 : 0) +
			(p[2] ? 1 : 0) + (p[3] ? 1 : 0) +
			(p[4] ? 1 : 0) + (p[5] ? 1 : 0) +
			(p[6] ? 1 : 0) + (p[7] ? 1 : 0);
		p += 8;
	}

#endif

	// tail
	for (; n > 0; n--) ans += (*p++) ? 1 : 0;

	return ans;
}


void vec_int32_set(int32_t *p, size_t n, int32_t val)
{
	for (; n > 0; n--) *p++ = val;
}


size_t vec_int32_count(int32_t *p, size_t n, int32_t val)
{
	size_t ans = 0;

#ifdef COREARRAY_SIMD_SSE2

	// header 1, 16-byte aligned
	size_t h = ((16 - ((size_t)p & 0x0F)) & 0x0F) >> 2;
	for (; (n > 0) && (h > 0); n--, h--)
		if (*p++ == val) ans++;

	// body, SSE2
	const __m128i mask = _mm_set1_epi32(val);
	__m128i sum = _mm_set1_epi32(0);

#   ifdef COREARRAY_SIMD_AVX2

	// header 2, 32-byte aligned
	if ((n >= 4) && ((size_t)p & 0x10))
	{
		__m128i c = _mm_cmpeq_epi32(_mm_load_si128((__m128i const*)p), mask);
		sum = _mm_sub_epi32(sum, c);
		n -= 4; p += 4;
	}

	const __m256i mask2 = _mm256_set1_epi32(val);
	__m256i sum2 = _mm256_castsi128_si256(sum);

	for (; n >= 8; n-=8, p+=8)
	{
		__m256i c = _mm256_cmpeq_epi32(_mm256_load_si256((__m256i const*)p), mask2);
		sum2 = _mm256_sub_epi32(sum2, c);
	}

	sum = _mm_add_epi32(
		_mm256_extracti128_si256(sum2, 0), _mm256_extracti128_si256(sum2, 1));

#   else

	for (; n >= 4; n-=4, p+=4)
	{
		__m128i c = _mm_cmpeq_epi32(_mm_load_si128((__m128i const*)p), mask);
		sum = _mm_sub_epi32(sum, c);
	}

#   endif

	int32_t array[4] __attribute__((aligned(16)));
	*((__m128i*)array) = sum;
	ans += array[0] + array[1] + array[2] + array[3];

#endif

	// tail
	for (; n > 0; n--) if (*p++ == val) ans++;

	return ans;
}


void vec_int32_count2(int32_t *p, size_t n, int32_t val1, int32_t val2,
	size_t *out_n1, size_t *out_n2)
{
	size_t n1=0, n2=0;

#ifdef COREARRAY_SIMD_SSE2

	// header 1, 16-byte aligned
	size_t h = ((16 - ((size_t)p & 0x0F)) & 0x0F) >> 2;
	for (; (n > 0) && (h > 0); n--, h--)
	{
		int32_t v = *p++;
		if (v == val1) n1++;
		if (v == val2) n2++;
	}

	// body, SSE2
	const __m128i mask1 = _mm_set1_epi32(val1);
	const __m128i mask2 = _mm_set1_epi32(val2);
	__m128i sum1, sum2;
	sum1 = sum2 = _mm_set1_epi32(0);

#   ifdef COREARRAY_SIMD_AVX2

	// header 2, 32-byte aligned
	if ((n >= 4) && ((size_t)p & 0x10))
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		__m128i c1 = _mm_cmpeq_epi32(v, mask1);
		sum1 = _mm_sub_epi32(sum1, c1);
		__m128i c2 = _mm_cmpeq_epi32(v, mask2);
		sum2 = _mm_sub_epi32(sum2, c2);
		n -= 4; p += 4;
	}

	const __m256i mask_1 = _mm256_set1_epi32(val1);
	const __m256i mask_2 = _mm256_set1_epi32(val2);
	__m256i sum_1 = _mm256_castsi128_si256(sum1);
	__m256i sum_2 = _mm256_castsi128_si256(sum2);

	for (; n >= 8; n-=8, p+=8)
	{
		__m256i v = _mm256_load_si256((__m256i const*)p);
		__m256i c1 = _mm256_cmpeq_epi32(v, mask_1);
		sum_1 = _mm256_sub_epi32(sum_1, c1);
		__m256i c2 = _mm256_cmpeq_epi32(v, mask_2);
		sum_2 = _mm256_sub_epi32(sum_2, c2);
	}

	sum1 = _mm_add_epi32(
		_mm256_extracti128_si256(sum_1, 0), _mm256_extracti128_si256(sum_1, 1));
	sum2 = _mm_add_epi32(
		_mm256_extracti128_si256(sum_2, 0), _mm256_extracti128_si256(sum_2, 1));

#   else

	for (; n >= 4; n-=4, p+=4)
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		__m128i c1 = _mm_cmpeq_epi32(v, mask1);
		sum1 = _mm_sub_epi32(sum1, c1);
		__m128i c2 = _mm_cmpeq_epi32(v, mask2);
		sum2 = _mm_sub_epi32(sum2, c2);
	}

#   endif

	int32_t a[4] __attribute__((aligned(16)));
	*((__m128i*)a) = sum1;
	n1 += a[0] + a[1] + a[2] + a[3];
	*((__m128i*)a) = sum2;
	n2 += a[0] + a[1] + a[2] + a[3];

#endif

	// tail
	for (; n > 0; n--)
	{
		int32_t v = *p++;
		if (v == val1) n1++;
		if (v == val2) n2++;
	}

	if (out_n1) *out_n1 = n1;
	if (out_n2) *out_n2 = n2;
}
