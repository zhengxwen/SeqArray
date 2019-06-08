// ===========================================================
//
// vectorization.h: compiler optimization with vectorization
//
// Copyright (C) 2016-2019    Xiuwen Zheng
//
// This file is part of SeqArray.
//
// SeqArray is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as
// published by the Free Software Foundation.
//
// SeqArray is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public
// License along with SeqArray.
// If not, see <http://www.gnu.org/licenses/>.

/**
 *	\file     vectorization.c
 *	\author   Xiuwen Zheng [zhengx@u.washington.edu]
 *	\version  1.0
 *	\date     2016-2019
 *	\brief    compiler optimization with vectorization
 *	\details
**/

#ifndef COREARRAY_COMPILER_OPTIMIZE_FLAG
#   define COREARRAY_COMPILER_OPTIMIZE_FLAG  3
#endif

#include "vectorization.h"


/// get the number of non-zero
size_t vec_i8_cnt_nonzero(const int8_t *p, size_t n)
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


/// get the number of non-zeros and the pointer to the first non-zero value
const int8_t *vec_i8_cnt_nonzero_ptr(const int8_t *p, size_t n, size_t *out_n)
{
#ifdef COREARRAY_SIMD_SSE2

#   ifdef COREARRAY_SIMD_AVX2

	const __m256i ZERO2 = _mm256_setzero_si256();
	for (; n >= 32; n-=32, p+=32)
	{
		__m256i v = _mm256_loadu_si256((__m256i const*)p);
		v = _mm256_cmpeq_epi8(v, ZERO2);
		if (_mm256_movemask_epi8(v) != -1) break;
	}

#   endif

	const __m128i ZERO = _mm_setzero_si128();
	for (; n >= 16; n-=16, p+=16)
	{
		__m128i v = _mm_loadu_si128((__m128i const*)p);
		v = _mm_cmpeq_epi8(v, ZERO);
		if (_mm_movemask_epi8(v) != 0xFFFF) break;
	}

#endif

	// tail
	for (; n > 0; n--, p++) if (*p) break;
	const int8_t *ans = p;

	n = vec_i8_cnt_nonzero(p, n);
	if (out_n) *out_n = n;
	return ans;
}



// ===========================================================
// functions for int8
// ===========================================================

size_t vec_i8_count(const char *p, size_t n, char val)
{
	size_t num = 0;

#ifdef COREARRAY_SIMD_SSE2

	// header 1, 16-byte aligned
	size_t h = (16 - ((size_t)p & 0x0F)) & 0x0F;
	for (; (n > 0) && (h > 0); n--, h--)
		if (*p++ == val) num++;

#   ifdef COREARRAY_SIMD_AVX2
	// body, AVX2
	const __m128i zeros = _mm_setzero_si128();
	const __m256i mask = _mm256_set1_epi8(val);
	__m256i sum = _mm256_setzero_si256();
	size_t offset = 0;

	// header 2, 32-byte aligned
	if ((n >= 16) && ((size_t)p & 0x10))
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		__m128i c1 = _mm_cmpeq_epi8(v, _mm256_castsi256_si128(mask));
		sum = MM_SET_M128(_mm_sub_epi8(zeros, c1), zeros);
		n -= 16; p += 16;
	}

	for (; n >= 128; n-=128)
	{
		__m256i v = _mm256_load_si256((__m256i const*)p); p += 32;
		sum = _mm256_sub_epi8(sum, _mm256_cmpeq_epi8(v, mask));
		v = _mm256_load_si256((__m256i const*)p); p += 32;
		sum = _mm256_sub_epi8(sum, _mm256_cmpeq_epi8(v, mask));
		v = _mm256_load_si256((__m256i const*)p); p += 32;
		sum = _mm256_sub_epi8(sum, _mm256_cmpeq_epi8(v, mask));
		v = _mm256_load_si256((__m256i const*)p); p += 32;
		sum = _mm256_sub_epi8(sum, _mm256_cmpeq_epi8(v, mask));
		offset += 4;
		if (offset >= 252)
		{
			num += vec_avx_sum_u8(sum);
			sum = _mm256_setzero_si256();
			offset = 0;
		}
	}
	for (; n >= 32; n-=32)
	{
		__m256i v = _mm256_load_si256((__m256i const*)p); p += 32;
		sum = _mm256_sub_epi8(sum, _mm256_cmpeq_epi8(v, mask));
		if ((++offset) >= 252)
		{
			num += vec_avx_sum_u8(sum);
			sum = _mm256_setzero_si256();
			offset = 0;
		}
	}
	if (n >= 16)
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		__m128i c1 = _mm_cmpeq_epi8(v, _mm256_castsi256_si128(mask));
		sum = _mm256_sub_epi8(sum, MM_SET_M128(zeros, c1));
		n -= 16; p += 16;
	}

	if (offset > 0)
		num += vec_avx_sum_u8(sum);

#   else
	// body, SSE2
	const __m128i mask = _mm_set1_epi8(val);
	__m128i sum = _mm_setzero_si128();
	size_t offset = 0;

	for (; n >= 64; n-=64)
	{
		__m128i v = _mm_load_si128((__m128i const*)p); p += 16;
		sum = _mm_sub_epi8(sum, _mm_cmpeq_epi8(v, mask));
		v = _mm_load_si128((__m128i const*)p); p += 16;
		sum = _mm_sub_epi8(sum, _mm_cmpeq_epi8(v, mask));
		v = _mm_load_si128((__m128i const*)p); p += 16;
		sum = _mm_sub_epi8(sum, _mm_cmpeq_epi8(v, mask));
		v = _mm_load_si128((__m128i const*)p); p += 16;
		sum = _mm_sub_epi8(sum, _mm_cmpeq_epi8(v, mask));
		offset += 4;
		if (offset >= 252)
		{
			num += vec_sum_u8(sum);
			sum = _mm_setzero_si128();
			offset = 0;
		}
	}
	for (; n >= 16; n-=16, p+=16)
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		sum = _mm_sub_epi8(sum, _mm_cmpeq_epi8(v, mask));
		if ((++offset) >= 252)
		{
			num += vec_sum_u8(sum);
			sum = _mm_setzero_si128();
			offset = 0;
		}
	}

	if (offset > 0)
		num += vec_sum_u8(sum);
#endif

#endif

	// tail
	for (; n > 0; n--)
		if (*p++ == val) num++;
	return num;
}


void vec_i8_count2(const char *p, size_t n, char val1, char val2,
	size_t *out_n1, size_t *out_n2)
{
	size_t n1 = 0, n2 = 0;

#ifdef COREARRAY_SIMD_SSE2

	// header 1, 16-byte aligned
	size_t h = (16 - ((size_t)p & 0x0F)) & 0x0F;
	for (; (n > 0) && (h > 0); n--, h--)
	{
		char v = *p++;
		if (v == val1) n1++;
		if (v == val2) n2++;
	}

#   ifdef COREARRAY_SIMD_AVX2
	// body, AVX2
	const __m128i zeros = _mm_setzero_si128();
	const __m256i mask1 = _mm256_set1_epi8(val1);
	const __m256i mask2 = _mm256_set1_epi8(val2);
	__m256i sum1, sum2;
	sum1 = sum2 = _mm256_setzero_si256();
	size_t offset = 0;

	// header 2, 32-byte aligned
	if ((n >= 16) && ((size_t)p & 0x10))
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		__m128i c1 = _mm_cmpeq_epi8(v, _mm256_castsi256_si128(mask1));
		sum1 = MM_SET_M128(_mm_sub_epi8(zeros, c1), zeros);
		__m128i c2 = _mm_cmpeq_epi8(v, _mm256_castsi256_si128(mask2));
		sum2 = MM_SET_M128(_mm_sub_epi8(zeros, c2), zeros);
		n -= 16; p += 16;
	}

	for (; n >= 32; n-=32, p+=32)
	{
		__m256i v = _mm256_load_si256((__m256i const*)p);
		sum1 = _mm256_sub_epi8(sum1, _mm256_cmpeq_epi8(v, mask1));
		sum2 = _mm256_sub_epi8(sum2, _mm256_cmpeq_epi8(v, mask2));
		if ((++offset) >= 252)
		{
			n1 += vec_avx_sum_u8(sum1);
			n2 += vec_avx_sum_u8(sum2);
			sum1 = sum2 = _mm256_setzero_si256();
			offset = 0;
		}
	}

	if (n >= 16)
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		__m128i c1 = _mm_cmpeq_epi8(v, _mm256_castsi256_si128(mask1));
		sum1 = _mm256_sub_epi8(sum1, MM_SET_M128(c1, zeros));
		__m128i c2 = _mm_cmpeq_epi8(v, _mm256_castsi256_si128(mask2));
		sum2 = _mm256_sub_epi8(sum2, MM_SET_M128(c2, zeros));
		n -= 16; p += 16;
	}

	if (offset > 0)
	{
		n1 += vec_avx_sum_u8(sum1);
		n2 += vec_avx_sum_u8(sum2);
	}

#   else
	// body, SSE2
	const __m128i mask1 = _mm_set1_epi8(val1);
	const __m128i mask2 = _mm_set1_epi8(val2);
	__m128i sum1, sum2;
	sum1 = sum2 = _mm_setzero_si128();
	size_t offset = 0;

	for (; n >= 16; n-=16, p+=16)
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		sum1 = _mm_sub_epi8(sum1, _mm_cmpeq_epi8(v, mask1));
		sum2 = _mm_sub_epi8(sum2, _mm_cmpeq_epi8(v, mask2));
		if ((++offset) >= 252)
		{
			n1 += vec_sum_u8(sum1);
			n2 += vec_sum_u8(sum2);
			sum1 = sum2 = _mm_setzero_si128();
			offset = 0;
		}
	}

	if (offset > 0)
	{
		n1 += vec_sum_u8(sum1);
		n2 += vec_sum_u8(sum2);
	}
#endif

#endif

	// tail
	for (; n > 0; n--)
	{
		char v = *p++;
		if (v == val1) n1++;
		if (v == val2) n2++;
	}

	if (out_n1) *out_n1 = n1;
	if (out_n2) *out_n2 = n2;
}


void vec_i8_count3(const char *p, size_t n, char val1, char val2, char val3,
	size_t *out_n1, size_t *out_n2, size_t *out_n3)
{
	size_t n1 = 0, n2 = 0, n3 = 0;

#ifdef COREARRAY_SIMD_SSE2

	// header 1, 16-byte aligned
	size_t h = (16 - ((size_t)p & 0x0F)) & 0x0F;
	for (; (n > 0) && (h > 0); n--, h--)
	{
		char v = *p++;
		if (v == val1) n1++;
		if (v == val2) n2++;
		if (v == val3) n3++;
	}

#   ifdef COREARRAY_SIMD_AVX2
	// body, AVX2
	const __m128i zeros = _mm_setzero_si128();
	const __m256i mask1 = _mm256_set1_epi8(val1);
	const __m256i mask2 = _mm256_set1_epi8(val2);
	const __m256i mask3 = _mm256_set1_epi8(val3);
	__m256i sum1, sum2, sum3;
	sum1 = sum2 = sum3 = _mm256_setzero_si256();
	size_t offset = 0;

	// header 2, 32-byte aligned
	if ((n >= 16) && ((size_t)p & 0x10))
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		__m128i c1 = _mm_cmpeq_epi8(v, _mm256_castsi256_si128(mask1));
		sum1 = MM_SET_M128(_mm_sub_epi8(zeros, c1), zeros);
		__m128i c2 = _mm_cmpeq_epi8(v, _mm256_castsi256_si128(mask2));
		sum2 = MM_SET_M128(_mm_sub_epi8(zeros, c2), zeros);
		__m128i c3 = _mm_cmpeq_epi8(v, _mm256_castsi256_si128(mask3));
		sum3 = MM_SET_M128(_mm_sub_epi8(zeros, c3), zeros);
		n -= 16; p += 16;
	}

	for (; n >= 32; n-=32, p+=32)
	{
		__m256i v = _mm256_load_si256((__m256i const*)p);
		sum1 = _mm256_sub_epi8(sum1, _mm256_cmpeq_epi8(v, mask1));
		sum2 = _mm256_sub_epi8(sum2, _mm256_cmpeq_epi8(v, mask2));
		sum3 = _mm256_sub_epi8(sum3, _mm256_cmpeq_epi8(v, mask3));
		if ((++offset) >= 252)
		{
			n1 += vec_avx_sum_u8(sum1);
			n2 += vec_avx_sum_u8(sum2);
			n3 += vec_avx_sum_u8(sum3);
			sum1 = sum2 = sum3 = _mm256_setzero_si256();
			offset = 0;
		}
	}

	if (n >= 16)
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		__m128i c1 = _mm_cmpeq_epi8(v, _mm256_castsi256_si128(mask1));
		sum1 = _mm256_sub_epi8(sum1, MM_SET_M128(c1, zeros));
		__m128i c2 = _mm_cmpeq_epi8(v, _mm256_castsi256_si128(mask2));
		sum2 = _mm256_sub_epi8(sum2, MM_SET_M128(c2, zeros));
		__m128i c3 = _mm_cmpeq_epi8(v, _mm256_castsi256_si128(mask3));
		sum3 = _mm256_sub_epi8(sum3, MM_SET_M128(c3, zeros));
		n -= 16; p += 16;
	}

	if (offset > 0)
	{
		n1 += vec_avx_sum_u8(sum1);
		n2 += vec_avx_sum_u8(sum2);
		n3 += vec_avx_sum_u8(sum3);
	}

#   else
	// body, SSE2
	const __m128i mask1 = _mm_set1_epi8(val1);
	const __m128i mask2 = _mm_set1_epi8(val2);
	const __m128i mask3 = _mm_set1_epi8(val3);
	__m128i sum1, sum2, sum3;
	sum1 = sum2 = sum3 = _mm_setzero_si128();
	size_t offset = 0;

	for (; n >= 16; n-=16, p+=16)
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		sum1 = _mm_sub_epi8(sum1, _mm_cmpeq_epi8(v, mask1));
		sum2 = _mm_sub_epi8(sum2, _mm_cmpeq_epi8(v, mask2));
		sum3 = _mm_sub_epi8(sum3, _mm_cmpeq_epi8(v, mask3));
		if ((++offset) >= 252)
		{
			n1 += vec_sum_u8(sum1);
			n2 += vec_sum_u8(sum2);
			n3 += vec_sum_u8(sum3);
			sum1 = sum2 = sum3 = _mm_setzero_si128();
			offset = 0;
		}
	}

	if (offset > 0)
	{
		n1 += vec_sum_u8(sum1);
		n2 += vec_sum_u8(sum2);
		n3 += vec_sum_u8(sum3);
	}
#endif

#endif

	// tail
	for (; n > 0; n--)
	{
		char v = *p++;
		if (v == val1) n1++;
		if (v == val2) n2++;
		if (v == val3) n3++;
	}

	if (out_n1) *out_n1 = n1;
	if (out_n2) *out_n2 = n2;
	if (out_n3) *out_n3 = n3;
}


void vec_i8_replace(int8_t *p, size_t n, int8_t val, int8_t substitute)
{
#ifdef COREARRAY_SIMD_SSE2

	// header 1, 16-byte aligned
	size_t h = (16 - ((size_t)p & 0x0F)) & 0x0F;
	for (; (n > 0) && (h > 0); n--, h--, p++)
		if (*p == val) *p = substitute;

	// body, SSE2
	const __m128i mask = _mm_set1_epi8(val);
	const __m128i sub  = _mm_set1_epi8(substitute);

#   ifdef COREARRAY_SIMD_AVX2

	// header 2, 32-byte aligned
	if ((n >= 16) && ((size_t)p & 0x10))
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		__m128i c = _mm_cmpeq_epi8(v, mask);
		if (_mm_movemask_epi8(c))
		{
			_mm_store_si128((__m128i *)p,
				_mm_or_si128(_mm_and_si128(c, sub), _mm_andnot_si128(c, v)));
		}
		n -= 16; p += 16;
	}

	const __m256i mask2 = _mm256_set1_epi8(val);
	const __m256i sub32 = _mm256_set1_epi8(substitute);
	const __m256i zero = _mm256_setzero_si256();
	const __m256i ones = _mm256_cmpeq_epi64(zero, zero);

	for (; n >= 32; n-=32, p+=32)
	{
		__m256i v = _mm256_load_si256((__m256i const*)p);
		__m256i c = _mm256_cmpeq_epi8(v, mask2);
		if (_mm256_movemask_epi8(c))
		{
			// TODO
			_mm256_store_si256((__m256i *)p,
				_mm256_or_si256(_mm256_and_si256(c, sub32),
				_mm256_andnot_si256(c, v)));
		}
	}

#   endif

	for (; n >= 16; n-=16, p+=16)
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		__m128i c = _mm_cmpeq_epi8(v, mask);
		if (_mm_movemask_epi8(c))
			_mm_maskmoveu_si128(sub, c, (char*)p);
	}

#endif

	// tail
	for (; n > 0; n--, p++)
		if (*p == val) *p = substitute;
}


void vec_i8_cnt_dosage2(const int8_t *p, int8_t *out, size_t n, int8_t val,
	int8_t missing, int8_t missing_substitute)
{
#ifdef COREARRAY_SIMD_SSE2

	// header 1, 16-byte aligned
	size_t h = (16 - ((size_t)out & 0x0F)) & 0x0F;
	for (; (n > 0) && (h > 0); n--, h--, p+=2)
	{
		*out ++ = ((p[0] == missing) || (p[1] == missing)) ?
			missing_substitute :
			(p[0]==val ? 1 : 0) + (p[1]==val ? 1 : 0);
	}

	// body, SSE2
	const __m128i val16  = _mm_set1_epi8(val);
	const __m128i miss16 = _mm_set1_epi8(missing);
	const __m128i sub16  = _mm_set1_epi8(missing_substitute);
	const __m128i mask   = _mm_set1_epi16(0x00FF);

#   ifdef COREARRAY_SIMD_AVX2

	// header 2, 32-byte aligned
	if ((n >= 16) && ((size_t)out & 0x10))
	{
		__m128i w1 = MM_LOADU_128((__m128i const*)p); p += 16;
		__m128i w2 = MM_LOADU_128((__m128i const*)p); p += 16;

		__m128i v1 = _mm_packus_epi16(_mm_and_si128(w1, mask), _mm_and_si128(w2, mask));
		__m128i v2 = _mm_packus_epi16(_mm_srli_epi16(w1, 8), _mm_srli_epi16(w2, 8));

		__m128i c = _mm_setzero_si128();
		c = _mm_sub_epi8(c, _mm_cmpeq_epi8(v1, val16));
		c = _mm_sub_epi8(c, _mm_cmpeq_epi8(v2, val16));

		w1 = _mm_cmpeq_epi8(v1, miss16);
		w2 = _mm_cmpeq_epi8(v2, miss16);
		__m128i w  = _mm_or_si128(w1, w2);
		c = _mm_or_si128(_mm_and_si128(w, sub16), _mm_andnot_si128(w, c));

		_mm_store_si128((__m128i *)out, c);
		n -= 16; out += 16;
	}

	const __m256i val32  = _mm256_set1_epi8(val);
	const __m256i miss32 = _mm256_set1_epi8(missing);
	const __m256i sub32  = _mm256_set1_epi8(missing_substitute);
	const __m256i mask2  = _mm256_set1_epi16(0x00FF);

	for (; n >= 32; n-=32)
	{
		__m256i w1 = MM_LOADU_256((__m256i const*)p); p += 32;
		__m256i w2 = MM_LOADU_256((__m256i const*)p); p += 32;

		__m256i v1 = _mm256_packus_epi16(_mm256_and_si256(w1, mask2), _mm256_and_si256(w2, mask2));
		__m256i v2 = _mm256_packus_epi16(_mm256_srli_epi16(w1, 8), _mm256_srli_epi16(w2, 8));

		__m256i c = _mm256_setzero_si256();
		c = _mm256_sub_epi8(c, _mm256_cmpeq_epi8(v1, val32));
		c = _mm256_sub_epi8(c, _mm256_cmpeq_epi8(v2, val32));

		w1 = _mm256_cmpeq_epi8(v1, miss32);
		w2 = _mm256_cmpeq_epi8(v2, miss32);
		__m256i w = _mm256_or_si256(w1, w2);
		c = _mm256_or_si256(_mm256_and_si256(w, sub32), _mm256_andnot_si256(w, c));

		c = _mm256_permute4x64_epi64(c, 0xD8);
		_mm256_store_si256((__m256i *)out, c);
		out += 32;
	}

#   endif

	// SSE2 only
	for (; n >= 16; n-=16)
	{
		__m128i w1 = MM_LOADU_128((__m128i const*)p); p += 16;
		__m128i w2 = MM_LOADU_128((__m128i const*)p); p += 16;

		__m128i v1 = _mm_packus_epi16(_mm_and_si128(w1, mask), _mm_and_si128(w2, mask));
		__m128i v2 = _mm_packus_epi16(_mm_srli_epi16(w1, 8), _mm_srli_epi16(w2, 8));

		__m128i c = _mm_setzero_si128();
		c = _mm_sub_epi8(c, _mm_cmpeq_epi8(v1, val16));
		c = _mm_sub_epi8(c, _mm_cmpeq_epi8(v2, val16));

		w1 = _mm_cmpeq_epi8(v1, miss16);
		w2 = _mm_cmpeq_epi8(v2, miss16);
		__m128i w = _mm_or_si128(w1, w2);
		c = _mm_or_si128(_mm_and_si128(w, sub16), _mm_andnot_si128(w, c));

		_mm_store_si128((__m128i *)out, c);
		out += 16;
	}

#endif

	// tail
	for (; n > 0; n--, p+=2)
	{
		*out ++ = ((p[0] == missing) || (p[1] == missing)) ?
			missing_substitute :
			(p[0]==val ? 1 : 0) + (p[1]==val ? 1 : 0);
	}
}


void vec_i8_cnt_dosage_alt2(const int8_t *p, int8_t *out, size_t n, int8_t val,
	int8_t missing, int8_t missing_substitute)
{
#ifdef COREARRAY_SIMD_SSE2

	// header 1, 16-byte aligned
	size_t h = (16 - ((size_t)out & 0x0F)) & 0x0F;
	for (; (n > 0) && (h > 0); n--, h--, p+=2)
	{
		*out ++ = ((p[0] == missing) || (p[1] == missing)) ?
			missing_substitute :
			(p[0]!=val ? 1 : 0) + (p[1]!=val ? 1 : 0);
	}

	// body, SSE2
	const __m128i val16  = _mm_set1_epi8(val);
	const __m128i miss16 = _mm_set1_epi8(missing);
	const __m128i sub16  = _mm_set1_epi8(missing_substitute);
	const __m128i two16  = _mm_set1_epi8(2);
	const __m128i mask   = _mm_set1_epi16(0x00FF);

#   ifdef COREARRAY_SIMD_AVX2

	// header 2, 32-byte aligned
	if ((n >= 16) && ((size_t)out & 0x10))
	{
		__m128i w1 = MM_LOADU_128((__m128i const*)p); p += 16;
		__m128i w2 = MM_LOADU_128((__m128i const*)p); p += 16;

		__m128i v1 = _mm_packus_epi16(_mm_and_si128(w1, mask), _mm_and_si128(w2, mask));
		__m128i v2 = _mm_packus_epi16(_mm_srli_epi16(w1, 8), _mm_srli_epi16(w2, 8));

		__m128i c = two16;
		c = _mm_add_epi8(c, _mm_cmpeq_epi8(v1, val16));
		c = _mm_add_epi8(c, _mm_cmpeq_epi8(v2, val16));

		w1 = _mm_cmpeq_epi8(v1, miss16);
		w2 = _mm_cmpeq_epi8(v2, miss16);
		__m128i w  = _mm_or_si128(w1, w2);
		c = _mm_or_si128(_mm_and_si128(w, sub16), _mm_andnot_si128(w, c));

		_mm_store_si128((__m128i *)out, c);
		n -= 16; out += 16;
	}

	const __m256i val32  = _mm256_set1_epi8(val);
	const __m256i miss32 = _mm256_set1_epi8(missing);
	const __m256i sub32  = _mm256_set1_epi8(missing_substitute);
	const __m256i two32  = _mm256_set1_epi8(2);
	const __m256i mask2  = _mm256_set1_epi16(0x00FF);

	for (; n >= 32; n-=32)
	{
		__m256i w1 = MM_LOADU_256((__m256i const*)p); p += 32;
		__m256i w2 = MM_LOADU_256((__m256i const*)p); p += 32;

		__m256i v1 = _mm256_packus_epi16(_mm256_and_si256(w1, mask2), _mm256_and_si256(w2, mask2));
		__m256i v2 = _mm256_packus_epi16(_mm256_srli_epi16(w1, 8), _mm256_srli_epi16(w2, 8));

		__m256i c = two32;
		c = _mm256_add_epi8(c, _mm256_cmpeq_epi8(v1, val32));
		c = _mm256_add_epi8(c, _mm256_cmpeq_epi8(v2, val32));

		w1 = _mm256_cmpeq_epi8(v1, miss32);
		w2 = _mm256_cmpeq_epi8(v2, miss32);
		__m256i w = _mm256_or_si256(w1, w2);
		c = _mm256_or_si256(_mm256_and_si256(w, sub32), _mm256_andnot_si256(w, c));

		c = _mm256_permute4x64_epi64(c, 0xD8);
		_mm256_store_si256((__m256i *)out, c);
		out += 32;
	}

#   endif

	// SSE2 only
	for (; n >= 16; n-=16)
	{
		__m128i w1 = MM_LOADU_128((__m128i const*)p); p += 16;
		__m128i w2 = MM_LOADU_128((__m128i const*)p); p += 16;

		__m128i v1 = _mm_packus_epi16(_mm_and_si128(w1, mask), _mm_and_si128(w2, mask));
		__m128i v2 = _mm_packus_epi16(_mm_srli_epi16(w1, 8), _mm_srli_epi16(w2, 8));

		__m128i c = two16;
		c = _mm_add_epi8(c, _mm_cmpeq_epi8(v1, val16));
		c = _mm_add_epi8(c, _mm_cmpeq_epi8(v2, val16));

		w1 = _mm_cmpeq_epi8(v1, miss16);
		w2 = _mm_cmpeq_epi8(v2, miss16);
		__m128i w = _mm_or_si128(w1, w2);
		c = _mm_or_si128(_mm_and_si128(w, sub16), _mm_andnot_si128(w, c));

		_mm_store_si128((__m128i *)out, c);
		out += 16;
	}

#endif

	// tail
	for (; n > 0; n--, p+=2)
	{
		*out ++ = ((p[0] == missing) || (p[1] == missing)) ?
			missing_substitute :
			(p[0]!=val ? 1 : 0) + (p[1]!=val ? 1 : 0);
	}
}



// ===========================================================
// functions for uint8
// ===========================================================

/// shifting *p right by 2 bits, assuming p is 2-byte aligned
void vec_u8_shr_b2(uint8_t *p, size_t n)
{
#ifdef COREARRAY_SIMD_SSE2

	// header 1, 16-byte aligned
	size_t h = ((16 - ((size_t)p & 0x0F)) & 0x0F);
	for (; (n > 0) && (h > 0); n--, h--)
		*p++ >>= 2;

	// body, SSE2
	const __m128i mask = _mm_set1_epi8(0x3F);
	for (; n >= 16; n-=16, p+=16)
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		v = _mm_srli_epi16(v, 2);
		_mm_store_si128((__m128i *)p, v & mask);
	}

#endif

	// tail
	for (; n > 0; n--) *p++ >>= 2;
}



// ===========================================================
// functions for int16
// ===========================================================

/// shifting *p right by 2 bits, assuming p is 2-byte aligned
void vec_i16_shr_b2(int16_t *p, size_t n)
{
#ifdef COREARRAY_SIMD_SSE2

	// header 1, 16-byte aligned
	size_t h = ((16 - ((size_t)p & 0x0F)) & 0x0F) >> 1;
	for (; (n > 0) && (h > 0); n--, h--)
		*p++ >>= 2;

	// body, SSE2
	for (; n >= 8; n-=8, p+=8)
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		_mm_store_si128((__m128i *)p, _mm_srli_epi16(v, 2));
	}

#endif

	// tail
	for (; n > 0; n--) *p++ >>= 2;
}



// ===========================================================
// functions for int32
// ===========================================================

/// count how many val in p, assuming p is 4-byte aligned
size_t vec_i32_count(const int32_t *p, size_t n, int32_t val)
{
	size_t ans = 0;

#ifdef __LP64__
	if (n > 2147483632) // 2^31 - 16
	{
		while (n > 0)
		{
			size_t nn = (n <= 2147483632) ? n : 2147483632;
			ans += vec_i32_count(p, nn, val);
			p += nn; n -= nn;
		}
		return ans;
	}
#endif

#ifdef COREARRAY_SIMD_SSE2

	// header 1, 16-byte aligned
	size_t h = ((16 - ((size_t)p & 0x0F)) & 0x0F) >> 2;
	for (; (n > 0) && (h > 0); n--, h--)
		if (*p++ == val) ans++;

#   ifdef COREARRAY_SIMD_AVX2

	// body, AVX2
	const __m128i zero = _mm_setzero_si128();
	const __m256i mask = _mm256_set1_epi32(val);
	__m256i sum = _mm256_setzero_si256();

	// header 2, 32-byte aligned
	if ((n >= 4) && ((size_t)p & 0x10))
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		__m128i c = _mm_cmpeq_epi32(v, _mm256_castsi256_si128(mask));
		sum = _mm256_sub_epi32(sum, MM_SET_M128(zero, c));
		n -= 4; p += 4;
	}
	for (; n >= 8; n-=8, p+=8)
	{
		__m256i c = _mm256_cmpeq_epi32(_mm256_load_si256((__m256i*)p), mask);
		sum = _mm256_sub_epi32(sum, c);
	}
	if (n >= 4)
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		__m128i c = _mm_cmpeq_epi32(v, _mm256_castsi256_si128(mask));
		sum = _mm256_sub_epi32(sum, MM_SET_M128(zero, c));
		n -= 4; p += 4;
	}
	ans += vec_avx_sum_i32(sum);

#   else

	// body, SSE2
	const __m128i mask = _mm_set1_epi32(val);
	__m128i sum = _mm_setzero_si128();
	for (; n >= 4; n-=4, p+=4)
	{
		__m128i c = _mm_cmpeq_epi32(_mm_load_si128((__m128i const*)p), mask);
		sum = _mm_sub_epi32(sum, c);
	}
	ans += vec_sum_i32(sum);

#   endif

#endif

	// tail
	for (; n > 0; n--) if (*p++ == val) ans++;

	return ans;
}


/// count how many val1 and val2 in p, assuming p is 4-byte aligned
void vec_i32_count2(const int32_t *p, size_t n, int32_t val1, int32_t val2,
	size_t *out_n1, size_t *out_n2)
{
	size_t n1 = 0, n2 = 0;

#ifdef __LP64__
	if (n > 2147483632) // 2^31 - 16
	{
		size_t m1 = 0, m2 = 0;
		while (n > 0)
		{
			size_t nn = (n <= 2147483632) ? n : 2147483632;
			vec_i32_count2(p, nn, val1, val2, &m1, &m2);
			n1 += m1; n2 += m2;
			p += nn; n -= nn;
		}
		if (out_n1) *out_n1 = n1;
		if (out_n2) *out_n2 = n2;
		return;
	}
#endif

#ifdef COREARRAY_SIMD_SSE2

	// header 1, 16-byte aligned
	size_t h = ((16 - ((size_t)p & 0x0F)) & 0x0F) >> 2;
	for (; (n > 0) && (h > 0); n--, h--)
	{
		int32_t v = *p++;
		if (v == val1) n1++;
		if (v == val2) n2++;
	}

#   ifdef COREARRAY_SIMD_AVX2

	// body, AVX2
	const __m128i zero  = _mm_setzero_si128();
	const __m256i mask1 = _mm256_set1_epi32(val1);
	const __m256i mask2 = _mm256_set1_epi32(val2);
	__m256i sum1, sum2;
	sum1 = sum2 = _mm256_setzero_si256();

	// header 2, 32-byte aligned
	if ((n >= 4) && ((size_t)p & 0x10))
	{
		__m128i v = _mm_load_si128((__m128i*)p);
		__m128i c1 = _mm_cmpeq_epi32(v, _mm256_castsi256_si128(mask1));
		sum1 = _mm256_sub_epi32(sum1, MM_SET_M128(zero, c1));
		__m128i c2 = _mm_cmpeq_epi32(v, _mm256_castsi256_si128(mask2));
		sum2 = _mm256_sub_epi32(sum2, MM_SET_M128(zero, c2));
		n -= 4; p += 4;
	}
	for (; n >= 8; n-=8, p+=8)
	{
		__m256i v = _mm256_load_si256((__m256i const*)p);
		sum1 = _mm256_sub_epi32(sum1, _mm256_cmpeq_epi32(v, mask1));
		sum2 = _mm256_sub_epi32(sum2, _mm256_cmpeq_epi32(v, mask2));
	}
	if (n >= 4)
	{
		__m128i v = _mm_load_si128((__m128i*)p);
		__m128i c1 = _mm_cmpeq_epi32(v, _mm256_castsi256_si128(mask1));
		sum1 = _mm256_sub_epi32(sum1, MM_SET_M128(zero, c1));
		__m128i c2 = _mm_cmpeq_epi32(v, _mm256_castsi256_si128(mask2));
		sum2 = _mm256_sub_epi32(sum2, MM_SET_M128(zero, c2));
		n -= 4; p += 4;
	}

	n1 += vec_avx_sum_i32(sum1);
	n2 += vec_avx_sum_i32(sum2);

#   else

	// body, SSE2
	const __m128i mask1 = _mm_set1_epi32(val1);
	const __m128i mask2 = _mm_set1_epi32(val2);
	__m128i sum1, sum2;
	sum1 = sum2 = _mm_setzero_si128();

	for (; n >= 4; n-=4, p+=4)
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		__m128i c1 = _mm_cmpeq_epi32(v, mask1);
		sum1 = _mm_sub_epi32(sum1, c1);
		__m128i c2 = _mm_cmpeq_epi32(v, mask2);
		sum2 = _mm_sub_epi32(sum2, c2);
	}

	n1 += vec_sum_i32(sum1);
	n2 += vec_sum_i32(sum2);

#   endif

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


/// count how many val1, val2 and val3 in p, assuming p is 4-byte aligned
void vec_i32_count3(const int32_t *p, size_t n, int32_t val1, int32_t val2,
	int32_t val3, size_t *out_n1, size_t *out_n2, size_t *out_n3)
{
	size_t n1 = 0, n2 = 0, n3 = 0;

#ifdef __LP64__
	if (n > 2147483632) // 2^31 - 16
	{
		size_t m1 = 0, m2 = 0, m3 = 0;
		while (n > 0)
		{
			size_t nn = (n <= 2147483632) ? n : 2147483632;
			vec_i32_count3(p, nn, val1, val2, val3, &m1, &m2, &m3);
			n1 += m1; n2 += m2; n3 += m3;
			p += nn; n -= nn;
		}
		if (out_n1) *out_n1 = n1;
		if (out_n2) *out_n2 = n2;
		if (out_n3) *out_n3 = n3;
		return;
	}
#endif

#ifdef COREARRAY_SIMD_SSE2

	// header 1, 16-byte aligned
	size_t h = ((16 - ((size_t)p & 0x0F)) & 0x0F) >> 2;
	for (; (n > 0) && (h > 0); n--, h--)
	{
		int32_t v = *p++;
		if (v == val1) n1++;
		if (v == val2) n2++;
		if (v == val3) n3++;
	}

#   ifdef COREARRAY_SIMD_AVX2

	// body, AVX2
	const __m128i zero  = _mm_setzero_si128();
	const __m256i mask1 = _mm256_set1_epi32(val1);
	const __m256i mask2 = _mm256_set1_epi32(val2);
	const __m256i mask3 = _mm256_set1_epi32(val3);
	__m256i sum1, sum2, sum3;
	sum1 = sum2 = sum3 = _mm256_setzero_si256();

	// header 2, 32-byte aligned
	if ((n >= 4) && ((size_t)p & 0x10))
	{
		__m128i v = _mm_load_si128((__m128i*)p);
		__m128i c1 = _mm_cmpeq_epi32(v, _mm256_castsi256_si128(mask1));
		sum1 = _mm256_sub_epi32(sum1, MM_SET_M128(zero, c1));
		__m128i c2 = _mm_cmpeq_epi32(v, _mm256_castsi256_si128(mask2));
		sum2 = _mm256_sub_epi32(sum2, MM_SET_M128(zero, c2));
		__m128i c3 = _mm_cmpeq_epi32(v, _mm256_castsi256_si128(mask3));
		sum3 = _mm256_sub_epi32(sum3, MM_SET_M128(zero, c3));
		n -= 4; p += 4;
	}
	for (; n >= 8; n-=8, p+=8)
	{
		__m256i v = _mm256_load_si256((__m256i const*)p);
		sum1 = _mm256_sub_epi32(sum1, _mm256_cmpeq_epi32(v, mask1));
		sum2 = _mm256_sub_epi32(sum2, _mm256_cmpeq_epi32(v, mask2));
		sum3 = _mm256_sub_epi32(sum3, _mm256_cmpeq_epi32(v, mask3));
	}
	if (n >= 4)
	{
		__m128i v = _mm_load_si128((__m128i*)p);
		__m128i c1 = _mm_cmpeq_epi32(v, _mm256_castsi256_si128(mask1));
		sum1 = _mm256_sub_epi32(sum1, MM_SET_M128(zero, c1));
		__m128i c2 = _mm_cmpeq_epi32(v, _mm256_castsi256_si128(mask2));
		sum2 = _mm256_sub_epi32(sum2, MM_SET_M128(zero, c2));
		__m128i c3 = _mm_cmpeq_epi32(v, _mm256_castsi256_si128(mask3));
		sum3 = _mm256_sub_epi32(sum3, MM_SET_M128(zero, c3));
		n -= 4; p += 4;
	}

	n1 += vec_avx_sum_i32(sum1);
	n2 += vec_avx_sum_i32(sum2);
	n3 += vec_avx_sum_i32(sum3);

#   else

	// body, SSE2
	const __m128i mask1 = _mm_set1_epi32(val1);
	const __m128i mask2 = _mm_set1_epi32(val2);
	const __m128i mask3 = _mm_set1_epi32(val3);
	__m128i sum1, sum2, sum3;
	sum1 = sum2 = sum3 = _mm_setzero_si128();

	for (; n >= 4; n-=4, p+=4)
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		__m128i c1 = _mm_cmpeq_epi32(v, mask1);
		sum1 = _mm_sub_epi32(sum1, c1);
		__m128i c2 = _mm_cmpeq_epi32(v, mask2);
		sum2 = _mm_sub_epi32(sum2, c2);
		__m128i c3 = _mm_cmpeq_epi32(v, mask3);
		sum3 = _mm_sub_epi32(sum3, c3);
	}

	n1 += vec_sum_i32(sum1);
	n2 += vec_sum_i32(sum2);
	n3 += vec_sum_i32(sum3);

#   endif

#endif

	// tail
	for (; n > 0; n--)
	{
		int32_t v = *p++;
		if (v == val1) n1++;
		if (v == val2) n2++;
		if (v == val3) n3++;
	}

	if (out_n1) *out_n1 = n1;
	if (out_n2) *out_n2 = n2;
	if (out_n3) *out_n3 = n3;
}


void vec_int32_set(int32_t *p, size_t n, int32_t val)
{
	for (; n > 0; n--) *p++ = val;
}


/// replace 'val' in the array of 'p' by 'substitute', assuming 'p' is 4-byte aligned
void vec_i32_replace(int32_t *p, size_t n, int32_t val, int32_t substitute)
{
#ifdef COREARRAY_SIMD_SSE2

	// header 1, 16-byte aligned
	size_t h = ((16 - ((size_t)p & 0x0F)) & 0x0F) >> 2;
	for (; (n > 0) && (h > 0); n--, h--, p++)
		if (*p == val) *p = substitute;

	// body, SSE2
	const __m128i mask = _mm_set1_epi32(val);
	const __m128i sub4 = _mm_set1_epi32(substitute);

#   ifdef COREARRAY_SIMD_AVX2

	// header 2, 32-byte aligned
	if ((n >= 4) && ((size_t)p & 0x10))
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		__m128i c = _mm_cmpeq_epi32(v, mask);
		if (_mm_movemask_epi8(c))
			_mm_maskstore_epi32(p, c, sub4);
		n -= 4; p += 4;
	}

	const __m256i mask2 = _mm256_set1_epi32(val);
	const __m256i sub8  = _mm256_set1_epi32(substitute);

	for (; n >= 8; n-=8, p+=8)
	{
		__m256i v = _mm256_load_si256((__m256i const*)p);
		__m256i c = _mm256_cmpeq_epi32(v, mask2);
		if (_mm256_movemask_epi8(c))
			_mm256_maskstore_epi32(p, c, sub8);
	}

#   endif

	for (; n >= 4; n-=4, p+=4)
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		__m128i c = _mm_cmpeq_epi32(v, mask);
		if (_mm_movemask_epi8(c))
			_mm_maskmoveu_si128(sub4, c, (char*)p);
	}

#endif

	// tail
	for (; n > 0; n--, p++)
		if (*p == val) *p = substitute;
}


/// assuming 'out' is 4-byte aligned, output (p[0]==val) + (p[1]==val) or missing_substitute
void vec_i32_cnt_dosage2(const int32_t *p, int32_t *out, size_t n, int32_t val,
	int32_t missing, int32_t missing_substitute)
{
#ifdef COREARRAY_SIMD_SSE2

	// header 1, 16-byte aligned
	size_t h = ((16 - ((size_t)out & 0x0F)) & 0x0F) >> 2;
	for (; (n > 0) && (h > 0); n--, h--, p+=2)
	{
		*out ++ = ((p[0] == missing) || (p[1] == missing)) ?
			missing_substitute :
			(p[0]==val ? 1 : 0) + (p[1]==val ? 1 : 0);
	}

	// body, SSE2
	const __m128i val4  = _mm_set1_epi32(val);
	const __m128i miss4 = _mm_set1_epi32(missing);
	const __m128i sub4  = _mm_set1_epi32(missing_substitute);

#   ifdef COREARRAY_SIMD_AVX2

	// header 2, 32-byte aligned
	if ((n >= 4) && ((size_t)out & 0x10))
	{
		__m128i v, w;

		v = MM_LOADU_128((__m128i const*)p); p += 4;
		__m128i v1 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,2,0));
		__m128i v2 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,3,1));

		v = MM_LOADU_128((__m128i const*)p); p += 4;
		__m128i w1 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,2,0));
		__m128i w2 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,3,1));

		v1 = _mm_unpacklo_epi64(v1, w1);
		v2 = _mm_unpacklo_epi64(v2, w2);

		__m128i c = _mm_setzero_si128();
		c = _mm_sub_epi32(c, _mm_cmpeq_epi32(v1, val4));
		c = _mm_sub_epi32(c, _mm_cmpeq_epi32(v2, val4));

		w1 = _mm_cmpeq_epi32(v1, miss4);
		w2 = _mm_cmpeq_epi32(v2, miss4);
		w = _mm_or_si128(w1, w2);
		c = _mm_or_si128(_mm_and_si128(w, sub4), _mm_andnot_si128(w, c));

		_mm_store_si128((__m128i *)out, c);
		n -= 4; out += 4;
	}

	const __m256i val8  = _mm256_set1_epi32(val);
	const __m256i miss8 = _mm256_set1_epi32(missing);
	const __m256i sub8  = _mm256_set1_epi32(missing_substitute);

	const __m256i shuffle1 = _mm256_set_epi32(0, 0, 0, 0, 6, 4, 2, 0);
	const __m256i shuffle2 = _mm256_set_epi32(0, 0, 0, 0, 7, 5, 3, 1);

	for (; n >= 8; n-=8)
	{
		__m256i v, w;

		v = MM_LOADU_256((__m256i const*)p); p += 8;
		__m256i v1 = _mm256_permutevar8x32_epi32(v, shuffle1);
		__m256i v2 = _mm256_permutevar8x32_epi32(v, shuffle2);

		v = MM_LOADU_256((__m256i const*)p); p += 8;
		__m256i w1 = _mm256_permutevar8x32_epi32(v, shuffle1);
		__m256i w2 = _mm256_permutevar8x32_epi32(v, shuffle2);

		v1 = _mm256_permute2f128_si256(v1, w1, 0x20);
		v2 = _mm256_permute2f128_si256(v2, w2, 0x20);

		__m256i c = _mm256_setzero_si256();
		c = _mm256_sub_epi32(c, _mm256_cmpeq_epi32(v1, val8));
		c = _mm256_sub_epi32(c, _mm256_cmpeq_epi32(v2, val8));

		w1 = _mm256_cmpeq_epi32(v1, miss8);
		w2 = _mm256_cmpeq_epi32(v2, miss8);
		w = _mm256_or_si256(w1, w2);
		c = _mm256_or_si256(_mm256_and_si256(w, sub8), _mm256_andnot_si256(w, c));

		_mm256_store_si256((__m256i *)out, c);
		out += 8;
	}

#   endif

	for (; n >= 4; n-=4)
	{
		__m128i v, w;

		v = MM_LOADU_128((__m128i const*)p); p += 4;
		__m128i v1 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,2,0));
		__m128i v2 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,3,1));

		v = MM_LOADU_128((__m128i const*)p); p += 4;
		__m128i w1 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,2,0));
		__m128i w2 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,3,1));

		v1 = _mm_unpacklo_epi64(v1, w1);
		v2 = _mm_unpacklo_epi64(v2, w2);

		__m128i c = _mm_setzero_si128();
		c = _mm_sub_epi32(c, _mm_cmpeq_epi32(v1, val4));
		c = _mm_sub_epi32(c, _mm_cmpeq_epi32(v2, val4));

		w1 = _mm_cmpeq_epi32(v1, miss4);
		w2 = _mm_cmpeq_epi32(v2, miss4);
		w = _mm_or_si128(w1, w2);
		c = _mm_or_si128(_mm_and_si128(w, sub4), _mm_andnot_si128(w, c));

		_mm_store_si128((__m128i *)out, c);
		out += 4;
	}

#endif

	// tail
	for (; n > 0; n--, p+=2)
	{
		*out ++ = ((p[0] == missing) || (p[1] == missing)) ?
			missing_substitute :
			(p[0]==val ? 1 : 0) + (p[1]==val ? 1 : 0);
	}
}


/// assuming 'out' is 4-byte aligned, output (p[0]!=val) + (p[1]!=val) or missing_substitute
void vec_i32_cnt_dosage_alt2(const int32_t *p, int32_t *out, size_t n, int32_t val,
	int32_t missing, int32_t missing_substitute)
{
#ifdef COREARRAY_SIMD_SSE2

	// header 1, 16-byte aligned
	size_t h = ((16 - ((size_t)out & 0x0F)) & 0x0F) >> 2;
	for (; (n > 0) && (h > 0); n--, h--, p+=2)
	{
		*out ++ = ((p[0] == missing) || (p[1] == missing)) ?
			missing_substitute :
			(p[0]!=val ? 1 : 0) + (p[1]!=val ? 1 : 0);
	}

	// body, SSE2
	const __m128i val4  = _mm_set1_epi32(val);
	const __m128i miss4 = _mm_set1_epi32(missing);
	const __m128i sub4  = _mm_set1_epi32(missing_substitute);
	const __m128i two4  = _mm_set1_epi32(2);

#   ifdef COREARRAY_SIMD_AVX2

	// header 2, 32-byte aligned
	if ((n >= 4) && ((size_t)out & 0x10))
	{
		__m128i v, w;

		v = MM_LOADU_128((__m128i const*)p); p += 4;
		__m128i v1 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,2,0));
		__m128i v2 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,3,1));

		v = MM_LOADU_128((__m128i const*)p); p += 4;
		__m128i w1 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,2,0));
		__m128i w2 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,3,1));

		v1 = _mm_unpacklo_epi64(v1, w1);
		v2 = _mm_unpacklo_epi64(v2, w2);

		__m128i c = two4;
		c = _mm_add_epi32(c, _mm_cmpeq_epi32(v1, val4));
		c = _mm_add_epi32(c, _mm_cmpeq_epi32(v2, val4));

		w1 = _mm_cmpeq_epi32(v1, miss4);
		w2 = _mm_cmpeq_epi32(v2, miss4);
		w = _mm_or_si128(w1, w2);
		c = _mm_or_si128(_mm_and_si128(w, sub4), _mm_andnot_si128(w, c));

		_mm_store_si128((__m128i *)out, c);
		n -= 4; out += 4;
	}

	const __m256i val8  = _mm256_set1_epi32(val);
	const __m256i miss8 = _mm256_set1_epi32(missing);
	const __m256i sub8  = _mm256_set1_epi32(missing_substitute);
	const __m256i two8  = _mm256_set1_epi32(2);

	const __m256i shuffle1 = _mm256_set_epi32(0, 0, 0, 0, 6, 4, 2, 0);
	const __m256i shuffle2 = _mm256_set_epi32(0, 0, 0, 0, 7, 5, 3, 1);

	for (; n >= 8; n-=8)
	{
		__m256i v, w;

		v = MM_LOADU_256((__m256i const*)p); p += 8;
		__m256i v1 = _mm256_permutevar8x32_epi32(v, shuffle1);
		__m256i v2 = _mm256_permutevar8x32_epi32(v, shuffle2);

		v = MM_LOADU_256((__m256i const*)p); p += 8;
		__m256i w1 = _mm256_permutevar8x32_epi32(v, shuffle1);
		__m256i w2 = _mm256_permutevar8x32_epi32(v, shuffle2);

		v1 = _mm256_permute2f128_si256(v1, w1, 0x20);
		v2 = _mm256_permute2f128_si256(v2, w2, 0x20);

		__m256i c = two8;
		c = _mm256_add_epi32(c, _mm256_cmpeq_epi32(v1, val8));
		c = _mm256_add_epi32(c, _mm256_cmpeq_epi32(v2, val8));

		w1 = _mm256_cmpeq_epi32(v1, miss8);
		w2 = _mm256_cmpeq_epi32(v2, miss8);
		w = _mm256_or_si256(w1, w2);
		c = _mm256_or_si256(_mm256_and_si256(w, sub8), _mm256_andnot_si256(w, c));

		_mm256_store_si256((__m256i *)out, c);
		out += 8;
	}

#   endif

	for (; n >= 4; n-=4)
	{
		__m128i v, w;

		v = MM_LOADU_128((__m128i const*)p); p += 4;
		__m128i v1 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,2,0));
		__m128i v2 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,3,1));

		v = MM_LOADU_128((__m128i const*)p); p += 4;
		__m128i w1 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,2,0));
		__m128i w2 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,3,1));

		v1 = _mm_unpacklo_epi64(v1, w1);
		v2 = _mm_unpacklo_epi64(v2, w2);

		__m128i c = two4;
		c = _mm_add_epi32(c, _mm_cmpeq_epi32(v1, val4));
		c = _mm_add_epi32(c, _mm_cmpeq_epi32(v2, val4));

		w1 = _mm_cmpeq_epi32(v1, miss4);
		w2 = _mm_cmpeq_epi32(v2, miss4);
		w = _mm_or_si128(w1, w2);
		c = _mm_or_si128(_mm_and_si128(w, sub4), _mm_andnot_si128(w, c));

		_mm_store_si128((__m128i *)out, c);
		out += 4;
	}

#endif

	// tail
	for (; n > 0; n--, p+=2)
	{
		*out ++ = ((p[0] == missing) || (p[1] == missing)) ?
			missing_substitute :
			(p[0]!=val ? 1 : 0) + (p[1]!=val ? 1 : 0);
	}
}


/// shifting *p right by 2 bits, assuming p is 2-byte aligned
void vec_i32_shr_b2(int32_t *p, size_t n)
{
#ifdef COREARRAY_SIMD_SSE2

	// header 1, 16-byte aligned
	size_t h = ((16 - ((size_t)p & 0x0F)) & 0x0F) >> 2;
	for (; (n > 0) && (h > 0); n--, h--)
		*p++ >>= 2;

	// body, SSE2
	for (; n >= 4; n-=4, p+=4)
	{
		__m128i v = _mm_load_si128((__m128i const*)p);
		_mm_store_si128((__m128i *)p, _mm_srli_epi32(v, 2));
	}

#endif

	// tail
	for (; n > 0; n--) *p++ >>= 2;
}


/// bounds checking, excluding NA_INTEGER
int vec_i32_bound_check(const int32_t *p, size_t n, int bound)
{
	#define NA_INTEGER  0x80000000

#ifdef COREARRAY_SIMD_AVX2
	__m256i NA8   = _mm256_set1_epi32(NA_INTEGER);
	__m256i ZERO8 = _mm256_setzero_si256();
	__m256i BND8  = _mm256_set1_epi32(bound+1);
	for (; n >= 8; n-=8)
	{
		__m256i i8 = _mm256_loadu_si256((__m256i const*)p);
		p += 8;
		__m256i m = _mm256_and_si256(_mm256_cmpgt_epi32(i8, ZERO8),
			_mm256_cmpgt_epi32(BND8, i8));
		m = _mm256_or_si256(m, _mm256_cmpeq_epi32(i8, NA8));
		if (_mm256_movemask_epi8(m) != -1)
			return 0;
	}
#endif
#ifdef COREARRAY_SIMD_SSE2
	__m128i NA   = _mm_set1_epi32(NA_INTEGER);
	__m128i ZERO = _mm_setzero_si128();
	__m128i BND  = _mm_set1_epi32(bound+1);
	for (; n >= 4; n-=4)
	{
		__m128i i4 = _mm_loadu_si128((__m128i const*)p);
		p += 4;
		__m128i m = _mm_and_si128(_mm_cmplt_epi32(ZERO, i4), _mm_cmplt_epi32(i4, BND));
		m = _mm_or_si128(m, _mm_cmpeq_epi32(i4, NA));
		if (_mm_movemask_epi8(m) != 0xFFFF)
			return 0;
	}
#endif
	for (; n > 0; n--)
	{
		int i = *p++;
		if ((i != NA_INTEGER) && ((i < 1) || (i > bound)))
			return 0;
	}
	return -1;
}



// ===========================================================
// functions for char
// ===========================================================

const char *vec_char_find_CRLF(const char *p, size_t n)
{
#ifdef COREARRAY_SIMD_SSE2

	// header 1, 16-byte aligned
	size_t h = (16 - ((size_t)p & 0x0F)) & 0x0F;
	for (; (n > 0) && (h > 0); n--, h--, p++)
		if (*p=='\n' || *p=='\r') return p;

	// body, SSE2
	const __m128i mask1 = _mm_set1_epi8('\n');
	const __m128i mask2 = _mm_set1_epi8('\r');

#   ifdef COREARRAY_SIMD_AVX2

	// header 2, 32-byte aligned
	if ((n >= 16) && ((size_t)p & 0x10))
	{
		__m128i v  = _mm_load_si128((__m128i const*)p);
		__m128i c1 = _mm_cmpeq_epi8(v, mask1);
		__m128i c2 = _mm_cmpeq_epi8(v, mask2);
		if (_mm_movemask_epi8(_mm_or_si128(c1, c2)))
			goto tail;
		n -= 16; p += 16;
	}

	// body, AVX2
	const __m256i mask3 = _mm256_set1_epi8('\n');
	const __m256i mask4 = _mm256_set1_epi8('\r');

	for (; n >= 32; n-=32, p+=32)
	{
		__m256i v = _mm256_load_si256((__m256i const*)p);
		__m256i c1 = _mm256_cmpeq_epi8(v, mask3);
		__m256i c2 = _mm256_cmpeq_epi8(v, mask4);
		if (_mm256_movemask_epi8(_mm256_or_si256(c1, c2)))
			goto tail;
	}

#endif

	for (; n >= 16; n-=16, p+=16)
	{
		__m128i v  = _mm_load_si128((__m128i const*)p);
		__m128i c1 = _mm_cmpeq_epi8(v, mask1);
		__m128i c2 = _mm_cmpeq_epi8(v, mask2);
		if (_mm_movemask_epi8(_mm_or_si128(c1, c2)))
			break;
	}

#ifdef COREARRAY_SIMD_AVX2
tail:
#endif

#endif

	// tail
	for (; n > 0; n--, p++)
		if (*p=='\n' || *p=='\r') break;

	return p;
}


/// find non-zero
COREARRAY_DLL_DEFAULT const int8_t *vec_bool_find_true(const int8_t *p,
	const int8_t *end)
{
#ifdef COREARRAY_SIMD_AVX
	for (; p+32 < end; p+=32)
	{
		__m256i v = _mm256_loadu_si256((__m256i const*)p);
		if (!_mm256_testz_si256(v, v)) break;
	}
#endif
#ifdef COREARRAY_SIMD_SSE2
	#ifndef COREARRAY_SIMD_SSE4_1
	__m128i zero = _mm_setzero_si128();
	#endif
	for (; p+16 < end; p+=16)
	{
		__m128i v = _mm_loadu_si128((__m128i const*)p);
	#ifdef COREARRAY_SIMD_SSE4_1
		if (!_mm_testz_si128(v, v)) break;
	#else
		if (_mm_movemask_epi8(_mm_cmpeq_epi8(v, zero)) != 0xFFFF) break;
	#endif
	}
#endif
	for (; p < end; p++) if (*p) break;
	return p;
}
