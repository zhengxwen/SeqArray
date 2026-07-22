// ===========================================================
//
// vectorization.cpp: compiler optimization with vectorization
//
// Copyright (C) 2016-2026    Xiuwen Zheng
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
 *	\date     2016-2026
 *	\brief    compiler optimization with vectorization
 *	\details
**/

#ifndef COREARRAY_COMPILER_OPTIMIZE_FLAG
#   define COREARRAY_COMPILER_OPTIMIZE_FLAG  3
#endif

#include "vectorization.h"
#include <Rdefines.h>


#ifdef __cplusplus
extern "C" {
#endif


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
		ans += 16 - vec_sum_u8(bit);
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

		__m128i lo = _mm256_castsi256_si128(bit);
		__m128i hi = _mm256_extracti128_si256(bit, 1);
		ans += 256
			- POPCNT_U64((uint64_t)_mm_cvtsi128_si64(lo))
			- POPCNT_U64((uint64_t)_mm_cvtsi128_si64(_mm_srli_si128(lo, 8)))
			- POPCNT_U64((uint64_t)_mm_cvtsi128_si64(hi))
			- POPCNT_U64((uint64_t)_mm_cvtsi128_si64(_mm_srli_si128(hi, 8)));
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
		// 'bit' packs 8 comparison results per byte (one bit each), so the number
		// of zero bytes is the popcount of 'bit', not the sum of its byte values.
		ans += 128
			- POPCNT_U64((uint64_t)_mm_cvtsi128_si64(bit))
			- POPCNT_U64((uint64_t)_mm_cvtsi128_si64(_mm_srli_si128(bit, 8)));
	}

	for (; n >= 16; n -= 16)
	{
		__m128i c = _mm_cmpeq_epi8(_mm_load_si128((__m128i const*)p), ZERO);
		__m128i bit = _mm_and_si128(c, ONES);
		p += 16;
		ans += 16 - vec_sum_u8(bit);
	}

#elif defined(COREARRAY_SIMD_NEON)

	// body, NEON
	const uint8x16_t ZERO = vdupq_n_u8(0);
	uint8x16_t sum = vdupq_n_u8(0);
	size_t offset = 0;

	for (; n >= 64; n -= 64)
	{
		uint8x16_t v = vld1q_u8((const uint8_t*)p); p += 16;
		sum = vsubq_u8(sum, vmvnq_u8(vceqq_u8(v, ZERO)));
		v = vld1q_u8((const uint8_t*)p); p += 16;
		sum = vsubq_u8(sum, vmvnq_u8(vceqq_u8(v, ZERO)));
		v = vld1q_u8((const uint8_t*)p); p += 16;
		sum = vsubq_u8(sum, vmvnq_u8(vceqq_u8(v, ZERO)));
		v = vld1q_u8((const uint8_t*)p); p += 16;
		sum = vsubq_u8(sum, vmvnq_u8(vceqq_u8(v, ZERO)));
		offset += 4;
		if (offset >= 252)
		{
			ans += vec_neon_sum_u8(sum);
			sum = vdupq_n_u8(0);
			offset = 0;
		}
	}
	for (; n >= 16; n -= 16, p += 16)
	{
		uint8x16_t v = vld1q_u8((const uint8_t*)p);
		sum = vsubq_u8(sum, vmvnq_u8(vceqq_u8(v, ZERO)));
		if ((++offset) >= 252)
		{
			ans += vec_neon_sum_u8(sum);
			sum = vdupq_n_u8(0);
			offset = 0;
		}
	}
	if (offset > 0)
		ans += vec_neon_sum_u8(sum);

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

#elif defined(COREARRAY_SIMD_NEON)

	for (; n >= 16; n-=16, p+=16)
	{
		uint8x16_t v = vld1q_u8((const uint8_t*)p);
		if (vmaxvq_u8(v) != 0) break;
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
		offset++;
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
		offset++;
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

#elif defined(COREARRAY_SIMD_NEON)

	// body, NEON
	const uint8x16_t mask = vdupq_n_u8((uint8_t)val);
	uint8x16_t sum = vdupq_n_u8(0);
	size_t offset = 0;

	for (; n >= 64; n-=64)
	{
		uint8x16_t v = vld1q_u8((const uint8_t*)p); p += 16;
		sum = vsubq_u8(sum, vreinterpretq_u8_s8(
			vreinterpretq_s8_u8(vceqq_u8(v, mask))));
		v = vld1q_u8((const uint8_t*)p); p += 16;
		sum = vsubq_u8(sum, vreinterpretq_u8_s8(
			vreinterpretq_s8_u8(vceqq_u8(v, mask))));
		v = vld1q_u8((const uint8_t*)p); p += 16;
		sum = vsubq_u8(sum, vreinterpretq_u8_s8(
			vreinterpretq_s8_u8(vceqq_u8(v, mask))));
		v = vld1q_u8((const uint8_t*)p); p += 16;
		sum = vsubq_u8(sum, vreinterpretq_u8_s8(
			vreinterpretq_s8_u8(vceqq_u8(v, mask))));
		offset += 4;
		if (offset >= 252)
		{
			num += vec_neon_sum_u8(sum);
			sum = vdupq_n_u8(0);
			offset = 0;
		}
	}
	for (; n >= 16; n-=16, p+=16)
	{
		uint8x16_t v = vld1q_u8((const uint8_t*)p);
		sum = vsubq_u8(sum, vreinterpretq_u8_s8(
			vreinterpretq_s8_u8(vceqq_u8(v, mask))));
		if ((++offset) >= 252)
		{
			num += vec_neon_sum_u8(sum);
			sum = vdupq_n_u8(0);
			offset = 0;
		}
	}
	if (offset > 0)
		num += vec_neon_sum_u8(sum);

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
		offset++;
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
		offset++;
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

#elif defined(COREARRAY_SIMD_NEON)

	// body, NEON
	const uint8x16_t mask1 = vdupq_n_u8((uint8_t)val1);
	const uint8x16_t mask2 = vdupq_n_u8((uint8_t)val2);
	uint8x16_t sum1 = vdupq_n_u8(0), sum2 = vdupq_n_u8(0);
	size_t offset = 0;

	for (; n >= 16; n-=16, p+=16)
	{
		uint8x16_t v = vld1q_u8((const uint8_t*)p);
		sum1 = vsubq_u8(sum1, vceqq_u8(v, mask1));
		sum2 = vsubq_u8(sum2, vceqq_u8(v, mask2));
		if ((++offset) >= 252)
		{
			n1 += vec_neon_sum_u8(sum1);
			n2 += vec_neon_sum_u8(sum2);
			sum1 = sum2 = vdupq_n_u8(0);
			offset = 0;
		}
	}
	if (offset > 0)
	{
		n1 += vec_neon_sum_u8(sum1);
		n2 += vec_neon_sum_u8(sum2);
	}

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
		offset++;
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
		offset++;
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

#elif defined(COREARRAY_SIMD_NEON)

	// body, NEON
	const uint8x16_t mask1 = vdupq_n_u8((uint8_t)val1);
	const uint8x16_t mask2 = vdupq_n_u8((uint8_t)val2);
	const uint8x16_t mask3 = vdupq_n_u8((uint8_t)val3);
	uint8x16_t sum1 = vdupq_n_u8(0), sum2 = vdupq_n_u8(0), sum3 = vdupq_n_u8(0);
	size_t offset = 0;

	for (; n >= 16; n-=16, p+=16)
	{
		uint8x16_t v = vld1q_u8((const uint8_t*)p);
		sum1 = vsubq_u8(sum1, vceqq_u8(v, mask1));
		sum2 = vsubq_u8(sum2, vceqq_u8(v, mask2));
		sum3 = vsubq_u8(sum3, vceqq_u8(v, mask3));
		if ((++offset) >= 252)
		{
			n1 += vec_neon_sum_u8(sum1);
			n2 += vec_neon_sum_u8(sum2);
			n3 += vec_neon_sum_u8(sum3);
			sum1 = sum2 = sum3 = vdupq_n_u8(0);
			offset = 0;
		}
	}
	if (offset > 0)
	{
		n1 += vec_neon_sum_u8(sum1);
		n2 += vec_neon_sum_u8(sum2);
		n3 += vec_neon_sum_u8(sum3);
	}

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

#elif defined(COREARRAY_SIMD_NEON)

	// body, NEON
	const uint8x16_t vmask = vdupq_n_u8((uint8_t)val);
	const uint8x16_t vsub  = vdupq_n_u8((uint8_t)substitute);

	for (; n >= 16; n-=16, p+=16)
	{
		uint8x16_t v = vld1q_u8((const uint8_t*)p);
		uint8x16_t c = vceqq_u8(v, vmask);
		if (vmaxvq_u8(c))
			vst1q_u8((uint8_t*)p, vbslq_u8(c, vsub, v));
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

#elif defined(COREARRAY_SIMD_NEON)

	// body, NEON
	const uint8x16_t val16  = vdupq_n_u8((uint8_t)val);
	const uint8x16_t miss16 = vdupq_n_u8((uint8_t)missing);
	const uint8x16_t sub16  = vdupq_n_u8((uint8_t)missing_substitute);

	for (; n >= 16; n-=16)
	{
		uint8x16_t w1 = vld1q_u8((const uint8_t*)p); p += 16;
		uint8x16_t w2 = vld1q_u8((const uint8_t*)p); p += 16;

		// deinterleave: v1 gets even bytes, v2 gets odd bytes
		uint8x16x2_t dz = vuzpq_u8(w1, w2);
		uint8x16_t v1 = dz.val[0];
		uint8x16_t v2 = dz.val[1];

		uint8x16_t c = vdupq_n_u8(0);
		c = vsubq_u8(c, vceqq_u8(v1, val16));
		c = vsubq_u8(c, vceqq_u8(v2, val16));

		uint8x16_t m1 = vceqq_u8(v1, miss16);
		uint8x16_t m2 = vceqq_u8(v2, miss16);
		uint8x16_t w  = vorrq_u8(m1, m2);
		c = vbslq_u8(w, sub16, c);

		vst1q_u8((uint8_t*)out, c);
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

#elif defined(COREARRAY_SIMD_NEON)

	// body, NEON
	const uint8x16_t val16  = vdupq_n_u8((uint8_t)val);
	const uint8x16_t miss16 = vdupq_n_u8((uint8_t)missing);
	const uint8x16_t sub16  = vdupq_n_u8((uint8_t)missing_substitute);
	const uint8x16_t two16  = vdupq_n_u8(2);

	for (; n >= 16; n-=16)
	{
		uint8x16_t w1 = vld1q_u8((const uint8_t*)p); p += 16;
		uint8x16_t w2 = vld1q_u8((const uint8_t*)p); p += 16;

		uint8x16x2_t dz = vuzpq_u8(w1, w2);
		uint8x16_t v1 = dz.val[0];
		uint8x16_t v2 = dz.val[1];

		uint8x16_t c = two16;
		c = vaddq_u8(c, vceqq_u8(v1, val16));
		c = vaddq_u8(c, vceqq_u8(v2, val16));

		uint8x16_t m1 = vceqq_u8(v1, miss16);
		uint8x16_t m2 = vceqq_u8(v2, miss16);
		uint8x16_t w  = vorrq_u8(m1, m2);
		c = vbslq_u8(w, sub16, c);

		vst1q_u8((uint8_t*)out, c);
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


/// output (p[0]!=val) + (p[1]!=val) allowing partial missing
void vec_i8_cnt_dosage_alt2_p(const int8_t *p, int8_t *out, size_t n,
	int8_t val, int8_t missing, int8_t missing_substitute)
{
#ifdef COREARRAY_SIMD_SSE2

	// body, SSE2
	const __m128i val16  = _mm_set1_epi8(val);
	const __m128i miss16 = _mm_set1_epi8(missing);
	const __m128i sub16  = _mm_set1_epi8(missing_substitute);
	const __m128i two16  = _mm_set1_epi8(2);
	const __m128i mask   = _mm_set1_epi16(0x00FF);
	// SSE2 only
	for (; n >= 16; n-=16)
	{
		__m128i w0 = MM_LOADU_128((__m128i const*)p); p += 16;
		__m128i w1 = MM_LOADU_128((__m128i const*)p); p += 16;
		__m128i v0 = _mm_packus_epi16(_mm_and_si128(w0, mask), _mm_and_si128(w1, mask));
		__m128i v1 = _mm_packus_epi16(_mm_srli_epi16(w0, 8), _mm_srli_epi16(w1, 8));

		__m128i b0 = _mm_cmpeq_epi8(v0, miss16);
		__m128i b1 = _mm_cmpeq_epi8(v1, miss16);
		__m128i bb = _mm_and_si128(b0, b1);

		__m128i c = two16;
		c = _mm_add_epi8(c, _mm_or_si128(b0, _mm_cmpeq_epi8(v0, val16)));
		c = _mm_add_epi8(c, _mm_or_si128(b1, _mm_cmpeq_epi8(v1, val16)));
		c = _mm_or_si128(_mm_and_si128(bb, sub16), _mm_andnot_si128(bb, c));

		_mm_storeu_si128((__m128i *)out, c);
		out += 16;
	}

#elif defined(COREARRAY_SIMD_NEON)

	// body, NEON
	const uint8x16_t val16  = vdupq_n_u8((uint8_t)val);
	const uint8x16_t miss16 = vdupq_n_u8((uint8_t)missing);
	const uint8x16_t sub16  = vdupq_n_u8((uint8_t)missing_substitute);
	const uint8x16_t two16  = vdupq_n_u8(2);

	for (; n >= 16; n-=16)
	{
		uint8x16_t w0 = vld1q_u8((const uint8_t*)p); p += 16;
		uint8x16_t w1 = vld1q_u8((const uint8_t*)p); p += 16;

		uint8x16x2_t dz = vuzpq_u8(w0, w1);
		uint8x16_t v0 = dz.val[0];
		uint8x16_t v1 = dz.val[1];

		uint8x16_t b0 = vceqq_u8(v0, miss16);
		uint8x16_t b1 = vceqq_u8(v1, miss16);
		uint8x16_t bb = vandq_u8(b0, b1);

		uint8x16_t c = two16;
		c = vaddq_u8(c, vorrq_u8(b0, vceqq_u8(v0, val16)));
		c = vaddq_u8(c, vorrq_u8(b1, vceqq_u8(v1, val16)));
		c = vbslq_u8(bb, sub16, c);

		vst1q_u8((uint8_t*)out, c);
		out += 16;
	}

#endif
	// tail
	for (; n > 0; n--, p+=2)
	{
		const bool b0 = p[0]==missing, b1 = p[1]==missing;
		*out ++ = (b0 && b1) ? missing_substitute :
			(p[0]!=val && !b0) + (p[1]!=val && !b1);
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

#elif defined(COREARRAY_SIMD_NEON)

	// body, NEON
	for (; n >= 16; n-=16, p+=16)
	{
		uint8x16_t v = vld1q_u8(p);
		vst1q_u8(p, vshrq_n_u8(v, 2));
	}

#endif

	// tail
	for (; n > 0; n--) *p++ >>= 2;
}


/// *p |= (*s) << nbit
void vec_u8_or_shl(uint8_t *p, size_t n, const uint8_t *s, const uint8_t nbit)
{
#ifdef COREARRAY_SIMD_SSE2
	if (n >= 16)
	{
		const __m128i mask = _mm_set1_epi8(((uint8_t)0xFF) >> nbit);
		for (; n >= 16; n-=16, p+=16, s+=16)
		{
			__m128i pv = _mm_loadu_si128((__m128i const*)p);
			__m128i sv = _mm_loadu_si128((__m128i const*)s);
			__m128i  v = _mm_slli_epi16(sv & mask, nbit);
			_mm_storeu_si128((__m128i*)p, pv | v);
		}
	}
#elif defined(COREARRAY_SIMD_NEON)
	if (n >= 16)
	{
		const uint8x16_t vmask = vdupq_n_u8(((uint8_t)0xFF) >> nbit);
		const int8x16_t vshift = vdupq_n_s8((int8_t)nbit);
		for (; n >= 16; n-=16, p+=16, s+=16)
		{
			uint8x16_t pv = vld1q_u8(p);
			uint8x16_t sv = vld1q_u8(s);
			uint8x16_t v = vshlq_u8(vandq_u8(sv, vmask), vshift);
			vst1q_u8(p, vorrq_u8(pv, v));
		}
	}
#endif
	// tail
	for (; n > 0; n--) *p++ |= (*s++) << nbit;
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
		_mm_store_si128((__m128i *)p, _mm_srai_epi16(v, 2));
	}

#elif defined(COREARRAY_SIMD_NEON)

	// body, NEON
	for (; n >= 8; n-=8, p+=8)
	{
		int16x8_t v = vld1q_s16(p);
		vst1q_s16(p, vshrq_n_s16(v, 2));
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

#elif defined(COREARRAY_SIMD_NEON)

	// body, NEON
	const int32x4_t vmask = vdupq_n_s32(val);
	int32x4_t sum = vdupq_n_s32(0);
	for (; n >= 4; n-=4, p+=4)
	{
		int32x4_t v = vld1q_s32(p);
		uint32x4_t c = vceqq_s32(v, vmask);
		sum = vsubq_s32(sum, vreinterpretq_s32_u32(c));
	}
	ans += vec_neon_sum_i32(sum);

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

#elif defined(COREARRAY_SIMD_NEON)

	// body, NEON
	const int32x4_t vmask1 = vdupq_n_s32(val1);
	const int32x4_t vmask2 = vdupq_n_s32(val2);
	int32x4_t sum1 = vdupq_n_s32(0), sum2 = vdupq_n_s32(0);

	for (; n >= 4; n-=4, p+=4)
	{
		int32x4_t v = vld1q_s32(p);
		sum1 = vsubq_s32(sum1, vreinterpretq_s32_u32(vceqq_s32(v, vmask1)));
		sum2 = vsubq_s32(sum2, vreinterpretq_s32_u32(vceqq_s32(v, vmask2)));
	}
	n1 += vec_neon_sum_i32(sum1);
	n2 += vec_neon_sum_i32(sum2);

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

#elif defined(COREARRAY_SIMD_NEON)

	// body, NEON
	const int32x4_t vmask1 = vdupq_n_s32(val1);
	const int32x4_t vmask2 = vdupq_n_s32(val2);
	const int32x4_t vmask3 = vdupq_n_s32(val3);
	int32x4_t sum1 = vdupq_n_s32(0), sum2 = vdupq_n_s32(0), sum3 = vdupq_n_s32(0);

	for (; n >= 4; n-=4, p+=4)
	{
		int32x4_t v = vld1q_s32(p);
		sum1 = vsubq_s32(sum1, vreinterpretq_s32_u32(vceqq_s32(v, vmask1)));
		sum2 = vsubq_s32(sum2, vreinterpretq_s32_u32(vceqq_s32(v, vmask2)));
		sum3 = vsubq_s32(sum3, vreinterpretq_s32_u32(vceqq_s32(v, vmask3)));
	}
	n1 += vec_neon_sum_i32(sum1);
	n2 += vec_neon_sum_i32(sum2);
	n3 += vec_neon_sum_i32(sum3);

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

#elif defined(COREARRAY_SIMD_NEON)

	// body, NEON
	const int32x4_t vmask = vdupq_n_s32(val);
	const int32x4_t vsub  = vdupq_n_s32(substitute);

	for (; n >= 4; n-=4, p+=4)
	{
		int32x4_t v = vld1q_s32(p);
		uint32x4_t c = vceqq_s32(v, vmask);
		if (vmaxvq_u32(c))
			vst1q_s32(p, vbslq_s32(c, vsub, v));
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

#elif defined(COREARRAY_SIMD_NEON)

	// body, NEON - deinterleave pairs: p[0],p[1] → v1[i],v2[i]
	const int32x4_t val4  = vdupq_n_s32(val);
	const int32x4_t miss4 = vdupq_n_s32(missing);
	const int32x4_t sub4  = vdupq_n_s32(missing_substitute);

	for (; n >= 4; n-=4)
	{
		int32x4x2_t pair = vld2q_s32(p); p += 8;
		int32x4_t v1 = pair.val[0];
		int32x4_t v2 = pair.val[1];

		int32x4_t c = vdupq_n_s32(0);
		c = vsubq_s32(c, vreinterpretq_s32_u32(vceqq_s32(v1, val4)));
		c = vsubq_s32(c, vreinterpretq_s32_u32(vceqq_s32(v2, val4)));

		uint32x4_t m1 = vceqq_s32(v1, miss4);
		uint32x4_t m2 = vceqq_s32(v2, miss4);
		uint32x4_t w  = vorrq_u32(m1, m2);
		c = vbslq_s32(w, sub4, c);

		vst1q_s32(out, c);
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

#elif defined(COREARRAY_SIMD_NEON)

	// body, NEON
	const int32x4_t val4  = vdupq_n_s32(val);
	const int32x4_t miss4 = vdupq_n_s32(missing);
	const int32x4_t sub4  = vdupq_n_s32(missing_substitute);
	const int32x4_t two4  = vdupq_n_s32(2);

	for (; n >= 4; n-=4)
	{
		int32x4x2_t pair = vld2q_s32(p); p += 8;
		int32x4_t v1 = pair.val[0];
		int32x4_t v2 = pair.val[1];

		int32x4_t c = two4;
		c = vaddq_s32(c, vreinterpretq_s32_u32(vceqq_s32(v1, val4)));
		c = vaddq_s32(c, vreinterpretq_s32_u32(vceqq_s32(v2, val4)));

		uint32x4_t m1 = vceqq_s32(v1, miss4);
		uint32x4_t m2 = vceqq_s32(v2, miss4);
		uint32x4_t w  = vorrq_u32(m1, m2);
		c = vbslq_s32(w, sub4, c);

		vst1q_s32(out, c);
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


/// output (p[0]!=val) + (p[1]!=val) allowing partial missing
COREARRAY_DLL_DEFAULT void vec_i32_cnt_dosage_alt2_p(const int32_t *p,
	int32_t *out, size_t n, int32_t val, int32_t missing,
	int32_t missing_substitute)
{
#ifdef COREARRAY_SIMD_SSE2

	// body, SSE2
	const __m128i val4  = _mm_set1_epi32(val);
	const __m128i miss4 = _mm_set1_epi32(missing);
	const __m128i sub4  = _mm_set1_epi32(missing_substitute);
	const __m128i two4  = _mm_set1_epi32(2);
	for (; n >= 4; n-=4)
	{
		__m128i v;
		v = MM_LOADU_128((__m128i const*)p); p += 4;
		__m128i v0 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,2,0));
		__m128i v1 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,3,1));

		v = MM_LOADU_128((__m128i const*)p); p += 4;
		__m128i w0 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,2,0));
		__m128i w1 = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,0,3,1));

		v0 = _mm_unpacklo_epi64(v0, w0);
		v1 = _mm_unpacklo_epi64(v1, w1);
		__m128i b0 = _mm_cmpeq_epi32(v0, miss4);
		__m128i b1 = _mm_cmpeq_epi32(v1, miss4);
		__m128i bb = _mm_and_si128(b0, b1);

		__m128i c = two4;
		c = _mm_add_epi32(c, _mm_or_si128(b0, _mm_cmpeq_epi32(v0, val4)));
		c = _mm_add_epi32(c, _mm_or_si128(b1, _mm_cmpeq_epi32(v1, val4)));
		c = _mm_or_si128(_mm_and_si128(bb, sub4), _mm_andnot_si128(bb, c));

		_mm_storeu_si128((__m128i *)out, c);
		out += 4;
	}

#elif defined(COREARRAY_SIMD_NEON)

	// body, NEON
	const int32x4_t val4  = vdupq_n_s32(val);
	const int32x4_t miss4 = vdupq_n_s32(missing);
	const int32x4_t sub4  = vdupq_n_s32(missing_substitute);
	const int32x4_t two4  = vdupq_n_s32(2);

	for (; n >= 4; n-=4)
	{
		int32x4x2_t pair = vld2q_s32(p); p += 8;
		int32x4_t v0 = pair.val[0];
		int32x4_t v1 = pair.val[1];

		uint32x4_t b0 = vceqq_s32(v0, miss4);
		uint32x4_t b1 = vceqq_s32(v1, miss4);
		uint32x4_t bb = vandq_u32(b0, b1);

		int32x4_t c = two4;
		c = vaddq_s32(c, vreinterpretq_s32_u32(
			vorrq_u32(b0, vceqq_s32(v0, val4))));
		c = vaddq_s32(c, vreinterpretq_s32_u32(
			vorrq_u32(b1, vceqq_s32(v1, val4))));
		c = vbslq_s32(bb, sub4, c);

		vst1q_s32(out, c);
		out += 4;
	}

#endif

	// tail
	for (; n > 0; n--, p+=2)
	{
		const bool b0 = p[0]==missing, b1 = p[1]==missing;
		*out ++ = (b0 && b1) ? missing_substitute :
			(p[0]!=val && !b0) + (p[1]!=val && !b1);
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
		_mm_store_si128((__m128i *)p, _mm_srai_epi32(v, 2));
	}

#elif defined(COREARRAY_SIMD_NEON)

	// body, NEON
	for (; n >= 4; n-=4, p+=4)
	{
		int32x4_t v = vld1q_s32(p);
		vst1q_s32(p, vshrq_n_s32(v, 2));
	}

#endif

	// tail
	for (; n > 0; n--) *p++ >>= 2;
}


/// bounds checking, excluding NA_INTEGER
int vec_i32_bound_check(const int32_t *p, size_t n, int bound)
{
#ifdef COREARRAY_SIMD_AVX2
	__m256i NA8   = _mm256_set1_epi32(NA_INTEGER);
	__m256i ZERO8 = _mm256_setzero_si256();
	__m256i BND8  = _mm256_set1_epi32(bound);
	for (; n >= 8; n-=8)
	{
		__m256i i8 = _mm256_loadu_si256((__m256i const*)p);
		p += 8;
		__m256i m = _mm256_andnot_si256(_mm256_cmpgt_epi32(i8, BND8),
			_mm256_cmpgt_epi32(i8, ZERO8));
		m = _mm256_or_si256(m, _mm256_cmpeq_epi32(i8, NA8));
		if (_mm256_movemask_epi8(m) != -1)
			return 0;
	}
#endif
#ifdef COREARRAY_SIMD_SSE2
	__m128i NA   = _mm_set1_epi32(NA_INTEGER);
	__m128i ZERO = _mm_setzero_si128();
	__m128i BND  = _mm_set1_epi32(bound);
	for (; n >= 4; n-=4)
	{
		__m128i i4 = _mm_loadu_si128((__m128i const*)p);
		p += 4;
		__m128i m = _mm_andnot_si128(_mm_cmpgt_epi32(i4, BND), _mm_cmpgt_epi32(i4, ZERO));
		m = _mm_or_si128(m, _mm_cmpeq_epi32(i4, NA));
		if (_mm_movemask_epi8(m) != 0xFFFF)
			return 0;
	}
#elif defined(COREARRAY_SIMD_NEON)
	const int32x4_t vNA   = vdupq_n_s32(NA_INTEGER);
	const int32x4_t vZERO = vdupq_n_s32(0);
	const int32x4_t vBND  = vdupq_n_s32(bound);
	for (; n >= 4; n-=4)
	{
		int32x4_t i4 = vld1q_s32(p);
		p += 4;
		// m = (i4 > 0) & !(i4 > bound) = valid and in-range
		uint32x4_t gt_zero = vcgtq_s32(i4, vZERO);
		uint32x4_t gt_bnd  = vcgtq_s32(i4, vBND);
		uint32x4_t m = vbicq_u32(gt_zero, gt_bnd);  // in-range: >0 and <=bound
		m = vorrq_u32(m, vceqq_s32(i4, vNA));       // or is NA
		if (vminvq_u32(m) == 0)
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


/// *p |= (*s) << nbit
COREARRAY_DLL_DEFAULT void vec_i32_or_shl(int32_t *p, size_t n,
	const int32_t *s, const uint8_t nbit)
{
#ifdef COREARRAY_SIMD_SSE2
	for (; n >= 4; n-=4, p+=4, s+=4)
	{
		__m128i pv = _mm_loadu_si128((__m128i const*)p);
		__m128i sv = _mm_loadu_si128((__m128i const*)s);
		__m128i  v = _mm_slli_epi32(sv, nbit);
		_mm_storeu_si128((__m128i*)p, pv | v);
	}
#elif defined(COREARRAY_SIMD_NEON)
	const int32x4_t vshift = vdupq_n_s32((int32_t)nbit);
	for (; n >= 4; n-=4, p+=4, s+=4)
	{
		int32x4_t pv = vld1q_s32(p);
		int32x4_t sv = vld1q_s32(s);
		int32x4_t v  = vshlq_s32(sv, vshift);
		vst1q_s32(p, vorrq_s32(pv, v));
	}
#endif
	// tail
	for (; n > 0; n--) *p++ |= (*s++) << nbit;
}


/// *p |= (*s) << nbit
COREARRAY_DLL_DEFAULT void vec_i32_or_shl2(int32_t *p, size_t n,
	const uint8_t *s, const uint8_t nbit)
{
#ifdef COREARRAY_SIMD_SSE2
	if (n >= 16)
	{
		const __m128i zero = _mm_setzero_si128();
		const __m128i mask = _mm_set1_epi8(((uint8_t)0xFF) >> nbit);
		for (; n >= 16; n-=16, s+=16)
		{
			__m128i sv4 = _mm_loadu_si128((__m128i const*)s) & mask;
			// low 8 uint8
			__m128i s_l = _mm_unpacklo_epi8(sv4, zero); // uint8 => uint16 (low)
			{ 	// 1st 4 int
				__m128i s = _mm_unpacklo_epi16(s_l, zero); // uint16 => uint32
				__m128i pv = _mm_loadu_si128((__m128i const*)p);
				_mm_storeu_si128((__m128i*)p, pv | _mm_slli_epi32(s, nbit));
				p += 4;
			}
			{ 	// 2nd 4 int
				__m128i s = _mm_unpackhi_epi16(s_l, zero); // uint16 => uint32
				__m128i pv = _mm_loadu_si128((__m128i const*)p);
				_mm_storeu_si128((__m128i*)p, pv | _mm_slli_epi32(s, nbit));
				p += 4;
			}
			// high 8 uint8
			__m128i s_h = _mm_unpackhi_epi8(sv4, zero); // uint8 => uint16 (high)
			{ 	// 1st 4 int
				__m128i s = _mm_unpacklo_epi16(s_h, zero); // uint16 => uint32
				__m128i pv = _mm_loadu_si128((__m128i const*)p);
				_mm_storeu_si128((__m128i*)p, pv | _mm_slli_epi32(s, nbit));
				p += 4;
			}
			{ 	// 2nd 4 int
				__m128i s = _mm_unpackhi_epi16(s_h, zero); // uint16 => uint32
				__m128i pv = _mm_loadu_si128((__m128i const*)p);
				_mm_storeu_si128((__m128i*)p, pv | _mm_slli_epi32(s, nbit));
				p += 4;
			}
		}
	}
#elif defined(COREARRAY_SIMD_NEON)
	if (n >= 16)
	{
		const int32x4_t vshift = vdupq_n_s32((int32_t)nbit);
		const uint8x16_t vmask = vdupq_n_u8(((uint8_t)0xFF) >> nbit);
		for (; n >= 16; n-=16, s+=16)
		{
			uint8x16_t sv16 = vandq_u8(vld1q_u8(s), vmask);
			// widen uint8 to uint16
			uint16x8_t s_l = vmovl_u8(vget_low_u8(sv16));
			uint16x8_t s_h = vmovl_u8(vget_high_u8(sv16));
			// widen uint16 to uint32 and shift-or
			{	// 1st 4 int
				int32x4_t sv = vreinterpretq_s32_u32(vmovl_u16(vget_low_u16(s_l)));
				int32x4_t pv = vld1q_s32(p);
				vst1q_s32(p, vorrq_s32(pv, vshlq_s32(sv, vshift)));
				p += 4;
			}
			{	// 2nd 4 int
				int32x4_t sv = vreinterpretq_s32_u32(vmovl_u16(vget_high_u16(s_l)));
				int32x4_t pv = vld1q_s32(p);
				vst1q_s32(p, vorrq_s32(pv, vshlq_s32(sv, vshift)));
				p += 4;
			}
			{	// 3rd 4 int
				int32x4_t sv = vreinterpretq_s32_u32(vmovl_u16(vget_low_u16(s_h)));
				int32x4_t pv = vld1q_s32(p);
				vst1q_s32(p, vorrq_s32(pv, vshlq_s32(sv, vshift)));
				p += 4;
			}
			{	// 4th 4 int
				int32x4_t sv = vreinterpretq_s32_u32(vmovl_u16(vget_high_u16(s_h)));
				int32x4_t pv = vld1q_s32(p);
				vst1q_s32(p, vorrq_s32(pv, vshlq_s32(sv, vshift)));
				p += 4;
			}
		}
	}
#endif
	// tail
	for (; n > 0; n--) *p++ |= (*s++) << nbit;
}


// ===========================================================
// functions for float64
// ===========================================================

/// get the number of finite numbers
size_t vec_f64_num_notfinite(const double p[], size_t n)
{
	size_t m = 0;
	for (size_t i=0; i < n; i++)
		if (!R_FINITE(p[i])) m ++;
	return m;
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

#ifdef COREARRAY_SIMD_AVX2

	// header 2, 32-byte aligned
	if ((n >= 16) && ((size_t)p & 0x10))
	{
		__m128i v  = _mm_load_si128((__m128i const*)p);
		__m128i c1 = _mm_cmpeq_epi8(v, mask1);
		__m128i c2 = _mm_cmpeq_epi8(v, mask2);
		if (_mm_movemask_epi8(_mm_or_si128(c1, c2)))
		{
			for (; n > 0; n--, p++)
				if (*p=='\n' || *p=='\r') break;
			return p;
		}
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
		{
			for (; n > 0; n--, p++)
				if (*p=='\n' || *p=='\r') break;
			return p;
		}
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

#elif defined(COREARRAY_SIMD_NEON)

	// body, NEON
	const uint8x16_t vmask1 = vdupq_n_u8('\n');
	const uint8x16_t vmask2 = vdupq_n_u8('\r');
	for (; n >= 16; n-=16, p+=16)
	{
		uint8x16_t v  = vld1q_u8((const uint8_t*)p);
		uint8x16_t c1 = vceqq_u8(v, vmask1);
		uint8x16_t c2 = vceqq_u8(v, vmask2);
		if (vmaxvq_u8(vorrq_u8(c1, c2)))
			break;
	}

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
#elif defined(COREARRAY_SIMD_NEON)
	for (; p+16 < end; p+=16)
	{
		uint8x16_t v = vld1q_u8((const uint8_t*)p);
		if (vmaxvq_u8(v) != 0) break;
	}
#endif
	for (; p < end; p++) if (*p) break;
	return p;
}


#ifdef __cplusplus
}
#endif
