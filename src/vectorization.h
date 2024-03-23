// ===========================================================
//
// vectorization.h: compiler optimization with vectorization
//
// Copyright (C) 2016-2024    Xiuwen Zheng
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
 *	\file     vectorization.h
 *	\author   Xiuwen Zheng [zhengx@u.washington.edu]
 *	\version  1.0
 *	\date     2016-2024
 *	\brief    compiler optimization with vectorization
 *	\details
**/


#ifndef _HEADER_COREARRAY_VECTORIZATION_
#define _HEADER_COREARRAY_VECTORIZATION_

#include <CoreDEF.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>


#if defined(COREARRAY_SIMD_SSE) && defined(COREARRAY_SIMD_SSE2)

#   include <xmmintrin.h>  // SSE
#   include <emmintrin.h>  // SSE2

#   if defined(COREARRAY_SIMD_SSE3)  // SSE3
#       include <pmmintrin.h>
#   endif

#   if defined(COREARRAY_SIMD_SSSE3)  // SSSE3
#       include <tmmintrin.h>
#   endif

#   if defined(COREARRAY_SIMD_SSE4_1)  // SSE_4_1
#       include <smmintrin.h>
#   endif

#   if defined(COREARRAY_SIMD_SSE4_2) || defined(__POPCNT__)
#       define COREARRAY_HARDWARE_POPCNT
#       include <nmmintrin.h>  // COREARRAY_SIMD_SSE4_2, for POPCNT
#   endif

#   if defined(COREARRAY_SIMD_AVX) || defined(COREARRAY_SIMD_AVX2)
#       include <immintrin.h>  // AVX, AVX2
#   endif

#endif



#ifdef __cplusplus

namespace Vectorization
{
	/// an aligned pointer
	struct COREARRAY_DLL_DEFAULT ALIGN_PTR
	{
	private:
		void *alloc_ptr, *base_ptr;

	public:
		ALIGN_PTR()
			{ alloc_ptr = base_ptr = NULL; }
		ALIGN_PTR(size_t n, size_t align)
			{
				alloc_ptr = base_ptr = NULL;
				reset(n, align);
			}
		~ALIGN_PTR()
			{
				if (alloc_ptr) free(alloc_ptr);
				alloc_ptr = base_ptr = NULL;
			}

		void reset(size_t n, size_t align)
		{
			if (n > 0)
			{
				if (align < 1) align = 1;
				alloc_ptr = realloc(alloc_ptr, n + align - 1);
				if (!alloc_ptr)
					throw "Insufficient memory.";
				size_t r = ((size_t)alloc_ptr) % align;
				base_ptr = r ? (void*)((char*)alloc_ptr + align - r) : alloc_ptr;
			} else {
				if (alloc_ptr) free(alloc_ptr);
				alloc_ptr = base_ptr = NULL;
			}
		}

		inline void *get() { return base_ptr; }
	};

	/// an aligned pointer with 16-byte alignment
	struct COREARRAY_DLL_DEFAULT ALIGN_PTR_SSE: public ALIGN_PTR
	{
	public:
		ALIGN_PTR_SSE(): ALIGN_PTR()
			{ }
		ALIGN_PTR_SSE(size_t n): ALIGN_PTR(n, 16u)
			{ }
		void reset(size_t n)
			{ ALIGN_PTR::reset(n, 16u); }
	};

	/// an aligned pointer with 32-byte alignment
	struct COREARRAY_DLL_DEFAULT ALIGN_PTR_AVX: public ALIGN_PTR
	{
	public:
		ALIGN_PTR_AVX(): ALIGN_PTR()
			{ }
		ALIGN_PTR_AVX(size_t n): ALIGN_PTR(n, 32u)
			{ }
		void reset(size_t n)
			{ ALIGN_PTR::reset(n, 32u); }
	};


	// auto pointer

#if defined(COREARRAY_SIMD_AVX)

	typedef ALIGN_PTR_AVX VEC_AUTO_PTR;

#elif defined(COREARRAY_SIMD_SSE)

	typedef ALIGN_PTR_SSE VEC_AUTO_PTR;

#else
	struct COREARRAY_DLL_DEFAULT VEC_AUTO_PTR: public ALIGN_PTR
	{
	public:
		VEC_AUTO_PTR(): ALIGN_PTR()
			{ }
		VEC_AUTO_PTR(size_t n): ALIGN_PTR(n, 1u)
			{ }
		void reset(size_t n)
			{ ALIGN_PTR::reset(n, 1u); }
	};
#endif
}

#endif



#ifdef __cplusplus
extern "C" {
#endif


#ifdef COREARRAY_HARDWARE_POPCNT

#   define POPCNT_U32(x)    _mm_popcnt_u32((uint32_t)(x))
#   ifdef COREARRAY_REGISTER_BIT32
#       define POPCNT_U64(x)    \\
#           _mm_popcnt_u32((uint32_t)(x)) + _mm_popcnt_u32((uint64_t)(x) >> 32)
#   else
#       define POPCNT_U64(x)    _mm_popcnt_u64((uint64_t)(x))
#   endif

#else

inline static int POPCNT_U32(uint32_t x)
{
	x = x - ((x >> 1) & 0x55555555);
	x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
	return (((x + (x >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24;
}

inline static int POPCNT_U64(uint64_t x)
{
	x -= ((x >> 1) & 0x5555555555555555LLU);
	x = (x & 0x3333333333333333LLU) + ((x >> 2) & 0x3333333333333333LLU);
	x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FLLU;
	return (x * 0x0101010101010101LLU) >> 56;
}

#endif



// ===========================================================

#ifdef COREARRAY_SIMD_SSE

#   ifdef COREARRAY_SIMD_SSE3
#       define MM_LOADU_128(p)    _mm_lddqu_si128((__m128i const*)(p))
#   else
#       define MM_LOADU_128(p)    _mm_loadu_si128((__m128i const*)(p))
#   endif

#   ifdef COREARRAY_SIMD_SSE2
#       define MM_BLEND_128(a, b, mask)  \
		    _mm_or_si128(_mm_and_si128(mask, a), _mm_andnot_si128(mask, b))
#   endif

#endif


#ifdef COREARRAY_SIMD_AVX

#   define MM_LOADU_256(p)    _mm256_loadu_si256((__m256i const *)(p))
#   define MM_SET_M128(v1, v0)    \
        _mm256_insertf128_si256(_mm256_castsi128_si256(v0), (v1), 1)

#   ifdef COREARRAY_SIMD_AVX2
#       define MM_BLEND_256(a, b, mask)  \
		    _mm256_or_si256(_mm256_and_si256(mask, a), _mm256_andnot_si256(mask, b))
#   endif

#endif



// ===========================================================
// Sum all elements in a SIMD register
// ===========================================================

#ifdef COREARRAY_SIMD_SSE2
	inline static double vec_sum_f64(__m128d s)
	{
		return _mm_cvtsd_f64(_mm_add_pd(s, _mm_shuffle_pd(s, s, 1)));
	}
	inline static int vec_sum_i32(__m128i s)
	{
		s = _mm_add_epi32(s, _mm_shuffle_epi32(s, _MM_SHUFFLE(1,0,3,2)));
		s = _mm_add_epi32(s, _mm_shuffle_epi32(s, _MM_SHUFFLE(0,0,0,1)));
		return _mm_cvtsi128_si32(s);
	}
	inline static int vec_sum_u8(__m128i s)
	{
		s = _mm_sad_epu8(s, _mm_setzero_si128());
		s = _mm_add_epi32(s, _mm_shuffle_epi32(s, 2));
		return _mm_cvtsi128_si32(s);
	}
#endif

#ifdef COREARRAY_SIMD_AVX
	inline static double vec_avx_sum_f64(__m256d s)
	{
		s = _mm256_add_pd(_mm256_permute_pd(s, 5), s);
		s = _mm256_add_pd(s, _mm256_permute2f128_pd(s, s, 1));
		return _mm_cvtsd_f64(_mm256_castpd256_pd128(s));
	}
#endif

#ifdef COREARRAY_SIMD_AVX2
	inline static int vec_avx_sum_i32(__m256i s)
	{
		s = _mm256_hadd_epi32(s, s);
		s = _mm256_add_epi32(s, _mm256_permute4x64_epi64(s, _MM_SHUFFLE(1,0,3,2)));
		__m128i a = _mm256_castsi256_si128(s);
		a = _mm_add_epi32(a, _mm_shuffle_epi32(a, _MM_SHUFFLE(0,0,0,1)));
		return _mm_cvtsi128_si32(a);
	}
	inline static int vec_avx_sum_u8(__m256i s)
	{
		s = _mm256_sad_epu8(s, _mm256_setzero_si256());
		s = _mm256_add_epi64(s, _mm256_permute4x64_epi64(s, _MM_SHUFFLE(1,0,3,2)));
		return _mm256_extract_epi32(s,0) + _mm256_extract_epi32(s,2);
	}
#endif



/// get the number of non-zeros
COREARRAY_DLL_DEFAULT size_t vec_i8_cnt_nonzero(const int8_t *p, size_t n);

/// get the number of non-zeros and the pointer to the first non-zero value
COREARRAY_DLL_DEFAULT const int8_t *vec_i8_cnt_nonzero_ptr(const int8_t *p,
	size_t n, size_t *out_n);



// ===========================================================
// functions for int8
// ===========================================================

/// count how many 'val' in 'p'
COREARRAY_DLL_DEFAULT size_t vec_i8_count(const char *p, size_t n, char val);

/// count how many val1 and val2 in p
COREARRAY_DLL_DEFAULT void vec_i8_count2(const char *p, size_t n,
	char val1, char val2, size_t *out_n1, size_t *out_n2);

/// count how many val1, val2 and val3 in p
COREARRAY_DLL_DEFAULT void vec_i8_count3(const char *p, size_t n,
	char val1, char val2, char val3, size_t *out_n1, size_t *out_n2,
	size_t *out_n3);

/// replace 'val' in the array of 'p' by 'substitute'
COREARRAY_DLL_DEFAULT void vec_i8_replace(int8_t *p, size_t n, int8_t val,
	int8_t substitute);

/// output (p[0]==val) + (p[1]==val) or missing_substitute
COREARRAY_DLL_DEFAULT void vec_i8_cnt_dosage2(const int8_t *p,
	int8_t *out, size_t n, int8_t val, int8_t missing, int8_t missing_substitute);

/// output (p[0]!=val) + (p[1]!=val) or missing_substitute
COREARRAY_DLL_DEFAULT void vec_i8_cnt_dosage_alt2(const int8_t *p,
	int8_t *out, size_t n, int8_t val, int8_t missing, int8_t missing_substitute);

/// output (p[0]!=val) + (p[1]!=val) allowing partial missing
COREARRAY_DLL_DEFAULT void vec_i8_cnt_dosage_alt2_p(const int8_t *p,
	int8_t *out, size_t n, int8_t val, int8_t missing, int8_t missing_substitute);



// ===========================================================
// functions for uint8
// ===========================================================

/// shifting *p right by 2 bits, assuming p is 2-byte aligned
COREARRAY_DLL_DEFAULT void vec_u8_shr_b2(uint8_t *p, size_t n);

/// *p |= (*s) << nbit
COREARRAY_DLL_DEFAULT void vec_u8_or_shl(uint8_t *p, size_t n,
	const uint8_t *s, const uint8_t nbit);



// ===========================================================
// functions for int16
// ===========================================================

/// shifting *p right by 2 bits, assuming p is 2-byte aligned
COREARRAY_DLL_DEFAULT void vec_i16_shr_b2(int16_t *p, size_t n);



// ===========================================================
// functions for int32
// ===========================================================

/// count how many val in p, assuming p is 4-byte aligned
COREARRAY_DLL_DEFAULT size_t vec_i32_count(const int32_t *p, size_t n, int32_t val);

/// count how many val1 and val2 in p, assuming p is 4-byte aligned
COREARRAY_DLL_DEFAULT void vec_i32_count2(const int32_t *p, size_t n,
	int32_t val1, int32_t val2, size_t *out_n1, size_t *out_n2);

/// count how many val1, val2 and val3 in p, assuming p is 4-byte aligned
COREARRAY_DLL_DEFAULT void vec_i32_count3(const int32_t *p, size_t n,
	int32_t val1, int32_t val2, int32_t val3, size_t *out_n1, size_t *out_n2,
	size_t *out_n3);

///
COREARRAY_DLL_DEFAULT void vec_int32_set(int32_t *p, size_t n, int32_t val);

/// replace 'val' in the array of 'p' by 'substitute', assuming 'p' is 4-byte aligned
COREARRAY_DLL_DEFAULT void vec_i32_replace(int32_t *p, size_t n, int32_t val,
	int32_t substitute);

/// assuming 'out' is 4-byte aligned, output (p[0]==val) + (p[1]==val) or missing_substitute
COREARRAY_DLL_DEFAULT void vec_i32_cnt_dosage2(const int32_t *p,
	int32_t *out, size_t n, int32_t val, int32_t missing, int32_t missing_substitute);

/// assuming 'out' is 4-byte aligned, output (p[0]!=val) + (p[1]!=val) or missing_substitute
COREARRAY_DLL_DEFAULT void vec_i32_cnt_dosage_alt2(const int32_t *p,
	int32_t *out, size_t n, int32_t val, int32_t missing, int32_t missing_substitute);

/// output (p[0]!=val) + (p[1]!=val) allowing partial missing
COREARRAY_DLL_DEFAULT void vec_i32_cnt_dosage_alt2_p(const int32_t *p,
	int32_t *out, size_t n, int32_t val, int32_t missing, int32_t missing_substitute);

/// shifting *p right by 2 bits, assuming p is 4-byte aligned
COREARRAY_DLL_DEFAULT void vec_i32_shr_b2(int32_t *p, size_t n);

/// bounds checking, return 0 if fails
COREARRAY_DLL_DEFAULT int vec_i32_bound_check(const int32_t *p, size_t n, int bound);

/// *p |= (*s) << nbit
COREARRAY_DLL_DEFAULT void vec_i32_or_shl(int32_t *p, size_t n,
	const int32_t *s, const uint8_t nbit);

/// *p |= (*s) << nbit
COREARRAY_DLL_DEFAULT void vec_i32_or_shl2(int32_t *p, size_t n,
	const uint8_t *s, const uint8_t nbit);



// ===========================================================
// functions for float64
// ===========================================================

/// get the number of finite numbers
COREARRAY_DLL_DEFAULT size_t vec_f64_num_notfinite(const double p[], size_t n);


// ===========================================================
// functions for searching characters
// ===========================================================

#define  VEC_BOOL_FIND_TRUE(p, end)    \
	(C_BOOL*)vec_bool_find_true((int8_t*)(p), (int8_t*)(end))

/// find CRLF character
COREARRAY_DLL_DEFAULT const char *vec_char_find_CRLF(const char *p, size_t n);

/// find non-zero
COREARRAY_DLL_DEFAULT const int8_t *vec_bool_find_true(const int8_t *p,
	const int8_t *end);


#ifdef __cplusplus
}
#endif

#endif /* _HEADER_COREARRAY_VECTORIZATION_ */
